"""Well-balanced first-order shallow water: HLL fluxes with Audusse hydrostatic reconstruction."""
from __future__ import annotations

import math

import numpy as np

from ..base import G, WaterSolver

H_DRY = 1e-5   # below this a cell has no velocity


def _hll(hL, hR, uL, uR, vL, vR):
    """HLL flux of (h, h u_n, h u_t) for 1-D normal states; transverse momentum upwinded on mass."""
    cL = np.sqrt(G * hL)
    cR = np.sqrt(G * hR)
    SL = np.minimum(uL - cL, uR - cR)
    SR = np.maximum(uL + cL, uR + cR)
    qL, qR = hL * uL, hR * uR
    FL1, FR1 = qL, qR
    FL2, FR2 = qL * uL + 0.5 * G * hL * hL, qR * uR + 0.5 * G * hR * hR
    den = SR - SL
    with np.errstate(divide="ignore", invalid="ignore"):
        F1 = (SR * FL1 - SL * FR1 + SL * SR * (hR - hL)) / den
        F2 = (SR * FL2 - SL * FR2 + SL * SR * (qR - qL)) / den
    F1 = np.where(SL >= 0, FL1, np.where(SR <= 0, FR1, F1))
    F2 = np.where(SL >= 0, FL2, np.where(SR <= 0, FR2, F2))
    dry = den <= 0
    F1 = np.where(dry, 0.0, F1)
    F2 = np.where(dry, 0.0, F2)
    F3 = F1 * np.where(F1 > 0, vL, vR)
    return F1, F2, F3


class SweHLL(WaterSolver):
    """Full 2-D shallow-water equations, first-order finite volume (HLL + hydrostatic reconstruction).

    Approximates: the depth-averaged Saint-Venant equations with momentum advection, so it carries
    bores, hydraulic jumps, supercritical jets and dam-break fronts. Fluxes are HLL on each face after
    Audusse et al. (2004) hydrostatic reconstruction (face bed ``max(zL, zR)``, depths cut to it, plus
    the matching pressure correction), which makes it well-balanced: a lake at rest on any bed stays
    at rest, and depths stay non-negative under CFL <= 0.5. Unsplit (both directions from the same
    state), so it keeps mirror symmetry. Manning friction is semi-implicit and uses |u| (coupled), and
    walls are reflecting. Mass moves through the shared donor-limited transfer, so volume is exact;
    ``stats["limited"]`` counts faces the limiter had to cut (should stay ~0).

    Cannot do: it is first order, so it smears fronts and damps waves over a few dozen cells; the
    transverse momentum is upwinded (HLL has no contact wave), so shear layers are diffusive. Wet/dry
    fronts are robust but thin films carry noisy velocity (cut at ``H_DRY``). Explicit with the full
    gravity-wave CFL: ``dt <= cfl dx / max(|u| + sqrt(g h))``, the tightest of all water solvers here.
    """

    name = "swe_hll"
    family = "shallow water (full, finite volume)"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, cfl=0.45):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.cfl = cfl
        self.hu = np.zeros(grid.shape)
        self.hv = np.zeros(grid.shape)

    def _uv(self):
        wet = self.h > H_DRY
        hs = np.where(wet, self.h, 1.0)
        u = np.where(wet, self.hu / hs, 0.0)
        v = np.where(wet, self.hv / hs, 0.0)
        return u, v

    def stable_dt(self) -> float:
        u, v = self._uv()
        c = np.sqrt(G * self.h)
        s = max(float((np.abs(u) + c).max()), float((np.abs(v) + c).max()), math.sqrt(G * 1e-3))
        return self.cfl * self.grid.dx / s

    def velocity(self):
        return self._uv()

    def _faces(self, h, z, un, ut, wall, axis):
        """Face states along ``axis`` (1 = x, 0 = y) including the outer faces, with reflecting walls.
        Returns flux arrays with one more entry along ``axis`` than cells, as the left-cell and
        right-cell variants (Audusse pressure correction differs per side)."""
        def pad(a, neg=False):
            # ghost = mirror of the edge cell (negated normal velocity)
            if axis == 1:
                lo, hi = a[:, :1], a[:, -1:]
            else:
                lo, hi = a[:1, :], a[-1:, :]
            if neg:
                lo, hi = -lo, -hi
            return np.concatenate([lo, a, hi], axis=axis)
        hp, zp, unp, utp = pad(h), pad(z), pad(un, True), pad(ut)
        wp = np.concatenate([np.ones_like(wall[:, :1] if axis == 1 else wall[:1, :]), wall,
                             np.ones_like(wall[:, :1] if axis == 1 else wall[:1, :])], axis=axis)
        sl = (slice(None), slice(None, -1)) if axis == 1 else (slice(None, -1), slice(None))
        sr = (slice(None), slice(1, None)) if axis == 1 else (slice(1, None), slice(None))
        hL, hR = hp[sl].copy(), hp[sr].copy()
        zL, zR = zp[sl].copy(), zp[sr].copy()
        uL, uR = unp[sl].copy(), unp[sr].copy()
        tL, tR = utp[sl].copy(), utp[sr].copy()
        wL, wR = wp[sl], wp[sr]
        # the outer ghosts are already mirrors; an inner wall side becomes the mirror of its neighbour
        inner_wL = wL.copy()
        inner_wR = wR.copy()
        if axis == 1:
            inner_wL[:, 0] = False
            inner_wR[:, -1] = False
        else:
            inner_wL[0, :] = False
            inner_wR[-1, :] = False
        mL = inner_wL & ~wR
        mR = inner_wR & ~wL
        hL = np.where(mL, hR, hL); zL = np.where(mL, zR, zL); uL = np.where(mL, -uR, uL); tL = np.where(mL, tR, tL)
        hR = np.where(mR, hL, hR); zR = np.where(mR, zL, zR); uR = np.where(mR, -uL, uR); tR = np.where(mR, tL, tR)
        both = wL & wR
        zs = np.maximum(zL, zR)
        hLs = np.maximum(hL + zL - zs, 0.0)
        hRs = np.maximum(hR + zR - zs, 0.0)
        F1, F2, F3 = _hll(hLs, hRs, uL, uR, tL, tR)
        F1 = np.where(both, 0.0, F1); F2 = np.where(both, 0.0, F2); F3 = np.where(both, 0.0, F3)
        F2_left = F2 + 0.5 * G * (hL * hL - hLs * hLs)     # as seen by the left cell
        F2_right = F2 + 0.5 * G * (hR * hR - hRs * hRs)    # as seen by the right cell
        return F1, F2_left, F2_right, F3

    def _step(self, dt):
        dx = self.grid.dx
        h, z, wall = self.h, self.z, self.wall
        u, v = self._uv()
        # x faces (ny, nx+1): normal u, transverse v
        X1, X2L, X2R, X3 = self._faces(h, z, u, v, wall, axis=1)
        # y faces (ny+1, nx): normal v, transverse u
        Y1, Y2L, Y2R, Y3 = self._faces(h, z, v, u, wall, axis=0)
        r = dt / dx
        hu = self.hu - r * (X2L[:, 1:] - X2R[:, :-1]) - r * (Y3[1:, :] - Y3[:-1, :])
        hv = self.hv - r * (Y2L[1:, :] - Y2R[:-1, :]) - r * (X3[:, 1:] - X3[:, :-1])
        # mass through the shared transfer (interior faces), volumes this step
        self.transfer(X1[:, 1:-1] * dx * dt, Y1[1:-1, :] * dx * dt)
        Fx, Fy = self._last_F
        self.record_flux(Fx, Fy, dt)
        self.stats["limited"] = self.stats.get("limited", 0) + self.limited_faces
        hn = self.h
        wet = (hn > H_DRY) & ~wall
        hs = np.where(wet, hn, 1.0)
        if self.manning > 0:
            sp = np.hypot(hu, hv) / hs
            fac = 1.0 / (1.0 + dt * G * self.manning ** 2 * sp / hs ** (4 / 3))
            hu, hv = hu * fac, hv * fac
        self.hu = np.where(wet, hu, 0.0)
        self.hv = np.where(wet, hv, 0.0)

    def apply_forcing(self, dt):
        super().apply_forcing(dt)
        dry = self.h <= H_DRY
        self.hu[dry] = 0.0
        self.hv[dry] = 0.0

    def blown_up(self) -> bool:
        return super().blown_up() or not (np.isfinite(self.hu).all() and np.isfinite(self.hv).all())


SOLVERS = (SweHLL,)
