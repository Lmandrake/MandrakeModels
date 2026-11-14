"""local_inertial_frcap with a ripple detector and surface diffusion ONLY where it fires (2026-11-14).

Item RIMMANDRAKE_TRY_RIPPLE_DETECTOR_TREAT_ONLY_WHERE_IT_RIPPLES; the owner's idea: "I also wonder if the ripples are
detectable and thus could be fixed with a local treatment that only happens when needed for efficiency."
Results and the detector census: D:\\ClaudeWork\\rimmandrake\\conductor\\water_ripples_20261114\\RESULT_DETECTOR.md.
Kept self-contained (the global dampers live in a sibling's file); the original engines are untouched.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from ..base import G
from .local_inertial import LocalInertialFrCap

WET = 1e-3


def wet_smooth(e, wet):
    """1-2-1 filter in both axes over wet cells only; a flat wet surface stays exactly flat."""
    a = np.where(wet, e, 0.0)
    w = wet.astype(float)
    for ax in (0, 1):
        lo = [slice(None)] * 2
        hi = [slice(None)] * 2
        lo[ax], hi[ax] = slice(None, -1), slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        pa = np.zeros_like(a)
        pw = np.zeros_like(w)
        pa[hi] += a[lo]
        pa[lo] += a[hi]
        pw[hi] += w[lo]
        pw[lo] += w[hi]
        a = np.where(wet, 2 * a + pa, 0.0)
        w = np.where(wet, 2 * w + pw, 0.0)
    return np.where(wet, a / np.where(w > 0, w, 1.0), e)


def lap4(e, wet):
    """|eta - mean of its wet 4-neighbours| (0 on dry cells and on cells with no wet neighbour)."""
    s = np.zeros_like(e)
    c = np.zeros_like(e)
    for ax in (0, 1):
        lo = [slice(None)] * 2
        hi = [slice(None)] * 2
        lo[ax], hi[ax] = slice(None, -1), slice(1, None)
        lo, hi = tuple(lo), tuple(hi)
        s[hi] += np.where(wet[lo], e[lo], 0.0)
        c[hi] += wet[lo]
        s[lo] += np.where(wet[hi], e[hi], 0.0)
        c[lo] += wet[hi]
    return np.where(wet & (c > 0), np.abs(e - s / np.maximum(c, 1)), 0.0)


def ripple_flags(s, rough=1.5e-4, slow_fr=0.2):
    """The detector: wet cells whose 3x3-mean surface roughness |eta - mean(4 neighbours)| exceeds ``rough`` (m) while
    the local flow is slow (cell Froude below ``slow_fr``). Census of it and five other candidates: RESULT_DETECTOR.md."""
    wet = (s.h > WET) & ~s.wall
    lp = ndimage.uniform_filter(lap4(s.eta, wet), 3, mode="constant")
    u, v = s.velocity()
    fr = np.hypot(u, v) / np.sqrt(G * np.maximum(s.h, WET))
    return wet & (lp > rough) & (fr < slow_fr)


class LocalInertialFrCapRipSense(LocalInertialFrCap):
    """local_inertial_frcap plus a mass-only diffusion of the 1-2-1-smoothed surface where a ripple detector fired.

    The detector (``ripple_flags``) runs every ``every`` steps; a flag holds ``hold`` refreshes after it last fired. The
    treatment strength is a weight in [0, 1]: the flags box-blurred ``taper`` cells each side and eased in and out by
    ``ramp`` per refresh, because hard patch edges and on/off switching made new ripples (measured 2026-11-14: reservoir
    roughness at 32.7 s rose 3.0e-4 to 3.8e-4 with hard edges). On a treated face both cells are wet and the face carries
    ``-kappa * weight * sqrt(g h_f) * d(smoothed eta)/dx`` of water, moved by the shared donor-limited transfer (volume exact,
    depth never negative) and NOT remembered as discharge. A flat surface gets exactly nothing (lake at rest).

    Cannot do: the detector is not a clean ripple test (it also fires on swell and on slow river water; census in
    RESULT_DETECTOR.md), so in moving scenes it treats most slow water (dam break 50%, river 96%, swell 84% of wet
    cell-steps). What keeps the pit and breach right is its slow-flow gate, not the ripple test; a damper faded by Froude
    everywhere (the sibling's local_inertial_frcap_smoothdiff_fr4) passes the same bars at about the same cost.
    """

    name = "local_inertial_frcap_ripsense"
    family = "shallow water (local inertial)"
    kappa = 0.5
    rough = 1.5e-4
    slow_fr = 0.2          # not fr_max: that name is the parent's face Froude cap
    every = 10             # steps between detector runs
    hold = 2
    taper = 2              # cells: the treatment weight is a box blur of the flags this wide each side (0 = hard edges)
    ramp = 0.25            # per detector run: the weight moves this share of the way to its new value (1 = switch)
    apply_every = 4        # damp every n-th step at n times kappa: same bars as every step, cost 1.03-1.10x vs 1.04-1.20x

    def __init__(self, *a, **kw):
        for k in ("kappa", "rough", "slow_fr", "every", "hold", "taper", "ramp", "apply_every"):
            if k in kw:
                setattr(self, k, kw.pop(k))
        super().__init__(*a, **kw)
        self._age = np.full(self.h.shape, 1 << 30)
        self.flag = np.zeros(self.h.shape, bool)
        self.weight = np.zeros(self.h.shape)
        self._fd = None
        self._box = None
        self._n_treat = self._n_wet = 0
        self.stats.update(treated_cell_steps=0, wet_cell_steps=0, detector_runs=0)

    def _step(self, dt):
        if self.nsteps % self.every == 0:
            self._detect()
        self._dt = dt
        self._fd = None
        self.stats["treated_cell_steps"] += self._n_treat
        self.stats["wet_cell_steps"] += self._n_wet
        super()._step(dt)

    def _detect(self):
        f = ripple_flags(self, self.rough, self.slow_fr)
        self._age = np.where(f, 0, self._age + 1)
        self.flag = self._age <= self.hold
        target = self.flag.astype(float)
        if self.taper:
            target = ndimage.uniform_filter(target, 2 * self.taper + 1, mode="constant")
        self.weight += self.ramp * (target - self.weight)
        self.weight[self.weight < 1e-3] = 0.0
        self.stats["detector_runs"] += 1
        # the treatment works on the bounding box of weighted cells plus 2 (the smoothing reads 1 beyond a treated face)
        wet = (self.h > WET) & ~self.wall
        on = (self.weight > 0) & wet
        self._n_treat, self._n_wet = int(on.sum()), int(wet.sum())
        if self._n_treat == 0:
            self._box = None
            return
        rows, cols = np.flatnonzero(on.any(1)), np.flatnonzero(on.any(0))
        ny, nx = self.h.shape
        self._box = (slice(max(rows[0] - 2, 0), min(rows[-1] + 3, ny)), slice(max(cols[0] - 2, 0), min(cols[-1] + 3, nx)))

    def _limit(self, nqx, nqy, hx, hy):
        nqx, nqy = super()._limit(nqx, nqy, hx, hy)
        self._q_kept = (nqx, nqy)
        if self._box is None or self.nsteps % self.apply_every:
            return nqx, nqy
        r, c = self._box
        h = self.h[r, c]
        wet = (h > WET) & ~self.wall[r, c]
        wgt = np.where(wet, self.weight[r, c], 0.0)
        fl = wgt > 0
        e = wet_smooth(h + self.z[r, c], wet)
        # faces inside the box: x-faces between box columns, y-faces between box rows
        hxb, hyb = hx[r, c.start:c.stop - 1], hy[r.start:r.stop - 1, c]
        tx = (hxb > 1e-4) & wet[:, 1:] & wet[:, :-1] & (fl[:, 1:] | fl[:, :-1])
        ty = (hyb > 1e-4) & wet[1:, :] & wet[:-1, :] & (fl[1:, :] | fl[:-1, :])
        kx = (-0.5 * self.kappa * self.apply_every) * (wgt[:, 1:] + wgt[:, :-1]) * np.sqrt(G * hxb) * (e[:, 1:] - e[:, :-1])
        ky = (-0.5 * self.kappa * self.apply_every) * (wgt[1:, :] + wgt[:-1, :]) * np.sqrt(G * hyb) * (e[1:, :] - e[:-1, :])
        self._fd = (np.where(tx, kx, 0.0), np.where(ty, ky, 0.0))
        return nqx, nqy

    def _fd_slices(self):
        r, c = self._box
        return (r, slice(c.start, c.stop - 1)), (slice(r.start, r.stop - 1), c)

    def transfer(self, Fx, Fy):
        if self._fd is not None:
            w = self.grid.dx * self._dt
            sx, sy = self._fd_slices()
            Fx[sx] += self._fd[0] * w            # in place: the parent hands over a fresh product
            Fy[sy] += self._fd[1] * w
        super().transfer(Fx, Fy)

    def record_flux(self, Fx, Fy, dt):
        # remember the scheme's discharge (scaled by the donor limiter like the moved volume), not the diffusion
        if self._fd is None:
            return super().record_flux(Fx, Fy, dt)
        w = self.grid.dx * dt
        sx, sy = self._fd_slices()
        if self.limited_faces == 0:            # nothing scaled: the moved volume is exactly (q + diffusion) * w
            super().record_flux(Fx, Fy, dt)
            self.qx[:, 1:-1][sx] -= self._fd[0]
            self.qy[1:-1, :][sy] -= self._fd[1]
            return
        qx, qy = self._q_kept
        fdx, fdy = np.zeros_like(qx), np.zeros_like(qy)
        fdx[sx], fdy[sy] = self._fd
        tx, ty = (qx + fdx) * w, (qy + fdy) * w
        with np.errstate(divide="ignore", invalid="ignore"):
            kx = np.where(tx != 0, Fx / tx, 1.0)
            ky = np.where(ty != 0, Fy / ty, 1.0)
        self.qx[:, 1:-1] = qx * kx
        self.qy[1:-1, :] = qy * ky

    def treated_fraction(self) -> float:
        """Treated cell-steps over wet cell-steps since the start."""
        return self.stats["treated_cell_steps"] / max(self.stats["wet_cell_steps"], 1)


SOLVERS = (LocalInertialFrCapRipSense,)
