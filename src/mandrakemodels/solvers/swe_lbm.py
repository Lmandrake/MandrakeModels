"""Shallow-water lattice Boltzmann (Zhou 2004), made volume-exact through the shared face transfer."""
from __future__ import annotations

import math

import numpy as np

from ..base import G, WaterSolver

H_DRY = 1e-5

# D2Q9 directions (ex, ey) in units of the lattice speed e; row index j is +y.
_D = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]])
_IDX = {(int(a), int(b)): n for n, (a, b) in enumerate(_D)}
_OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])
_AXIAL = np.array([False, True, True, True, True, False, False, False, False])
_WP = np.array([0.0, 1 / 6, 1 / 6, 1 / 6, 1 / 6, 1 / 24, 1 / 24, 1 / 24, 1 / 24])   # pressure weights
_WS = 2 * _WP                                                                      # bed-slope weights


def _add(arr, M, dj, di, s):
    """arr[j + dj, i + di] += s * M[j, i] for every (j, i) whose target lies inside arr."""
    ny, nx = M.shape
    j0, j1 = max(0, -dj), min(ny, arr.shape[0] - dj)
    i0, i1 = max(0, -di), min(nx, arr.shape[1] - di)
    if j1 > j0 and i1 > i0:
        arr[j0 + dj:j1 + dj, i0 + di:i1 + di] += s * M[j0:j1, i0:i1]


class SweLBM(WaterSolver):
    """Shallow-water lattice Boltzmann, D2Q9 (Zhou 2004), single-relaxation BGK.

    Approximates: the depth-averaged shallow-water equations. The nine populations carry depth
    (their sum is h) and momentum (h u); the equilibrium carries the g h^2 / 2 pressure; the bed slope
    enters as a centred force (h averaged over the two cells a population links, times the bed step),
    weighted so a lake at rest on any bed, in two dimensions, is reproduced exactly (the populations
    are then precisely the equilibrium); Manning friction is a semi-implicit force. The lattice speed
    is e = dx / dt, so each step streams a population one cell. ``stable_dt`` keeps
    e >= 2.5 (sqrt(g h) + |u|)_max, with hysteresis so e is not changed every step (a changed e
    resets the populations to equilibrium at the current h and u). Viscosity is
    nu = e dx (2 tau - 1) / 6.

    Volume is exact: the depth each population carries across a cell face is summed per face (a
    diagonal population is routed half x-then-y, half y-then-x) and moved only through
    ``transfer``; the source part of the flux is counted once, which makes this the same as the
    textbook LBM arrival sum. Walls, blocked shores (a dry cell whose bed is above the water surface
    of its neighbour) and the domain edge are half-way bounce-back. At the start of every step the
    populations are rescaled to the current depth (so forcing, digging and the donor limiter keep
    them consistent).

    Cannot do: it is a low-Froude scheme (|u| well below e, bores and jets near Froude 1 smear and
    ring, there is no shock capturing); dry fronts are crude (a thin film carries noisy populations,
    no hydrostatic reconstruction, so a flood front over a dry sloping bed can show negative
    populations that the donor limiter then cuts); there is numerical viscosity from tau; cost is nine
    distributions per cell and a diagonal-routing step, several times a local-inertial step.
    """

    name = "swe_lbm"
    family = "shallow water (lattice Boltzmann)"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, tau=0.8, e_factor=2.5):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.tau = float(tau)
        self.e_factor = float(e_factor)
        self._u = np.zeros(grid.shape)
        self._v = np.zeros(grid.shape)
        self.e = self._e_needed()
        self.f = self._feq(self.h, np.zeros(grid.shape), np.zeros(grid.shape), self.e)
        self.f[:, self.wall] = 0.0
        self._e_f = self.e
        self.stats["limited"] = 0

    # ------------------------------------------------------------------ lattice helpers
    @staticmethod
    def _feq(h, u, v, e):
        eu = e * (_D[:, 0, None, None] * u + _D[:, 1, None, None] * v)
        usq = u * u + v * v
        w2 = np.where(_AXIAL, 1 / 6, 1 / 24)[:, None, None]
        c = np.where(_AXIAL, 1 / 3, 1 / 12)[:, None, None]
        c2 = np.where(_AXIAL, 1 / 2, 1 / 8)[:, None, None]
        fe = (w2 * G * h * h / (e * e) + c * h * eu / (e * e) + c2 * h * eu * eu / e ** 4
              - w2 * h * usq / (e * e))
        fe[0] = h - 5 * G * h * h / (6 * e * e) - 2 * h * usq / (3 * e * e)
        return fe

    def _wave_speed_max(self):
        hh = np.where(self.wall, 0.0, self.h)
        c = np.sqrt(G * hh)
        wet = hh > H_DRY
        s = float((np.abs(self._u) + np.abs(self._v) + c)[wet].max()) if wet.any() else 0.0
        return max(s, math.sqrt(G * 1e-3))

    def _e_needed(self):
        return self.e_factor * self._wave_speed_max()

    def stable_dt(self) -> float:
        need = self._e_needed()
        if need > self.e or need < 0.5 * self.e:      # hysteresis: keep e (and the populations) steady
            self.e = need
        return self.grid.dx / self.e

    # ------------------------------------------------------------------ step
    def _moments(self, e):
        h = self.f.sum(0)
        hu = e * (_D[:, 0, None, None] * self.f).sum(0)
        hv = e * (_D[:, 1, None, None] * self.f).sum(0)
        wet = (self.h > H_DRY) & ~self.wall
        hs = np.where(wet, np.maximum(h, 1e-300), 1.0)
        u = np.where(wet, hu / hs, 0.0)
        v = np.where(wet, hv / hs, 0.0)
        sp = np.hypot(u, v)
        cap = 0.4 * e
        k = np.where(sp > cap, cap / np.maximum(sp, 1e-300), 1.0)
        return u * k, v * k

    def _step(self, dt):
        ny, nx = self.grid.shape
        dx = self.grid.dx
        e = dx / dt
        wall = self.wall
        h, z = self.h, self.z
        wet = (h > H_DRY) & ~wall
        # 1. make the populations consistent with h (forcing, events, limiter), and with e
        s = self.f.sum(0)
        ok = wet & (s > 1e-12)
        self.f = np.where(ok, self.f * (h / np.where(ok, s, 1.0)), 0.0)
        u, v = self._moments(self._e_f)
        if abs(e - self._e_f) > 1e-9 * e:
            self.f = self._feq(h, u, v, e)
            self.f[:, ~wet] = 0.0
            self._e_f = e
        else:
            self.f = np.where(wet, self.f, self._feq(h, 0 * u, 0 * v, e))
            self.f[:, wall] = 0.0
        u, v = self._moments(e)
        self._u, self._v = u, v
        # 2. collision (BGK) with friction force
        feq = self._feq(h, u, v, e)
        fb = self.f - (self.f - feq) / self.tau
        if self.manning > 0:
            sp = np.hypot(u, v)
            hs = np.where(wet, h, 1.0)
            fac = 1.0 / (1.0 + dt * G * self.manning ** 2 * sp / hs ** (4 / 3))
            dhu, dhv = (fac - 1.0) * h * u, (fac - 1.0) * h * v
            fb = fb + (e * (_D[:, 0, None, None] * dhu + _D[:, 1, None, None] * dhv)) / (6 * e * e)
        # 3. bed-slope source, blocked directions, per-direction moved depth
        eta = z + h
        hp = np.pad(h, 1)
        zp = np.pad(z, 1)
        wp = np.pad(wall, 1, constant_values=True)
        B = np.zeros_like(fb)             # full post-collision population (streams into the next f)
        Mflux = np.zeros_like(fb)         # depth that carries mass across faces (source counted once)
        allowed = np.zeros(fb.shape, bool)
        for k in range(1, 9):
            ex, ey = int(_D[k, 0]), int(_D[k, 1])
            sl = (slice(1 + ey, 1 + ey + ny), slice(1 + ex, 1 + ex + nx))
            hn, zn, wn = hp[sl], zp[sl], wp[sl]
            blocked = wall | wn | (h < H_DRY) | ((hn < H_DRY) & (zn > eta))
            if ex and ey:
                blocked = blocked | wp[1:1 + ny, 1 + ex:1 + ex + nx] | wp[1 + ey:1 + ey + ny, 1:1 + nx]
                blocked = blocked | (hp[1:1 + ny, 1 + ex:1 + ex + nx] < H_DRY) | (hp[1 + ey:1 + ey + ny, 1:1 + nx] < H_DRY)
            allowed[k] = ~blocked
            src = -_WS[k] * G / (e * e) * 0.5 * (h + hn) * (zn - z)
            src = np.where(blocked, 0.0, src)
            B[k] = fb[k] + src
            Mflux[k] = fb[k] + 0.5 * src
        B[0] = fb[0]
        # 4. face fluxes (volumes) -> the shared transfer
        FX = np.zeros((ny + 1, nx + 1))
        FY = np.zeros((ny + 1, nx + 1))
        for k in range(1, 9):
            ex, ey = int(_D[k, 0]), int(_D[k, 1])
            Mk = np.where(allowed[k], Mflux[k], 0.0)
            if not ex and not ey:
                continue
            if ey == 0:                         # axial x: face right of col i (ex=+1) or left (ex=-1)
                _add(FX, Mk, 0, 1 if ex > 0 else 0, ex)
            elif ex == 0:
                _add(FY, Mk, 1 if ey > 0 else 0, 0, ey)
            else:                               # diagonal: half x-then-y, half y-then-x
                half = 0.5 * Mk
                xi = 1 if ex > 0 else 0
                yj = 1 if ey > 0 else 0
                _add(FX, half, 0, xi, ex)                    # x leg at the source row
                _add(FX, half, ey, xi, ex)                   # x leg at the target row
                # y legs: at the target column (x first) and at the source column (y first)
                tcol = ex
                _add(FY, half, yj, tcol, ey)
                _add(FY, half, yj, 0, ey)
        area = dx * dx
        Fx = FX[:ny, 1:nx] * area
        Fy = FY[1:ny, :nx] * area
        self.transfer(Fx, Fy)
        self.record_flux(*self._last_F, dt)
        self.stats["limited"] += self.limited_faces
        # 5. stream the full populations (bounce back where blocked)
        fn = np.zeros_like(fb)
        fn[0] = B[0]
        for k in range(1, 9):
            ex, ey = int(_D[k, 0]), int(_D[k, 1])
            fn[k] += np.roll(np.where(allowed[k], B[k], 0.0), (ey, ex), axis=(0, 1))
            rest = np.where(allowed[k], 0.0, B[k])
            if ex and ey:                   # blocked diagonal: free-slip (specular) reflection at a flat wall
                kx = _IDX[(ex, 0)]          # axial neighbours of the diagonal
                ky = _IDX[(0, ey)]
                refl_y = ~allowed[ky] & allowed[kx] & allowed[_IDX[(-ex, 0)]]   # straight wall normal to y: keep ex, flip ey
                refl_x = ~allowed[kx] & allowed[ky] & allowed[_IDX[(0, -ey)]]   # straight wall normal to x: flip ex, keep ey
                fn[_IDX[(ex, -ey)]] += np.where(refl_y, rest, 0.0)
                fn[_IDX[(-ex, ey)]] += np.where(refl_x, rest, 0.0)
                fn[_OPP[k]] += np.where(refl_x | refl_y, 0.0, rest)
            else:
                fn[_OPP[k]] += rest
        fn[:, wall] = 0.0
        self.f = fn
        self.stats["max_u"] = float(np.hypot(u, v).max())
        self.stats["e"] = e

    def blown_up(self) -> bool:
        return super().blown_up() or not np.isfinite(self.f).all()

    def velocity(self):
        return self._u.copy(), self._v.copy()


SOLVERS = (SweLBM,)
