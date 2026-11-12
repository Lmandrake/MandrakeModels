"""Weighted cellular automaton for flood inundation (WCA2D)."""
from __future__ import annotations

import numpy as np

from ..base import G, WaterSolver


class WCA2D(WaterSolver):
    """Weighted cellular automaton in the style of WCA2D (Guidolin et al. 2016, Env. Modelling &
    Software 84).

    RECONSTRUCTED FROM THE PAPER'S DESCRIPTION, NOT CHECKED AGAINST THE CADDIES SOURCE: details such as
    where the tolerance enters, what "previous outflow" means in the cap, and the dt rule may differ
    from the reference code.

    Approximates: diffusion-like flood spreading. Per cell, level differences dz_i to the (four,
    von Neumann) lower neighbours give volumes dV_i = A dz_i; weights w_i = dV_i/(dV_min + dV_tot),
    the cell keeping w_0 = dV_min/(dV_min + dV_tot). Total outflow is limited by the cell's water, by
    I_M / w_M (I_M = v d dx dt through the steepest face with v = min(sqrt(g d), Manning velocity)) and
    by dV_min + the cell's previous outflow, and (an addition, for stability: without it the explicit
    update amplifies roundoff into a mirror-asymmetric mess) by ``max_frac`` (0.2) of dV_tot; it is shared in proportion to the weights. Differences
    below ``tau`` are ignored. dt adapts: dt = alpha dx / max(v), clamped to [dt_min, dt_max]. The
    update is synchronous and goes through ``transfer``, so volume is exact and depth >= 0.

    Cannot do: inertia, waves, bores or momentum (the slope-driven velocity is a limiter, not a state);
    the tolerance leaves a residual surface slope of up to ``tau`` per cell; timing is only as good as
    the velocity limiter.
    """

    name = "wca2d"
    family = "cellular automaton (weighted)"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, tau=1e-3, alpha=0.5,
                 dt_min=1e-3, dt_max=1.0, max_frac=0.2):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.tau, self.alpha, self.dt_min, self.dt_max = tau, alpha, dt_min, dt_max
        self.max_frac = max_frac
        self._prev_out = np.zeros(grid.shape)

    def _face_vel(self):
        """Velocity limiter and depth on interior faces, from the current state."""
        dEx, dEy = self.face_surface_diff()
        hx, hy = self.face_flow_depth()
        dx, n = self.grid.dx, self.manning
        out = []
        for dE, hf in ((dEx, hx), (dEy, hy)):
            S = np.abs(dE) / dx
            v = np.minimum(np.sqrt(G * hf), hf ** (2 / 3) * np.sqrt(S) / n)
            v = np.where(np.abs(dE) > self.tau, v, 0.0)
            out.append((dE, hf, v))
        return out

    def stable_dt(self) -> float:
        (_, _, vx), (_, _, vy) = self._face_vel()
        vmax = max(float(vx.max(initial=0.0)), float(vy.max(initial=0.0)))
        if vmax <= 1e-12:
            return self.dt_max
        return float(np.clip(self.alpha * self.grid.dx / vmax, self.dt_min, self.dt_max))

    def _step(self, dt):
        dx, A = self.grid.dx, self.grid.area
        (dEx, hx, vx), (dEy, hy, vy) = self._face_vel()
        ny, nx = self.grid.shape
        # per cell: drop toward E, W, N(+y), S(-y) neighbour; 0 where not lower, dry, closed or < tau
        def drop(d, hf, sign):
            dz = np.where(hf > 1e-9, sign * d, 0.0)
            return np.where(dz > self.tau, dz, 0.0)
        dz = np.zeros((4, ny, nx))
        dz[0, :, :-1] = drop(dEx, hx, -1.0)      # to +x: eta_i - eta_{i+1} = -dEx
        dz[1, :, 1:] = drop(dEx, hx, 1.0)        # to -x
        dz[2, :-1, :] = drop(dEy, hy, -1.0)
        dz[3, 1:, :] = drop(dEy, hy, 1.0)
        vel = np.zeros((4, ny, nx))
        vel[0, :, :-1], vel[1, :, 1:] = vx, vx
        vel[2, :-1, :], vel[3, 1:, :] = vy, vy
        dep = np.zeros((4, ny, nx))
        dep[0, :, :-1], dep[1, :, 1:] = hx, hx
        dep[2, :-1, :], dep[3, 1:, :] = hy, hy
        dV = dz * A
        dVtot = (dV[0] + dV[1]) + (dV[2] + dV[3])      # pairwise: bitwise mirror-symmetric
        lower = dV > 0
        dVmin = np.where(lower, dV, np.inf).min(0)
        has = lower.any(0)
        dVmin = np.where(has, dVmin, 0.0)
        den = np.where(has, dVmin + dVtot, 1.0)
        w = dV / den
        # velocity-volume limit through every lower face, I_i / w_i (the steepest face, M, normally binds)
        with np.errstate(divide="ignore", invalid="ignore"):
            Ilim = np.where(lower, vel * dep * dx * dt / np.where(w > 0, w, 1.0), np.inf).min(0)
        Ilim = np.where(has, Ilim, 0.0)
        T = np.minimum(np.minimum(self.h * A, Ilim), dVmin + self._prev_out)
        T = np.minimum(T, self.max_frac * dVtot)     # never give more than a stable share of the head
        T = np.where(has, np.maximum(T, 0.0), 0.0)
        share = dV / np.where(dVtot > 0, dVtot, 1.0)
        out = T[None] * share
        Fx = out[0, :, :-1] - out[1, :, 1:]
        Fy = out[2, :-1, :] - out[3, 1:, :]
        self.transfer(Fx, Fy)
        Fx, Fy = self._last_F
        self._prev_out = T
        self.record_flux(Fx, Fy, dt)


SOLVERS = (WCA2D,)
