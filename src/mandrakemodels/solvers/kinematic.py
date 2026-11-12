"""Kinematic wave: flow down the bed slope."""
from __future__ import annotations

import numpy as np

from ..base import WaterSolver


class KinematicWave(WaterSolver):
    """Kinematic-wave overland flow on a grid (D8-free, four faces).

    Approximates: Manning flow driven by the BED slope only, q = h^(5/3) sqrt(S0) / n per unit width,
    where S0 is the magnitude of the downhill bed-slope vector and the discharge is shared among the
    downhill neighbours in proportion to their slope. Uses the donor cell's depth. Explicit; dt from
    the CFL condition on the kinematic celerity (5/3) v. Volume moves only through ``transfer``.

    Cannot do: backwater (a full pond does not hold back its inflow), flow up a bed slope, or keep a
    lake at rest on an uneven bed: the water just runs down to local bed minima and sits there as a
    set of puddles whose surfaces differ (the common test xfails it by name). Flat bed regions have no
    flow at all. Good only for steep, fast-draining overland flow and rain-runoff routing.
    """

    name = "kinematic_wave"
    family = "kinematic wave"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, alpha=0.5, dt_max=10.0):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.alpha, self.dt_max = alpha, dt_max

    def _slopes(self):
        """Downhill bed slopes from each cell to E, W, N, S neighbours (>= 0, open faces only)."""
        dx = self.grid.dx
        ny, nx = self.grid.shape
        z = self.z
        s = np.zeros((4, ny, nx))
        dzx = (z[:, :-1] - z[:, 1:]) / dx      # positive: bed falls toward +x
        dzy = (z[:-1, :] - z[1:, :]) / dx
        s[0, :, :-1] = np.where(self.open_x, np.maximum(dzx, 0), 0)
        s[1, :, 1:] = np.where(self.open_x, np.maximum(-dzx, 0), 0)
        s[2, :-1, :] = np.where(self.open_y, np.maximum(dzy, 0), 0)
        s[3, 1:, :] = np.where(self.open_y, np.maximum(-dzy, 0), 0)
        return s

    def _speed(self, h, s):
        S0 = np.sqrt((s ** 2).sum(0))
        return h ** (2 / 3) * np.sqrt(S0) / self.manning, S0

    def stable_dt(self) -> float:
        v, _ = self._speed(self.h, self._slopes())
        vmax = float(v.max())
        if vmax < 1e-12:
            return self.dt_max
        return min(self.dt_max, self.alpha * self.grid.dx / (5.0 / 3.0 * vmax))

    def _step(self, dt):
        dx = self.grid.dx
        s = self._slopes()
        v, S0 = self._speed(self.h, s)
        q = v * self.h                              # discharge per unit width
        # proportional to slope: weights s_i / sum s_i
        tot = s.sum(0)
        share = s / np.where(tot > 0, tot, 1.0)
        out = q[None] * share * dx * dt
        Fx = out[0, :, :-1] - out[1, :, 1:]
        Fy = out[2, :-1, :] - out[3, 1:, :]
        self.transfer(Fx, Fy)
        self.record_flux(*self._last_F, dt)


SOLVERS = (KinematicWave,)
