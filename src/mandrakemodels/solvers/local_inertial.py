"""Local-inertial shallow water (Bates, Horritt & Fewtrell 2010; de Almeida et al. 2012)."""
from __future__ import annotations

import numpy as np

from ..base import G, WaterSolver


class LocalInertial(WaterSolver):
    """Local-inertial shallow water on faces (LISFLOOD-FP "acceleration" kind).

    Approximates: the shallow-water momentum equation with local acceleration, hydrostatic pressure
    gradient and Manning bed friction, dropping the advective (u . grad u) term. Each face keeps a
    signed discharge with memory; the de Almeida (2012) theta weighting (``theta`` < 1 mixes in the
    neighbouring faces) damps the checkerboard mode; friction is semi-implicit and uses the coupled
    vector |q| (not axis-wise, avoiding the 2^(1/4) 45-degree error). Water moves through the shared
    donor-limited transfer, so volume is exact and depth never negative.

    Cannot do: fast jets, bores and hydraulic jumps correctly (no momentum advection); it is valid for
    Froude < ~0.5 (de Almeida & Bates 2013); wet/dry fronts recede crudely (Cozzolino 2019). Explicit:
    dt <= alpha dx / sqrt(g h_max). Default alpha 0.6: with theta 0.8 and alpha 0.7 a lake at rest on
    this rough bed grew roundoff into metre-scale sloshing (measured 2026-11-12 at 64, 128, 250); 0.6
    held it to 1e-15 at all three.
    """

    name = "local_inertial"
    family = "shallow water (local inertial)"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, theta=0.8, alpha=0.6):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.theta = theta
        self.alpha = alpha

    def stable_dt(self) -> float:
        return self.wave_dt(self.alpha)

    def _step(self, dt):
        dx, n2, th = self.grid.dx, self.manning ** 2, self.theta
        dEx, dEy = self.face_surface_diff()
        hx, hy = self.face_flow_depth()
        qx = self.qx[:, 1:-1]          # interior faces, previous step
        qy = self.qy[1:-1, :]
        # theta-weighted previous discharge (de Almeida 2012), neighbours along the flow direction
        qxp = np.pad(qx, ((0, 0), (1, 1)))
        qyp = np.pad(qy, ((1, 1), (0, 0)))
        qx_t = th * qx + 0.5 * (1 - th) * (qxp[:, :-2] + qxp[:, 2:])
        qy_t = th * qy + 0.5 * (1 - th) * (qyp[:-2, :] + qyp[2:, :])
        # coupled |q|: the cross-component averaged onto each face
        qyc = 0.5 * (self.qy[1:, :] + self.qy[:-1, :])        # cell-centred qy (ny, nx)
        qxc = 0.5 * (self.qx[:, 1:] + self.qx[:, :-1])        # cell-centred qx
        qy_on_x = 0.5 * (qyc[:, 1:] + qyc[:, :-1])
        qx_on_y = 0.5 * (qxc[1:, :] + qxc[:-1, :])
        magx = np.hypot(qx, qy_on_x)
        magy = np.hypot(qy, qx_on_y)
        with np.errstate(divide="ignore", invalid="ignore"):
            wx = hx > 1e-4
            wy = hy > 1e-4
            hx_s = np.where(wx, hx, 1.0)
            hy_s = np.where(wy, hy, 1.0)
            nqx = (qx_t - G * hx * dt * dEx / dx) / (1 + G * dt * n2 * magx / hx_s ** (7 / 3))
            nqy = (qy_t - G * hy * dt * dEy / dx) / (1 + G * dt * n2 * magy / hy_s ** (7 / 3))
        nqx = np.where(wx, nqx, 0.0)
        nqy = np.where(wy, nqy, 0.0)
        self.transfer(nqx * dx * dt, nqy * dx * dt)
        Fx, Fy = self._last_F
        self.record_flux(Fx, Fy, dt)


SOLVERS = (LocalInertial,)
