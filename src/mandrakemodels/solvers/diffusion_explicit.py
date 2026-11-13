"""Explicit diffusion-wave (zero-inertia) shallow water.

Sources: zero-inertia / diffusion wave of Hromadka & Yen (1986); the explicit 2-D raster form and its
quadratic time-step limit, and the slope floor, follow Hunter, Horritt, Bates, Wilson & Werner (2005),
Adv. Water Resour. 28. Own implementation.
"""
from __future__ import annotations

import numpy as np

from ..base import WaterSolver


def face_conductance(solver, eta, s_min):
    """Diffusion-wave conductance K = h_f^(5/3) / (n sqrt(max(|S|, s_min))) (m^2/s) on interior x and y
    faces for the surface ``eta``, with the coupled-vector slope magnitude |S| = hypot(Sx, Sy) (the
    cross component taken from the cell-centred gradient averaged onto the face). Zero on closed or
    dry faces. Flux per unit width is q = -K dEta/dx; volume per second is K (eta_i - eta_j)."""
    z, dx, n = solver.z, solver.grid.dx, solver.manning
    hx = np.maximum(np.maximum(eta[:, 1:], eta[:, :-1]) - np.maximum(z[:, 1:], z[:, :-1]), 0.0)
    hy = np.maximum(np.maximum(eta[1:, :], eta[:-1, :]) - np.maximum(z[1:, :], z[:-1, :]), 0.0)
    hx = hx * solver.open_x
    hy = hy * solver.open_y
    dEx = eta[:, 1:] - eta[:, :-1]
    dEy = eta[1:, :] - eta[:-1, :]
    # cell-centred gradients (one-sided at the edges), then averaged onto the faces
    gx = np.zeros_like(eta)
    gx[:, 1:-1] = 0.5 * (dEx[:, 1:] + dEx[:, :-1])
    gx[:, 0], gx[:, -1] = dEx[:, 0], dEx[:, -1]
    gy = np.zeros_like(eta)
    gy[1:-1, :] = 0.5 * (dEy[1:, :] + dEy[:-1, :])
    gy[0, :], gy[-1, :] = dEy[0, :], dEy[-1, :]
    cx = 0.5 * (gy[:, 1:] + gy[:, :-1])       # y-gradient on x-faces
    cy = 0.5 * (gx[1:, :] + gx[:-1, :])       # x-gradient on y-faces
    Sx = np.hypot(dEx / dx, cx / dx)
    Sy = np.hypot(dEy / dx, cy / dx)
    Kx = hx ** (5 / 3) / (n * np.sqrt(np.maximum(Sx, s_min)))
    Ky = hy ** (5 / 3) / (n * np.sqrt(np.maximum(Sy, s_min)))
    return Kx, Ky


class DiffusionWaveExplicit(WaterSolver):
    """Explicit zero-inertia (diffusion) wave: q = -(h_f^(5/3)/n) S / sqrt(max(|S|, S_min)).

    Approximates: the shallow-water momentum balance with only the pressure gradient and Manning
    friction (inertia dropped); S is the SURFACE slope on the face, h_f the Bates face flow depth, and
    |S| the coupled-vector slope magnitude (cross component from neighbouring faces). Backwater and a
    lake at rest are handled (flow follows the surface). Volume moves only through ``transfer``.

    Cannot do: waves, inertia, bores, supercritical flow, or fast flow that outruns the surface
    gradient. THE FAMOUS WEAKNESS: it is a nonlinear diffusion with D = h_f^(5/3)/(n sqrt(max(|S|,
    S_min))), which grows without bound as the slope flattens, so the explicit limit dt <= dx^2/(4
    D_max) collapses toward 0 as S_min -> 0 (S_min = 1e-4 gives D ~ 100 x deeper-than-1m values and
    sub-millisecond steps on flat water). S_min is therefore a stability knob, not just a regulariser.
    """

    name = "diffusion_wave_explicit"
    family = "diffusion wave"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, s_min=1e-4,
                 safety=0.9, dt_max=10.0):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.s_min, self.safety, self.dt_max = s_min, safety, dt_max

    def stable_dt(self) -> float:
        Kx, Ky = face_conductance(self, self.eta, self.s_min)
        D = max(float(Kx.max(initial=0.0)), float(Ky.max(initial=0.0)))
        if D <= 0:
            return self.dt_max
        return min(self.dt_max, self.safety * self.grid.dx ** 2 / (4.0 * D))

    def _step(self, dt):
        Kx, Ky = face_conductance(self, self.eta, self.s_min)
        dEx, dEy = self.face_surface_diff()
        self.transfer(-Kx * dEx * dt, -Ky * dEy * dt)
        self.record_flux(*self._last_F, dt)


SOLVERS = (DiffusionWaveExplicit,)
