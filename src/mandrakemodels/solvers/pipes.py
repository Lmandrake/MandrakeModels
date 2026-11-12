"""Virtual pipes (O'Brien & Hodgins 1995; Mei, Decaudin & Hu 2007)."""
from __future__ import annotations

import math

import numpy as np

from ..base import G, WaterSolver


class VirtualPipes(WaterSolver):
    """Virtual-pipe water: every cell holds four outflow pipes with memory (the graphics workhorse).

    Approximates: each face is a pipe of cross-section ``A = dx^2`` and length ``dx``; the outflow in
    a pipe accelerates with the surface drop, ``f += dt * A * g * d_eta / dx`` (m^3/s), clamped at zero
    (pipes only push out), and all four outflows are scaled by ``K = min(1, V / (dt * sum f))`` so a
    cell never gives more than it holds (Mei 2007). The outflow memory gives momentum, so it sloshes,
    makes waves and can climb a little uphill. An optional semi-implicit Manning friction (using the
    donor depth) damps it; with ``manning=0`` it rings for a long time. Volume moves through the shared
    transfer, so it is exact.

    Cannot do: the pipe cross-section does not depend on depth, so a 1 cm film and a 3 m channel
    accelerate alike: shallow sheets run far too fast and wave speed is depth-independent
    (sqrt(g dx), not sqrt(g h)). It is not a discretisation of the shallow-water equations; timing and
    bores are not physical. Explicit: dt <= ~0.5 sqrt(dx / g) whatever the depth.
    """

    name = "pipes"
    family = "virtual pipes (graphics)"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, alpha=0.4):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.alpha = alpha
        ny, nx = grid.shape
        # outflows (m^3/s) on interior faces: +x out of the left cell, -x out of the right cell, etc.
        self.f_r = np.zeros((ny, nx - 1))
        self.f_l = np.zeros((ny, nx - 1))
        self.f_u = np.zeros((ny - 1, nx))
        self.f_d = np.zeros((ny - 1, nx))

    def stable_dt(self) -> float:
        return self.alpha * math.sqrt(self.grid.dx / G)

    def _step(self, dt):
        dx, A = self.grid.dx, self.grid.area
        dEx, dEy = self.face_surface_diff()           # eta[right] - eta[left]
        hx, hy = self.face_flow_depth()
        k = dt * dx * G                               # dt * A * g / l with A = dx^2, l = dx
        fr = np.maximum(self.f_r - k * dEx, 0.0)      # left cell -> right, driven by left being higher
        fl = np.maximum(self.f_l + k * dEx, 0.0)
        fu = np.maximum(self.f_u - k * dEy, 0.0)
        fd = np.maximum(self.f_d + k * dEy, 0.0)
        if self.manning > 0:
            n2 = self.manning ** 2
            with np.errstate(divide="ignore", invalid="ignore"):
                def fric(f, hf):
                    hs = np.where(hf > 1e-4, hf, 1.0)
                    u = f / (dx * hs)
                    return np.where(hf > 1e-4, f / (1 + dt * G * n2 * u / hs ** (4 / 3)), 0.0)
                fr, fl, fu, fd = fric(fr, hx), fric(fl, hx), fric(fu, hy), fric(fd, hy)
        fr = np.where(self.open_x, fr, 0.0)
        fl = np.where(self.open_x, fl, 0.0)
        fu = np.where(self.open_y, fu, 0.0)
        fd = np.where(self.open_y, fd, 0.0)
        # Mei's K: scale every outflow of a cell so it never sends more than it holds
        out = np.zeros(self.grid.shape)
        out[:, :-1] += fr
        out[:, 1:] += fl
        out[:-1, :] += fu
        out[1:, :] += fd
        V = self.h * A
        with np.errstate(divide="ignore", invalid="ignore"):
            K = np.where(out * dt > V, V / (out * dt), 1.0)
        K = np.where(out > 0, K, 1.0)
        fr = fr * K[:, :-1]
        fl = fl * K[:, 1:]
        fu = fu * K[:-1, :]
        fd = fd * K[1:, :]
        self.f_r, self.f_l, self.f_u, self.f_d = fr, fl, fu, fd
        self.transfer((fr - fl) * dt, (fu - fd) * dt)
        Fx, Fy = self._last_F
        self.record_flux(Fx, Fy, dt)


SOLVERS = (VirtualPipes,)
