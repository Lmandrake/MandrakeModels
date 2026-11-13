"""Linear scalar wave equation (leapfrog / Verlet) on a cell grid.

Sources: second-order wave equation, leapfrog / Stormer-Verlet time stepping (Hairer, Lubich & Wanner
2003, Acta Numerica 12); CFL limit c dt/dx <= 1/sqrt(2) in 2-D. Own implementation.
"""
from __future__ import annotations

import math

import numpy as np

from ..base import Solver


class LinearWave(Solver):
    """Damped linear wave equation u_tt = c^2 lap(u) - damping * u_t, leapfrog (Stormer-Verlet).

    Approximates: small-amplitude surface or membrane waves (ripples on a pond) with a constant or
    per-cell speed c (c^2 = g h over a bed). Obstacle cells are clamped to 0 (a reflecting rigid
    boundary for the field); the outer edges are closed (zero normal gradient, Neumann), so waves
    reflect without phase flip. Drops from ``scene.init["drops"]`` (t, x_frac, y_frac) add a small
    Gaussian bump when the clock passes their time (deterministic).

    Cannot do: nonlinear steepening, breaking, dispersion (all frequencies travel at c, apart from the
    scheme's own grid dispersion at short wavelengths), wet/dry fronts, or true mass transport. The
    field is not conserved (damping removes energy; drops and obstacle clamping add or remove it).
    Explicit: dt <= dx / (c sqrt 2).
    """

    name = "linear_wave"
    family = "linear wave (leapfrog)"
    kind = "wave"
    conserves = False

    def __init__(self, grid, c=10.0, damping=0.0, obstacle=None, drops=(), drop_sigma=1.5,
                 drop_amp=1.0, u0=None, cfl=0.9):
        super().__init__(grid)
        self.c = np.broadcast_to(np.asarray(c, dtype=np.float64), grid.shape).copy()
        self.c2 = self.c ** 2
        self.damping = float(damping)
        self.obstacle = np.zeros(grid.shape, bool) if obstacle is None else np.array(obstacle, bool)
        self.drops = sorted([tuple(map(float, d)) for d in drops])
        self._next_drop = 0
        self.drop_sigma = drop_sigma
        self.drop_amp = drop_amp
        self.cfl = cfl
        self.u = np.zeros(grid.shape) if u0 is None else np.array(u0, dtype=np.float64)
        self.u[self.obstacle] = 0.0
        self.u_old = self.u.copy()
        self._dt_last = None
        self._fire_drops()

    @classmethod
    def from_scene(cls, scene, **kw):
        p = dict(scene.params)
        p.update(kw)
        p.pop("steps", None)
        return cls(scene.grid, obstacle=scene.init.get("obstacle"), drops=scene.init.get("drops", ()),
                   c=p.pop("c", 10.0), damping=p.pop("damping", 0.0), **p)

    def stable_dt(self) -> float:
        return self.cfl * self.grid.dx / (float(self.c.max()) * math.sqrt(2.0))

    def _fire_drops(self):
        ny, nx = self.grid.shape
        while self._next_drop < len(self.drops) and self.drops[self._next_drop][0] <= self.t + 1e-12:
            _, fx, fy = self.drops[self._next_drop]
            self._next_drop += 1
            X = np.arange(nx)[None, :] + 0.5
            Y = np.arange(ny)[:, None] + 0.5
            r2 = (X - fx * nx) ** 2 + (Y - fy * ny) ** 2
            bump = self.drop_amp * np.exp(-r2 / (2 * self.drop_sigma ** 2))
            bump[self.obstacle] = 0.0
            self.u += bump
            self.u_old += bump          # a displacement at rest (no initial velocity)

    @staticmethod
    def _lap(u):
        p = np.pad(u, 1, mode="edge")   # Neumann: ghost = edge value
        return p[1:-1, :-2] + p[1:-1, 2:] + p[:-2, 1:-1] + p[2:, 1:-1] - 4 * u

    def _step(self, dt):
        self._fire_drops()
        dx2 = self.grid.dx ** 2
        a = self.damping * dt / 2
        un = (2 * self.u - self.u_old * (1 - a) + dt * dt * self.c2 * self._lap(self.u) / dx2) / (1 + a)
        un[self.obstacle] = 0.0
        self.u_old, self.u = self.u, un
        self._dt_last = dt
        self.stats["energy"] = self.energy(dt)
        self.stats["max_abs"] = float(np.abs(self.u).max())

    def energy(self, dt=None) -> float:
        """Discrete leapfrog energy 0.5 sum[(u_t)^2 + c^2 <grad u_n, grad u_{n-1}>] * dx^2 (conserved
        exactly by the undamped scheme with free edges; obstacle cells are held at 0)."""
        dt = dt or self._dt_last
        if not dt:
            return 0.0
        ut = (self.u - self.u_old) / dt
        gx = (self.u[:, 1:] - self.u[:, :-1]) * (self.u_old[:, 1:] - self.u_old[:, :-1])
        gy = (self.u[1:, :] - self.u[:-1, :]) * (self.u_old[1:, :] - self.u_old[:-1, :])
        cx = np.minimum(self.c2[:, 1:], self.c2[:, :-1])
        cy = np.minimum(self.c2[1:, :], self.c2[:-1, :])
        pot = (cx * gx).sum() + (cy * gy).sum()
        return float(0.5 * ((ut ** 2).sum() * self.grid.dx ** 2 + pot))

    def blown_up(self) -> bool:
        return (not np.isfinite(self.u).all()) or float(np.abs(self.u).max()) > 1e6

    def snapshot(self) -> dict:
        return {"height": self.u.copy(), "obstacle": self.obstacle}


SOLVERS = (LinearWave,)
