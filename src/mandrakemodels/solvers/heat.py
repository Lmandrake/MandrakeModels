"""Heat / diffusion of a scalar on a cell grid: explicit FTCS and Peaceman-Rachford ADI."""
from __future__ import annotations

import math

import numpy as np

from ..base import Solver


class _HeatBase(Solver):
    kind = "scalar"
    conserves = True

    def __init__(self, grid, field, D=1.0):
        super().__init__(grid)
        self.u = np.array(field, dtype=np.float64)
        self.D = float(D)
        self.total0 = self.total()

    @classmethod
    def from_scene(cls, scene, **kw):
        p = dict(scene.params)
        p.update(kw)
        return cls(scene.grid, scene.init["field"], D=p.get("D", 1.0))

    def total(self) -> float:
        return float(self.u.sum() * self.grid.area)

    def blown_up(self) -> bool:
        return (not np.isfinite(self.u).all()) or float(np.abs(self.u).max()) > 1e6

    def _book(self):
        tot = self.total()
        self.stats["total"] = tot
        self.stats["total_rel_error"] = (tot - self.total0) / (abs(self.total0) or 1.0)
        self.stats["max"] = float(self.u.max())
        self.stats["min"] = float(self.u.min())

    def diagnostics(self) -> dict:
        self._book()
        return super().diagnostics()

    def snapshot(self) -> dict:
        return {"field": self.u.copy()}


def _lap1(u, axis):
    """1-D Neumann (insulated) second difference along ``axis`` in finite-volume form: sum of face
    differences, zero flux through the outer faces."""
    d = np.diff(u, axis=axis)
    out = np.zeros_like(u)
    if axis == 1:
        out[:, :-1] += d
        out[:, 1:] -= d
    else:
        out[:-1, :] += d
        out[1:, :] -= d
    return out


class HeatExplicit(_HeatBase):
    """Forward-time centred-space (FTCS) heat equation u_t = D lap(u), finite-volume form.

    Approximates: linear diffusion with constant D and insulated walls (zero flux through every outer
    face). The update moves quantity through shared faces with opposite signs, so the total is
    conserved to round-off.

    Cannot do: large steps (conditionally stable, dt <= dx^2 / 4D), variable D, sources, or sharp-front
    sensitivity (it is first order in time). Past the bound it blows up (checkerboard growth).
    """

    name = "heat_explicit"
    family = "heat (explicit FTCS)"

    def stable_dt(self) -> float:
        return 0.9 * self.grid.dx ** 2 / (4 * self.D)

    def _step(self, dt):
        r = self.D * dt / self.grid.dx ** 2
        self.u += r * (_lap1(self.u, 0) + _lap1(self.u, 1))
        self._book()


class _Thomas:
    """Batched tridiagonal solver for (I - s L) x = d with the 1-D Neumann Laplacian L. Coefficients
    are constant along a line, so the elimination is factored once; the sweep loops over the n
    positions along the line and is vectorised across every other line (no loop over rows)."""

    def __init__(self, n, s):
        self.n = n
        diag = np.full(n, 1 + 2 * s)
        diag[0] -= s
        diag[-1] -= s
        off = -s
        cp = np.zeros(n)          # modified super-diagonal
        den = np.zeros(n)         # pivots
        den[0] = diag[0]
        cp[0] = off / den[0]
        for i in range(1, n):
            den[i] = diag[i] - off * cp[i - 1]
            cp[i] = off / den[i]
        self.cp, self.den, self.off = cp, den, off

    def solve(self, d, axis):
        """Solve along ``axis`` of the 2-D array ``d``."""
        x = np.moveaxis(d, axis, 0).copy()      # (n, batch)
        n = self.n
        x[0] /= self.den[0]
        for i in range(1, n):
            x[i] = (x[i] - self.off * x[i - 1]) / self.den[i]
        for i in range(n - 2, -1, -1):
            x[i] -= self.cp[i] * x[i + 1]
        return np.moveaxis(x, 0, axis)


class HeatADI(_HeatBase):
    """Peaceman-Rachford alternating-direction implicit heat equation, finite-volume form.

    Approximates: the same linear diffusion u_t = D lap(u) with insulated walls, second order in time,
    unconditionally stable. Each step is two half steps (implicit in x with explicit y, then the
    reverse), solved with a vectorised batched Thomas algorithm. Every matrix has unit column sums, so
    the total is conserved to round-off.

    Cannot do: positivity at large steps (Peaceman-Rachford rings on sharp fronts when D dt / dx^2 is
    large; the solution is stable and conservative but can undershoot below zero), variable D, or
    sources. ``dt_hint`` is 10x the explicit limit, a sensible step for animation.
    """

    name = "heat_adi"
    family = "heat (ADI Peaceman-Rachford)"

    def __init__(self, grid, field, D=1.0):
        super().__init__(grid, field, D)
        self._cache = {}
        self.dt_hint = 10 * 0.9 * grid.dx ** 2 / (4 * self.D)

    def stable_dt(self) -> float:
        return math.inf

    def _solvers(self, r):
        key = round(r, 14)
        if key not in self._cache:
            self._cache = {key: (_Thomas(self.grid.nx, r / 2), _Thomas(self.grid.ny, r / 2))}
        return self._cache[key]

    def _step(self, dt):
        r = self.D * dt / self.grid.dx ** 2
        sx, sy = self._solvers(r)
        u = self.u
        # (I - r/2 Lx) u* = (I + r/2 Ly) u
        u = sx.solve(u + 0.5 * r * _lap1(u, 0), axis=1)
        # (I - r/2 Ly) u+ = (I + r/2 Lx) u*
        u = sy.solve(u + 0.5 * r * _lap1(u, 1), axis=0)
        self.u = u
        self._book()


SOLVERS = (HeatExplicit, HeatADI)
