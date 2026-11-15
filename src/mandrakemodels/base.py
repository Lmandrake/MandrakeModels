"""Common solver interface on 2-D gridded cells.

Conventions
-----------
* Arrays are ``(ny, nx)``, row ``j`` is y, column ``i`` is x. ``dx`` is the cell width in metres.
* Water solvers carry bed ``z`` (m), true depth ``h`` (m, >= 0) and a ``wall`` mask (no water, no flux).
  Surface ``eta = z + h``.
* Face fluxes: ``qx`` is ``(ny, nx+1)``, ``qy`` is ``(ny+1, nx)``, signed volume flux per unit width
  (m^2/s), positive toward +x / +y. The outer faces are closed (zero); every external volume enters or
  leaves through a booked ``Forcing`` term, so ``volume_error()`` is exact bookkeeping.
* A core step must move water only through shared faces (one transfer, opposite signs). The
  ``transfer`` helper does that with a donor limiter, so depth never goes negative and nothing leaks.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

G = 9.81
EPS_H = 1e-6  # depth below which a cell counts as dry (m)


def even_dt(dt: float, left: float) -> float:
    """The stepping rule for water solvers: split the time ``left`` to the next target (checkpoint, frame, event) into
    equal steps no longer than ``dt``, never full steps plus a short clipped tail. Stated once, with its evidence, at
    ``EvenDt`` in ``csharp/LiTiming/LocalInertialDL.cs``."""
    return left / max(1.0, math.ceil(left / dt - 1e-9))


@dataclass
class Grid:
    ny: int
    nx: int
    dx: float = 1.0

    @property
    def shape(self):
        return (self.ny, self.nx)

    @property
    def area(self) -> float:
        return self.dx * self.dx

    def coords(self):
        """Cell-centre coordinates (X, Y) in metres."""
        x = (np.arange(self.nx) + 0.5) * self.dx
        y = (np.arange(self.ny) + 0.5) * self.dx
        return np.meshgrid(x, y)


class Solver:
    """Base for every model. Subclasses implement ``_step(dt)``; the docstring of each subclass says what
    physics it approximates and what it cannot do."""

    name = "base"
    family = ""
    kind = "abstract"      # water | wave | scalar | flow
    conserves = False      # does the core step conserve the primary quantity exactly?

    def __init__(self, grid: Grid):
        self.grid = grid
        self.t = 0.0
        self.nsteps = 0
        self.stats: dict = {}

    def stable_dt(self) -> float:
        """Largest dt the scheme's stability bound allows for the current state (inf if unconditional)."""
        return math.inf

    def step(self, dt: float) -> None:
        self._step(dt)
        self.t += dt
        self.nsteps += 1

    def _step(self, dt: float) -> None:  # pragma: no cover
        raise NotImplementedError

    def blown_up(self) -> bool:
        return False

    def diagnostics(self) -> dict:
        return {"t": self.t, "steps": self.nsteps, **self.stats}

    def snapshot(self) -> dict:
        return {}


# --------------------------------------------------------------------------------------------- forcing

@dataclass
class Ledger:
    """External volume (m^3) booked by kind. ``net`` is in minus out."""
    terms: dict = field(default_factory=dict)
    gross: float = 0.0          # sum of |every booking|, so an in/out exchange does not cancel

    def book(self, kind: str, vol: float) -> None:
        self.terms[kind] = self.terms.get(kind, 0.0) + float(vol)
        self.gross += abs(float(vol))

    @property
    def net(self) -> float:
        return sum(self.terms.values())

    @property
    def throughput(self) -> float:
        return self.gross


@dataclass
class Source:
    """Inflow of ``Q`` m^3/s spread evenly over ``mask`` (``Q`` may be a function of t)."""
    mask: np.ndarray
    Q: float | Callable[[float], float]


@dataclass
class Sink:
    """Drains up to ``rate`` m/s of depth from ``mask`` (``math.inf`` empties it every step)."""
    mask: np.ndarray
    rate: float = math.inf


@dataclass
class FixedLevel:
    """Holds the surface at ``level`` (m, or function of t) on ``mask``: a level-setter, ocean edge, or
    wave maker. Adds or removes whatever it takes and books it."""
    mask: np.ndarray
    level: float | Callable[[float], float]


@dataclass
class Groundwater:
    """Seepage up from a water table ``table`` (m, array or scalar): rate ``k * (table - eta)`` m/s
    where the surface is below the table, capped so it never overshoots the table."""
    table: float | np.ndarray
    k: float = 1e-4
    mask: np.ndarray | None = None


@dataclass
class Forcing:
    sources: list = field(default_factory=list)
    sinks: list = field(default_factory=list)
    levels: list = field(default_factory=list)
    rain: float | Callable[[float], float] = 0.0          # m/s on every open cell
    evaporation: float = 0.0                             # m/s on wet cells
    groundwater: Groundwater | None = None


def _val(v, t):
    return v(t) if callable(v) else v


class WaterSolver(Solver):
    """Base for depth-on-bed water models. Holds bed, depth, walls, face fluxes and the forcing ledger."""

    kind = "water"
    conserves = True

    def __init__(self, grid: Grid, z: np.ndarray, h: np.ndarray | None = None,
                 wall: np.ndarray | None = None, forcing: Forcing | None = None,
                 manning: float = 0.03):
        super().__init__(grid)
        self.z = np.array(z, dtype=np.float64)
        self.h = np.zeros(grid.shape) if h is None else np.array(h, dtype=np.float64)
        self.wall = np.zeros(grid.shape, bool) if wall is None else np.array(wall, bool)
        self.h[self.wall] = 0.0
        self.forcing = forcing or Forcing()
        self.manning = manning
        self.qx = np.zeros((grid.ny, grid.nx + 1))
        self.qy = np.zeros((grid.ny + 1, grid.nx))
        self.ledger = Ledger()
        self._set_open()
        self.V0 = self.volume()
        self.min_h_seen = float(self.h.min())
        self.limited_faces = 0

    @classmethod
    def from_scene(cls, scene, **kw):
        import copy
        p = dict(scene.params)
        p.update(kw)
        return cls(scene.grid, z=scene.init["z"], h=scene.init["h"], wall=scene.init["wall"],
                   forcing=copy.deepcopy(scene.forcing), **p)

    # -- geometry
    def _set_open(self):
        o = ~self.wall
        self.open_x = o[:, 1:] & o[:, :-1]      # interior x-faces (ny, nx-1)
        self.open_y = o[1:, :] & o[:-1, :]      # interior y-faces (ny-1, nx)

    @property
    def eta(self) -> np.ndarray:
        return self.z + self.h

    def volume(self) -> float:
        return float(self.h[~self.wall].sum() * self.grid.area)

    def volume_error(self) -> float:
        return self.volume() - self.V0 - self.ledger.net

    def blown_up(self) -> bool:
        return (not np.isfinite(self.h).all()) or self.h.min() < -1e-6 or self.h.max() > 1e4

    # -- the one conservative transfer
    def transfer(self, Fx: np.ndarray, Fy: np.ndarray) -> None:
        """Move volume through interior faces. ``Fx (ny, nx-1)``: m^3 this step, positive toward +x;
        ``Fy (ny-1, nx)`` toward +y. Each donor's total outflow is scaled down to what it holds (donor
        limiter), then every face is applied once with opposite signs. Exact, and depth stays >= 0."""
        Fx = np.where(self.open_x, Fx, 0.0)
        Fy = np.where(self.open_y, Fy, 0.0)
        V = self.h * self.grid.area
        out = np.zeros_like(V)
        out[:, :-1] += np.maximum(Fx, 0)
        out[:, 1:] += np.maximum(-Fx, 0)
        out[:-1, :] += np.maximum(Fy, 0)
        out[1:, :] += np.maximum(-Fy, 0)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            s = np.where(out > V, V / out, 1.0)
        s = np.where(out > 0, s, 1.0)
        sx = np.where(Fx > 0, s[:, :-1], s[:, 1:])
        sy = np.where(Fy > 0, s[:-1, :], s[1:, :])
        self.limited_faces = int((sx < 1).sum() + (sy < 1).sum())
        Fx = Fx * sx
        Fy = Fy * sy
        dV = np.zeros_like(V)
        dV[:, :-1] -= Fx
        dV[:, 1:] += Fx
        dV[:-1, :] -= Fy
        dV[1:, :] += Fy
        self.h = np.maximum((V + dV) / self.grid.area, 0.0)
        self._last_F = (Fx, Fy)

    def record_flux(self, Fx, Fy, dt):
        """Store the per-width flux of this step (for velocity and current) from transferred volumes."""
        w = self.grid.dx * dt
        self.qx[:, 1:-1] = Fx / w
        self.qy[1:-1, :] = Fy / w

    # -- velocity from face fluxes
    def velocity(self):
        hc = np.maximum(self.h, 1e-3)
        u = 0.5 * (self.qx[:, 1:] + self.qx[:, :-1]) / hc
        v = 0.5 * (self.qy[1:, :] + self.qy[:-1, :]) / hc
        dry = self.h < 1e-3
        u[dry] = 0
        v[dry] = 0
        return u, v

    # -- forcing, after each core step
    def apply_forcing(self, dt: float) -> None:
        f, A, t = self.forcing, self.grid.area, self.t + dt
        o = ~self.wall
        for s in f.sources:
            m = s.mask & o
            n = int(m.sum())
            if n:
                q = _val(s.Q, t) * dt
                self.h[m] += q / (A * n)
                self.ledger.book("source", q)
        r = _val(f.rain, self.t + 0.5 * dt)      # midpoint: an end-time sample gave the step ending at a cutoff none
        if r:
            self.h[o] += r * dt
            self.ledger.book("rain", r * dt * A * o.sum())
        gw = f.groundwater
        if gw is not None:
            m = o if gw.mask is None else (gw.mask & o)
            deficit = np.where(m, np.maximum(gw.table - self.eta, 0.0), 0.0)
            add = np.minimum(gw.k * deficit * dt, deficit)
            self.h += add
            self.ledger.book("groundwater", add.sum() * A)
        if f.evaporation:
            take = np.minimum(self.h, f.evaporation * dt)
            take[self.wall] = 0
            self.h -= take
            self.ledger.book("evaporation", -take.sum() * A)
        for s in f.sinks:
            take = np.where(s.mask & o, np.minimum(self.h, s.rate * dt), 0.0)
            self.h -= take
            self.ledger.book("sink", -take.sum() * A)
        for lv in f.levels:
            m = lv.mask & o
            target = np.maximum(_val(lv.level, t) - self.z[m], 0.0)
            d = (target - self.h[m]) * A
            self.ledger.book("level", d[d > 0].sum())        # in and out booked apart (both count as throughput)
            self.ledger.book("level", d[d < 0].sum())
            self.h[m] = target
        self.min_h_seen = min(self.min_h_seen, float(self.h.min()))

    def step(self, dt: float) -> None:
        self._step(dt)
        self.apply_forcing(dt)
        self.t += dt
        self.nsteps += 1

    # -- terrain edits that keep every cell's water volume
    def dig(self, mask: np.ndarray, depth: float) -> None:
        """Lower the bed by ``depth`` on ``mask``; the water in each cell stays (its surface drops)."""
        self.z[mask] -= depth

    def raise_bed(self, mask: np.ndarray, height: float) -> None:
        """Raise the bed on ``mask``; water in each cell stays and is lifted (it then flows away)."""
        self.z[mask] += height

    def diagnostics(self) -> dict:
        d = super().diagnostics()
        d.update(volume=self.volume(), volume_error=self.volume_error(),
                 throughput=self.ledger.throughput, min_h=self.min_h_seen,
                 max_h=float(self.h.max()), wet=int((self.h > 1e-3).sum()))
        return d

    def snapshot(self) -> dict:
        u, v = self.velocity()
        return {"z": self.z.copy(), "h": self.h.copy(), "wall": self.wall, "speed": np.hypot(u, v)}

    # -- shared helpers
    def face_surface_diff(self):
        """eta differences on interior faces: dEx = eta[i+1]-eta[i] (ny, nx-1), dEy (ny-1, nx)."""
        e = self.eta
        return e[:, 1:] - e[:, :-1], e[1:, :] - e[:-1, :]

    def face_flow_depth(self):
        """Bates (2010) flow depth on faces: max(eta) - max(z), >= 0, zero across walls."""
        e, z = self.eta, self.z
        hx = np.maximum(np.maximum(e[:, 1:], e[:, :-1]) - np.maximum(z[:, 1:], z[:, :-1]), 0.0)
        hy = np.maximum(np.maximum(e[1:, :], e[:-1, :]) - np.maximum(z[1:, :], z[:-1, :]), 0.0)
        return hx * self.open_x, hy * self.open_y

    def wave_dt(self, cfl: float = 0.7) -> float:
        hmax = float(self.h.max())
        return cfl * self.grid.dx / math.sqrt(G * max(hmax, 1e-3))
