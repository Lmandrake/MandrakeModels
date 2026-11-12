"""Scenes as data. ``build(name, n, seed)`` returns a ``Scene``: grid, initial arrays, forcing, timed events,
duration and settle tolerance. Same inputs give the same scene (seeded). Water scenes use dx = 1.414 m
(a RimWorld square is 2 m^2); the map is ``n x n`` cells, so bigger sizes are bigger places.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .base import FixedLevel, Forcing, Grid, Groundwater, Sink, Source

DX = 1.414


@dataclass
class Scene:
    name: str
    kind: str                      # water | wave | scalar | flow
    grid: Grid
    duration: float                # simulated seconds
    init: dict                     # arrays the solver starts from
    forcing: Forcing = field(default_factory=Forcing)
    events: list = field(default_factory=list)   # (t, method name, kwargs)
    params: dict = field(default_factory=dict)   # scene-level physics knobs (manning, c, nu, ...)
    settle_tol: float = 1e-4       # water: max |dh/dt| m/s counted as settled
    still: bool = False            # True: nothing should move (lake at rest)
    note: str = ""


def smooth_noise(n: int, ny: int | None = None, scale: float = 0.08, seed: int = 0) -> np.ndarray:
    """Zero-mean, unit-std smooth random field; ``scale`` is the feature size as a fraction of the map."""
    ny = ny or n
    rng = np.random.default_rng(seed)
    w = rng.standard_normal((ny, n))
    ky = np.fft.fftfreq(ny)[:, None]
    kx = np.fft.fftfreq(n)[None, :]
    k2 = kx * kx + ky * ky
    sigma = scale * max(n, ny)                      # Gaussian smoothing width in cells
    filt = np.exp(-2 * math.pi ** 2 * sigma ** 2 * k2)
    f = np.real(np.fft.ifft2(np.fft.fft2(w) * filt))
    f -= f.mean()
    return f / (f.std() + 1e-12)


def _water(name, n, duration, z, h, wall=None, forcing=None, events=None, still=False, note="",
           settle_tol=1e-4, manning=0.03):
    g = Grid(n, n, DX)
    wall = np.zeros((n, n), bool) if wall is None else wall
    h = np.where(wall, 0.0, np.maximum(h, 0.0))
    return Scene(name, "water", g, duration, {"z": z, "h": h, "wall": wall}, forcing or Forcing(),
                 events or [], {"manning": manning}, settle_tol, still, note)


def _channel(n, center_frac, width, depth, meander=0.0, seed=0):
    """Boolean channel mask along x and its centreline; returns (mask, dist to centreline in cells)."""
    j = np.arange(n)[:, None]
    i = np.arange(n)[None, :]
    c = center_frac * n + meander * n * np.sin(2 * math.pi * i / n * 1.5 + seed)
    d = np.abs(j - c)
    return d <= width / 2, d


# ---------------------------------------------------------------------------------------------- water

def lake_at_rest(n, seed=0):
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    z = 3.0 * (X ** 2 + Y ** 2) + 0.6 * smooth_noise(n, scale=0.05, seed=seed)
    z += 2.5 * np.exp(-((X - 0.3) ** 2 + (Y + 0.2) ** 2) / 0.02)      # an island
    level = 1.5
    return _water("lake_at_rest", n, 60.0, z, level - z, still=True, settle_tol=1e-6,
                  note="closed lake on a rough bed with an island: nothing may move")


def dam_break(n, seed=0):
    z = np.zeros((n, n))
    wall = np.zeros((n, n), bool)
    x0 = n // 3
    wall[:, x0:x0 + 2] = True
    gap = slice(n // 2 - n // 8, n // 2 + n // 8)
    wall[gap, x0:x0 + 2] = False
    h = np.where(np.arange(n)[None, :] < x0, 2.0, 0.3) * np.ones((n, 1))
    return _water("dam_break", n, 40.0, z, h, wall, note="partial dam breach, 2 m onto 0.3 m wet bed")


def _river_bed(n, seed, slope=1e-3, width=None, meander=0.08, center=0.5):
    width = width or max(6, n // 18)
    X = np.arange(n)[None, :] * DX
    bank = 2.0 + 0.4 * smooth_noise(n, scale=0.06, seed=seed)
    chan, d = _channel(n, center, width, 1.5, meander, seed)
    bottom = 0.0 + 0.3 * (d / (width / 2)) ** 2
    z = np.where(chan, bottom, bank + 0.02 * np.maximum(d - width / 2, 0)) - slope * X
    return z, chan, width


def river(n, seed=0):
    z, chan, w = _river_bed(n, seed)
    inlet = chan & (np.arange(n)[None, :] < 2)
    outlet = chan & (np.arange(n)[None, :] >= n - 2)
    h = np.where(chan, 1.0 - (z - z.min(axis=0, keepdims=True)), 0.0)
    f = Forcing(sources=[Source(inlet, 0.6 * w * DX)], sinks=[Sink(outlet)])
    return _water("river", n, 120.0, z, h, forcing=f,
                  note="meandering sloped channel, inflow source upstream, free-outfall sink downstream")


def lake_levelset(n, seed=0):
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    z = 2.0 * (X ** 2 + Y ** 2) + 0.3 * smooth_noise(n, scale=0.05, seed=seed)
    chan = (np.abs(Y) < max(3, n // 40) / n * 2) & (X > 0)
    z = np.where(chan, np.minimum(z, 0.2 + 0.5 * X), z)            # outlet channel to the east edge
    edge = chan & (np.arange(n)[None, :] >= n - 2)
    h = 1.2 - z
    h[~(np.hypot(X, Y) < 0.8)] = 0
    f = Forcing(levels=[FixedLevel(edge, 0.9)])
    return _water("lake_levelset", n, 120.0, z, h, forcing=f,
                  note="lake at 1.2 m drains through a channel to an edge level-setter held at 0.9 m")


def ocean_edge(n, seed=0):
    X = np.arange(n)[None, :] / n * np.ones((n, 1))
    z = -6.0 + 10.0 * X + 0.4 * smooth_noise(n, scale=0.05, seed=seed)
    chan, d = _channel(n, 0.55, max(5, n // 25), 1.2, 0.05, seed)
    z = np.where(chan & (X > 0.45), np.minimum(z, 0.2 + 2.0 * (X - 0.45)), z)
    sea = np.arange(n)[None, :] < 3
    inlet = chan & (np.arange(n)[None, :] >= n - 2)
    h = np.maximum(0.0 - z, 0.0)
    f = Forcing(levels=[FixedLevel(sea & np.ones((n, 1), bool), 0.0)],
                sources=[Source(inlet, 0.5 * max(5, n // 25) * DX)])
    return _water("ocean_edge", n, 120.0, z, h, forcing=f,
                  note="beach sloping into a sea held at 0 m on the west edge; river enters from the east")


def pit_dig(n, seed=0):
    z, chan, w = _river_bed(n, seed, meander=0.0, center=0.35)
    J, I = np.mgrid[0:n, 0:n]
    pit = (np.abs(J - 0.72 * n) < n * 0.12) & (np.abs(I - 0.5 * n) < n * 0.15)
    canal = (np.abs(I - 0.5 * n) < max(2, w // 3)) & (J > 0.35 * n) & (J < 0.72 * n)
    inlet = chan & (I < 2)
    outlet = chan & (I >= n - 2)
    h = np.where(chan, 1.0 - (z - z.min(axis=0, keepdims=True)), 0.0)
    f = Forcing(sources=[Source(inlet, 0.6 * w * DX)], sinks=[Sink(outlet)])
    zc = z.copy()
    ev = [(5.0, "dig", {"mask": pit, "depth": 3.0}),
          (15.0, "dig", {"mask": canal & (zc > 0.5), "depth": 2.0})]
    return _water("pit_dig", n, 150.0, z, h, forcing=f, events=ev,
                  note="river; a 3 m pit dug at t=5 s, a canal from the river dug at t=15 s fills it")


def moat_evap(n, seed=0):
    J, I = np.mgrid[0:n, 0:n] / n
    r = np.maximum(np.abs(J - 0.5), np.abs(I - 0.5))
    moat = (r > 0.25) & (r < 0.33)
    z = np.where(moat, -1.0, 1.0) + 0.1 * smooth_noise(n, scale=0.05, seed=seed)
    h = np.where(moat, 0.5 - z, 0.0)
    h = np.maximum(h, 0)
    f = Forcing(evaporation=2e-3)
    return _water("moat_evap", n, 120.0, z, h, forcing=f,
                  note="closed square moat, evaporation 2 mm/s (exaggerated so it shows in a short run)")


def rain_groundwater(n, seed=0):
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    z = 1.5 * (X ** 2 + Y ** 2) + 0.8 * smooth_noise(n, scale=0.07, seed=seed)
    dur = 120.0
    f = Forcing(rain=lambda t: 4e-3 if t < 30 else 0.0, groundwater=Groundwater(table=0.2, k=2e-2))
    return _water("rain_groundwater", n, dur, z, np.zeros((n, n)), forcing=f,
                  note="dry rough basin: 4 mm/s rain for 30 s, groundwater table at 0.2 m seeps up after")


def swell(n, seed=0):
    X = np.arange(n)[None, :] / n * np.ones((n, 1))
    z = -8.0 + 10.0 * X + 0.3 * smooth_noise(n, scale=0.05, seed=seed)
    wall = np.zeros((n, n), bool)
    J = np.arange(n)[:, None] * np.ones((1, n))
    bw = (np.abs(X - 0.45) < 0.012) & ((J < 0.35 * n) | ((J > 0.45 * n) & (J < 0.75 * n)))
    wall |= bw
    sea = np.arange(n)[None, :] < 2
    h = np.maximum(0.0 - z, 0.0)
    f = Forcing(levels=[FixedLevel(sea & np.ones((n, 1), bool), lambda t: 0.3 * math.sin(2 * math.pi * t / 8.0))])
    return _water("swell", n, 60.0, z, h, wall, forcing=f, settle_tol=0.0,
                  note="sea edge rises and falls 0.3 m every 8 s onto a beach with two breakwaters")


# ---------------------------------------------------------------------------------- wave / scalar / flow

def ripples(n, seed=0):
    rng = np.random.default_rng(seed)
    X, Y = np.meshgrid(np.linspace(-1, 1, n), np.linspace(-1, 1, n))
    obstacle = (np.hypot(X - 0.25, Y + 0.1) < 0.18) | ((np.abs(X + 0.4) < 0.03) & (np.abs(Y) > 0.15))
    drops = [(float(t), float(x), float(y)) for t, x, y in
             zip(np.sort(rng.uniform(0, 15, 8)), rng.uniform(0.1, 0.9, 8), rng.uniform(0.1, 0.9, 8))]
    return Scene("ripples", "wave", Grid(n, n, 1.0), 20.0, {"obstacle": obstacle, "drops": drops},
                 params={"c": 10.0, "damping": 0.05}, note="pond with an island and a slotted wall; 8 drops")


def gaussian(n, seed=0):
    g = Grid(n, n, 1.0)
    X, Y = g.coords()
    L = n * g.dx
    sigma0 = 0.05 * L
    u0 = np.exp(-((X - L / 2) ** 2 + (Y - L / 2) ** 2) / (2 * sigma0 ** 2))
    D = 1.0
    T = 0.5 * (0.1 * L) ** 2 / D      # sigma grows to sqrt(sigma0^2 + 2 D T)
    return Scene("gaussian", "scalar", g, T, {"field": u0}, params={"D": D, "sigma0": sigma0},
                 note="Gaussian blob diffusing; exact solution known (insulated walls far away)")


def dye(n, seed=0):
    g = Grid(n, n, 1.0)
    rng = np.random.default_rng(seed)
    X, Y = g.coords()
    u = np.zeros(g.shape)
    for _ in range(12):
        cx, cy, s = rng.uniform(0.15, 0.85) * n, rng.uniform(0.15, 0.85) * n, rng.uniform(0.01, 0.04) * n
        u += rng.uniform(0.5, 1.0) * ((np.hypot(X - cx, Y - cy) < s * 1.5))
    D = 1.0
    return Scene("dye", "scalar", g, 0.5 * (0.06 * n) ** 2 / D, {"field": u}, params={"D": D},
                 note="12 sharp dye discs diffusing with insulated walls (total dye must be conserved)")


def smoke_plume(n, seed=0):
    g = Grid(n, n, 1.0)
    X, Y = np.meshgrid(np.linspace(0, 1, n), np.linspace(0, 1, n))
    obstacle = (np.hypot(X - 0.5, Y - 0.55) < 0.07) | ((np.abs(Y - 0.3) < 0.015) & (np.abs(X - 0.3) < 0.12))
    src = (np.hypot(X - 0.5, Y - 0.06) < 0.05)
    return Scene("smoke_plume", "flow", g, 6.0, {"obstacle": obstacle, "source": src},
                 params={"buoyancy": 1.0, "inflow": 0.0, "nu": 0.0, "steps": 3 * n},
                 note="buoyant smoke rising from a vent past a disc and a shelf (closed box)")


def cylinder_wake(n, seed=0):
    ny = max(16, n // 2)
    g = Grid(ny, n, 1.0)
    X, Y = np.meshgrid(np.arange(n), np.arange(ny))
    obstacle = np.hypot(X - n * 0.2, Y - ny * 0.5 - 0.5) < ny * 0.09
    src = np.zeros((ny, n), bool)
    src[(np.arange(ny) % max(4, ny // 8)) < 2, 0:2] = True
    return Scene("cylinder_wake", "flow", g, float(4 * n), {"obstacle": obstacle, "source": src},
                 params={"buoyancy": 0.0, "inflow": 0.08, "Re": 150.0, "steps": 8 * n},
                 note="channel flow past a cylinder (Re ~150): Karman vortex street; dye streaks from the inlet")


WATER = ["lake_at_rest", "dam_break", "river", "lake_levelset", "ocean_edge", "pit_dig", "moat_evap",
         "rain_groundwater", "swell"]
OTHER = ["ripples", "gaussian", "dye", "smoke_plume", "cylinder_wake"]
ALL = WATER + OTHER


def build(name: str, n: int, seed: int = 0) -> Scene:
    return globals()[name](n, seed)
