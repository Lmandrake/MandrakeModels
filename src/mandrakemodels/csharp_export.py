"""Export water scenes and a local-inertial reference run for the C# console port (stage C).

    python -m mandrakemodels.csharp_export --out DIR --sizes 128,250 --scenes river,lake_at_rest,pit_dig,dam_break,swell \
        [--solver local_inertial|local_inertial_frcap|local_inertial_frcap07]

Writes ``DIR/<scene>_<n>.lis`` (little-endian binary, read by ``csharp/LiTiming/SceneFile.cs``): the scene's bed, depth,
walls, forcing and events, then the reference solver's depth at a few checkpoint times. The file records the reference's face Froude cap
(0 = none), so the C# side runs the same scheme. Both sides choose dt by the
same rule (bakeoff ``choose_dt`` plus landing exactly on each checkpoint), so the comparison is step for step.
Callable forcing is encoded by kind: the swell level is a sine (amplitude, period); any other callable is refused.
"""
from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path

import numpy as np

from . import scenes
from .bakeoff import WATER_DT_CAP
from .solvers.local_inertial import SOLVERS

REFERENCES = {c.name: c for c in SOLVERS}
MAGIC = b"LIS2"          # LIS1 plus the face Froude cap after alpha
CHECKPOINTS = 4          # evenly spaced over the duration, the last at the end


def _mask(f, m):
    f.write(np.ascontiguousarray(m, dtype=np.uint8).tobytes())


def _arr(f, a):
    f.write(np.ascontiguousarray(a, dtype="<f8").tobytes())


def _level_kind(sc, lv):
    if not callable(lv.level):
        return 0, float(lv.level), 0.0
    if sc.name in ("swell", "stress_half"):
        return 1, 0.3, 8.0                     # scenes.swell and stress_half: 0.3 * sin(2 pi t / 8)
    raise SystemExit(f"{sc.name}: callable level not encodable")


def fr_max_of(cls) -> float:
    return float(getattr(cls, "fr_max", 0.0))


def reference(sc, times, cls):
    """Run the reference solver with the bakeoff dt rule, landing on each checkpoint time; return h at each."""
    s = cls.from_scene(sc)
    events = sorted(sc.events, key=lambda e: e[0])
    out, steps = [], 0
    for tc in times:
        while s.t < tc - 1e-9:
            while events and s.t >= events[0][0] - 1e-9:
                _, meth, kw = events.pop(0)
                getattr(s, meth)(**kw)
            dt = min(s.stable_dt(), WATER_DT_CAP, tc - s.t)
            if events:
                nxt = events[0][0] - s.t
                if nxt > 1e-9:
                    dt = min(dt, nxt)
            s.step(dt)
            steps += 1
        out.append((s.t, s.h.copy(), steps, s.volume_error()))
    return out


# Load scenes for timing only (no reference run). Wet targets per size are the plan's census figures (section 5.1).
CENSUS_RIVER_WET = {250: 13300, 325: 17300, 400: 21200}


def census_river(n, seed=0):
    """A river whose wet area matches the plan's river-map census for this size (one wide meandering channel)."""
    target = CENSUS_RIVER_WET.get(n, int(13300 * n / 250))
    w = max(6, round(target / n))
    z, chan, w = scenes._river_bed(n, seed, width=w)
    inlet = chan & (np.arange(n)[None, :] < 2)
    outlet = chan & (np.arange(n)[None, :] >= n - 2)
    h = np.where(chan, 1.0 - (z - z.min(axis=0, keepdims=True)), 0.0)
    f = scenes.Forcing(sources=[scenes.Source(inlet, 0.6 * w * scenes.DX)], sinks=[scenes.Sink(outlet)])
    return scenes._water("census_river", n, 120.0, z, h, forcing=f, note=f"river sized to {target} wet squares")


def stress_half(n, seed=0):
    """Half the map under water up to 10 m deep, the west edge swelling 0.3 m every 8 s: the plan's 50%-wet stress."""
    X = np.arange(n)[None, :] / n * np.ones((n, 1))
    z = -10.0 + 20.0 * X + 0.3 * scenes.smooth_noise(n, scale=0.05, seed=seed)
    sea = np.arange(n)[None, :] < 2
    h = np.maximum(0.0 - z, 0.0)
    f = scenes.Forcing(levels=[scenes.FixedLevel(sea & np.ones((n, 1), bool), lambda t: 0.3 * math.sin(2 * math.pi * t / 8.0))])
    return scenes._water("stress_half", n, 60.0, z, h, forcing=f, settle_tol=0.0, note="50% wet, 10 m deep, swelling edge")


LOCAL = {"census_river": census_river, "stress_half": stress_half}


def export(name, n, out_dir: Path, seed=0, ref=True, solver="local_inertial"):
    cls = REFERENCES[solver]
    sc = LOCAL[name](n, seed) if name in LOCAL else scenes.build(name, n, seed)
    if callable(sc.forcing.rain) or sc.forcing.groundwater is not None and not np.isscalar(sc.forcing.groundwater.table):
        raise SystemExit(f"{name}: rain/groundwater form not encodable")
    p = {"manning": 0.03, "theta": 0.8, "alpha": 0.6, **sc.params}
    times = [sc.duration * (k + 1) / CHECKPOINTS for k in range(CHECKPOINTS)]
    ref = reference(sc, times, cls) if ref else []
    path = out_dir / f"{name}_{n}.lis"
    g, fo = sc.grid, sc.forcing
    with path.open("wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<ii", g.ny, g.nx))
        f.write(struct.pack("<6d", g.dx, p["manning"], p["theta"], p["alpha"], fr_max_of(cls), sc.duration))
        _arr(f, sc.init["z"]); _arr(f, sc.init["h"]); _mask(f, sc.init["wall"])
        f.write(struct.pack("<dd", float(fo.evaporation), float(fo.rain or 0.0)))
        f.write(struct.pack("<i", len(fo.sources)))
        for s in fo.sources:
            if callable(s.Q):
                raise SystemExit(f"{name}: callable source not encodable")
            _mask(f, s.mask); f.write(struct.pack("<d", float(s.Q)))
        f.write(struct.pack("<i", len(fo.sinks)))
        for s in fo.sinks:
            _mask(f, s.mask); f.write(struct.pack("<d", float(s.rate)))
        f.write(struct.pack("<i", len(fo.levels)))
        for lv in fo.levels:
            _mask(f, lv.mask); f.write(struct.pack("<idd", *_level_kind(sc, lv)))
        gw = fo.groundwater
        f.write(struct.pack("<i", 0 if gw is None else 1))
        if gw is not None:
            f.write(struct.pack("<dd", float(gw.table), float(gw.k)))
        f.write(struct.pack("<i", len(sc.events)))
        for t, meth, kw in sorted(sc.events, key=lambda e: e[0]):
            if meth != "dig":
                raise SystemExit(f"{name}: event {meth} not encodable")
            depth = np.zeros(g.shape)
            depth[kw["mask"]] = kw["depth"]
            f.write(struct.pack("<d", float(t))); _mask(f, kw["mask"]); _arr(f, depth)
        f.write(struct.pack("<i", len(ref)))
        for t, h, steps, verr in ref:
            f.write(struct.pack("<did", t, steps, verr)); _arr(f, h)
    return path, ref[-1][2] if ref else 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sizes", default="128,250")
    ap.add_argument("--scenes", default="river,lake_at_rest,pit_dig,dam_break,swell")
    ap.add_argument("--no-ref", action="store_true", help="scene only, no Python reference checkpoints (timing loads)")
    ap.add_argument("--solver", default="local_inertial", choices=sorted(REFERENCES),
                    help="reference solver; its face Froude cap is written into the scene (timing loads too)")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for n in [int(x) for x in a.sizes.split(",")]:
        for name in a.scenes.split(","):
            path, steps = export(name, n, out, ref=not a.no_ref, solver=a.solver)
            print(f"{path.name}: {steps} reference steps", flush=True)


if __name__ == "__main__":
    main()
