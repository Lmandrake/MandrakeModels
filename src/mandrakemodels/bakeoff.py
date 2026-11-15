"""The bakeoff: every applicable solver x scene x size. Resumable; each phase appends JSON lines.

    python -m mandrakemodels.bakeoff quality --sizes 64,128,250 --workers 10 --budget 480
    python -m mandrakemodels.bakeoff timing  --sizes 64,128,250,325,400,500 --budget 480
    python -m mandrakemodels.bakeoff maxdt   --budget 480
    python -m mandrakemodels.bakeoff probe                     # numba unit cost for the C# guess
    python -m mandrakemodels.bakeoff report                    # results.csv, results.md, gallery.html

Phases
------
quality  full-duration runs (wall-capped per run, flagged ``capped``); settle, volume error, blow-ups,
         still-water drift; renders a GIF and stills when ``n == --frame-size``. Runs in parallel
         processes, so its ms/step is NOT used for timing.
timing   sequential, one process, single-threaded: warm-up then ``--steps`` timed steps per run.
maxdt    largest fixed dt that survives 20 s of the 128 dam break (water), 64 gaussian (scalar),
         64 ripples (wave), by doubling then bisection.
probe    numba kernels of three solvers, to turn numpy ms into a compiled-code guess.
Output goes under ``paths.out_root()/bakeoff/<run>`` (default run name ``main``).
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import json
import math
import statistics
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from . import render, scenes
from .base import even_dt
from .paths import out_root
from .solvers import REGISTRY, for_kind

SIZES = [64, 128, 250, 325, 400, 500]
TICKS_PER_SIM_SECOND = 60          # sibling plan's assumption: one river-second per 60 game ticks
WATER_DT_CAP = 5.0


def run_dir(name: str) -> Path:
    d = out_root() / "bakeoff" / name
    d.mkdir(parents=True, exist_ok=True)
    return d


def tasks(sizes, only_solvers=None, only_scenes=None):
    out = []
    for kind, names in (("water", scenes.WATER), ("other", scenes.OTHER)):
        for sc_name in names:
            if only_scenes and sc_name not in only_scenes:
                continue
            probe = scenes.build(sc_name, 32)
            for cls in for_kind(probe.kind):
                if only_solvers and cls.name not in only_solvers:
                    continue
                for n in sizes:
                    out.append((cls.name, sc_name, n))
    return out


def done_keys(path: Path) -> set:
    if not path.exists():
        return set()
    keys = set()
    for line in path.read_text().splitlines():
        try:
            r = json.loads(line)
            keys.add((r["solver"], r["scene"], r["n"]))
        except Exception:
            pass
    return keys


def _num(d):
    o = {}
    for k, v in d.items():
        if isinstance(v, (bool, np.bool_)):
            o[k] = bool(v)
        elif isinstance(v, (int, float, np.integer, np.floating)):
            o[k] = float(v) if not isinstance(v, (int, np.integer)) else int(v)
    return o


def choose_dt(s, sc, events_left):
    if sc.kind == "water":
        end = sc.duration
        if events_left and events_left[0][0] - s.t > 1e-9:
            end = min(end, events_left[0][0])
        return max(even_dt(min(s.stable_dt(), WATER_DT_CAP), end - s.t), 1e-6)
    if sc.kind == "flow":
        return s.stable_dt()
    return min(s.stable_dt(), sc.duration / 50.0, sc.duration - s.t)


def n_steps_total(sc):
    return int(sc.params.get("steps", 0)) if sc.kind == "flow" else None


# ------------------------------------------------------------------------------------------- quality

def quality_run(solver: str, scene: str, n: int, max_wall: float, frame_size: int, out: Path):
    cls = REGISTRY[solver]
    sc = scenes.build(scene, n)
    s = cls.from_scene(sc)
    events = sorted(sc.events, key=lambda e: e[0])
    last_event = events[-1][0] if events else 0.0
    want_frames = n == frame_size
    frames, stills = [], []
    style = None
    flow_steps = n_steps_total(sc)
    n_frames = 48
    t_frames = np.linspace(0, sc.duration, n_frames) if flow_steps is None else None
    s_frames = np.linspace(0, flow_steps, n_frames).astype(int) if flow_steps else None
    fi = 0

    def grab():
        nonlocal style
        snap = s.snapshot()
        if style is None:
            style = render.Style.from_snapshot(sc.kind, snap)
        return render.frame(snap, style, sc.grid.dx)

    rec = dict(solver=solver, scene=scene, n=n, kind=sc.kind, family=cls.family, conserves=cls.conserves)
    h0 = s.h.copy() if sc.kind == "water" else None
    last_unsettled_t, last_rate = 0.0, math.inf
    blow_t = None
    t0 = time.perf_counter()
    step_time = 0.0
    capped = False
    err = None
    try:
        while True:
            if flow_steps is not None:
                if s.nsteps >= flow_steps:
                    break
            elif s.t >= sc.duration - 1e-9:
                break
            while events and s.t >= events[0][0] - 1e-9:
                _, meth, kw = events.pop(0)
                getattr(s, meth)(**kw)
            if want_frames:
                if flow_steps is None and fi < n_frames and s.t >= t_frames[fi] - 1e-9:
                    frames.append(grab()); fi += 1
                elif flow_steps is not None and fi < n_frames and s.nsteps >= s_frames[fi]:
                    frames.append(grab()); fi += 1
            dt = choose_dt(s, sc, events)
            hb = s.h.copy() if sc.kind == "water" else None
            a = time.perf_counter()
            s.step(dt)
            step_time += time.perf_counter() - a
            if s.blown_up():
                blow_t = s.t
                break
            if sc.kind == "water":
                rate = float(np.abs(s.h - hb).max()) / dt
                last_rate = rate
                if rate >= sc.settle_tol or s.t <= last_event:
                    last_unsettled_t = s.t
                    rec["_settle_steps"] = s.nsteps
            if time.perf_counter() - t0 > max_wall:
                capped = True
                break
    except Exception as e:  # a solver that throws is reported, not dropped
        err = f"{type(e).__name__}: {e}"
        traceback.print_exc()
    rec.update(steps=s.nsteps, t_end=s.t, wall_s=time.perf_counter() - t0,
               ms_per_step_parallel=1000 * step_time / max(s.nsteps, 1), capped=capped,
               blew_up=blow_t is not None, blow_t=blow_t, error=err)
    if sc.kind == "water":
        rec["settle_tol"] = sc.settle_tol
        rec["final_rate"] = last_rate if math.isfinite(last_rate) else None
        settled = (not capped) and blow_t is None and err is None and last_rate < sc.settle_tol
        rec["settled"] = bool(settled)
        rec["settle_t"] = last_unsettled_t if settled else None
        rec["settle_steps"] = rec.pop("_settle_steps", 0) if settled else None
        rec.pop("_settle_steps", None)
        rec["vol_err_rel"] = abs(s.volume_error()) / max(s.V0, s.ledger.throughput, 1e-12)
        rec["min_h"] = float(s.min_h_seen)
        rec["max_h"] = float(np.nanmax(s.h))
        rec["wet_frac"] = float((s.h > 1e-3).mean())
        if sc.still:
            rec["still_drift_m"] = float(np.abs(s.h - h0).max())
        u, v = s.velocity()
        rec["max_speed"] = float(np.nanmax(np.hypot(u, v)))
    rec["diag"] = _num(s.diagnostics())
    if want_frames and err is None:
        frames.append(grab())
        d = out / "frames" / scene
        d.mkdir(parents=True, exist_ok=True)
        stem = f"{solver}_{n}"
        render.save_gif(frames, d / f"{stem}.gif", fps=12, max_px=320)
        for k, idx in (("a", len(frames) // 3), ("b", 2 * len(frames) // 3), ("c", len(frames) - 1)):
            render.save_png(frames[idx], d / f"{stem}_{k}.png", max_px=480)
        rec["gif"] = f"frames/{scene}/{stem}.gif"
    return rec


def _quality_worker(args):
    solver, scene, n, max_wall, frame_size, out = args
    try:
        return quality_run(solver, scene, n, max_wall, frame_size, Path(out))
    except Exception as e:
        return dict(solver=solver, scene=scene, n=n, error=f"{type(e).__name__}: {e}", blew_up=False)


def phase_quality(a):
    import multiprocessing as mp
    out = run_dir(a.run)
    path = out / "quality.jsonl"
    done = done_keys(path)
    todo = [t for t in tasks(a.sizes, a.solvers, a.scenes) if t not in done]
    todo.sort(key=lambda t: t[2])
    print(f"quality: {len(done)} done, {len(todo)} to do", flush=True)
    t0 = time.time()
    ctx = mp.get_context("fork")
    with ctx.Pool(a.workers, maxtasksperchild=4) as pool:
        it = pool.imap_unordered(_quality_worker,
                                 [(s, sc, n, a.max_wall, a.frame_size, str(out)) for s, sc, n in todo])
        k = 0
        while True:
            left = a.budget - (time.time() - t0)
            if left <= 0:
                print("quality: budget spent, stopping (resume to continue)", flush=True)
                pool.terminate()
                break
            try:
                r = it.next(timeout=left)
            except StopIteration:
                break
            except mp.TimeoutError:
                print("quality: budget spent, stopping (resume to continue)", flush=True)
                pool.terminate()
                break
            with path.open("a") as f:
                f.write(json.dumps(r) + "\n")
            k += 1
            if k % 10 == 0:
                print(f"quality: {k}/{len(todo)} at {time.time() - t0:.0f}s", flush=True)
    print(f"quality: wrote {k} this call", flush=True)


# -------------------------------------------------------------------------------------------- timing

def timing_run(solver, scene, n, steps, warm, max_wall):
    cls = REGISTRY[solver]
    sc = scenes.build(scene, n)
    s = cls.from_scene(sc)
    events = sorted(sc.events, key=lambda e: e[0])
    rec = dict(solver=solver, scene=scene, n=n, kind=sc.kind, cells=sc.grid.ny * sc.grid.nx)
    times, dts = [], []
    t0 = time.perf_counter()
    for i in range(warm + steps):
        while events and s.t >= events[0][0] - 1e-9:
            _, meth, kw = events.pop(0)
            getattr(s, meth)(**kw)
        a = time.perf_counter()        # stable_dt is part of the step's cost (wca2d, kinematic recompute velocities)
        dt = choose_dt(s, sc, events) if sc.kind != "water" else min(s.stable_dt(), WATER_DT_CAP)
        s.step(dt)
        b = time.perf_counter() - a
        if i >= warm:
            times.append(b)
            dts.append(dt)
        if s.blown_up() or time.perf_counter() - t0 > max_wall:
            break
    rec.update(timed_steps=len(times), blew_up=bool(s.blown_up()))
    if times:
        ms = 1000 * statistics.median(times)
        dtm = statistics.mean(dts)
        rec.update(ms_per_step=ms, ms_per_step_mean=1000 * statistics.mean(times), dt_mean=dtm,
                   ns_per_cell_step=1e6 * ms / rec["cells"])
        if sc.kind == "water":
            rec["py_ms_per_tick"] = ms / dtm / TICKS_PER_SIM_SECOND
    if "iters" in s.stats:
        rec["iters"] = s.stats["iters"]
    return rec


def phase_timing(a):
    out = run_dir(a.run)
    path = out / "timing.jsonl"
    done = done_keys(path)
    todo = [t for t in tasks(a.sizes, a.solvers, a.scenes) if t not in done]
    todo.sort(key=lambda t: t[2])
    print(f"timing: {len(done)} done, {len(todo)} to do", flush=True)
    t0 = time.time()
    k = 0
    for solver, scene, n in todo:
        if time.time() - t0 > a.budget:
            print("timing: budget spent (resume to continue)", flush=True)
            break
        try:
            r = timing_run(solver, scene, n, a.steps, a.warm, a.max_wall)
        except Exception as e:
            r = dict(solver=solver, scene=scene, n=n, error=f"{type(e).__name__}: {e}")
        with path.open("a") as f:
            f.write(json.dumps(r) + "\n")
        k += 1
    print(f"timing: wrote {k} this call", flush=True)


# --------------------------------------------------------------------------------------------- maxdt

MAXDT_SCENE = {"water": ("dam_break", 128, 20.0), "scalar": ("gaussian", 64, None),
               "wave": ("ripples", 64, 10.0)}


def survives(cls, kind, dt):
    sc_name, n, T = MAXDT_SCENE[kind]
    sc = scenes.build(sc_name, n)
    T = T or sc.duration
    s = cls.from_scene(sc)
    ref = None
    if kind == "water":
        ref = float(s.h.max())
    elif kind == "scalar":
        ref = float(s.snapshot()["field"].max())
    steps = int(math.ceil(T / dt - 1e-9))
    if steps > 40000:
        return None
    for k in range(steps):
        s.step(min(dt, T - k * dt))           # the last step is clipped: exactly T seconds are simulated
        if s.blown_up():
            return False
    snap = s.snapshot()
    if kind == "water":
        return bool(np.isfinite(s.h).all() and s.h.max() <= 1.25 * ref)
    if kind == "scalar":
        f = snap["field"]
        return bool(np.isfinite(f).all() and f.max() <= ref * 1.0001 and f.min() >= -1e-6 * ref)
    u = snap["height"]
    return bool(np.isfinite(u).all() and np.abs(u).max() < 10.0)


def max_dt(cls, kind, own=math.inf, lo=1e-3, hi_cap=64.0):
    """Largest fixed dt that survives, searched from the solver's own stable dt (0.5 s when it has
    none): halve until it survives, double until it fails, then bisect 6 times in log space."""
    if cls.name in ("ca_game",):
        return None, "n/a: fixed pass"
    t = min(max(own, 0.005), 2.0) if math.isfinite(own) else 0.5
    while True:
        r = survives(cls, kind, t)
        if r:
            break
        if r is None or t / 2 < lo:
            return None, f"< {t:g} (fails at the smallest dt tried)" if r is not None else \
                f"not searched: < {t:g} needs over 40000 steps"
        t /= 2
    good, bad = t, None
    while t < hi_cap:
        t = min(2 * t, hi_cap)
        r = survives(cls, kind, t)
        if r is None:
            break
        if not r:
            bad = t
            break
        good = t
    if cls.name == "kinematic_wave":
        return good, "meaningless here: the dam-break bed is flat, so kinematic flow never starts"
    if bad is None:
        return good, f">= {good:g} (cap)"
    for _ in range(6):
        mid = math.sqrt(good * bad)
        if survives(cls, kind, mid):
            good = mid
        else:
            bad = mid
    return good, ""


def phase_maxdt(a):
    out = run_dir(a.run)
    path = out / "maxdt.jsonl"
    done = {r.split('"solver": "')[1].split('"')[0] for r in path.read_text().splitlines()} if path.exists() else set()
    t0 = time.time()
    for cls in sorted(REGISTRY.values(), key=lambda c: c.name):
        if cls.kind not in MAXDT_SCENE or cls.name in done:
            continue
        if a.solvers and cls.name not in a.solvers:
            continue
        if time.time() - t0 > a.budget:
            print("maxdt: budget spent (resume)", flush=True)
            break
        a0 = time.time()
        sc_name, n, _ = MAXDT_SCENE[cls.kind]
        s = cls.from_scene(scenes.build(sc_name, n))
        own = s.stable_dt()
        try:
            v, note = max_dt(cls, cls.kind, own)
        except Exception as e:
            v, note = None, f"error {type(e).__name__}: {e}"
        r = dict(solver=cls.name, kind=cls.kind, scene=sc_name, n=n, max_dt=v, note=note,
                 own_stable_dt=None if math.isinf(own) else own, wall_s=time.time() - a0)
        with path.open("a") as f:
            f.write(json.dumps(r) + "\n")
        print(f"maxdt {cls.name}: {v} {note} ({time.time() - a0:.0f}s)", flush=True)


# --------------------------------------------------------------------------------------------- probe

def phase_probe(a):
    from . import probe_numba
    out = run_dir(a.run)
    res = probe_numba.run(n=a.probe_n)
    (out / "probe.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


# ------------------------------------------------------------------------------------------ accuracy

def breach_volume(cls, n=128, T=20.0):
    """m^3 that crosses the dam line of ``dam_break`` in T seconds (each solver at its own dt)."""
    sc = scenes.build("dam_break", n)
    c = int(np.where(sc.init["wall"].any(0))[0].max()) + 1
    s = cls.from_scene(sc)
    r0 = float(s.h[:, c:].sum())
    while s.t < T - 1e-9:
        s.step(even_dt(min(s.stable_dt(), WATER_DT_CAP), T - s.t))
        if s.blown_up():
            return None
    return (float(s.h[:, c:].sum()) - r0) * sc.grid.dx ** 2


def phase_accuracy(a):
    """Breach throughput of every water solver against swe_hll, the solver that matches Stoker's
    dam break to 0.1% (tests/test_water_main.py). A ratio far from 1 means wrong flood timing."""
    out = run_dir(a.run)
    res = {"scene": "dam_break", "n": 128, "T": 20.0, "reference": "swe_hll", "m3": {}}
    if a.solvers and (out / "accuracy.json").exists():     # re-run only the named rows, keep the others
        res = json.loads((out / "accuracy.json").read_text())
        res.pop("ratio", None)
    for cls in sorted(REGISTRY.values(), key=lambda c: c.name):
        if cls.kind != "water" or (a.solvers and cls.name not in a.solvers):
            continue
        t0 = time.time()
        try:
            res["m3"][cls.name] = breach_volume(cls)
        except Exception as e:
            res["m3"][cls.name] = None
            res.setdefault("errors", {})[cls.name] = f"{type(e).__name__}: {e}"
        print(f"accuracy {cls.name}: {res['m3'][cls.name]} ({time.time() - t0:.0f}s)", flush=True)
    ref = res["m3"].get("swe_hll")
    res["ratio"] = {k: (v / ref if (v is not None and ref) else None) for k, v in res["m3"].items()}
    (out / "accuracy.json").write_text(json.dumps(res, indent=1))


# -------------------------------------------------------------------------------------------- report

def load(path):
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines():
        try:
            rows.append(json.loads(line))
        except Exception:
            pass
    return rows


MONO_LO, MONO_HI = 1.5, 4.0     # guess: Unity Mono JIT vs numba/LLVM on the same loop


def phase_report(a):
    from . import report
    report.write(run_dir(a.run), Path(a.curate) if a.curate else None)


def main(argv=None):
    p = argparse.ArgumentParser(prog="mandrakemodels.bakeoff")
    p.add_argument("phase", choices=["quality", "timing", "maxdt", "probe", "accuracy", "report"])
    p.add_argument("--run", default="main")
    p.add_argument("--sizes", default=",".join(map(str, SIZES)))
    p.add_argument("--solvers", default="")
    p.add_argument("--scenes", default="")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--budget", type=float, default=480.0, help="wall seconds for this call")
    p.add_argument("--max-wall", type=float, default=150.0, help="wall cap per run (quality/timing)")
    p.add_argument("--frame-size", type=int, default=250)
    p.add_argument("--steps", type=int, default=30)
    p.add_argument("--warm", type=int, default=3)
    p.add_argument("--probe-n", type=int, default=400)
    p.add_argument("--curate", default="", help="report: copy small tables/images into this dir")
    a = p.parse_args(argv)
    a.sizes = [int(x) for x in a.sizes.split(",") if x]
    a.solvers = [x for x in a.solvers.split(",") if x] or None
    a.scenes = [x for x in a.scenes.split(",") if x] or None
    {"quality": phase_quality, "timing": phase_timing, "maxdt": phase_maxdt, "probe": phase_probe,
     "accuracy": phase_accuracy, "report": phase_report}[a.phase](a)


if __name__ == "__main__":
    main()
