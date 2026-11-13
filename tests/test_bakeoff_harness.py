"""Harness fixes from the stage-2 review (RIMMANDRAKE_HARDEN_FAST_WATER_SOLVERS_VIA_GPT_PIPES_RING)."""
import math
import time

from mandrakemodels import bakeoff
from mandrakemodels.solvers import REGISTRY


def test_timing_includes_stable_dt(monkeypatch):
    """The timed region of a water step covers stable_dt too (wca2d and kinematic recompute velocities there)."""
    base = REGISTRY["local_inertial"]

    class Slow(base):
        name = "slow_dt_probe"

        def stable_dt(self):
            time.sleep(0.01)
            return super().stable_dt()

    monkeypatch.setitem(REGISTRY, Slow.name, Slow)
    r = bakeoff.timing_run(Slow.name, "lake_at_rest", 32, steps=5, warm=1, max_wall=30)
    assert r["ms_per_step"] >= 10.0, r


def test_maxdt_search_respects_its_cap(monkeypatch):
    """The doubling search never reports a dt above hi_cap (it reported wca2d at 81.7 s past a 64 s cap)."""
    monkeypatch.setattr(bakeoff, "survives", lambda cls, kind, dt: True)
    v, note = bakeoff.max_dt(REGISTRY["wca2d"], "water", own=0.16, hi_cap=64.0)
    assert v <= 64.0, (v, note)


def test_maxdt_survival_run_stops_at_T(monkeypatch):
    """survives() simulates exactly T seconds: the last step is clipped, never ceil(T/dt) full steps."""
    seen = []
    cls = REGISTRY["local_inertial"]

    class Spy(cls):
        def step(self, dt):
            seen.append(dt)
            return super().step(dt)

    bakeoff.survives(Spy, "water", 7.0)
    assert math.isclose(sum(seen), bakeoff.MAXDT_SCENE["water"][2], rel_tol=1e-9), sum(seen)


def test_report_projects_no_csharp_cost_from_the_simplified_numba_probe(tmp_path):
    """The numba probe's kernels are different solvers (pipes without friction), so no C# ms/tick is derived."""
    import json
    from mandrakemodels import report
    (tmp_path / "probe.json").write_text(json.dumps({"n": 32, "pipes": {"numpy_ms": 1.0, "numba_ms": 0.1, "ratio": 10.0,
                                                                         "numba_ns_per_cell": 1.0}}))
    (tmp_path / "timing.jsonl").write_text(json.dumps(dict(solver="pipes", scene="river", n=325, kind="water",
                                                           ms_per_step=9.0, dt_mean=0.15, py_ms_per_tick=1.0)) + "\n")
    rows, _, _ = report.merged(tmp_path)
    assert not any(k.startswith("cs_") for k in rows[0]), rows[0]
    assert "C# ms/tick guess" not in report.md_headline(report.headline(rows, 325), 325)


def test_ledger_throughput_counts_an_exchange_both_ways():
    """A level-setter that adds 5 m^3 then takes 5 m^3 has 10 m^3 of throughput, not 0 (|net per kind| cancelled)."""
    from mandrakemodels.base import Ledger
    L = Ledger()
    L.book("level", 5.0)
    L.book("level", -5.0)
    assert L.throughput == 10.0, L.throughput


def test_rain_reaches_the_step_that_ends_at_its_cutoff():
    """Rain while t < 1 s: the step from 0.5 to 1.0 s gets its rain (it was sampled at the step's end, so got none)."""
    import numpy as np
    from mandrakemodels.base import Forcing, Grid
    g = Grid(2, 2, 1.0)
    s = REGISTRY["local_inertial"](g, np.zeros((2, 2)), np.zeros((2, 2)),
                                   forcing=Forcing(rain=lambda t: 1e-3 if t < 1.0 else 0.0))
    s.t = 0.5
    s.step(0.5)
    assert abs(s.h.mean() - 5e-4) < 1e-12, s.h


def test_water_headline_leaves_out_kinematic_wave():
    """kinematic_wave is out of the water headline (bed-slope only; it blew up in ocean_edge and swell)."""
    from mandrakemodels import report
    rows = [dict(solver=k, scene="river", n=325, kind="water", ms_per_step=1.0) for k in ("kinematic_wave", "pipes")]
    assert [r["solver"] for r in report.headline(rows, 325)] == ["pipes"]
