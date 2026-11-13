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
