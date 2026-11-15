"""Regression tests for the stage-2 review fixes (results/stage2.md).

Each defect is pinned twice where a solver fix became a new competitor: the original solver is a strict xfail (the
defect is still there, unchanged on purpose), the new competitor must pass. Harness fixes are plain tests.
"""
import numpy as np
import pytest

from mandrakemodels import bakeoff, scenes
from mandrakemodels.base import Grid
from mandrakemodels.solvers import REGISTRY

CAP = 5.0


def _run(s, sc, T):
    events = sorted(sc.events, key=lambda e: e[0])
    while s.t < T - 1e-9:
        while events and s.t >= events[0][0] - 1e-9:
            _, m, kw = events.pop(0)
            getattr(s, m)(**kw)
        dt = min(s.stable_dt(), CAP, T - s.t)
        if events:
            dt = min(dt, max(events[0][0] - s.t, 1e-6))
        s.step(dt)
        assert not s.blown_up()
    return s


def _lap_rough(s):
    """mean |eta - mean of its 4 neighbours| over interior cells whose 4 neighbours are wet and open (m)."""
    e, ok = s.eta, (s.h > 0.05) & ~s.wall
    c = ok[1:-1, 1:-1] & ok[2:, 1:-1] & ok[:-2, 1:-1] & ok[1:-1, 2:] & ok[1:-1, :-2]
    r = np.abs(e[1:-1, 1:-1] - 0.25 * (e[2:, 1:-1] + e[:-2, 1:-1] + e[1:-1, 2:] + e[1:-1, :-2]))
    return float(r[c].mean())


_REF = {}


def _ref(key, fn):
    if key not in _REF:
        _REF[key] = fn()
    return _REF[key]


def _pit_depth(name, n):
    sc = scenes.build("pit_dig", n)
    pit = sc.events[0][2]["mask"]
    s = _run(REGISTRY[name].from_scene(sc), sc, 150.0)
    return float(s.h[pit].mean())


DEFECT = pytest.mark.xfail(strict=True, reason="the original keeps its defect; the fix is a new competitor")


@pytest.mark.parametrize("name", [pytest.param("local_inertial", marks=DEFECT), "local_inertial_frcap"])
def test_breach_is_not_early(name):
    """Dam-break breach volume in 20 s within 15% of swe_hll (local_inertial: 1.48x, no momentum advection)."""
    ref = _ref("breach", lambda: bakeoff.breach_volume(REGISTRY["swe_hll"]))
    ratio = bakeoff.breach_volume(REGISTRY[name]) / ref
    assert 0.85 <= ratio <= 1.15, ratio


@pytest.mark.parametrize("name", [pytest.param("local_inertial", marks=DEFECT), "local_inertial_frcap"])
def test_pit_fills_on_time_at_128(name):
    """Pit depth at 150 s of pit_dig, 128 cells, within 25% of swe_hll (local_inertial: 0.72 m vs 0.42)."""
    ref = _ref("pit128", lambda: _pit_depth("swe_hll", 128))
    d = _pit_depth(name, 128)
    assert abs(d / ref - 1) <= 0.25, (d, ref)


def test_frcap07_caps_less_than_frcap():
    """The 0.7 cap sits between the 0.544 cap and no cap on the breach (monotone in the cap)."""
    v = {k: bakeoff.breach_volume(REGISTRY[k]) for k in ("local_inertial_frcap", "local_inertial_frcap07",
                                                           "local_inertial")}
    assert v["local_inertial_frcap"] < v["local_inertial_frcap07"] < v["local_inertial"], v


@pytest.mark.parametrize("name", [pytest.param("pipes", marks=DEFECT), "pipes_theta"])
def test_draining_lake_has_no_checkerboard(name):
    """lake_levelset at 120 s, 128 cells: surface roughness below 2e-4 m (pipes: 1.8e-3 m, a checkerboard)."""
    sc = scenes.build("lake_levelset", 128)
    s = _run(REGISTRY[name].from_scene(sc), sc, 120.0)
    assert _lap_rough(s) < 2e-4, _lap_rough(s)


@pytest.mark.parametrize("name", [pytest.param("pipes", marks=DEFECT), "pipes_theta"])
def test_frictionless_pipe_cannot_push_water_over_a_dry_face(name):
    """manning 0: water at eta 1 beside a dry step with bed 2 and a preloaded outflow must not climb the step."""
    g = Grid(1, 2, 1.0)
    z = np.array([[0.0, 2.0]])
    h = np.array([[1.0, 0.0]])
    s = REGISTRY[name](g, z, h, manning=0.0)
    s.f_r[:] = 50.0
    s.step(0.05)
    assert s.h[0, 1] == 0.0, s.h


def test_ca_game_step_scales_with_dt():
    """A ca_game step of dt 0.25 moves half the volume of a 0.5 s pass (it ignored dt and moved a whole pass)."""
    g = Grid(1, 2, 1.0)
    z = np.zeros((1, 2))
    moved = {}
    for dt in (0.5, 0.25):
        s = REGISTRY["ca_game"](g, z.copy(), np.array([[1.0, 0.9]]))
        s.step(dt)
        moved[dt] = 1.0 - s.h[0, 0]
    assert abs(moved[0.25] / moved[0.5] - 0.5) < 1e-9, moved
