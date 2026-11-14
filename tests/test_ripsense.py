"""local_inertial_frcap_ripsense: a ripple detector plus surface diffusion only where it fires (2026-11-14).

Item RIMMANDRAKE_TRY_RIPPLE_DETECTOR_TREAT_ONLY_WHERE_IT_RIPPLES. Pinned numbers are the swe_hll referee's, measured
2026-11-14: reservoir roughness at 32.7 s 1.84e-4 m, pit_dig 128 mean pit depth at 150 s 0.418 m.
"""
import numpy as np

from mandrakemodels import bakeoff, scenes
from mandrakemodels.solvers import REGISTRY
from mandrakemodels.solvers.local_inertial_ripsense import LocalInertialFrCapRipSense as RS, ripple_flags

SWE_ROUGH_32 = 1.84e-4
SWE_PIT_128 = 0.418


def _run(s, sc, T):
    events = sorted(sc.events, key=lambda e: e[0])
    while s.t < T - 1e-9:
        while events and s.t >= events[0][0] - 1e-9:
            _, m, kw = events.pop(0)
            getattr(s, m)(**kw)
        dt = min(s.stable_dt(), 5.0, T - s.t)
        if events:
            dt = min(dt, max(events[0][0] - s.t, 1e-6))
        s.step(dt)
        assert not s.blown_up()
    return s


def _reservoir_rough(s, c):
    e = s.h + s.z
    nb = 0.25 * (np.roll(e, 1, 0) + np.roll(e, -1, 0) + np.roll(e, 1, 1) + np.roll(e, -1, 1))
    m = np.zeros_like(e, bool)
    m[2:-2, 2:c - 2] = True
    m &= s.h > 1e-3
    return float(np.abs(e - nb)[m].mean())


def test_zero_kappa_is_the_chosen_engine():
    # the parent's Froude cap is fr_max; a gate of the same name once silently replaced it (breach fell to 0.46x)
    C = type("C", (RS,), {"kappa": 0.0})
    assert C.fr_max == REGISTRY["local_inertial_frcap"].fr_max
    a, b = bakeoff.breach_volume(C), bakeoff.breach_volume(REGISTRY["local_inertial_frcap"])
    assert abs(a - b) <= 1e-9 * b


def test_lake_at_rest_flags_nothing_and_stays_still():
    sc = scenes.build("lake_at_rest", 64)
    s = RS.from_scene(sc)
    e0 = s.eta.copy()
    _run(s, sc, 60.0)
    assert s.stats["treated_cell_steps"] == 0
    assert np.abs(s.eta - e0)[sc.init["h"] > 0].max() < 1e-12


def test_reservoir_ripples_flagged_damped_and_books_exact():
    sc = scenes.build("dam_break", 128)
    c = int(np.where(sc.init["wall"].any(0))[0].max())
    raw = _run(REGISTRY["local_inertial_frcap"].from_scene(sc), sc, 32.7)
    flags = ripple_flags(raw)
    res = np.zeros_like(flags)
    res[2:-2, 2:c - 3] = True
    assert flags[res].mean() > 0.25                     # 39% measured
    s = _run(RS.from_scene(sc), sc, 32.7)
    r = _reservoir_rough(s, c)
    assert r <= 1.3 * SWE_ROUGH_32                      # 2.06e-4 measured (every=10)
    assert r < 0.8 * _reservoir_rough(raw, c)           # frcap 3.05e-4
    assert abs(s.volume_error()) / s.V0 < 1e-12


def test_pit_not_flooded_by_the_damping():
    # damping everywhere (the same diffusion, no detector) filled the 128 pit 20% deeper than swe_hll
    sc = scenes.build("pit_dig", 128)
    pit = sc.events[0][2]["mask"]
    s = _run(RS.from_scene(sc), sc, 150.0)
    assert abs(s.h[pit].mean() / SWE_PIT_128 - 1) <= 0.10   # 1.06 measured
