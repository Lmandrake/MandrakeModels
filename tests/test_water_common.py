"""Properties every water solver is held to (parametrised over the registry)."""
import math

import numpy as np
import pytest

from mandrakemodels import scenes
from mandrakemodels.solvers import REGISTRY, for_kind

WATER = sorted(c.name for c in for_kind("water"))
# Models that by construction cannot keep a lake still on an uneven bed (they flow down the BED slope).
NOT_WELL_BALANCED = {"kinematic_wave"}


def run(cls, scene, t_end, dt_cap=1.0):
    s = cls.from_scene(scene)
    while s.t < t_end - 1e-9:
        dt = min(s.stable_dt(), dt_cap, t_end - s.t)
        s.step(dt)
        assert not s.blown_up(), f"{cls.name} blew up at t={s.t}"
    return s


@pytest.mark.parametrize("name", WATER)
def test_lake_at_rest_stays_still(name):
    cls = REGISTRY[name]
    sc = scenes.build("lake_at_rest", 48)
    s = run(cls, sc, 10.0)
    eta = s.eta[s.h > 1e-3]
    drift = np.abs(s.h - sc.init["h"]).max()
    if name in NOT_WELL_BALANCED:
        pytest.xfail(f"{name} is not well-balanced by design (max drift {drift:.3f} m)")
    assert drift < 1e-3, f"lake moved {drift:.2e} m"
    assert eta.max() - eta.min() < 1e-3


@pytest.mark.parametrize("name", WATER)
def test_closed_dam_break_conserves_and_is_symmetric(name):
    cls = REGISTRY[name]
    sc = scenes.build("dam_break", 48)
    s = run(cls, sc, 5.0)
    rel = abs(s.volume_error()) / s.V0
    assert rel < 1e-9, f"volume error {rel:.2e}"
    assert s.min_h_seen >= -1e-12
    asym = np.abs(s.h - s.h[::-1, :]).max()
    assert asym < 1e-8, f"y-mirror asymmetry {asym:.2e} m"
    moved = np.abs(s.h - sc.init["h"]).max()
    if name in NOT_WELL_BALANCED:
        # characterisation: kinematic flow follows the bed slope only, and this bed is flat, so the
        # dam cannot break at all (a documented limit, pinned so a change shows up)
        assert moved < 1e-9, f"{name} moved {moved:.2e} m on a flat bed"
        return
    assert moved > 1e-2, "nothing flowed"


@pytest.mark.parametrize("name", WATER)
def test_forcing_ledger_balances(name):
    """River with source and sink, plus rain/evaporation: volume change equals the booked terms."""
    cls = REGISTRY[name]
    sc = scenes.build("river", 48)
    sc.forcing.rain = 1e-4
    sc.forcing.evaporation = 5e-5
    s = run(cls, sc, 5.0)
    assert abs(s.volume_error()) <= 1e-9 * max(s.ledger.throughput, s.V0)
    assert s.ledger.terms["source"] > 0 and s.ledger.terms["rain"] > 0


@pytest.mark.parametrize("name", WATER)
def test_dig_event_keeps_volume(name):
    cls = REGISTRY[name]
    sc = scenes.build("pit_dig", 48)
    s = cls.from_scene(sc)
    for t, meth, kw in sc.events:
        getattr(s, meth)(**kw)
    s.step(min(s.stable_dt(), 0.5))
    assert abs(s.volume_error()) <= 1e-9 * max(s.ledger.throughput, s.V0)



# Pit dig at 128 (the size where the old scene's sill kept every pit dry). Bars, measured 2026-11-12 on the fixed scene:
# momentum and implicit solvers fill the pit to 0.17-0.77 m mean depth in 150 s; the two cellular automata spread one cell per
# pass, so their front is late (ca_game 0.013 m, wca2d 0.076 m): for them the bar is only that the water arrives.
PIT_SLOW_FRONT = {"ca_game", "wca2d"}
# Too slow to run 150 simulated seconds at 128 inside a test (explicit diffusion: ~1e-4 s steps; CG: ~6 min wall). UNTESTED.
PIT_TOO_SLOW = {"diffusion_wave_explicit", "diffusion_wave_implicit_cg"}


@pytest.mark.parametrize("name", WATER)
def test_pit_dig_pit_actually_floods(name):
    """The dug canal must carry river water into the dug pit (the old scene's sill left every solver's pit dry)."""
    if name in PIT_TOO_SLOW:
        pytest.skip(f"{name}: too slow to run the scene at 128 in a test; its pit fill is untested")
    cls = REGISTRY[name]
    sc = scenes.build("pit_dig", 128)
    pit = sc.events[0][2]["mask"]
    s = cls.from_scene(sc)
    events = sorted(sc.events, key=lambda e: e[0])
    while s.t < sc.duration - 1e-9:
        while events and s.t >= events[0][0] - 1e-9:
            _, meth, kw = events.pop(0)
            getattr(s, meth)(**kw)
        nxt = events[0][0] - s.t if events else 1e9
        s.step(max(min(s.stable_dt(), 5.0, sc.duration - s.t, nxt), 1e-6))
        assert not s.blown_up(), f"{name} blew up at t={s.t}"
    depth, volume = s.h[pit].mean(), s.h[pit].sum() * scenes.DX ** 2
    if name in NOT_WELL_BALANCED:
        pytest.xfail(f"{name}: no backwater, the river runs as a film down its bed and never turns into the canal "
                     f"(pit {volume:.1f} m3)")
    if name in PIT_SLOW_FRONT:
        assert volume > 5.0, f"{name}: pit holds {volume:.2f} m3 after {sc.duration:.0f} s"
    else:
        assert depth > 0.1, f"{name}: pit mean depth {depth:.3f} m after {sc.duration:.0f} s"
