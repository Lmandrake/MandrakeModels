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
