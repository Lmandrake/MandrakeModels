"""Ripple dampers for local_inertial_frcap (2026-11-14): exact books, a lake at rest stays at rest, and the chosen
damper settles the dam-break reservoir where the original does not."""
import numpy as np
import pytest

from mandrakemodels import scenes
from mandrakemodels.solvers import REGISTRY

DAMPERS = [n for n in REGISTRY if n.startswith("local_inertial_frcap_") and n != "local_inertial_frcap07"]


def _run(s, T):
    while s.t < T - 1e-9:
        s.step(min(s.stable_dt(), 5.0, T - s.t))
    return s


@pytest.mark.parametrize("name", DAMPERS)
def test_lake_at_rest_and_books(name):
    sc = scenes.build("lake_at_rest", 64)
    s = REGISTRY[name].from_scene(sc)
    e0 = s.eta.copy()
    _run(s, 30.0)
    assert np.abs(s.eta - e0)[sc.init["h"] > 0].max() < 1e-12
    assert abs(s.volume_error()) < 1e-12 * s.V0


def _rough(name, T=32.7):
    sc = scenes.build("dam_break", 128)
    s = _run(REGISTRY[name].from_scene(sc), T)
    c = int(np.where(sc.init["wall"].any(0))[0].max())
    e = s.eta
    nb = 0.25 * (np.roll(e, 1, 0) + np.roll(e, -1, 0) + np.roll(e, 1, 1) + np.roll(e, -1, 1))
    m = np.zeros_like(e, bool)
    m[2:-2, 2:c - 2] = True
    m &= s.h > 1e-3
    return float(np.abs(e - nb)[m].mean()), abs(s.volume_error()) / s.V0


def test_smoothdiff_fr_settles_the_reservoir():
    r_cap, _ = _rough("local_inertial_frcap")
    r_fix, err = _rough("local_inertial_frcap_smoothdiff_fr")
    r_ref, _ = _rough("swe_hll")
    assert r_cap > 1.5 * r_ref            # the defect: frcap's reservoir is rough (1.65x measured)
    assert r_fix < 1.3 * r_ref            # the fix: within 1.3x of the full shallow-water referee (1.16x measured)
    assert err < 1e-12
