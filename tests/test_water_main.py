"""Analytic and solver-specific checks for pipes, swe_hll and local_inertial."""
import math

import numpy as np
import pytest
from scipy.optimize import brentq

from mandrakemodels.base import G, Grid
from mandrakemodels.solvers import REGISTRY


def stoker(hl, hr):
    """Middle depth and shock speed of the wet-bed dam break (Stoker 1957)."""
    cl = math.sqrt(G * hl)

    def f(hm):
        um_rare = 2 * (cl - math.sqrt(G * hm))
        um_shock = (hm - hr) * math.sqrt(0.5 * G * (hm + hr) / (hm * hr))
        return um_rare - um_shock
    hm = brentq(f, hr * 1.0001, hl * 0.9999)
    um = 2 * (cl - math.sqrt(G * hm))
    s = hm * um / (hm - hr)
    return hm, um, s


def strip_dam(cls, n=400, dx=0.25, hl=2.0, hr=0.3, t_end=3.0, **kw):
    g = Grid(3, n, dx)
    z = np.zeros(g.shape)
    h = np.where(np.arange(n)[None, :] < n // 2, hl, hr) * np.ones(g.shape)
    s = cls(g, z, h, manning=0.0, **kw)
    while s.t < t_end - 1e-12:
        s.step(min(s.stable_dt(), t_end - s.t))
    return s


def test_swe_hll_matches_stoker_dam_break():
    hl, hr, T = 2.0, 0.3, 3.0
    hm, um, sp = stoker(hl, hr)
    s = strip_dam(REGISTRY["swe_hll"], t_end=T)
    x = (np.arange(400) + 0.5) * 0.25 - 50.0
    prof = s.h[1]
    # plateau between the rarefaction tail and the shock
    tail = (um - math.sqrt(G * hm)) * T
    shock = sp * T
    mid = (x > tail + 1.0) & (x < shock - 1.5)
    assert mid.sum() > 5
    assert abs(prof[mid].mean() - hm) / hm < 0.02
    # shock position: where depth crosses halfway between hm and hr
    xs = x[np.argmax((prof < 0.5 * (hm + hr)) & (x > 0))]
    assert abs(xs - shock) < 1.0, (xs, shock)
    assert s.stats["limited"] == 0


@pytest.mark.parametrize("name", ["local_inertial", "pipes"])
def test_approximate_solvers_miss_the_stoker_bore(name):
    """Characterisation, not a pass mark. Froude ~1 here, outside both models' validity. Measured
    2026-11-12 at t = 3 s (Stoker plateau 0.915 m, shock at x = 12.9 m; swe_hll 0.914 m, 12.9 m):
    local_inertial plateau 1.265 m (+38%), front 9.9 m; pipes 0.848 m (-7%), front 4.9 m.
    The test pins that they stay finite, conservative and lag the true shock, so a change shows up."""
    hm, _, sp = stoker(2.0, 0.3)
    s = strip_dam(REGISTRY[name], t_end=3.0)
    x = (np.arange(400) + 0.5) * 0.25 - 50.0
    prof = s.h[1]
    xs = x[np.argmax((prof < 0.5 * (hm + 0.3)) & (x > 0))]
    assert xs < sp * 3.0 - 1.0
    assert abs(s.volume_error()) / s.V0 < 1e-12


def test_swe_hll_dry_bed_front_advances_and_stays_positive():
    """Ritter dry-bed dam break: front speed 2 sqrt(g h0); first order lags it, but it must advance."""
    s = strip_dam(REGISTRY["swe_hll"], hr=0.0, t_end=2.0)
    x = (np.arange(400) + 0.5) * 0.25 - 50.0
    front = x[(s.h[1] > 1e-3)].max()
    exact = 2 * math.sqrt(G * 2.0) * 2.0
    assert 0.6 * exact < front < 1.05 * exact, (front, exact)
    assert s.min_h_seen >= 0


def test_pipes_wave_speed_is_depth_independent():
    """The documented flaw: pipes propagate a disturbance at a speed that does not grow with depth."""
    def arrival(depth):
        g = Grid(3, 200, 1.414)
        z = np.zeros(g.shape)
        h = np.full(g.shape, depth)
        h[:, :5] += 0.01 * depth
        s = REGISTRY["pipes"](g, z, h, manning=0.0)
        base = h.copy()
        while s.t < 60:
            s.step(s.stable_dt())
            if abs(s.h[1, 100] - base[1, 100]) > 1e-4 * depth:
                return s.t
        return math.inf
    t1, t4 = arrival(0.25), arrival(4.0)
    # shallow-water physics would make the deep case 4x faster; pipes are within 30% of each other
    assert abs(t1 - t4) / t4 < 0.3, (t1, t4)
