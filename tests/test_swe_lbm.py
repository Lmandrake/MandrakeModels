"""SweLBM (Zhou shallow-water lattice Boltzmann): wet-bed dam break vs Stoker, and mass conservation."""
import math

import numpy as np
from scipy.optimize import brentq

from mandrakemodels import scenes
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


def test_swe_lbm_dam_break_plateau_matches_stoker():
    hl, hr, T = 2.0, 0.3, 3.0
    hm, um, sp = stoker(hl, hr)
    s = strip_dam(REGISTRY["swe_lbm"], t_end=T)
    x = (np.arange(400) + 0.5) * 0.25 - 50.0
    tail = (um - math.sqrt(G * hm)) * T
    shock = sp * T
    mid = (x > tail + 1.0) & (x < shock - 1.5)
    assert mid.sum() > 5
    assert abs(s.h[1][mid].mean() - hm) / hm < 0.10


def test_swe_lbm_conserves_mass_on_river():
    sc = scenes.build("river", 48)
    s = REGISTRY["swe_lbm"].from_scene(sc)
    while s.t < 60.0 - 1e-9:
        s.step(min(s.stable_dt(), 1.0, 60.0 - s.t))
        assert not s.blown_up()
    assert abs(s.volume_error()) <= 1e-9 * max(s.ledger.throughput, s.V0)
