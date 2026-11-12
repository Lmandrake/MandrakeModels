"""Family A water solvers: CA, WCA2D, kinematic, diffusion wave (explicit, CG, ADI)."""
import math

import numpy as np
import pytest

from mandrakemodels import scenes
from mandrakemodels.base import FixedLevel, Forcing, Grid, Sink, Source
from mandrakemodels.solvers import REGISTRY

DIFFUSION = ["diffusion_wave_explicit", "diffusion_wave_implicit_cg", "diffusion_wave_implicit_adi"]


@pytest.mark.parametrize("name", DIFFUSION)
def test_steady_wide_channel_matches_manning(name):
    ny, nx, dx, S, n = 6, 48, 1.414, 0.01, 0.03
    q = 0.2                                              # m^2/s per unit width
    hn = (q * n / math.sqrt(S)) ** 0.6
    g = Grid(ny, nx, dx)
    z = -S * (np.arange(nx)[None, :] + 0.5) * dx * np.ones((ny, 1))
    src = np.zeros((ny, nx), bool)
    src[:, 0] = True
    snk = np.zeros((ny, nx), bool)
    snk[:, -1] = True
    f = Forcing(sources=[Source(src, q * ny * dx)], sinks=[Sink(snk)])
    s = REGISTRY[name](g, z, h=np.full((ny, nx), hn), forcing=f, manning=n)
    t_end = 1500.0
    while s.t < t_end:
        s.step(min(s.stable_dt(), 0.5, t_end - s.t))
    assert not s.blown_up()
    im = nx // 2
    qmid = s.qx[:, im].mean()
    hmid = s.h[:, im].mean()
    manning_q = hmid ** (5 / 3) * math.sqrt(S) / n
    assert abs(qmid - manning_q) / manning_q < 0.15, (qmid, manning_q)
    assert abs(qmid - q) / q < 0.05       # steady: through-flow equals the inflow


@pytest.mark.parametrize("name", ["diffusion_wave_implicit_cg", "diffusion_wave_implicit_adi"])
def test_implicit_large_step_bounded_and_conservative(name):
    sc = scenes.build("dam_break", 48)
    s = REGISTRY[name].from_scene(sc)
    assert s.stable_dt() == 5.0
    for _ in range(8):
        s.step(5.0)
    assert not s.blown_up()
    assert abs(s.volume_error()) / s.V0 < 1e-9
    assert s.h.max() <= sc.init["h"].max() + 1e-9 and s.h.min() >= 0
    assert np.abs(s.h - sc.init["h"]).max() > 1e-2


def test_cg_reports_iterations():
    s = REGISTRY["diffusion_wave_implicit_cg"].from_scene(scenes.build("dam_break", 32))
    s.step(5.0)
    assert s.stats["iters"] > 0


@pytest.mark.parametrize("name", ["ca_game", "wca2d"])
def test_closed_two_level_box_levels_out(name):
    n = 12
    g = Grid(n, n, 1.0)
    z = np.zeros((n, n))
    h = np.where(np.arange(n)[None, :] < n // 2, 1.0, 0.2) * np.ones((n, 1))
    kw = {"tau": 1e-5} if name == "wca2d" else {}
    s = REGISTRY[name](g, z, h=h, **kw)
    V0 = s.volume()
    for _ in range(20000):
        s.step(s.stable_dt())
        if s.h.max() - s.h.min() < 5e-4:
            break
    assert s.h.max() - s.h.min() < 1e-3, s.h.max() - s.h.min()
    assert abs(s.volume() - V0) < 1e-9 * V0
    assert abs(s.h.mean() - 0.6) < 1e-12


def test_kinematic_conserves_on_river():
    sc = scenes.build("river", 48)
    s = REGISTRY["kinematic_wave"].from_scene(sc)
    while s.t < 40.0:
        s.step(min(s.stable_dt(), 1.0))
        assert not s.blown_up()
    assert abs(s.volume_error()) <= 1e-9 * max(s.ledger.throughput, s.V0)
    assert s.ledger.terms["source"] > 0


def test_kinematic_flat_bed_no_flow():
    g = Grid(8, 8, 1.0)
    s = REGISTRY["kinematic_wave"](g, np.zeros((8, 8)), h=np.full((8, 8), 0.5))
    s.step(1.0)
    assert np.abs(s.h - 0.5).max() == 0
