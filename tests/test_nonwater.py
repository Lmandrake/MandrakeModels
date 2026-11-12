"""Tests for the wave, heat, stable-fluids and LBM solvers (grids <= 64)."""
import math

import numpy as np
import pytest

from mandrakemodels import scenes
from mandrakemodels.base import Grid
from mandrakemodels.solvers import REGISTRY
from mandrakemodels.solvers.heat import HeatADI, HeatExplicit
from mandrakemodels.solvers.lbm_d2q9 import LBMD2Q9
from mandrakemodels.solvers.linear_wave import LinearWave
from mandrakemodels.solvers.stable_fluids import StableFluids


def run_to(solver, T, dt):
    nsteps = int(math.ceil(T / dt - 1e-9))
    h = T / nsteps
    for _ in range(nsteps):
        solver.step(h)
    return nsteps


def gaussian_exact_peak(sc, T):
    s0 = sc.params["sigma0"]
    return s0 ** 2 / (s0 ** 2 + 2 * sc.params["D"] * T)


def test_registry_has_all():
    for n in ("linear_wave", "heat_explicit", "heat_adi", "stable_fluids", "lbm_d2q9"):
        assert n in REGISTRY


@pytest.mark.parametrize("cls,mult", [(HeatExplicit, 1.0), (HeatADI, 1.0), (HeatADI, 10.0)])
def test_heat_matches_exact_gaussian(cls, mult):
    sc = scenes.build("gaussian", 64)
    s = cls.from_scene(sc)
    base_dt = 0.9 * sc.grid.dx ** 2 / (4 * sc.params["D"])
    run_to(s, sc.duration, base_dt * mult)
    assert not s.blown_up()
    peak = s.u.max()
    exact = gaussian_exact_peak(sc, sc.duration)
    assert abs(peak - exact) / exact < 0.02, (peak, exact)
    # whole profile along the centre row
    X, Y = sc.grid.coords()
    L = 64 * sc.grid.dx
    sig2 = sc.params["sigma0"] ** 2 + 2 * sc.params["D"] * sc.duration
    ex = exact * np.exp(-((X - L / 2) ** 2 + (Y - L / 2) ** 2) / (2 * sig2))
    assert np.abs(s.u - ex).max() / exact < 0.03


@pytest.mark.parametrize("cls,mult", [(HeatExplicit, 1.0), (HeatADI, 1.0), (HeatADI, 10.0), (HeatADI, 50.0)])
def test_heat_conserves_total_on_dye(cls, mult):
    sc = scenes.build("dye", 64)
    s = cls.from_scene(sc)
    t0 = s.total()
    base_dt = 0.9 * sc.grid.dx ** 2 / (4 * sc.params["D"])
    run_to(s, sc.duration, base_dt * mult)
    assert abs(s.total() - t0) / t0 < 1e-10
    assert not s.blown_up()


def test_adi_stable_at_10x_explicit_dt_explicit_is_not():
    sc = scenes.build("dye", 48)
    dt = 10 * 0.9 * sc.grid.dx ** 2 / (4 * sc.params["D"])
    adi, ex = HeatADI.from_scene(sc), HeatExplicit.from_scene(sc)
    assert adi.stable_dt() == math.inf
    for _ in range(40):
        adi.step(dt)
        ex.step(dt)
    assert not adi.blown_up() and np.abs(adi.u).max() < 2.0
    assert ex.blown_up() or np.abs(ex.u).max() > 10     # the explicit scheme is unstable there


def test_linear_wave_energy_and_symmetry():
    n = 48
    s = LinearWave(Grid(n, n, 1.0), c=10.0, damping=0.0, drops=[(0.0, 0.5, 0.5)])
    dt = s.stable_dt()
    assert dt == pytest.approx(0.9 / (10.0 * math.sqrt(2)))
    s.step(dt)
    e1 = s.stats["energy"]
    for _ in range(300):
        s.step(dt)
    e2 = s.stats["energy"]
    assert e1 > 0 and abs(e2 - e1) / e1 < 0.05, (e1, e2)
    assert np.abs(s.u - s.u[:, ::-1]).max() < 1e-12
    assert np.abs(s.u - s.u[::-1, :]).max() < 1e-12
    assert not s.blown_up()


def test_linear_wave_scene_runs_with_obstacles_and_drops():
    sc = scenes.build("ripples", 48)
    s = LinearWave.from_scene(sc)
    dt = s.stable_dt()
    for _ in range(int(sc.duration / dt)):
        s.step(dt)
    snap = s.snapshot()
    assert not s.blown_up()
    assert (snap["height"][snap["obstacle"]] == 0).all()
    assert s._next_drop == len(sc.init["drops"])


def test_stable_fluids_projection_reduces_divergence():
    sc = scenes.build("smoke_plume", 48)
    s = StableFluids.from_scene(sc)
    for _ in range(30):
        s.step(s.stable_dt())
    assert not s.blown_up()
    before, after = s.stats["div_before"], s.stats["max_divergence"]
    assert before > 1e-3 and after < 1e-3 * before, (before, after)
    # smoke has risen
    Y = np.arange(48)[:, None]
    assert (s.dye * Y).sum() / s.dye.sum() > 0.06 * 48 + 2


def test_stable_fluids_cylinder_wake_runs():
    sc = scenes.build("cylinder_wake", 48)
    s = StableFluids.from_scene(sc)
    for _ in range(120):
        s.step(1.0)
    assert not s.blown_up() and s.stats["max_divergence"] < 1e-5
    snap = s.snapshot()
    assert set(snap) == {"dye", "speed", "vorticity", "obstacle"}


def test_lbm_conserves_mass_in_closed_box():
    sc = scenes.build("smoke_plume", 32)
    s = LBMD2Q9.from_scene(sc)
    m0 = s.mass()
    for _ in range(80):
        s.step(1.0)
    assert abs(s.mass() - m0) / m0 < 1e-10
    assert not s.blown_up()
    assert s.stats["max_speed"] > 1e-3      # the plume is actually moving


def test_lbm_poiseuille_profile():
    ny, nx, g, tau = 24, 8, 2e-5, 0.8
    s = LBMD2Q9(Grid(ny, nx, 1.0), tau=tau, force=(g, 0.0), periodic_x=True)
    nu = (tau - 0.5) / 3
    for _ in range(14000):
        s.step(1.0)
    y = np.arange(1, ny - 1) + 0.0          # fluid rows; walls half a cell outside the end rows
    y0, y1 = 0.5, ny - 1.5
    exact = g / (2 * nu) * (y - y0) * (y1 - y)
    u = s.u[1:-1, nx // 2]
    assert np.abs(u - exact).max() / exact.max() < 0.03, (u.max(), exact.max())
    assert abs(s.mass() - s.mass0) / s.mass0 < 1e-10


def test_lbm_cylinder_wake_stable_and_sheds():
    sc = scenes.build("cylinder_wake", 64)
    s = LBMD2Q9.from_scene(sc)
    ny = sc.grid.ny
    probe = []
    for _ in range(sc.params["steps"]):
        s.step(1.0)
        probe.append(s.v[ny // 2 + 1, int(64 * 0.4)])
    assert not s.blown_up()
    assert s.tau >= 0.51
    assert np.abs(probe[-128:]).max() > 0.02      # transverse velocity behind the cylinder: the street has started
