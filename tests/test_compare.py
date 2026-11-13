"""The side-by-side tool: one initial grid for every solver, common sample times, files written."""
import json

import numpy as np
import pytest

from mandrakemodels import compare, scenes


def test_identical_start_and_outputs(tmp_path):
    compare.main(["--scene", "dam_break", "--n", "24", "--frames", "4", "--max-wall", "20",
                  "--solvers", "swe_hll,local_inertial,pipes,ca_game", "--out", str(tmp_path)])
    d = tmp_path / "dam_break_24"
    meta = json.loads((d / "meta.json").read_text())
    assert meta["ref"] == "swe_hll"
    assert set(meta["final"]) == {"swe_hll", "local_inertial", "pipes", "ca_game"}
    assert meta["final"]["swe_hll"]["rms_vs_ref"] == 0.0
    for f in meta["final"].values():
        assert f["status"] == "ok" and abs(f["t"] - meta["samples"][-1]) < 1e-9
    init = np.load(d / "init.npz")
    sc = scenes.build("dam_break", 24)
    assert np.array_equal(init["h"], sc.init["h"]) and np.array_equal(init["z"], sc.init["z"])
    assert (d / "state.gif").exists() and (d / "diff.gif").exists()
    html = (tmp_path / "compare.html").read_text()
    assert "dam_break" in html and "data:image/webp;base64" in html


def test_user_grid_from_npy(tmp_path):
    z = np.zeros((20, 20))
    z[:, 10:] = -1.0
    h = np.where(z < 0, 1.0, 0.0)
    np.save(tmp_path / "bed.npy", z)
    np.save(tmp_path / "depth.npy", h)
    compare.main(["--bed", str(tmp_path / "bed.npy"), "--depth", str(tmp_path / "depth.npy"), "--name", "mine",
                  "--duration", "2", "--frames", "2", "--solvers", "swe_hll,pipes", "--out", str(tmp_path)])
    meta = json.loads((tmp_path / "mine_20" / "meta.json").read_text())
    assert meta["final"]["pipes"]["vol_err"] < 1e-9


def test_external_landlab_conserves():
    pytest.importorskip("landlab")
    from mandrakemodels import external
    sc = scenes.build("dam_break", 16)
    s = external.LandlabOverland.from_scene(sc)
    for _ in range(5):
        s.step(min(s.stable_dt(), 0.5))
    assert abs(s.volume_error()) < 1e-9 * s.V0
    c = int(np.where(sc.init["wall"].any(0))[0].max()) + 1
    assert s.h[:, c:].sum() > sc.init["h"][:, c:].sum() + 0.3   # water crossed the breach


def test_external_anuga_agrees_with_swe_hll():
    """ANUGA on the identical dam_break grid moves about the same water past the wall as swe_hll."""
    pytest.importorskip("anuga")
    import copy
    from mandrakemodels import external, scenes
    from mandrakemodels.solvers import REGISTRY
    sc = scenes.build("dam_break", 24)
    cols = np.where(sc.init["wall"].any(0))[0]
    right = {}
    for key, cls in (("ours", REGISTRY["swe_hll"]), ("anuga", external.AnugaSWE)):
        s = cls.from_scene(copy.deepcopy(sc))
        while s.t < 4.0 - 1e-9:
            s.step(min(s.stable_dt(), 4.0 - s.t))
        assert abs(s.volume_error()) < 1e-6 * s.h.sum()
        right[key] = s.h[:, cols[-1] + 1:].sum()
    assert right["ours"] > 0
    assert abs(right["anuga"] - right["ours"]) < 0.1 * right["ours"]
