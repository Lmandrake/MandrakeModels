"""Cross-check adapters: open-source solvers run from the same initial grid as ours (optional extras).

Nothing here is a core dependency. Each adapter imports its library lazily; ``available()`` lists the ones
whose library imports, and ``compare --external`` adds them as extra tiles. No code from these libraries is
copied into this repository; they are called. Licences (checked from each package's metadata, 2026-11-12):

==========================  =========================================  ==================  ==============
adapter                     library and method                         licence             scene kind
==========================  =========================================  ==================  ==============
ext_landlab_overland        Landlab ``OverlandFlow`` (de Almeida,       MIT                 water
                            Bates, Freer & Souvignet 2012)
ext_landlab_kinwave         Landlab ``KinematicWaveRengers``            MIT                 water
                            (kinematic wave, Manning)
ext_landlab_lindiff         Landlab ``LinearDiffusionOverlandFlow-     MIT                 water
                            Router`` (linearised diffusion wave)
ext_phiflow_smoke           PhiFlow (Holl et al.) semi-Lagrangian       MIT                 flow
                            advection + pressure projection
==========================  =========================================  ==================  ==============

How an adapter works (water): the library holds its own grid, padded by one closed ring so every one of our
cells is a core node. Each ``_step`` pushes our bed and depth in, advances the library by ``dt`` (the
library may sub-step internally), and pulls the depth back. Our forcing, ledger and events then run as for
any ``WaterSolver``, so volume error is measured the same way. Walls become closed nodes.
"""
from __future__ import annotations

import numpy as np

from .base import Solver, WaterSolver


def _pad(a, v):
    return np.pad(np.asarray(a, dtype=np.float64), 1, constant_values=v)


class _LandlabWater(WaterSolver):
    family = "external"
    conserves = False       # measured, not assumed
    library = "landlab"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, **_):
        super().__init__(grid, z, h=h, wall=wall, forcing=forcing, manning=manning)
        from landlab import RasterModelGrid
        ny, nx = grid.shape
        self.lg = RasterModelGrid((ny + 2, nx + 2), xy_spacing=grid.dx)
        self.lg.add_field("topographic__elevation", _pad(self.z, 0.0).ravel(), at="node")
        self.lg.add_field("surface_water__depth", _pad(self.h, 0.0).ravel(), at="node")
        closed = np.pad(self.wall, 1, constant_values=True).ravel()
        self.lg.status_at_node[closed] = self.lg.BC_NODE_IS_CLOSED
        self._make()

    def _make(self):  # pragma: no cover
        raise NotImplementedError

    def _push(self):
        self.lg.at_node["topographic__elevation"][:] = _pad(self.z, 0.0).ravel()
        self.lg.at_node["surface_water__depth"][:] = _pad(self.h, 0.0).ravel()

    def _pull(self):
        ny, nx = self.grid.shape
        h = self.lg.at_node["surface_water__depth"].reshape(ny + 2, nx + 2)[1:-1, 1:-1]
        self.h[:] = np.where(self.wall, 0.0, h)

    def _step(self, dt):
        self._push()
        self._advance(dt)
        self._pull()


class LandlabOverland(_LandlabWater):
    """Landlab ``OverlandFlow``: de Almeida, Bates, Freer & Souvignet (2012) local-inertial scheme with the
    theta-weighted face discharge. Same family as our ``local_inertial``; the cross-check for it, ``pipes``
    and ``swe_hll``. Library: Landlab (Hobley et al. 2017; Barnhart et al. 2020), MIT."""
    name = "ext_landlab_overland"

    def _make(self):
        from landlab.components import OverlandFlow
        self.of = OverlandFlow(self.lg, mannings_n=self.manning, steep_slopes=True, h_init=1e-5)

    def stable_dt(self):
        self._push()
        try:
            return float(min(self.of.calc_time_step(), 5.0))
        except Exception as e:  # Landlab raises NoWaterError on a dry grid: nothing limits dt there
            if type(e).__name__ != "NoWaterError":
                raise
            return 5.0

    def _advance(self, dt):
        self.of.overland_flow(dt=dt)


class LandlabKinwave(_LandlabWater):
    """Landlab ``KinematicWaveRengers``: kinematic wave with Manning friction, explicit, sub-stepped under a
    Courant limit (Rengers et al. 2016). Cross-check for our ``kinematic_wave``. Library: Landlab, MIT."""
    name = "ext_landlab_kinwave"

    def _make(self):
        from landlab.components import KinematicWaveRengers
        self.kw = KinematicWaveRengers(self.lg, mannings_n=self.manning, dt_max=1.0)

    def stable_dt(self):
        return 1.0

    def _advance(self, dt):
        # np.float64: the library calls .clip() on min(its dt, ours), which fails on a python float
        self.kw.run_one_step(np.float64(dt), rainfall_intensity=0.0)


class LandlabLinDiff(_LandlabWater):
    """Landlab ``LinearDiffusionOverlandFlowRouter``: linearised diffusion-wave overland flow, explicit with
    CFL sub-stepping. Cross-check for our diffusion-wave solvers. Library: Landlab, MIT."""
    name = "ext_landlab_lindiff"

    def _make(self):
        from landlab.components import LinearDiffusionOverlandFlowRouter
        self.ld = LinearDiffusionOverlandFlowRouter(self.lg, roughness=self.manning, rain_rate=0.0)

    def stable_dt(self):
        self._push()
        return float(min(self.ld._cfl_time_step(), 5.0))   # one library iteration per step

    def _advance(self, dt):
        self.ld.run_one_step(dt)


class AnugaSWE(WaterSolver):
    """ANUGA: finite-volume shallow-water equations on unstructured triangles, central-upwind fluxes
    (Kurganov, Noelle & Petrova 2001) with wetting and drying (Roberts et al. 2015, ANUGA user manual;
    Davies & Roberts 2015 for the 'DE' algorithms). Our grid becomes a rectangular-cross mesh, four
    triangles per cell; walls become a bed 10 m above the highest surface; reflective outer boundary.
    ANUGA sub-steps under its own CFL. Library: ANUGA (Geoscience Australia / ANU), Apache-2.0; needs an
    MPI library at import (the ``mpich`` wheel serves). The cross-check for ``swe_hll``."""
    name = "ext_anuga"
    family = "external"
    conserves = False
    library = "anuga"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, **_):
        super().__init__(grid, z, h=h, wall=wall, forcing=forcing, manning=manning)
        import anuga
        ny, nx = grid.shape
        dx = grid.dx
        d = anuga.rectangular_cross_domain(nx, ny, len1=nx * dx, len2=ny * dx)
        d.set_flow_algorithm("DE0")
        d.set_quantity("friction", manning)
        d.set_name("ext_anuga")
        d.set_store(False)
        d.set_quantities_to_be_monitored(None) if hasattr(d, "set_quantities_to_be_monitored") else None
        c = d.centroid_coordinates
        self.cj = np.clip((c[:, 1] / dx).astype(int), 0, ny - 1)
        self.ci = np.clip((c[:, 0] / dx).astype(int), 0, nx - 1)
        self.cell = self.cj * nx + self.ci
        self.count = np.bincount(self.cell, minlength=ny * nx).astype(float)
        br = anuga.Reflective_boundary(d)
        d.set_boundary({tag: br for tag in d.get_boundary_tags()})
        self.d = d
        self._pz = self._ph = None
        self._push()

    def _bed(self):
        top = float(np.max(self.z + self.h)) + 10.0
        return np.where(self.wall, top, self.z)

    def _push(self):
        if self._pz is None or not np.array_equal(self._pz, self.z):
            self.d.set_quantity("elevation", self._bed()[self.cj, self.ci], location="centroids")
            self._pz = self.z.copy()
            self._ph = None
        if self._ph is None or not np.array_equal(self._ph, self.h):
            el = self.d.quantities["elevation"].centroid_values
            self.d.set_quantity("stage", el + self.h[self.cj, self.ci], location="centroids")

    def _pull(self):
        q = self.d.quantities
        ht = np.maximum(q["stage"].centroid_values - q["elevation"].centroid_values, 0.0)
        h = (np.bincount(self.cell, ht, minlength=self.count.size) / self.count).reshape(self.grid.shape)
        self.h[:] = np.where(self.wall, 0.0, h)
        self._ph = self.h.copy()

    def stable_dt(self):
        return 1.0

    def _step(self, dt):
        self._push()
        for _ in self.d.evolve(yieldstep=dt, finaltime=self.d.get_time() + dt):
            pass
        self._pull()


class PhiFlowSmoke(Solver):
    """PhiFlow buoyant smoke: semi-Lagrangian advection of dye and staggered velocity, buoyancy
    ``buoyancy * gain * dye`` upward, pressure projection by conjugate gradient with box obstacles
    (one unit box per obstacle cell). The same recipe as our ``stable_fluids`` (Stam 1999), written with
    PhiFlow's operators (Holl, Koltun & Thuerey 2020, "Learning to control PDEs with differentiable
    physics"). Closed box, dt in cell units. Library: PhiFlow / phiml, MIT. Only the smoke_plume scene."""
    name = "ext_phiflow_smoke"
    kind = "flow"
    family = "external"
    library = "phi"

    @classmethod
    def applies(cls, sc):
        return sc.name == "smoke_plume"

    def __init__(self, grid, obstacle=None, source=None, buoyancy=1.0, gain=0.04, **_):
        super().__init__(grid)
        from phi import flow as pf
        self.pf = pf
        ny, nx = grid.shape
        self.obstacle = np.zeros(grid.shape, bool) if obstacle is None else np.array(obstacle, bool)
        src = np.zeros(grid.shape, bool) if source is None else np.array(source, bool)
        self.force = float(buoyancy) * float(gain)
        b = pf.Box(x=nx, y=ny)
        jj, ii = np.nonzero(self.obstacle)
        self.obs = ()
        if len(ii):
            lo = pf.math.tensor(np.stack([ii, jj], -1).astype(float), pf.instance("o"), pf.channel(vector="x,y"))
            self.obs = (pf.Obstacle(pf.Box(lo, lo + 1)),)
        self.src = pf.math.tensor(src.T.astype(float), pf.spatial("x,y"))
        self.keep = pf.math.tensor((~self.obstacle).T.astype(float), pf.spatial("x,y"))
        self.dye_g = pf.CenteredGrid(0, pf.ZERO_GRADIENT, b, x=nx, y=ny)
        self.vel = pf.StaggeredGrid(0, 0, b, x=nx, y=ny)
        self.p = None
        self._dye_np()

    @classmethod
    def from_scene(cls, scene, **kw):
        p = {k: v for k, v in scene.params.items() if k in ("buoyancy", "gain")}
        p.update(kw)
        return cls(scene.grid, obstacle=scene.init.get("obstacle"), source=scene.init.get("source"), **p)

    def stable_dt(self):
        return 1.0

    def _dye_np(self):
        self.dye = np.asarray(self.dye_g.values.numpy("y,x"), dtype=np.float64)

    def _step(self, dt):
        pf, m = self.pf, self.pf.math
        d = pf.advect.semi_lagrangian(self.dye_g, self.vel, dt)
        d = d.with_values(m.where(self.src > 0.5, 1.0, d.values) * self.keep)
        buoy = pf.resample(d * (0, self.force * dt), to=self.vel)
        v = pf.advect.semi_lagrangian(self.vel, self.vel, dt) + buoy
        solve = pf.Solve("CG", 1e-3, 1e-4, x0=self.p, max_iterations=2000, suppress=(m.NotConverged, m.Diverged))
        self.vel, self.p = pf.fluid.make_incompressible(v, self.obs, solve)
        self.dye_g = d
        self._dye_np()

    def blown_up(self):
        return not np.isfinite(self.dye).all() or float(np.abs(self.dye).max()) > 1e3

    def snapshot(self):
        c = self.vel.at_centers()
        uv = np.asarray(c.values.numpy("y,x,vector"), dtype=np.float64)
        uc, vc = uv[..., 0], uv[..., 1]
        vort = np.gradient(vc, axis=1) - np.gradient(uc, axis=0)
        vort[self.obstacle] = 0.0
        return {"dye": self.dye.copy(), "speed": np.hypot(uc, vc), "vorticity": vort, "obstacle": self.obstacle}


SOLVERS = (LandlabOverland, LandlabKinwave, LandlabLinDiff, AnugaSWE, PhiFlowSmoke)


def available():
    """(name -> class for every adapter whose library imports, name -> reason for the rest)."""
    out, notes = {}, {}
    for cls in SOLVERS:
        try:
            __import__(cls.library)
            out[cls.name] = cls
        except Exception as e:  # noqa: BLE001
            notes[cls.name] = f"{type(e).__name__}: {e}"
    return out, notes
