"""Lattice Boltzmann D2Q9 (BGK / TRT) incompressible flow with bounce-back obstacles and passive dye.

Sources: D2Q9 BGK of Qian, d'Humieres & Lallemand (1992), Europhys. Lett. 17; Guo, Zheng & Shi (2002)
forcing; two-relaxation-time with magic 1/4 family of Ginzburg (2005; 3/16 choice per Ginzburg et al.
2008); textbook: Kruger et al. (2017) "The Lattice Boltzmann Method". Own implementation.
"""
from __future__ import annotations

import numpy as np

from ..base import Solver
from .stable_fluids import bilinear

_E = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]])
_W = np.array([4 / 9] + [1 / 9] * 4 + [1 / 36] * 4)
_OPP = np.array([0, 3, 4, 1, 2, 7, 8, 5, 6])


class LBMD2Q9(Solver):
    """D2Q9 lattice Boltzmann, single-relaxation BGK or two-relaxation (TRT, magic parameter 3/16).

    Approximates: weakly compressible incompressible Navier-Stokes (Mach << 1), with viscosity
    nu = (tau - 1/2) / 3 in lattice units, 2nd-order accurate, full-way bounce-back at solid nodes
    (obstacles and closed walls; the wall sits half a cell inside the solid, mass is conserved
    exactly in a closed box). Optional body force (Guo scheme, BGK only) and buoyancy
    ``buoyancy * gain * dye`` upward, a passive dye advected by the lattice velocity (semi-Lagrangian
    bilinear, vectorised) and injected at the source. With ``inflow`` > 0: equilibrium velocity inlet
    on the left column, zero-gradient outlet on the right, no-slip top and bottom walls. tau is set from
    Re = U D / nu with D the obstacle diameter and is floored at 0.51 (so the effective Re can be
    lower than requested; ``effective_Re`` reports it). BGK is stable at cylinder_wake's Re (checked);
    ``collision="trt"`` is available but blew up at the inlet corner there, so it is not the default.

    Cannot do: compressible or high-Mach flow (keep |u| < ~0.15 lattice units; ``substeps`` runs
    several lattice steps per ``step`` to turn lattice speed into scene speed instead), free surfaces,
    thermal flow (buoyancy is a body force on the dye, not temperature), curved walls to better than
    first order (staircase boundary), accurate pressure at the inlet/outlet.
    """

    name = "lbm_d2q9"
    family = "lattice Boltzmann D2Q9"
    kind = "flow"
    conserves = True

    def __init__(self, grid, obstacle=None, source=None, inflow=0.0, Re=None, nu=None, tau=None,
                 buoyancy=0.0, gain=2e-3, force=(0.0, 0.0), collision=None, substeps=1,
                 periodic_x=False, periodic_y=False, dye_decay=0.0):
        super().__init__(grid)
        ny, nx = grid.shape
        self.obstacle = np.zeros(grid.shape, bool) if obstacle is None else np.array(obstacle, bool)
        self.source = np.zeros(grid.shape, bool) if source is None else np.array(source, bool)
        self.inflow = float(inflow)
        self.open_flow = self.inflow > 0
        self.buoyancy, self.gain = float(buoyancy), float(gain)
        self.force = np.array(force, dtype=np.float64)
        self.substeps = int(substeps)
        self.dye_decay = dye_decay
        diam = 2.0 * np.sqrt(self.obstacle.sum() / np.pi) if self.obstacle.any() else float(ny)
        if tau is None:
            if nu is None:
                nu = self.inflow * diam / Re if (Re and self.inflow > 0) else 0.05
            tau = max(0.51, 3 * nu + 0.5)
        self.tau = float(tau)
        self.nu = (self.tau - 0.5) / 3
        self.effective_Re = self.inflow * diam / self.nu if self.inflow > 0 else None
        self.collision = collision or "bgk"
        if self.collision == "trt":
            self.tau_m = 0.5 + (3 / 16) / (self.tau - 0.5)
        # solid = obstacle plus closed outer walls
        solid = self.obstacle.copy()
        if not periodic_y and not self.open_flow:
            solid[0, :] = solid[-1, :] = True
        if not periodic_x and not self.open_flow:
            solid[:, 0] = solid[:, -1] = True
        if not periodic_y and self.open_flow:
            solid[0, :] = solid[-1, :] = True
        self.solid = solid
        self.fluid = ~solid
        self.f = self._feq(np.ones(grid.shape), np.zeros(grid.shape), np.zeros(grid.shape))
        if self.open_flow:
            self.f = self._feq(np.ones(grid.shape), np.where(self.fluid, self.inflow, 0.0), np.zeros(grid.shape))
        self.dye = np.zeros(grid.shape)
        self.dye[self.source] = 1.0
        self.u = np.zeros(grid.shape)
        self.v = np.zeros(grid.shape)
        self.mass0 = self.mass()
        self._moments()

    @classmethod
    def from_scene(cls, scene, **kw):
        p = dict(scene.params)
        p.update(kw)
        p.pop("steps", None)
        inflow = p.pop("inflow", 0.0)
        buoy = p.pop("buoyancy", 0.0)
        Re = p.pop("Re", None)
        nu = p.pop("nu", 0.0) or None
        if inflow == 0 and buoy and "substeps" not in p:
            p["substeps"] = 4          # smoke: lattice speed is capped by Mach, so run 4 lattice steps per step
        if inflow > 0 and "substeps" not in p:
            p["substeps"] = 3          # wake: at tau 0.51 the street needs ~1000 lattice steps to develop
        return cls(scene.grid, obstacle=scene.init.get("obstacle"), source=scene.init.get("source"),
                   inflow=inflow, Re=Re, nu=nu, buoyancy=buoy, **p)

    def stable_dt(self) -> float:
        return 1.0

    # ------------------------------------------------------------------ helpers
    @staticmethod
    def _feq(rho, u, v):
        eu = _E[:, 0, None, None] * u + _E[:, 1, None, None] * v
        usq = u * u + v * v
        return _W[:, None, None] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * usq)

    def mass(self) -> float:
        return float(self.f.sum())

    def _moments(self):
        rho = self.f.sum(0)
        mx = (self.f * _E[:, 0, None, None]).sum(0)
        my = (self.f * _E[:, 1, None, None]).sum(0)
        Fx, Fy = self._body_force()
        self.rho = rho
        self.u = np.where(self.fluid, (mx + 0.5 * Fx) / rho, 0.0)
        self.v = np.where(self.fluid, (my + 0.5 * Fy) / rho, 0.0)

    def _body_force(self):
        Fx = np.full(self.grid.shape, self.force[0])
        Fy = np.full(self.grid.shape, self.force[1])
        if self.buoyancy:
            Fy = Fy + self.buoyancy * self.gain * self.dye
        return Fx, Fy

    # ------------------------------------------------------------------ one lattice step
    def _lattice_step(self):
        f = self.f
        self._moments()
        rho, u, v = self.rho, self.u, self.v
        feq = self._feq(rho, u, v)
        if self.collision == "trt":
            fo, feo = f[_OPP], feq[_OPP]
            fp, fm = 0.5 * (f + fo), 0.5 * (f - fo)
            ep, em = 0.5 * (feq + feo), 0.5 * (feq - feo)
            fn = f - (fp - ep) / self.tau - (fm - em) / self.tau_m
        else:
            fn = f - (f - feq) / self.tau
            Fx, Fy = self._body_force()
            if np.any(Fx) or np.any(Fy):
                eu = _E[:, 0, None, None] * u + _E[:, 1, None, None] * v
                ef = _E[:, 0, None, None] * Fx + _E[:, 1, None, None] * Fy
                uf = u * Fx + v * Fy
                fn += (1 - 0.5 / self.tau) * _W[:, None, None] * (3 * ef - 3 * uf + 9 * eu * ef)
        fn = np.where(self.fluid, fn, f)
        # stream (periodic via roll; walls are solid nodes, inlet/outlet fixed below)
        out = np.empty_like(fn)
        for k in range(9):
            out[k] = np.roll(fn[k], (_E[k, 1], _E[k, 0]), axis=(0, 1))
        # full-way bounce-back at solid nodes
        out[:, self.solid] = out[_OPP][:, self.solid]
        if self.open_flow:
            fl0 = self.fluid[:, 0]
            ny = self.grid.ny
            eq_in = self._feq(np.ones((ny, 1)), np.where(fl0, self.inflow, 0.0)[:, None], np.zeros((ny, 1)))[:, :, 0]
            out[:, :, 0] = np.where(fl0, eq_in, out[:, :, 0])
            out[:, :, -1] = np.where(self.fluid[:, -1], out[:, :, -2], out[:, :, -1])
        self.f = out

    def _advect_dye(self):
        ny, nx = self.grid.shape
        J, I = np.mgrid[0:ny, 0:nx]
        x0, y0 = I - self.u, J - self.v
        self.dye = bilinear(self.dye, x0, y0)
        if self.dye_decay:
            self.dye *= 1 - self.dye_decay
        self.dye[self.source] = 1.0
        self.dye[self.solid] = 0.0

    def _step(self, dt):
        for _ in range(self.substeps):
            self._lattice_step()
            self._moments()
            self._advect_dye()
        speed = np.hypot(self.u, self.v)
        self.stats["total"] = self.mass()
        self.stats["mass_rel_error"] = (self.mass() - self.mass0) / self.mass0
        self.stats["max_speed"] = float(speed.max())
        self.stats["energy"] = float(0.5 * (self.rho * speed ** 2).sum())
        self.stats["tau"] = self.tau

    def blown_up(self) -> bool:
        return (not np.isfinite(self.f).all()) or float(np.abs(self.f).max()) > 1e3 \
            or float(np.hypot(self.u, self.v).max()) > 0.7

    def snapshot(self) -> dict:
        vort = np.gradient(self.v, axis=1) - np.gradient(self.u, axis=0)
        vort[self.solid] = 0.0
        return {"dye": self.dye.copy(), "speed": np.hypot(self.u, self.v), "vorticity": vort,
                "obstacle": self.obstacle}


SOLVERS = (LBMD2Q9,)
