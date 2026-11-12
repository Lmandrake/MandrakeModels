"""Stam (1999) "Stable Fluids": semi-Lagrangian incompressible Navier-Stokes on a MAC grid."""
from __future__ import annotations

import numpy as np

from ..base import Solver


def bilinear(f, xi, yi):
    """Sample ``f`` (ny, nx) at fractional index coordinates (xi along x/columns, yi along y/rows),
    clamped to the array. Vectorised; no loops."""
    ny, nx = f.shape
    xi = np.clip(xi, 0.0, nx - 1.0)
    yi = np.clip(yi, 0.0, ny - 1.0)
    i0 = np.minimum(xi.astype(np.intp), nx - 2) if nx > 1 else np.zeros(xi.shape, np.intp)
    j0 = np.minimum(yi.astype(np.intp), ny - 2) if ny > 1 else np.zeros(yi.shape, np.intp)
    fx = xi - i0
    fy = yi - j0
    i1 = np.minimum(i0 + 1, nx - 1)
    j1 = np.minimum(j0 + 1, ny - 1)
    return ((1 - fy) * ((1 - fx) * f[j0, i0] + fx * f[j0, i1])
            + fy * ((1 - fx) * f[j1, i0] + fx * f[j1, i1]))


class StableFluids(Solver):
    """Stam-style stable fluids: incompressible Navier-Stokes, semi-Lagrangian advection, pressure
    projection, on a staggered (MAC) grid in cell units (dx = 1 per cell, dt = 1 per step).

    Approximates: smoke and air flow. Buoyancy ``buoyancy * gain * dye`` pushes dye upward (+y); the
    velocity and the dye are advected by a vectorised bilinear backtrace (RK2 midpoint); optional
    viscosity ``nu`` (implicit, Jacobi); optional vorticity confinement ``confinement`` to re-inject
    swirl the advection smears. Pressure projection uses a matrix-free, masked **conjugate-gradient**
    Poisson solve (warm-started), so divergence after projection is at solver tolerance, with solid
    cells and closed edges as Neumann, and Dirichlet p = 0 at the outflow edge. With ``inflow`` > 0 the
    left edge is a velocity inlet (``inflow * speed_scale``), the right edge an outflow, the top and
    bottom free-slip walls; otherwise the box is closed. Obstacles are solid cells with no flow
    through them (weak no-slip: solid cells read as zero velocity when the backtrace samples them).

    Cannot do: it is strongly numerically dissipative (semi-Lagrangian) so vortex streets at Re ~ 150
    are damped and weak, there is no real Reynolds number (``nu`` only adds to the numerical viscosity),
    no compressibility, no turbulence model, no exact conservation of anything (conserves = False).
    Velocity units are cells per step; ``stable_dt`` is a fixed 1.
    """

    name = "stable_fluids"
    family = "stable fluids (Stam 1999)"
    kind = "flow"
    conserves = False

    def __init__(self, grid, obstacle=None, source=None, buoyancy=1.0, inflow=0.0, nu=0.0,
                 gain=0.04, speed_scale=4.0, confinement=0.0, dye_decay=0.0, cg_iters=200, cg_tol=1e-6):
        super().__init__(grid)
        ny, nx = grid.shape
        self.obstacle = np.zeros(grid.shape, bool) if obstacle is None else np.array(obstacle, bool)
        self.source = np.zeros(grid.shape, bool) if source is None else np.array(source, bool)
        self.buoyancy, self.nu, self.gain = float(buoyancy), float(nu), float(gain)
        self.inflow = float(inflow) * (speed_scale if inflow > 0 else 1.0)
        self.confinement, self.dye_decay = float(confinement), float(dye_decay)
        self.cg_iters, self.cg_tol = int(cg_iters), float(cg_tol)
        self.open_out = self.inflow > 0
        fl = ~self.obstacle
        # faces open to pressure (both neighbours fluid); outer faces closed
        self.ox = fl[:, 1:] & fl[:, :-1]
        self.oy = fl[1:, :] & fl[:-1, :]
        # faces where velocity is held at 0 (touching a solid cell, or outer wall)
        self.ufix = np.ones((ny, nx + 1), bool)
        self.ufix[:, 1:-1] = ~self.ox
        self.vfix = np.ones((ny + 1, nx), bool)
        self.vfix[1:-1, :] = ~self.oy
        self.out_rows = np.zeros(ny, bool)
        if self.open_out:
            self.out_rows = fl[:, -1].copy()
            self.ufix[:, -1] = ~self.out_rows
            self.ufix[:, 0] = True            # inlet value set explicitly
        self.u = np.zeros((ny, nx + 1))
        self.v = np.zeros((ny + 1, nx))
        if self.open_out:
            self.u[:, :] = self.inflow
        self.dye = np.zeros(grid.shape)
        self.p = np.zeros(grid.shape)
        self._fl = fl
        self._bc()
        self.dye[self.source] = 1.0
        self.project()

    @classmethod
    def from_scene(cls, scene, **kw):
        p = dict(scene.params)
        p.update(kw)
        p.pop("steps", None)
        p.pop("Re", None)
        return cls(scene.grid, obstacle=scene.init.get("obstacle"), source=scene.init.get("source"), **p)

    def stable_dt(self) -> float:
        return 1.0

    # ------------------------------------------------------------------ boundary conditions
    def _bc(self):
        self.u[self.ufix] = 0.0
        self.v[self.vfix] = 0.0
        if self.open_out:
            self.u[:, 0] = np.where(self._fl[:, 0], self.inflow, 0.0)
            self.v[0, :] = 0.0
            self.v[-1, :] = 0.0

    # ------------------------------------------------------------------ advection
    def _advect(self, f, dt, xs, ys, ox, oy):
        """Advect array f whose samples sit at x = i*sx + ox, y = j*sy + oy (cell units): RK2 backtrace."""
        ny, nx = f.shape
        J, I = np.mgrid[0:ny, 0:nx]
        x = I + ox
        y = J + oy
        u1, v1 = self._vel_at(x, y)
        xm, ym = x - 0.5 * dt * u1, y - 0.5 * dt * v1
        u2, v2 = self._vel_at(xm, ym)
        x0, y0 = x - dt * u2, y - dt * v2
        return bilinear(f, x0 - ox, y0 - oy)

    def _vel_at(self, x, y):
        """Velocity at physical position (x, y) in cell units; cell centres sit at i + 0.5."""
        u = bilinear(self.u, x, y - 0.5)
        v = bilinear(self.v, x - 0.5, y)
        return u, v

    def _centre_velocity(self):
        return 0.5 * (self.u[:, 1:] + self.u[:, :-1]), 0.5 * (self.v[1:, :] + self.v[:-1, :])

    # ------------------------------------------------------------------ forces
    def _forces(self, dt):
        if self.buoyancy:
            d = self.dye
            fv = 0.5 * (d[1:, :] + d[:-1, :])
            self.v[1:-1, :] += dt * self.buoyancy * self.gain * fv
        if self.confinement:
            uc, vc = self._centre_velocity()
            w = np.gradient(vc, axis=1) - np.gradient(uc, axis=0)
            aw = np.abs(w)
            gx, gy = np.gradient(aw, axis=1), np.gradient(aw, axis=0)
            m = np.hypot(gx, gy) + 1e-9
            nxn, nyn = gx / m, gy / m
            fx = self.confinement * nyn * w
            fy = -self.confinement * nxn * w
            self.u[:, 1:-1] += dt * 0.5 * (fx[:, 1:] + fx[:, :-1])
            self.v[1:-1, :] += dt * 0.5 * (fy[1:, :] + fy[:-1, :])

    def _diffuse(self, f, a, fixed):
        """Implicit diffusion by Jacobi iteration; fixed faces hold their value."""
        f0 = f.copy()
        for _ in range(20):
            p = np.pad(f, 1, mode="edge")
            f = (f0 + a * (p[1:-1, :-2] + p[1:-1, 2:] + p[:-2, 1:-1] + p[2:, 1:-1])) / (1 + 4 * a)
            f[fixed] = 0.0 if not self.open_out else f[fixed]
        return f

    # ------------------------------------------------------------------ projection (masked CG)
    def _A(self, p):
        fx = (p[:, :-1] - p[:, 1:]) * self.ox
        fy = (p[:-1, :] - p[1:, :]) * self.oy
        out = np.zeros_like(p)
        out[:, :-1] += fx
        out[:, 1:] -= fx
        out[:-1, :] += fy
        out[1:, :] -= fy
        if self.open_out:
            out[:, -1] += p[:, -1] * self.out_rows       # ghost cell p = 0 beyond the outflow edge
        return np.where(self._fl, out, 0.0)

    def divergence(self):
        d = (self.u[:, 1:] - self.u[:, :-1]) + (self.v[1:, :] - self.v[:-1, :])
        return np.where(self._fl, d, 0.0)

    def project(self):
        div = self.divergence()
        self.stats["div_before"] = float(np.abs(div).max())
        b = -div          # div_new = div + A p  =>  A p = -div
        if not self.open_out:
            b[self._fl] -= b[self._fl].mean()
        x = np.where(self._fl, self.p, 0.0)
        r = b - self._A(x)
        d = r.copy()
        rs = float((r * r).sum())
        bn = max(float(np.sqrt((b * b).sum())), 1e-30)
        it = 0
        while it < self.cg_iters and np.sqrt(rs) > self.cg_tol * bn:
            Ad = self._A(d)
            dAd = float((d * Ad).sum())
            if dAd <= 0:
                break
            a = rs / dAd
            x += a * d
            r -= a * Ad
            rs_new = float((r * r).sum())
            d = r + (rs_new / rs) * d
            rs = rs_new
            it += 1
        self.p = x
        # u -= grad p on open faces (and across the outflow edge)
        self.u[:, 1:-1] -= (x[:, 1:] - x[:, :-1]) * self.ox
        self.v[1:-1, :] -= (x[1:, :] - x[:-1, :]) * self.oy
        if self.open_out:
            self.u[:, -1] -= (0.0 - x[:, -1]) * self.out_rows
        self._bc()
        self.stats["iterations"] = it
        div = self.divergence()
        self.stats["max_divergence"] = float(np.abs(div).max())

    # ------------------------------------------------------------------ step
    def _step(self, dt):
        self._forces(dt)
        un = self._advect(self.u, dt, 1, 1, 0.0, 0.5)
        vn = self._advect(self.v, dt, 1, 1, 0.5, 0.0)
        self.u, self.v = un, vn
        self._bc()
        if self.nu > 0:
            a = self.nu * dt
            self.u = self._diffuse(self.u, a, self.ufix)
            self.v = self._diffuse(self.v, a, self.vfix)
            self._bc()
        self.project()
        # dye: advect by the projected field, then inject
        J, I = np.mgrid[0:self.grid.ny, 0:self.grid.nx]
        x, y = I + 0.5, J + 0.5
        u1, v1 = self._vel_at(x, y)
        u2, v2 = self._vel_at(x - 0.5 * dt * u1, y - 0.5 * dt * v1)
        self.dye = bilinear(self.dye, x - dt * u2 - 0.5, y - dt * v2 - 0.5)
        if self.dye_decay:
            self.dye *= 1 - self.dye_decay * dt
        self.dye[self.source] = 1.0
        self.dye[self.obstacle] = 0.0
        uc, vc = self._centre_velocity()
        sp2 = uc * uc + vc * vc
        self.stats["energy"] = float(0.5 * sp2.sum())
        self.stats["total"] = float(self.dye.sum())
        self.stats["max_speed"] = float(np.sqrt(sp2.max()))

    def blown_up(self) -> bool:
        return (not (np.isfinite(self.u).all() and np.isfinite(self.v).all() and np.isfinite(self.dye).all())
                or float(np.abs(self.u).max()) > 1e3 or float(self.dye.max()) > 1e3)

    def snapshot(self) -> dict:
        uc, vc = self._centre_velocity()
        vort = np.gradient(vc, axis=1) - np.gradient(uc, axis=0)
        vort[self.obstacle] = 0.0
        return {"dye": self.dye.copy(), "speed": np.hypot(uc, vc), "vorticity": vort, "obstacle": self.obstacle}


SOLVERS = (StableFluids,)
