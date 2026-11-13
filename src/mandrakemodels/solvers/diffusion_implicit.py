"""Implicit diffusion-wave solvers (backward Euler on the surface, conservative corrector).

Sources: implicit diffusion-wave overland flow (e.g. Bates & De Roo 2000 discuss the explicit limit;
implicit Picard/Newton forms as in Panday & Huyakorn 2004 / ParFlow's overland module, Kollet & Maxwell
2006). Solvers: scipy conjugate gradient (Hestenes & Stiefel 1952) and Peaceman-Rachford ADI (1955).
Own implementation.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ..base import WaterSolver
from .diffusion_explicit import face_conductance


def _thomas(a, b, c, d):
    """Batched tridiagonal solve along axis 1: rows are independent systems
    a[:,i] x[i-1] + b[:,i] x[i] + c[:,i] x[i+1] = d[:,i] (a[:,0], c[:,-1] unused)."""
    n = b.shape[1]
    cp = np.empty_like(b)
    dp = np.empty_like(b)
    cp[:, 0] = c[:, 0] / b[:, 0]
    dp[:, 0] = d[:, 0] / b[:, 0]
    for i in range(1, n):
        m = b[:, i] - a[:, i] * cp[:, i - 1]
        cp[:, i] = c[:, i] / m
        dp[:, i] = (d[:, i] - a[:, i] * dp[:, i - 1]) / m
    x = np.empty_like(b)
    x[:, -1] = dp[:, -1]
    for i in range(n - 2, -1, -1):
        x[:, i] = dp[:, i] - cp[:, i] * x[:, i + 1]
    return x


class _ImplicitDiffusion(WaterSolver):
    """Shared parts: face conductance lagged at the old state and the conservative corrector."""

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, s_min=1e-4,
                 max_dt=5.0, picard=2, max_slices=400):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.s_min, self.max_dt, self.picard = s_min, max_dt, picard
        self.max_slices = max_slices

    def stable_dt(self) -> float:
        return self.max_dt

    def _solve(self, Kx, Ky, eta0, dt):  # pragma: no cover
        """Return (eta_new, Fx, Fy) with F the signed face volumes (m^3, +x/+y)."""
        raise NotImplementedError

    def _step(self, dt):
        eta0 = self.eta
        eta = eta0
        for _ in range(max(1, self.picard)):
            Kx, Ky = face_conductance(self, eta, self.s_min)
            eta, Fx, Fy = self._solve(Kx, Ky, eta0, dt)
        # conservative corrector: water moves as face volumes, h is never assigned from eta
        self._corrector(Fx, Fy, dt)

    def _corrector(self, Fx, Fy, dt):
        """Apply the solved face volumes through ``transfer``. A cell's through-flow can exceed its own
        water even when its net change is fine (it is refilled by other faces in the same solve), and
        ``transfer``'s donor limiter would then distort the balance. So the volumes are applied in m equal
        slices, m chosen so no cell gives more than it holds per slice (capped at ``max_slices``)."""
        Fx = np.where(self.open_x, Fx, 0.0)
        Fy = np.where(self.open_y, Fy, 0.0)
        out = np.zeros(self.grid.shape)
        out[:, :-1] += np.maximum(Fx, 0)
        out[:, 1:] += np.maximum(-Fx, 0)
        out[:-1, :] += np.maximum(Fy, 0)
        out[1:, :] += np.maximum(-Fy, 0)
        V = self.h * self.grid.area
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.where(out > 1e-9 * self.grid.area, out / np.maximum(V, 1e-300), 0.0)
        m = int(min(max(1, math.ceil(float(ratio.max()))), self.max_slices))
        self.stats["slices"] = m
        tx = np.zeros_like(Fx)
        ty = np.zeros_like(Fy)
        for _ in range(m):
            self.transfer(Fx / m, Fy / m)
            tx += self._last_F[0]
            ty += self._last_F[1]
        self._last_F = (tx, ty)
        self.record_flux(tx, ty, dt)


class DiffusionWaveImplicitCG(_ImplicitDiffusion):
    """Implicit diffusion wave solved with conjugate gradients.

    Approximates: the same zero-inertia (diffusion) wave as ``diffusion_wave_explicit`` with Manning
    friction and surface slope, advanced by backward Euler: (A/dt + L(K)) eta_new = (A/dt) eta_old,
    with the face conductance K lagged at the old state (``picard`` fixed-point passes, default 2),
    closed or dry faces K = 0. The sparse symmetric positive-definite system is built vectorised
    (COO) and solved by scipy CG with a Jacobi preconditioner (iterations in ``stats['iters']``; it
    solves for the change in eta, which is the same system, so convergence tolerance is relative to
    the actual update). A conservative corrector then moves F = K (eta_i - eta_j) dt through
    ``transfer`` (never assigning h from eta), so volume is exact and depth >= 0.

    Cannot do: waves, inertia, bores; accuracy at large dt. ``max_dt`` (default 5 s) is a large
    stable step, NOT an accurate one: the timing of transients is smeared (backward Euler damps) and
    the lagged conductance and donor limiter make wet/dry fronts crude.
    """

    name = "diffusion_wave_implicit_cg"
    family = "diffusion wave (implicit)"

    def __init__(self, *a, rtol=1e-13, **kw):
        super().__init__(*a, **kw)
        self.rtol = rtol

    def _solve(self, Kx, Ky, eta0, dt):
        ny, nx = self.grid.shape
        A = self.grid.area
        idx = np.arange(ny * nx).reshape(ny, nx)
        ia, ib = idx[:, :-1].ravel(), idx[:, 1:].ravel()
        ja, jb = idx[:-1, :].ravel(), idx[1:, :].ravel()
        kx, ky = Kx.ravel(), Ky.ravel()
        rows = np.concatenate([ia, ib, ja, jb])
        cols = np.concatenate([ib, ia, jb, ja])
        vals = -np.concatenate([kx, kx, ky, ky])
        diag = np.full(ny * nx, A / dt)
        diag += np.bincount(np.concatenate([ia, ib, ja, jb]),
                            weights=np.concatenate([kx, kx, ky, ky]), minlength=ny * nx)
        M = sp.coo_matrix((np.concatenate([vals, diag]),
                           (np.concatenate([rows, idx.ravel()]), np.concatenate([cols, idx.ravel()]))),
                          shape=(ny * nx, ny * nx)).tocsr()
        e0 = eta0.ravel()
        # residual of the old state: -L eta_old = flux divergence; solve for the change d
        rhs = np.zeros(ny * nx)
        fx = kx * (e0[ia] - e0[ib])
        fy = ky * (e0[ja] - e0[jb])
        rhs -= np.bincount(ia, weights=fx, minlength=ny * nx)
        rhs += np.bincount(ib, weights=fx, minlength=ny * nx)
        rhs -= np.bincount(ja, weights=fy, minlength=ny * nx)
        rhs += np.bincount(jb, weights=fy, minlength=ny * nx)
        # (rhs is -L eta_old; units m^3/s)
        count = [0]

        def cb(_):
            count[0] += 1
        Minv = sp.diags(1.0 / diag)
        if np.abs(rhs).max() == 0:
            d = np.zeros(ny * nx)
        else:
            d, info = spla.cg(M, rhs, x0=np.zeros(ny * nx), rtol=self.rtol, atol=0.0, M=Minv,
                              maxiter=2000, callback=cb)
        self.stats["iters"] = count[0]
        eta_new = (e0 + d).reshape(ny, nx)
        Fx = Kx * (eta_new[:, :-1] - eta_new[:, 1:]) * dt
        Fy = Ky * (eta_new[:-1, :] - eta_new[1:, :]) * dt
        return eta_new, Fx, Fy


class DiffusionWaveImplicitADI(_ImplicitDiffusion):
    """Implicit diffusion wave solved by directional splitting (ADI-style, Lie splitting).

    Approximates: the same zero-inertia wave as the CG variant, but each step does an implicit x sweep
    then an implicit y sweep, each a batched tridiagonal (Thomas) solve in numpy over all rows at once:
    (A/dt + Lx(K)) eta* = (A/dt) eta_old, then (A/dt + Ly(K)) eta_new = (A/dt) eta*, K lagged at the old
    state (``picard`` passes). Each sweep is a conservative backward-Euler step in one direction. The
    corrector moves the x-sweep face volumes and the y-sweep face volumes together through
    ``transfer`` (h is never assigned from eta), so volume is exact and depth >= 0.

    Cannot do: waves, inertia, bores. The splitting adds an O(dt) error (first order, and it
    treats the two directions in sequence, so diagonal flow is slightly anisotropic at large dt).
    ``max_dt`` (default 5 s) is a large stable step, NOT an accurate one.
    """

    name = "diffusion_wave_implicit_adi"
    family = "diffusion wave (implicit)"

    @staticmethod
    def _sweep(K, e, c, axis):
        """One backward-Euler tridiagonal sweep along ``axis`` with face conductances K (c = A/dt)."""
        if axis == 0:
            e, K = e.T, K.T
        n_, m_ = e.shape
        Kp = np.zeros((n_, m_ + 1))
        Kp[:, 1:-1] = K
        lo, hi = Kp[:, :-1], Kp[:, 1:]            # conductance to the previous / next cell
        b = c + lo + hi
        a = -lo
        cc = -hi
        x = _thomas(a, b, cc, c * e)
        return x.T if axis == 0 else x

    def _solve(self, Kx, Ky, eta0, dt):
        c = self.grid.area / dt
        e1 = self._sweep(Kx, eta0, c, 1)
        Fx = Kx * (e1[:, :-1] - e1[:, 1:]) * dt
        e2 = self._sweep(Ky, e1, c, 0)
        Fy = Ky * (e2[:-1, :] - e2[1:, :]) * dt
        return e2, Fx, Fy


SOLVERS = (DiffusionWaveImplicitCG, DiffusionWaveImplicitADI)
