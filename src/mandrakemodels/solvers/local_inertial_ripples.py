"""Ripple dampers for ``local_inertial_frcap``: each is a named competitor beside the original (2026-11-14).

Why: ``local_inertial_frcap`` (the author's pick, 2026-11-14) leaves criss-cross stripes in the dam-break reservoir
at 16 s and 33 s: wave fronts reflecting between the walls, which nothing in the scheme damps. The Froude cap never binds
in that slow water. Each class here adds one damper and is otherwise ``local_inertial_frcap`` unchanged. Every damper
acts on the face discharge before the cap and before the shared donor-limited transfer, so volume stays exact.
"""
from __future__ import annotations

import numpy as np

from ..base import G, WaterSolver

WaterSolverTransfer = WaterSolver.transfer
WaterSolverRecord = WaterSolver.record_flux
from .local_inertial import LocalInertialFrCap


def _face_lap(q, wet, axis_along):
    """Free-slip 5-point Laplacian (in units of 1/dx^2) of face discharges ``q`` over wet faces ``wet``.

    A dry or closed neighbour counts as the face itself (no flux of momentum into banks or walls), so the sum of q over
    a channel cross-section is kept: momentum is spread, not removed. ``axis_along`` is unused (both axes are mixed)."""
    out = np.zeros_like(q)
    for ax, sh in ((0, 1), (0, -1), (1, 1), (1, -1)):
        nb = np.roll(q, sh, axis=ax)
        nw = np.roll(wet, sh, axis=ax)
        # np.roll wraps: the wrapped row/column is the array edge, treat it as closed
        if ax == 0:
            edge = (slice(0, 1) if sh == 1 else slice(-1, None), slice(None))
        else:
            edge = (slice(None), slice(0, 1) if sh == 1 else slice(-1, None))
        nw = nw.copy()
        nw[edge] = False
        out += np.where(nw, nb, q) - q
    return np.where(wet, out, 0.0)


class _Damped(LocalInertialFrCap):
    family = "shallow water (local inertial)"

    def _limit(self, nqx, nqy, hx, hy):
        nqx, nqy = self._damp(nqx, nqy, hx, hy)
        return super()._limit(nqx, nqy, hx, hy)

    def _damp(self, nqx, nqy, hx, hy):  # pragma: no cover
        raise NotImplementedError


class LocalInertialFrCapDivDamp(_Damped):
    """local_inertial_frcap plus divergence damping (a Kelvin-Voigt term): each face discharge is pushed against the
    gradient of the surface's rate of change, ``dq/dt -= nu * d(dh/dt)/dx``, with ``nu = kappa * dx * sqrt(g h_f)``.

    Why: this is the gravity-wave damper of atmospheric models. A wave of wavenumber k loses amplitude at nu k^2 / 2, so
    short crossing fronts die and long swell barely feels it; where the surface is not changing (a lake at rest, a river
    at steady state) the term is exactly zero. dh/dt is last step's depth change over last step's dt (depth, not surface,
    so a dig or raise of the bed is not read as a rate).

    Cannot do: it also slows any changing surface, so it acts during the river's spin-up and on every advancing front;
    explicit, and it shares the stability margin with theta 0.8 at alpha 0.6: kappa 0.25 grew the reservoir's
    roughness 60-fold in 90 s (measured 2026-11-14); kappa stays small.
    """

    name = "local_inertial_frcap_divdamp"
    kappa = 0.04

    def __init__(self, *a, kappa=None, **kw):
        super().__init__(*a, **kw)
        if kappa is not None:
            self.kappa = kappa
        self._hprev = None
        self._rate = np.zeros(self.grid.shape)

    def _step(self, dt):
        if self._hprev is not None:
            self._rate = (self.h - self._hprev) / self._dtprev
        self._hprev = self.h.copy()
        self._dtprev = dt
        self._dt = dt
        super()._step(dt)

    def _damp(self, nqx, nqy, hx, hy):
        dx, dt, r = self.grid.dx, self._dt, self._rate
        wet = self.h > 1e-3
        wx = (hx > 1e-4) & wet[:, 1:] & wet[:, :-1]
        wy = (hy > 1e-4) & wet[1:, :] & wet[:-1, :]
        nux = self.kappa * dx * np.sqrt(G * hx)
        nuy = self.kappa * dx * np.sqrt(G * hy)
        nqx = nqx - np.where(wx, nux * dt * (r[:, 1:] - r[:, :-1]) / dx, 0.0)
        nqy = nqy - np.where(wy, nuy * dt * (r[1:, :] - r[:-1, :]) / dx, 0.0)
        return nqx, nqy


class LocalInertialFrCapDivDampSoft(LocalInertialFrCapDivDamp):
    """local_inertial_frcap_divdamp at half the strength (kappa 0.02)."""

    name = "local_inertial_frcap_divdamp_soft"
    kappa = 0.02


class LocalInertialFrCapVisc(_Damped):
    """local_inertial_frcap plus a small eddy viscosity on the discharge: ``dq/dt += nu * lap(q)`` over wet faces, both
    axes, free-slip at banks and walls (a dry neighbour counts as the face itself), ``nu = kappa * dx * sqrt(g h_f)``.

    Why: a Laplacian damps every short wave in q, crossing fronts included. Free-slip is the difference from the
    rejected cross-axis theta mixing, which pulled discharge into dry bank faces and so dragged the river (-43%).

    Cannot do: it smooths every shear, so a river's fast core slows a little toward its edges; long swell loses a
    little too (rate nu k^2).
    """

    name = "local_inertial_frcap_visc"
    kappa = 0.03

    def __init__(self, *a, kappa=None, **kw):
        super().__init__(*a, **kw)
        if kappa is not None:
            self.kappa = kappa

    def _step(self, dt):
        self._dt = dt
        super()._step(dt)

    def _damp(self, nqx, nqy, hx, hy):
        dx, dt = self.grid.dx, self._dt
        wx, wy = hx > 1e-4, hy > 1e-4
        nqx = nqx + self.kappa * np.sqrt(G * hx) * dt / dx * _face_lap(nqx, wx, 1)
        nqy = nqy + self.kappa * np.sqrt(G * hy) * dt / dx * _face_lap(nqy, wy, 0)
        return nqx, nqy


class LocalInertialFrCapHyper(_Damped):
    """local_inertial_frcap plus hyperviscosity on the discharge: ``q -= beta * lap(lap(q))`` per step (lap in units of
    1/dx^2, free-slip as in ``local_inertial_frcap_visc``).

    Why: a biharmonic damps at k^4, so a 4-cell front loses far more than a 50-cell swell (about 150x less at k ratio
    12); the usual cure for grid-scale noise that must leave the long waves alone.

    Cannot do: per step, so the damping per second grows with the step count (a smaller dt damps more); explicit, beta
    must stay under 1/32.
    """

    name = "local_inertial_frcap_hyper"
    beta = 0.0025

    def __init__(self, *a, beta=None, **kw):
        super().__init__(*a, **kw)
        if beta is not None:
            self.beta = beta

    def _damp(self, nqx, nqy, hx, hy):
        wx, wy = hx > 1e-4, hy > 1e-4
        nqx = nqx - self.beta * _face_lap(_face_lap(nqx, wx, 1), wx, 1)
        nqy = nqy - self.beta * _face_lap(_face_lap(nqy, wy, 0), wy, 0)
        return nqx, nqy


class LocalInertialFrCapDeepTheta(LocalInertialFrCap):
    """local_inertial_frcap with de Almeida's theta lowered toward ``theta_lo`` only in deep, slow water.

    Each step, before the scheme's own theta-0.8 mixing, the remembered face discharge is mixed again with its two
    along-axis neighbours by ``(theta - theta_eff) / 2``, where ``theta_eff`` falls from 0.8 to ``theta_lo`` as the face
    Froude number falls below ``fr0`` and its depth rises past ``h0``. Fast or shallow faces keep theta 0.8.

    Why: lower theta damps more; doing it only where the flow is deep and slow aims it at the reservoir and spares the
    river and the breach.

    Cannot do: mixing is along each axis only (the original theta's direction), so a front running diagonally is damped
    less than one running along an axis.
    """

    name = "local_inertial_frcap_deeptheta"
    theta_lo = 0.7
    fr0 = 0.2
    h0 = 0.5

    def _step(self, dt):
        hx, hy = self.face_flow_depth()
        for q, hf, ax in ((self.qx[:, 1:-1], hx, 1), (self.qy[1:-1, :], hy, 0)):
            ok = hf > 1e-4
            fr = np.where(ok, np.abs(q) / np.where(ok, hf * np.sqrt(G * np.maximum(hf, 1e-4)), 1.0), 0.0)
            s = (1 - np.minimum(fr / self.fr0, 1.0)) * np.minimum(hf / self.h0, 1.0) * ok
            w = 0.5 * (self.theta - self.theta_lo) * s
            if ax == 1:
                p = np.pad(q, ((0, 0), (1, 1)))
                nb = p[:, :-2] + p[:, 2:]
            else:
                p = np.pad(q, ((1, 1), (0, 0)))
                nb = p[:-2, :] + p[2:, :]
            q[...] = q + w * (nb - 2 * q)
        super()._step(dt)


class LocalInertialFrCapSponge(_Damped):
    """local_inertial_frcap with wall-reflection damping: a face's discharge decays at ``gamma`` per second, scaled by
    how close the face is to a wall or the map edge along its own axis (1 at the wall, 0 from ``width`` cells out).

    Why: the stripes are reflections; a wall that absorbs part of each normally incident wave stops feeding them.

    Cannot do: real walls reflect, so a pool beside a wall settles faster than it should and swell against a breakwater
    reflects less; a slow flow along an axis into a wall (a canal ending at a wall) is slowed near the end.
    """

    name = "local_inertial_frcap_sponge"
    gamma = 0.5
    width = 4

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        closed = self.wall
        ny, nx = self.grid.shape
        # distance (cells) from each cell to the nearest wall or edge cell along x and along y
        def dist(axis):
            d = np.full((ny, nx), self.width + 1.0)
            n = nx if axis == 1 else ny
            idx = np.arange(n)
            for k in range(n):
                sl = (slice(None), k) if axis == 1 else (k, slice(None))
                lim = np.minimum(idx[k] + 1, n - idx[k])   # distance to the map edge, counted to the edge face
                d[sl] = np.minimum(d[sl], lim)
            for k in range(n):
                sl = (slice(None), k) if axis == 1 else (k, slice(None))
                wk = closed[sl]
                if not wk.any():
                    continue
                for off in range(-self.width, self.width + 1):
                    j = k + off
                    if 0 <= j < n:
                        tl = (slice(None), j) if axis == 1 else (j, slice(None))
                        d[tl] = np.where(wk, np.minimum(d[tl], abs(off)), d[tl])
            return d
        dxc, dyc = dist(1), dist(0)
        fx = np.minimum(dxc[:, 1:], dxc[:, :-1])
        fy = np.minimum(dyc[1:, :], dyc[:-1, :])
        self._wx = np.clip(1 - (fx - 1) / self.width, 0, 1)
        self._wy = np.clip(1 - (fy - 1) / self.width, 0, 1)

    def _step(self, dt):
        self._dt = dt
        super()._step(dt)

    def _damp(self, nqx, nqy, hx, hy):
        return nqx * np.exp(-self.gamma * self._dt * self._wx), nqy * np.exp(-self.gamma * self._dt * self._wy)


class LocalInertialFrCapEtaDiff(_Damped):
    """local_inertial_frcap plus a mass-only surface diffusion: each wet face also carries ``-nu * d(eta)/dx`` of water,
    ``nu = kappa * dx * sqrt(g h_f)``, moved by the same transfer but NOT remembered as discharge.

    Why: this is the diffusion the first-order HLL referee carries in its mass flux (about c dx / 2), the reason its
    reservoir settles. It reads the surface slope, so a lake at rest (flat surface) gets exactly nothing, and it fades
    with face depth (sqrt(h)), so shallow flow barely feels it. Not remembering it keeps it from acting as extra gravity.

    Cannot do: it adds a little downhill conveyance on any sloped surface (river: about nu * slope per metre width, tiny);
    explicit, and it shares the stability margin with theta 0.8 at alpha 0.6, so kappa stays small (0.04).
    """

    name = "local_inertial_frcap_etadiff"
    kappa = 0.04

    def __init__(self, *a, kappa=None, **kw):
        super().__init__(*a, **kw)
        if kappa is not None:
            self.kappa = kappa
        self._fd = None

    def _step(self, dt):
        self._dt = dt
        self._fd = None
        super()._step(dt)

    def _damp(self, nqx, nqy, hx, hy):
        dx = self.grid.dx
        dEx, dEy = self.face_surface_diff()
        wx, wy = hx > 1e-4, hy > 1e-4
        self._fd = (np.where(wx, -self.kappa * np.sqrt(G * hx) * dEx, 0.0),
                    np.where(wy, -self.kappa * np.sqrt(G * hy) * dEy, 0.0))
        return nqx, nqy

    def _limit(self, nqx, nqy, hx, hy):
        nqx, nqy = self._damp(nqx, nqy, hx, hy)
        nqx, nqy = LocalInertialFrCap._limit(self, nqx, nqy, hx, hy)
        self._q_kept = (nqx, nqy)
        return nqx, nqy

    def transfer(self, Fx, Fy):
        if self._fd is not None:
            w = self.grid.dx * self._dt
            Fx, Fy = Fx + self._fd[0] * w, Fy + self._fd[1] * w
        super().transfer(Fx, Fy)

    def record_flux(self, Fx, Fy, dt):
        # remember the discharge the scheme computed (scaled by the donor limiter like the moved volume), not the diffusion
        w = self.grid.dx * dt
        if self._fd is None:
            return super().record_flux(Fx, Fy, dt)
        qx, qy = self._q_kept
        tx, ty = qx * w + self._fd[0] * w, qy * w + self._fd[1] * w
        with np.errstate(divide="ignore", invalid="ignore"):
            sx = np.where(tx != 0, Fx / tx, 1.0)
            sy = np.where(ty != 0, Fy / ty, 1.0)
        self.qx[:, 1:-1] = qx * sx
        self.qy[1:-1, :] = qy * sy


def _normal_121(d, wet_face, axis):
    """1-2-1 filter (sum, not yet /4) of face surface differences ``d`` along the face-normal ``axis``, faces without
    flow depth counted as zero. In 1-D this is exactly the difference of a 1-2-1-filtered surface, so the checkerboard
    along the flow is removed; a lake at rest has zero difference on every face with flow depth, so it gets nothing."""
    d *= wet_face
    f = d * 2
    if axis == 1:
        f[:, 1:] += d[:, :-1]; f[:, :-1] += d[:, 1:]
    else:
        f[1:] += d[:-1]; f[:-1] += d[1:]
    return f


class LocalInertialFrCapSmoothDiff(LocalInertialFrCapEtaDiff):
    """local_inertial_frcap_etadiff, but each face's surface difference is first smoothed 1-2-1 along the flow axis.

    Why: plain surface diffusion is held to kappa ~0.04 because it also damps the grid's checkerboard mode, which uses up
    the stability margin theta 0.8 and alpha 0.6 leave (kappa 0.08 moved the 250 lake at rest 0.49 m, measured
    2026-11-14). The filter removes the checkerboard along the flow, so the diffusion is zero there and peaks on
    4-cell waves, the crossing fronts' scale; a much larger kappa then stays stable. Faces without flow depth count as
    zero, so a lake at rest (zero difference on every face with flow depth) gets exactly nothing. (A full 2-D wet-only
    smoothing of the surface did the same job, measured 2026-11-14, but cost 1.3x frcap's step at 250 in numpy.)

    Cannot do: as etadiff (a little extra downhill conveyance; untapered it filled the dug pit 15% early at 128);
    2-cell noise across the flow is left to theta, as before.
    """

    name = "local_inertial_frcap_smoothdiff"
    kappa = 0.3
    fr0 = None
    every = 1          # apply the damper on every ``every``-th step at ``every`` times the strength (cost / every)

    def _limit(self, nqx, nqy, hx, hy):
        # K: the damper's volume per (m/s of wave speed) per (m of surface difference) this step; c_K = K * sqrt(g h_f)
        # is shared by the damper and the Froude cap, so the strength costs no extra pass
        K = 0.25 * self.kappa * self.every * self.grid.dx * self._dt
        gk = G * K * K
        cx, cy = np.sqrt(gk * hx), np.sqrt(gk * hy)
        mx, my = cx * hx, cy * hy
        mx *= self.fr_max / K; my *= self.fr_max / K     # the frcap limit, fr_max * h_f * sqrt(g h_f)
        ax, ay = np.abs(nqx), np.abs(nqy)
        self.stats["capped"] += int((ax > mx).sum() + (ay > my).sum())
        self._fd = None
        if self.nsteps % self.every == 0:
            self._fd = self._diffusion(hx, hy, cx, ax, ay, cy, K)
        nqx, nqy = np.clip(nqx, -mx, mx), np.clip(nqy, -my, my)
        self._q_kept = (nqx, nqy)
        return nqx, nqy

    def _diffusion(self, hx, hy, cx, ax, ay, cy, K):
        """Volume (m^3, this step) each face carries DOWN the smoothed surface, sign flipped (transfer subtracts it).
        ``cx``, ``cy`` are K * c; ``ax``, ``ay`` are |q| (overwritten).

        Per face: K * c * (smoothed surface difference) * taper, with c * taper = max(c - |q| / (fr0 h), 0), that is
        c * max(1 - Fr / fr0, 0), folded so it costs one division. A face with no flow depth has c = 0."""
        e = self.z + self.h
        out = []
        for c, hf, d, t, axis in ((cx, hx, e[:, 1:] - e[:, :-1], ax, 1), (cy, hy, e[1:, :] - e[:-1, :], ay, 0)):
            f = _normal_121(d, hf > 0, axis)
            if self.fr0 is not None:
                np.divide(t, hf, out=t, where=hf > 0)
                t *= K / self.fr0
                np.subtract(c, t, out=t)
                np.maximum(t, 0.0, out=t)
                f *= t
            else:
                f *= c
            out.append(f)
        return tuple(out)

    def transfer(self, Fx, Fy):
        if self._fd is not None:              # Fx, Fy are fresh arrays from the step: subtract in place (downhill)
            Fx -= self._fd[0]
            Fy -= self._fd[1]
        WaterSolverTransfer(self, Fx, Fy)

    def record_flux(self, Fx, Fy, dt):
        # remember only the scheme's discharge, scaled as the donor limiter scaled the moved volume
        if self._fd is None:
            return WaterSolverRecord(self, Fx, Fy, dt)
        if not self.limited_faces:
            self.qx[:, 1:-1], self.qy[1:-1, :] = self._q_kept
            return
        w = self.grid.dx * dt
        qx, qy = self._q_kept
        tx, ty = qx * w - self._fd[0], qy * w - self._fd[1]
        with np.errstate(divide="ignore", invalid="ignore"):
            self.qx[:, 1:-1] = qx * np.where(tx != 0, Fx / tx, 1.0)
            self.qy[1:-1, :] = qy * np.where(ty != 0, Fy / ty, 1.0)


class LocalInertialFrCapSmoothDiffFr(LocalInertialFrCapSmoothDiff):
    """local_inertial_frcap_smoothdiff that fades out in fast flow: the diffusion on a face is scaled by
    ``max(0, 1 - Fr / fr0)``, Fr the face's remembered |q| over ``h_f sqrt(g h_f)``.

    Why: untapered, the filtered diffusion also adds conveyance down a canal, and the dug pit filled 15% early at 128
    (kappa 0.3, measured 2026-11-14). The crossing fronts live in slow water, canals and breaches run near the 0.544
    cap, so a Froude taper keeps the damper on the first and off the second.

    Cannot do: as smoothdiff; a slow but steep sheet flow still gets the diffusion.
    """

    name = "local_inertial_frcap_smoothdiff_fr"
    kappa = 0.4
    fr0 = 0.2


class LocalInertialFrCapSmoothDiffFr2(LocalInertialFrCapSmoothDiffFr):
    """local_inertial_frcap_smoothdiff_fr applied every second step at twice the strength: the same damping per second
    for half the extra cost."""

    name = "local_inertial_frcap_smoothdiff_fr2"
    every = 2


class LocalInertialFrCapSmoothDiffFr4(LocalInertialFrCapSmoothDiffFr):
    """local_inertial_frcap_smoothdiff_fr applied every fourth step at four times the strength. Same reservoir, breach,
    pit, river, lake and swell numbers as every step (measured 2026-11-14), at a quarter of the extra cost."""

    name = "local_inertial_frcap_smoothdiff_fr4"
    every = 4


SOLVERS = (LocalInertialFrCapSmoothDiffFr4, LocalInertialFrCapSmoothDiffFr2, LocalInertialFrCapSmoothDiffFr, LocalInertialFrCapSmoothDiff, LocalInertialFrCapEtaDiff, LocalInertialFrCapDivDamp, LocalInertialFrCapDivDampSoft, LocalInertialFrCapVisc, LocalInertialFrCapHyper,
           LocalInertialFrCapDeepTheta, LocalInertialFrCapSponge)
