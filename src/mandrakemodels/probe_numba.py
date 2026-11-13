"""Compiled unit cost: plain-loop numba kernels of three solvers, timed beside their numpy versions.

The kernels are SIMPLIFIED solvers, not the numpy arithmetic: the pipes kernel has no friction and the
local inertial kernel has no theta weighting or coupled-|q| friction (stage-2 review, 2026-11-13). Their
numpy/numba ratio therefore prices different arithmetic and is no basis for a C# cost; the report no
longer projects one (measured C#: results/stagec.md). Kernels: local inertial face update plus
conservative divergence with a donor clamp; virtual pipes outflow update plus Mei scaling; 5-point FTCS
heat; single thread, float64.
"""
from __future__ import annotations

import math
import time

import numpy as np

try:
    import numba
    from numba import njit
except ImportError:  # pragma: no cover
    numba = None

G = 9.81


if numba is not None:
    @njit(cache=False, fastmath=False)
    def li_step(z, h, qx, qy, dx, dt, n2):
        ny, nx = h.shape
        for j in range(ny):
            for i in range(1, nx):
                el, er = z[j, i - 1] + h[j, i - 1], z[j, i] + h[j, i]
                hf = max(el, er) - max(z[j, i - 1], z[j, i])
                if hf > 1e-4:
                    q = qx[j, i]
                    qx[j, i] = (q - G * hf * dt * (er - el) / dx) / (1 + G * dt * n2 * abs(q) / hf ** (7 / 3))
                else:
                    qx[j, i] = 0.0
        for j in range(1, ny):
            for i in range(nx):
                el, er = z[j - 1, i] + h[j - 1, i], z[j, i] + h[j, i]
                hf = max(el, er) - max(z[j - 1, i], z[j, i])
                if hf > 1e-4:
                    q = qy[j, i]
                    qy[j, i] = (q - G * hf * dt * (er - el) / dx) / (1 + G * dt * n2 * abs(q) / hf ** (7 / 3))
                else:
                    qy[j, i] = 0.0
        r = dt / dx
        for j in range(ny):
            for i in range(nx):
                out = (max(qx[j, i + 1], 0.0) + max(-qx[j, i], 0.0) + max(qy[j + 1, i], 0.0)
                       + max(-qy[j, i], 0.0)) * r
                if out > h[j, i] and out > 0:
                    s = h[j, i] / out
                    qx[j, i + 1] = qx[j, i + 1] * s if qx[j, i + 1] > 0 else qx[j, i + 1]
                    qx[j, i] = qx[j, i] * s if qx[j, i] < 0 else qx[j, i]
                    qy[j + 1, i] = qy[j + 1, i] * s if qy[j + 1, i] > 0 else qy[j + 1, i]
                    qy[j, i] = qy[j, i] * s if qy[j, i] < 0 else qy[j, i]
        for j in range(ny):
            for i in range(nx):
                h[j, i] -= r * (qx[j, i + 1] - qx[j, i] + qy[j + 1, i] - qy[j, i])

    @njit(cache=False)
    def pipes_step(z, h, fl, fr, fd, fu, dx, dt):
        ny, nx = h.shape
        k = dt * dx * G
        A = dx * dx
        for j in range(ny):
            for i in range(nx):
                e = z[j, i] + h[j, i]
                fl[j, i] = max(fl[j, i] + k * (e - z[j, i - 1] - h[j, i - 1]), 0.0) if i > 0 else 0.0
                fr[j, i] = max(fr[j, i] + k * (e - z[j, i + 1] - h[j, i + 1]), 0.0) if i < nx - 1 else 0.0
                fd[j, i] = max(fd[j, i] + k * (e - z[j - 1, i] - h[j - 1, i]), 0.0) if j > 0 else 0.0
                fu[j, i] = max(fu[j, i] + k * (e - z[j + 1, i] - h[j + 1, i]), 0.0) if j < ny - 1 else 0.0
                tot = (fl[j, i] + fr[j, i] + fd[j, i] + fu[j, i]) * dt
                if tot > h[j, i] * A and tot > 0:
                    s = h[j, i] * A / tot
                    fl[j, i] *= s; fr[j, i] *= s; fd[j, i] *= s; fu[j, i] *= s
        for j in range(ny):
            for i in range(nx):
                inflow = 0.0
                if i > 0:
                    inflow += fr[j, i - 1]
                if i < nx - 1:
                    inflow += fl[j, i + 1]
                if j > 0:
                    inflow += fu[j - 1, i]
                if j < ny - 1:
                    inflow += fd[j + 1, i]
                h[j, i] += dt * (inflow - fl[j, i] - fr[j, i] - fd[j, i] - fu[j, i]) / A

    @njit(cache=False)
    def heat_step(u, un, a):
        ny, nx = u.shape
        for j in range(ny):
            for i in range(nx):
                c = u[j, i]
                s = 0.0
                if i > 0:
                    s += u[j, i - 1] - c
                if i < nx - 1:
                    s += u[j, i + 1] - c
                if j > 0:
                    s += u[j - 1, i] - c
                if j < ny - 1:
                    s += u[j + 1, i] - c
                un[j, i] = c + a * s
        u[:, :] = un


def _time(fn, reps):
    fn()
    ts = []
    for _ in range(reps):
        a = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - a)
    ts.sort()
    return ts[len(ts) // 2]


def run(n: int = 400, reps: int = 30) -> dict:
    from . import scenes
    from .solvers import REGISTRY
    if numba is None:
        return {"error": "numba not installed (pip install mandrakemodels[fast])"}
    sc = scenes.build("river", n)
    cells = n * n
    res = {"n": n, "cells": cells, "numba": numba.__version__, "basis":
           "median of %d steps, single thread; numba kernels are SIMPLIFIED (pipes: no friction; local inertial: no theta, no coupled friction)" % reps}
    z = sc.init["z"].astype(np.float64)
    h0 = sc.init["h"].astype(np.float64) + 0.2
    dx = sc.grid.dx
    out = {}
    # local inertial
    qx = np.zeros((n, n + 1)); qy = np.zeros((n + 1, n)); h = h0.copy()
    t_nb = _time(lambda: li_step(z, h, qx, qy, dx, 0.05, 0.03 ** 2), reps)
    s = REGISTRY["local_inertial"](sc.grid, z, h0.copy(), manning=0.03)
    t_np = _time(lambda: s.step(0.05), reps)
    out["local_inertial"] = (t_np, t_nb)
    # pipes
    fl, fr, fd, fu = (np.zeros((n, n)) for _ in range(4)); h = h0.copy()
    t_nb = _time(lambda: pipes_step(z, h, fl, fr, fd, fu, dx, 0.05), reps)
    s = REGISTRY["pipes"](sc.grid, z, h0.copy(), manning=0.03)
    t_np = _time(lambda: s.step(0.05), reps)
    out["pipes"] = (t_np, t_nb)
    # heat
    if "heat_explicit" in REGISTRY:
        g = scenes.build("dye", n)
        u = g.init["field"].astype(np.float64).copy(); un = np.empty_like(u)
        t_nb = _time(lambda: heat_step(u, un, 0.2), reps)
        s = REGISTRY["heat_explicit"].from_scene(g)
        t_np = _time(lambda: s.step(0.2), reps)
        out["heat_explicit"] = (t_np, t_nb)
    for k, (a, b) in out.items():
        res[k] = {"numpy_ms": 1000 * a, "numba_ms": 1000 * b, "ratio": a / b,
                  "numba_ns_per_cell": 1e9 * b / cells}
    ratios = [v["ratio"] for k, v in res.items() if isinstance(v, dict)]
    res["ratio_median"] = float(np.median(ratios))
    return res
