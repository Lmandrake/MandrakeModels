"""Solver registry. Every module in ``MODULES`` defines one or more ``Solver`` subclasses listed in its
``SOLVERS`` tuple. ``REGISTRY`` maps ``name -> class``."""
from __future__ import annotations

import importlib

MODULES = [
    "ca_game", "wca2d", "pipes", "local_inertial", "local_inertial_ripples", "kinematic", "diffusion_explicit",
    "diffusion_implicit", "swe_hll", "swe_lbm", "linear_wave", "heat", "stable_fluids", "lbm_d2q9",
]

REGISTRY: dict = {}
MISSING: list = []
for _m in MODULES:
    try:
        _mod = importlib.import_module(f"{__name__}.{_m}")
    except ModuleNotFoundError as e:
        if e.name == f"{__name__}.{_m}":
            MISSING.append(_m)
            continue
        raise
    for _cls in _mod.SOLVERS:
        REGISTRY[_cls.name] = _cls


def for_kind(kind: str):
    return [c for c in REGISTRY.values() if c.kind == kind]
