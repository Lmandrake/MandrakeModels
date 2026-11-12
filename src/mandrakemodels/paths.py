"""Where scratch, renders and bakeoff output go.

Order: ``MANDRAKEMODELS_OUT`` if set; else ``<agent storage root>/mandrakemodels`` where the root comes
from ``AGENT_STORAGE_ROOT`` or Lodestar's ``agent_storage.py`` helper. Fails by name if neither resolves
(never falls back into the repo).
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

_HELPER = Path.home() / "dev" / "Lodestar" / "bin" / "lib" / "agent_storage.py"


def out_root() -> Path:
    env = os.environ.get("MANDRAKEMODELS_OUT")
    if env:
        p = Path(env)
    else:
        root = os.environ.get("AGENT_STORAGE_ROOT")
        if not root and _HELPER.exists():
            root = subprocess.run(["python3", str(_HELPER)], capture_output=True, text=True,
                                  check=True).stdout.strip()
        if not root:
            raise RuntimeError("no output root: set MANDRAKEMODELS_OUT or AGENT_STORAGE_ROOT")
        p = Path(root) / "mandrakemodels"
    p.mkdir(parents=True, exist_ok=True)
    return p
