"""Falling/levelling cellular automaton, the water model most games use."""
from __future__ import annotations

import numpy as np

from ..base import WaterSolver


class GameCA(WaterSolver):
    """Game-style water cellular automaton (Minecraft/Terraria/Noita school), written on faces.

    Approximates: nothing physical. Each pass, every cell offers 1/5 of its surface-height difference
    to each lower neighbour (so one cell cannot give away more than 4/5 of a difference and never
    overshoots). The update is synchronous (Jacobi: all faces from the same old state), hence
    mirror-symmetric. Water crosses a face only where the surface differs and the donor stands above
    the receiver's bed, capped at the available face depth (Bates flow depth), so a still lake on a
    rough bed stays still. Volume moves only through ``transfer``: exact, never negative.

    Cannot do: time. One step is one pass; ``pass_seconds`` (default 0.5) only labels how long a pass
    is so fluxes and velocities have units. Spreading speed is one cell per pass whatever the depth or
    slope (no momentum, no friction, no waves, no overshoot), and levelling is diffusive: a long
    thin pond needs ~N^2 passes.
    """

    name = "ca_game"
    family = "cellular automaton (game)"

    def __init__(self, grid, z, h=None, wall=None, forcing=None, manning=0.03, pass_seconds=0.5,
                 share=0.2):
        super().__init__(grid, z, h, wall, forcing, manning)
        self.pass_seconds = pass_seconds
        self.share = share

    def stable_dt(self) -> float:
        return self.pass_seconds

    def _step(self, dt):
        A = self.grid.area
        dEx, dEy = self.face_surface_diff()      # eta[right]-eta[left]; positive: flows toward -x
        hx, hy = self.face_flow_depth()
        # volume per pass: share of the level difference, capped at the face depth (water above the
        # higher bed); sign toward +x/+y is minus the difference
        Fx = -np.sign(dEx) * np.minimum(np.abs(dEx) * self.share, hx) * A
        Fy = -np.sign(dEy) * np.minimum(np.abs(dEy) * self.share, hy) * A
        self.transfer(Fx, Fy)
        self.record_flux(*self._last_F, self.pass_seconds)


SOLVERS = (GameCA,)
