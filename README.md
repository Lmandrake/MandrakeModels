# MandrakeModels

Gridded-cell physics models in Python: water (cellular automata, virtual pipes,
local-inertial, diffusion-wave, kinematic-wave, full shallow-water), waves,
diffusion, incompressible smoke (stable fluids) and lattice Boltzmann, all behind
one `Solver` interface on 2-D grids, numpy-vectorised.

Uses: a bakeoff of fluid solvers for game water (RimMandrake's FlowWorks is the
first customer), visualizations, animated backgrounds.

- Design: `design/specs/2026-11-12-fluid-solver-bakeoff.md`
- Library: `src/mandrakemodels/`
- Bakeoff: `uv run python -m mandrakemodels.bakeoff` (writes under the agent
  storage root, `mandrakemodels/`; override with `MANDRAKEMODELS_OUT`)
- Tests: `uv run pytest`
- Curated results: `results/`
