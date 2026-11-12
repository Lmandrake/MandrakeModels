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
- Curated results: `results/` (tables, numba probe, contact sheets); the full run with GIFs and
  `gallery.html` stays under the output root, `bakeoff/<run>/`.

## Running

```sh
export UV_PROJECT_ENVIRONMENT=/mnt/d/ClaudeWork/mandrakemodels/venv   # keep the venv off ~/dev
uv sync --extra fast                                     # numba for the probe
uv run python -m mandrakemodels.bakeoff probe   --run main
uv run python -m mandrakemodels.bakeoff maxdt   --run main
uv run python -m mandrakemodels.bakeoff timing  --run main --sizes 64,128,250,325,400,500
uv run python -m mandrakemodels.bakeoff quality --run main --sizes 64,128,250 --frame-size 250
uv run python -m mandrakemodels.bakeoff report  --run main --curate results
```

Every phase is resumable: it skips (solver, scene, size) keys already in its jsonl and stops at
`--budget` seconds, so call it again until it reports nothing left to do. Timing runs one thread,
in sequence; quality runs in parallel processes, so its ms/step is not a timing.

## Solvers

`python -c "from mandrakemodels.solvers import REGISTRY; print(sorted(REGISTRY))"`. Each class
docstring states the physics it approximates and what it cannot do.
