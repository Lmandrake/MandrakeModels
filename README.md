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

## Provenance

Every solver here is our own implementation written from the published method; no code was copied
from another project, so nothing in `src/` carries a third-party licence. Citations are from memory
and have not been re-checked against the papers.

| solver | method | source |
|---|---|---|
| `ca_game` | falling/levelling cellular automaton | game practice (Noita GDC 2019 talk; Terraria/Minecraft liquid rules), no paper |
| `wca2d` | weighted cellular automaton | Guidolin et al. 2016, Env. Modelling & Software 84 |
| `pipes` | virtual pipes | O'Brien & Hodgins 1995; Mei, Decaudin & Hu 2007 |
| `local_inertial` | local-inertial SWE | Bates, Horritt & Fewtrell 2010; de Almeida et al. 2012 |
| `kinematic_wave` | kinematic wave, Manning | Lighthill & Whitham 1955; Singh 1996 |
| `diffusion_wave_explicit` | zero-inertia wave, explicit | Hromadka & Yen 1986; Hunter et al. 2005 |
| `diffusion_wave_implicit_cg` / `_adi` | zero-inertia wave, backward Euler | Panday & Huyakorn 2004; Kollet & Maxwell 2006; CG (Hestenes & Stiefel 1952); ADI (Peaceman & Rachford 1955) |
| `swe_hll` | full SWE, HLL + hydrostatic reconstruction | Harten, Lax & van Leer 1983; Toro 2009; Audusse et al. 2004 |
| `swe_lbm` | shallow-water lattice Boltzmann | Zhou 2004 |
| `linear_wave` | damped wave equation, leapfrog | Hairer, Lubich & Wanner 2003 |
| `heat_explicit` / `heat_adi` | FTCS / Peaceman-Rachford ADI | Peaceman & Rachford 1955; Thomas 1949 |
| `stable_fluids` | semi-Lagrangian + projection | Stam 1999 |
| `lbm_d2q9` | D2Q9 BGK / TRT | Qian, d'Humieres & Lallemand 1992; Guo et al. 2002; Ginzburg 2005 |

Runtime dependencies: numpy (BSD-3), scipy (BSD-3), Pillow (MIT-CMU); optional numba (BSD-2).

### Cross-check libraries (optional, `uv sync --extra crosscheck`)

`src/mandrakemodels/external.py` calls these; it copies nothing from them. None is a core dependency.

| library | licence (from package metadata) | used for | status |
|---|---|---|---|
| Landlab 2.11 (`OverlandFlow`, `KinematicWaveRengers`, `LinearDiffusionOverlandFlowRouter`) | MIT | water tiles `ext_landlab_*` | runs; see `results/crosscheck.md` |
| ANUGA 4.1.0 (DE0 finite-volume SWE on triangles) | Apache-2.0 | water tile `ext_anuga` | runs (needs the `mpich` wheel for its MPI import); see `results/crosscheck.md` |
| PhiFlow 3.4 / phiml 1.16 | MIT | smoke tile `ext_phiflow_smoke` | runs on `smoke_plume`; see `results/crosscheck.md` |
| pylbm 0.12 | BSD-3-Clause | LBM | installed; no adapter: its geometry is built from analytic shapes, not cell masks |
| Clawpack/PyClaw 5.14 | BSD-3-Clause (PyPI) | - | not installed: PyPI has only a macOS-arm64 wheel and an sdist, and this machine has no Fortran compiler |
| lbmpy 2.1 / pystencils 2.1 | GPL-3.0-or-later | - | not used: GPL; never imported by this package |

## Side-by-side comparison

```sh
uv run python -m mandrakemodels.compare --scene dam_break --n 128              # all water solvers
uv run python -m mandrakemodels.compare --scene all --n 96 --external           # every scene, plus open-source tiles
uv run python -m mandrakemodels.compare --bed my_bed.npy --depth my_depth.npy --name mine --n 128
```

One initial grid is built once and copied to every solver; each runs on its own stable time step and
is sampled at common times (120 frames by default, over the active part of the event), so every tile
shows the same simulated moment. Water is drawn as water (`waterview.py`): the bed as shaded relief
in earth tones, water in blues that darken with depth, the water surface lit from the north-west so
wave fronts and bores show as light and dark bands, a pale shoreline, walls as dark stone. One look
(terrain range, depth scale, shading gain) per scene, taken from the reference run. Dam break and
swell get a side-view strip of the water surface along the row marked with yellow triangles.

Output under the output root: `compare/index.html` (water scenes first, other physics after) and one
self-contained page per scene, `compare/<scene>_<n>.html` (play/pause, frame step, time slider, speed,
"what the water does" / "difference from reference" toggle, results table); both work from `file://`.
`compare/<scene>_<n>/` holds `state.mp4`, `state.webm`/`diff.webm`, small GIFs, three stills,
`init.npz` (the exact starting grid) and `meta.json`; the stills are also copied to `compare/stills/`.
A solver that blew up, crashed or hit the wall-time cap keeps its tile, frozen at its last good frame,
under a red (blew up) or orange (stopped) banner.
