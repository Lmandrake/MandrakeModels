# MandrakeModels

Physics models on 2-D grids of cells, written in Python and vectorised with numpy. The package covers:

- water: cellular automata, virtual pipes, local-inertial, diffusion-wave, kinematic-wave, full shallow-water
  and shallow-water lattice Boltzmann
- linear waves
- heat and diffusion
- incompressible smoke (Stam's stable fluids)
- D2Q9 lattice Boltzmann flow

Every model sits behind one `Solver` interface. On top of the solvers there are three tools:

- a **bakeoff** that runs every applicable solver on every scene at several grid sizes and tabulates cost,
  stability, volume conservation and settling
- a **side-by-side comparison** that starts every solver from the same grid and renders the runs as videos
  and HTML pages
- a small **C# port** of the local-inertial solver, used to time it on the Mono runtime a game uses

## Why

The project began as a bakeoff to choose a water model for a game. The first customer is FlowWorks, a water
system for the author's RimWorld mod, which needs rivers, lakes, coasts and dug channels on maps of up to
about 500 x 500 squares. The question it answers is which model is cheap enough per game tick and still
behaves like water: a still lake stays still, volume is exact, floods arrive at roughly the right time, and
the model does not blow up. The solvers are also useful on their own for visualisations and animated
backgrounds.

The design spec is `design/specs/2026-11-12-fluid-solver-bakeoff.md`.

## Quick start

You need [uv](https://docs.astral.sh/uv/) and Python 3.11 to 3.13 (`.python-version` pins 3.12). Rendering
videos in `compare` needs `ffmpeg` on `PATH`. Without it the videos and the per-scene HTML pages are
skipped; the GIFs and stills are still written.

```sh
git clone https://github.com/Lmandrake/MandrakeModels.git
cd MandrakeModels
uv sync                       # core: numpy, scipy, Pillow, matplotlib; pytest in the dev group
uv sync --extra fast          # optional: numba, used only by the bakeoff "probe" phase
uv sync --extra crosscheck    # optional: Landlab, PhiFlow, ANUGA, mpich, for the open-source cross-check tiles

export MANDRAKEMODELS_OUT=$PWD/out     # where every run writes (out/ is git-ignored); see "Where output goes"

uv run pytest -q                                                             # about 4 minutes on a laptop
uv run python -c "from mandrakemodels.solvers import REGISTRY; print(sorted(REGISTRY))"
uv run python -m mandrakemodels.compare --scene dam_break --n 32 --frames 8   # a small comparison, under a minute
```

Open `$MANDRAKEMODELS_OUT/compare/index.html` in a browser. The pages work from `file://` and need no server.

The extras are separate `uv sync` calls here for clarity. Each `uv sync` installs only the extras it is
given and removes the others, so to keep both use `uv sync --extra fast --extra crosscheck`.

If you would rather keep the virtual environment outside the checkout, set `UV_PROJECT_ENVIRONMENT` to a
directory of your choice before `uv sync`. uv reads this variable; the project itself does not.

### A bakeoff run

The bakeoff is long: the full matrix at sizes up to 500 takes hours of CPU. Every phase is resumable, so you
can run it in pieces. This is the sequence used to make `results/`:

```sh
uv run python -m mandrakemodels.bakeoff probe    --run main          # needs --extra fast
uv run python -m mandrakemodels.bakeoff maxdt    --run main
uv run python -m mandrakemodels.bakeoff timing   --run main --sizes 64,128,250,325,400,500
uv run python -m mandrakemodels.bakeoff quality  --run main --sizes 64,128,250 --frame-size 250
uv run python -m mandrakemodels.bakeoff accuracy --run main
uv run python -m mandrakemodels.bakeoff report   --run main --curate results
```

Each phase skips the (solver, scene, size) keys it has already written and stops after `--budget` seconds
(480 by default), so run it again until it reports nothing left to do. `timing` runs on one thread, one run at
a time. `quality` runs in parallel worker processes (`--workers`, 8 by default), so its ms/step is not used as
a timing. `--solvers` and `--scenes` (comma lists) narrow any phase. Full details: `docs/bakeoff.md`.

## Solvers

29 solvers are registered. List them with the one-liner above. The class docstring of each one says what
physics it approximates and what it cannot do; `docs/solvers.md` collects them.

| kind | solvers |
|---|---|
| water | `ca_game`, `wca2d`, `pipes`, `pipes_theta`, `local_inertial`, `local_inertial_frcap`, `local_inertial_frcap07`, eleven `local_inertial_frcap_*` ripple-damper variants, `kinematic_wave`, `diffusion_wave_explicit`, `diffusion_wave_implicit_cg`, `diffusion_wave_implicit_adi`, `swe_hll`, `swe_lbm` |
| wave | `linear_wave` |
| scalar | `heat_explicit`, `heat_adi` |
| flow | `stable_fluids`, `lbm_d2q9` |

`swe_hll` is the reference for water: it matches Stoker's analytic dam break in the tests, and ANUGA, an
independent open-source solver, agrees with it (`results/crosscheck.md`).

`local_inertial_frcap_ripsense` (in `solvers/local_inertial_ripsense.py`) has tests but is not in the
registry, because its module is not listed in `solvers/__init__.py`. The bakeoff and `compare` do not see it.

## Provenance

Every solver here is our own implementation written from the published method; no code was copied from
another project, so nothing in `src/` carries a third-party licence. Citations are from memory and have not
been re-checked against the papers.

| solver | method | source |
|---|---|---|
| `ca_game` | falling/levelling cellular automaton | game practice (Noita GDC 2019 talk; Terraria/Minecraft liquid rules), no paper |
| `wca2d` | weighted cellular automaton | Guidolin et al. 2016, Env. Modelling & Software 84 |
| `pipes`, `pipes_theta` | virtual pipes | O'Brien & Hodgins 1995; Mei, Decaudin & Hu 2007; theta smoothing after de Almeida et al. 2012 |
| `local_inertial` and the `local_inertial_frcap*` variants | local-inertial SWE | Bates, Horritt & Fewtrell 2010; de Almeida et al. 2012; Froude cap and ripple dampers are our additions (see `docs/solvers.md`) |
| `kinematic_wave` | kinematic wave, Manning | Lighthill & Whitham 1955; Singh 1996 |
| `diffusion_wave_explicit` | zero-inertia wave, explicit | Hromadka & Yen 1986; Hunter et al. 2005 |
| `diffusion_wave_implicit_cg` / `_adi` | zero-inertia wave, backward Euler | Panday & Huyakorn 2004; Kollet & Maxwell 2006; CG (Hestenes & Stiefel 1952); ADI (Peaceman & Rachford 1955) |
| `swe_hll` | full SWE, HLL + hydrostatic reconstruction | Harten, Lax & van Leer 1983; Toro 2009; Audusse et al. 2004 |
| `swe_lbm` | shallow-water lattice Boltzmann | Zhou 2004 |
| `linear_wave` | damped wave equation, leapfrog | Hairer, Lubich & Wanner 2003 |
| `heat_explicit` / `heat_adi` | FTCS / Peaceman-Rachford ADI | Peaceman & Rachford 1955; Thomas 1949 |
| `stable_fluids` | semi-Lagrangian + projection | Stam 1999 |
| `lbm_d2q9` | D2Q9 BGK / TRT | Qian, d'Humieres & Lallemand 1992; Guo et al. 2002; Ginzburg 2005 |

Runtime dependencies: numpy (BSD-3), scipy (BSD-3), Pillow (MIT-CMU), matplotlib (PSF-style, used for
colormaps); optional numba (BSD-2).

### Cross-check libraries (optional, `uv sync --extra crosscheck`)

`src/mandrakemodels/external.py` calls these; it copies nothing from them. None is a core dependency.

| library | licence (from package metadata) | used for | status |
|---|---|---|---|
| Landlab 2.11 (`OverlandFlow`, `KinematicWaveRengers`, `LinearDiffusionOverlandFlowRouter`) | MIT | water tiles `ext_landlab_*` | runs; see `results/crosscheck.md` |
| ANUGA 4.1.0 (DE0 finite-volume SWE on triangles) | Apache-2.0 | water tile `ext_anuga` | runs (needs the `mpich` wheel for its MPI import); see `results/crosscheck.md` |
| PhiFlow 3.4 / phiml 1.16 | MIT | smoke tile `ext_phiflow_smoke` | runs on `smoke_plume`; see `results/crosscheck.md` |
| pylbm 0.12 | BSD-3-Clause | LBM | not adapted: its geometry is built from analytic shapes, not cell masks (not in the `crosscheck` extra) |
| Clawpack/PyClaw 5.14 | BSD-3-Clause (PyPI) | - | not used: PyPI had only a macOS-arm64 wheel and an sdist, and the machine the cross-check ran on had no Fortran compiler |
| lbmpy 2.1 / pystencils 2.1 | GPL-3.0-or-later | - | not used: GPL; never imported by this package |

## Side-by-side comparison

```sh
uv run python -m mandrakemodels.compare --scene dam_break --n 128              # all water solvers
uv run python -m mandrakemodels.compare --scene all --n 96 --external           # every scene, plus open-source tiles
uv run python -m mandrakemodels.compare --bed my_bed.npy --depth my_depth.npy --name mine --n 128
```

The scene is built once and a copy goes to every solver. Each solver runs on its own stable time step and is
sampled at common times: 120 frames by default, over the active part of the event. Every tile therefore shows
the same simulated moment.

Water is drawn as water (`waterview.py`):

- the bed as shaded relief in earth tones
- water in blues that darken with depth
- the water surface lit from the north-west, so wave fronts and bores show as light and dark bands
- a pale shoreline, and walls as dark stone

Each scene has one look (terrain range, depth scale, shading gain), taken from the reference run. Dam break
and swell also get a side-view strip of the water surface along the row marked with yellow triangles.

Output under `<output root>/compare/`:

- `index.html`: water scenes first, other physics after
- `<scene>_<n>.html`: one self-contained page per scene, with play/pause, frame step, a time slider, speed,
  a toggle between "what the water does" and "difference from reference", and a results table. Both pages work
  from `file://`.
- `<scene>_<n>/`: `state.mp4`, `state.webm`, `diff.webm`, small `state.gif`/`diff.gif`, `poster.webp`,
  three state stills and three difference stills (`.jpg`), `init.npz` (the exact starting grid; water only)
  and `meta.json`. The state stills are also copied to `compare/stills/`.

A solver that blew up, crashed or hit the wall-time cap keeps its tile, frozen at its last good frame, under
a red (blew up) or orange (stopped) banner. Full option list: `docs/compare.md`.

## Where output goes

Nothing a run makes is written into the repository unless you ask for it (`report --curate results`).
`src/mandrakemodels/paths.py` picks the output root in this order:

1. `MANDRAKEMODELS_OUT`, if set: used as is.
2. Otherwise `$AGENT_STORAGE_ROOT/mandrakemodels`, if `AGENT_STORAGE_ROOT` is set.
3. Otherwise, if the author's own storage helper exists at `~/dev/Lodestar/bin/lib/agent_storage.py`, the
   root it prints, plus `/mandrakemodels`. On any other machine that file does not exist and this step is
   skipped.
4. Otherwise it stops with `RuntimeError: no output root: set MANDRAKEMODELS_OUT or AGENT_STORAGE_ROOT`.
   It never falls back to a folder inside the repository.

So on a fresh clone, set `MANDRAKEMODELS_OUT` (for example to `$PWD/out`, which `.gitignore` already
covers). `compare --out DIR` overrides the root for one comparison. The bakeoff writes to
`<root>/bakeoff/<run>/` and the comparison to `<root>/compare/`.

## Repository layout

```
src/mandrakemodels/
  base.py          Solver and WaterSolver base classes, Grid, forcing (sources, sinks, levels, rain, ...)
  solvers/         one module per solver family; __init__.py builds REGISTRY
  scenes.py        the 14 test scenes, built from (name, n, seed)
  bakeoff.py       the bakeoff phases (CLI)
  report.py        merges bakeoff output into results.csv, results.md, contact sheets, gallery.html
  compare.py       the side-by-side tool (CLI)
  waterview.py     water drawn as water, for compare
  render.py        frames, GIFs and stills for the bakeoff
  external.py      adapters for Landlab, ANUGA and PhiFlow (optional)
  probe_numba.py   simplified numba kernels for the bakeoff probe phase
  csharp_export.py writes scenes and Python reference runs as binary files for the C# port
  paths.py         where output goes
tests/             pytest suite (analytic checks, conservation, lake at rest, harness)
results/           curated results: tables, contact sheets, small compare outputs, write-ups (results/README.md)
csharp/            the C# port of local_inertial for runtime timing (csharp/README.md)
design/specs/      the design spec
tickets/           the work item this repository was built under
docs/              the documentation below
```

## Documentation

- `docs/architecture.md`: the Solver interface, grid and face conventions, units, the conservative transfer,
  forcing and the ledger, the registry
- `docs/solvers.md`: each solver's physics, stability limit and limits
- `docs/scenes.md`: every scene
- `docs/bakeoff.md`: phases, resumability, output files
- `docs/compare.md`: the side-by-side tool and its options
- `docs/csharp.md`: the C# port (with `csharp/README.md`)
- `results/README.md`: what each curated result file is

## Licence

MIT. See `LICENSE`.
