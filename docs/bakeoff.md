# The bakeoff

`python -m mandrakemodels.bakeoff PHASE [options]` (`src/mandrakemodels/bakeoff.py`) runs every applicable
registered solver on every scene of its kind, at each size. External adapters are not included. Each phase
writes into `<output root>/bakeoff/<run>/`; see the README for how the root is chosen.

The bakeoff sets `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `NUMBA_NUM_THREADS` to 1,
unless they are already set, so numpy and numba run single-threaded.

## Options

| option | default | meaning |
|---|---|---|
| `--run` | `main` | run name; the output folder is `bakeoff/<run>/` |
| `--sizes` | `64,128,250,325,400,500` | grid sizes n (n x n) |
| `--solvers`, `--scenes` | all | comma lists that narrow the task list |
| `--budget` | 480 | wall seconds for this call; the phase stops cleanly when it is spent |
| `--max-wall` | 150 | wall cap per run (quality, timing); a quality run that hits it is flagged `capped` |
| `--workers` | 8 | quality: parallel worker processes |
| `--frame-size` | 250 | quality: render a GIF and stills only for runs at this n |
| `--steps`, `--warm` | 30, 3 | timing: timed steps after warm-up steps |
| `--probe-n` | 400 | probe: grid size |
| `--curate` | none | report: copy the small tables and contact sheets into this directory |

## Resumability

`quality` and `timing` append one JSON line per (solver, scene, n) to `quality.jsonl` or `timing.jsonl`. On
start they read the keys already present and skip them. `maxdt` does the same per solver. A phase that runs
out of `--budget` stops and says so; call it again to continue.

The report keeps the last line for each key. To redo a run, delete its line, or start under a new `--run`.

## Phases

### `quality`

Each run covers the scene's full duration (flow scenes: their step count), with scene events applied on time.
Water steps follow `even_dt(min(stable_dt, 5 s), time to the next event or the end)`.

Each line records:

- steps, simulated end time, wall seconds
- `blew_up` / `blow_t`, `error`, `capped`
- for water:
  - `settled` and `settle_t`: the last time |dh/dt| was at or above the scene's `settle_tol`, or still within
    an event
  - `vol_err_rel`: |volume error| over the larger of the starting volume and the throughput
  - `min_h`, `max_h`, `wet_frac`, `max_speed`
  - `still_drift_m` (only for `lake_at_rest`)
- `diag`: the solver's diagnostics

At `n == --frame-size` it also saves `frames/<scene>/<solver>_<n>.gif` and three PNG stills (`_a`, `_b`,
`_c`). Runs go to a process pool (`fork` context), smallest sizes first. Because they compete for cores, the
recorded `ms_per_step_parallel` is not used as a timing.

### `timing`

One process, one run at a time, single-threaded. Each run does `--warm` steps and then `--steps` timed steps.
The time of `stable_dt()` counts as part of each step, because some solvers recompute velocities in it.

Each line records:

- `ms_per_step` (median), `ms_per_step_mean`
- `dt_mean`
- `ns_per_cell_step`
- for water, `py_ms_per_tick`: ms/step / dt / 60, which assumes one simulated second per 60 game ticks
- for the CG solver, `iters`

### `maxdt`

For each solver, the largest fixed dt that survives on a reference scene:

- water: `dam_break` at 128 for 20 s
- scalar: `gaussian` at 64 for its duration
- wave: `ripples` at 64 for 10 s

Flow solvers are skipped.

"Survives" means no blow-up and no NaN, and:

- water: max depth never above 1.25x the start
- scalar: no new extreme
- wave: |u| < 10

The search starts at the solver's own `stable_dt` (0.5 s if it has none), halves until a run survives,
doubles until one fails (capped at 64 s), then bisects 6 times in log space. A dt that would need more than
40000 steps is not searched. `ca_game` is skipped (a fixed pass). `kinematic_wave`'s value is marked
meaningless, because the dam-break bed is flat.

The generated `results.md` says it plainly: bounded is not accurate. A limiter can keep a solver
bounded at a dt where its wave speeds are already wrong.

Output: `maxdt.jsonl`, one line per solver with `max_dt`, `note`, `own_stable_dt` and `wall_s`.

### `probe`

Needs `numba` (`uv sync --extra fast`). It times plain-loop numba kernels against their numpy versions on the
`river` scene at `--probe-n`. The kernels are a local-inertial face update, a virtual-pipes update and 5-point
FTCS heat. Output: `probe.json`.

The kernels are simplified solvers: pipes has no friction, and local inertial has no theta weighting or
coupled friction. The ratio therefore compares different arithmetic, and the report does not project a C#
cost from it. Measured C# cost is in `results/stagec.md`.

### `accuracy`

Every water solver runs `dam_break` at 128 for 20 s at its own stable dt. The phase measures the m^3 that
crosses the dam line and divides by `swe_hll`'s. A ratio well below 1 means floods arrive late; above 1, early.
Output: `accuracy.json`. With `--solvers` it re-runs only those rows and keeps the others.

### `report`

Merges what exists into the following files (missing phases are left out):

- `results.csv`: one row per (solver, scene, n), quality and timing fields joined
- `results.md`, with these sections:
  - the probe table
  - a water headline at 325 and at the largest timed size: the median over water scenes, with
    `kinematic_wave` left out
  - a scaling table
  - the maxdt table
  - the flood-timing table
  - a per-scene table at 250
- `sheets/<scene>_<n>.png`: contact sheets of every solver's mid and final frame
- `gallery.html`: the GIFs by scene

With `--curate DIR` it copies `results.csv`, `results.md`, `probe.json`, `maxdt.jsonl` and downscaled sheets
into `DIR`. `results/` in this repository was made that way.

## Cost

The full matrix (29 solvers, 14 scenes, sizes up to 500) takes hours. `diffusion_wave_explicit` alone hits
the wall cap on every water scene at large sizes, because of its sub-millisecond step on flat water. To try
the harness cheaply, narrow it, for example:

```sh
uv run python -m mandrakemodels.bakeoff timing  --run try --sizes 32 --solvers swe_hll,local_inertial --scenes dam_break
uv run python -m mandrakemodels.bakeoff report  --run try
```

The commands shown in this document were not all run for it. `report` (on an empty run), `probe --probe-n 32`
and `maxdt --solvers heat_adi,ca_game` were checked to run and write their files. The narrowed `timing` call
above was not run.
