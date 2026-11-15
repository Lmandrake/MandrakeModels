# Side-by-side comparison

`python -m mandrakemodels.compare [options]` (`src/mandrakemodels/compare.py`) starts every applicable solver
from one initial grid and renders them as one tiled video per scene, with an HTML player.

## What makes the comparison fair

- The scene is built once. Each solver gets a deep copy of the same `Scene`: the same bed, water, walls,
  forcing, events and seed. Any difference between tiles is the solver.
- Water, wave and scalar solvers are sampled at the same simulated times (`--frames` + 1 samples over the
  duration).
- Between samples, each solver takes equal steps no longer than its own `stable_dt` (water: also no longer
  than `--dt-cap`, 5 s by default) and stops at each scene event.
- Flow scenes are sampled at equal step counts instead. `stable_fluids` and `lbm_d2q9` do not share a time
  unit, so their tiles compare shape, not speed.
- The reference is `swe_hll` (water), `heat_explicit` (scalar), `linear_wave` (wave) or `stable_fluids`
  (flow), unless `--ref` names another. The reference runs first and sets the colour scales for every tile.
- A solver that throws, blows up or exceeds `--max-wall` keeps its last good frame:
  - "blown up" or "error" gets a red banner
  - "stopped" (the wall cap) gets an orange banner
  - it is never dropped from the comparison

`dam_break` is shown for its first 36 s (scaled by n/128 above 128), since the settled tail is not
interesting. `--duration` overrides any scene's length.

## Options

| option | default | meaning |
|---|---|---|
| `--scene` | `dam_break` | a scene name, a comma list, `water` or `all` |
| `--n` | 128 | grid size |
| `--solvers` | every solver of the scene's kind | comma list (with `--external`, the `ext_*` tiles are added too) |
| `--exclude` | | comma list to leave out |
| `--ref` | per kind, as above | reference for the difference view |
| `--frames` | 120 | frames (sample intervals) |
| `--duration` | scene's | simulated seconds |
| `--max-wall` | 45 | wall seconds per solver per scene |
| `--tile` | at least 320 px, a whole multiple of n | tile size |
| `--seed` | 0 | scene seed |
| `--out` | `<output root>/compare` | output directory |
| `--html-only` | | rebuild the HTML pages from saved runs without running anything |
| `--dt-cap` | 5.0 | largest water dt, in s (it matters for the implicit solvers) |
| `--external` | | add tiles from Landlab, ANUGA and PhiFlow if they import (`uv sync --extra crosscheck`); the ones that do not import are listed with the reason |

### Your own starting grid (water)

| option | meaning |
|---|---|
| `--init file.npz` | arrays `z`, `h` and optionally `wall` (every water run saves one as `init.npz`) |
| `--bed file` | bed elevation, `.npy` in metres or `.png` with `--bed-range lo,hi` (metres for black and white) |
| `--depth file` / `--depth-range` | water depth, the same formats |
| `--level L` | fill to surface level L instead of `--depth` |
| `--wall file` | wall mask, `.npy` or `.png` (white = wall) |
| `--dx`, `--manning`, `--name` | cell size (default 1.414 m), Manning n (0.03), scene name (`custom`) |

A custom scene has no forcing and lasts 60 s unless `--duration` is given. PNG row 0 is taken as north and
flipped to array row 0 = south.

## Outputs

Under the output directory:

- `index.html`: cards for every saved scene, water first, each with a poster image
- `<scene>_<n>.html`: a self-contained page with the videos embedded as base64, so it works from `file://`.
  It has play/pause, frame step, a time slider, speed, a toggle between "what the water does" and "difference
  from reference", and a results table:
  - status, time reached and steps
  - ms/step
  - volume error
  - RMS against the reference
- `<scene>_<n>/` holds:
  - `state.mp4` (H.264), `state.webm` and `diff.webm` (VP9)
  - `state.gif`, `diff.gif`: about 40 frames
  - `poster.webp`
  - three state stills and three difference stills, `<scene>_<n>_{state,diff}_t<time>.jpg`
  - `init.npz` (water)
  - `meta.json`: scene, solvers, sample times, per-solver final status and errors
- `stills/`: copies of the state stills

The videos need `ffmpeg` on `PATH`. Without it no videos are written, and `index.html` lists only scenes that
have both `.webm` files. The GIFs, stills and `meta.json` are still written.

The difference view shows solver minus reference in the primary field (water depth, wave height, scalar
field, dye): red is more than the reference, blue is less.

## Checked

`compare --scene dam_break --n 32 --frames 8` was run for this document with `MANDRAKEMODELS_OUT` set. It
took about 23 s of wall time on a laptop and wrote every file listed above. `tests/test_compare.py` checks
the shared start, the outputs and a user grid from `.npy`.
