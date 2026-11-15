# Scenes

`scenes.build(name, n, seed=0)` returns a `Scene`. The same inputs always give the same scene: random terrain
comes from `smooth_noise`, a seeded smooth field. A `Scene` holds:

- `name` and `kind`
- `grid`
- `duration`: simulated seconds
- `init`: the starting arrays
- `forcing`: for water scenes
- `events`: a list of `(t, method name, kwargs)` called on the solver when its clock passes `t`
- `params`: scene-level physics knobs such as `manning`, `c`, `D` or `buoyancy`
- `settle_tol`: max |dh/dt| in m/s that counts as settled
- `still`: True if nothing should move
- `note`: a one-line description

Water scenes are `n x n` cells at `dx = 1.414 m`. The lists `scenes.WATER`, `scenes.OTHER` and `scenes.ALL`
name them. The numbers below are read from `scenes.py`.

## Water

| scene | duration | what it is | forcing and events | what it tests |
|---|---|---|---|---|
| `lake_at_rest` | 60 s | closed bowl on a rough bed, with an island; surface at 1.5 m everywhere | none; `still = True`, `settle_tol = 1e-6` | nothing may move: the well-balanced test |
| `dam_break` | 40 s | flat bed; a 2-cell-thick dam at x = n/3 with a gap of n/4 in the middle; 2 m of water behind it onto 0.3 m in front | none | flood front speed and breach throughput (wet-bed dam break) |
| `river` | 120 s | meandering channel on a 1e-3 slope between noisy banks | a source at the upstream channel end (`0.6 * width * dx` m^3/s), a free-outfall sink at the downstream end | steady flow, inflow and outflow bookkeeping |
| `lake_levelset` | 120 s | bowl lake at 1.2 m with an outlet channel to the east edge | `FixedLevel` at 0.9 m on the channel's edge cells | draining to a held level |
| `ocean_edge` | 120 s | beach rising from -6 m in the west, with a river channel entering from the east | sea held at 0 m on the west 3 columns; a river source at the east end | a level-setter coast plus a river |
| `pit_dig` | 150 s | straight river at y = 0.35 n | source and sink as `river`; at t = 5 s a 3 m pit is dug at y = 0.6 n, beside the river; at t = 15 s a canal is dug from the river to the pit, its bed falling 0.3 m toward the pit | terrain edits and filling a dug hole |
| `moat_evap` | 120 s | closed square moat whose floor tilts 0.8 m west to east | evaporation of 15 mm/s (exaggerated on purpose, so it dries in a short run) | evaporation and a receding shoreline |
| `rain_groundwater` | 120 s | dry rough basin | 4 mm/s rain for the first 30 s; groundwater table at 0.2 m, `k = 2e-2`, seeping up after | rain, runoff and seepage |
| `swell` | 60 s | beach from -8 m with two breakwater walls at x = 0.45 n | the west 2 columns held at `0.3 sin(2 pi t / 8)` m; `settle_tol = 0` (never settles) | waves on a beach and behind breakwaters |

`compare` shows only the first 36 s of `dam_break` (scaled with n/128 above 128 cells), to keep the event and
skip the settled tail.

## Wave, scalar and flow

| scene | kind | grid | duration | what it is |
|---|---|---|---|---|
| `ripples` | wave | n x n, dx 1 | 20 s | pond with a round island and a slotted wall; 8 drops at seeded times within the first 15 s; `c = 10`, `damping = 0.05` |
| `gaussian` | scalar | n x n, dx 1 | `0.5 (0.1 L)^2 / D` | Gaussian blob (sigma 0.05 L) diffusing with D = 1; the exact solution is known |
| `dye` | scalar | n x n, dx 1 | `0.5 (0.06 n)^2 / D` | 12 sharp dye discs diffusing between insulated walls; total dye must be conserved |
| `smoke_plume` | flow | n x n, dx 1 | `3 n` steps | buoyant smoke rising from a vent past a disc and a shelf, in a closed box |
| `cylinder_wake` | flow | (n/2, at least 16) x n | `8 n` steps | channel flow past a cylinder at Re about 150, with dye streaks from the inlet: a Karman vortex street |

Flow scenes set `params["steps"]`, and the tools run them for that many steps rather than for `duration`.

## Scenes used only by the C# export

`csharp_export.py` defines two more water scenes that are not in `scenes.ALL`:

- `census_river`: a river sized to a target wet area per map size
- `stress_half`: half the map under water up to 10 m deep, with a swelling west edge

Both are timing loads for the C# port.

## Your own grid

`compare` builds a water scene from your own arrays: `--init file.npz` (with `z`, `h` and optionally `wall`),
or `--bed` with `--depth` or `--level`. See `compare.md`.
