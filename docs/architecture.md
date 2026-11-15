# Architecture

Everything a solver shares lives in `src/mandrakemodels/base.py`. Scenes (`scenes.py`) are plain data that a
solver builds itself from. The tools (`bakeoff.py`, `compare.py`) only use the interface described here.

## Grid and array conventions

- `Grid(ny, nx, dx)`: `ny` rows, `nx` columns, square cells of width `dx` metres. `grid.area` is `dx^2`, and
  `grid.coords()` returns cell-centre coordinates `(X, Y)` in metres.
- Arrays are `(ny, nx)`. Row `j` is y and column `i` is x. Row 0 is the south edge: `compare` flips PNG input
  so that image row 0 (north) lands on the last array row.
- Water scenes use `dx = 1.414 m` (`scenes.DX`), so one cell is 2 m^2, the area of one RimWorld square. The
  map is `n x n` cells, so a bigger `n` is a bigger place, not a finer grid.
- Wave, scalar and flow scenes use `dx = 1.0`.

## Units

SI throughout for water: metres, seconds, m^3. `G = 9.81`. A cell with depth below `EPS_H = 1e-6 m` counts
as dry, and several solvers and diagnostics also use a 1 mm wet threshold (`h > 1e-3`).

The exceptions:

- `stable_fluids` works in cell units: `dx = 1` per cell, `dt = 1` per step, velocities in cells per step.
- `lbm_d2q9` works in lattice units, with `stable_dt() = 1`.

The two flow solvers therefore do not share a time unit. `compare` samples flow scenes by step count, not by
time.

## The `Solver` interface

```python
class Solver:
    name: str          # registry key, e.g. "swe_hll"
    family: str        # human-readable family
    kind: str          # "water" | "wave" | "scalar" | "flow"
    conserves: bool    # does the core step conserve the primary quantity exactly?

    def __init__(self, grid): ...          # sets t = 0, nsteps = 0, stats = {}
    @classmethod
    def from_scene(cls, scene, **kw): ...  # build from a scenes.Scene (subclasses define it)
    def stable_dt(self) -> float: ...      # largest dt the scheme's own stability bound allows now (inf if none)
    def step(self, dt): ...                # _step(dt), then t += dt, nsteps += 1
    def _step(self, dt): ...               # the scheme; subclasses implement it
    def blown_up(self) -> bool: ...
    def diagnostics(self) -> dict: ...     # t, steps, plus stats
    def snapshot(self) -> dict: ...        # arrays for rendering
```

The caller chooses `dt`. A solver never sub-steps on its own to meet `stable_dt` (the external adapters are
the exception: their libraries sub-step internally). `stable_dt` is "the step the scheme claims is stable",
not "an accurate step". The implicit solvers return a large fixed value (`max_dt`, 5 s by default), and the
bakeoff's `maxdt` phase measures how far past `stable_dt` each solver really survives.

`snapshot()` keys by kind:

| kind | keys |
|---|---|
| water | `z`, `h`, `wall`, `speed` |
| wave | `height`, `obstacle` |
| scalar | `field` |
| flow | `dye`, `speed`, `vorticity`, `obstacle` (and more for `lbm_d2q9`) |

## `WaterSolver`

Every water model subclasses `WaterSolver`. It holds:

- bed `z` (m), true depth `h` (m, never negative) and a boolean `wall` mask (no water, no flux). The surface
  is `eta = z + h`.
- face fluxes `qx` of shape `(ny, nx+1)` and `qy` of shape `(ny+1, nx)`: signed volume flux per unit width
  (m^2/s), positive toward +x / +y. These are for velocity and display. A solver with memory (local
  inertial, pipes) keeps its own discharge arrays as well.
- `manning` (default 0.03), `forcing`, and a `Ledger`.

The outer faces of the map are closed. Water enters or leaves only through a booked forcing term.

### The one conservative transfer

`WaterSolver.transfer(Fx, Fy)` is the only way the core step moves water.

- `Fx` has shape `(ny, nx-1)` and `Fy` has shape `(ny-1, nx)`: m^3 to move this step through each interior
  face.
- Faces touching a wall are zeroed.
- Each donor cell's total outflow is scaled down to what it holds (the donor limiter). The count of faces the
  limiter cut is in `limited_faces`.
- Every face is then applied once, with opposite signs to its two cells.

So depth never goes negative and nothing leaks. `volume_error()` (volume now, minus starting volume, minus
the ledger's net) is exact bookkeeping up to round-off, and the tests hold every water solver to it.

Shared helpers:

- `face_surface_diff()`: `eta[i+1] - eta[i]` on interior faces
- `face_flow_depth()`: the Bates (2010) face flow depth `max(eta) - max(z)`, at least 0 and zero across walls
- `wave_dt(cfl)`: `cfl * dx / sqrt(g * h_max)`
- `velocity()`: cell-centre `u, v` from the face fluxes, zero where `h < 1 mm`
- `dig(mask, depth)` and `raise_bed(mask, height)`: change the bed and keep each cell's water volume (scene
  events call these)

### Forcing and the ledger

`WaterSolver.step(dt)` runs `_step(dt)` and then `apply_forcing(dt)`, in this order:

| term | dataclass | effect, booked as |
|---|---|---|
| sources | `Source(mask, Q)` | `Q` m^3/s (number or function of t) spread evenly over the open cells of `mask`; `source` |
| rain | `Forcing.rain` | m/s on every open cell (number or function of t, sampled at the step midpoint); `rain` |
| groundwater | `Groundwater(table, k, mask)` | seepage `k * (table - eta)` m/s where the surface is below the table, capped so it never overshoots; `groundwater` |
| evaporation | `Forcing.evaporation` | m/s removed from wet cells, never more than is there; `evaporation` |
| sinks | `Sink(mask, rate)` | up to `rate` m/s of depth drained from `mask` (`math.inf` empties it every step); `sink` |
| fixed levels | `FixedLevel(mask, level)` | the surface on `mask` set to `level` (number or function of t): a level-setter, ocean edge or wave maker. Water added and water removed are booked separately; `level` |

`Ledger.net` is in minus out. `Ledger.gross` (also `throughput`) is the sum of the absolute value of every
booking, so a level-setter that adds and removes water does not cancel out. The bakeoff's relative volume
error divides by the larger of the starting volume and the throughput.

## Stepping rule: equal steps

`base.even_dt(dt, left)` splits the time `left` to the next target (a frame, checkpoint or event) into equal
steps no longer than `dt`, instead of full steps plus a short tail. The bakeoff, `compare` and
`csharp_export` all step this way. The code comments give the reason: on `dam_break`, a short tail step after
full ones grew `local_inertial`'s reservoir roughness about 500-fold by 33 s. The C# port states the same rule
at `LocalInertialDL.EvenDt`.

## Non-water solvers

`linear_wave`, `heat_explicit`, `heat_adi`, `stable_fluids` and `lbm_d2q9` subclass `Solver` directly. Each
defines its own `from_scene`, which reads `scene.init` (for example `field`, `obstacle`, `source`, `drops`)
and `scene.params` (for example `D`, `c`, `damping`, `buoyancy`, `inflow`, `Re`). Heat solvers expose
`total()`, which `compare` uses for their conservation error.

## The registry

`solvers/__init__.py` imports each module named in `MODULES`. Every module ends with a `SOLVERS` tuple of its
classes, and each class is entered in `REGISTRY` under its `name`.

- `for_kind(kind)` returns the classes of one kind.
- A module that fails to import because the module file itself is missing is listed in `MISSING`. Any other
  import error is raised.
- `local_inertial_ripsense` has a module and tests but is not in `MODULES`, so it is not registered.

To add a solver: write a subclass with `name`, `family` and `kind`, implement `_step` and `stable_dt`
(water solvers move water only through `transfer`), add it to its module's `SOLVERS`, and add the module to
`MODULES` if it is new. `tests/test_water_common.py` is parametrised over every registered water solver, so a
new water solver is tested for exact volume and a lake at rest from that point on. Solvers that cannot keep a
lake at rest by construction are listed there as expected failures.

## External adapters

`external.py` wraps Landlab, ANUGA and PhiFlow as `WaterSolver` / `Solver` subclasses named `ext_*`. Each
step pushes our bed and depth into the library, advances it by `dt`, and pulls the depth back. Our forcing and
ledger then run as for any water solver, so volume error is measured the same way. They are marked
`conserves = False` ("measured, not assumed"). `external.available()` returns the adapters whose library
imports, and `compare --external` adds them as tiles. They are never in `REGISTRY`, so the bakeoff does not
run them.
