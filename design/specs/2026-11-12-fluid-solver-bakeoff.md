# Fluid solver bakeoff: design (2026-11-12)

Ticket `MANDRAKEMODELS_BUILD_FLUID_BAKEOFF_ONE_INTERFACE_MANY_WATERS`. Owner (2026-11-12): a Python bakeoff of
fluid solvers on gridded cells, "useful for visualizations, backgrounds, and a variety of other realizations",
first customer RimMandrake FlowWorks (rivers, lakes, oceans, maybe waves, on big maps).
Inputs: the FlowWorks literature review (2026-11-12) and RimMandrake's bakeoff plan
(`design/plans/2026-11-12-flowworks-water-solver-bakeoff.md` in RimMandrake; read, not edited).

## 1. Library shape

- `mandrakemodels.grid.Grid`: `ny, nx, dx` (metres). Cell-centred arrays `(ny, nx)`; face arrays
  `qx (ny, nx+1)`, `qy (ny+1, nx)` hold signed volume flux per unit width (m^2/s). Outer faces are walls.
- `mandrakemodels.base.Solver`: `step(dt)`, `stable_dt()`, `diagnostics()`, `snapshot()` (render fields),
  class attributes `name`, `family`, `kind` (`water`, `wave`, `scalar`, `flow`), `conserves` (bool),
  docstring = physics approximated and what it cannot do.
- `mandrakemodels.base.WaterSolver(Solver)`: state `z` (bed, m), `h` (true depth, m), `wall` (bool);
  boundaries and forcing applied by the base class after each core step, every external volume booked in
  a `Ledger` so `volume_error = V - V0 - (in - out)` is exact bookkeeping, not a guess:
  - `Source(mask, Q)` m^3/s spread over the mask; `Sink(mask, rate)` drains up to `rate` m/s (inf = all);
  - `FixedLevel(mask, level_fn(t))` holds the surface (lake level-setter, ocean edge, swell);
  - `rain` m/s, `evaporation` m/s on wet cells, `Groundwater(table, k)` seeps up while the surface is
    below the table;
  - events at time t: `dig(mask, depth)` / `fill(mask, height)` change the bed and keep each cell's
    water volume (conservative by construction).
- numpy-vectorised; numba used only for the compiled unit-cost probe (C# projection basis), optional.
- Deterministic: all randomness from `numpy.random.default_rng(seed)` in the scene builders.

## 2. Solvers (one module each)

Water (depth on a bed): game cellular automaton (falling/levelling, FlowWorks `diff/5` kind);
WCA2D weighted CA (Guidolin 2016); virtual pipes (Mei 2007 outflow pipes, O'Brien 1995 lineage);
local-inertial (Bates 2010, de Almeida 2012 theta-weighting, coupled vector friction);
kinematic wave (bed slope, Manning); explicit diffusion wave (regularised slope, adaptive dt);
implicit diffusion wave (design-7 kind: backward Euler on the surface, CG solve, conservative
corrector) and its ADI line-solve variant; full 2-D shallow water (HLL finite volume, Audusse hydrostatic
reconstruction, positivity, semi-implicit Manning).
Waves: linear wave equation heightfield (leapfrog, damping, c^2 = g h optional).
Scalar: heat/diffusion explicit (FTCS) and implicit ADI (Peaceman-Rachford).
Flow (incompressible / mesoscopic): stable fluids (Stam 1999 semi-Lagrangian, pressure projection);
lattice Boltzmann D2Q9 BGK; shallow-water LBM (Zhou 2004) if time allows.
Not solvers, cut with reason: subgrid bathymetry and fill-spill-merge (orthogonal accelerations the
review says to measure separately); SPH (not gridded).

## 3. Scenes (data: a `Scene` = bed, initial depth, walls, boundaries, events, duration, render hints)

Water: `lake_at_rest` (rough bed, closed, must stay still), `dam_break` (wet bed, Stoker-type),
`river` (sloped channel, source upstream, sink downstream), `lake_levelset` (lake with an edge level
setter), `ocean_edge` (sloping beach, sea held at a level on one edge), `pit_dig` (river beside a pit dug
mid-run through a canal), `moat_evap` (closed moat, evaporation), `rain_groundwater` (dry basin, rain then
groundwater refill), `swell` (ocean edge level oscillating, waves on a beach).
Wave: `ripples` (drops on a pond with an island). Scalar: `gaussian` (analytic diffusion check), `dye`.
Flow: `smoke_plume` (buoyant smoke past obstacles), `cylinder_wake` (LBM Karman street).
Sizes: 64, 128, 250, 325, 400, 500 cells a side (RimWorld: 250 common, 325 largest normal, 400 dev max).
dx = 1.414 m for water scenes (a RimWorld square is 2 m^2, RimMandrake plan section 2).

## 4. Harness and metrics

`python -m mandrakemodels.bakeoff` runs every applicable solver x scene x size. Per run: ms per step
(median, after warm-up), steps and ms per simulated second, steps to settle (max |dh/dt| < tol), volume
error (absolute and relative to throughput), min depth, blow-up (NaN/inf/negative beyond tol),
max stable dt (bisection on the dam break at 128). Big sizes run a capped number of steps; a capped run
says so. Output: `results.csv`, `results.md`, per run 4 frames + one GIF, under the agent storage root
`mandrakemodels/bakeoff/` (`MANDRAKEMODELS_OUT` overrides). Curated tables and a few images are copied
into `results/` in the repo.

Timing is sequential (one process, numpy single-threaded) so runs do not contend; quality runs may go in
parallel. Projected C# ms per tick is a labelled guess: numpy ns/cell-step x (numba/numpy ratio measured on
probe kernels) x a Mono factor of 2-4 [guess], at 1 simulated second per 60 ticks.

## 5. Proof

pytest per solver: conservation (where the physics conserves), still water stays still, mirror symmetry,
analytic checks (Gaussian diffusion against the exact solution, LBM Poiseuille profile, linear-wave
energy). Then the full bakeoff at all sizes, and every solver's frames viewed. A solver that cannot be
made stable is reported as such.

## 6. Gallery

`gallery.html` beside the results: per scene, one row of GIFs (one per solver) at 250, plus the table.
