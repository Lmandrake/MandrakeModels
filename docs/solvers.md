# Solvers

One section per module in `src/mandrakemodels/solvers/`, condensed from each class docstring and its
`stable_dt`. The docstrings are the primary source; read them for detail. Constructor keywords are given with
their defaults. Figures quoted as "measured" come from the docstrings and were measured by the author on the
bakeoff scenes; they were not re-run for this document.

Common to every water solver (see `architecture.md`): water moves only through the shared donor-limited
`transfer`, so volume is exact and depth never goes negative. Manning friction defaults to `manning = 0.03`.

## Water

### `ca_game`: game-style cellular automaton (`ca_game.py`)

- **Physics.** None, by intent. This is the falling and levelling rule most games use. In each pass, every
  cell offers 1/5 (`share = 0.2`) of its surface-height difference to each lower neighbour, capped at the
  Bates face flow depth. The update is synchronous, so it is mirror-symmetric, and a still lake on a rough bed
  stays still.
- **Step.** One step is one pass. `stable_dt()` returns `pass_seconds` (0.5), which only labels a pass with a
  duration so fluxes have units. A shorter `dt` moves that fraction of a pass.
- **Cannot do.** Time: water spreads one cell per pass whatever the depth or slope. There is no momentum, no
  friction and no waves. Levelling is diffusive: a long thin pond needs about N^2 passes. The `maxdt` phase
  skips it.

### `wca2d`: weighted cellular automaton (`wca2d.py`)

- **Physics.** Diffusion-like flood spreading in the style of WCA2D (Guidolin et al. 2016). Outflow is shared
  among lower neighbours by volume weights and limited three ways: by the cell's water, by a velocity limiter
  using `min(sqrt(g d), Manning velocity)`, and by the previous outflow. An extra cap, `max_frac = 0.2`, was
  added for stability. Level differences below `tau = 1e-3` m are ignored.
- **Step.** `dt = alpha dx / max(v)`, with `alpha = 0.5`, clamped to `[dt_min, dt_max] = [1e-3, 1.0]`.
- **Cannot do.** Inertia, waves, bores or momentum. It leaves a residual surface slope of up to `tau` per
  cell.
- **Caveat.** It was reconstructed from the paper's description and not checked against the CADDIES
  reference code, so the details may differ.

### `pipes`, `pipes_theta`: virtual pipes (`pipes.py`)

- **Physics.** The graphics workhorse (O'Brien & Hodgins 1995; Mei et al. 2007). Each cell has four outflow
  pipes with memory. A pipe accelerates with the surface drop, `f += dt A g d_eta / dx`, with `A = dx^2`, and
  is clamped at zero. All outflows are scaled so a cell never gives more than it holds. The memory gives
  sloshing and waves. Optional semi-implicit Manning friction damps it.
- **Step.** `alpha sqrt(dx / g)`, with `alpha = 0.4`, whatever the depth.
- **Cannot do.** The pipe cross-section does not depend on depth, so a 1 cm film and a 3 m channel accelerate
  alike. Thin sheets run far too fast, and the wave speed is `sqrt(g dx)`, not `sqrt(g h)`. It is not a
  discretisation of the shallow-water equations, so its timing and bores are not physical.
- **`pipes_theta`.** Adds de Almeida theta-0.8 smoothing of each pipe memory along its axis, and applies the
  dry-face gate even when `manning = 0`. It is visual only. It removes the checkerboard that plain pipes leave
  on a draining lake (lake roughness 70x smaller, measured), and it also damps swell.

### `local_inertial`, `local_inertial_frcap`, `local_inertial_frcap07` (`local_inertial.py`)

- **Physics.** Local-inertial shallow water, the LISFLOOD-FP "acceleration" scheme (Bates et al. 2010; de
  Almeida et al. 2012). The momentum equation keeps local acceleration, the hydrostatic pressure gradient and
  semi-implicit Manning friction (using the coupled vector |q|), and drops advection. Each face keeps a
  discharge with memory. De Almeida's theta weighting (`theta = 0.8`) damps the checkerboard mode.
- **Step.** `alpha dx / sqrt(g h_max)`, with `alpha = 0.6`. The docstring records that alpha 0.7 with theta
  0.8 let a lake at rest grow round-off into metre-scale sloshing (measured at 64, 128 and 250), while 0.6
  held it to 1e-15.
- **Cannot do.** Fast jets, bores and hydraulic jumps, since there is no momentum advection. It is valid for
  Froude numbers below about 0.5. Wet/dry fronts recede crudely.
- **`local_inertial_frcap`.** Caps every face discharge at `fr_max h_f sqrt(g h_f)`, with
  `fr_max = (2/3)^1.5 = 0.544`, the broad-crested-weir (critical-flow) discharge. Without the cap, flow
  through the dam-break gap reached face Froude 1.36 and breach throughput was 1.48x `swe_hll`'s. The cost is
  that it also caps genuinely fast reaches: the river delivers about 11% less to its outfall than `swe_hll`
  (measured at 128). This is the engine the C# port implements.
- **`local_inertial_frcap07`.** The same cap at Froude 0.7: breach 1.28x `swe_hll`, and more of the river's
  fast reach kept.

### Ripple dampers: `local_inertial_frcap_*` (`local_inertial_ripples.py`)

`local_inertial_frcap` leaves criss-cross stripes in the dam-break reservoir: wave fronts reflecting between
the walls, which nothing in the scheme damps. Each class here is `local_inertial_frcap` plus one damper. The
damper acts on the face discharge before the cap and before the transfer, so volume stays exact.

| name | damper | known cost (from the docstring) |
|---|---|---|
| `local_inertial_frcap_divdamp` | divergence damping, `dq/dt -= nu d(dh/dt)/dx`, `nu = kappa dx sqrt(g h_f)`, kappa 0.04; zero where the surface is not changing | slows every changing surface; kappa 0.25 went unstable |
| `local_inertial_frcap_divdamp_soft` | the same at kappa 0.02 | |
| `local_inertial_frcap_visc` | eddy viscosity `nu lap(q)` over wet faces, free-slip at banks, kappa 0.03 | smooths every shear; long swell loses a little |
| `local_inertial_frcap_hyper` | hyperviscosity `q -= beta lap(lap(q))` per step, beta 0.0025 | damping per second grows as dt shrinks; beta must stay under 1/32 |
| `local_inertial_frcap_deeptheta` | theta lowered toward 0.7 only on deep (over 0.5 m), slow (Froude under 0.2) faces | along-axis only, so diagonal fronts are damped less |
| `local_inertial_frcap_sponge` | discharge decays at 0.5/s near walls and map edges (4 cells) | walls reflect less than real walls do |
| `local_inertial_frcap_etadiff` | mass-only surface diffusion `-nu d(eta)/dx`, kappa 0.04, not remembered as discharge; zero on a lake at rest | a little extra downhill conveyance |
| `local_inertial_frcap_smoothdiff` | `etadiff` on a 1-2-1-smoothed surface difference, so kappa 0.3 stays stable | filled the dug pit 15% early at 128 |
| `local_inertial_frcap_smoothdiff_fr` | `smoothdiff` faded out in fast flow, `max(0, 1 - Fr/0.2)`, kappa 0.4 | a slow but steep sheet flow still gets it |
| `local_inertial_frcap_smoothdiff_fr2` / `_fr4` | `smoothdiff_fr` applied every 2nd / 4th step at 2x / 4x strength | `_fr4` gives the same scene numbers as every step (measured) at a quarter of the extra cost |

The latest commit message names `smoothdiff_fr` as the damper that settles the reservoir. Per the
`local_inertial_ripsense` docstring, `_fr4` passes the same bars as the ripple detector at about the same
cost.

### `local_inertial_frcap_ripsense` (`local_inertial_ripsense.py`, not registered)

`local_inertial_frcap` plus a ripple detector, with mass-only diffusion of the smoothed surface applied only
where the detector fires. The treatment weight is blurred in space and eased in time, because hard patch edges
made new ripples.

The docstring reports that the detector is not a clean ripple test: it also fires on swell and slow river
water, and treated 50% (dam break), 96% (river) and 84% (swell) of wet cell-steps. The module is tested
(`tests/test_ripsense.py`) but not listed in `solvers/__init__.py`, so the bakeoff and `compare` do not run
it. Its docstring also cites a result file by a path on the author's machine, which is not in this
repository.

### `kinematic_wave` (`kinematic.py`)

- **Physics.** Manning flow driven by the bed slope only: `q = h^(5/3) sqrt(S0) / n`, shared among downhill
  neighbours by slope, using the donor's depth.
- **Step.** CFL on the kinematic celerity `(5/3) v`, with `alpha = 0.5` and `dt_max = 10`.
- **Cannot do.** Backwater, flow up a bed slope, or a lake at rest on an uneven bed: water runs to local bed
  minima and sits as puddles whose surfaces differ. There is no flow at all on a flat bed. It is good only for
  steep, fast-draining overland flow. The tests mark its lake-at-rest failure as expected, and the bakeoff
  report leaves it out of the water headline. Its `maxdt` value is meaningless, since the dam-break bed is
  flat.

### `diffusion_wave_explicit` (`diffusion_explicit.py`)

- **Physics.** Zero-inertia (diffusion) wave: `q = -(h_f^(5/3)/n) S / sqrt(max(|S|, S_min))`, with S the
  surface slope on the face (the coupled vector magnitude). It handles backwater and a lake at rest.
- **Step.** `dt <= 0.9 dx^2 / (4 D_max)`, with `D = h_f^(5/3) / (n sqrt(max(|S|, S_min)))`, capped at 10 s.
- **Cannot do.** Waves, inertia, bores or supercritical flow. Its known weakness is that D grows without bound
  as the slope flattens, so the explicit step collapses on flat water. With `S_min = 1e-4`, steps are below a
  millisecond, and in the bakeoff it hits the wall-time cap on every water scene at 325. `S_min` is a
  stability knob, not only a regulariser.

### `diffusion_wave_implicit_cg`, `diffusion_wave_implicit_adi` (`diffusion_implicit.py`)

- **Physics.** The same zero-inertia wave, advanced by backward Euler with the face conductance lagged at the
  old state (`picard = 2` fixed-point passes). A conservative corrector then moves `K (eta_i - eta_j) dt`
  through `transfer`; h is never assigned from eta.
  - CG: the sparse SPD system is solved with scipy conjugate gradient and a Jacobi preconditioner, with
    `rtol = 1e-13`.
  - ADI: an implicit x sweep, then an implicit y sweep (Lie splitting), each a batched Thomas solve.
- **Step.** `stable_dt()` returns `max_dt = 5 s`. That is a large stable step, not an accurate one.
- **Cannot do.** Waves, inertia or bores. Backward Euler smears transient timing, and wet/dry fronts are
  crude. ADI adds an O(dt) splitting error and a slight anisotropy on diagonal flow at large dt.

### `swe_hll`: full shallow water, the reference (`swe_hll.py`)

- **Physics.** The 2-D Saint-Venant equations with momentum advection, as first-order finite volume: HLL
  fluxes after Audusse hydrostatic reconstruction, which makes it well-balanced (a lake at rest stays at rest
  on any bed). It is unsplit, so it is mirror-symmetric. Manning friction is semi-implicit, and walls reflect.
  It carries bores, hydraulic jumps, jets and dam-break fronts.
- **Step.** `cfl dx / max(|u| + sqrt(g h))`, with `cfl = 0.45`: the tightest step of the water solvers.
- **Cannot do.** It is first order, so it smears fronts and damps waves over a few dozen cells. Shear layers
  are diffusive, because HLL has no contact wave. Thin films carry noisy velocity.
- **Checked against.** Stoker's analytic wet-bed dam break (`tests/test_water_main.py`) and ANUGA
  (`results/crosscheck.md`).

### `swe_lbm`: shallow-water lattice Boltzmann (`swe_lbm.py`)

- **Physics.** Zhou's (2004) D2Q9 shallow-water LBM with BGK collision (`tau = 0.8`). The bed slope enters as
  a centred force, weighted so a lake at rest is exact. Manning friction is a semi-implicit force. The depth
  each population carries across a face is summed per face and moved through `transfer`, so volume is exact.
  Walls, blocked shores and the map edge use half-way bounce-back.
- **Step.** The lattice speed `e = dx/dt` is kept at least `2.5 (sqrt(g h) + |u|)_max`, with hysteresis. A
  change of `e` resets the populations to equilibrium.
- **Cannot do.** It is a low-Froude scheme: bores and jets near Froude 1 smear and ring, and there is no shock
  capturing. Dry fronts are crude. There is numerical viscosity from tau, and each step costs several
  local-inertial steps.
- **Checked against.** Stoker's dam break (`tests/test_swe_lbm.py`).

## Wave

### `linear_wave` (`linear_wave.py`)

- **Physics.** `u_tt = c^2 lap(u) - damping u_t`, leapfrog (Stormer-Verlet): small-amplitude ripples.
  Obstacles are clamped to 0, and the outer edges are Neumann (reflecting without a phase flip). Scene drops
  add a Gaussian bump at their time.
- **Step.** `cfl dx / (c_max sqrt 2)`, with `cfl = 0.9`.
- **Cannot do.** Nonlinear steepening, breaking, dispersion, wet/dry fronts or mass transport. The field is
  not conserved.

## Scalar

### `heat_explicit`, `heat_adi` (`heat.py`)

- **Physics.** `u_t = D lap(u)` with constant D and insulated walls, in finite-volume form. The total is
  conserved to round-off.
- **Step.**
  - `heat_explicit` (FTCS): `0.9 dx^2 / (4 D)`. Past that it blows up as a growing checkerboard.
  - `heat_adi` (Peaceman-Rachford): second order and unconditionally stable, with `stable_dt = inf`.
    `dt_hint` is 10x the explicit limit.
- **Cannot do.** Variable D or sources. ADI rings and can undershoot below zero on sharp fronts at large
  `D dt / dx^2`.

## Flow

### `stable_fluids` (`stable_fluids.py`)

- **Physics.** Stam (1999) on a staggered MAC grid in cell units. Advection is a semi-Lagrangian RK2 bilinear
  backtrace. The pressure projection is a masked conjugate-gradient Poisson solve. Options: implicit viscosity
  `nu`, vorticity confinement, buoyancy on the dye, and an inflow inlet/outflow mode.
- **Step.** A fixed `dt = 1` (one step).
- **Cannot do.** It is strongly numerically dissipative, so vortex streets at Re about 150 are weak, and there
  is no real Reynolds number. It has no compressibility, no turbulence model, and conserves nothing exactly
  (`conserves = False`).

### `lbm_d2q9` (`lbm_d2q9.py`)

- **Physics.** D2Q9 lattice Boltzmann with BGK (default) or TRT collision, for weakly compressible flow at
  Mach well below 1. Viscosity is `(tau - 1/2)/3` in lattice units, and full-way bounce-back is used at solids.
  Options: a Guo body force, buoyancy on a passive dye, and an inflow mode (equilibrium inlet, zero-gradient
  outlet). tau is set from Re and floored at 0.51, so the effective Re can be lower than requested
  (`effective_Re`). TRT blew up at the inlet corner on `cylinder_wake`, so it is not the default.
- **Step.** `stable_dt = 1`, one lattice step. `substeps` runs several lattice steps per `step`.
- **Cannot do.** High-Mach flow (keep |u| under about 0.15 lattice units), free surfaces or thermal flow.
  Curved walls are staircase (first order), and pressure at the inlet and outlet is not accurate.
