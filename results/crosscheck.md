# Cross-check against open-source solvers

Same initial grid (`compare --scene all --n 64 --external`, 2026-11-12), same clock, reference `swe_hll`
(water) / `stable_fluids` (smoke). Every number [measured] on a shared, loaded machine (load 7-34 on
14 cores). rms = RMS of depth (m) or dye against the reference over the frames both reached;
ms/sim-s = wall ms per simulated second (external libraries sub-step internally, so ms/step is not
comparable). Wall cap 10-20 s per solver per scene; capped runs are compared over what they reached.

| scene (64^2) | ANUGA rms | Landlab OverlandFlow rms | Landlab kinwave rms | swe_hll ms/sim-s | ANUGA | OverlandFlow | kinwave |
|---|---|---|---|---|---|---|---|
| dam_break | 0.039 | 0.21 | 0.42 | 6.9 | 90 | 263 | 0.3 |
| lake_at_rest | 0.000 | 0.044 | 3.4 | 5.1 | 76 | 522 | 4471 |
| lake_levelset | 0.004 | 0.040 | 2.7 | 4.2 | 47 | 579 | 2106 |
| river | 0.026 | 0.035 | 0.71 | 5.2 | 43 | 248 | 845 |
| ocean_edge | 0.007 | 0.019 | 4.5 | 10.1 | 101 | 640 | 3871 |
| swell | 0.099 | 0.28 | 3.9 | 9.8 | 141 | 921 | 35 |
| pit_dig | 0.009 | 0.16 | 0.89 | 4.7 | 39 | 457 | 1650 |
| moat_evap | 0.000 | 0.051 | 2.2 | 3.7 | 32 | 375 | 2328 |
| rain_groundwater | 0.006 | 0.079 | 1.0 | 4.2 | 77 | 316 | 425 |

Volume: ANUGA and OverlandFlow conserve to <=1e-6 relative except OverlandFlow on river (9e-4,
inflow/outflow edges). Landlab kinwave loses or gains 10-200% (it routes down the bed slope only, with no
pressure term, and does not honour our closed edges); our own `kinematic_wave` shares the physics
limit (drains lakes, blows up in pits) but conserves.

Dam break, 48^2, water moved past the wall by t=10 s: swe_hll 675.0, ANUGA 676.5 (0.2%), OverlandFlow
606 (90%). ANUGA (Kurganov-type central-upwind finite volume on triangles) is the independent
confirmation of `swe_hll`; OverlandFlow (de Almeida local inertial) shows the same checkerboard on deep
water as our `pipes`.

Smoke: PhiFlow 3.4 agrees at 48^2 (dye centroid row 7.2 vs 7.7 after 40 steps) but not at 64^2: rms
0.17, spurious speeds of 2.2 against 0.2 by step 12 even with CG to 1e-6 (scipy-direct reports
Diverged). Probably our adapter (obstacles passed as an instance-batched Box set, source on the bottom
wall) rather than PhiFlow [guess ~60%]; owed. Its speed: ~200-260 ms/step against ~6-10 for ours.

Not run: Clawpack/PyClaw 5.14 (PyPI has only a macOS-arm64 wheel and an sdist; no Fortran compiler
here); pylbm 0.12 installed, no adapter (geometry from analytic shapes, not cell masks); lbmpy
(GPL-3.0-or-later, not used). Landlab `LinearDiffusionOverlandFlowRouter` runs but its CFL step is
4e-5 s at 2 m depth (built for thin sheet flow), so it is excluded from the all-scene run.
