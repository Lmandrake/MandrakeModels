# Fluid solver bakeoff: results

Tags: every number here is [measured] by this harness. numpy ms/tick = numpy ms/step x steps per tick at the solver's own dt, with one simulated second per 60 ticks (sibling plan's assumption). No C# cost is projected from these: the numba probe's kernels are simplified solvers, so its ratio is not a C# projection; measured C# cost (local_inertial on Unity Mono) is in results/stagec.md. Full-grid cost: no sleeping or active sets.

### numba unit cost probe (simplified kernels: NOT the bakeoff's arithmetic, not a C# projection)

| kernel | numpy ms/step | numba ms/step | ratio | numba ns/cell |
|---|---|---|---|---|
| local_inertial | 15.52 | 2.91 | 5.3 | 18.2 |
| pipes | 15.24 | 0.39 | 39.4 | 2.4 |
| heat_explicit | 0.66 | 0.07 | 8.9 | 0.5 |

Grid 400^2, median of 30 steps, single thread. The pipes kernel has no friction and the local_inertial kernel has no theta weighting or coupled-|q| friction (probe_numba.py), so these ratios compare different solvers.

### Water solvers at 325 x 325 (median over water scenes)

| solver | ms/step numpy | dt s | numpy ms/tick | worst vol err | lake drift m | settled | blow-ups | errors | capped |
|---|---|---|---|---|---|---|---|---|---|
| ca_game | 3.58 | 0.5 | 0.119 | 3.1e-15 | 4.4e-16 | 1/9 | 0 | 0 | 0 |
| diffusion_wave_explicit | 9.25 | 4.25e-05 | 3.63e+03 | 7.9e-16 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| diffusion_wave_implicit_adi | 905 | 5 | 3.02 | 6.9e-16 | 5.9e-11 | 1/9 | 0 | 0 | 0 |
| diffusion_wave_implicit_cg | 2.77e+03 | 5 | 9.22 | 7.3e-16 | 2.2e-16 | 1/9 | 0 | 0 | 6 |
| local_inertial | 8.27 | 0.192 | 0.744 | 2.6e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| local_inertial_frcap | 8.34 | 0.192 | 0.745 | 6.6e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| local_inertial_frcap07 | 8.02 | 0.192 | 0.747 | 3.0e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| pipes | 9.46 | 0.152 | 1.04 | 6.7e-15 | 3.6e-15 | 1/9 | 0 | 0 | 0 |
| pipes_theta | 10.3 | 0.152 | 1.14 | 6.7e-15 | 4.4e-16 | 1/9 | 0 | 0 | 0 |
| swe_hll | 18.4 | 0.135 | 2.53 | 5.9e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| swe_lbm | 56 | 0.0965 | 14.6 | 4.5e-15 | 4.4e-16 | 1/9 | 0 | 0 | 7 |
| wca2d | 15.6 | 0.251 | 1.62 | 3.3e-15 | 2.2e-16 | 1/9 | 0 | 0 | 0 |

Left out: kinematic_wave (bed-slope only: no flat-bed flow, no backwater; blows up under level-setters).

### Water solvers at 500 x 500 (median over water scenes)

| solver | ms/step numpy | dt s | numpy ms/tick | worst vol err | lake drift m | settled | blow-ups | errors | capped |
|---|---|---|---|---|---|---|---|---|---|
| ca_game | 9.73 | 0.5 | 0.324 | 1.6e-15 | 4.4e-16 | 1/9 | 0 | 0 | 0 |
| diffusion_wave_explicit | 24.9 | 4.25e-05 | 9.47e+03 | 1.8e-15 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| diffusion_wave_implicit_adi | 2.48e+03 | 5 | 8.27 | 5.0e-16 | 3.0e-11 | 1/9 | 0 | 0 | 4 |
| diffusion_wave_implicit_cg | 9.64e+03 | 5 | 32.1 | 2.6e-16 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| local_inertial | 20.5 | 0.192 | 1.65 | 2.5e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| local_inertial_frcap | 20.4 | 0.192 | 1.64 | 4.7e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| local_inertial_frcap07 | 20.2 | 0.192 | 1.68 | 6.2e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| pipes | 23 | 0.152 | 2.53 | 6.7e-15 | 2.7e-15 | 1/9 | 0 | 0 | 0 |
| pipes_theta | 26.1 | 0.152 | 2.86 | 6.7e-15 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| swe_hll | 48.6 | 0.135 | 5.86 | 1.2e-15 | 8.9e-16 | 1/9 | 0 | 0 | 4 |
| swe_lbm | 140 | 0.118 | 19.8 | 1.6e-15 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| wca2d | 38.3 | 0.254 | 2.55 | 1.8e-15 | 4.4e-16 | 1/9 | 0 | 0 | 1 |

Left out: kinematic_wave (bed-slope only: no flat-bed flow, no backwater; blows up under level-setters).

### numpy ms/step by size (median over scenes; timing phase, sequential, one thread)

| solver | kind | 64 | 128 | 250 | 325 | 400 | 500 |
|---|---|---|---|---|---|---|---|
| lbm_d2q9 | flow | 3.22 | 13.1 | 51.8 | 141 | 143 | 224 |
| stable_fluids | flow | 7.22 | 25.6 | 87.4 | 245 | 246 | 481 |
| heat_adi | scalar | 0.253 | 0.709 | 1.35 | 2.34 | 2.56 | 3.87 |
| heat_explicit | scalar | 0.0224 | 0.0775 | 0.241 | 0.521 | 0.685 | 1.23 |
| ca_game | water | 0.133 | 0.389 | 1.78 | 3.58 | 5.24 | 9.73 |
| diffusion_wave_explicit | water | 0.375 | 1.31 | 5.09 | 9.25 | 13.8 | 24.9 |
| diffusion_wave_implicit_adi | water | 1.66 | 9.01 | 415 | 905 | 1.45e+03 | 2.48e+03 |
| diffusion_wave_implicit_cg | water | 30.4 | 141 | 852 | 2.77e+03 | 6.53e+03 | 9.64e+03 |
| kinematic_wave | water | 0.26 | 0.85 | 3.96 | 8.05 | 11.6 | 21 |
| local_inertial | water | 0.324 | 0.941 | 4.37 | 8.27 | 11.4 | 20.5 |
| local_inertial_frcap | water | 0.345 | 1.07 | 4.25 | 8.34 | 11.7 | 20.4 |
| local_inertial_frcap07 | water | 0.376 | 1.11 | 4.15 | 8.02 | 11.6 | 20.2 |
| pipes | water | 0.327 | 1.08 | 4.45 | 9.46 | 12.8 | 23 |
| pipes_theta | water | 0.427 | 1.31 | 5.22 | 10.3 | 14.5 | 26.1 |
| swe_hll | water | 0.657 | 2.07 | 9.12 | 18.4 | 27 | 48.6 |
| swe_lbm | water | 2.02 | 6.47 | 29 | 56 | 83.7 | 140 |
| wca2d | water | 0.47 | 2.31 | 7.62 | 15.6 | 21.7 | 38.3 |
| linear_wave | wave | 0.0667 | 0.234 | 0.762 | 3.95 | 3.47 | 4.58 |

### Largest fixed dt that survives 20 s of the 128 dam break (water) / gaussian / ripples

Survives means no NaN, no blow-up and depth never above 1.25x its start (scalar: no new extreme; wave: bounded) over the whole run. Bounded is not accurate: a limiter can keep a solver bounded at a dt where its wave speeds are already wrong.

| solver | max dt s | solver's own stable_dt s | note |
|---|---|---|---|
| ca_game |  | 0.5 | n/a: fixed pass |
| diffusion_wave_explicit | 0.00403 | 4.25e-05 |  |
| diffusion_wave_implicit_adi | 64 | 5 | >= 64 (cap) |
| diffusion_wave_implicit_cg | 64 | 5 | >= 64 (cap) |
| heat_adi | 64 |  | >= 64 (cap) |
| heat_explicit | 0.391 | 0.225 |  |
| kinematic_wave | 64 | 10 | meaningless here: the dam-break bed is flat, so kinematic flow never starts |
| linear_wave | 0.0709 | 0.0636 |  |
| local_inertial | 0.207 | 0.192 |  |
| local_inertial_frcap | 0.207 | 0.192 |  |
| local_inertial_frcap07 | 0.207 | 0.192 |  |
| pipes | 0.197 | 0.152 |  |
| pipes_theta | 0.197 | 0.152 |  |
| swe_hll | 0.171 | 0.144 |  |
| swe_lbm | 0.294 | 0.128 |  |
| wca2d | 64 | 0.16 | >= 64 (cap) |

### Flood timing: water through the dam breach in 20 s (dam_break, 128^2)

Reference swe_hll (matches Stoker's analytic dam break to 0.1%, tests/test_water_main.py). Ratio well below 1: floods arrive late; above 1: early.

| solver | m^3 through breach | ratio to reference |
|---|---|---|
| diffusion_wave_implicit_cg | 7049 | 2.73 |
| pipes | 4113 | 1.59 |
| pipes_theta | 4023 | 1.56 |
| diffusion_wave_implicit_adi | 4007 | 1.55 |
| local_inertial | 3818 | 1.48 |
| local_inertial_frcap07 | 3294 | 1.28 |
| local_inertial_frcap | 2820 | 1.09 |
| swe_hll | 2583 | 1.00 |
| swe_lbm | 1595 | 0.62 |
| wca2d | 230 | 0.09 |
| ca_game | 89 | 0.03 |
| kinematic_wave | 0 | 0.00 |

### Every run at 250 x 250

| scene | solver | ms/step | dt s | steps | sim s | settled at s | final max dh/dt m/s | vol err | min h | blew up | capped | error |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cylinder_wake | lbm_d2q9 | 14.3 | 1 | 2000 | 2e+03 |  |  |  |  | no | no |  |
| cylinder_wake | stable_fluids | 63.5 | 1 | 1350 | 1.35e+03 |  |  |  |  | no | yes |  |
| dam_break | ca_game | 1.87 | 0.5 | 80 | 40 |  | 5.2e-03 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | diffusion_wave_explicit | 5.34 | 4.25e-05 | 10468 | 0.445 |  | 1.6e+00 | 1.9e-15 | 0.0e+00 | no | yes |  |
| dam_break | diffusion_wave_implicit_adi | 26.1 | 5 | 8 | 40 |  | 1.0e-02 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | diffusion_wave_implicit_cg | 852 | 5 | 8 | 40 |  | 1.1e-02 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | kinematic_wave | 4.27 | 5 | 8 | 40 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | local_inertial | 4.69 | 0.192 | 209 | 40 |  | 4.4e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | local_inertial_frcap | 4.22 | 0.192 | 209 | 40 |  | 3.0e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | local_inertial_frcap07 | 4.23 | 0.192 | 209 | 40 |  | 4.1e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | pipes | 5.36 | 0.152 | 264 | 40 |  | 3.0e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | pipes_theta | 6.21 | 0.152 | 264 | 40 |  | 1.0e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | swe_hll | 9.1 | 0.109 | 395 | 40 |  | 1.5e-01 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | swe_lbm | 28.9 | 0.0965 | 1527 | 19.8 |  | 1.5e-01 | 1.4e-16 | 0.0e+00 | no | yes |  |
| dam_break | wca2d | 8.11 | 0.16 | 250 | 40 |  | 5.4e-02 | 0.0e+00 | 0.0e+00 | no | no |  |
| dye | heat_adi | 1.23 | 2.25 | 50 | 112 |  |  |  |  | no | no |  |
| dye | heat_explicit | 0.247 | 0.225 | 500 | 112 |  |  |  |  | no | no |  |
| gaussian | heat_adi | 1.47 | 6.25 | 50 | 312 |  |  |  |  | no | no |  |
| gaussian | heat_explicit | 0.235 | 0.225 | 1389 | 312 |  |  |  |  | no | no |  |
| lake_at_rest | ca_game | 1.78 | 0.5 | 120 | 60 | 0 | 1.8e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | diffusion_wave_explicit | 5.18 | 2.4e-05 | 10878 | 0.261 |  | 1.8e-11 | 0.0e+00 | 0.0e+00 | no | yes |  |
| lake_at_rest | diffusion_wave_implicit_adi | 13.7 | 5 | 12 | 60 | 0 | 1.3e-11 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | diffusion_wave_implicit_cg | 660 | 5 | 12 | 60 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | kinematic_wave | 4.02 | 0.0254 | 1245 | 60 |  | 6.9e+00 | 4.9e-16 | 0.0e+00 | no | no |  |
| lake_at_rest | local_inertial | 4.78 | 0.161 | 372 | 60 | 0 | 6.2e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | local_inertial_frcap | 3.97 | 0.161 | 372 | 60 | 0 | 6.2e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | local_inertial_frcap07 | 3.94 | 0.161 | 372 | 60 | 0 | 6.2e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | pipes | 4.45 | 0.152 | 396 | 60 | 0 | 3.1e-14 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | pipes_theta | 4.97 | 0.152 | 396 | 60 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | swe_hll | 12.3 | 0.121 | 496 | 60 | 0 | 8.6e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | swe_lbm | 30.7 | 0.108 | 558 | 60 | 0 | 4.8e-18 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | wca2d | 8.02 | 1 | 60 | 60 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_levelset | ca_game | 1.75 | 0.5 | 240 | 120 |  | 1.2e-03 | 6.7e-18 | 0.0e+00 | no | no |  |
| lake_levelset | diffusion_wave_explicit | 5.03 | 4.88e-05 | 9891 | 0.483 |  | 6.1e-01 | 3.3e-16 | 0.0e+00 | no | yes |  |
| lake_levelset | diffusion_wave_implicit_adi | 415 | 5 | 24 | 120 |  | 8.1e-04 | 1.5e-17 | 0.0e+00 | no | no |  |
| lake_levelset | diffusion_wave_implicit_cg | 957 | 5 | 24 | 120 |  | 2.9e-04 | 3.1e-16 | 0.0e+00 | no | no |  |
| lake_levelset | kinematic_wave | 3.91 | 0.0459 | 1493 | 120 |  | 4.6e+00 | 1.2e-16 | 0.0e+00 | no | no |  |
| lake_levelset | local_inertial | 4.07 | 0.2 | 599 | 120 |  | 6.5e-03 | 2.0e-16 | 0.0e+00 | no | no |  |
| lake_levelset | local_inertial_frcap | 4.25 | 0.2 | 599 | 120 |  | 3.4e-03 | 2.7e-16 | 0.0e+00 | no | no |  |
| lake_levelset | local_inertial_frcap07 | 4.06 | 0.2 | 599 | 120 |  | 3.8e-03 | 1.3e-16 | 0.0e+00 | no | no |  |
| lake_levelset | pipes | 4.41 | 0.152 | 791 | 120 |  | 2.0e-02 | 1.5e-16 | 0.0e+00 | no | no |  |
| lake_levelset | pipes_theta | 5.22 | 0.152 | 791 | 120 |  | 3.0e-02 | 7.3e-17 | 0.0e+00 | no | no |  |
| lake_levelset | swe_hll | 8.77 | 0.15 | 802 | 120 |  | 3.5e-03 | 5.3e-17 | 0.0e+00 | no | no |  |
| lake_levelset | swe_lbm | 29.9 | 0.0264 | 1446 | 18.9 |  | 1.2e-02 | 1.9e-16 | 0.0e+00 | no | yes |  |
| lake_levelset | wca2d | 7.44 | 0.291 | 390 | 120 |  | 3.6e-03 | 8.3e-18 | 0.0e+00 | no | no |  |
| moat_evap | ca_game | 1.63 | 0.5 | 240 | 120 |  | 1.5e-02 | 5.0e-16 | 0.0e+00 | no | no |  |
| moat_evap | diffusion_wave_explicit | 5.41 | 3.39e-05 | 14244 | 0.484 |  | 1.5e-02 | 2.2e-15 | 0.0e+00 | no | yes |  |
| moat_evap | diffusion_wave_implicit_adi | 12.5 | 5 | 24 | 120 |  | 1.5e-02 | 5.0e-16 | 0.0e+00 | no | no |  |
| moat_evap | diffusion_wave_implicit_cg | 577 | 5 | 24 | 120 |  | 1.5e-02 | 5.0e-16 | 0.0e+00 | no | no |  |
| moat_evap | kinematic_wave | 3.77 | 0.0467 | 824 | 120 |  | 1.5e-02 | 4.1e-16 | 0.0e+00 | no | no |  |
| moat_evap | local_inertial | 4.37 | 0.181 | 514 | 120 |  | 1.5e-02 | 9.9e-16 | 0.0e+00 | no | no |  |
| moat_evap | local_inertial_frcap | 4 | 0.181 | 514 | 120 |  | 1.5e-02 | 9.9e-16 | 0.0e+00 | no | no |  |
| moat_evap | local_inertial_frcap07 | 4.15 | 0.181 | 514 | 120 |  | 1.5e-02 | 9.9e-16 | 0.0e+00 | no | no |  |
| moat_evap | pipes | 4.33 | 0.152 | 791 | 120 |  | 1.5e-02 | 2.4e-14 | 0.0e+00 | no | no |  |
| moat_evap | pipes_theta | 5 | 0.152 | 791 | 120 |  | 1.5e-02 | 2.4e-14 | 0.0e+00 | no | no |  |
| moat_evap | swe_hll | 9.63 | 0.135 | 685 | 120 |  | 1.5e-02 | 8.3e-16 | 0.0e+00 | no | no |  |
| moat_evap | swe_lbm | 28.4 | 0.119 | 985 | 120 |  | 1.5e-02 | 1.1e-14 | 0.0e+00 | no | no |  |
| moat_evap | wca2d | 7.09 | 1 | 120 | 120 |  | 1.5e-02 | 3.6e-15 | 0.0e+00 | no | no |  |
| ocean_edge | ca_game | 1.81 | 0.5 | 240 | 120 |  | 2.0e-02 | 1.1e-16 | 0.0e+00 | no | no |  |
| ocean_edge | diffusion_wave_explicit | 5.09 | 5.49e-06 | 10142 | 0.0557 |  | 1.8e-01 | 2.8e-17 | 0.0e+00 | no | yes |  |
| ocean_edge | diffusion_wave_implicit_adi | 440 | 5 | 24 | 120 |  | 1.8e-01 | 2.7e-17 | 0.0e+00 | no | no |  |
| ocean_edge | diffusion_wave_implicit_cg | 801 | 5 | 24 | 120 |  | 2.3e-01 | 4.1e-17 | 0.0e+00 | no | no |  |
| ocean_edge | kinematic_wave | 4.41 | 0.0131 | 160 | 1.83 |  | 5.4e+03 | 1.6e-17 | 0.0e+00 | yes (1.83 s) | no |  |
| ocean_edge | local_inertial | 3.9 | 0.104 | 1158 | 120 |  | 6.5e-02 | 2.8e-17 | 0.0e+00 | no | no |  |
| ocean_edge | local_inertial_frcap | 4.62 | 0.104 | 1158 | 120 |  | 9.1e-02 | 9.8e-17 | 0.0e+00 | no | no |  |
| ocean_edge | local_inertial_frcap07 | 5.19 | 0.104 | 1158 | 120 |  | 4.3e-02 | 2.8e-17 | 0.0e+00 | no | no |  |
| ocean_edge | pipes | 4.59 | 0.152 | 791 | 120 |  | 5.5e-03 | 3.4e-18 | 0.0e+00 | no | no |  |
| ocean_edge | pipes_theta | 5.32 | 0.152 | 791 | 120 |  | 1.8e-02 | 5.6e-17 | 0.0e+00 | no | no |  |
| ocean_edge | swe_hll | 10.1 | 0.0777 | 1544 | 120 |  | 6.5e-03 | 3.0e-17 | 0.0e+00 | no | no |  |
| ocean_edge | swe_lbm | 30.9 | 0.0388 | 1360 | 16.9 |  | 3.3e-01 | 4.8e-17 | 0.0e+00 | no | yes |  |
| ocean_edge | wca2d | 7.72 | 0.251 | 823 | 120 |  | 3.0e+00 | 2.7e-17 | 0.0e+00 | no | no |  |
| pit_dig | ca_game | 1.55 | 0.5 | 300 | 150 |  | 1.7e-02 | 6.6e-16 | 0.0e+00 | no | no |  |
| pit_dig | diffusion_wave_explicit | 4.73 | 0.000223 | 17907 | 2.9 |  | 3.0e-02 | 4.2e-16 | 0.0e+00 | no | yes |  |
| pit_dig | diffusion_wave_implicit_adi | 453 | 5 | 30 | 150 |  | 1.5e-01 | 9.2e-16 | 0.0e+00 | no | no |  |
| pit_dig | diffusion_wave_implicit_cg | 782 | 5 | 30 | 150 |  | 2.2e-01 | 4.4e-16 | 0.0e+00 | no | no |  |
| pit_dig | kinematic_wave | 3.54 | 0.115 | 1152 | 150 |  | 1.1e-03 | 4.4e-15 | 0.0e+00 | no | no |  |
| pit_dig | local_inertial | 3.98 | 0.25 | 570 | 150 |  | 5.7e-02 | 1.8e-16 | 0.0e+00 | no | no |  |
| pit_dig | local_inertial_frcap | 3.58 | 0.25 | 554 | 150 |  | 1.7e-02 | 2.6e-16 | 0.0e+00 | no | no |  |
| pit_dig | local_inertial_frcap07 | 3.88 | 0.25 | 556 | 150 |  | 8.3e-03 | 4.0e-16 | 0.0e+00 | no | no |  |
| pit_dig | pipes | 4.14 | 0.152 | 988 | 150 |  | 2.5e-02 | 6.8e-15 | 0.0e+00 | no | no |  |
| pit_dig | pipes_theta | 4.55 | 0.152 | 988 | 150 |  | 3.6e-03 | 5.5e-15 | 0.0e+00 | no | no |  |
| pit_dig | swe_hll | 8.24 | 0.162 | 1016 | 150 |  | 5.8e-03 | 1.1e-16 | 0.0e+00 | no | no |  |
| pit_dig | swe_lbm | 27.7 | 0.151 | 2815 | 42.3 |  | 7.7e-02 | 1.3e-16 | 0.0e+00 | no | yes |  |
| pit_dig | wca2d | 7.22 | 0.184 | 1222 | 150 |  | 8.0e+00 | 2.9e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | ca_game | 1.79 | 0.5 | 240 | 120 |  | 1.4e-02 | 1.1e-15 | 0.0e+00 | no | no |  |
| rain_groundwater | diffusion_wave_explicit | 5.58 | 0.00217 | 13017 | 10.7 |  | 6.0e-02 | 3.4e-15 | 0.0e+00 | no | yes |  |
| rain_groundwater | diffusion_wave_implicit_adi | 460 | 5 | 24 | 120 |  | 6.3e-02 | 4.2e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | diffusion_wave_implicit_cg | 1.18e+03 | 5 | 24 | 120 |  | 2.8e-02 | 0.0e+00 | 0.0e+00 | no | no |  |
| rain_groundwater | kinematic_wave | 4.3 | 0.165 | 2246 | 120 |  | 7.2e+01 | 2.9e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | local_inertial | 4.37 | 0.36 | 510 | 120 |  | 2.2e-02 | 2.8e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | local_inertial_frcap | 4.57 | 0.361 | 504 | 120 |  | 8.1e-03 | 5.4e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | local_inertial_frcap07 | 4.66 | 0.36 | 507 | 120 |  | 8.6e-03 | 5.5e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | pipes | 4.92 | 0.152 | 791 | 120 |  | 8.9e-03 | 1.7e-15 | 0.0e+00 | no | no |  |
| rain_groundwater | pipes_theta | 5.24 | 0.152 | 791 | 120 |  | 9.8e-03 | 3.1e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | swe_hll | 9.12 | 0.247 | 710 | 120 |  | 7.9e-03 | 5.6e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | swe_lbm | 28.3 | 0.191 | 2171 | 36.6 |  | 1.2e-02 | 7.6e-16 | 0.0e+00 | no | yes |  |
| rain_groundwater | wca2d | 8.38 | 0.415 | 509 | 120 |  | 1.7e-01 | 9.0e-16 | 0.0e+00 | no | no |  |
| ripples | linear_wave | 0.762 | 0.0636 | 315 | 20 |  |  |  |  | no | no |  |
| river | ca_game | 1.83 | 0.5 | 240 | 120 |  | 2.0e-02 | 1.1e-15 | 0.0e+00 | no | no |  |
| river | diffusion_wave_explicit | 4.89 | 0.000199 | 10374 | 2.81 |  | 2.8e-02 | 1.7e-16 | 0.0e+00 | no | yes |  |
| river | diffusion_wave_implicit_adi | 454 | 5 | 24 | 120 |  | 1.1e-02 | 8.4e-16 | 0.0e+00 | no | no |  |
| river | diffusion_wave_implicit_cg | 869 | 5 | 24 | 120 |  | 2.4e-04 | 1.9e-16 | 0.0e+00 | no | no |  |
| river | kinematic_wave | 3.74 | 0.0857 | 992 | 120 |  | 2.7e+00 | 2.3e-16 | 0.0e+00 | no | no |  |
| river | local_inertial | 4.44 | 0.242 | 490 | 120 |  | 2.1e-02 | 7.8e-17 | 0.0e+00 | no | no |  |
| river | local_inertial_frcap | 4.36 | 0.242 | 490 | 120 |  | 2.3e-02 | 3.9e-17 | 0.0e+00 | no | no |  |
| river | local_inertial_frcap07 | 4.18 | 0.242 | 490 | 120 |  | 2.2e-02 | 2.3e-16 | 0.0e+00 | no | no |  |
| river | pipes | 4.36 | 0.152 | 791 | 120 |  | 2.9e-02 | 3.6e-15 | 0.0e+00 | no | no |  |
| river | pipes_theta | 4.84 | 0.152 | 791 | 120 |  | 8.9e-02 | 3.6e-15 | 0.0e+00 | no | no |  |
| river | swe_hll | 8.8 | 0.162 | 780 | 120 |  | 2.0e-03 | 1.9e-16 | 0.0e+00 | no | no |  |
| river | swe_lbm | 29 | 0.145 | 1010 | 120 |  | 1.5e-02 | 1.6e-16 | 0.0e+00 | no | no |  |
| river | wca2d | 7.38 | 0.183 | 973 | 120 |  | 1.3e-01 | 4.7e-16 | 0.0e+00 | no | no |  |
| smoke_plume | lbm_d2q9 | 89.4 | 1 | 141 | 141 |  |  |  |  | yes (141 s) | no |  |
| smoke_plume | stable_fluids | 111 | 1 | 520 | 520 |  |  |  |  | no | yes |  |
| swell | ca_game | 1.56 | 0.5 | 120 | 60 |  | 2.3e-01 | 1.0e-17 | 0.0e+00 | no | no |  |
| swell | diffusion_wave_explicit | 4.99 | 3.82e-06 | 9487 | 0.0362 |  | 2.4e-01 | 3.4e-15 | 0.0e+00 | no | yes |  |
| swell | diffusion_wave_implicit_adi | 13.6 | 5 | 12 | 60 |  | 4.2e-02 | 1.4e-16 | 0.0e+00 | no | no |  |
| swell | diffusion_wave_implicit_cg | 1.09e+03 | 5 | 12 | 60 |  | 4.2e-02 | 1.1e-16 | 0.0e+00 | no | no |  |
| swell | kinematic_wave | 3.96 | 0.0115 | 146 | 1.49 |  | 7.7e+03 | 3.7e-17 | 0.0e+00 | yes (1.49 s) | no |  |
| swell | local_inertial | 3.89 | 0.0922 | 649 | 60 |  | 4.3e-01 | 1.5e-16 | 0.0e+00 | no | no |  |
| swell | local_inertial_frcap | 4.27 | 0.0922 | 649 | 60 |  | 4.3e-01 | 1.5e-16 | 0.0e+00 | no | no |  |
| swell | local_inertial_frcap07 | 4.05 | 0.0922 | 649 | 60 |  | 4.3e-01 | 1.5e-16 | 0.0e+00 | no | no |  |
| swell | pipes | 4.58 | 0.152 | 396 | 60 |  | 4.2e-01 | 2.0e-17 | 0.0e+00 | no | no |  |
| swell | pipes_theta | 5.3 | 0.152 | 396 | 60 |  | 2.4e-01 | 5.2e-17 | 0.0e+00 | no | no |  |
| swell | swe_hll | 9.37 | 0.0683 | 883 | 60 |  | 2.4e-01 | 4.0e-17 | 0.0e+00 | no | no |  |
| swell | swe_lbm | 30.8 | 0.0611 | 1327 | 60 |  | 2.4e-01 | 2.5e-17 | 0.0e+00 | no | no |  |
| swell | wca2d | 7.62 | 0.0768 | 765 | 60 |  | 4.6e-01 | 5.3e-17 | 0.0e+00 | no | no |  |
