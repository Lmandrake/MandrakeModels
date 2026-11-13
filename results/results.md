# Fluid solver bakeoff: results

Tags: every number here is [measured] by this harness unless the column says guess. C# ms/tick is a GUESS: numpy ms/step / (numpy/numba ratio from the probe, or the probes' median where marked) x Mono factor 1.5-4.0 [guess], x steps per tick at the solver's own dt, with one simulated second per 60 ticks (sibling plan's assumption). Full-grid cost: no sleeping or active sets, so a real C# port with sleeping cells costs less on mostly dry or settled maps.

### Compiled unit cost probe (numba, single thread)

| kernel | numpy ms/step | numba ms/step | ratio | numba ns/cell |
|---|---|---|---|---|
| local_inertial | 15.52 | 2.91 | 5.3 | 18.2 |
| pipes | 15.24 | 0.39 | 39.4 | 2.4 |
| heat_explicit | 0.66 | 0.07 | 8.9 | 0.5 |

Grid 400^2, median of 30 steps, single thread; numba kernels are plain loops of the same arithmetic.

### Water solvers at 325 x 325 (median over water scenes)

| solver | ms/step numpy | dt s | numpy ms/tick | C# ms/tick guess | fits 0.5 ms/tick | worst vol err | lake drift m | settled | blow-ups | errors | capped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ca_game | 3.73 | 0.5 | 0.124 | 0.021-0.056 | yes | 3.1e-15 | 4.4e-16 | 1/9 | 0 | 0 | 0 |
| diffusion_wave_explicit | 6.49 | 4.25e-05 | 2.76e+03 | 4.6e+02-1.2e+03 | no | 1.4e-15 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| diffusion_wave_implicit_adi | 905 | 5 | 3.02 | 0.51-1.4 | no | 6.9e-16 | 5.9e-11 | 1/9 | 0 | 0 | 0 |
| diffusion_wave_implicit_cg | 2.77e+03 | 5 | 9.22 | 1.5-4.1 | no | 7.3e-16 | 2.2e-16 | 1/9 | 0 | 0 | 6 |
| kinematic_wave | 7.28 | 0.0576 | 2.11 | 0.35-0.94 | maybe | 2.6e-15 | 1.0e+04 | 1/9 | 5 | 0 | 0 |
| local_inertial | 8.24 | 0.192 | 0.731 | 0.21-0.55 | maybe | 2.6e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| pipes | 9.12 | 0.152 | 1 | 0.038-0.1 | yes | 6.7e-15 | 3.6e-15 | 1/9 | 0 | 0 | 0 |
| swe_hll | 23.5 | 0.135 | 2.79 | 0.47-1.3 | maybe | 5.9e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| swe_lbm | 71.8 | 0.0965 | 10.8 | 1.8-4.8 | no | 4.5e-15 | 4.4e-16 | 1/9 | 0 | 0 | 7 |
| wca2d | 12.8 | 0.251 | 0.862 | 0.14-0.39 | yes | 3.3e-15 | 2.2e-16 | 1/9 | 0 | 0 | 0 |

### Water solvers at 500 x 500 (median over water scenes)

| solver | ms/step numpy | dt s | numpy ms/tick | C# ms/tick guess | fits 0.5 ms/tick | worst vol err | lake drift m | settled | blow-ups | errors | capped |
|---|---|---|---|---|---|---|---|---|---|---|---|
| ca_game | 13.5 | 0.5 | 0.45 | 0.076-0.2 | yes | 1.6e-15 | 4.4e-16 | 1/9 | 0 | 0 | 0 |
| diffusion_wave_explicit | 17.4 | 4.25e-05 | 7.49e+03 | 1.3e+03-3.4e+03 | no | 1.8e-15 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| diffusion_wave_implicit_adi | 2.48e+03 | 5 | 8.27 | 1.4-3.7 | no | 4.1e-16 | 3.0e-11 | 1/9 | 0 | 0 | 5 |
| diffusion_wave_implicit_cg | 9.64e+03 | 5 | 32.1 | 5.4-14 | no | 2.6e-16 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| kinematic_wave | 17.9 | 0.0777 | 3.7 | 0.62-1.7 | no | 9.9e-16 | 1.0e+04 | 1/9 | 6 | 0 | 0 |
| local_inertial | 21.9 | 0.192 | 1.98 | 0.56-1.5 | no | 6.3e-16 | 8.9e-16 | 1/9 | 0 | 0 | 0 |
| pipes | 26.8 | 0.152 | 2.94 | 0.11-0.3 | yes | 6.7e-15 | 2.7e-15 | 1/9 | 0 | 0 | 0 |
| swe_hll | 58 | 0.135 | 7.16 | 1.2-3.2 | no | 1.2e-15 | 8.9e-16 | 1/9 | 0 | 0 | 5 |
| swe_lbm | 165 | 0.118 | 24.7 | 4.2-11 | no | 1.6e-15 | 4.4e-16 | 0/9 | 0 | 0 | 9 |
| wca2d | 40.2 | 0.254 | 2.72 | 0.46-1.2 | maybe | 1.8e-15 | 4.4e-16 | 1/9 | 0 | 0 | 1 |

### numpy ms/step by size (median over scenes; timing phase, sequential, one thread)

| solver | kind | 64 | 128 | 250 | 325 | 400 | 500 |
|---|---|---|---|---|---|---|---|
| lbm_d2q9 | flow | 3.22 | 13.1 | 51.8 | 141 | 143 | 224 |
| stable_fluids | flow | 7.22 | 25.6 | 87.4 | 245 | 246 | 481 |
| heat_adi | scalar | 0.253 | 0.709 | 1.35 | 2.34 | 2.56 | 3.87 |
| heat_explicit | scalar | 0.0224 | 0.0775 | 0.241 | 0.521 | 0.685 | 1.23 |
| ca_game | water | 0.125 | 0.382 | 1.84 | 3.73 | 8.79 | 13.5 |
| diffusion_wave_explicit | water | 0.235 | 0.851 | 3.38 | 6.49 | 14.4 | 17.4 |
| diffusion_wave_implicit_adi | water | 1.66 | 9.01 | 415 | 905 | 1.45e+03 | 2.48e+03 |
| diffusion_wave_implicit_cg | water | 30.4 | 141 | 852 | 2.77e+03 | 6.53e+03 | 9.64e+03 |
| kinematic_wave | water | 0.18 | 0.635 | 3.05 | 7.28 | 14.5 | 17.9 |
| local_inertial | water | 0.307 | 1.1 | 4.6 | 8.24 | 19.2 | 21.9 |
| pipes | water | 0.317 | 1.16 | 4.93 | 9.12 | 21.1 | 26.8 |
| swe_hll | water | 0.546 | 2.25 | 9.76 | 23.5 | 40.9 | 58 |
| swe_lbm | water | 1.95 | 6.85 | 32.2 | 71.8 | 116 | 165 |
| wca2d | water | 0.368 | 1.86 | 7.15 | 12.8 | 26.9 | 40.2 |
| linear_wave | wave | 0.0667 | 0.234 | 0.762 | 3.95 | 3.47 | 4.58 |

### Largest fixed dt that survives 20 s of the 128 dam break (water) / gaussian / ripples

Survives means no NaN, no blow-up and depth never above 1.25x its start (scalar: no new extreme; wave: bounded) over the whole run. Bounded is not accurate: a limiter can keep a solver bounded at a dt where its wave speeds are already wrong.

| solver | max dt s | solver's own stable_dt s | note |
|---|---|---|---|
| ca_game |  | 0.5 | n/a: fixed pass |
| diffusion_wave_explicit | 0.00373 | 4.25e-05 |  |
| diffusion_wave_implicit_adi | 64 | 5 | >= 64 (cap) |
| diffusion_wave_implicit_cg | 64 | 5 | >= 64 (cap) |
| heat_adi | 56.2 |  |  |
| heat_explicit | 0.387 | 0.225 |  |
| kinematic_wave | 64 | 10 | meaningless here: the dam-break bed is flat, so kinematic flow never starts |
| linear_wave | 0.0709 | 0.0636 |  |
| local_inertial | 0.207 | 0.192 |  |
| pipes | 0.197 | 0.152 |  |
| swe_hll | 0.171 | 0.144 |  |
| swe_lbm | 0.294 | 0.128 |  |
| wca2d | 81.7 | 0.16 | >= 81.7222 (cap) |

### Flood timing: water through the dam breach in 20 s (dam_break, 128^2)

Reference swe_hll (matches Stoker's analytic dam break to 0.1%, tests/test_water_main.py). Ratio well below 1: floods arrive late; above 1: early.

| solver | m^3 through breach | ratio to reference |
|---|---|---|
| diffusion_wave_implicit_cg | 7049 | 2.73 |
| pipes | 4113 | 1.59 |
| diffusion_wave_implicit_adi | 4007 | 1.55 |
| local_inertial | 3818 | 1.48 |
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
| dam_break | ca_game | 1.68 | 0.5 | 80 | 40 |  | 5.2e-03 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | diffusion_wave_explicit | 3.28 | 4.25e-05 | 10468 | 0.445 |  | 1.6e+00 | 1.9e-15 | 0.0e+00 | no | yes |  |
| dam_break | diffusion_wave_implicit_adi | 26.1 | 5 | 8 | 40 |  | 1.0e-02 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | diffusion_wave_implicit_cg | 852 | 5 | 8 | 40 |  | 1.1e-02 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | kinematic_wave | 3.07 | 5 | 8 | 40 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | local_inertial | 4.62 | 0.192 | 209 | 40 |  | 4.4e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | pipes | 5 | 0.152 | 264 | 40 |  | 3.0e-01 | 1.4e-16 | 0.0e+00 | no | no |  |
| dam_break | swe_hll | 8.99 | 0.109 | 395 | 40 |  | 1.5e-01 | 0.0e+00 | 0.0e+00 | no | no |  |
| dam_break | swe_lbm | 30 | 0.0965 | 1527 | 19.8 |  | 1.5e-01 | 1.4e-16 | 0.0e+00 | no | yes |  |
| dam_break | wca2d | 7.45 | 0.16 | 250 | 40 |  | 5.4e-02 | 0.0e+00 | 0.0e+00 | no | no |  |
| dye | heat_adi | 1.23 | 2.25 | 50 | 112 |  |  |  |  | no | no |  |
| dye | heat_explicit | 0.247 | 0.225 | 500 | 112 |  |  |  |  | no | no |  |
| gaussian | heat_adi | 1.47 | 6.25 | 50 | 312 |  |  |  |  | no | no |  |
| gaussian | heat_explicit | 0.235 | 0.225 | 1389 | 312 |  |  |  |  | no | no |  |
| lake_at_rest | ca_game | 1.85 | 0.5 | 120 | 60 | 0 | 1.8e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | diffusion_wave_explicit | 4.3 | 2.4e-05 | 10878 | 0.261 |  | 1.8e-11 | 0.0e+00 | 0.0e+00 | no | yes |  |
| lake_at_rest | diffusion_wave_implicit_adi | 13.7 | 5 | 12 | 60 | 0 | 1.3e-11 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | diffusion_wave_implicit_cg | 660 | 5 | 12 | 60 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | kinematic_wave | 3.05 | 0.0254 | 1245 | 60 |  | 6.9e+00 | 4.9e-16 | 0.0e+00 | no | no |  |
| lake_at_rest | local_inertial | 5.18 | 0.161 | 372 | 60 | 0 | 6.2e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | pipes | 5.26 | 0.152 | 396 | 60 | 0 | 3.1e-14 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | swe_hll | 10.2 | 0.121 | 496 | 60 | 0 | 8.6e-15 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | swe_lbm | 32.2 | 0.108 | 558 | 60 | 0 | 4.8e-18 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_at_rest | wca2d | 6.46 | 1 | 60 | 60 | 0 | 0.0e+00 | 0.0e+00 | 0.0e+00 | no | no |  |
| lake_levelset | ca_game | 1.84 | 0.5 | 240 | 120 |  | 1.2e-03 | 6.7e-18 | 0.0e+00 | no | no |  |
| lake_levelset | diffusion_wave_explicit | 3.38 | 4.88e-05 | 9891 | 0.483 |  | 6.1e-01 | 3.3e-16 | 0.0e+00 | no | yes |  |
| lake_levelset | diffusion_wave_implicit_adi | 415 | 5 | 24 | 120 |  | 8.1e-04 | 1.5e-17 | 0.0e+00 | no | no |  |
| lake_levelset | diffusion_wave_implicit_cg | 957 | 5 | 24 | 120 |  | 2.9e-04 | 3.1e-16 | 0.0e+00 | no | no |  |
| lake_levelset | kinematic_wave | 3 | 0.0459 | 1493 | 120 |  | 4.6e+00 | 1.2e-16 | 0.0e+00 | no | no |  |
| lake_levelset | local_inertial | 3.92 | 0.2 | 599 | 120 |  | 6.5e-03 | 2.0e-16 | 0.0e+00 | no | no |  |
| lake_levelset | pipes | 4.92 | 0.152 | 791 | 120 |  | 2.0e-02 | 1.5e-16 | 0.0e+00 | no | no |  |
| lake_levelset | swe_hll | 9.1 | 0.15 | 802 | 120 |  | 3.5e-03 | 5.3e-17 | 0.0e+00 | no | no |  |
| lake_levelset | swe_lbm | 34.1 | 0.0264 | 1446 | 18.9 |  | 1.2e-02 | 1.9e-16 | 0.0e+00 | no | yes |  |
| lake_levelset | wca2d | 7.15 | 0.291 | 390 | 120 |  | 3.6e-03 | 8.3e-18 | 0.0e+00 | no | no |  |
| moat_evap | ca_game | 2.75 | 0.5 | 240 | 120 |  | 1.5e-02 | 5.0e-16 | 0.0e+00 | no | no |  |
| moat_evap | diffusion_wave_explicit | 4.04 | 3.39e-05 | 14244 | 0.484 |  | 1.5e-02 | 2.2e-15 | 0.0e+00 | no | yes |  |
| moat_evap | diffusion_wave_implicit_adi | 12.5 | 5 | 24 | 120 |  | 1.5e-02 | 5.0e-16 | 0.0e+00 | no | no |  |
| moat_evap | diffusion_wave_implicit_cg | 577 | 5 | 24 | 120 |  | 1.5e-02 | 5.0e-16 | 0.0e+00 | no | no |  |
| moat_evap | kinematic_wave | 3.43 | 0.0467 | 824 | 120 |  | 1.5e-02 | 4.1e-16 | 0.0e+00 | no | no |  |
| moat_evap | local_inertial | 6.84 | 0.181 | 514 | 120 |  | 1.5e-02 | 9.9e-16 | 0.0e+00 | no | no |  |
| moat_evap | pipes | 6 | 0.152 | 791 | 120 |  | 1.5e-02 | 2.4e-14 | 0.0e+00 | no | no |  |
| moat_evap | swe_hll | 9.76 | 0.135 | 685 | 120 |  | 1.5e-02 | 8.3e-16 | 0.0e+00 | no | no |  |
| moat_evap | swe_lbm | 30.9 | 0.119 | 985 | 120 |  | 1.5e-02 | 1.1e-14 | 0.0e+00 | no | no |  |
| moat_evap | wca2d | 8.23 | 1 | 120 | 120 |  | 1.5e-02 | 3.6e-15 | 0.0e+00 | no | no |  |
| ocean_edge | ca_game | 1.84 | 0.5 | 240 | 120 |  | 2.0e-02 | 1.1e-16 | 0.0e+00 | no | no |  |
| ocean_edge | diffusion_wave_explicit | 3.32 | 5.49e-06 | 10142 | 0.0557 |  | 1.8e-01 | 2.8e-17 | 0.0e+00 | no | yes |  |
| ocean_edge | diffusion_wave_implicit_adi | 440 | 5 | 24 | 120 |  | 1.8e-01 | 2.7e-17 | 0.0e+00 | no | no |  |
| ocean_edge | diffusion_wave_implicit_cg | 801 | 5 | 24 | 120 |  | 2.3e-01 | 4.1e-17 | 0.0e+00 | no | no |  |
| ocean_edge | kinematic_wave | 3.01 | 0.0131 | 160 | 1.83 |  | 5.4e+03 | 1.6e-17 | 0.0e+00 | yes (1.83 s) | no |  |
| ocean_edge | local_inertial | 4.6 | 0.104 | 1158 | 120 |  | 6.5e-02 | 2.8e-17 | 0.0e+00 | no | no |  |
| ocean_edge | pipes | 4.93 | 0.152 | 791 | 120 |  | 5.5e-03 | 3.4e-18 | 0.0e+00 | no | no |  |
| ocean_edge | swe_hll | 10.3 | 0.0777 | 1544 | 120 |  | 6.5e-03 | 3.0e-17 | 0.0e+00 | no | no |  |
| ocean_edge | swe_lbm | 36.6 | 0.0388 | 1360 | 16.9 |  | 3.3e-01 | 4.8e-17 | 0.0e+00 | no | yes |  |
| ocean_edge | wca2d | 6.58 | 0.251 | 823 | 120 |  | 3.0e+00 | 2.7e-17 | 0.0e+00 | no | no |  |
| pit_dig | ca_game | 1.94 | 0.5 | 300 | 150 |  | 1.7e-02 | 6.6e-16 | 0.0e+00 | no | no |  |
| pit_dig | diffusion_wave_explicit | 3.07 | 0.000223 | 17907 | 2.9 |  | 3.0e-02 | 4.2e-16 | 0.0e+00 | no | yes |  |
| pit_dig | diffusion_wave_implicit_adi | 453 | 5 | 30 | 150 |  | 1.5e-01 | 9.2e-16 | 0.0e+00 | no | no |  |
| pit_dig | diffusion_wave_implicit_cg | 782 | 5 | 30 | 150 |  | 2.2e-01 | 4.4e-16 | 0.0e+00 | no | no |  |
| pit_dig | kinematic_wave | 2.91 | 0.115 | 1152 | 150 |  | 1.1e-03 | 4.4e-15 | 0.0e+00 | no | no |  |
| pit_dig | local_inertial | 4.57 | 0.25 | 570 | 150 |  | 5.7e-02 | 1.8e-16 | 0.0e+00 | no | no |  |
| pit_dig | pipes | 4.57 | 0.152 | 988 | 150 |  | 2.5e-02 | 6.8e-15 | 0.0e+00 | no | no |  |
| pit_dig | swe_hll | 12.1 | 0.162 | 1016 | 150 |  | 5.8e-03 | 1.1e-16 | 0.0e+00 | no | no |  |
| pit_dig | swe_lbm | 33.1 | 0.151 | 2815 | 42.3 |  | 7.7e-02 | 1.3e-16 | 0.0e+00 | no | yes |  |
| pit_dig | wca2d | 7.82 | 0.184 | 1222 | 150 |  | 8.0e+00 | 2.9e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | ca_game | 2.38 | 0.5 | 240 | 120 |  | 1.3e-02 | 7.6e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | diffusion_wave_explicit | 4.54 | 0.00217 | 9865 | 9.89 |  | 6.7e-02 | 3.2e-15 | 0.0e+00 | no | yes |  |
| rain_groundwater | diffusion_wave_implicit_adi | 460 | 5 | 24 | 120 |  | 5.2e-02 | 5.9e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | diffusion_wave_implicit_cg | 1.18e+03 | 5 | 24 | 120 |  | 3.8e-02 | 4.5e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | kinematic_wave | 3.74 | 0.165 | 2246 | 120 |  | 7.2e+01 | 5.9e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | local_inertial | 5.82 | 0.36 | 510 | 120 |  | 2.2e-02 | 2.8e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | pipes | 6.79 | 0.152 | 791 | 120 |  | 9.5e-03 | 1.2e-15 | 0.0e+00 | no | no |  |
| rain_groundwater | swe_hll | 9.07 | 0.247 | 710 | 120 |  | 7.9e-03 | 5.6e-16 | 0.0e+00 | no | no |  |
| rain_groundwater | swe_lbm | 31.5 | 0.191 | 1868 | 36 |  | 1.4e-02 | 3.1e-16 | 0.0e+00 | no | yes |  |
| rain_groundwater | wca2d | 8.17 | 0.415 | 508 | 120 |  | 2.6e-02 | 6.4e-16 | 0.0e+00 | no | no |  |
| ripples | linear_wave | 0.762 | 0.0636 | 315 | 20 |  |  |  |  | no | no |  |
| river | ca_game | 1.83 | 0.5 | 240 | 120 |  | 2.0e-02 | 1.1e-15 | 0.0e+00 | no | no |  |
| river | diffusion_wave_explicit | 3.22 | 0.000199 | 10374 | 2.81 |  | 2.8e-02 | 1.7e-16 | 0.0e+00 | no | yes |  |
| river | diffusion_wave_implicit_adi | 454 | 5 | 24 | 120 |  | 1.1e-02 | 8.4e-16 | 0.0e+00 | no | no |  |
| river | diffusion_wave_implicit_cg | 869 | 5 | 24 | 120 |  | 2.4e-04 | 1.9e-16 | 0.0e+00 | no | no |  |
| river | kinematic_wave | 3.11 | 0.0857 | 992 | 120 |  | 2.7e+00 | 2.3e-16 | 0.0e+00 | no | no |  |
| river | local_inertial | 4.55 | 0.242 | 490 | 120 |  | 2.1e-02 | 7.8e-17 | 0.0e+00 | no | no |  |
| river | pipes | 4.65 | 0.152 | 791 | 120 |  | 2.9e-02 | 3.6e-15 | 0.0e+00 | no | no |  |
| river | swe_hll | 9.23 | 0.162 | 780 | 120 |  | 2.0e-03 | 1.9e-16 | 0.0e+00 | no | no |  |
| river | swe_lbm | 29 | 0.145 | 1010 | 120 |  | 1.5e-02 | 1.6e-16 | 0.0e+00 | no | no |  |
| river | wca2d | 6.08 | 0.183 | 973 | 120 |  | 1.3e-01 | 4.7e-16 | 0.0e+00 | no | no |  |
| smoke_plume | lbm_d2q9 | 89.4 | 1 | 141 | 141 |  |  |  |  | yes (141 s) | no |  |
| smoke_plume | stable_fluids | 111 | 1 | 520 | 520 |  |  |  |  | no | yes |  |
| swell | ca_game | 1.75 | 0.5 | 120 | 60 |  | 2.3e-01 | 1.0e-17 | 0.0e+00 | no | no |  |
| swell | diffusion_wave_explicit | 3.54 | 3.82e-06 | 9487 | 0.0362 |  | 2.4e-01 | 3.4e-15 | 0.0e+00 | no | yes |  |
| swell | diffusion_wave_implicit_adi | 13.6 | 5 | 12 | 60 |  | 4.2e-02 | 1.4e-16 | 0.0e+00 | no | no |  |
| swell | diffusion_wave_implicit_cg | 1.09e+03 | 5 | 12 | 60 |  | 4.2e-02 | 1.1e-16 | 0.0e+00 | no | no |  |
| swell | kinematic_wave | 3.01 | 0.0115 | 146 | 1.49 |  | 7.7e+03 | 3.7e-17 | 0.0e+00 | yes (1.49 s) | no |  |
| swell | local_inertial | 4.19 | 0.0922 | 649 | 60 |  | 4.3e-01 | 1.5e-16 | 0.0e+00 | no | no |  |
| swell | pipes | 4.84 | 0.152 | 396 | 60 |  | 4.2e-01 | 2.0e-17 | 0.0e+00 | no | no |  |
| swell | swe_hll | 10.5 | 0.0683 | 883 | 60 |  | 2.4e-01 | 4.0e-17 | 0.0e+00 | no | no |  |
| swell | swe_lbm | 33.5 | 0.0611 | 1327 | 60 |  | 2.4e-01 | 2.5e-17 | 0.0e+00 | no | no |  |
| swell | wca2d | 6.76 | 0.0768 | 765 | 60 |  | 4.6e-01 | 5.3e-17 | 0.0e+00 | no | no |  |
