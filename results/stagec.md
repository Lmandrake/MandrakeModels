# Stage C: local-inertial in C# on integer decilitres (2026-11-13)

Code: `csharp/` (see `csharp/README.md`). Run: `bash csharp/run_stagec.sh`; raw output
`D:\ClaudeWork\mandrakemodels\stagec\out\` (timing_final.csv, fidelity_*.jsonl, stability_*.jsonl, fuzz_*.jsonl), drift probes
`D:\ClaudeWork\mandrakemodels\stagec\drift\`. All figures below are [measured] unless tagged.

## Runtimes

- **Unity Mono 6.13** ("Visual Studio built mono"): RimWorld's own `mono-2.0-bdwgc.dll` and `RimWorldWin64_Data\Managed` class
  library, embedded from Windows Python by ctypes (`csharp/MonoHost/monohost.py`), net472 build. Same JIT and class library the
  game runs; not inside Unity, so the game's own main-thread work is not competing for cache [guess ~70% that in-game cost is
  within +30% of these figures].
- .NET Framework 4.8 (Windows CLR, same net472 exe) and .NET 8.0.29 (net8.0 build), for comparison.
- Machine: Archmagi, 14 cores, shared (load 3 to 11; RimWorldWin64.exe was running on the desktop throughout). Single thread,
  process priority High. Unity Mono matrix run twice; each row reports the worse repeat, repeats agreed within 5%.

## Port and tolerances (fixed before the first run)

Water is an int per square (dL); each face keeps a double discharge and a sub-dL rounding residue carried to the next step;
integer donor limiter; each face moves whole dL once with opposite signs. Wet-only stepping: 16x16 blocks holding water, dilated
by one block.

A scene passes when, at all 4 checkpoints: relative L1 depth difference against Python <= 2%; wet-area IoU (h > 1 mm) >= 0.98;
books exactly 0 dL; wet-only bitwise equal to full grid; lake at rest moves < 1 mm; pit_dig pit mean depth within 2%; dam_break
volume east of the dam within 2%. Result: **40 of 40 rows pass on Unity Mono and on .NET 8** (river, lake_at_rest, pit_dig,
dam_break, swell at 128 and 250); worst relative L1 9.6e-5, min IoU 0.9999, books 0 dL, wet-only identical. Checker proven:
Manning x1.3 turns 16 rows red, theta 1.0 turns 32 red.

## Findings

1. **A recompute must split into equal substeps.** Splitting each 0.2 s recompute as full CFL steps plus a short remainder
   (0.17 s, 0.03 s, ...) made the lake at rest slosh 2.1 m within 20 water-s, in Python float64 too (drift probe: ragged 2.1 m,
   equal 4e-16 m), so it is the scheme, not the integers. With k = ceil(0.2 / CFL dt) equal substeps the C# lake holds 0.2 mm
   over 10,000 s. Equal substeps cost more steps than CFL needs (2 per recompute at 1 m depth where 1.2 would do).
2. **Unity Mono is 2.1 to 4.9 times slower than .NET 8** on this kernel (typically 3); .NET Framework 4.8 is within 25% of .NET 8.
   Profiling: the face pass costs ~55 ns per wet square-step under Mono against ~18 under .NET 8; removing bounds checks
   (pointers) and edge tests, and a cheaper cube root, gained only ~5%.
3. **The recompute is a spike.** Its whole cost lands in one tick unless spread over the 12: on the realistic 325 river map
   1.9 ms in that tick (fine), on the 50%-wet stress 17 to 20 ms (a whole frame at speed 1).
4. **Fuzz found one real bug:** forcing kept writing into a square walled over a source or level-setter, and demolishing it
   then deleted that water. Fixed (forcing skips walled squares); a planted fault (a wall deletes its water) turns 145 of 200
   sequences red.

## Bars (plan section 5.4; Unity Mono, wet-only, mean ms per tick)

| bar | 250 | 325 (binding) | 400 (measured) |
|---|---|---|---|
| realistic river map (census_river, 13,178 / 17,121 / 21,094 wet) under 0.5 ms | 0.123 PASS | **0.163 PASS** (p95 0.186) | 0.206 PASS |
| 50%-wet stress (up to 10 m deep, swelling edge) under 2 ms | 0.855 PASS | **1.452 PASS** | 2.154 FAIL |
| settled under 0.05 ms | no sleeping is built: a settled lake costs what a moving one does (lake_at_rest 0.573 / 0.859 / 1.376) | **FAIL** | FAIL |
| books 0 dL over 100,000 recomputes and 1,000 fuzz sequences; no square below empty | PASS (river 250: 100,000 recomputes = 20,000 water-s; fuzz 1,000 x 150 recomputes with 50,037 digs, earth dams, walls and demolitions, on river 128 under Mono and pit_dig 128 under .NET 8) | | |
| lake at rest under 1 mm over 10,000 river-s, 0 dL | PASS (250: 0.20 mm, 0 dL, Unity Mono and .NET 8) | | |
| long runs stay finite (swell, dam_break, pit_dig, stress_half 2,000 water-s at 250) | PASS, books 0 dL | | |

### Unity Mono 6.13 (RimWorld's own runtime), wet-only stepping

ms per tick = the 12-tick recompute's cost spread over its 12 ticks; the recompute itself lands in one tick unless split.

| scene | n | wet squares | substeps per recompute | ms/tick mean | median | p95 | ms per recompute median | p95 | max | repeat spread | books dL |
|---|---|---|---|---|---|---|---|---|---|---|---|
| census_river | 250 | 13178 | 1.0 | 0.123 | 0.122 | 0.136 | 1.46 | 1.64 | 4.39 | 0% | 0 |
| census_river | 325 | 17121 | 1.0 | 0.163 | 0.161 | 0.186 | 1.93 | 2.24 | 2.59 | 2% | 0 |
| census_river | 400 | 21094 | 1.0 | 0.206 | 0.202 | 0.233 | 2.43 | 2.80 | 3.32 | 3% | 0 |
| river | 250 | 3224 | 1.0 | 0.044 | 0.043 | 0.045 | 0.52 | 0.54 | 0.73 | 1% | 0 |
| river | 325 | 5860 | 1.0 | 0.071 | 0.070 | 0.076 | 0.84 | 0.91 | 1.24 | 1% | 0 |
| river | 400 | 8758 | 1.0 | 0.100 | 0.099 | 0.103 | 1.19 | 1.24 | 1.57 | 1% | 0 |
| pit_dig | 250 | 4781 | 1.0 | 0.051 | 0.054 | 0.059 | 0.64 | 0.71 | 0.77 | 0% | 0 |
| pit_dig | 325 | 7537 | 1.0 | 0.080 | 0.083 | 0.096 | 0.99 | 1.15 | 1.20 | 0% | 0 |
| pit_dig | 400 | 11079 | 1.0 | 0.110 | 0.111 | 0.130 | 1.33 | 1.55 | 4.46 | 1% | 0 |
| lake_at_rest | 250 | 24873 | 2.0 | 0.573 | 0.565 | 0.630 | 6.78 | 7.56 | 9.66 | 4% | 0 |
| lake_at_rest | 325 | 38553 | 2.0 | 0.859 | 0.857 | 0.892 | 10.28 | 10.70 | 12.75 | 0% | 0 |
| lake_at_rest | 400 | 63314 | 2.0 | 1.376 | 1.371 | 1.439 | 16.45 | 17.27 | 17.65 | 0% | 0 |
| dam_break | 250 | 62124 | 2.0 | 1.079 | 1.077 | 1.165 | 12.93 | 13.98 | 15.76 | 0% | 0 |
| dam_break | 325 | 105135 | 2.0 | 1.866 | 1.847 | 2.014 | 22.17 | 24.16 | 25.94 | 1% | 0 |
| dam_break | 400 | 159400 | 2.0 | 2.845 | 2.820 | 2.991 | 33.84 | 35.90 | 37.72 | 2% | 0 |
| swell | 250 | 48564 | 3.0 | 1.313 | 1.293 | 1.446 | 15.51 | 17.35 | 23.24 | 0% | 0 |
| swell | 325 | 83686 | 3.0 | 2.316 | 2.301 | 2.546 | 27.61 | 30.55 | 35.70 | 0% | 0 |
| swell | 400 | 126113 | 3.0 | 3.564 | 3.565 | 3.857 | 42.78 | 46.28 | 49.30 | 1% | 0 |
| stress_half | 250 | 32462 | 3.0 | 0.855 | 0.832 | 0.958 | 9.98 | 11.49 | 12.55 | 0% | 0 |
| stress_half | 325 | 53853 | 3.0 | 1.452 | 1.400 | 1.662 | 16.80 | 19.95 | 25.96 | 1% | 0 |
| stress_half | 400 | 79656 | 3.0 | 2.154 | 2.105 | 2.415 | 25.26 | 28.98 | 32.18 | 1% | 0 |

### Full grid against wet-only, and the three runtimes (mean ms per tick)

| scene | n | Mono full | Mono wet | .NET Fx 4.8 wet | .NET 8 wet | Mono / .NET 8 |
|---|---|---|---|---|---|---|
| census_river | 250 | 0.166 | 0.123 | 0.043 | 0.035 | 3.5 |
| census_river | 325 | 0.244 | 0.163 | 0.058 | 0.045 | 3.6 |
| census_river | 400 | 0.345 | 0.206 | 0.072 | 0.057 | 3.6 |
| river | 250 | 0.099 | 0.044 | 0.022 | 0.016 | 2.8 |
| river | 325 | 0.164 | 0.071 | 0.032 | 0.024 | 2.9 |
| river | 400 | 0.252 | 0.100 | 0.044 | 0.034 | 2.9 |
| pit_dig | 250 | 0.104 | 0.051 | 0.022 | 0.017 | 3.0 |
| pit_dig | 325 | 0.175 | 0.080 | 0.034 | 0.026 | 3.0 |
| pit_dig | 400 | 0.264 | 0.110 | 0.043 | 0.035 | 3.2 |
| lake_at_rest | 250 | 0.609 | 0.573 | 0.283 | 0.260 | 2.2 |
| lake_at_rest | 325 | 0.993 | 0.859 | 0.434 | 0.404 | 2.1 |
| lake_at_rest | 400 | 1.517 | 1.376 | 0.700 | 0.653 | 2.1 |
| dam_break | 250 | 1.081 | 1.079 | 0.253 | 0.214 | 5.0 |
| dam_break | 325 | 1.846 | 1.866 | 0.441 | 0.381 | 4.9 |
| dam_break | 400 | 2.828 | 2.845 | 0.675 | 0.569 | 5.0 |
| swell | 250 | 1.337 | 1.313 | 0.503 | 0.453 | 2.9 |
| swell | 325 | 2.345 | 2.316 | 0.920 | 0.899 | 2.6 |
| swell | 400 | 3.619 | 3.564 | 1.461 | 1.342 | 2.7 |
| stress_half | 250 | 0.941 | 0.855 | 0.297 | 0.272 | 3.1 |
| stress_half | 325 | 1.615 | 1.452 | 0.519 | 0.464 | 3.1 |
| stress_half | 400 | 2.405 | 2.154 | 0.786 | 0.712 | 3.0 |

fidelity_mono.jsonl: 40/40 rows pass; worst relative L1 9.63e-05; min IoU 0.9999; books [0]; wet-only identical True

fidelity_net8.jsonl: 40/40 rows pass; worst relative L1 9.63e-05; min IoU 0.9999; books [0]; wet-only identical True

stability_mono.jsonl:
- lake_at_rest_250.lis: 10000 water-s, 50000 recomputes, 100000 steps, worst books 0 dL, min square 0 dL, max depth 2.82 m, lake surface moved 0.20 mm: PASS
- river_250.lis: 20000 water-s, 100000 recomputes, 100000 steps, worst books 0 dL, min square 0 dL, max depth 1.26 m: PASS

stability_net8.jsonl:
- lake_at_rest_250.lis: 10000 water-s, 50000 recomputes, 100000 steps, worst books 0 dL, min square 0 dL, max depth 2.82 m, lake surface moved 0.20 mm: PASS
- river_250.lis: 20000 water-s, 100000 recomputes, 100000 steps, worst books 0 dL, min square 0 dL, max depth 1.26 m: PASS
- swell_250.lis: 2000 water-s, 10000 recomputes, 30000 steps, worst books 0 dL, min square 0 dL, max depth 8.74 m: PASS
- dam_break_250.lis: 2000 water-s, 10000 recomputes, 10207 steps, worst books 0 dL, min square 0 dL, max depth 2.00 m: PASS
- pit_dig_250.lis: 2000 water-s, 10000 recomputes, 10000 steps, worst books 0 dL, min square 0 dL, max depth 1.50 m: PASS
- stress_half_250.lis: 2000 water-s, 10000 recomputes, 30000 steps, worst books 0 dL, min square 0 dL, max depth 10.66 m: PASS

fuzz_mono.jsonl: 1000 sequences, 50037 edits, 150000 recompute checks, 0 failures (Mono 6.13.0 (Visual Studio built mono))

fuzz_net8_pit.jsonl: 1000 sequences, 50037 edits, 150000 recompute checks, 0 failures (.NET 8.0.29)
