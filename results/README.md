# Curated results

Small outputs copied into the repository from full runs. The full runs, with every GIF, frame and
`gallery.html`, stay under the output root and are not committed.

All of these were measured on the author's machines; `crosscheck.md` and `stagec.md` name a shared 14-core workstation.
Absolute timings depend on the machine and its load. The ratios and pass/fail findings are the useful part.

## Bakeoff (written by `bakeoff report --curate results`)

| file | what it is |
|---|---|
| `results.md` | the generated report, with these sections: the numba probe table; water headlines at 325 and 500 (median over water scenes: numpy ms/step, dt, ms per game tick, worst volume error, lake drift, settled/blow-up/error/capped counts); numpy ms/step by size; largest surviving fixed dt; flood timing through the dam breach against `swe_hll`; and every run at 250 |
| `results.csv` | one row per (solver, scene, n): 756 rows, 18 solvers, sizes 64 to 500. Columns are documented in `docs/bakeoff.md` (quality and timing fields) |
| `maxdt.jsonl` | the `maxdt` phase: largest fixed dt each solver survives on `dam_break` 128 / `gaussian` 64 / `ripples` 64, and its own `stable_dt` |
| `probe.json` | the numba probe at 400^2. The kernels are simplified solvers, so the ratio is not a C# cost estimate |
| `sheets/<scene>_250.png` | contact sheets: every solver's mid and final frame at 250, one sheet per scene |

These predate the ripple-damper variants (`local_inertial_frcap_*` beyond `frcap` and `frcap07`), which were
added on 2026-11-14 and have not been through the bakeoff. Their figures are in their docstrings and tests.

## Comparison outputs (`compare/`)

Five scenes at 64^2 from `compare --scene all --n 64 --external` (2026-11-12): `dam_break`, `river`,
`ocean_edge`, `swell` and `smoke_plume`. Each folder has:

- `state.mp4`: every solver tiled
- `final_state.webp`, `final_diff.webp`: the last frame, state and difference from the reference
- `meta.json`: solvers, sample times, per-solver final status, volume error and RMS against the reference

These are a curated subset (one video and two stills per scene) kept small for the repository, not the full `compare` output set described in
`docs/compare.md`. `final_state.webp` and `final_diff.webp` are not written by the current `compare.py`.

## Write-ups

| file | what it is |
|---|---|
| `crosscheck.md` | our solvers against ANUGA, Landlab and PhiFlow from the same grid at 64^2. ANUGA agrees with `swe_hll` (depth RMS at most 0.04 m on 8 of 9 water scenes) |
| `stage2.md` | stage 2 (2026-11-13): fixes from an external review, the new competitors `local_inertial_frcap`, `local_inertial_frcap07` and `pipes_theta`, a Froude-cap sweep, and notes on what each scene looks like on screen |
| `stagec.md` | stage C (2026-11-13 and 14): the C# integer port, its fidelity against Python, cost per game tick on Unity Mono, .NET Framework 4.8 and .NET 8, the bars it was held to, and the equal-steps finding |

`stage2.md` and `stagec.md` are working records, left as written. They cite paths on the author's machines
(`D:\ClaudeWork\...`) and work items in his private tracker (names such as `RIMMANDRAKE_...`). Those files and
items are not part of this repository; the write-ups carry the figures that matter.
