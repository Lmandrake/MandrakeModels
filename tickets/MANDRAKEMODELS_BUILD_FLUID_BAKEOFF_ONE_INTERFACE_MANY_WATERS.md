---
id: MANDRAKEMODELS_BUILD_FLUID_BAKEOFF_ONE_INTERFACE_MANY_WATERS
title: Python bakeoff of gridded fluid solvers (water, waves, diffusion, Navier-Stokes, LBM) at RimWorld map sizes
project: mandrakemodels
type: build
status: open
owner: lukas
created: 2026-11-12T14:14:01-08:00
updated: 2026-11-12T17:20:00-08:00
effort: 6h
source: manual
check: "{\"reason\":\"bakeoff results table and gallery not yet produced\",\"review_by\":\"2026-11-26\",\"type\":\"gap\"}"
---

## Log
- 2026-11-12T14:14 archmagi/ticket: new (manual)
- 2026-11-12T17:20 bench/builder: bakeoff landed (594 runs, 64-500, results/results.md), compare tool and
  cross-check (ANUGA, Landlab, PhiFlow; results/crosscheck.md). 76 tests pass.
  OWED: PhiFlow smoke adapter disagrees at 64^2 (spurious speeds; likely adapter obstacles); pylbm adapter;
  re-time on an idle machine (timing taken at load 7-34); a C# port to measure the guessed column;
  compare at 128+ for the externals (slow: 40-900 ms per simulated second at 64).
NEXT: owner looks at D:\ClaudeWork\mandrakemodels\compare\compare.html and picks the water family.
