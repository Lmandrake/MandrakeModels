---
id: MANDRAKEMODELS_BUILD_FLUID_BAKEOFF_ONE_INTERFACE_MANY_WATERS
title: Python bakeoff of gridded fluid solvers (water, waves, diffusion, Navier-Stokes, LBM) at RimWorld map sizes
project: mandrakemodels
type: build
status: open
owner: lukas
created: 2026-11-12T14:14:01-08:00
updated: 2026-11-16T07:16:44-08:00
effort: 6h
source: manual
check: "{\"reason\":\"bakeoff results table and gallery not yet produced\",\"review_by\":\"2026-11-26\",\"type\":\"gap\"}"
---

2026-11-16: repo made PUBLIC (MIT) with full docs (README, docs/, results/README.md); private paths/tracker ids stripped. Still owed: register ripsense, refresh results/ with ripple dampers.

## Log
- 2026-11-12T14:14 ticket: new (manual)
- 2026-11-12T17:20 bench/builder: bakeoff landed (594 runs, 64-500, results/results.md), compare tool and
  cross-check (ANUGA, Landlab, PhiFlow; results/crosscheck.md). 76 tests pass.
  OWED: PhiFlow smoke adapter disagrees at 64^2 (spurious speeds; likely adapter obstacles); pylbm adapter;
  re-time on an idle machine (timing taken at load 7-34); a C# port to measure the guessed column;
  compare at 128+ for the externals (slow: 40-900 ms per simulated second at 64).
- 2026-11-13 bench/builder: viewer rework after the owner said the compare images "don't look like anything
  we'd see in a water-surface simulation". DONE: waterview.py (earth hillshade, blue depth, surface-slope
  light/dark bands, shoreline, stone walls), compare.py (water scenes first, plain names + code name +
  volume error + ms/step, red/orange banner frozen at last good frame, common clock, side-view strip for
  dam_break/swell, 120 frames, tiles 384 px at n=128), index.html + one page per scene (embedded webm).
  Old compare.html deleted. Re-rendered at 128 with --external: dam_break, river, lake_levelset,
  ocean_edge, pit_dig, moat_evap. Stills: <output root>/compare/stills/ (18).
  Looked at: dam_break, river, lake, ocean_edge, pit_dig (all read as water); moat_evap not yet looked at.
  DEFECT FOUND: pit_dig at n>=128 never floods the pit (the canal is cut a fixed 2 m into a bank that
  climbs away from the river, leaving a sill). Fix drafted, not landed: cut the canal to a level bed
  (z = 0.3) and stop it at the pit; then re-run the bakeoff pit_dig rows (they measured the dry-pit scene).
NEXT, in order: fix pit_dig; render rain_groundwater, swell, lake_at_rest, then other physics, one at a
  time (run1.sh in the output root); look at moat + 3 more; optional dam_break n=250
  showcase; replace results/compare/*_64 with the new meta.json + stills; owner looks at index.html.
- 2026-11-12 bench/detached: pit_dig fixed (canal bed falls from the river bed
  to a smaller, nearer pit; 03ba39e) and moat_evap made to dry visibly (tilted floor, 15 mm/s; fad3ee1), each with a test
  that fails on the old scene (7/7) and passes now; pit_dig and moat_evap bakeoff rows re-run (6840ddf and this commit);
  rendered at 128 and looked at: all nine water scenes (index.html now lists them). Viewing notes kept outside the repo. Still owed: other-physics pages, PhiFlow
  adapter at 64, pylbm adapter, idle re-time, C# measurement.
- 2026-11-16T07:16 mt-400028/ticket: note
