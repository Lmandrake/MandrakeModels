# Stage 2: review fixes, new competitors, rerun (2026-11-13)

Input: a skeptical external review of stage 1 (a ranked fix list, each claim checked before acting on it). Every figure is
[measured] unless tagged. The probe scripts and per-scene viewing notes are working files kept outside this repository.

## New competitors (originals unchanged beside them)

- `local_inertial_frcap`: every face discharge capped at the broad-crested-weir value, i.e. face Froude <= (2/3)^1.5 = 0.544.
- `local_inertial_frcap07`: the same cap at Froude 0.7.
- `pipes_theta`: pipes with theta 0.8 smoothing of each one-way memory and a dry-face gate at any manning (visual only:
  hydraulics as pipes).

| solver | breach vs swe_hll (dam_break 128, 20 s) | pit depth 150 s, 64 / 128 / 250 (swe_hll 1.12 / 0.42 / 0.26) | river to outfall, 128 (swe_hll 774 m^3) | ms/tick numpy 250 / 325 |
|---|---|---|---|---|
| swe_hll (referee) | 1.00 | 1.12 / 0.42 / 0.26 | 774 | 1.19 / 2.53 |
| local_inertial | 1.48 | 1.55 / 0.72 / 0.49 | 779 | 0.40 / 0.74 |
| **local_inertial_frcap** | **1.09** | **1.30 / 0.44 / 0.26** | 683 (-12%) | 0.37 / 0.75 |
| local_inertial_frcap07 | 1.28 | not run / 0.50 / 0.33 | 718 (-7%) | 0.37 / 0.75 |
| pipes | 1.59 | 1.42 / 0.77 / 0.58 | 841 at 7.8 m/s max (REVIEW) | 0.49 / 1.04 |
| pipes_theta | 1.56 | not measured; early like pipes on screen | 1266 held vs 841 (REVIEW) | 0.57 / 1.14 |

Cap sweep at 128 (`probes\out_frsweep_stage2.json`): breach 1.09 / 1.17 / 1.23 / 1.28 and river 683 / 698 / 709 / 718 m^3 at
Froude 0.544 / 0.6 / 0.65 / 0.7. **No cap holds both the breach within 15% and the river within 5%**: the review's hope of
one value near 0.6 doing both [guess ~55% then] is refuted. The river deficit closes only with no cap.

## On screen (LOOKED2.md)

- **pipes_theta removes the pipes ringing in every scene that had it** (dam_break, river, lake_levelset, ocean_edge, pit_dig,
  rain_groundwater, swell). Floods stay early (it is pipes underneath).
- **The Froude cap fixes local_inertial's timing on screen** (pit_dig front and stage match swe_hll) but **not the stripes
  in the dam-break reservoir**: those are wall reflections that local_inertial never damps (upstream surface roughness at
  32.7 s: swe_hll 1.8e-4 m, local_inertial 4.1e-4, frcap 3.1e-4). Mixing each face with its cross-axis neighbours as well
  cut it to 2.3e-4 but cost the river 43% and made the 250 lake at rest move 2.9 cm: rejected (`probes2\stripes.py`,
  `out_cross_safety.json`).

## Harness fixes (rows re-run)

- e157c92: stable_dt is timed with the step (wca2d's ms/tick at 325 rose 0.86 to 1.62); the max-dt search stops at its cap
  (wca2d's 81.7 s was the overshoot; now >= 64 cap) and simulates exactly T. Every water timing row and every max-dt row re-run.
- 24e4485: no C# ms/tick is projected from the numba probe any more; its kernels are simplified solvers (pipes without
  friction, local_inertial without theta), so the 39x pipes ratio priced different arithmetic. Measured C#: `stagec.md`.
- b774658: ledger throughput sums every booking's size (a level-setter's in and out no longer cancel); rain is sampled at
  the step midpoint. rain_groundwater rows re-run at all sizes; level-scene rows at 250 and up keep the older, smaller
  throughput (their vol err reads larger than it is; every one is below 1e-14).
- ec367a1: ca_game moves the fraction dt/pass of a pass per step. compare.py clips dt to frame times, so the compare pages
  ran the CA about 1.5x fast; bakeoff rows (dt always 0.5 s) are unchanged.
- 2d0ba1f: kinematic_wave out of the water headline (bed-slope only); its rows stay in the csv.

## Not done from the review, and why

- Rank 7, second half (timing windows after the last event): not done; timing still samples steps 4 to 33 from t = 0.
- Rank 8 (WCA consistent weights): not entered, as the review advised (late and noisy even repaired).
- Rank 9 (river inlet as a boundary-face flux): not done; the CA inlet "spray" is still on screen in river and pit_dig.
- Section 6 boundary proposals (free outfall, Flather swell edge): not done.
