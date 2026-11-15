#!/usr/bin/env bash
# Stage C full run on Archmagi (WSL): build, fidelity, timing on three runtimes, stability, fuzz.
#   bash csharp/run_stagec.sh            (from the MandrakeModels root)
#   REF=local_inertial bash csharp/run_stagec.sh                  the uncapped reference (stage C as first measured)
#   STEPS=fidelity,timing SIZES=250,325 bash csharp/run_stagec.sh  a subset; timing on the sizes exported
#   steps: evensteps (every driver takes equal steps, one short interval per edit at most: LocalInertialDL.EvenDt), fidelity (40 rows vs Python; a sleeping twin within the same tolerances), timing (modes full, wet, sleep),
#   spread (each recompute over its 12 ticks, bitwise against whole), settled (sleeping lake bar), wake (a dig into
#   sleeping water flows as in water that never slept), stability, fuzz
#   LITIMING_MUTATE=frcap1 STEPS=fidelity OUT=... bash csharp/run_stagec.sh   a planted fault (passed to Windows too)
# REF names the Python reference solver (default local_inertial_frcap, the chosen engine). Its face Froude cap travels in
# the scene files, so the C# runs the same scheme. Scenes come from `python -m mandrakemodels.csharp_export --solver $REF`
# into $STAGEC/scenes_$REF (with checkpoints) and $STAGEC/load_$REF (timing loads), exported here when absent.
# Output: $STAGEC/out_$REF (or $OUT).
set -uo pipefail
STAGEC=${STAGEC:-/mnt/d/ClaudeWork/mandrakemodels/stagec}
REF=${REF:-local_inertial_frcap}
STEPS=${STEPS:-evensteps,fidelity,timing,spread,settled,wake,stability,fuzz}
SIZES=${SIZES:-250,325,400}
W=$(wslpath -w "$STAGEC")
DOTNET=${DOTNET:-/mnt/c/Users/Mandrake/.dotnet/dotnet.exe}
export WSLENV=${WSLENV:+$WSLENV:}LITIMING_MUTATE
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SC="$STAGEC/scenes_$REF"; LD="$STAGEC/load_$REF"
WSC="$W\\scenes_$REF"; WLD="$W\\load_$REF"
has() { [[ ",$STEPS," == *",$1,"* ]]; }
if [ ! -d "$SC" ]; then
  ( cd "$ROOT" && uv run python -m mandrakemodels.csharp_export --out "$SC" --sizes 128,250 --solver "$REF" ) || exit 1
fi
if [ ! -d "$LD" ]; then
  ( cd "$ROOT" && uv run python -m mandrakemodels.csharp_export --out "$LD" --sizes "$SIZES" --no-ref --solver "$REF" \
      --scenes river,lake_at_rest,pit_dig,dam_break,swell,census_river,stress_half ) || exit 1
fi
rsync -a --delete "$(dirname "$0")/" "$STAGEC/src/"
( cd "$STAGEC/src/LiTiming" && "$DOTNET" build -c Release --nologo -v quiet ) || exit 1
EXE8="$W\\src\\LiTiming\\bin\\Release\\net8.0\\LiTiming.dll"
EXE4="$W\\src\\LiTiming\\bin\\Release\\net472\\LiTiming.exe"
HOST="$W\\src\\MonoHost\\monohost.py"
mono() { python.exe "$HOST" "$EXE4" "$@"; }
net8() { "$DOTNET" "$EXE8" "$@"; }
netfx() { "$STAGEC/src/LiTiming/bin/Release/net472/LiTiming.exe" "$@"; }
O=${OUT:-$STAGEC/out_$REF}
mkdir -p "$O"; WO=$(wslpath -w "$O")
step() { echo "$(date '+%H:%M') $*" >> "$O/run_stagec.milestones"; }
rc=0
step "start ref=$REF steps=$STEPS mutate=${LITIMING_MUTATE:-none}"
if has evensteps; then
  : > "$O/evensteps_net8.log"
  for f in dam_break_128 pit_dig_128 swell_128; do
    net8 evensteps "$WSC\\$f.lis" >> "$O/evensteps_net8.log" 2>&1 || rc=1
  done; step "evensteps net8: $(grep -c PASS "$O/evensteps_net8.log") of 3 PASS"
fi
if has fidelity; then
  mono fidelity "$WSC" "$WO\\fidelity_mono.jsonl" > "$O/fidelity_mono.log" 2>&1 || rc=1; step "fidelity mono: $(tail -1 "$O/fidelity_mono.log")"
  net8 fidelity "$WSC" "$WO\\fidelity_net8.jsonl" > "$O/fidelity_net8.log" 2>&1 || rc=1; step "fidelity net8: $(tail -1 "$O/fidelity_net8.log")"
fi
if has timing; then
  rm -f "$O/timing_final.csv"
  for rep in 1 2; do
    mono timing "$WLD" "$WO\\timing_final.csv" --modes full,wet,sleep --label "unity_mono_6.13 r$rep" > "$O/timing_mono_r$rep.log" 2>&1; step "timing mono r$rep"
  done
  net8 timing "$WLD" "$WO\\timing_final.csv" --label "net8 r1" > "$O/timing_net8_r1.log" 2>&1; step "timing net8"
  netfx timing "$WLD" "$WO\\timing_final.csv" --label "netfx48 r1" > "$O/timing_netfx_r1.log" 2>&1; step "timing netfx"
fi
if has spread; then
  rm -f "$O/spread_mono.csv"
  mono spread "$WLD" "$WO\\spread_mono.csv" --modes wet,sleep > "$O/spread_mono.log" 2>&1 || rc=1; step "spread mono: $(tail -1 "$O/spread_mono.log")"
fi
if has settled; then
  rm -f "$O/settled_mono.csv"
  mono settled "$WLD" "$WO\\settled_mono.csv" > "$O/settled_mono.log" 2>&1 || rc=1; step "settled mono: $(tail -1 "$O/settled_mono.log")"
fi
if has wake; then
  rm -f "$O/wake_mono.jsonl"
  for f in lake_at_rest_250 lake_at_rest_325; do
    mono wake "$WLD\\$f.lis" "$WO\\wake_mono.jsonl" > "$O/wake_mono_$f.log" 2>&1 || rc=1
  done; step "wake mono: $(tail -1 "$O/wake_mono_lake_at_rest_325.log")"
fi
if has stability; then
  rm -f "$O/stability_mono.jsonl" "$O/stability_net8.jsonl"
  mono stability "$WLD" "$WO\\stability_mono.jsonl" --only lake_at_rest > "$O/stability_mono.log" 2>&1 || rc=1
  mono stability "$WLD" "$WO\\stability_mono.jsonl" --only river > "$O/stability_mono_river.log" 2>&1 || rc=1; step "stability mono (lake, river)"
  net8 stability "$WLD" "$WO\\stability_net8.jsonl" > "$O/stability_net8.log" 2>&1 || rc=1; step "stability net8: $(tail -1 "$O/stability_net8.log")"
  net8 stability "$WLD" "$WO\\stability_net8.jsonl" --modes sleep > "$O/stability_net8_sleep.log" 2>&1 || rc=1; step "stability net8 sleep: $(tail -1 "$O/stability_net8_sleep.log")"
fi
if has fuzz; then
  mono fuzz "$WSC\\river_128.lis" "$WO\\fuzz_mono.jsonl" --seqs 1000 > "$O/fuzz_mono.log" 2>&1 || rc=1; step "fuzz mono: $(tail -1 "$O/fuzz_mono.log")"
  net8 fuzz "$WSC\\pit_dig_128.lis" "$WO\\fuzz_net8_pit.jsonl" --seqs 1000 > "$O/fuzz_net8_pit.log" 2>&1 || rc=1; step "fuzz net8 pit: $(tail -1 "$O/fuzz_net8_pit.log")"
fi
step "done rc=$rc"
exit $rc
