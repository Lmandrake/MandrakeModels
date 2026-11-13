#!/usr/bin/env bash
# Stage C full run on Archmagi (WSL): build, fidelity, timing on three runtimes, stability, fuzz.
#   bash csharp/run_stagec.sh            (from the MandrakeModels root)
# Scenes come from `python -m mandrakemodels.csharp_export` (see csharp/README.md). Output: $STAGEC/out.
set -uo pipefail
STAGEC=${STAGEC:-/mnt/d/ClaudeWork/mandrakemodels/stagec}
W=$(wslpath -w "$STAGEC")
DOTNET=${DOTNET:-/mnt/c/Users/Mandrake/.dotnet/dotnet.exe}
rsync -a --delete "$(dirname "$0")/" "$STAGEC/src/"
( cd "$STAGEC/src/LiTiming" && "$DOTNET" build -c Release --nologo -v quiet ) || exit 1
EXE8="$W\\src\\LiTiming\\bin\\Release\\net8.0\\LiTiming.dll"
EXE4="$W\\src\\LiTiming\\bin\\Release\\net472\\LiTiming.exe"
HOST="$W\\src\\MonoHost\\monohost.py"
mono() { python.exe "$HOST" "$EXE4" "$@"; }
net8() { "$DOTNET" "$EXE8" "$@"; }
netfx() { "$STAGEC/src/LiTiming/bin/Release/net472/LiTiming.exe" "$@"; }
O="$STAGEC/out"; WO="$W\\out"
mkdir -p "$O"
step() { echo "$(date '+%H:%M') $*" >> "$O/run_stagec.milestones"; }
step start
mono fidelity "$W\\scenes" "$WO\\fidelity_mono.jsonl" > "$O/fidelity_mono.log" 2>&1; step "fidelity mono: $(tail -1 "$O/fidelity_mono.log")"
net8 fidelity "$W\\scenes" "$WO\\fidelity_net8.jsonl" > "$O/fidelity_net8.log" 2>&1; step "fidelity net8: $(tail -1 "$O/fidelity_net8.log")"
rm -f "$O/timing_final.csv"
for rep in 1 2; do
  mono timing "$W\\load" "$WO\\timing_final.csv" --label "unity_mono_6.13 r$rep" > "$O/timing_mono_r$rep.log" 2>&1; step "timing mono r$rep"
done
net8 timing "$W\\load" "$WO\\timing_final.csv" --label "net8 r1" > "$O/timing_net8_r1.log" 2>&1; step "timing net8"
netfx timing "$W\\load" "$WO\\timing_final.csv" --label "netfx48 r1" > "$O/timing_netfx_r1.log" 2>&1; step "timing netfx"
rm -f "$O/stability_mono.jsonl" "$O/stability_net8.jsonl"
mono stability "$W\\load" "$WO\\stability_mono.jsonl" --only lake_at_rest > "$O/stability_mono.log" 2>&1
mono stability "$W\\load" "$WO\\stability_mono.jsonl" --only river > "$O/stability_mono_river.log" 2>&1; step "stability mono (lake, river)"
net8 stability "$W\\load" "$WO\\stability_net8.jsonl" > "$O/stability_net8.log" 2>&1; step "stability net8: $(tail -1 "$O/stability_net8.log")"
mono fuzz "$W\\scenes\\river_128.lis" "$WO\\fuzz_mono.jsonl" --seqs 1000 > "$O/fuzz_mono.log" 2>&1; step "fuzz mono: $(tail -1 "$O/fuzz_mono.log")"
net8 fuzz "$W\\scenes\\pit_dig_128.lis" "$WO\\fuzz_net8_pit.jsonl" --seqs 1000 > "$O/fuzz_net8_pit.log" 2>&1; step "fuzz net8 pit: $(tail -1 "$O/fuzz_net8_pit.log")"
step done
