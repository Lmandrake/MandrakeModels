"""Turn stage C output into results/stagec.md tables.

    python3 csharp/summarize_stagec.py $STAGEC/out > results/stagec_tables.md

Unity-Mono rows: the WORSE of the repeats (higher mean ms/tick) is reported, and the spread between repeats is shown.
"""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

out = Path(sys.argv[1])
rows = list(csv.DictReader(open(out / "timing_final.csv")))
by = defaultdict(list)
for r in rows:
    rt = r["runtime"].split(" r")[0]
    by[(rt, r["scene"], int(r["n"]), r["mode"])].append(r)


def pick(rt, scene, n, mode):
    rs = by.get((rt, scene, n, mode), [])
    if not rs:
        return None, 0.0
    worst = max(rs, key=lambda r: float(r["mean_ms_tick"]))
    means = [float(r["mean_ms_tick"]) for r in rs]
    spread = (max(means) - min(means)) / max(means) if len(means) > 1 else 0.0
    return worst, spread


MONO = "unity_mono_6.13"
SCENES = ["census_river", "river", "pit_dig", "lake_at_rest", "dam_break", "swell", "stress_half"]
print("### Unity Mono 6.13 (RimWorld's own runtime), wet-only stepping\n")
print("ms per tick = the 12-tick recompute's cost spread over its 12 ticks; the recompute itself lands in one tick unless split.\n")
print("| scene | n | wet squares | substeps per recompute | ms/tick mean | median | p95 | ms per recompute median | p95 | max | repeat spread | books dL |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
for sc in SCENES:
    for n in (250, 325, 400):
        r, sp = pick(MONO, sc, n, "wet")
        if r is None:
            continue
        sub = int(r["substeps"]) / int(r["recomputes"])
        print(f"| {sc} | {n} | {r['mean_wet_cells']} | {sub:.1f} | {float(r['mean_ms_tick']):.3f} | {float(r['med_ms_tick']):.3f} | "
              f"{float(r['p95_ms_tick']):.3f} | {float(r['med_ms_recompute']):.2f} | {float(r['p95_ms_recompute']):.2f} | "
              f"{float(r['max_ms_recompute']):.2f} | {sp:.0%} | {r['books_dl']} |")

print("\n### Full grid against wet-only, and the three runtimes (mean ms per tick)\n")
print("| scene | n | Mono full | Mono wet | .NET Fx 4.8 wet | .NET 8 wet | Mono / .NET 8 |")
print("|---|---|---|---|---|---|---|")
for sc in SCENES:
    for n in (250, 325, 400):
        mf, _ = pick(MONO, sc, n, "full")
        mw, _ = pick(MONO, sc, n, "wet")
        fx, _ = pick("netfx48", sc, n, "wet")
        n8, _ = pick("net8", sc, n, "wet")
        if mw is None:
            continue
        g = lambda r: f"{float(r['mean_ms_tick']):.3f}" if r else "n/a"
        ratio = float(mw["mean_ms_tick"]) / float(n8["mean_ms_tick"]) if n8 else float("nan")
        print(f"| {sc} | {n} | {g(mf)} | {g(mw)} | {g(fx)} | {g(n8)} | {ratio:.1f} |")

for name in ("fidelity_mono.jsonl", "fidelity_net8.jsonl"):
    p = out / name
    if p.exists():
        rs = [json.loads(l) for l in open(p)]
        ok = sum(r["pass"] for r in rs)
        worst = max(r["l1rel"] for r in rs)
        print(f"\n{name}: {ok}/{len(rs)} rows pass; worst relative L1 {worst:.2e}; min IoU {min(r['iou'] for r in rs):.4f}; "
              f"books {sorted(set(r['books_dl'] for r in rs))}; wet-only identical {all(r['wet_only_identical'] for r in rs)}")
for name in ("stability_mono.jsonl", "stability_net8.jsonl"):
    p = out / name
    if p.exists():
        print(f"\n{name}:")
        for l in open(p):
            r = json.loads(l)
            extra = f", lake surface moved {r['max_surface_move_m'] * 1000:.2f} mm" if "max_surface_move_m" in r else ""
            print(f"- {r['run']}: {r['water_s']:.0f} water-s, {r['recomputes']} recomputes, {r['substeps']} steps, worst books "
                  f"{r['worst_books_dl']} dL, min square {r['min_v']} dL, max depth {r['max_depth_m']:.2f} m{extra}: "
                  f"{'PASS' if r['pass'] else 'FAIL'}")
for name in ("fuzz_mono.jsonl", "fuzz_net8_pit.jsonl"):
    p = out / name
    if p.exists():
        last = json.loads(open(p).read().splitlines()[-1])
        print(f"\n{name}: {last['seqs']} sequences, {last['edits']} edits, {last['recompute_checks']} recompute checks, "
              f"{last['fails']} failures ({last['runtime']})")
