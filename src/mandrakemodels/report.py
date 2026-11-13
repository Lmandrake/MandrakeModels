"""Merge the bakeoff phases into results.csv, results.md, contact sheets and gallery.html."""
from __future__ import annotations

import csv
import html
import json
import math
import shutil
import statistics
from pathlib import Path

from PIL import Image, ImageDraw

from . import scenes

STRESS_TICK_MS = 2.0


def _load(p: Path):
    if not p.exists():
        return []
    out = []
    for line in p.read_text().splitlines():
        try:
            out.append(json.loads(line))
        except Exception:
            pass
    return out


def _latest(rows):
    d = {}
    for r in rows:
        d[(r.get("solver"), r.get("scene"), r.get("n"))] = r
    return d


def _f(v, fmt="{:.3g}"):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return ""
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return fmt.format(v)
    return str(v)


def merged(out: Path):
    q = _latest(_load(out / "quality.jsonl"))
    t = _latest(_load(out / "timing.jsonl"))
    md = {r["solver"]: r for r in _load(out / "maxdt.jsonl")}
    probe = json.loads((out / "probe.json").read_text()) if (out / "probe.json").exists() else {}
    keys = sorted(set(q) | set(t), key=lambda k: (k[1], k[0], k[2]))
    rows = []
    for k in keys:
        a, b = q.get(k, {}), t.get(k, {})
        solver, scene, n = k
        r = dict(solver=solver, scene=scene, n=n, kind=a.get("kind") or b.get("kind"))
        r["ms_per_step"] = b.get("ms_per_step")
        r["dt_mean"] = b.get("dt_mean")
        r["py_ms_per_tick"] = b.get("py_ms_per_tick")
        # No C# projection from the numba probe: its kernels are simplified solvers (pipes without friction,
        # local_inertial without theta weighting), so their numpy/numba ratio prices different arithmetic.
        # Measured C# cost: results/stagec.md.
        for f in ("final_rate", "steps", "t_end", "capped", "blew_up", "blow_t", "error", "settled", "settle_t",
                  "settle_steps", "vol_err_rel", "min_h", "max_h", "wet_frac", "still_drift_m",
                  "max_speed", "gif"):
            r[f] = a.get(f)
        if not r["error"]:
            r["error"] = b.get("error")
        if b.get("blew_up"):
            r["blew_up_timing"] = True
        r["max_dt_dam128"] = md.get(solver, {}).get("max_dt")
        r["diag"] = a.get("diag", {})
        rows.append(r)
    return rows, md, probe


CSV_FIELDS = ["solver", "scene", "n", "kind", "ms_per_step", "dt_mean", "py_ms_per_tick", "steps", "t_end", "capped", "blew_up", "blow_t", "error",
              "settled", "settle_t", "settle_steps", "final_rate", "vol_err_rel", "min_h", "max_h", "wet_frac",
              "still_drift_m", "max_speed", "max_dt_dam128"]


def headline(rows, n):
    """Per water solver at size n: median over water scenes."""
    by = {}
    for r in rows:
        if r["n"] == n and r["kind"] == "water":
            by.setdefault(r["solver"], []).append(r)
    out = []
    for s, rs in sorted(by.items()):
        def med(f):
            v = [r[f] for r in rs if r.get(f) is not None]
            return statistics.median(v) if v else None
        def mx(f):
            v = [r[f] for r in rs if r.get(f) is not None]
            return max(v) if v else None
        out.append(dict(
            solver=s, scenes=len(rs), ms_per_step=med("ms_per_step"), dt=med("dt_mean"),
            py_tick=med("py_ms_per_tick"), vol=mx("vol_err_rel"),
            blow=sum(1 for r in rs if r.get("blew_up") or r.get("blew_up_timing")),
            err=sum(1 for r in rs if r.get("error")),
            settled=sum(1 for r in rs if r.get("settled")),
            capped=sum(1 for r in rs if r.get("capped")),
            ran=sum(1 for r in rs if r.get("steps") is not None),
            drift=next((r["still_drift_m"] for r in rs if r.get("still_drift_m") is not None), None)))
    return out


def md_headline(h, n):
    L = [f"### Water solvers at {n} x {n} (median over water scenes)", "",
         "| solver | ms/step numpy | dt s | numpy ms/tick | worst vol err | "
         "lake drift m | settled | blow-ups | errors | capped |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for r in h:
        L.append("| {solver} | {ms} | {dt} | {pt} | {vol} | {drift} | {st}/{ran} | {bl} | {er} | {cp} |".format(
            solver=r["solver"], ms=_f(r["ms_per_step"]), dt=_f(r["dt"]), pt=_f(r["py_tick"]), vol=_f(r["vol"], "{:.1e}"),
            drift=_f(r["drift"], "{:.1e}"), st=r["settled"], ran=r["ran"], bl=r["blow"], er=r["err"],
            cp=r["capped"]))
    return "\n".join(L)


def md_scene(rows, n):
    L = [f"### Every run at {n} x {n}", "",
         "| scene | solver | ms/step | dt s | steps | sim s | settled at s | final max dh/dt m/s | vol err | min h | blew up | capped | error |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["n"] != n:
            continue
        L.append(f"| {r['scene']} | {r['solver']} | {_f(r['ms_per_step'])} | {_f(r['dt_mean'])} | {_f(r['steps'])} | "
                 f"{_f(r['t_end'])} | {_f(r['settle_t'])} | {_f(r.get('final_rate'), '{:.1e}')} | {_f(r['vol_err_rel'], '{:.1e}')} | {_f(r['min_h'], '{:.1e}')} | "
                 f"{_f(r['blew_up'])}{' (' + _f(r['blow_t']) + ' s)' if r.get('blow_t') else ''} | {_f(r['capped'])} | "
                 f"{html.escape(str(r['error']))[:60] if r['error'] else ''} |")
    return "\n".join(L)


def md_scaling(rows, sizes):
    by = {}
    for r in rows:
        if r.get("ms_per_step") is not None:
            by.setdefault((r["solver"], r["kind"]), {}).setdefault(r["n"], []).append(r["ms_per_step"])
    L = ["### numpy ms/step by size (median over scenes; timing phase, sequential, one thread)", "",
         "| solver | kind | " + " | ".join(str(s) for s in sizes) + " |", "|---|---|" + "---|" * len(sizes)]
    for (s, k), d in sorted(by.items(), key=lambda x: (x[0][1], x[0][0])):
        L.append(f"| {s} | {k} | " + " | ".join(_f(statistics.median(d[n])) if n in d else "" for n in sizes) + " |")
    return "\n".join(L)


def contact_sheets(out: Path, rows, n):
    """One PNG per scene: every solver's mid and final frame, labelled."""
    sheets = {}
    by = {}
    for r in rows:
        if r["n"] == n and r.get("gif"):
            by.setdefault(r["scene"], []).append(r)
    for scene, rs in by.items():
        rs.sort(key=lambda r: r["solver"])
        cell = 200
        W, H = cell * len(rs), 2 * cell + 34
        sheet = Image.new("RGB", (W, H), (24, 24, 24))
        d = ImageDraw.Draw(sheet)
        for i, r in enumerate(rs):
            stem = out / "frames" / scene / f"{r['solver']}_{n}"
            for row, k in enumerate(("b", "c")):
                p = Path(f"{stem}_{k}.png")
                if p.exists():
                    im = Image.open(p).convert("RGB")
                    im.thumbnail((cell, cell))
                    sheet.paste(im, (i * cell, 30 + row * cell))
            tag = r["solver"] + (" BLEW UP" if r.get("blew_up") else "")
            d.text((i * cell + 4, 4), tag, fill=(255, 120, 120) if r.get("blew_up") else (230, 230, 230))
        p = out / "sheets" / f"{scene}_{n}.png"
        p.parent.mkdir(exist_ok=True)
        sheet.save(p)
        sheets[scene] = p
    return sheets


def gallery(out: Path, rows, n, md_text):
    by = {}
    for r in rows:
        if r["n"] == n and r.get("gif"):
            by.setdefault(r["scene"], []).append(r)
    notes = {}
    for name in scenes.ALL:
        try:
            notes[name] = scenes.build(name, 32).note
        except Exception:
            notes[name] = ""
    H = ["<!doctype html><meta charset=utf-8><title>MandrakeModels fluid bakeoff</title>",
         "<style>body{background:#161616;color:#ddd;font:14px system-ui;margin:24px}"
         "h2{margin-top:36px;border-bottom:1px solid #444}.row{display:flex;flex-wrap:wrap;gap:14px}"
         ".card{width:330px}.card img{width:320px;image-rendering:pixelated;border:1px solid #333}"
         ".m{font:12px monospace;color:#aaa}.bad{color:#f77}table{border-collapse:collapse;font:12px monospace}"
         "td,th{border:1px solid #333;padding:2px 6px}</style>",
         f"<h1>Fluid solver bakeoff: gallery at {n} x {n}</h1>",
         "<p>Each card loops the whole scene duration in 48 frames. Numbers: numpy ms/step (timing phase), "
         "relative volume error, when it settled. Full tables in <a href=results.md>results.md</a> and "
         "<a href=results.csv>results.csv</a>.</p>"]
    order = [s for s in scenes.ALL if s in by]
    H.append("<p>" + " | ".join(f"<a href=#{s}>{s}</a>" for s in order) + "</p>")
    for s in order:
        H.append(f"<h2 id={s}>{s}</h2><p>{html.escape(notes.get(s, ''))}</p><div class=row>")
        for r in sorted(by[s], key=lambda r: r["solver"]):
            bad = r.get("blew_up") or r.get("error")
            m = (f"{_f(r.get('ms_per_step'))} ms/step, dt {_f(r.get('dt_mean'))} s<br>"
                 f"vol err {_f(r.get('vol_err_rel'), '{:.1e}')}, settled {_f(r.get('settle_t'))} s")
            if r.get("blew_up"):
                m += f"<br><b class=bad>BLEW UP at t={_f(r.get('blow_t'))} s</b>"
            if r.get("capped"):
                m += "<br>capped (wall budget)"
            H.append(f"<div class=card><b class={'bad' if bad else ''}>{r['solver']}</b><br>"
                     f"<img src='{r['gif']}' loading=lazy><div class=m>{m}</div></div>")
        H.append("</div>")
    (out / "gallery.html").write_text("\n".join(H))


def write(out: Path, curate: Path | None = None):
    rows, md, probe = merged(out)
    with (out / "results.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, CSV_FIELDS, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in CSV_FIELDS})
    sizes = sorted({r["n"] for r in rows})
    big = max((r["n"] for r in rows if r.get("ms_per_step") is not None), default=None)
    parts = ["# Fluid solver bakeoff: results", "",
             "Tags: every number here is [measured] by this harness. numpy ms/tick = numpy ms/step x steps per "
             "tick at the solver's own dt, with one simulated second per 60 ticks (sibling plan's assumption). "
             "No C# cost is projected from these: the numba probe's kernels are simplified solvers, so its ratio "
             "is not a C# projection; measured C# cost (local_inertial on Unity Mono) is in results/stagec.md. "
             "Full-grid cost: no sleeping or active sets.", ""]
    if probe:
        parts += ["### numba unit cost probe (simplified kernels: NOT the bakeoff's arithmetic, not a C# projection)", "",
                  "| kernel | numpy ms/step | numba ms/step | ratio | numba ns/cell |", "|---|---|---|---|---|"]
        for k, v in probe.items():
            if isinstance(v, dict):
                parts.append(f"| {k} | {v['numpy_ms']:.2f} | {v['numba_ms']:.2f} | {v['ratio']:.1f} | {v['numba_ns_per_cell']:.1f} |")
        parts += ["", f"Grid {probe.get('n')}^2, median of 30 steps, single thread. The pipes kernel has no "
                  "friction and the local_inertial kernel has no theta weighting or coupled-|q| friction "
                  "(probe_numba.py), so these ratios compare different solvers.", ""]
    for n in (325, big):
        if n in sizes:
            parts += [md_headline(headline(rows, n), n), ""]
    parts += [md_scaling(rows, sizes), ""]
    if md:
        parts += ["### Largest fixed dt that survives 20 s of the 128 dam break (water) / gaussian / ripples", "",
                  "Survives means no NaN, no blow-up and depth never above 1.25x its start (scalar: no new "
                  "extreme; wave: bounded) over the whole run. Bounded is not accurate: a limiter can keep a "
                  "solver bounded at a dt where its wave speeds are already wrong.", "",
                  "| solver | max dt s | solver's own stable_dt s | note |", "|---|---|---|---|"]
        for s, r in sorted(md.items()):
            parts.append(f"| {s} | {_f(r.get('max_dt'))} | {_f(r.get('own_stable_dt'))} | {r.get('note', '')} |")
        parts.append("")
    acc = json.loads((out / "accuracy.json").read_text()) if (out / "accuracy.json").exists() else None
    if acc:
        parts += [f"### Flood timing: water through the dam breach in {acc['T']:g} s ({acc['scene']}, {acc['n']}^2)", "",
                  f"Reference {acc['reference']} (matches Stoker's analytic dam break to 0.1%, tests/test_water_main.py). "
                  "Ratio well below 1: floods arrive late; above 1: early.", "",
                  "| solver | m^3 through breach | ratio to reference |", "|---|---|---|"]
        for k, v in sorted(acc["m3"].items(), key=lambda kv: -(acc["ratio"].get(kv[0]) or 0)):
            parts.append(f"| {k} | {_f(v, '{:.0f}')} | {_f(acc['ratio'].get(k), '{:.2f}')} |")
        parts.append("")
    frame_n = 250 if 250 in sizes else (sizes[0] if sizes else None)
    if frame_n:
        parts += [md_scene(rows, frame_n), ""]
    text = "\n".join(parts)
    (out / "results.md").write_text(text)
    sheets = contact_sheets(out, rows, frame_n) if frame_n else {}
    gallery(out, rows, frame_n, text)
    if curate:
        curate.mkdir(parents=True, exist_ok=True)
        for f in ("results.csv", "results.md", "probe.json", "maxdt.jsonl"):
            if (out / f).exists():
                shutil.copy(out / f, curate / f)
        sd = curate / "sheets"
        sd.mkdir(exist_ok=True)
        for scene, p in sheets.items():
            im = Image.open(p)
            im.thumbnail((1400, 500))
            im.save(sd / p.name, optimize=True)
    print(f"report: {len(rows)} rows -> {out / 'results.md'}, {out / 'gallery.html'}")
