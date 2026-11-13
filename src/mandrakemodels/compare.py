"""Side-by-side comparison: one initial grid, every applicable solver, one simulated-time clock.

    python -m mandrakemodels.compare --scene dam_break --n 128 [--solvers a,b,c] [--ref swe_hll]
    python -m mandrakemodels.compare --scene all --n 96            # every scene, one HTML with a picker
    python -m mandrakemodels.compare --bed bed.png --bed-range 0,4 --level 1.5 --duration 60
    python -m mandrakemodels.compare --init my_grid.npz            # z, h, wall arrays (from a past run)
    python -m mandrakemodels.compare --html-only                    # rebuild the HTML from saved runs

The scene is built ONCE; every solver gets a deep copy of the same Scene (same bed, water, walls, forcing,
events, seed), so any difference between tiles is the solver. All solvers are sampled at the same
simulated times (water, wave, scalar); each takes its own stable dt between samples. Flow scenes are
sampled at the same step counts instead, because the two flow solvers do not share a time unit (lattice
steps vs seconds): their tiles compare shape, not speed.

Each tile shows the solver, the time, its volume error (relative; water and scalar) and its mean
ms/step on this machine (numpy, single process, sequential: comparable within a run). A solver that
blows up, throws, or runs out of its wall budget keeps its last good frame and is labelled in red; it is
never dropped. The difference view shows ``solver field - reference field`` (water: depth, m) on one symmetric
scale per scene.

Outputs go to ``<out root>/compare/<scene>_<n>/``: ``state.gif``, ``diff.gif``, ``state.mp4`` (if ffmpeg
is on PATH), ``init.npz`` (the exact starting grid, editable and re-loadable with ``--init``) and
``meta.json``. ``<out root>/compare/compare.html`` embeds every saved run as WebP frames and works from
file:// with no server (play/pause, time slider, state/difference toggle, scene picker).
"""
from __future__ import annotations

import argparse
import base64
import copy
import io
import json
import math
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
from matplotlib import colormaps
from PIL import Image, ImageDraw, ImageFont

from . import render, scenes
from .base import Grid
from .paths import out_root
from .solvers import REGISTRY

REF = {"water": "swe_hll", "scalar": "heat_explicit", "wave": "linear_wave", "flow": "stable_fluids"}
WATER_DT_CAP = 5.0   # s; implicit solvers have no stability limit, so their dt is min(cap, sample gap)
CLASSES = dict(REGISTRY)   # ours, plus open-source adapters when --external is given
STATUS_COLOR = {"ok": (235, 235, 235), "blown up": (255, 80, 80), "error": (255, 80, 80),
                "stopped": (255, 170, 60)}


# ---------------------------------------------------------------------------------- initial grid

def _load_grid(path: str, rng_lo_hi: str | None, name: str) -> np.ndarray:
    p = Path(path)
    if p.suffix == ".npy":
        return np.load(p).astype(np.float64)
    im = np.asarray(Image.open(p).convert("F"), dtype=np.float64) / 255.0
    im = im[::-1]                                   # image row 0 is north; array row 0 is south
    if rng_lo_hi:
        lo, hi = (float(v) for v in rng_lo_hi.split(","))
        im = lo + (hi - lo) * im
    elif name != "wall":
        raise SystemExit(f"--{name} is a PNG: give --{name}-range lo,hi (metres for black,white)")
    return im


def custom_scene(a) -> scenes.Scene:
    """A water scene from user arrays: ``--init npz`` or ``--bed`` with ``--depth`` or ``--level``."""
    if a.init:
        d = np.load(a.init)
        z, h = d["z"].astype(np.float64), d["h"].astype(np.float64)
        wall = d["wall"].astype(bool) if "wall" in d.files else np.zeros(z.shape, bool)
    else:
        z = _load_grid(a.bed, a.bed_range, "bed")
        if a.depth:
            h = _load_grid(a.depth, a.depth_range, "depth")
        elif a.level is not None:
            h = np.maximum(a.level - z, 0.0)
        else:
            raise SystemExit("give --depth or --level with --bed")
        wall = _load_grid(a.wall, None, "wall") > 0.5 if a.wall else np.zeros(z.shape, bool)
    if not (z.shape == h.shape == wall.shape):
        raise SystemExit(f"shape mismatch: bed {z.shape}, depth {h.shape}, wall {wall.shape}")
    ny, nx = z.shape
    h = np.where(wall, 0.0, np.maximum(h, 0.0))
    name = a.name or "custom"
    return scenes.Scene(name, "water", Grid(ny, nx, a.dx), a.duration or 60.0,
                        {"z": z, "h": h, "wall": wall}, params={"manning": a.manning},
                        note="user grid" + (f" from {a.init}" if a.init else f" bed {a.bed}"))


# ---------------------------------------------------------------------------------- running

def _primary(kind, snap):
    key = {"water": "h", "wave": "height", "scalar": "field", "flow": "dye"}[kind]
    return np.asarray(snap[key], dtype=np.float32)


def _total(s, kind):
    if kind == "water":
        return None
    if kind == "scalar" and hasattr(s, "total"):
        return s.total()
    return None


def run_one(cls, sc, samples, max_wall, store_snaps, style):
    """Advance one solver through the common sample points. Returns per-frame records."""
    s = cls.from_scene(copy.deepcopy(sc))
    events = sorted(copy.deepcopy(sc.events), key=lambda e: e[0])
    flow = sc.kind == "flow"
    tot0 = _total(s, sc.kind)
    frames, status, note = [], "ok", ""
    step_wall, t_start = 0.0, time.perf_counter()

    def record():
        snap = s.snapshot()
        r = {"t": s.t, "steps": s.nsteps, "status": status, "note": note,
             "ms": 1000 * step_wall / s.nsteps if s.nsteps else None, "field": _primary(sc.kind, snap)}
        if sc.kind == "water":
            r["vol_err"] = abs(s.volume_error()) / max(s.V0, s.ledger.throughput, 1e-12)
        elif tot0:
            r["vol_err"] = abs(_total(s, sc.kind) - tot0) / abs(tot0)
        else:
            r["vol_err"] = None
        if store_snaps:
            r["snap"] = {k: (np.asarray(v, np.float32) if isinstance(v, np.ndarray) and v.dtype != bool else v)
                         for k, v in snap.items()}
        else:
            r["rgb"] = render.frame(snap, style, sc.grid.dx).astype(np.float32)
        return r

    for target in samples:
        while status == "ok":
            done = (s.nsteps >= target) if flow else (s.t >= target - 1e-9)
            if done:
                break
            while events and s.t >= events[0][0] - 1e-9:
                _, meth, kw = events.pop(0)
                getattr(s, meth)(**kw)
            try:
                if flow:
                    dt = s.stable_dt()
                else:
                    dt = min(s.stable_dt(), target - s.t)
                    if sc.kind == "water":
                        dt = min(dt, WATER_DT_CAP)
                    if events and events[0][0] - s.t > 1e-9:
                        dt = min(dt, events[0][0] - s.t)
                    dt = max(dt, 1e-6)
                a = time.perf_counter()
                s.step(dt)
                step_wall += time.perf_counter() - a
            except Exception as e:  # reported, not hidden
                status, note = "error", f"{type(e).__name__}: {str(e)[:60]}"
                break
            if s.blown_up() or (s.nsteps % 25 == 0 and not np.isfinite(_primary(sc.kind, s.snapshot())).all()):
                status, note = "blown up", f"at t={s.t:.3g}"
                break
            if time.perf_counter() - t_start > max_wall:
                status, note = "stopped", f"wall cap {max_wall:.0f}s at t={s.t:.3g}"
                break
        if status != "ok" and frames and frames[-1]["status"] == "ok":
            frames.append(record())      # one frame showing the failed state
        elif status != "ok" and frames:
            prev = dict(frames[-1]); prev.update(status=status, note=note)
            frames.append(prev)
        else:
            frames.append(record())
    while len(frames) < len(samples):
        prev = dict(frames[-1]); prev.update(status=status, note=note)
        frames.append(prev)
    return frames


def run_scene(sc, names, ref, n_frames, max_wall):
    flow = sc.kind == "flow"
    if flow:
        total = int(sc.params.get("steps", 200))
        samples = np.linspace(0, total, n_frames + 1).round().astype(int).tolist()
    else:
        samples = np.linspace(0, sc.duration, n_frames + 1).tolist()
    order = [ref] + [n for n in names if n != ref] if ref in names else list(names)
    style = None
    if sc.kind == "water":
        probe = CLASSES[order[0]].from_scene(copy.deepcopy(sc))
        style = render.Style.from_snapshot("water", probe.snapshot())
    store = sc.kind != "water"
    runs = {}
    for name in order:
        a = time.perf_counter()
        runs[name] = run_one(CLASSES[name], sc, samples, max_wall, store, style)
        last = runs[name][-1]
        print(f"  {name:28s} t={last['t']:9.3f} steps={last['steps']:7d} {last['status']:9s} "
              f"{time.perf_counter() - a:6.1f}s {last['note']}", flush=True)
    if store:                                   # global style from the reference's whole run
        snaps = [f["snap"] for f in runs[order[0]]]
        st = _global_style(sc.kind, snaps)
        for name in order:
            for f in runs[name]:
                f["rgb"] = render.frame(f["snap"], st, sc.grid.dx).astype(np.float32)
                del f["snap"]
    return order, samples, runs


def _global_style(kind, snaps):
    if kind == "wave":
        st = render.Style("wave")
        st.extra["amp"] = max(float(np.abs(s["height"]).max()) for s in snaps) or 1.0
        return st
    if kind == "scalar":
        lo = min(float(s["field"].min()) for s in snaps)
        hi = max(float(s["field"].max()) for s in snaps)
        return render.Style("scalar", lo=lo, hi=hi + 1e-12, cmap="inferno")
    hi = max(float(np.percentile(s["dye"], 99.5)) for s in snaps)
    return render.Style("flow", lo=0.0, hi=max(hi, 1e-9), cmap="magma")


# ---------------------------------------------------------------------------------- drawing

def _font(px):
    try:
        return ImageFont.load_default(size=px)
    except TypeError:  # Pillow < 10.1
        return ImageFont.load_default()


def _tile(rgb, tile_px, lines, color):
    im = render.to_image(np.nan_to_num(rgb, nan=1.0), max_px=tile_px)
    head = 34
    out = Image.new("RGB", (tile_px, tile_px + head), (24, 24, 24))
    out.paste(im, ((tile_px - im.width) // 2, head + (tile_px - im.height) // 2))
    d = ImageDraw.Draw(out)
    f = _font(12)
    d.text((4, 2), lines[0], fill=(255, 255, 255), font=f)
    d.text((4, 18), lines[1], fill=color, font=f)
    if color != STATUS_COLOR["ok"]:
        d.rectangle([0, 0, tile_px - 1, tile_px + head - 1], outline=color, width=3)
    return out


def _sheet(tiles, cols, title):
    w, h = tiles[0].size
    rows = math.ceil(len(tiles) / cols)
    top = 22
    out = Image.new("RGB", (cols * w + (cols - 1) * 4, top + rows * (h + 4)), (12, 12, 12))
    d = ImageDraw.Draw(out)
    d.text((4, 4), title, fill=(255, 255, 255), font=_font(13))
    for k, t in enumerate(tiles):
        out.paste(t, ((k % cols) * (w + 4), top + (k // cols) * (h + 4)))
    return out


def _fmt(f, flow):
    t = f"step {f['steps']}" if flow else f"t={f['t']:.1f}s"
    ve = "-" if f["vol_err"] is None else f"{f['vol_err']:.0e}"
    ms = "-" if f["ms"] is None else f"{f['ms']:.2f}"
    return f"{t}  vol {ve}  {ms} ms/step"


def compose(sc, order, samples, runs, tile_px):
    flow = sc.kind == "flow"
    cols = min(len(order), 5 if len(order) > 8 else 4)
    ref = order[0]
    diffs = {}
    for name in order[1:]:
        diffs[name] = [f["field"] - rf["field"] for f, rf in zip(runs[name], runs[ref])]
    allabs = np.concatenate([np.abs(d[np.isfinite(d)]).ravel() for v in diffs.values() for d in v] or [np.zeros(1)])
    scale = float(np.percentile(allabs, 99.5)) if allabs.size else 1.0
    scale = scale if scale > 0 else 1.0
    cmap = colormaps["RdBu_r"]
    state_frames, diff_frames = [], []
    unit = {"water": "m depth", "wave": "height", "scalar": "value", "flow": "dye"}[sc.kind]
    for k in range(len(samples)):
        clock = f"step {samples[k]}" if flow else f"t = {samples[k]:.1f} s"
        tiles, dtiles = [], []
        for name in order:
            f = runs[name][k]
            color = STATUS_COLOR[f["status"]]
            l2 = _fmt(f, flow) if f["status"] == "ok" else f"{f['status'].upper()} {f['note']}"
            tiles.append(_tile(f["rgb"], tile_px, [name + ("  (ref)" if name == ref else ""), l2], color))
            if name == ref:
                rgb = f["rgb"]
                dtiles.append(_tile(rgb, tile_px, [name + "  (reference)", _fmt(f, flow)], STATUS_COLOR["ok"]))
            else:
                d = diffs[name][k]
                img = cmap(np.clip(0.5 + 0.5 * np.nan_to_num(d / scale), 0, 1))[..., :3][::-1]
                rms = float(np.sqrt(np.nanmean(d.astype(np.float64) ** 2)))
                l2d = f"rms diff {rms:.3g}" if f["status"] == "ok" else f"{f['status'].upper()} {f['note']}"
                dtiles.append(_tile(img, tile_px, [name + " minus ref", l2d], color))
        head = f"{sc.name}  n={sc.grid.nx}x{sc.grid.ny}  {clock}"
        state_frames.append(_sheet(tiles, cols, head + "  (same colour scale on every tile)"))
        diff_frames.append(_sheet(dtiles, cols, head + f"  solver minus reference, {unit}, "
                                                       f"blue/red = +/-{scale:.3g}"))
    return state_frames, diff_frames, scale


# ---------------------------------------------------------------------------------- outputs

def rms_vs_ref(ref_frames, frames):
    """RMS of (solver - reference) primary field over the frames both reached, in the field's units
    (water: depth in m). None when the solver never got past the start."""
    k = min(len(ref_frames), len(frames))
    sq = [np.mean((np.asarray(frames[i]["field"], float) - np.asarray(ref_frames[i]["field"], float)) ** 2)
          for i in range(1, k) if frames[i]["status"] == "ok" or i < k - 1]
    sq = [v for v in sq if np.isfinite(v)]
    return float(np.sqrt(np.mean(sq))) if sq else None


def save_outputs(sc, order, samples, runs, state_frames, diff_frames, scale, d: Path, fps=8):
    d.mkdir(parents=True, exist_ok=True)
    fr = d / "frames"
    if fr.exists():
        shutil.rmtree(fr)
    fr.mkdir()
    for k, (a, b) in enumerate(zip(state_frames, diff_frames)):
        a.save(fr / f"state_{k:03d}.webp", quality=72, method=4)
        b.save(fr / f"diff_{k:03d}.webp", quality=72, method=4)
        a.save(fr / f"state_{k:03d}.png")
    for nm, frames in (("state", state_frames), ("diff", diff_frames)):
        ims = [f.convert("P", palette=Image.ADAPTIVE, colors=192) for f in frames]
        ims[0].save(d / f"{nm}.gif", save_all=True, append_images=ims[1:] + [ims[-1]] * fps,
                    duration=int(1000 / fps), loop=0, optimize=True)
    mp4 = None
    if shutil.which("ffmpeg"):
        mp4 = d / "state.mp4"
        r = subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i",
                            str(fr / "state_%03d.png"), "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                            "-pix_fmt", "yuv420p", "-c:v", "libx264", "-crf", "23", str(mp4)],
                           capture_output=True, text=True)
        if r.returncode:
            print("ffmpeg failed:", r.stderr[-300:])
            mp4 = None
    for p in fr.glob("state_*.png"):
        p.unlink()
    final = {n: {k: runs[n][-1][k] for k in ("t", "steps", "status", "note", "ms", "vol_err")} for n in order}
    for n in order:
        final[n]["rms_vs_ref"] = rms_vs_ref(runs[order[0]], runs[n])
    meta = {"scene": sc.name, "kind": sc.kind, "n": [sc.grid.ny, sc.grid.nx], "dx": sc.grid.dx,
            "duration": sc.duration, "note": sc.note, "ref": order[0], "solvers": order,
            "samples": samples, "diff_scale": scale, "frames": len(state_frames),
            "final": final,
            "made": time.strftime("%Y-%m-%dT%H:%M:%S")}
    (d / "meta.json").write_text(json.dumps(meta, indent=1, default=float))
    return mp4


HTML_WIDTH, HTML_QUALITY, HTML_FRAMES = 760, 45, 20   # keeps the 14-scene player under ~10 MB


def _html_enc(p: Path) -> str:
    """Frame re-encoded small for the player; the full-size frames stay on disk next to it."""
    import io
    im = Image.open(p)
    if im.width > HTML_WIDTH:
        im = im.resize((HTML_WIDTH, round(im.height * HTML_WIDTH / im.width)), Image.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, "WEBP", quality=HTML_QUALITY, method=4)
    return "data:image/webp;base64," + base64.b64encode(buf.getvalue()).decode()


def build_html(root: Path) -> Path:
    runs = []
    for m in sorted(root.glob("*/meta.json")):
        meta = json.loads(m.read_text())
        fr = m.parent / "frames"
        st = sorted(fr.glob("state_*.webp"))
        df = sorted(fr.glob("diff_*.webp"))
        if not st:
            continue
        keep = sorted({round(i * (len(st) - 1) / max(1, HTML_FRAMES - 1)) for i in range(min(len(st), HTML_FRAMES))})
        samples = meta["samples"]
        runs.append({"id": m.parent.name, "label": f"{meta['scene']} {meta['n'][1]}x{meta['n'][0]}",
                     "note": meta.get("note", ""), "kind": meta["kind"],
                     "samples": [samples[i] for i in keep if i < len(samples)],
                     "state": [_html_enc(st[i]) for i in keep], "diff": [_html_enc(df[i]) for i in keep if i < len(df)],
                     "final": meta["final"]})
    page = HTML.replace("__DATA__", json.dumps(runs))
    out = root / "compare.html"
    out.write_text(page)
    return out


HTML = """<!doctype html><html><head><meta charset="utf-8"><title>MandrakeModels: same grid, every solver</title>
<style>body{background:#111;color:#ddd;font:14px system-ui,sans-serif;margin:12px}
select,button{font-size:14px;margin-right:8px}#img{max-width:100%;display:block;margin-top:8px}
#bar{display:flex;align-items:center;gap:8px;flex-wrap:wrap}input[type=range]{width:360px}
table{border-collapse:collapse;margin-top:8px;font-size:12px}td,th{border:1px solid #333;padding:2px 6px}
.bad{color:#f66}.stop{color:#fa4}</style></head><body>
<div id="bar"><select id="scene"></select><button id="play">pause</button>
<label><input type="radio" name="v" value="state" checked>state</label>
<label><input type="radio" name="v" value="diff">difference from reference</label>
<input type="range" id="t" min="0" value="0"><span id="lab"></span>
<label>fps <input type="number" id="fps" value="8" min="1" max="30" style="width:3em"></label></div>
<div id="note"></div><img id="img"><table id="tab"></table>
<script>
const R=__DATA__;let cur=0,k=0,playing=true,view='state',timer=null;
const $=id=>document.getElementById(id);
R.forEach((r,i)=>{const o=document.createElement('option');o.value=i;o.textContent=r.label;$('scene').appendChild(o)});
function table(){const r=R[cur];let h='<tr><th>solver</th><th>status</th><th>end t</th><th>steps</th><th>ms/step</th><th>vol err</th><th>rms vs ref</th></tr>';
for(const [n,f] of Object.entries(r.final)){const c=f.status=='ok'?'':(f.status=='stopped'?'stop':'bad');
h+=`<tr class="${c}"><td>${n}</td><td>${f.status} ${f.note||''}</td><td>${(+f.t).toFixed(2)}</td><td>${f.steps}</td>`+
`<td>${f.ms==null?'-':(+f.ms).toFixed(3)}</td><td>${f.vol_err==null?'-':(+f.vol_err).toExponential(1)}</td>`+
`<td>${f.rms_vs_ref==null?'-':(+f.rms_vs_ref).toPrecision(3)}</td></tr>`}
$('tab').innerHTML=h}
function show(){const r=R[cur];const fr=r[view];k=Math.min(k,fr.length-1);$('img').src=fr[k];$('t').value=k;
$('lab').textContent=(r.kind=='flow'?'step ':'t = ')+r.samples[k]+(r.kind=='flow'?'':' s')}
function pick(i){cur=+i;k=0;const r=R[cur];$('t').max=r.state.length-1;$('note').textContent=r.note;table();show()}
function tick(){if(playing){k=(k+1)%R[cur][view].length;show()}}
function loop(){clearInterval(timer);timer=setInterval(tick,1000/Math.max(1,+$('fps').value))}
$('scene').onchange=e=>pick(e.target.value);$('t').oninput=e=>{k=+e.target.value;show()};
$('play').onclick=()=>{playing=!playing;$('play').textContent=playing?'pause':'play'};$('fps').onchange=loop;
document.querySelectorAll('input[name=v]').forEach(x=>x.onchange=e=>{view=e.target.value;show()});
if(R.length){pick(0);loop()}else{$('note').textContent='no saved runs'}
</script></body></html>"""


# ---------------------------------------------------------------------------------- cli

def applicable(sc, only=None, exclude=()):
    def fits(c):
        return c.kind == sc.kind and getattr(c, "applies", lambda _sc: True)(sc)
    names = [n for n, c in CLASSES.items() if fits(c)]
    if only:
        bad = [n for n in only if n not in CLASSES]
        if bad:
            raise SystemExit(f"unknown solvers {bad}; known: {sorted(CLASSES)}")
        names = [n for n in only if fits(CLASSES[n])] + [n for n in names if n.startswith("ext_")
                                                          and n not in only]
    return [n for n in names if n not in exclude]


def main(argv=None):
    global WATER_DT_CAP
    p = argparse.ArgumentParser(prog="mandrakemodels.compare", description=__doc__.split("\n")[0])
    p.add_argument("--scene", default="dam_break", help="scene name, comma list, 'all' or 'water'")
    p.add_argument("--n", type=int, default=128)
    p.add_argument("--solvers", default="", help="comma list (default: every solver of the scene's kind)")
    p.add_argument("--exclude", default="", help="comma list of solvers to leave out")
    p.add_argument("--ref", default="", help="reference solver for the difference view")
    p.add_argument("--frames", type=int, default=40)
    p.add_argument("--duration", type=float, default=None, help="override the scene's simulated seconds")
    p.add_argument("--max-wall", type=float, default=60.0, help="wall seconds per solver per scene")
    p.add_argument("--tile", type=int, default=200, help="tile size in pixels")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="", help="output root (default <agent storage>/mandrakemodels/compare)")
    p.add_argument("--html-only", action="store_true")
    p.add_argument("--dt-cap", type=float, default=WATER_DT_CAP,
                   help="largest water dt in s (matters for the implicit solvers, which have no limit)")
    p.add_argument("--external", action="store_true",
                   help="add tiles from open-source libraries that import (see external.py)")
    g = p.add_argument_group("your own starting grid (water)")
    g.add_argument("--init", help=".npz with z, h[, wall] (every run saves one as init.npz)")
    g.add_argument("--bed", help="bed elevation .npy (m) or .png (needs --bed-range)")
    g.add_argument("--bed-range", help="lo,hi metres for black,white of a PNG bed")
    g.add_argument("--depth", help="water depth .npy (m) or .png (needs --depth-range)")
    g.add_argument("--depth-range")
    g.add_argument("--level", type=float, help="fill to this surface level instead of --depth")
    g.add_argument("--wall", help="wall mask .npy or .png (white = wall)")
    g.add_argument("--dx", type=float, default=scenes.DX)
    g.add_argument("--manning", type=float, default=0.03)
    g.add_argument("--name", default="")
    a = p.parse_args(argv)
    WATER_DT_CAP = a.dt_cap
    if a.external:
        from . import external
        got, notes = external.available()
        CLASSES.update(got)
        for k, v in notes.items():
            print(f"external {k}: not available ({v})")
    root = Path(a.out) if a.out else out_root() / "compare"
    root.mkdir(parents=True, exist_ok=True)
    if not a.html_only:
        if a.init or a.bed:
            built = [custom_scene(a)]
        else:
            names = {"all": scenes.ALL, "water": scenes.WATER}.get(a.scene, a.scene.split(","))
            built = []
            for nm in names:
                sc = scenes.build(nm, a.n, a.seed)
                if a.duration:
                    sc.duration = a.duration
                built.append(sc)
        only = [s for s in a.solvers.split(",") if s]
        exclude = [s for s in a.exclude.split(",") if s]
        for sc in built:
            names = applicable(sc, only, exclude)
            if not names:
                print(f"{sc.name}: no applicable solver, skipped")
                continue
            ref = a.ref if a.ref in names else (REF[sc.kind] if REF[sc.kind] in names else names[0])
            print(f"{sc.name} ({sc.kind}, {sc.grid.nx}x{sc.grid.ny}): {len(names)} solvers, ref {ref}",
                  flush=True)
            order, samples, runs = run_scene(sc, names, ref, a.frames, a.max_wall)
            st, df, scale = compose(sc, order, samples, runs, a.tile)
            d = root / f"{sc.name}_{sc.grid.nx}"
            if sc.kind == "water":
                np.savez_compressed(d.mkdir(parents=True, exist_ok=True) or d / "init.npz",
                                    z=sc.init["z"], h=sc.init["h"], wall=sc.init["wall"])
            mp4 = save_outputs(sc, order, samples, runs, st, df, scale, d)
            print(f"  -> {d}/state.gif{' + state.mp4' if mp4 else ''}", flush=True)
    html = build_html(root)
    print(f"html: {html} ({html.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
