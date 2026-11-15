"""Side-by-side comparison: one initial grid, every applicable solver, one simulated-time clock.

    python -m mandrakemodels.compare --scene dam_break --external          # one water scene, 128 cells
    python -m mandrakemodels.compare --scene water --external              # every water scene, one by one
    python -m mandrakemodels.compare --scene dam_break --n 250 --tile 500  # a bigger showcase run
    python -m mandrakemodels.compare --bed bed.png --bed-range 0,4 --level 1.5 --duration 60
    python -m mandrakemodels.compare --init my_grid.npz            # z, h, wall arrays (from a past run)
    python -m mandrakemodels.compare --html-only                    # rebuild the pages from saved runs

The scene is built ONCE; every solver gets a deep copy of the same Scene (same bed, water, walls, forcing,
events, seed), so any difference between tiles is the solver. All solvers are sampled at the same
simulated times (water, wave, scalar); each takes its own stable dt between samples. Flow scenes are
sampled at the same step counts instead, because the two flow solvers do not share a time unit (lattice
steps vs seconds): their tiles compare shape, not speed.

Water is drawn as water (``waterview``): shaded-relief terrain, water in blues (darker = deeper), the
water surface lit so wave fronts show as light and dark bands, a pale shoreline, dark walls; one look per
scene taken from the reference run, so every tile shares colours and scales. Dam break and swell add a
side view of the water surface along one grid row under each tile. Each tile names the solver in plain
words, then its code name, volume error and mean ms/step on this machine. A solver that blows up, throws,
or runs out of its wall budget keeps its last good frame under a red banner; it is never dropped. The
difference view shows solver depth minus reference depth (red = more water than the reference).

Outputs go to ``<out root>/compare/<scene>_<n>/``: ``state.mp4`` and ``state.webm``/``diff.webm``,
small ``state.gif``/``diff.gif``, three stills, ``init.npz`` (the exact starting grid, editable and
re-loadable with ``--init``) and ``meta.json``. ``<out root>/compare/index.html`` lists the scenes
(water first) and links one self-contained page per scene (``<scene>_<n>.html``, the videos embedded,
works from file:// with no server). Stills for checking also go to ``<out root>/compare/stills/``.
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
import tempfile
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

from . import render, scenes, waterview
from .base import Grid
from .paths import out_root
from .solvers import REGISTRY

REF = {"water": "swe_hll", "scalar": "heat_explicit", "wave": "linear_wave", "flow": "stable_fluids"}
WATER_DT_CAP = 5.0   # s; implicit solvers have no stability limit, so their dt is min(cap, sample gap)
CLASSES = dict(REGISTRY)   # ours, plus open-source adapters when --external is given
BAD = (220, 40, 40)
STOP = (225, 120, 20)

PLAIN = {
    "swe_hll": "shallow water equations (HLL finite volume)",
    "pipes": "virtual pipes",
    "pipes_theta": "virtual pipes, theta-smoothed (visual only)",
    "ca_game": "game cellular automaton",
    "wca2d": "weighted cellular automaton (WCA2D)",
    "local_inertial": "local inertial (LISFLOOD-FP)",
    "local_inertial_frcap": "local inertial, Froude cap 0.544 (weir)",
    "local_inertial_frcap07": "local inertial, Froude cap 0.7",
    "kinematic_wave": "kinematic wave",
    "diffusion_wave_explicit": "diffusion wave, explicit",
    "diffusion_wave_implicit_cg": "diffusion wave, implicit (CG)",
    "diffusion_wave_implicit_adi": "diffusion wave, implicit (ADI)",
    "swe_lbm": "lattice Boltzmann shallow water",
    "ext_anuga": "ANUGA (open-source shallow water)",
    "ext_landlab_overland": "Landlab overland flow (open source)",
    "ext_landlab_kinwave": "Landlab kinematic wave (open source)",
    "ext_landlab_lindiff": "Landlab linear diffusion (open source)",
    "ext_phiflow_smoke": "PhiFlow smoke (open source)",
    "linear_wave": "linear wave equation",
    "heat_explicit": "heat equation, explicit",
    "heat_adi": "heat equation, implicit (ADI)",
    "stable_fluids": "stable fluids (Stam)",
    "lbm_d2q9": "lattice Boltzmann (D2Q9)",
}

# the order a person meets the scenes in: recognisable water first, other physics last
ORDER = ["dam_break", "river", "lake_levelset", "ocean_edge", "pit_dig", "moat_evap", "rain_groundwater",
         "swell", "lake_at_rest", "ripples", "gaussian", "dye", "smoke_plume", "cylinder_wake"]
TITLE = {
    "dam_break": "Dam break: a breach opens in a dam holding 2 m of water",
    "river": "River: a meandering channel fed upstream, draining downstream",
    "lake_levelset": "Lake draining through an outlet channel",
    "ocean_edge": "Ocean edge: a river runs down a beach into the sea",
    "pit_dig": "Pit dig: a pit and a canal are dug beside a river, the canal floods the pit",
    "moat_evap": "Moat drying out under (exaggerated) evaporation",
    "rain_groundwater": "Rain on a dry basin, then groundwater seeping up",
    "swell": "Swell: the sea rises and falls onto a beach behind two breakwaters",
    "lake_at_rest": "Still lake with an island (nothing should move)",
}
# simulated seconds worth watching, as a fraction of the domain crossing at 128 cells (scaled with n)
WATCH = {"dam_break": 36.0}
SECTION_ROW = {"dam_break": 0.5, "swell": 0.40}   # side view along this fraction of the height


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


def run_one(cls, sc, samples, max_wall):
    """Advance one solver through the common sample points. Returns per-frame records."""
    s = cls.from_scene(copy.deepcopy(sc))
    events = sorted(copy.deepcopy(sc.events), key=lambda e: e[0])
    flow = sc.kind == "flow"
    tot0 = _total(s, sc.kind)
    frames, status, note = [], "ok", ""
    step_wall, t_start = 0.0, time.perf_counter()
    last_z = [None]

    def record():
        snap = s.snapshot()
        r = {"t": s.t, "steps": s.nsteps, "status": status, "note": note,
             "ms": 1000 * step_wall / s.nsteps if s.nsteps else None, "field": _primary(sc.kind, snap)}
        if sc.kind == "water":
            r["vol_err"] = abs(s.volume_error()) / max(s.V0, s.ledger.throughput, 1e-12)
            z = np.asarray(snap["z"], np.float32)
            if last_z[0] is not None and np.array_equal(z, last_z[0]):
                z = last_z[0]                       # share the bed between frames until it is dug
            last_z[0] = z
            r["snap"] = {"z": z, "h": r["field"], "wall": np.asarray(snap["wall"], bool)}
        else:
            r["vol_err"] = abs(_total(s, sc.kind) - tot0) / abs(tot0) if tot0 else None
            r["snap"] = {k: (np.asarray(v, np.float32) if isinstance(v, np.ndarray) and v.dtype != bool else v)
                         for k, v in snap.items()}
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
                    gap = target - s.t
                    dt = min(s.stable_dt(), gap)
                    if sc.kind == "water":
                        dt = min(dt, WATER_DT_CAP)
                    if dt < gap:
                        # equal steps to the frame: a short tail step after full ones (0.19 s then 0.11 s on dam_break)
                        # grew local_inertial's reservoir roughness 500-fold by 33 s (measured 2026-11-14)
                        dt = gap / math.ceil(gap / dt - 1e-9)
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
                status, note = "blown up", f"at t={s.t:.3g} s"
                break
            if time.perf_counter() - t_start > max_wall:
                status, note = "stopped", f"hit the {max_wall:.0f} s wall-time cap at t={s.t:.3g} s"
                break
        if status != "ok" and frames and frames[-1]["status"] == "ok" and status != "blown up":
            frames.append(record())      # one frame showing where it got to
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
    runs = {}
    for name in order:
        a = time.perf_counter()
        runs[name] = run_one(CLASSES[name], sc, samples, max_wall)
        last = runs[name][-1]
        print(f"  {name:28s} t={last['t']:9.3f} steps={last['steps']:7d} {last['status']:9s} "
              f"{time.perf_counter() - a:6.1f}s {last['note']}", flush=True)
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

_FONTS = {}


def _font(px, bold=False):
    key = (px, bold)
    if key not in _FONTS:
        name = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
        try:
            _FONTS[key] = ImageFont.truetype(name, px)
        except OSError:
            try:
                _FONTS[key] = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/" + name, px)
            except OSError:
                _FONTS[key] = ImageFont.load_default(size=px)
    return _FONTS[key]


def _u8(rgb):
    return Image.fromarray((np.clip(np.nan_to_num(rgb, nan=1.0), 0, 1) * 255).astype(np.uint8))


def _fit(d, text, width, px, bold=False, least=9):
    """Largest font (down to ``least``) that fits ``text`` in ``width``; cut with an ellipsis below that."""
    while px > least and d.textlength(text, font=_font(px, bold)) > width:
        px -= 1
    f = _font(px, bold)
    while d.textlength(text, font=f) > width and len(text) > 4:
        text = text[:-2].rstrip() + "\u2026"
    return text, f


def _wrap(d, text, width, f):
    lines, cur = [], ""
    for w in text.split(" "):
        nxt = (cur + " " + w).strip()
        if cur and d.textlength(nxt, font=f) > width:
            lines.append(cur)
            cur = w
        else:
            cur = nxt
    return lines + [cur]


def plain(name):
    return PLAIN.get(name, name.replace("_", " "))


def _tile(img, strip, tile_px, f, name, is_ref, clock, mark_row=None):
    """One solver: plain name, code name + numbers, the picture, a banner if it failed, the side view."""
    head, sh = 44, (strip.height + 4 if strip is not None else 0)
    out = Image.new("RGB", (tile_px, head + tile_px + sh), (20, 20, 22))
    if img.size != (tile_px, tile_px):
        img = img.resize((tile_px, tile_px), Image.NEAREST if img.width < tile_px else Image.LANCZOS)
    out.paste(img, (0, head))
    if strip is not None:
        out.paste(strip, (0, head + tile_px + 4))
    d = ImageDraw.Draw(out)
    t, fn = _fit(d, plain(name) + ("  [reference]" if is_ref else ""), tile_px - 10, 15, True)
    d.text((6, 3), t, fill=(255, 255, 255), font=fn)
    ve = "-" if f["vol_err"] is None else f"{f['vol_err']:.0e}"
    ms = "-" if f["ms"] is None else f"{f['ms']:.2f}"
    t, fn = _fit(d, f"{name}   volume error {ve}   {ms} ms/step", tile_px - 10, 12)
    d.text((6, 25), t, fill=(185, 185, 190), font=fn)
    tw = d.textlength(clock, font=_font(13, True))
    yb = head + tile_px
    d.rectangle([tile_px - tw - 12, yb - 24, tile_px - 4, yb - 4], fill=(0, 0, 0))
    d.text((tile_px - tw - 8, yb - 22), clock, fill=(255, 255, 255), font=_font(13, True))
    if mark_row is not None:                    # where the side view is cut
        y = head + mark_row
        for x0, sgn in ((0, 1), (tile_px - 1, -1)):
            d.polygon([(x0, y - 7), (x0, y + 7), (x0 + sgn * 10, y)], fill=(255, 235, 80))
    if f["status"] != "ok":
        col = STOP if f["status"] == "stopped" else BAD
        msg = {"blown up": f"BLEW UP {f['note']}: frozen at last good frame",
               "stopped": f"STOPPED at t={f['t']:.1f} s: too slow (wall-time cap)",
               "error": f"CRASHED: {f['note']}"}[f["status"]]
        d.rectangle([0, head, tile_px - 1, head + 26], fill=col)
        t, fn = _fit(d, msg, tile_px - 12, 13, True)
        d.text((6, head + 5), t, fill=(255, 255, 255), font=fn)
        d.rectangle([0, 0, tile_px - 1, out.height - 1], outline=col, width=4)
    return out


def _sheet(tiles, cols, title, legend, clock):
    w, h = tiles[0].size
    rows = math.ceil(len(tiles) / cols)
    gap = 6
    W = cols * w + (cols - 1) * gap
    probe = ImageDraw.Draw(Image.new("RGB", (8, 8)))
    cw = probe.textlength(clock, font=_font(26, True)) + 24
    lines = _wrap(probe, legend, W - 16, _font(13))[:3]
    top = 40 + 18 * len(lines) + 6
    out = Image.new("RGB", (W, top + rows * (h + gap)), (10, 10, 12))
    d = ImageDraw.Draw(out)
    t, fn = _fit(d, title, W - cw - 16, 20, True)
    d.text((8, 6), t, fill=(255, 255, 255), font=fn)
    for i, ln in enumerate(lines):
        d.text((8, 38 + 18 * i), ln, fill=(200, 200, 205), font=_font(13))
    d.text((W - cw + 12, 6), clock, fill=(255, 255, 255), font=_font(26, True))
    for k, t in enumerate(tiles):
        out.paste(t, ((k % cols) * (w + gap), top + (k // cols) * (h + gap)))
    return out


def _clock(sc, samples, k):
    return f"step {samples[k]}" if sc.kind == "flow" else f"t = {samples[k]:.1f} s"


class _Video:
    """Frames piped straight into ffmpeg (no temporary images)."""

    def __init__(self, path, size, fps, codec):
        w, h = size
        enc = (["-c:v", "libx264", "-crf", "22", "-preset", "medium", "-g", str(fps)] if codec == "h264" else
               ["-c:v", "libvpx-vp9", "-crf", "33", "-b:v", "0", "-g", str(fps), "-row-mt", "1",
                "-deadline", "good", "-cpu-used", "4"])
        self.p = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
                                   "-s", f"{w}x{h}", "-framerate", str(fps), "-i", "-",
                                   "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-pix_fmt", "yuv420p", *enc, str(path)],
                                  stdin=subprocess.PIPE)
        self.path = path

    def add(self, im):
        self.p.stdin.write(im.tobytes())

    def close(self):
        self.p.stdin.close()
        return self.p.wait() == 0


def render_outputs(sc, order, samples, runs, tile_px, d: Path, stills_dir: Path | None, fps=10):
    """Draw every frame (state and difference), stream them into videos, keep small GIFs and stills."""
    d.mkdir(parents=True, exist_ok=True)
    water = sc.kind == "water"
    ref = order[0]
    nfr = len(samples)
    up = max(1, tile_px // max(sc.grid.nx, sc.grid.ny))
    look = sec = None
    if water:
        look = waterview.Look.from_frames([f["snap"] for f in runs[ref]], sc.grid.dx, up)
        if sc.name in SECTION_ROW:
            sec = waterview.Section.from_frames([f["snap"] for f in runs[ref]],
                                                int(SECTION_ROW[sc.name] * sc.grid.ny))
    else:
        st = _global_style(sc.kind, [f["snap"] for f in runs[ref]])
    diffs = {n: [f["field"] - rf["field"] for f, rf in zip(runs[n], runs[ref])] for n in order[1:]}
    allabs = np.concatenate([np.abs(x[np.isfinite(x)]).ravel() for v in diffs.values() for x in v] or [np.zeros(1)])
    scale = float(np.percentile(allabs, 99.5)) if allabs.size else 1.0
    scale = scale if scale > 0 else 1.0
    unit = {"water": "m of depth", "wave": "height", "scalar": "value", "flow": "dye"}[sc.kind]
    refname = plain(ref)
    title = TITLE.get(sc.name, sc.name.replace("_", " ")) + f"   ({sc.grid.nx} x {sc.grid.ny} cells)"
    if water:
        leg_s = ("Same starting grid for every solver, one clock.  Blue = water, darker = deeper; light and "
                 "dark bands = wave fronts; green/brown = dry land; dark grey = wall."
                 + ("  Strip: side view along the yellow marks." if sec else ""))
    else:
        leg_s = "Same starting grid for every solver, one clock, one colour scale on every tile."
    leg_d = (f"Each solver minus the reference ({refname}).  RED = more {unit.split(' of ')[-1]} than the "
             f"reference here, BLUE = less, grey = the same.  Full colour = {scale:.2g} {unit}.")
    cols = min(len(order), 5 if len(order) > 8 else 4)
    mark = None
    if sec is not None:
        mark = int(((sc.grid.ny - 1 - sec.row) + 0.5) * up * tile_px / (up * sc.grid.ny))
    vids, gif_s, gif_d = {}, [], []
    stills = sorted({round(nfr * x) for x in (0.12, 0.45, 0.9)} | ({nfr - 1} if nfr < 4 else set()))
    stills = [min(k, nfr - 1) for k in stills]
    gif_every = max(1, nfr // 40)
    for k in range(nfr):
        clock = _clock(sc, samples, k)
        ts, td = [], []
        for name in order:
            f = runs[name][k]
            sn = f["snap"]
            if water:
                img = _u8(waterview.water_view(sn["z"], sn["h"], sn["wall"], look, sc.grid.dx))
                strip = _u8(waterview.section_strip(sn["z"], sn["h"], sn["wall"], sec, img.width)) if sec else None
            else:
                img, strip = _u8(render.frame(sn, st, sc.grid.dx)), None
            ts.append(_tile(img, strip, tile_px, f, name, name == ref, clock, mark))
            if name == ref:
                td.append(ts[-1])
                continue
            if water:
                dimg = _u8(waterview.diff_view(diffs[name][k], scale, sn["z"], sn["wall"], look, sc.grid.dx))
            else:
                from matplotlib import colormaps
                dimg = _u8(colormaps["RdBu_r"](np.clip(0.5 + 0.5 * np.nan_to_num(diffs[name][k] / scale), 0, 1))
                           [..., :3][::-1])
            td.append(_tile(dimg, strip, tile_px, f, name, False, clock, mark))
        a = _sheet(ts, cols, title, leg_s, clock)
        b = _sheet(td, cols, title + "  - difference from reference", leg_d, clock)
        if k == 0:
            for key, ext, codec in (("state_mp4", "state.mp4", "h264"), ("state_webm", "state.webm", "vp9"),
                                    ("diff_webm", "diff.webm", "vp9")):
                vids[key] = _Video(d / ext, a.size, fps, codec) if shutil.which("ffmpeg") else None
        for key, im in (("state_mp4", a), ("state_webm", a), ("diff_webm", b)):
            if vids.get(key):
                vids[key].add(im)
        if k % gif_every == 0 or k == nfr - 1:
            sm = (640, round(a.height * 640 / a.width))
            gif_s.append(a.resize(sm, Image.LANCZOS).convert("P", palette=Image.ADAPTIVE, colors=128))
            gif_d.append(b.resize(sm, Image.LANCZOS).convert("P", palette=Image.ADAPTIVE, colors=128))
        if k in stills:
            for dd in [d] + ([stills_dir] if stills_dir else []):
                dd.mkdir(parents=True, exist_ok=True)
                a.save(dd / f"{d.name}_state_t{samples[k]:06.1f}.jpg".replace("/", ""), quality=90)
            b.save(d / f"{d.name}_diff_t{samples[k]:06.1f}.jpg", quality=88)
        if k == min(nfr - 1, round(nfr * 0.45)):
            a.resize((960, round(a.height * 960 / a.width)), Image.LANCZOS).save(d / "poster.webp", quality=70)
    ok = {key: v.close() for key, v in vids.items() if v}
    for nm, ims in (("state", gif_s), ("diff", gif_d)):
        ims[0].save(d / f"{nm}.gif", save_all=True, append_images=ims[1:] + [ims[-1]] * 6,
                    duration=int(1000 * gif_every / fps), loop=0, optimize=True)
    return scale, fps, {key: v for key, v in ok.items()}


# ---------------------------------------------------------------------------------- outputs

def rms_vs_ref(ref_frames, frames):
    """RMS of (solver - reference) primary field over the frames both reached, in the field's units
    (water: depth in m). None when the solver never got past the start."""
    k = min(len(ref_frames), len(frames))
    sq = [np.mean((np.asarray(frames[i]["field"], float) - np.asarray(ref_frames[i]["field"], float)) ** 2)
          for i in range(1, k) if frames[i]["status"] == "ok" or i < k - 1]
    sq = [v for v in sq if np.isfinite(v)]
    return float(np.sqrt(np.mean(sq))) if sq else None


def save_meta(sc, order, samples, runs, scale, fps, videos, d: Path):
    final = {n: {k: runs[n][-1][k] for k in ("t", "steps", "status", "note", "ms", "vol_err")} for n in order}
    for n in order:
        final[n]["rms_vs_ref"] = rms_vs_ref(runs[order[0]], runs[n])
        final[n]["plain"] = plain(n)
    meta = {"scene": sc.name, "title": TITLE.get(sc.name, sc.name), "kind": sc.kind,
            "n": [sc.grid.ny, sc.grid.nx], "dx": sc.grid.dx, "duration": sc.duration, "note": sc.note,
            "ref": order[0], "solvers": order, "samples": samples, "diff_scale": scale,
            "frames": len(samples), "fps": fps, "videos": videos, "final": final,
            "made": time.strftime("%Y-%m-%dT%H:%M:%S")}
    (d / "meta.json").write_text(json.dumps(meta, indent=1, default=float))


def _rank(meta):
    s = meta["scene"]
    return (ORDER.index(s) if s in ORDER else len(ORDER), meta["n"][1])


def _b64(p: Path) -> str:
    return base64.b64encode(p.read_bytes()).decode()


def build_html(root: Path) -> Path:
    """index.html (water scenes first, other physics last) plus one self-contained page per saved run."""
    metas = []
    for m in root.glob("*/meta.json"):
        meta = json.loads(m.read_text())
        if (m.parent / "state.webm").exists() and (m.parent / "diff.webm").exists():
            meta["id"] = m.parent.name
            metas.append(meta)
    metas.sort(key=_rank)
    for meta in metas:
        d = root / meta["id"]
        data = {k: meta[k] for k in ("id", "title", "kind", "n", "note", "ref", "solvers", "samples",
                                     "diff_scale", "fps", "final")}
        data["state"], data["diff"] = _b64(d / "state.webm"), _b64(d / "diff.webm")
        (root / f"{meta['id']}.html").write_text(PAGE.replace("__DATA__", json.dumps(data, default=float)))
    cards = {"water": [], "other": []}
    for meta in metas:
        poster = d = root / meta["id"] / "poster.webp"
        img = f'<img src="data:image/webp;base64,{_b64(poster)}">' if poster.exists() else ""
        bad = [v["plain"] for v in meta["final"].values() if v["status"] != "ok"]
        span = (f"{meta['samples'][-1]:.0f} s simulated" if meta["kind"] != "flow"
                else f"{meta['samples'][-1]} steps")
        card = (f'<a class="card" href="{meta["id"]}.html">{img}<b>{meta["title"]}</b>'
                f'<span>{meta["n"][1]} x {meta["n"][0]} cells, {len(meta["solvers"])} solvers, {span}, '
                f'{meta["frames"]} frames.{" Failed or stopped: " + ", ".join(bad) + "." if bad else ""}</span></a>')
        cards["water" if meta["kind"] == "water" else "other"].append(card)
    page = INDEX.replace("__WATER__", "\n".join(cards["water"]) or "<p>no water runs saved</p>")
    page = page.replace("__OTHER__", "\n".join(cards["other"]) or "<p>none saved</p>")
    out = root / "index.html"
    out.write_text(page)
    return out


INDEX = """<!doctype html><html><head><meta charset="utf-8"><title>Same water, every solver</title>
<style>body{background:#111;color:#ddd;font:15px system-ui,sans-serif;margin:16px 24px;max-width:1500px}
h1{font-size:22px;margin:4px 0}h2{font-size:18px;margin:22px 0 8px;color:#fff}p{max-width:1000px;line-height:1.4}
.grid{display:flex;flex-wrap:wrap;gap:14px}.card{display:flex;flex-direction:column;width:460px;
background:#1c1c20;border:1px solid #333;border-radius:6px;padding:8px;color:#ddd;text-decoration:none}
.card:hover{border-color:#6af}.card img{width:100%;border-radius:3px;margin-bottom:6px}.card b{color:#fff}
.card span{font-size:13px;color:#aaa;margin-top:4px}.other .card{width:300px}</style></head><body>
<h1>Same water, every solver</h1>
<p>Each page starts every water solver from the <b>identical grid</b> (same ground, same water, same walls),
advances them all on <b>one clock</b>, and shows them side by side, so any difference you see is the solver.
Open a scene and press play. The picture is drawn the same way in every tile: ground as shaded relief,
water in blue (darker = deeper), wave fronts as light and dark bands, walls dark grey.</p>
<h2>Water</h2><div class="grid">__WATER__</div>
<h2>Other physics (not water: heat, waves on a membrane, smoke, dye, flow past a cylinder)</h2>
<div class="grid other">__OTHER__</div></body></html>"""


PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>compare</title>
<style>body{background:#111;color:#ddd;font:15px system-ui,sans-serif;margin:10px 16px}
a{color:#8bf}button,select{font-size:15px}#bar{display:flex;align-items:center;gap:12px;flex-wrap:wrap;
margin:6px 0}input[type=range]{width:520px}video{max-width:100%;display:block;background:#000}
#lab{font:bold 18px system-ui;min-width:110px}#legend{margin:4px 0 6px;color:#ccc}
table{border-collapse:collapse;margin-top:10px;font-size:13px}td,th{border:1px solid #333;padding:3px 7px}
.bad{color:#f66}.stop{color:#fa4}</style></head><body>
<div><a href="index.html">&larr; all scenes</a> &nbsp; <b id="title"></b></div>
<div id="note"></div>
<div id="bar"><button id="play">play</button><button id="back">&lsaquo; frame</button><button id="fwd">frame &rsaquo;</button>
<input type="range" id="t" min="0" value="0"><span id="lab"></span>
<label>speed <select id="rate"><option>0.25</option><option>0.5</option><option selected>1</option><option>2</option></select></label>
<label><input type="radio" name="v" value="state" checked> what the water does</label>
<label><input type="radio" name="v" value="diff"> difference from reference</label></div>
<div id="legend"></div>
<video id="vs" muted loop playsinline preload="auto"></video><video id="vd" muted loop playsinline preload="auto" style="display:none"></video>
<table id="tab"></table>
<script>
const D=__DATA__;const $=id=>document.getElementById(id);
function blob(b){const s=atob(b),u=new Uint8Array(s.length);for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i);
return URL.createObjectURL(new Blob([u],{type:'video/webm'}))}
$('vs').src=blob(D.state);$('vd').src=blob(D.diff);let v=$('vs'),view='state';const N=D.samples.length;
document.title=D.title;$('title').textContent=D.title+'  ('+D.n[1]+' x '+D.n[0]+' cells)';$('note').textContent='Set-up: '+D.note;
$('t').max=N-1;const unit=D.kind=='water'?' m of depth':'';
const LEG={state:'Every tile started from the same grid; all share one clock and one colour scale. Tiles with a red or orange banner blew up or were stopped and stay frozen at their last good frame.',
diff:'Each tile is that solver minus the reference ('+D.final[D.ref].plain+'). RED = more water than the reference at that spot, BLUE = less, grey = the same. Full colour = '+(+D.diff_scale).toPrecision(2)+unit+'. The reference tile itself is shown as it is.'};
function frame(){return Math.min(N-1,Math.max(0,Math.floor(v.currentTime*D.fps+1e-3)))}
function lab(k){$('lab').textContent=D.kind=='flow'?'step '+D.samples[k]:'t = '+(+D.samples[k]).toFixed(1)+' s'}
function seek(k){k=Math.min(N-1,Math.max(0,k));v.currentTime=(k+0.5)/D.fps;$('t').value=k;lab(k)}
function upd(){const k=frame();$('t').value=k;lab(k);if(!v.paused)requestAnimationFrame(upd)}
$('play').onclick=()=>{if(v.paused){v.play();$('play').textContent='pause';requestAnimationFrame(upd)}else{v.pause();$('play').textContent='play'}};
$('t').oninput=e=>{v.pause();$('play').textContent='play';seek(+e.target.value)};
$('back').onclick=()=>{v.pause();$('play').textContent='play';seek(frame()-1)};
$('fwd').onclick=()=>{v.pause();$('play').textContent='play';seek(frame()+1)};
$('rate').onchange=e=>{$('vs').playbackRate=$('vd').playbackRate=+e.target.value};
document.querySelectorAll('input[name=v]').forEach(x=>x.onchange=e=>{const k=frame(),was=!v.paused;v.pause();
view=e.target.value;$('vs').style.display=view=='state'?'':'none';$('vd').style.display=view=='diff'?'':'none';
v=view=='state'?$('vs'):$('vd');$('legend').textContent=LEG[view];seek(k);if(was){v.play();requestAnimationFrame(upd)}});
$('legend').textContent=LEG.state;
let h='<tr><th>solver</th><th>code name</th><th>status</th><th>reached t</th><th>steps</th><th>ms/step</th><th>volume error</th><th>rms vs reference</th></tr>';
for(const n of D.solvers){const f=D.final[n];const c=f.status=='ok'?'':(f.status=='stopped'?'stop':'bad');
h+=`<tr class="${c}"><td>${f.plain}${n==D.ref?' (reference)':''}</td><td>${n}</td><td>${f.status} ${f.note||''}</td><td>${(+f.t).toFixed(1)}</td><td>${f.steps}</td>`+
`<td>${f.ms==null?'-':(+f.ms).toFixed(3)}</td><td>${f.vol_err==null?'-':(+f.vol_err).toExponential(1)}</td>`+
`<td>${f.rms_vs_ref==null?'-':(+f.rms_vs_ref).toPrecision(3)}</td></tr>`}
$('tab').innerHTML=h;
$('vs').addEventListener('loadeddata',()=>seek(Math.round(N*0.12)),{once:true});
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
    p.add_argument("--frames", type=int, default=120)
    p.add_argument("--duration", type=float, default=None, help="override the scene's simulated seconds")
    p.add_argument("--max-wall", type=float, default=45.0, help="wall seconds per solver per scene")
    p.add_argument("--tile", type=int, default=0, help="tile size in pixels (default: >=320, a whole multiple of n)")
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
            names = {"all": [s for s in ORDER if s in scenes.ALL],
                     "water": [s for s in ORDER if s in scenes.WATER]}.get(a.scene, a.scene.split(","))
            built = []
            for nm in names:
                sc = scenes.build(nm, a.n, a.seed)
                if a.duration:
                    sc.duration = a.duration
                elif nm in WATCH:                   # the event, not the settled tail
                    sc.duration = WATCH[nm] * max(1.0, a.n / 128)
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
            n = max(sc.grid.nx, sc.grid.ny)
            tile = a.tile or n * math.ceil(320 / n)
            d = root / f"{sc.name}_{sc.grid.nx}"
            d.mkdir(parents=True, exist_ok=True)
            if sc.kind == "water":
                np.savez_compressed(d / "init.npz", z=sc.init["z"], h=sc.init["h"], wall=sc.init["wall"])
            for old in d.glob("*.jpg"):
                old.unlink()
            a0 = time.perf_counter()
            scale, fps, vids = render_outputs(sc, order, samples, runs, tile, d, root / "stills")
            save_meta(sc, order, samples, runs, scale, fps, vids, d)
            del runs
            print(f"  -> {d} (drawn in {time.perf_counter() - a0:.0f}s; videos {vids})", flush=True)
    html = build_html(root)
    pages = sum(p.stat().st_size for p in root.glob("*.html")) / 1e6
    print(f"html: {html} (+ one page per scene, {pages:.1f} MB in all)")


if __name__ == "__main__":
    main()
