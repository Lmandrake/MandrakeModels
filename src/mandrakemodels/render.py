"""Turn solver snapshots into RGB frames, PNG stills and GIFs (numpy colouring, PIL output).

Water: hillshaded bed (sand to grass by height), water blended over it by depth (shallow turquoise to
deep navy), white foam where it runs fast, walls dark. Wave: shaded heightfield in a blue-white
diverging ramp. Scalar and flow: a matplotlib colormap on a fixed range taken from the first frame so a
run's frames compare. Pass the same ``Style`` to every frame of one run.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from matplotlib import colormaps
from PIL import Image


def _lerp(a, b, t):
    t = np.clip(t, 0, 1)[..., None]
    return a * (1 - t) + b * t


def _hillshade(z, dx=1.0, az=315.0, alt=45.0):
    gy, gx = np.gradient(z, dx)
    slope = np.arctan(np.hypot(gx, gy))
    aspect = np.arctan2(-gx, gy)
    a, al = np.radians(az), np.radians(alt)
    s = np.sin(al) * np.cos(slope) + np.cos(al) * np.sin(slope) * np.cos(a - aspect)
    return np.clip(s, 0, 1)


@dataclass
class Style:
    kind: str
    lo: float = 0.0
    hi: float = 1.0
    zlo: float = 0.0
    zhi: float = 1.0
    depth_ref: float = 2.0
    speed_ref: float = 2.0
    cmap: str = "magma"
    extra: dict = field(default_factory=dict)

    @classmethod
    def from_snapshot(cls, kind: str, snap: dict) -> "Style":
        if kind == "water":
            z = snap["z"]
            open_ = ~snap["wall"]
            zz = z[open_] if open_.any() else z
            hmax = float(np.percentile(snap["h"][snap["h"] > 1e-3], 95)) if (snap["h"] > 1e-3).any() else 1.0
            return cls(kind, zlo=float(zz.min()), zhi=float(zz.max()) + 1e-9, depth_ref=max(hmax, 0.5))
        if kind == "wave":
            return cls(kind, hi=1.0)
        if kind == "scalar":
            f = snap["field"]
            return cls(kind, lo=float(f.min()), hi=float(f.max()) + 1e-12, cmap="inferno")
        if kind == "flow":
            return cls(kind, lo=0.0, hi=max(float(snap["dye"].max()), 1e-9), cmap="magma")
        raise ValueError(kind)


SAND = np.array([0.80, 0.72, 0.52])
GRASS = np.array([0.36, 0.50, 0.26])
ROCK = np.array([0.55, 0.53, 0.50])
SHALLOW = np.array([0.30, 0.75, 0.78])
DEEP = np.array([0.04, 0.12, 0.35])
FOAM = np.array([0.95, 0.97, 1.0])
WALL = np.array([0.20, 0.19, 0.18])


def water_frame(snap: dict, st: Style, dx: float = 1.414) -> np.ndarray:
    z, h, wall, sp = snap["z"], snap["h"], snap["wall"], snap.get("speed")
    zt = (z - st.zlo) / (st.zhi - st.zlo)
    land = np.where(zt[..., None] < 0.55, _lerp(SAND, GRASS, zt / 0.55), _lerp(GRASS, ROCK, (zt - 0.55) / 0.45))
    shade = 0.55 + 0.6 * _hillshade(z, dx)
    land = land * shade[..., None]
    d = np.clip(h / st.depth_ref, 0, 1)
    water = _lerp(SHALLOW, DEEP, d ** 0.7)
    a = np.clip(h / 0.08, 0, 1) * 0.92          # thin films are see-through
    img = land * (1 - a[..., None]) + water * a[..., None]
    if sp is not None:
        f = np.clip((sp - 0.3) / st.speed_ref, 0, 1) * (h > 1e-3)
        img = img * (1 - 0.8 * f[..., None]) + FOAM * (0.8 * f[..., None])
    img[wall] = WALL
    return np.clip(img, 0, 1)


def wave_frame(snap: dict, st: Style) -> np.ndarray:
    u = snap["height"]
    m = float(np.abs(u).max()) or 1.0
    st.extra.setdefault("amp", m)
    amp = st.extra["amp"]
    shade = 0.6 + 0.5 * _hillshade(u / amp * 6.0)
    t = np.clip(0.5 + 0.5 * u / amp, 0, 1)
    img = _lerp(DEEP, np.array([0.85, 0.93, 1.0]), t) * shade[..., None]
    if "obstacle" in snap and snap["obstacle"] is not None:
        img[snap["obstacle"]] = WALL
    return np.clip(img, 0, 1)


def scalar_frame(snap: dict, st: Style) -> np.ndarray:
    f = (snap["field"] - st.lo) / (st.hi - st.lo)
    return colormaps[st.cmap](np.clip(f, 0, 1))[..., :3]


def flow_frame(snap: dict, st: Style) -> np.ndarray:
    d = np.clip(snap["dye"] / st.hi, 0, 1)
    img = colormaps[st.cmap](d ** 0.6)[..., :3]
    if snap.get("vorticity") is not None:
        w = snap["vorticity"]
        wr = st.extra.setdefault("wref", float(np.percentile(np.abs(w), 99)) or 1.0)
        tint = np.clip(w / (wr + 1e-12), -1, 1)
        img = img * 0.85 + 0.15 * colormaps["coolwarm"](0.5 + 0.5 * tint)[..., :3]
    if snap.get("obstacle") is not None:
        img[snap["obstacle"]] = np.array([0.6, 0.6, 0.6])
    return np.clip(img, 0, 1)


def frame(snap: dict, st: Style, dx: float = 1.414) -> np.ndarray:
    """RGB float image, row 0 at the top = largest y (north up)."""
    if st.kind == "water":
        img = water_frame(snap, st, dx)
    elif st.kind == "wave":
        img = wave_frame(snap, st)
    elif st.kind == "scalar":
        img = scalar_frame(snap, st)
    else:
        img = flow_frame(snap, st)
    return img[::-1]


def to_image(img: np.ndarray, max_px: int = 320) -> Image.Image:
    im = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8))
    s = max_px / max(im.size)
    if s < 1:
        im = im.resize((max(1, round(im.width * s)), max(1, round(im.height * s))), Image.LANCZOS)
    elif s >= 2:
        k = int(s)
        im = im.resize((im.width * k, im.height * k), Image.NEAREST)
    return im


def save_png(img: np.ndarray, path, max_px: int = 320) -> None:
    to_image(img, max_px).save(path)


def save_gif(frames: list, path, fps: int = 12, max_px: int = 320) -> None:
    ims = [to_image(f, max_px).convert("P", palette=Image.ADAPTIVE, colors=128) for f in frames]
    if ims:
        ims[0].save(path, save_all=True, append_images=ims[1:], duration=int(1000 / fps), loop=0,
                    optimize=True)
