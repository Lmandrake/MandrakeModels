"""Draw a water state the way a person expects to see water: for the side-by-side viewer.

Terrain is shaded relief of the bed (earth tones, lit from the north-west). Water lies over it in blues:
deeper is darker and more opaque, thin water lets the bed show through. The water surface is lit from
the same side using the slope of the surface height ``eta = z + h``, so wave fronts, bores and ripples
show as light and dark bands, and very steep fronts get a white crest. A pale line marks the shoreline.
Walls are dark stone. Dry land stays land.

Every tile of one scene uses one ``Look`` (terrain range, depth scale, surface-shading gain), taken
from the reference solver's whole run, so a colour means the same depth in every tile.
The grid is upsampled bilinearly to the display size, the same light smoothing for every solver.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage

WET = 0.01                       # m; thinner water is drawn as dry
EARTH = np.array([[0.47, 0.42, 0.33],   # low ground: damp mud
                  [0.74, 0.67, 0.48],   # sand
                  [0.52, 0.58, 0.33],   # grass
                  [0.44, 0.48, 0.30],   # darker grass
                  [0.58, 0.54, 0.48]])  # high: stony
SHALLOW = np.array([0.36, 0.66, 0.80])
DEEP = np.array([0.03, 0.14, 0.36])
CREST = np.array([0.93, 0.96, 0.98])
SHORE = np.array([0.90, 0.88, 0.78])
STONE = np.array([0.24, 0.23, 0.22])


@dataclass
class Look:
    zlo: float
    zhi: float
    depth_ref: float        # m; water this deep is ~63 % of the way to the deep colour
    slope_k: float          # surface shading = tanh(slope_k * slope toward the light)
    up: int                 # display pixels per cell
    relief: float = 3.0     # vertical exaggeration of the bed's hillshade

    @classmethod
    def from_frames(cls, frames, dx, up):
        """``frames``: list of dicts with z, h, wall from the reference run."""
        z0, wall = frames[0]["z"], frames[0]["wall"]
        zs = z0[~wall] if (~wall).any() else z0.ravel()
        zlo, zhi = float(np.percentile(zs, 1)), float(np.percentile(zs, 99))
        if zhi - zlo < 0.5:
            zlo, zhi = zlo - 0.5, zhi + 0.5
        wet = np.concatenate([f["h"][f["h"] > WET].ravel() for f in frames] or [np.ones(1)])
        dref = float(np.percentile(wet, 60)) if wet.size else 1.0
        sl = np.concatenate([np.abs(_surface_slope(f["z"], f["h"], f["wall"], dx)[2]).ravel()
                             for f in frames[1:]] or [np.zeros(1)])
        sl = sl[sl > 0]
        p = float(np.percentile(sl, 98)) if sl.size else 0.0
        k = 1.0 / max(p, 2e-3)          # flat lakes: no amplified rounding noise
        return cls(zlo, zhi, max(dref, 0.05), k, up)


def _surface_slope(z, h, wall, dx):
    """Slope of eta toward the light (north-west), only across faces where both cells are wet."""
    eta = z + h
    wet = (h > WET) & ~wall
    gx = np.zeros_like(eta)
    gy = np.zeros_like(eta)
    cx = np.zeros_like(eta)
    cy = np.zeros_like(eta)
    d = np.where(wet[:, 1:] & wet[:, :-1], eta[:, 1:] - eta[:, :-1], 0.0)
    m = (wet[:, 1:] & wet[:, :-1]).astype(float)
    gx[:, 1:] += d; gx[:, :-1] += d; cx[:, 1:] += m; cx[:, :-1] += m
    d = np.where(wet[1:, :] & wet[:-1, :], eta[1:, :] - eta[:-1, :], 0.0)
    m = (wet[1:, :] & wet[:-1, :]).astype(float)
    gy[1:, :] += d; gy[:-1, :] += d; cy[1:, :] += m; cy[:-1, :] += m
    gx = gx / np.maximum(cx, 1) / dx
    gy = gy / np.maximum(cy, 1) / dx
    # array row 0 is south; light from the north-west: lit where the surface faces -x and +y
    toward = (-gx + gy) / np.sqrt(2.0)
    # wet-weighted smoothing over ~1 cell: removes the checker a staircase bank leaves in the gradient
    w = wet.astype(float)
    toward = ndimage.gaussian_filter(toward * w, 0.8) / np.maximum(ndimage.gaussian_filter(w, 0.8), 1e-6) * w
    return gx, gy, toward


def _up(a, f):
    if f == 1:
        return a
    return ndimage.zoom(a, f, order=1, mode="nearest", grid_mode=True)


def _hillshade(z, dx):
    gy, gx = np.gradient(z, dx)
    nx, ny, nz = -gx, -gy, np.ones_like(z)
    n = np.sqrt(nx ** 2 + ny ** 2 + 1)
    lx, ly, lz = -0.5, 0.5, 0.707          # north-west, 45 degrees up
    return np.clip((nx * lx + ny * ly + nz * lz) / n, 0, 1)


def terrain_rgb(z, wall, look, dx):
    t = np.clip((z - look.zlo) / (look.zhi - look.zlo), 0, 1) * (len(EARTH) - 1)
    i = np.minimum(t.astype(int), len(EARTH) - 2)
    f = (t - i)[..., None]
    col = EARTH[i] * (1 - f) + EARTH[i + 1] * f
    hs = _hillshade(np.where(wall, z, z) * look.relief, dx)
    hs0 = _hillshade(np.zeros((3, 3)), 1.0)[1, 1]
    shade = 0.80 + 0.55 * (hs - hs0)       # flat ground keeps its colour
    return np.clip(col * shade[..., None], 0, 1)


def water_view(z, h, wall, look, dx):
    """RGB image (display size, north up) of one water state."""
    z = np.asarray(z, float)
    h = np.where(wall, 0.0, np.maximum(np.asarray(h, float), 0.0))
    _, _, toward = _surface_slope(z, h, wall, dx)
    up = look.up
    zu, hu, tu = _up(z, up), np.maximum(_up(h, up), 0.0), _up(toward, up)
    wu = np.kron(wall, np.ones((up, up), bool)) if up > 1 else wall
    img = terrain_rgb(zu, wu, look, dx / up)
    wet = (hu > WET) & ~wu
    f = (1 - np.exp(-hu / look.depth_ref))[..., None]
    col = SHALLOW * (1 - f) + DEEP * f
    s = np.tanh(look.slope_k * tu)[..., None]
    col = np.where(s > 0, col + (1 - col) * 0.55 * s, col * (1 + 0.55 * s))
    crest = np.clip(np.abs(look.slope_k * tu) - 1.6, 0, 1)[..., None] * 0.7
    col = col * (1 - crest) + CREST * crest
    a = (np.clip(hu / 0.12, 0, 1) * (0.62 + 0.36 * f[..., 0]))[..., None]
    img = np.where(wet[..., None], img * (1 - a) + col * a, img)
    edge = wet & ~ndimage.binary_erosion(wet, iterations=1, border_value=1)
    img[edge] = img[edge] * 0.3 + SHORE * 0.7
    img[wu] = STONE
    return np.clip(img[::-1], 0, 1)


def diff_view(d, scale, z, wall, look, dx):
    """Solver depth minus reference depth over grey relief: red = more water, blue = less."""
    from matplotlib import colormaps
    up = look.up
    zu = _up(np.asarray(z, float), up)
    wu = np.kron(wall, np.ones((up, up), bool)) if up > 1 else wall
    grey = terrain_rgb(zu, wu, look, dx / up).mean(-1, keepdims=True) * np.ones(3)
    grey = 0.55 + 0.4 * grey
    du = _up(np.nan_to_num(np.asarray(d, float)), up)
    c = colormaps["RdBu_r"](np.clip(0.5 + 0.5 * du / scale, 0, 1))[..., :3]
    a = np.clip(np.abs(du) / (0.1 * scale), 0, 1)[..., None]
    img = grey * (1 - a) + c * a
    img[wu] = STONE
    return np.clip(img[::-1], 0, 1)


@dataclass
class Section:
    row: int
    lo: float
    hi: float

    @classmethod
    def from_frames(cls, frames, row):
        z = frames[0]["z"][row]
        wall = frames[0]["wall"][row]
        etas = [f["z"][row] + f["h"][row] for f in frames]
        wet = [e[(f["h"][row] > WET) & ~wall] for e, f in zip(etas, frames)]
        wet = np.concatenate([w for w in wet if w.size] or [z])
        elo, ehi = float(wet.min()), float(wet.max())
        er = max(ehi - elo, 0.2)
        zlo = float(z[~wall].min()) if (~wall).any() else float(z.min())
        lo = max(zlo - 0.1 * er, elo - 2 * er - 0.3)
        return cls(row, lo, ehi + 0.6 * er + 0.2)


def section_strip(z, h, wall, sec, width, height=96):
    """Side view along one grid row: bed brown, water blue, walls dark; same vertical range per scene."""
    zr = np.asarray(z[sec.row], float)
    hr = np.where(wall[sec.row], 0.0, np.maximum(np.asarray(h[sec.row], float), 0.0))
    wr = wall[sec.row]
    n = zr.size
    xs = (np.arange(width) + 0.5) / width * n - 0.5
    zc = np.interp(xs, np.arange(n), zr)
    ec = np.interp(xs, np.arange(n), zr + hr)
    wc = wr[np.clip(np.round(xs).astype(int), 0, n - 1)]
    hc = ec - zc
    ys = sec.hi - (np.arange(height) + 0.5) / height * (sec.hi - sec.lo)     # row 0 at the top
    Y = ys[:, None]
    img = np.empty((height, width, 3))
    img[:] = (0.86, 0.90, 0.95)                                             # sky
    water = (Y <= ec[None, :]) & (Y > zc[None, :]) & (hc[None, :] > WET)
    img[water] = (0.16, 0.42, 0.72)
    img[Y <= zc[None, :]] = (0.55, 0.43, 0.30)
    img[:, wc] = STONE
    line = np.abs(Y - ec[None, :]) < (sec.hi - sec.lo) / height * 0.8
    img[line & (hc[None, :] > WET) & ~wc[None, :]] = (0.05, 0.20, 0.45)
    return np.clip(img, 0, 1)
