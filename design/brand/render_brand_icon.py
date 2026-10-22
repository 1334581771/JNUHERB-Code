#!/usr/bin/env python
"""Rasterise the user-designed JNUHERB brand SVG (ui-redesign-proposal/jnuherb-brand-icon.svg).

Dedicated renderer for this file's feature set: path data (M/L/C/Z, implicit
negative separators), <use> instances with translate/rotate/scale transforms,
userSpaceOnUse linear gradients evaluated through inverse affine maps.
"""
import re
import numpy as np
from pathlib import Path
from PIL import Image, ImageDraw

SRC = Path("/mnt/d/herbtar/jnuherb/ui-redesign-proposal/jnuherb-brand-icon.svg")
OUT = Path("/mnt/d/herbtar/jnuherb/src/jnuherb/data/assets/jnuherb-brand-icon.png")

S = 4            # supersample: px per viewBox unit
SIZE = 256       # viewBox is 256x256
SZ = SIZE * S

# ---------------------------------------------------------------- path data
def points_from_path(d):
    tokens = re.findall(r"([MmLlCcZz])|(-?\d*\.?\d+(?:e-?\d+)?)", d)
    cur, start, pts = None, None, []
    i = 0
    nums = [t[1] for t in tokens if t[1] != ""]
    cmds = [t[0] for t in tokens if t[0]]
    j = 0
    out = []
    def isnum(s):
        return re.fullmatch(r"-?\d*\.?\d+(?:e-?\d+)?", s) is not None
    # simpler: sequential scan
    seq = [t for t in tokens]
    k = 0
    n = len(seq)
    last_cmd = None
    cur = (0.0, 0.0)
    start = None
    while k < n:
        cmd, num = seq[k]
        if cmd:
            c = cmd
            k += 1
        else:
            c = last_cmd if last_cmd and last_cmd.upper() in "MLC" and last_cmd.upper() != "Z" else None
            if c is None:
                raise ValueError(f"number without command at token {k}: {num}")
        def take(cnt):
            nonlocal k
            vals = []
            for _ in range(cnt):
                cc, vv = seq[k]
                assert not cc, f"expected number got {cc}"
                vals.append(float(vv)); k += 1
            return vals
        if c in "Mm":
            x, y = take(2)
            if c == "m":
                x += cur[0]; y += cur[1]
            cur = (x, y); start = cur
            out.append(cur)
            last_cmd = "L" if c == "M" else "l"
            if c == "m": last_cmd = "l"
        elif c in "Ll":
            x, y = take(2)
            if c == "l":
                x += cur[0]; y += cur[1]
            cur = (x, y)
            out.append(cur)
            last_cmd = c
        elif c in "Cc":
            x1, y1, x2, y2, x, y = take(6)
            if c == "c":
                x1 += cur[0]; y1 += cur[1]; x2 += cur[0]; y2 += cur[1]
                x += cur[0]; y += cur[1]
            p0 = cur
            for t in np.linspace(0, 1, 24)[1:]:
                mt = 1 - t
                px = mt**3*p0[0] + 3*mt*mt*t*x1 + 3*mt*t*t*x2 + t**3*x
                py = mt**3*p0[1] + 3*mt*mt*t*y1 + 3*mt*t*t*y2 + t**3*y
                out.append((px, py))
            cur = (x, y)
            last_cmd = c
        elif c in "Zz":
            if start:
                out.append(start)
            cur = start
            last_cmd = None
        else:
            raise ValueError(f"unsupported command {c}")
    return out

# ---------------------------------------------------------------- gradients
def hex2rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

BLADE_STOPS = [(0.0, "#355f36"), (0.36, "#608843"), (0.53, "#a3be6b"),
               (0.70, "#7c9f4d"), (1.0, "#b1c77c")]
STEM_STOPS = [(0.0, "#577643"), (0.5, "#789a45"), (1.0, "#a2bb61")]
B_P1 = np.array([-13.0, -6.0]); B_P2 = np.array([17.0, -38.0])
S_P1 = np.array([92.0, 220.0]); S_P2 = np.array([161.0, 43.0])

def grad_color(stops, p1, p2, pt):
    d = p2 - p1
    t = float(np.dot(np.asarray(pt, float) - p1, d) / np.dot(d, d))
    t = min(1.0, max(0.0, t))
    for (t0, c0), (t1, c1) in zip(stops, stops[1:]):
        if t0 <= t <= t1:
            a, b = hex2rgb(c0), hex2rgb(c1)
            f = 0 if t1 == t0 else (t - t0) / (t1 - t0)
            return tuple(int(a[i] + (b[i] - a[i]) * f) for i in range(3))
    return hex2rgb(stops[-1][1])

# ---------------------------------------------------------------- affine util
def parse_use_transform(s):
    """translate(tx ty) rotate(a) scale(sx [sy]) -> 3x3 forward matrix."""
    M = np.eye(3)
    for fn, args in re.findall(r"(translate|rotate|scale)\(([^)]*)\)", s):
        vals = [float(v) for v in re.findall(r"-?\d*\.?\d+(?:e-?\d+)?", args)]
        if fn == "translate":
            T = np.array([[1, 0, vals[0]], [0, 1, vals[1]], [0, 0, 1.0]])
        elif fn == "rotate":
            a = np.radians(vals[0]); c, sn = np.cos(a), np.sin(a)
            T = np.array([[c, -sn, 0], [sn, c, 0], [0, 0, 1.0]])
        else:
            sx = vals[0]; sy = vals[1] if len(vals) > 1 else vals[0]
            T = np.array([[sx, 0, 0], [0, sy, 0], [0, 0, 1.0]])
        M = M @ T                      # SVG applies left-to-right
    return M

def apply(M, pts):
    P = np.asarray(pts, float).T
    H = np.vstack([P, np.ones((1, len(pts)))])
    R = M @ H
    return [(R[0, i], R[1, i]) for i in range(len(pts))]

# ---------------------------------------------------------------- render
txt = SRC.read_text("utf-8")

canvas = Image.new("RGBA", (SZ, SZ), (0, 0, 0, 0))
d = ImageDraw.Draw(canvas)

def S_(v):
    return v * S

def draw_polyline(pts, color, width, opacity=1.0):
    p = [(S_(x), S_(y)) for x, y in pts]
    if len(p) < 2:
        return
    ov = Image.new("RGBA", (SZ, SZ), (0, 0, 0, 0))
    ImageDraw.Draw(ov).line(p, fill=color + (int(255 * opacity),),
                            width=max(1, int(S_(width))), joint="curve")
    canvas.alpha_composite(ov)

# stem + twigs (gradient by segment midpoint)
stem_d = re.search(r'<path d="(M73 205[^"]+)" stroke="url\(#stem\)" stroke-width="4"', txt).group(1)
draw_polyline(points_from_path(stem_d), (0, 0, 0), 4)  # placeholder replaced below
canvas = Image.new("RGBA", (SZ, SZ), (0, 0, 0, 0))
d = ImageDraw.Draw(canvas)

stem_pts = points_from_path(stem_d)
for i in range(len(stem_pts) - 1):
    col = grad_color(STEM_STOPS, S_P1, S_P2,
                     ((stem_pts[i][0] + stem_pts[i + 1][0]) / 2,
                      (stem_pts[i][1] + stem_pts[i + 1][1]) / 2))
    draw_polyline(stem_pts[i:i + 2], col, 4)

twig_block = re.search(r'<g stroke="url\(#stem\)" stroke-width="2"[^>]*>(.*?)</g>', txt, re.S).group(1)
for m in re.finditer(r'<path d="([^"]+)"/>', twig_block):
    pts = points_from_path(m.group(1))
    for i in range(len(pts) - 1):
        col = grad_color(STEM_STOPS, S_P1, S_P2,
                         ((pts[i][0] + pts[i + 1][0]) / 2,
                          (pts[i][1] + pts[i + 1][1]) / 2))
        draw_polyline(pts[i:i + 2], col, 2)

# leaves: <use href="#naturalLeaf" transform="..."/>
blade = points_from_path(re.search(r'<path id="bladeShape" d="([^"]+)"/>', txt).group(1))
vein = points_from_path("M0-1C-2-20 2-40 1-62")
side_block = re.search(r'<g stroke="#c0d18f"[^>]*>(.*?)</g>', txt, re.S).group(1)
sides = [points_from_path(m.group(1)) for m in re.finditer(r'<path d="([^"]+)"/>', side_block)]
highlight = points_from_path("M-12-35C-12-44-5-54 0-61")

for m in re.finditer(r'<use href="#naturalLeaf" transform="([^"]+)"/>', txt):
    M = parse_use_transform(m.group(1))
    Minv = np.linalg.inv(M)

    # fill: gradient evaluated in local coords via inverse map
    fb = apply(M, blade)
    xs = [p[0] for p in fb]; ys = [p[1] for p in fb]
    x0, x1 = int(min(xs) * S) - 2, int(max(xs) * S) + 3
    y0, y1 = int(min(ys) * S) - 2, int(max(ys) * S) + 3
    gx, gy = np.meshgrid(np.arange(x0, x1) / S, np.arange(y0, y1) / S)
    H = np.vstack([gx.ravel(), gy.ravel(), np.ones(gx.size)])
    L = Minv @ H
    lx, ly = L[0].reshape(gx.shape), L[1].reshape(gx.shape)
    dvec = B_P2 - B_P1
    t = ((lx - B_P1[0]) * dvec[0] + (ly - B_P1[1]) * dvec[1]) / (dvec @ dvec)
    t = np.clip(t, 0, 1)
    arr = np.zeros((gx.shape[0], gx.shape[1], 4), dtype=np.uint8)
    stops = [(o, np.array(hex2rgb(c))) for o, c in BLADE_STOPS]
    ts = t.ravel()
    cols = np.empty((ts.size, 3))
    for i in range(len(stops) - 1):
        (o0, c0), (o1, c1) = stops[i], stops[i + 1]
        msk = (ts >= o0) & (ts <= o1)
        f = np.where(o1 > o0, (ts[msk] - o0) / (o1 - o0), 0)
        cols[msk] = c0[None, :] + (c1 - c0)[None, :] * f[:, None]
    cols[ts <= stops[0][0]] = stops[0][1]
    cols[ts >= stops[-1][0]] = stops[-1][1]
    arr[..., :3] = cols.reshape(gx.shape + (3,)).astype(np.uint8)
    arr[..., 3] = 255
    patch = Image.fromarray(arr, "RGBA")
    mask = Image.new("L", patch.size, 0)
    ImageDraw.Draw(mask).polygon(
        [(S_(x) - x0, S_(y) - y0) for x, y in fb], fill=255)
    patch.putalpha(mask)
    canvas.alpha_composite(patch, (x0, y0))

    # strokes (width scaled by mean scale factor)
    sx = np.linalg.norm(M[:2, 0]); sy = np.linalg.norm(M[:2, 1])
    wfac = (sx + sy) / 2
    draw_polyline(apply(M, vein), hex2rgb("#d0d9a0"), 0.85 * wfac, 0.75)
    for sp in sides:
        draw_polyline(apply(M, sp), hex2rgb("#c0d18f"), 0.55 * wfac, 0.55)
    draw_polyline(apply(M, highlight), hex2rgb("#d9dfb4"), 0.65 * wfac, 0.5)

canvas = canvas.resize((SIZE, SIZE), Image.LANCZOS)
canvas.save(OUT, optimize=True)
print("saved", OUT, OUT.stat().st_size // 1024, "KB")

# white-bg preview
bg = Image.new("RGBA", (SIZE, SIZE), (255, 255, 255, 255))
bg.paste(canvas, (0, 0), canvas)
bg.convert("RGB").save("/tmp/brand_icon_preview.png")
