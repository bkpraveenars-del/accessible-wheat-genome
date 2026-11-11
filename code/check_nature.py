"""Audit an SVG against Nature figure specs. python check_nature.py fig.svg"""
import re, sys, itertools, numpy as np
from lxml import etree

PT_PER_MM = 72 / 25.4
OK_WIDTHS = {89: "single", 120: "1.5-col (verify journal)", 183: "double"}
FAIL, WARN = [], []

def to_pt(v):
    m = re.match(r"([\d.]+)\s*(pt|mm|px|in|cm)?", v); x, u = float(m[1]), m[2] or "px"
    return x * {"pt": 1, "px": .75, "mm": PT_PER_MM, "in": 72, "cm": 10 * PT_PER_MM}[u]

def hex2rgb(h):
    h = h.lstrip("#"); h = "".join(c * 2 for c in h) if len(h) == 3 else h
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]) / 255

def lin(c): return np.where(c <= .04045, c / 12.92, ((c + .055) / 1.055) ** 2.4)
def lab(rgb):
    x = lin(rgb) @ np.array([[.4124, .3576, .1805], [.2126, .7152, .0722], [.0193, .1192, .9505]]).T
    x = x / np.array([.95047, 1, 1.08883]); f = np.where(x > .008856, np.cbrt(x), 7.787 * x + 16 / 116)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])
# Machado et al. 2009 CVD matrices, severity 1.0
CVD = {"deuteranopia": np.array([[.367322, .860646, -.227968], [.280085, .672501, .047413], [-.011820, .042940, .968881]]),
       "protanopia": np.array([[.152286, 1.052583, -.204868], [.114503, .786281, .099216], [-.003882, -.048116, 1.051998]])}
def sim(rgb, M):
    l = np.clip(M @ lin(rgb), 0, 1); return np.where(l <= .0031308, 12.92 * l, 1.055 * l ** (1 / 2.4) - .055)
def dE(a, b): return float(np.linalg.norm(lab(a) - lab(b)))

def style_of(el):
    d = {}
    for k in ("fill", "stroke", "stroke-width", "font-size", "font-family", "font-weight"):
        if el.get(k): d[k] = el.get(k)
    for part in (el.get("style") or "").split(";"):
        if ":" in part:
            k, v = part.split(":", 1); d[k.strip()] = v.strip()
    if "font" in d:
        m = re.search(r"([\d.]+)px", d["font"]); d.setdefault("font-size", m[1] + "px" if m else None)
        if "bold" in d["font"] or "700" in d["font"]: d["font-weight"] = "bold"
        d.setdefault("font-family", d["font"].split("px", 1)[-1])
    return d

def num(v): return float(re.match(r"[\d.]+", v)[0])

def main(path):
    root = etree.parse(path).getroot(); ns = "{http://www.w3.org/2000/svg}"
    w, h = to_pt(root.get("width")), to_pt(root.get("height"))
    vb = [float(v) for v in root.get("viewBox", f"0 0 {w} {h}").replace(",", " ").split()]
    scale = w / vb[2]                                   # pt per user unit
    wmm, hmm = w / PT_PER_MM, h / PT_PER_MM
    if not [k for k in OK_WIDTHS if abs(wmm - k) <= 1.5]: WARN.append(f"width {wmm:.1f} mm is not 89/120/183 mm")
    if hmm > 171: FAIL.append(f"height {hmm:.1f} mm > 170 mm max")
    if root.findall(f".//{ns}linearGradient") or root.findall(f".//{ns}radialGradient"): FAIL.append("gradients present")
    if root.findall(f".//{ns}filter"): FAIL.append("filters (shadows/blur) present")
    if root.findall(f".//{ns}pattern"): WARN.append("pattern fills present (Nature avoids patterns)")
    glyphs = [p for p in root.iter(f"{ns}path") if re.match(r".*-[0-9a-f]{2,}$", p.get("id") or "")]
    if glyphs: FAIL.append(f"{len(glyphs)} glyph paths: text is outlined, must be live <text>")
    texts = list(root.iter(f"{ns}text")); n_lab = 0
    if not texts: FAIL.append("no <text> elements: text outlined or missing")
    for t in texts:
        s = {}
        for anc in list(reversed(list(t.iterancestors()))) + [t]: s.update({k: v for k, v in style_of(anc).items() if v})
        content = re.sub(r"\s+", "", "".join(t.itertext()))
        fs = s.get("font-size") or next((style_of(c).get("font-size") for c in t.iter() if style_of(c).get("font-size")), None)
        size = num(fs) * scale if fs else None
        fam = (s.get("font-family") or "").lower()
        bold = s.get("font-weight") in ("bold", "700")
        if fam and not re.search(r"arial|helvetica|sans", fam): FAIL.append(f"font '{fam}' on '{content[:20]}'")
        if size is None: WARN.append(f"no font-size on '{content[:20]}'"); continue
        if bold and len(content) == 1 and content.islower():
            n_lab += 1
            if abs(size - 8) > .3: FAIL.append(f"panel label '{content}' is {size:.1f} pt (need 8)")
        elif not (4.95 <= size <= 7.05): FAIL.append(f"text '{content[:25]}' is {size:.2f} pt (5-7 pt)")
        fill = s.get("fill", "#000000").lower()
        if fill.startswith("#") and fill not in ("#000", "#000000", "#fff", "#ffffff", "#000000ff"):
            FAIL.append(f"coloured text {fill} on '{content[:20]}' (use black + keyline)")
    for ts in root.iter(f"{ns}tspan"):
        fs = style_of(ts).get("font-size")
        if fs and num(fs) * scale < 4.95:
            txt = "".join(ts.itertext()).strip()[:10]
            FAIL.append(f"sub/superscript '{txt}' is {num(fs) * scale:.2f} pt (<5 pt): call nf.apply(), base text 6 pt")
    if n_lab == 0: WARN.append("no 8 pt bold lowercase panel labels found")
    widths = []
    for el in root.iter():
        st = style_of(el); sw = st.get("stroke-width")
        if sw and st.get("stroke", "") not in ("none", ""):
            try: widths.append(num(sw) * scale)
            except Exception: pass
    thin = [x for x in widths if 0 < x < .25]
    if thin: WARN.append(f"{len(thin)} strokes thinner than 0.25 pt (may vanish in print)")
    cols = set()
    for el in root.iter():
        if el.tag == f"{ns}text": continue
        st = style_of(el)
        for k in ("fill", "stroke"):
            v = (st.get(k) or "").lower()
            if re.fullmatch(r"#[0-9a-f]{6}|#[0-9a-f]{3}", v): cols.add(v if len(v) == 7 else "#" + "".join(c * 2 for c in v[1:]))
    data_cols = [c for c in cols if not (np.ptp(hex2rgb(c)) < .06)]   # drop greys/black/white
    if len(data_cols) <= 40:
        for a, b in itertools.combinations(data_cols, 2):
            A, B = hex2rgb(a), hex2rgb(b)
            if dE(A, B) < 20: continue
            for name, M in CVD.items():
                if dE(sim(A, M), sim(B, M)) < 8:
                    FAIL.append(f"{a} vs {b} indistinguishable under {name}"); break
    else: WARN.append(f"{len(data_cols)} colours (colormap?): CVD pair test skipped; use viridis/cividis")
    print(f"{path}: {wmm:.1f} x {hmm:.1f} mm, {len(texts)} text, {n_lab} panel labels, {len(data_cols)} data colours")
    for m in dict.fromkeys(FAIL): print("  FAIL", m)
    for m in dict.fromkeys(WARN): print("  WARN", m)
    print("  RESULT:", "FAIL" if FAIL else "PASS")
    return 1 if FAIL else 0

if __name__ == "__main__": sys.exit(main(sys.argv[1]))
