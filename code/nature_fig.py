"""Nature-compliant matplotlib helpers. Usage: import nature_fig as nf; nf.apply()."""
import re, matplotlib as mpl, matplotlib.pyplot as plt
from cycler import cycler

MM = 1 / 25.4
WIDTHS = {"single": 89, "onehalf": 120, "double": 183}   # mm; 120 is journal-specific
MAX_H = 170                                             # mm, Nature max height
WONG = dict(black="#000000", orange="#E69F00", skyblue="#56B4E9", green="#009E73",
            yellow="#F0E442", blue="#0072B2", vermillion="#D55E00", purple="#CC79A7")
CYCLE = [WONG[k] for k in ("blue", "orange", "green", "vermillion", "skyblue", "purple", "black")]
GREY = "#7F7F7F"

def apply():
    # mathtext shrinks scripts to 70 % (4.2 pt on 6 pt text, below Nature's 5 pt minimum); 0.84 keeps them >= 5 pt
    import matplotlib._mathtext as _mt; _mt.SHRINK_FACTOR = 0.84
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 6, "axes.labelsize": 6, "axes.titlesize": 6,
        "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
        "mathtext.default": "it", "mathtext.fontset": "custom",
        "mathtext.it": "sans:italic", "mathtext.rm": "sans", "mathtext.cal": "sans",
        "mathtext.bf": "sans:bold", "mathtext.sf": "sans",
        "text.color": "black", "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black", "axes.edgecolor": "black",
        "axes.linewidth": 0.5, "lines.linewidth": 1.0, "lines.markersize": 3,
        "patch.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2, "ytick.major.size": 2, "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "legend.frameon": False, "legend.handlelength": 1.2, "legend.borderaxespad": 0.2,
        "axes.prop_cycle": cycler(color=CYCLE), "image.cmap": "viridis",
        "figure.facecolor": "white", "savefig.facecolor": "white", "savefig.transparent": False,
        "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 450, "savefig.bbox": None, "figure.dpi": 150,
    })

def sup(x):
    """Upright superscript for units: 'kg ha' + sup('−1'). Variables stay $italic$."""
    return r"$^{\mathregular{%s}}$" % x

def figure(width="double", height_mm=100, **kw):
    w = WIDTHS[width] if isinstance(width, str) else width
    if height_mm > MAX_H:
        raise ValueError(f"height {height_mm} mm exceeds Nature max {MAX_H} mm")
    return plt.figure(figsize=(w * MM, height_mm * MM), **kw)

def label_panels(fig, axes, letters="abcdefghijklmnop", dx_mm=-6, dy_mm=1.5):
    """8 pt bold lowercase upright labels at the top-left of each panel."""
    fw, fh = fig.get_size_inches() / MM
    for ax, l in zip(axes, letters):
        b = ax.get_position()
        fig.text(b.x0 + dx_mm / fw, b.y1 + dy_mm / fh, l, fontsize=8, fontweight="bold",
                 ha="left", va="bottom")

def save(fig, stem, family="Arial"):
    """Writes stem.svg (live text), stem.pdf (Type 42), stem_preview.png (300 dpi)."""
    fig.savefig(stem + ".svg"); fig.savefig(stem + ".pdf"); fig.savefig(stem + "_preview.png", dpi=300)
    s = open(stem + ".svg", encoding="utf8").read()
    s = re.sub(r"font-family:[^;\"]*", f"font-family: {family}, Helvetica, sans-serif", s)
    open(stem + ".svg", "w", encoding="utf8").write(s)
