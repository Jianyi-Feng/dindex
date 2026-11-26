"""One house style for every figure of the giant-virus manuscript (Nature production spec).

Rules encoded here, so that figures cannot drift apart again:
  * width 183 mm (double column) or 89 mm; height <= 170 mm
  * Arial, all text 5-7 pt at final size; panel labels 8 pt bold lowercase
  * open axes (left + bottom spines), 0.5 pt lines, ticks out, no gridlines, no boxed frames
  * EVERY HUE HAS ONE MEANING (V61, 'vivid' theme): Imitervirales blue, Algavirales vermilion (Okabe-Ito);
    the 2014-2015 heatwave purple (#AA3377) with a pale purple band; hosts (18S) green (#1B8A5A);
    temperature as a covariate in viridis. Everything else (models, size fractions, basins) is slate/navy by marker shape.
  * exponents as x10^n, never 'e-15'; P values italic P; one decimal style per panel
"""
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import glob, os
# Arial is resolved from the usual Linux/macOS/Windows font directories. Without Arial the style falls
# back to DejaVu Sans and figures/v73/qa_figures.py will report a non-Arial font.
for _d in ('/usr/share/fonts/truetype/msttcorefonts', '/usr/local/share/fonts', '/Library/Fonts',
           os.path.expanduser('~/.fonts'), os.path.expanduser('~/.local/share/fonts'), 'C:/Windows/Fonts'):
    for _f in glob.glob(os.path.join(_d, '[Aa]rial*.ttf')):
        fm.fontManager.addfont(_f)
MM = 1 / 25.4
IMI, ALG = '#0072B2', '#D55E00'
PRE, EVT, BAND = '#9A9A9A', '#1A1A1A', '#F6E3D3'
INK, MUTED, LIGHT = '#1A1A1A', '#6B6B6B', '#D9D9D9'
# size fractions and Earth-system models are told apart by marker shape and grey level, never by hue
FRAC = {'viral': dict(marker='o', color='#1A1A1A', label='Viral'), 'pico': dict(marker='s', color='#7A7A7A', label='Pico'), 'micro': dict(marker='^', color='#B4B4B4', label='Micro/nano')}
ESM = {'CESM2-WACCM': dict(marker='o', ls='-'), 'IPSL-CM6A-LR': dict(marker='s', ls=(0, (4, 1.5))), 'NorESM2-MM': dict(marker='^', ls=(0, (1, 1))), 'UKESM1-0-LL': dict(marker='D', ls=(0, (4, 1.5, 1, 1.5)))}
GREY3 = ['#4D4D4D', '#C8C8C8', '#1A1A1A']
RC = {'font.family': 'sans-serif', 'font.sans-serif': ['Arial'], 'font.size': 6.5, 'axes.labelsize': 6.5, 'axes.titlesize': 7,
      'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.titlepad': 3, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
      'legend.fontsize': 6, 'legend.frameon': False, 'legend.handlelength': 1.4, 'legend.handletextpad': 0.5, 'legend.borderaxespad': 0.2,
      'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2.2, 'ytick.major.size': 2.2,
      'xtick.direction': 'out', 'ytick.direction': 'out', 'axes.spines.top': False, 'axes.spines.right': False,
      'lines.linewidth': 0.9, 'lines.markersize': 3, 'patch.linewidth': 0.4, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
      'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular',
      'figure.dpi': 150, 'savefig.dpi': 600, 'axes.unicode_minus': True}
def cmap_orders():
    """diverging map for D and for changes in D: Algavirales vermilion <- white -> Imitervirales blue"""
    return mpl.colors.LinearSegmentedColormap.from_list('alg_imi', [ALG, '#F7F7F7', IMI])
CMAP_T = 'Greys'      # temperature (and any other covariate) is shown in greys
HOST = INK            # host (18S) series: neutral in the classic theme
# ---- optional theme (GV_THEME=vivid): every hue still has ONE meaning, but four more things get a hue of their own.
#      orders blue / vermilion (unchanged); the 2014-2015 event purple; hosts green; temperature a sequential map; neutrals cool slate.
import os as _os
THEME = _os.environ.get('GV_THEME', 'vivid')   # V61: the vivid theme is the house scheme; GV_THEME=classic restores the grey one
N_DARK, N_MID, N_LIGHT = '#4D4D4D', '#9A9A9A', '#DDDDDD'
if THEME == 'vivid':
    EVT, BAND = '#AA3377', '#F4E3EC'
    PRE, INK, MUTED, LIGHT = '#98A2B3', '#1F2A44', '#5F6B7A', '#D5DAE3'
    HOST = '#1B8A5A'
    FRAC = {'viral': dict(marker='o', color='#1F2A44', label='Viral'), 'pico': dict(marker='s', color='#6C7A90', label='Pico'), 'micro': dict(marker='^', color='#AAB4C3', label='Micro/nano')}
    GREY3 = ['#5F6B7A', '#D5DAE3', '#1F2A44']
    CMAP_T = 'viridis'
    N_DARK, N_MID, N_LIGHT = '#4A5568', '#98A2B3', '#D5DAE3'
SCALE = 1.18          # global text scale for the 2026-11 size pass; per-figure overrides call use(scale)
def fs(x):
    """Scale a font size written against the original 6.5 pt body style."""
    return x * SCALE
def use(scale=None):
    global SCALE
    if scale is not None: SCALE = scale
    rc = dict(RC)
    for _k in ('font.size', 'axes.labelsize', 'axes.titlesize', 'xtick.labelsize', 'ytick.labelsize', 'legend.fontsize'):
        rc[_k] = RC[_k] * SCALE
    mpl.rcParams.update(rc)
    import matplotlib._mathtext as _mt
    _mt.SHRINK_FACTOR = 0.85          # sub/superscripts at 6.5 pt body -> 5.2 pt, above Nature's 5 pt floor (default 0.7 gives 4.55)
def figure(w_mm=183, h_mm=120): return plt.figure(figsize=(w_mm * MM, h_mm * MM))
_LABELS = {}
def label(ax, s, dx=None, dy=1.02, x_mm=None):
    """8 pt bold lowercase panel label. Placed at save time, just left of the y-axis label text,
    so long axis labels (or larger type) can never collide with the letter."""
    _LABELS.setdefault(id(ax.figure), []).append((ax, s))
def _place_labels(fig):
    items = _LABELS.pop(id(fig), [])
    if not items:
        return
    for ax, s in items:
        bb = ax.get_position()
        fig.text(bb.x0 - 4.0 / 183.0, bb.y1 + 0.012, s, fontsize=fs(8), fontweight='bold', va='bottom', ha='left')
def raincloud(ax, data, pos, color, width=0.34, seed=0, cloud=True, box=True, ms=2.2, side=1):
    """Raincloud: half-violin (the cloud) on one side, a slim box plot on the centre line, and the raw observations (the rain)
    jittered on the other side. `side` = +1 puts the cloud to the right. Everything is drawn in one colour."""
    import numpy as np
    from scipy.stats import gaussian_kde
    v = np.asarray(data, float); v = v[np.isfinite(v)]; rng = np.random.default_rng(seed); gap = 0.05 * width / 0.34
    if cloud and len(v) >= 10 and np.ptp(v) > 0:
        yy = np.linspace(v.min(), v.max(), 200); dens = gaussian_kde(v)(yy); dens = dens / dens.max() * (width - gap)
        ax.fill_betweenx(yy, pos + side * gap, pos + side * (gap + dens), color=color, alpha=0.28, lw=0); ax.plot(pos + side * (gap + dens), yy, color=color, lw=0.6)
    if len(v):
        ax.scatter(pos - side * (gap + rng.uniform(0, width - 2 * gap, len(v))), v, s=ms, color=color, alpha=0.5, lw=0, zorder=2)
    if box and len(v) >= 5:
        q1, md, q3 = np.percentile(v, [25, 50, 75]); iqr = q3 - q1; lo = v[v >= q1 - 1.5 * iqr].min(); hi = v[v <= q3 + 1.5 * iqr].max(); bw = 0.11 * width / 0.34
        ax.plot([pos, pos], [lo, hi], color=color, lw=0.7, zorder=3); ax.add_patch(mpl.patches.Rectangle((pos - bw / 2, q1), bw, iqr, facecolor='white', edgecolor=color, lw=0.8, zorder=4)); ax.plot([pos - bw / 2, pos + bw / 2], [md, md], color=color, lw=1.3, zorder=5, solid_capstyle='butt')
    elif len(v): ax.plot([pos - 0.06, pos + 0.06], [np.median(v)] * 2, color=color, lw=1.3, zorder=5)
def sci(p):
    """1.2e-15 -> '1.2 x 10^-15' in mathtext"""
    m, e = f'{p:.1e}'.split('e'); return rf'{m} × 10$^{{{int(e)}}}$'
def save(fig, stem):
    _place_labels(fig)
    fig.savefig(stem + '.pdf'); fig.savefig(stem + '.png', dpi=300)
