# -*- coding: utf-8 -*-
# 基于 big115（1.15× 字号版）重做的修改；原文件未动。
"""Main and Extended Data figures for the  D-index manuscript ( builder +  Fig. 5, ED Figs 4-5).

Nature specification: 180 mm width, Arial-metric sans serif (Liberation Sans), 6-8 pt text,
8 pt bold lower-case panel letters, 0.5 pt axes, Okabe-Ito colour-blind-safe palette,
vector PDF (TrueType fonts embedded) + 600-dpi PNG + per-figure source data (CSV).

Inputs: ../raw_tables (supplied tables) and analysis (outputs of dindex_analyses.py).
Run with a Python that has cartopy (job venv) for the Tara map.
Optional genome tables, placed in ../raw_tables:
  GVMAG_traits_long.csv    genome_id, order, trait, present
  GVMAG_microdiversity.csv genome_id, order, size_fraction, mean_pi
When absent, Fig. 4b,c use labelled crops of the previous raster panels.
"""
from pathlib import Path
import sys
import io, zipfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
import statsmodels.formula.api as smf
from scipy import stats
from PIL import Image

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
sys.path.insert(0, str(PKG / "tools"))
import paths as P  # noqa: E402  (repository layout is defined in tools/paths.py)

DATA = P.RAW
AN = P.MODEL_STATS
A29 = P.MODEL_STATS
MR = P.GENOMES                     # marker-recovery tables live in genomes_and_microdiversity/
SRC = P.outdir(P.OUT_SOURCE_DATA)
OUTFIG = P.outdir(P.OUT)
MM = 1 / 25.4
RNG = np.random.default_rng(7)

# ------------------------------------------------------------------ style
plt.rcParams.update({
    # real Arial for text AND mathtext (was Liberation Sans + DejaVu Sans mixed)
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial", "mathtext.tt": "Arial", "mathtext.cal": "Arial:italic",
    "mathtext.default": "regular",
    "font.size": 7.705, "axes.titlesize": 8.28, "axes.labelsize": 7.705, "xtick.labelsize": 7.13,
    "ytick.labelsize": 7.13, "legend.fontsize": 7.13, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.9, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "axes.titlepad": 3, "axes.labelpad": 2, "savefig.dpi": 600,
})
# ---------------------------------------------------------------- 小提琴垫片（第二轮反馈）
# 只在这三段里把 boxplot 临时替换为小提琴图；其余箱线图面板保持原样。
import matplotlib.axes as _ma
_ORIG_BOXPLOT = _ma.Axes.boxplot


def _violin_boxplot(self, data, positions=None, widths=None, **kw):
    import numpy as _np
    # 取色：优先系列色（boxprops.edgecolor），medianprops 的黑色只用于中位线
    col = "k"
    for key, field in (("boxprops", "edgecolor"), ("boxprops", "color"),
                       ("whiskerprops", "color"), ("whiskerprops", "edgecolor")):
        d_ = kw.get(key) or {}
        if d_.get(field):
            col = d_[field]
            break
    # 规整输入：Series/ndarray -> 单元素列表；已是 list/tuple（含长度不等的 ragged 列表）原样保留
    if isinstance(data, (list, tuple)):
        data = list(data)
    else:
        data = [data]
    pos = positions if positions is not None else range(1, len(data) + 1)
    w = widths if isinstance(widths, (int, float)) else 0.5
    for vals, p in zip(data, pos):
        v = _np.asarray(vals, float).ravel()
        v = v[_np.isfinite(v)]
        if v.size < 2 or _np.ptp(v) == 0:
            continue
        vp = self.violinplot([v], positions=[p], widths=w, showextrema=False,
                             showmedians=True, bw_method=0.35)
        lo, hi = float(_np.min(v)), float(_np.max(v))
        for b in vp["bodies"]:
            # KDE 尾巴会伸出数据范围，按数据极值裁掉（避免长尖刺）
            for path in b.get_paths():
                xy = path.vertices
                xy[:, 1] = _np.clip(xy[:, 1], lo, hi)
            b.set_facecolor("white"); b.set_edgecolor(col); b.set_linewidth(0.55); b.set_alpha(1)
        vp["cmedians"].set_color("k"); vp["cmedians"].set_linewidth(0.8)
    return None


def violin_panels_on():
    _ma.Axes.boxplot = _violin_boxplot


def violin_panels_off():
    # 作者要求全局小提琴：保持启用（不再还原）
    pass


# 作者要求：本图集内所有箱线图一律改为小提琴图
violin_panels_on()

I_C, A_C = "#0072B2", "#D55E00"          # Imitervirales blue, Algavirales vermillion
I_L, A_L = "#9CC3E0", "#F0B48A"
FR_C = {"viral": "#000000", "pico": "#CC79A7", "micro": "#009E73"}
FR_M = {"viral": "o", "pico": "s", "micro": "^"}
FR_LAB = {"viral": "Viral (0.02–0.2 µm)", "pico": "Pico (0.2–1.2 µm)", "micro": "Micro/nano (1.2–80 µm)"}
INK = "#333333"
GREY, LGREY, HW = "#7F7F7F", "#D9D9D9", "#FBE3CF"
ESM_C = {"CESM2-WACCM": "#E69F00", "IPSL-CM6A-LR": "#56B4E9", "NorESM2-MM": "#009E73", "UKESM1-0-LL": "#CC79A7"}
D_CMAP = LinearSegmentedColormap.from_list("D", [A_C, "#F7F7F7", I_C])
NORM = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
MONTHS = list("JFMAMJJASOND")
MHW0, MHW1 = pd.Timestamp("2014-01-01"), pd.Timestamp("2016-08-31")
KM = 111.32


def letter(ax, s, x=-0.2, y=1.08):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9.2, fontweight="bold", va="bottom", ha="left")


def title(ax, s):
    ax.set_title(s, loc="left", fontweight="bold")


def save(fig, name, sources):
    _fill_canvas_width(fig)
    _tighten_row_gaps(fig, gap_mm=1.2 if name in ("Fig1_derive", "Fig3_validate") else 3.0)
    _crop_canvas_height(fig)
    fig.savefig(OUTFIG / f"{name}.pdf")
    fig.savefig(OUTFIG / f"{name}.svg")
    fig.savefig(OUTFIG / f"{name}.png", dpi=600)
    plt.close(fig)
    with pd.ExcelWriter(SRC / f"{name}_source_data.xlsx") as xw:
        for k, df in sources.items():
            df.to_excel(xw, sheet_name=k[:31], index=False)
    print("wrote", name)



def fmt_p(p):
    if p < 1e-3:
        e = int(np.floor(np.log10(p)))
        return f"$P$ = {p / 10**e:.1f}e{e}"
    return f"$P$ = {p:.2g}"


# ------------------------------------------------------------------ data
spot = pd.read_csv(DATA / "SPOT_sample_Dindex.csv", parse_dates=["date"])
spot["frac"] = spot.size_fraction_ascii.map({"0.02-0.2um": "viral", "0.2-1.2um": "pico", "1.2-80um": "micro"})
spot["I"], spot["A"] = spot.Imitervirales_relative_abundance, spot.Algavirales_relative_abundance
tara = pd.read_csv(DATA / "Tara_observed_Dindex.csv")
tara = tara[tara.D_index_qc == "ok"].copy()
val = pd.read_csv(DATA / "Dindex_model_validation.csv")
heat = pd.read_csv(DATA / "SPOT_heatwave_Dindex_stats.csv")
latyr = pd.read_csv(DATA / "CMIP6_Dindex_latitude_year.csv")
front = pd.read_csv(DATA / "Dindex_front_trajectory.csv")
front = front[front.ESM != "ENSEMBLE-mean"].copy()
front["abslat"] = front.D0_front_latitude.abs()
hs = pd.read_csv(AN / "headline_statistics.csv")
V = pd.read_csv(AN / "inherited__spot_viral_heatwave_residuals.csv", parse_dates=["date"])
cv = pd.read_csv(A29 / "tara_cross_validation.csv")
imp = pd.read_csv(A29 / "tara_oof_permutation_importance.csv")
stn = pd.read_csv(AN / "inherited__tara_station_means.csv")
per_esm = pd.read_csv(AN / "inherited__cmip6_front_per_esm.csv")


def H(name):
    r = hs[hs.statistic == name]
    assert len(r) == 1, name
    return r.iloc[0]


# ================================================================== Fig. 1 DERIVE
def _crop_canvas_height(fig, pad_mm=1.0, min_frac=0.35):
    """Cut the blank bands off the top and bottom of the canvas.

    The page width is fixed by the journal, so the horizontal pass has to grow
    the panels into the margin; the page height is free, so here the paper is
    simply taken away instead.  The canvas is resized and every axes rectangle
    is remapped so that each panel keeps its exact absolute size and aspect
    ratio -- nothing is stretched, the drawing is only translated.  Font sizes,
    line widths and marker sizes are untouched.  If the result would not fit,
    the original canvas and positions are restored.

    Positions are snapshotted for the same reason as in _fill_canvas_width:
    ax.twinx() gives two Axes one shared rectangle.
    """
    try:
        W, H = fig.get_figwidth(), fig.get_figheight()
        pad = pad_mm / 25.4
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        y0, y1 = bb.y0 / H, bb.y1 / H
        span = y1 - y0
        newH = span * H + 2 * pad
        if span <= 0.02 or newH >= H - 0.002 or newH < min_frac * H:
            print("    [crop] left as is: %.1f mm tall, top %.2f mm, bottom %.2f mm"
                  % (H * 25.4, y0 * H * 25.4, (1 - y1) * H * 25.4))
            return fig
        snap = [(ax, ax.get_position().frozen()) for ax in fig.axes]
        fig.set_size_inches(W, newH)
        k = H / newH
        for ax, p in snap:
            ax.set_position([p.x0, (p.y0 - y0) * k + pad / newH, p.width, p.height * k])
        fig.canvas.draw()
        bb2 = fig.get_tightbbox(fig.canvas.get_renderer())
        if bb2.y0 < -0.0005 or bb2.y1 > newH + 0.0005:
            fig.set_size_inches(W, H)
            for ax, p in snap:
                ax.set_position(p)
            fig.canvas.draw()
            print("    [crop] REVERTED: content would not fit")
            return fig
        print("    [crop] height %.1f -> %.1f mm   top %.2f mm   bottom %.2f mm"
              % (H * 25.4, newH * 25.4, bb2.y0 * 25.4, (newH - bb2.y1) * 25.4))
    except Exception as exc:              # cosmetics must never break a render
        print("    [crop] skipped:", type(exc).__name__, exc)
    return fig


def _tighten_row_gaps(fig, gap_mm=3.0, min_gap_mm=5.0):
    """Pull vertically stacked rows of panels together.

    Row spacing was hand-set as hspace -- a fraction of the panel height -- and
    on these figures it reserves far more room than the x labels of one row and
    the titles of the next row actually occupy.  Measure where each panel's ink
    really starts and ends, group the panels into rows by that ink, and
    translate the lower rows up until only gap_mm is left between two rows.
    Panels keep their exact size and aspect ratio; only their y position moves.
    A gap below min_gap_mm is left alone on purpose: that is the deliberate
    offset between a panel and its own inset colour bar, not row spacing.
    Nothing is ever moved apart, and a translation that would leave two rows
    touching is discarded.

    Positions are snapshotted before anything is written, for the same reason
    as in _fill_canvas_width: ax.twinx() gives two Axes one shared rectangle.
    """
    def measure():
        fig.canvas.draw()
        R = fig.canvas.get_renderer()
        raw = []
        inv = fig.dpi_scale_trans.inverted()      # Axes.get_tightbbox is in pixels
        for ax in fig.axes:
            b = ax.get_tightbbox(R).transformed(inv)
            raw.append((b.y0, b.y1, ax))
        raw.sort(key=lambda t: t[0])
        out = []
        for y0, y1, ax in raw:
            if out and y0 <= out[-1][1] + 0.5 / 25.4:
                out[-1][1] = max(out[-1][1], y1)
                out[-1][2].append(ax)
            else:
                out.append([y0, y1, [ax]])
        out.sort(key=lambda b: -b[1])
        return out

    try:
        H = fig.get_figheight()
        bands = measure()
        if len(bands) < 2:
            print("    [rows] single row -- nothing to tighten")
            return fig
        gaps = lambda bs: [(bs[k][0] - bs[k + 1][1]) * 25.4 for k in range(len(bs) - 1)]
        before = gaps(bands)
        cum, shift_of = 0.0, {}
        for k in range(len(bands) - 1):
            if before[k] > min_gap_mm:
                cum += max(0.0, before[k] - gap_mm)
            for ax in bands[k + 1][2]:
                shift_of[id(ax)] = cum
        if cum < 0.2:
            print("    [rows] already tight: gaps %s mm"
                  % " ".join("%.1f" % g for g in before))
            return fig
        if cum > 0.25 * H * 25.4:
            print("    [rows] REFUSED: %.1f mm shift is implausible on a %.1f mm canvas"
                  % (cum, H * 25.4))
            return fig
        snap = [(ax, ax.get_position().frozen()) for ax in fig.axes]
        for ax, p in snap:                       # absolute values, from the snapshot
            ax.set_position([p.x0, p.y0 + shift_of.get(id(ax), 0.0) / 25.4 / H,
                             p.width, p.height])
        after = gaps(measure())
        if len(after) != len(before) or any(g < 0.5 for g in after):
            for ax, p in snap:
                ax.set_position(p)
            fig.canvas.draw()
            print("    [rows] REVERTED: rows would touch")
            return fig
        print("    [rows] gaps %s -> %s mm   (%.1f mm reclaimed)"
              % (" ".join("%.1f" % g for g in before),
                 " ".join("%.1f" % g for g in after), cum))
    except Exception as exc:              # cosmetics must never break a render
        print("    [rows] skipped:", type(exc).__name__, exc)
    return fig


def _fill_canvas_width(fig, iters=8, m=0.005):
    """Grow the layout sideways until the drawn content fills 99 % of the canvas.

    A blank band at the left or right edge of a fixed-width figure is dead
    space: it shrinks every panel without buying room for anything.  Measure the
    real content bbox and affinely remap every axes rectangle about it so the
    grid expands into the margin.  The content is aimed at [m, 1-m] -- 0.9 mm of
    slack at each edge -- rather than exactly at the page edge: touching the
    edge is a knife-edge that font-metric rounding can push a fraction of a
    millimetre over, and clipping is far worse than this 0.5 % inset, which also
    leaves room to nudge elements by hand later.  Font sizes are untouched.

    Every position is snapshotted before anything is written: ax.twinx() gives
    two Axes one shared rectangle, so re-reading a position that an earlier line
    already rewrote would scale that panel twice and push it off the page.
    """
    try:
        W = fig.get_figwidth()
        hist = []
        for _ in range(iters):
            fig.canvas.draw()
            bb = fig.get_tightbbox(fig.canvas.get_renderer())
            x0, x1 = bb.x0 / W, bb.x1 / W
            span = x1 - x0
            hist.append((span, x0, [(ax, ax.get_position().frozen()) for ax in fig.axes]))
            if span <= 0.02 or span > 1.0:
                break                                  # degenerate / already over the edge
            if span >= 1.0 - 2 * m - 0.002 and x0 >= m - 0.004 and x1 <= 1.0 - m + 0.004:
                break                                  # measured, on target, nothing pending
            s = (1.0 - 2 * m) / span
            for ax, p in hist[-1][2]:
                ax.set_position([m + (p.x0 - x0) * s, p.y0, p.width * s, p.height])
        fit = [h for h in hist if h[0] <= 1.0 and h[1] >= -0.0005 and h[1] + h[0] <= 1.0005]
        pick = max(fit, key=lambda h: h[0]) if fit else min(hist, key=lambda h: h[0])
        if pick[2] is not hist[-1][2]:                 # keep the widest layout that still fits
            for ax, p in pick[2]:
                ax.set_position(p)
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        print("    [fill] width used {:5.1f}%  blank left {:4.2f} mm  blank right {:4.2f} mm"
              .format((bb.x1 - bb.x0) / W * 100, bb.x0 / W * 180, (1 - bb.x1 / W) * 180))
        print("           trace   " + "   ".join(
            "{:.1f}%[{:+.2f},{:+.2f}]".format(h[0] * 100, h[1] * 180, (1 - h[1] - h[0]) * 180)
            for h in hist))
    except Exception as exc:              # cosmetics must never break a render
        print("    [fill] skipped:", type(exc).__name__, exc)
    return fig


fig = plt.figure(figsize=(180 * MM, 112 * MM))
gs = fig.add_gridspec(2, 3, left=0.10, right=0.92, bottom=0.09, top=0.92, wspace=0.42, hspace=0.55)

ax = fig.add_subplot(gs[0, 0]); letter(ax, "a"); title(ax, "SPOT metagenomes (n = 167)")
ax.axvspan(MHW0, MHW1, color=HW, lw=0, zorder=0)
ax.text(MHW0 + (MHW1 - MHW0) / 2, 3.05, "Marine\nheatwave", ha="center", va="bottom", fontsize=5.75, color="#A0522D")
for yi, f in enumerate(["micro", "pico", "viral"][::-1]):
    g = spot[spot.frac == f]
    ax.scatter(g.date, np.full(len(g), yi + 0.5), s=4, marker=FR_M[f], color=FR_C[f], lw=0)
    ax.text(g.date.max() + pd.Timedelta(days=90), yi + 0.5, f"{len(g)}", va="center", fontsize=5.75, color=GREY)
ax.set_yticks([0.5, 1.5, 2.5], ["Viral", "Pico", "Micro/nano"])
ax.set_ylim(0, 3.6); ax.set_xlim(pd.Timestamp("2009-01-01"), pd.Timestamp("2020-06-01"))
ax.xaxis.set_major_locator(matplotlib.dates.YearLocator(2)); ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)

src1 = {"a_samples": spot[["sample_id", "date", "frac", "I", "A", "D_index", "SST", "SST_anomaly"]]}
for k, (f, pos) in enumerate([("viral", (0, 1)), ("pico", (0, 2)), ("micro", (1, 0))]):
    ax = fig.add_subplot(gs[pos]); letter(ax, "bcd"[k]); title(ax, FR_LAB[f])
    g = spot[spot.frac == f]
    s = g.groupby("month").agg(Im=("I", "median"), Iq1=("I", lambda x: x.quantile(.25)), Iq3=("I", lambda x: x.quantile(.75)),
                               Am=("A", "median"), Aq1=("A", lambda x: x.quantile(.25)), Aq3=("A", lambda x: x.quantile(.75)),
                               n=("I", "size")).reindex(range(1, 13)).reset_index()
    for c, lo, hi, med, lab in [(I_C, "Iq1", "Iq3", "Im", "Imitervirales"), (A_C, "Aq1", "Aq3", "Am", "Algavirales")]:
        ax.fill_between(s.month, s[lo], s[hi], color=c, alpha=0.18, lw=0)
        ax.plot(s.month, s[med], color=c, marker="o", ms=2, lw=0.9, label=lab)
    ax.set_xticks(range(1, 13), MONTHS); ax.set_ylim(0, 1); ax.set_xlim(0.5, 12.5)
    ax.set_ylabel("Relative abundance"); ax.set_xlabel("Month")
    if k == 0:
        ax.legend(loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0), handlelength=1.2, columnspacing=0.8)
    src1[f"{'bcd'[k]}_{f}_monthly"] = s

ax = fig.add_subplot(gs[1, 1]); letter(ax, "e"); title(ax, "D is associated with SST")
g = spot.dropna(subset=["SST"])
m = smf.ols("D_index ~ SST + C(frac)", data=g).fit(cov_type="HC3")
for f in ["viral", "pico", "micro"]:
    gg = g[g.frac == f]
    ax.scatter(gg.SST, gg.D_index, s=5, marker=FR_M[f], facecolor="none", edgecolor=FR_C[f], lw=0.4, alpha=0.8)
    xs = np.linspace(gg.SST.min(), gg.SST.max(), 50)
    ax.plot(xs, m.predict(pd.DataFrame({"SST": xs, "frac": f})), color=FR_C[f], lw=0.9)
ax.axhline(0, color=LGREY, lw=0.5, ls="--", zorder=0)
b = H("SPOT dD/dSST, fraction-adjusted (HC3)")
ax.text(0.03, 0.97, f"{b.value:.2f} D/°C (95% CI {b.ci_low:.2f}–{b.ci_high:.2f})\n$r$ = {H('SPOT D vs SST Pearson r').value:.2f}, n = {int(b.n)}",
        bbox=dict(facecolor="white", alpha=1, edgecolor="none", pad=1),
        transform=ax.transAxes, va="top", fontsize=6.325)
ax.set_xlabel("SST (°C)"); ax.set_ylabel("D-index"); ax.set_ylim(-1.05, 1.05)
ax.legend(handles=[Line2D([], [], marker=FR_M[f], ls="-", color=FR_C[f], mfc="none", ms=3, lw=0.8, label=l)
                   for f, l in [("viral", "Viral"), ("pico", "Pico"), ("micro", "Micro/nano")]],
          loc="lower right", handlelength=1.5, fontsize=5.75, frameon=True, facecolor="white", edgecolor="none", framealpha=1)

ax = fig.add_subplot(gs[1, 2]); letter(ax, "f"); title(ax, "The index: D = (I − A)/(I + A)")
for f in ["viral", "pico", "micro"]:
    s = spot[spot.frac == f].groupby("month").D_index.mean().reindex(range(1, 13))
    ax.plot(s.index, s.values, color=FR_C[f], lw=0.5, alpha=0.6)
allm = spot.groupby("month").D_index.agg(["mean", "sem", "size"]).reindex(range(1, 13))
tcrit = stats.t.ppf(0.975, allm["size"] - 1)
ax.fill_between(allm.index, allm["mean"] - tcrit * allm["sem"], allm["mean"] + tcrit * allm["sem"], color=GREY, alpha=0.2, lw=0)
ax.plot(allm.index, allm["mean"], color="k", lw=1.3, marker="o", ms=2)
ax.axhline(0, color=LGREY, lw=0.5, ls="--", zorder=0)
ax.set_ylim(-1.05, 1.05); ax.set_xticks(range(1, 13), MONTHS); ax.set_xlim(0.5, 12.5)
ax.set_ylabel("D-index"); ax.set_xlabel("Month")
ax.text(12.4, 0.93, "Imitervirales-dominated", ha="right", va="top", color=I_C, fontsize=5.75)
ax.text(12.4, -0.93, "Algavirales-dominated", ha="right", va="bottom", color=A_C, fontsize=5.75)
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
clim = spot.drop_duplicates("env_month").groupby("month").SST.mean()
ax2.plot(clim.index, clim.values, color=GREY, ls="--", lw=0.8)
ax2.set_ylabel("SST climatology (°C)", color=GREY); ax2.tick_params(axis="y", colors=GREY)
ax2.set_ylim(12, 24)
src1["e_regression"] = pd.DataFrame({"term": m.params.index, "coef": m.params.values, "P": m.pvalues.values})
src1["f_monthly_D"] = allm.reset_index().assign(SST_climatology=clim.reindex(range(1, 13)).values)
save(fig, "Fig1_derive", src1)

# ================================================================== Fig. 2 TEST
fig = plt.figure(figsize=(180 * MM, 105 * MM))
gs = fig.add_gridspec(2, 7, left=0.09, right=0.955, bottom=0.11, top=0.92, wspace=1.6, hspace=0.6)
v = V[V.date <= "2015-12-31"].copy()
pre, hw = v[v.date < MHW0], v[v.date >= MHW0]

ax = fig.add_subplot(gs[0, :4]); letter(ax, "a", x=-0.1); title(ax, "Northeast Pacific marine heatwave at SPOT")
ax.axvspan(MHW0, pd.Timestamp("2015-12-31"), color=HW, lw=0, zorder=0)
ax.bar(v.date, v.SST_anomaly, width=25, color=np.where(v.SST_anomaly > 0, A_C, I_C), lw=0)
ax.axhline(0, color="k", lw=0.4, ls="--")
ax.set_ylabel("SST anomaly (°C)"); ax.set_ylim(-3.3, 3.8)
ax2 = ax.twinx(); ax2.spines["right"].set_visible(True)
ax2.plot(v.date, v.MEI, color=GREY, lw=0.8); ax2.text(1.0, 1.02, "MEI", transform=ax2.transAxes, ha="center", va="bottom", color=GREY, fontsize=6.9); ax2.tick_params(axis="y", colors=GREY); ax2.set_ylim(-2.5, 3)
ax.xaxis.set_major_locator(matplotlib.dates.YearLocator()); ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.text(0.01, 0.97, f"Mean anomaly: {pre.SST_anomaly.mean():+.1f} °C (2009–2013), {hw.SST_anomaly.mean():+.1f} °C (2014–2015)",
        transform=ax.transAxes, va="top", fontsize=5.75)

ax = fig.add_subplot(gs[0, 4:]); letter(ax, "b", x=-0.14); title(ax, "Viral-fraction D against its\npre-event seasonal cycle")
ax.axvspan(MHW0, pd.Timestamp("2015-12-31"), color=HW, lw=0, zorder=0)
full_t = pd.date_range(v.date.min(), v.date.max(), freq="MS")
harm = smf.ols("D_index ~ s + c", data=pre).fit()
curve = harm.predict(pd.DataFrame({"s": np.sin(2 * np.pi * full_t.month / 12), "c": np.cos(2 * np.pi * full_t.month / 12)}))
ax.plot(full_t, curve, color=GREY, ls="--", lw=0.8, label="Expected (pre-2014 harmonic fit)")
ax.plot(v.date, v.D_index, color="k", lw=0.5)
ax.scatter(v.date, v.D_index, s=5, c=np.where(v.date >= MHW0, A_C, "k"), lw=0, zorder=3, label="Observed")
ax.axhline(0, color=LGREY, lw=0.5, ls="--", zorder=0)
ax.set_ylim(-1.05, 1.05); ax.set_ylabel("D-index")
ax.xaxis.set_major_locator(matplotlib.dates.YearLocator()); ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%Y"))
ax.legend(loc="upper left", ncol=1, handlelength=1.5, borderaxespad=0.2)

ax = fig.add_subplot(gs[1, :2]); letter(ax, "c"); title(ax, "Seasonal\nresiduals")
for xi, (vals, col) in enumerate([(pre.resid, GREY), (hw.resid, A_C)]):
    bp = ax.boxplot(vals, positions=[xi], widths=0.5, showfliers=False, patch_artist=True,
                    medianprops=dict(color="k", lw=0.8), boxprops=dict(lw=0.5, facecolor="white"),
                    whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    ax.scatter(xi + RNG.normal(0, 0.07, len(vals)), vals, s=4, color=col, lw=0, alpha=0.8, zorder=3)
ax.axhline(0, color=LGREY, lw=0.5, ls="--", zorder=0)
ax.set_xticks([0, 1], [f"2009–2013\n(n = {len(pre)})", f"2014–2015\n(n = {len(hw)})"])
ax.set_ylabel("Observed − expected D")

ax = fig.add_subplot(gs[1, 2:5]); letter(ax, "d", x=-0.42); title(ax, "Heatwave-associated change of D")
sup = heat.set_index("contrast")
dm = H("matched-month D difference, median")
diffs = (hw.D_index - hw.month.map(pre.groupby("month").D_index.mean())).dropna().to_numpy()
mm_boot = [np.median(RNG.choice(diffs, len(diffs))) for _ in range(5000)]
yb = H("season-adjusted mean shift, year-block bootstrap CI")
forest = pd.DataFrame([
    ("Year-level mean", yb.value, yb.ci_low, yb.ci_high, f"$P_{{exact}}$ = {H('year-label permutation P (one-sided, exact)').value:.3f}", "k"),
    ("Sample-level median\n(season-adjusted)", sup.loc["pre_vs_heatwave_2014_2015", "delta_D_season_adjusted"], sup.loc["pre_vs_heatwave_2014_2015", "delta_D_season_adj_CI_low"], sup.loc["pre_vs_heatwave_2014_2015", "delta_D_season_adj_CI_high"], f"$P_{{BH}}$ = {sup.loc['pre_vs_heatwave_2014_2015', 'P_season_adjusted_BH']:.4f}", INK),
    ("Same-month\ncomparison", dm.value, np.percentile(mm_boot, 2.5), np.percentile(mm_boot, 97.5), f"$P_{{bootstrap}}$ = {dm.P:.4f}", GREY),
    ("2014 only", sup.loc["pre_vs_2014", "delta_D_season_adjusted"], sup.loc["pre_vs_2014", "delta_D_season_adj_CI_low"], sup.loc["pre_vs_2014", "delta_D_season_adj_CI_high"], f"$P_{{BH}}$ = {sup.loc['pre_vs_2014', 'P_season_adjusted_BH']:.4f}", INK),
    ("2015 only", sup.loc["pre_vs_2015", "delta_D_season_adjusted"], sup.loc["pre_vs_2015", "delta_D_season_adj_CI_low"], sup.loc["pre_vs_2015", "delta_D_season_adj_CI_high"], f"$P_{{BH}}$ = {sup.loc['pre_vs_2015', 'P_season_adjusted_BH']:.4f}", GREY),
    ("Unadjusted median", sup.loc["pre_vs_heatwave_2014_2015", "delta_D"], sup.loc["pre_vs_heatwave_2014_2015", "CI_low"], sup.loc["pre_vs_heatwave_2014_2015", "CI_high"], f"$P_{{BH}}$ = {sup.loc['pre_vs_heatwave_2014_2015', 'P_adjusted_BH']:.2f}", GREY),
], columns=["estimate", "dD", "lo", "hi", "p", "col"])
forest["col"] = "#000000"   # 作者要求：森林图所有行统一黑色
yy = np.arange(len(forest))[::-1]
for y_, r in zip(yy, forest.itertuples()):
    ax.plot([r.lo, r.hi], [y_, y_], color=r.col, lw=0.8)
    ax.plot(r.dD, y_, "o", ms=3, color=r.col, lw=0)
    ax.text(1.32, y_, r.p, va="center", fontsize=6.67, color="k")
ax.axvline(0, color=LGREY, lw=0.6, ls="--", zorder=0)
ax.set_yticks(yy, forest.estimate, fontsize=5.75); ax.set_xlim(-0.2, 1.72); ax.set_xticks([0, 0.4, 0.8, 1.2]); ax.spines["bottom"].set_bounds(-0.2, 1.25); ax.set_ylim(-0.6, len(forest) - 0.4)
ax.set_xlabel("ΔD, 2014–2015 minus 2009–2013 (95% CI)", x=0.4)

# [panel-width FIG2-D] d: 14 % off the left, 4 % off the right.  Measured in
# the rendered PDF: d's widest y tick label ran from x=111 to x=176 pt while
# panel c ends at x=136 pt, i.e. 8.6 mm inside c; the panel letter 'd' sat
# 10.1 mm inside c.  Moving d's left edge right by 24.8 pt widens the gap to
# 24.9 mm, which is what the 22.8 mm label needs.
_p = ax.get_position(); ax.set_position([_p.x0 + 0.14 * _p.width, _p.y0, 0.82 * _p.width, _p.height])
ax = fig.add_subplot(gs[1, 5:]); letter(ax, "e", x=-0.45); title(ax, "Observed ΔD exceeds the\nSST-calibrated expectation")
cal = smf.ols("D_index ~ SST", data=pre).fit()
pr = cal.get_prediction(hw).summary_frame(); pr.index = hw.index
exp_mean = pr["mean"].mean() - cal.fittedvalues.mean()
exp_boot, obs_boot = [], []
for _ in range(3000):
    bp_ = pre.sample(len(pre), replace=True, random_state=RNG.integers(1e9)); bh_ = hw.sample(len(hw), replace=True, random_state=RNG.integers(1e9))
    c_ = np.polyfit(bp_.SST, bp_.D_index, 1)
    exp_boot.append(c_[0] * (bh_.SST.mean() - bp_.SST.mean()))
    obs_boot.append(bh_.D_index.mean() - bp_.D_index.mean())
obs_mean = hw.D_index.mean() - pre.D_index.mean()
bars = pd.DataFrame({"quantity": ["Expected from\nSST calibration", "Observed"], "dD": [exp_mean, obs_mean],
                     "lo": [np.percentile(exp_boot, 2.5), np.percentile(obs_boot, 2.5)], "hi": [np.percentile(exp_boot, 97.5), np.percentile(obs_boot, 97.5)]})
ax.bar([0, 1], bars.dD, color=[GREY, A_C], width=0.6)
ax.errorbar([0, 1], bars.dD, yerr=[bars.dD - bars.lo, bars.hi - bars.dD], fmt="none", ecolor="k", elinewidth=0.6, capsize=1.5)
ax.set_xticks([0, 1], bars.quantity, fontsize=5.75)
ax.set_ylabel("Mean ΔD during 2014–2015")
ax.text(0.03, 0.97, f"ΔSST = +{hw.SST.mean() - pre.SST.mean():.1f} °C\npre-event slope {cal.params.SST:.3f} D/°C",
        transform=ax.transAxes, ha="left", va="top", fontsize=5.75)
ax.set_ylim(0, 0.75)
save(fig, "Fig2_test", {"a_b_viral_series": v, "d_forest": forest.drop(columns="col"), "e_expected_vs_observed": bars})

# ================================================================== Fig. 3 VALIDATE
import cartopy.crs as ccrs
import cartopy.feature as cfeature
fig = plt.figure(figsize=(180 * MM, 124 * MM))
gs = fig.add_gridspec(2, 3, left=0.095, right=0.94, bottom=0.16, top=0.93, wspace=0.85, hspace=0.45, height_ratios=[1, 0.78], width_ratios=[1, 1, 1])
ax = fig.add_subplot(gs[0, :2], projection=ccrs.Robinson(central_longitude=0))
ax.set_global(); ax.add_feature(cfeature.LAND, facecolor="#E6E6E6", lw=0); ax.coastlines(lw=0.25, color="#9E9E9E"); ax.spines["geo"].set_linewidth(0.4)
for lat in (-35, 35):
    ax.plot(np.linspace(-180, 180, 100), np.full(100, lat), transform=ccrs.PlateCarree(), color=GREY, lw=0.4, ls=":")
sc = ax.scatter(stn.lon, stn.lat, c=stn.D, cmap=D_CMAP, norm=NORM, s=11, edgecolor="k", lw=0.25, transform=ccrs.PlateCarree(), zorder=3)
ax.text(-0.06, 1.04, "a", transform=ax.transAxes, fontsize=9.2, fontweight="bold")
ax.set_title(f"Tara Oceans: station-mean D ({len(stn)} stations, {len(tara)} metagenomes)", loc="left", fontweight="bold")
cax = ax.inset_axes([1.01, 0.12, 0.015, 0.76]); cb = fig.colorbar(sc, cax=cax); cb.set_label("D-index"); cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=5.75)

ax = fig.add_subplot(gs[0, 2]); letter(ax, "b"); title(ax, "Global temperature\ncalibration")
sc2 = ax.scatter(stn["T"], stn.D, c=stn.lat.abs(), cmap="viridis_r", vmin=0, vmax=80, s=7, lw=0.2, edgecolor="k")
EPS = 1e-3
tara["L"] = np.log((tara.I_relative + EPS) / (tara.A_relative + EPS)); tara["femto"] = (tara.size_fraction == "Femto-size").astype(int)
lt = smf.ols("L ~ temperature", data=tara).fit()
xs = np.linspace(-2, 31, 100)
ax.plot(xs, np.tanh(lt.predict(pd.DataFrame({"temperature": xs})) / 2), color="k", lw=1)
ax.axhline(0, color=LGREY, lw=0.5, ls="--", zorder=0)
r_ = H("Tara station-mean D vs temperature Pearson r")
ax.text(0.03, 0.97, f"$r$ = {r_.value:.2f}, {fmt_p(r_.P)}", transform=ax.transAxes, va="top", fontsize=6.325,
        bbox=dict(facecolor="white", alpha=1, edgecolor="none", pad=1))
ax.set_xlabel("Temperature (°C)"); ax.set_ylabel("Station-mean D"); ax.set_ylim(-1.05, 1.05)
cax = ax.inset_axes([0.62, 0.20, 0.33, 0.035]); cb = fig.colorbar(sc2, cax=cax, orientation="horizontal")
cb.set_label("|Latitude| (°)", fontsize=5.75, labelpad=1); cb.ax.tick_params(labelsize=6.325, length=1.5, pad=1); cb.outline.set_linewidth(0.3)

BAS = {"Arctic Ocean": "#56B4E9", "Southern Ocean": "#000000", "Atlantic Ocean": "#E69F00", "Pacific Ocean": "#009E73",
       "Indian Ocean": "#CC79A7", "Mediterranean Sea": "#D55E00", "Red Sea": "#F0E442"}
oof14 = pd.read_csv(A29 / "nestedcv_oof_7pred.csv")
_t = tara.reset_index(drop=True)
assert len(_t) == len(oof14) and (_t.sample_id.values == oof14.sample_id.values).all()
vv = pd.DataFrame({"sample_id": oof14.sample_id, "station_id": oof14.station_id, "ocean_basin": _t.ocean_basin, "temperature": _t.temperature,
                   "observed_D": oof14.observed_D, "predicted_D": oof14.oof_D})
vv["residual_D"] = vv.observed_D - vv.predicted_D
ax = fig.add_subplot(gs[1, 0]); letter(ax, "c"); title(ax, "Station-blocked out-of-fold prediction")
for bname, gb in vv.groupby("ocean_basin"):
    ax.scatter(gb.observed_D, gb.predicted_D, s=5, color=BAS[bname], lw=0.2, edgecolor="k", label=bname.replace(" Ocean", ""))
ax.plot([-1, 1], [-1, 1], color=GREY, ls="--", lw=0.6)
ens14 = pd.read_csv(A29 / "ensemble_validation.csv").iloc[0]
pr2 = pd.Series({"value": ens14.R2}); rm = pd.Series({"value": ens14.RMSE})
ax.text(0.03, 0.97, f"$R^2$ = {pr2.value:.2f}, RMSE = {rm.value:.2f}\nn = 196 (83 stations)", transform=ax.transAxes, va="top", fontsize=6.325)
ax.set_xlabel("Observed D"); ax.set_ylabel("Predicted D"); ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=5.2, frameon=False, handletextpad=0.25, columnspacing=0.7, labelspacing=0.18)

ax = fig.add_subplot(gs[1, 1]); letter(ax, "d"); title(ax, "Geographic validation\n(pooled basin holdout)")
_p = ax.get_position(); ax.set_position([_p.x0, 0.09, _p.width, _p.y1 - 0.09])   # stretch into the space the panel-c legend does not use
order = ["temperature-only D", "ridge log-ratio (7 predictors)", "random forest log-ratio (7 predictors)", "random forest D (7 predictors)"]
labs = ["Temperature +\nfraction (D)", "Ridge, 7 env\n+ fraction",
        "Random forest, 7 env\n+ fraction (log-ratio)", "Random forest, 7 env\n+ fraction (D)"]
yy = np.arange(len(order))[::-1]
rows_d = []
for y_, mname, lab in zip(yy, order, labs):
    s = cv[(cv.scheme == "5-fold station-blocked") & (cv.model == mname)].R2
    b_ = cv[(cv.scheme == "leave-one-basin-out") & (cv.model == mname)].R2.iloc[0]
    ax.plot([s.min(), s.max()], [y_ + 0.12] * 2, color=I_C, lw=0.8); ax.plot(s.mean(), y_ + 0.12, "o", ms=3, color=I_C, lw=0)
    ax.plot(b_, y_ - 0.12, "D", ms=3, color=A_C, lw=0)
    rows_d.append(dict(model=lab, station_blocked_R2_mean=s.mean(), station_blocked_R2_min=s.min(), station_blocked_R2_max=s.max(), leave_one_basin_out_R2=b_))
ax.plot(pr2.value, len(order), "o", color="k", ms=3)
ax.set_yticks(list(yy) + [len(order)], labs + ["Ensemble\n(14 models)"], fontsize=5.0)
_r2min = min([0.0] + [r["station_blocked_R2_min"] for r in rows_d])
ax.axvline(0, color="k", lw=0.4, ls="--", zorder=0)
ax.set_xlim(min(_r2min - 0.03, -0.03), 0.75); ax.set_ylim(-0.5, len(order) + 1.6)
ax.set_xlabel("Out-of-sample R² for D")
ax.legend(handles=[Line2D([], [], color=I_C, marker="o", ms=3, lw=0.8, label="Station-blocked (20 partitions)"),
                   Line2D([], [], color=A_C, marker="D", ms=3, lw=0, label="Leave one basin out (pooled)")],
          loc="upper left", fontsize=6.67, bbox_to_anchor=(-0.02, 1.0))

ax = fig.add_subplot(gs[1, 2]); letter(ax, "e"); title(ax, "Permutation importance")
_p = ax.get_position(); ax.set_position([_p.x0, 0.09, _p.width, _p.y1 - 0.09])   # stretch into the space the panel-c legend does not use
NAMES = {"temperature": "Temperature", "log_chl": "Chlorophyll a", "femto": "Size fraction", "log_CO3": "Carbonate", "salinity": "Salinity",
         "log_Fe": "Iron", "log_PO4": "Phosphate", "log_Si": "Silicate", "log_NO3": "Nitrate"}
ii = imp.sort_values("mean")
xerr = np.vstack([ii["mean"].to_numpy() - ii["ci_low"].to_numpy(), ii["ci_high"].to_numpy() - ii["mean"].to_numpy()])
ax.barh(range(len(ii)), ii["mean"], xerr=xerr, color=[I_C if f == "temperature" else "#BFBFBF" for f in ii.feature], edgecolor="k", lw=0.4,
        error_kw=dict(lw=0.5, capsize=1), height=0.65)
ax.set_yticks(range(len(ii)), [NAMES[f] for f in ii.feature], fontsize=7.13)
ax.axvline(0, color="k", lw=0.4, ls="--"); ax.set_xlabel("ΔR² after shuffling each predictor")
save(fig, "Fig3_validate", {"a_b_station_means": stn, "c_oof_predictions": vv, "d_model_transfer": pd.DataFrame(rows_d), "e_permutation_importance": imp})

# ================================================================== Fig. 4 INTERPRET
traits_file, micro_file = P.SUPP_DATA / "Supplementary_Data_7_GVMAG_functional_traits.csv", P.SUPP_DATA / "Supplementary_Data_8_GVMAG_microdiversity.csv"
fig = plt.figure(figsize=(180 * MM, 82 * MM))
gs = fig.add_gridspec(1, 3, left=0.15, right=0.97, bottom=0.19, top=0.85, wspace=0.34, width_ratios=[1.25, 0.9, 1.45])
ax = fig.add_subplot(gs[0]); letter(ax, "a", x=-0.28, y=1.02)
tl = pd.read_csv(traits_file)
tl = tl[tl.order.isin(["Imitervirales", "Algavirales"])].copy()
nI = int((tl.order == "Imitervirales").sum()); nA = int((tl.order == "Algavirales").sum())
trait_defs = [
    ("Heat-shock\nproteins", "any_HSP"),
    ("DNA repair", "pw_DNA_repair"),
    ("Type-I rhodopsin", "type_I_rhodopsin"),
    ("Glycolysis", "pw_glycolysis"),
    ("TCA cycle", "pw_TCA"),
]
tr_rows = []
for lab, col in trait_defs:
    kI = int(tl.loc[tl.order == "Imitervirales", col].sum())
    kA = int(tl.loc[tl.order == "Algavirales", col].sum())
    tr_rows.append((lab, kI, kA))
tr = pd.DataFrame(tr_rows, columns=["trait", "kI", "kA"])

def wilson(k, n, z=1.96):
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(c - h, 0), min(c + h, 1)

tr["pI"], tr["pA"] = tr.kI / nI, tr.kA / nA
tr["fisher_P"] = [stats.fisher_exact([[r.kI, nI-r.kI], [r.kA, nA-r.kA]])[1] for r in tr.itertuples()]
from statsmodels.stats.multitest import multipletests
tr["P_BH"] = multipletests(tr.fisher_P, method="fdr_bh")[1]
yy = np.arange(len(tr))[::-1]
for y_, r in zip(yy, tr.itertuples()):
    for k, n, col, dy in [(r.kA, nA, A_C, -0.14), (r.kI, nI, I_C, 0.14)]:
        lo, hi = wilson(k, n)
        ax.plot([lo*100, hi*100], [y_+dy]*2, color=col, lw=0.9)
        ax.plot(k/n*100, y_+dy, "o", ms=3.2, color=col, lw=0)
    # 作者要求：行末百分比注释去掉（数值在 source data 与图注中给出）
ax.set_yticks(yy, tr.trait); ax.set_xlim(0, 122); ax.set_xticks(range(0, 101, 20)); ax.spines["bottom"].set_bounds(0, 100)
ax.set_xlabel("Genomes encoding trait (%, 95% CI)")
ax.legend(handles=[Line2D([], [], color=I_C, marker="o", ms=3, label=f"Imitervirales ({nI} genomes)"),
                   Line2D([], [], color=A_C, marker="o", ms=3, label=f"Algavirales ({nA})")],
          loc="lower left", bbox_to_anchor=(-0.02, 1.04), ncol=2, fontsize=5.75)
_pmax = tr.P_BH.max(); _pe = int(np.floor(np.log10(_pmax)))
ax.text(1.0, 1.02, f"Fisher exact, all P_BH ≤ {_pmax / 10 ** _pe:.0f}e{_pe}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.67)
src4 = {"a_trait_prevalence": tr}

# Pathway richness per genome from the curated wide table
violin_panels_on()
ax = fig.add_subplot(gs[1]); letter(ax, "b", x=-0.14, y=1.06); pw = tl[["genome_id", "order", "n_pathways_present"]].dropna().copy()
valsA = pw.loc[pw.order == "Algavirales", "n_pathways_present"].to_numpy()
valsI = pw.loc[pw.order == "Imitervirales", "n_pathways_present"].to_numpy()
for xi, (vals, col) in enumerate([(valsA, A_C), (valsI, I_C)]):
    bp = ax.boxplot(vals, positions=[xi], widths=0.48, showfliers=False, patch_artist=True,
                    medianprops=dict(color="k", lw=0.8), boxprops=dict(facecolor="white", edgecolor=col, lw=0.8),
                    whiskerprops=dict(color=col, lw=0.6), capprops=dict(color=col, lw=0.6))
    ax.scatter(xi + RNG.normal(0, .05, len(vals)), vals, s=2.8, facecolor="none", edgecolor=col, lw=0.35, alpha=0.7)
up = stats.mannwhitneyu(valsI, valsA, alternative="two-sided")
ax.text(0.5, 0.97, f"median {np.median(valsA):.0f} vs {np.median(valsI):.0f}\n{fmt_p(up.pvalue)}", transform=ax.transAxes,
        ha="center", va="top", fontsize=5.75)
ax.set_xticks([0, 1], ["Algavirales", "Imitervirales"]); ax.set_ylabel("Number of pathway categories per genome")
src4["b_pathway_richness"] = pw

# Per-genome median nucleotide diversity to avoid sample-level pseudoreplication
violin_panels_on()
ax = fig.add_subplot(gs[2]); letter(ax, "c", x=-0.14, y=1.06); title(ax, "Population microdiversity")
md0 = pd.read_csv(micro_file)
md0 = md0[md0.order.isin(["Imitervirales", "Algavirales"])].copy()
md = (md0.groupby(["genome_id", "order", "size_fraction"], as_index=False)
      .mean_nucleotide_diversity_pi.median()
      .rename(columns={"mean_nucleotide_diversity_pi": "median_pi"}))
frac_order = ["0.02–0.2 µm", "0.2–1.2 µm", "1.2–80 µm"]
frac_labels = ["Viral", "Pico", "Micro/nano"]
for k, f in enumerate(frac_order):
    Avals = md[(md.size_fraction == f) & (md.order == "Algavirales")].median_pi.to_numpy()
    Ivals = md[(md.size_fraction == f) & (md.order == "Imitervirales")].median_pi.to_numpy()
    for dx, vals, col in [(-0.18, Avals, A_C), (0.18, Ivals, I_C)]:
        ax.boxplot(vals, positions=[k+dx], widths=0.28, showfliers=False, patch_artist=True,
                   medianprops=dict(color="k", lw=0.7), boxprops=dict(facecolor="white", edgecolor=col, lw=0.7),
                   whiskerprops=dict(color=col, lw=0.55), capprops=dict(color=col, lw=0.55))
        ax.scatter(k+dx+RNG.normal(0,0.035,len(vals)), vals, s=2.2, facecolor="none", edgecolor=col, lw=0.3, alpha=0.6)
    if len(Avals) and len(Ivals):
        pv = stats.mannwhitneyu(Ivals, Avals, alternative="two-sided").pvalue
        ytxt = max(np.nanpercentile(np.r_[Avals,Ivals], 95), 1e-6)
        frac_labels[k] = frac_labels[k] + "\n" + fmt_p(pv)
ax.set_xticks(range(3), frac_labels, fontsize=6.21); ax.set_ylabel("Nucleotide diversity (π) per genome")
_zz = np.linspace(-1, 1, 40); _ww = np.exp(-4 * _zz ** 2)
_vs = np.concatenate([_ww * 0.5, -_ww[::-1] * 0.5]); _vy = np.concatenate([_zz, _zz[::-1]])
def _violin_handle(c):
    return Polygon(np.c_[_vs, _vy], closed=True, facecolor="white", edgecolor=c, lw=0.7)
ax.legend(handles=[_violin_handle(A_C), _violin_handle(I_C)], labels=["Algavirales", "Imitervirales"],
          loc="upper left", bbox_to_anchor=(0.01, 0.99), fontsize=5.4, frameon=False,
          handlelength=1.1, handleheight=0.9, handletextpad=0.3, labelspacing=0.25)
src4["c_microdiversity"] = md
violin_panels_off()
save(fig, "Fig4_interpret", src4)

# ================================================================== Extended Data Fig. 1: index properties
fig, axs = plt.subplots(1, 4, figsize=(180 * MM, 58 * MM), gridspec_kw=dict(wspace=0.95, left=0.085, right=0.925, bottom=0.26, top=0.85, width_ratios=[1, 1, 1, 1]))
ax = axs[0]; letter(ax, "a", x=-0.3); title(ax, "D is a bounded log-ratio")
xs = np.linspace(-7, 7, 200); ax.plot(xs, np.tanh(xs / 2), color="k", lw=1)
ax.scatter(np.log(spot.I / spot.A), spot.D_index, s=3, color="#7570B3", lw=0, label="SPOT", zorder=3)
ax.scatter(np.log((tara.I_relative + EPS) / (tara.A_relative + EPS)), tara.D_index, s=3, color="#E6AB02", lw=0, label="Tara", zorder=2)
ax.set_xlabel("ln(I/A)"); ax.set_ylabel("D"); ax.legend(loc="upper left")
violin_panels_on()
ax = axs[1]; letter(ax, "b", x=-0.3); title(ax, "PolB is dominated by two orders")
for xi, f in enumerate(["viral", "pico", "micro"]):
    vals = spot[spot.frac == f].I_plus_A
    ax.boxplot(vals, positions=[xi], widths=0.5, showfliers=False, medianprops=dict(color="k"))
    ax.scatter(xi + RNG.normal(0, 0.06, len(vals)), vals, s=2, color=FR_C[f], lw=0)
ax.set_xticks([0, 1, 2], ["Viral", "Pico", "Micro/\nnano"]); ax.set_ylabel("(I + A) share of PolB"); ax.set_ylim(0.35, 1.02)
ax = axs[2]; letter(ax, "c", x=-0.3); title(ax, "Seasonal vs anomaly sensitivity")
coef = [("SPOT seasonal SST", "SPOT dD/dSST_clim, fraction-adjusted (HC3)"), ("SPOT SST anomaly", "SPOT dD/dSST_anomaly, fraction-adjusted (HC3)"),
        ("SPOT viral, pre-2014", "viral pre-event dD/dSST"), ("Tara, 13–24 °C", "Tara dD/dT within SPOT SST range 13-24 degC"),
        ("Tara, all", "Tara dD/dT (station-clustered)")]
yy = np.arange(len(coef))[::-1]
for y_, (lab, key) in zip(yy, coef):
    r = H(key); col = "#E6AB02" if "Tara" in lab else "#7570B3"
    ax.plot([r.ci_low, r.ci_high], [y_, y_], color=col, lw=0.8); ax.plot(r.value, y_, "o", ms=3, color=col, lw=0)
ax.set_yticks(yy, [c[0] for c in coef], fontsize=5.75); ax.axvline(0, color=LGREY, lw=0.5, ls="--"); ax.set_xlabel("dD/dT (D/°C, 95% CI)")
ax = axs[3]; letter(ax, "d", x=-0.3); zp = [("SPOT viral", "SST at D=0 (viral)"), ("SPOT pico", "SST at D=0 (pico)"), ("SPOT micro/nano", "SST at D=0 (micro)"),
      ("Tara <0.22 µm", "Tara temperature at D=0, femto"), ("Tara 0.22–3 µm", "Tara temperature at D=0, pico")]
yy = np.arange(len(zp))[::-1]
for y_, (lab, key) in zip(yy, zp):
    r = H(key); col = "#E6AB02" if "Tara" in lab else "#7570B3"
    ax.plot([r.ci_low, min(r.ci_high, 32)], [y_, y_], color=col, lw=0.8); ax.plot(r.value, y_, "o", ms=3, color=col, lw=0)
ax.set_yticks(yy, [z[0] for z in zp], fontsize=5.75); ax.set_xlabel("Temperature at D = 0\n(°C, 95% CI)"); ax.set_xlim(5, 32)
violin_panels_off()
save(fig, "ExtendedDataFig1_index_properties", {"c_slopes": hs[hs.statistic.isin([c[1] for c in coef])], "d_zero_points": hs[hs.statistic.isin([z[1] for z in zp])]})

# ================================================================== Extended Data Fig. 2: heatwave diagnostics
fig, axs = plt.subplots(1, 3, figsize=(180 * MM, 80 * MM), gridspec_kw=dict(wspace=0.28, left=0.09, right=0.955, bottom=0.22, top=0.86))
ax = axs[0]; letter(ax, "a", x=-0.25)
title(ax, "Year-by-year residual shift")
for per, g_, col in [("2009–2013", pre, GREY), ("2014–2015", hw, A_C)]:
    ax.scatter(g_.SST_anomaly, g_.resid, s=6, color=col, lw=0, label=per)
ax.axhline(0, color=LGREY, lw=0.5, ls="--"); ax.set_xlabel("SST anomaly (°C)"); ax.set_ylabel("Observed − expected D"); ax.legend(loc="upper left")
ma = H("viral 2009-2015 resid ~ SST_anomaly + MEI: SST_anomaly"); me = H("viral 2009-2015 resid ~ SST_anomaly + MEI: MEI")
ax.set_title(f"Residual D vs SST anomaly (viral)\nSST anomaly {ma.value:.2f} D/°C ({fmt_p(ma.P)}); MEI P = {me.P:.2f}", loc="left", fontweight="bold", fontsize=6.325)
ax = axs[1]; letter(ax, "b", x=-0.25); title(ax, "Year-label permutation test")
allyrs = sorted(v.year.unique()); import itertools
perm = [(combo, v[v.year.isin(combo)].resid.mean() - v[~v.year.isin(combo)].resid.mean()) for combo in itertools.combinations(allyrs, 2)]
pdf_ = pd.DataFrame(perm, columns=["years", "shift"]).sort_values("shift")
cols = [A_C if set(c) == {2014, 2015} else LGREY for c in pdf_.years]
ax.barh(range(len(pdf_)), pdf_["shift"], color=cols, height=0.7)
ax.set_yticks(range(len(pdf_)), [f"{a}+{b}" for a, b in pdf_.years], fontsize=5.98); ax.axvline(0, color="k", lw=0.4, ls="--")
ax.set_xlabel("Mean residual shift for that pair of years")
ax = axs[2]; letter(ax, "c", x=-0.25); title(ax, "Same-month differences")
hw2 = hw.assign(diff=hw.D_index - hw.month.map(pre.groupby("month").D_index.mean())).dropna(subset=["diff"])
ax.bar(range(len(hw2)), hw2["diff"], color=np.where(hw2["diff"] > 0, A_C, I_C))
ax.set_xticks(range(len(hw2)), hw2.date.dt.strftime("%b %y"), rotation=90, fontsize=6.325); ax.axhline(0, color="k", lw=0.4, ls="--")
ax.set_ylabel("D minus 2009–2013 mean of same month")
save(fig, "ExtendedDataFig2_heatwave", {"b_permutations": pdf_.assign(years=pdf_.years.astype(str)), "c_same_month": hw2[["date", "month", "D_index", "diff"]]})

# ================================================================== Extended Data Fig. 3: global robustness
fig, axs = plt.subplots(1, 4, figsize=(180 * MM, 72 * MM), gridspec_kw=dict(wspace=0.82, left=0.085, right=0.945, bottom=0.30, top=0.86))
ax = axs[0]; letter(ax, "a", x=-0.35); title(ax, "Chlorophyll gradient")
sc3 = ax.scatter(np.log10(stn.chl), stn.D, c=stn["T"], cmap="viridis", s=6, lw=0.2, edgecolor="k")
cax3 = ax.inset_axes([0.585, 0.835, 0.385, 0.024]); cb = fig.colorbar(sc3, cax=cax3, orientation="horizontal")
cb.set_label("T (°C)", fontsize=5.0, labelpad=1)
cb.ax.xaxis.label.set_bbox(dict(facecolor="white", alpha=1, edgecolor="none", pad=1))
cb.ax.tick_params(labelsize=5.98, length=1.5, pad=1); cb.outline.set_linewidth(0.3)
rc = H("Tara station-mean D vs chlorophyll a Pearson r")
ax.text(0.97, 0.97, f"$r$ = {rc.value:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=5.75)
ax.set_xlabel("log10 chlorophyll a (mg/m³)"); ax.set_ylabel("Station-mean D")
ax = axs[1]; letter(ax, "b", x=-0.35); title(ax, "Held-out basin error")
lob = cv[cv.scheme.str.startswith("LOBO")].copy(); lob["basin"] = lob.scheme.str.replace("LOBO held-out basin: ", "").str.replace(" Ocean", "")
for k, (mname, col) in enumerate([("temperature-only log-ratio", "k"), ("random forest D (7 predictors)", A_C)]):
    s = lob[lob.model == mname].reset_index(drop=True)
    ax.barh(np.arange(len(s)) + (0.2 if k == 0 else -0.2), s.RMSE, height=0.38, color=col, label="Temperature only (D)" if k == 0 else "Random forest")
ax.set_yticks(range(len(s)), [f"{b.replace(' Sea', '')}\n({int(n)})" for b, n in zip(s.basin, s.n)], fontsize=6.67); ax.set_xlabel("RMSE (D units)"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=5.4, frameon=False)
ax = axs[2]; letter(ax, "c", x=-0.35); title(ax, "Residual D by basin")
bas = [b for b in BAS if b in set(vv.ocean_basin)]
ax.boxplot([vv[vv.ocean_basin == b].residual_D for b in bas], widths=0.5, showfliers=False, medianprops=dict(color="k"))
ax.axhline(0, color=LGREY, lw=0.5, ls="--"); ax.set_xticks(range(1, len(bas) + 1), [b.replace(" Ocean", "").replace(" Sea", "") for b in bas], rotation=90, fontsize=5.75)
ax.set_ylabel("Observed − predicted D")
ax = axs[3]; letter(ax, "d", x=-0.35); title(ax, "Tara calibration\nat SPOT")
tr_ = pd.read_csv(AN / "inherited__spot_transfer_from_tara.csv")
for f in ["viral", "pico"]:
    g_ = tr_[tr_.frac == f]; ax.scatter(g_.D_tara, g_.D_index, s=5, marker=FR_M[f], facecolor="none", edgecolor=FR_C[f], lw=0.4, label=f.capitalize())
ax.plot([-1, 1], [-1, 1], color=GREY, ls="--", lw=0.5)
rt = H("SPOT observed D vs D predicted from Tara temperature calibration, r"); bt = H("transfer bias (predicted - observed D), mean")
ax.text(0.03, 0.97, f"$r$ = {rt.value:.2f}; bias {bt.value:+.2f}", transform=ax.transAxes, va="top", fontsize=5.75)
ax.set_xlabel("D predicted from Tara\n(SST only)"); ax.set_ylabel("Observed SPOT D"); ax.legend(loc="upper left", bbox_to_anchor=(0.01, 0.93), fontsize=6.325, handletextpad=0.4, borderpad=0.2, labelspacing=0.15, frameon=True, facecolor="white", edgecolor="none", framealpha=1)
save(fig, "ExtendedDataFig3_global_robustness", {"b_lobo": lob, "d_transfer": tr_})

# ================================================================== Fig. 5 PROJECT (: gridded maps + robustness)
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# A29 is set in the header
maps = np.load(P.CMIP6 / "cmip6_maps_2deg.npz")
D0m, D1m = np.nanmean(maps["D2000"], 0), np.nanmean(maps["D2090"], 0)
dstack = maps["D2090"] - maps["D2000"]
common = np.all(np.isfinite(dstack), 0)
agree = common & (np.all(dstack > 0, 0) | np.all(dstack < 0, 0))
dmean = np.where(common, np.nanmean(dstack, 0), np.nan)
lon2 = -179 + 2 * np.arange(180); lat2 = -89 + 2 * np.arange(90)
LON, LAT = np.meshgrid(lon2, lat2)
rates = pd.read_csv(A29 / "cmip6_front_rates_by_projection.csv")
band = pd.read_csv(A29 / "cmip6_agreement_by_band.csv")
traj = pd.read_csv(A29 / "cmip6_front_trajectories_by_projection.csv")
traj = traj[traj.projection == "final ensemble (7 predictors)"]
gst = pd.read_csv(A29 / "cmip6_statistics.csv")
def G(name):
    r = gst[gst.statistic == name]; assert len(r) == 1, name; return r.iloc[0]

fig = plt.figure(figsize=(180 * MM, 112 * MM))
gs = fig.add_gridspec(2, 3, left=0.08, right=0.96, bottom=0.10, top=0.95, wspace=0.3, hspace=0.42, height_ratios=[0.62, 1])
proj = ccrs.Robinson(central_longitude=-150)
map_panels = [("a", "Four-model mean D, 2000s", D0m, D_CMAP, NORM), ("b", "Four-model mean D, 2090s (SSP5-8.5)", D1m, D_CMAP, NORM)]
for k, (lt, ttl, arr, cmap, norm) in enumerate(map_panels):
    ax = fig.add_subplot(gs[0, k], projection=proj); ax.set_global()
    pm = ax.pcolormesh(np.r_[lon2 - 1, 181], np.r_[lat2 - 1, 91], arr, cmap=cmap, norm=norm, transform=ccrs.PlateCarree(), rasterized=True)
    ax.contour(LON, LAT, np.where(np.isfinite(arr), arr, np.nan), levels=[0], colors="k", linewidths=0.6, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#D9D9D9", lw=0, zorder=2); ax.coastlines(lw=0.2, color="#8C8C8C", zorder=3)
    ax.spines["geo"].set_linewidth(0.4)
    ax.text(-0.02, 1.04, lt, transform=ax.transAxes, fontsize=9.2, fontweight="bold"); ax.set_title(ttl, loc="left", fontweight="bold", x=0.06)
cax = ax.inset_axes([-0.95, -0.12, 0.8, 0.045]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
cb.set_label("D (black line, D = 0)", fontsize=5.75, labelpad=1); cb.ax.tick_params(labelsize=6.325, length=1.5, pad=1); cb.outline.set_linewidth(0.3)
ax = fig.add_subplot(gs[0, 2], projection=proj); ax.set_global()
DCM = LinearSegmentedColormap.from_list("dD", ["#7570B3", "#F7F7F7", "#1B7837"])
pm2 = ax.pcolormesh(np.r_[lon2 - 1, 181], np.r_[lat2 - 1, 91], dmean, cmap=DCM, norm=TwoSlopeNorm(vmin=-0.4, vcenter=0, vmax=0.4),
                    transform=ccrs.PlateCarree(), rasterized=True)
st = (~agree) & common
ax.scatter(LON[st][::2], LAT[st][::2], s=0.25, color="#4D4D4D", lw=0, transform=ccrs.PlateCarree(), zorder=4)
ax.add_feature(cfeature.LAND, facecolor="#D9D9D9", lw=0, zorder=2); ax.coastlines(lw=0.2, color="#8C8C8C", zorder=3); ax.spines["geo"].set_linewidth(0.4)
ax.text(-0.02, 1.04, "c", transform=ax.transAxes, fontsize=9.2, fontweight="bold"); ax.set_title("ΔD, 2090s − 2000s", loc="left", fontweight="bold", x=0.06)
cax = ax.inset_axes([0.1, -0.12, 0.8, 0.045]); cb = fig.colorbar(pm2, cax=cax, orientation="horizontal", extend="both")
cb.set_label("ΔD (four-model mean); dots: models disagree on sign", fontsize=5.75, labelpad=1); cb.ax.tick_params(labelsize=6.325, length=1.5, pad=1); cb.outline.set_linewidth(0.3)

ax = fig.add_subplot(gs[1, 0]); letter(ax, "d", x=-0.2); title(ax, "D = 0 fronts by Earth-system model")
for esm, ge in traj.groupby("ESM"):
    for hem, ls, mk in [("North", "-", "o"), ("South", "--", "s")]:
        g_ = ge[ge.hemisphere == hem]
        ax.plot(g_.year, g_.front, ls=ls, marker=mk, ms=1.8, lw=0.8, color=ESM_C[esm], label=esm if hem == "North" else None)
ax.set_xlabel("Decade"); ax.set_ylabel("Front latitude (° N solid, ° S dashed)")
ax.legend(loc="upper left", fontsize=6.325, handlelength=1.4, labelspacing=0.2, borderaxespad=0.2)

ax = fig.add_subplot(gs[1, 1]); letter(ax, "e", x=-0.42); title(ax, "Poleward rate across projection designs")
PROJ = [("final ensemble (7 predictors)", "Final ensemble"),
        ("final ensemble, AOA cells only", "Final ensemble,\nsupport domain only"),
        ("random forest, 7 predictors", "Random forest"),
        ("temperature-only", "Temperature only")]
yy = np.arange(len(PROJ))[::-1]
for y_, (key, lab) in zip(yy, PROJ):
    for hem, col, dy in [("North", I_C, 0.14), ("South", A_C, -0.14)]:
        g_ = rates[(rates.projection == key) & (rates.hemisphere == hem)]
        ax.plot([g_.km_per_decade.min(), g_.km_per_decade.max()], [y_ + dy] * 2, color=col, lw=0.7)
        ax.scatter(g_.km_per_decade, np.full(len(g_), y_ + dy), s=5, color=[ESM_C[e] for e in g_.ESM], lw=0, zorder=3)
        ax.plot(g_.km_per_decade.mean(), y_ + dy, "|", color=col, ms=6, mew=1.4, zorder=4)
ax.axvline(0, color="k", lw=0.4, ls="--")
ax.set_yticks(yy, [p[1] for p in PROJ], fontsize=6.67)
for _t in ax.get_yticklabels():
    _t.set_bbox(dict(facecolor="white", edgecolor="none", pad=0.6))
ax.set_xlim(-5, 150); ax.set_xticks(range(0, 121, 20)); ax.spines["bottom"].set_bounds(-5, 120); ax.set_xlabel("Rate (km per decade)")
ax.legend(handles=[Line2D([], [], color=I_C, marker="|", ms=6, mew=1.4, lw=0.7, label="North"),
                   Line2D([], [], color=A_C, marker="|", ms=6, mew=1.4, lw=0.7, label="South")], loc="upper right", fontsize=6.325)

# [panel-width FIG5-E] e: 20 % off the left, 6 % off the right.  Measured:
# e's y tick label "support domain only" ran from x=136 to x=197 pt while
# panel d ends at x=162 pt, i.e. 9.0 mm inside d.
_p = ax.get_position(); ax.set_position([_p.x0 + 0.20 * _p.width, _p.y0, 0.74 * _p.width, _p.height])
ax = fig.add_subplot(gs[1, 2]); letter(ax, "f", x=-0.2); title(ax, "Model agreement by latitude")
b_ = band[band.hemisphere == "both"].reset_index(drop=True)
labs_b = {"|lat| < 30": "< 30°", "30-45": "30–45°", "45-60": "45–60°", "> 60": "> 60°"}
ax.bar(range(len(b_)), b_.share_sign_agreement * 100, color=["#4D4D4D", "#1A1A1A", "#8C8C8C", "#D0D0D0"], width=0.65, lw=0)
for i, r in b_.iterrows():
    ax.text(i, r.share_sign_agreement * 100 + 2, f"{r.share_sign_agreement * 100:.0f}%", ha="center", fontsize=5.75)
ax.set_xticks(range(len(b_)), [labs_b[x] for x in b_.band]); ax.set_ylim(0, 105)
ax.set_xlabel("Absolute latitude band"); ax.set_ylabel("Ocean area where all four\nmodels agree on sign of ΔD (%)")
ax.text(1, 8, "fronts", ha="center", color="white", fontsize=5.75, fontweight="bold")
# [panel-width FIG5-F] f: 5 % off each side, small, to keep e and f balanced.
_p = ax.get_position(); ax.set_position([_p.x0 + 0.05 * _p.width, _p.y0, 0.90 * _p.width, _p.height])
save(fig, "Fig5_project", {"a_b_c_maps_2deg_4model_mean": pd.DataFrame({"lat": LAT.ravel(), "lon": LON.ravel(), "D_2000s": D0m.ravel(), "D_2090s": D1m.ravel(),
                                                                        "dD": dmean.ravel(), "all_models_agree": agree.ravel()}).dropna(subset=["D_2000s"]),
                           "d_front_trajectories": traj, "e_rates_by_projection": rates, "f_agreement_by_band": band})

# ================================================================== Extended Data Fig. 5: projection robustness
aoa_maps = np.load(P.CMIP6 / "cmip6_aoa_2deg.npz")
units = pd.read_csv(A29 / "cmip6_vs_tara_predictor_units.csv")
aoa_tab = pd.read_csv(A29 / "cmip6_area_of_applicability.csv")
fig = plt.figure(figsize=(180 * MM, 98 * MM))
gs = fig.add_gridspec(2, 3, left=0.09, right=0.96, bottom=0.11, top=0.94, wspace=0.35, hspace=0.5, height_ratios=[0.8, 1])
for k, (key, ttl) in enumerate([("aoa2000", "Calibration support, 2000s"), ("aoa2090", "Calibration support, 2090s")]):
    ax = fig.add_subplot(gs[0, k], projection=proj); ax.set_global()
    pm = ax.pcolormesh(np.r_[lon2 - 1, 181], np.r_[lat2 - 1, 91], aoa_maps[key] * 100, cmap="Greys_r", vmin=0, vmax=100, transform=ccrs.PlateCarree(), rasterized=True)
    ax.contour(LON, LAT, D0m if k == 0 else D1m, levels=[0], colors=[I_C], linewidths=0.6, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#F2E6D9", lw=0, zorder=2); ax.coastlines(lw=0.2, color="#8C8C8C", zorder=3); ax.spines["geo"].set_linewidth(0.4)
    ax.text(-0.02, 1.04, "ab"[k], transform=ax.transAxes, fontsize=9.2, fontweight="bold"); ax.set_title(ttl, loc="left", fontweight="bold", x=0.06)
cax = ax.inset_axes([-0.9, -0.12, 0.8, 0.045]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
cb.set_label("Grid cells within area of applicability (%, four-model mean); blue, D = 0", fontsize=6.67, labelpad=1); cb.ax.tick_params(labelsize=6.325, length=1.5, pad=1)
ax = fig.add_subplot(gs[0, 2]); letter(ax, "c", x=-0.22, y=1.02); title(ax, "Share of ocean within calibration domain")
spec = "unweighted, 7 predictors (no CO3)"
cats = [("share_all", "All ocean"), ("share_near_front", "Near fronts\n(|D| < 0.15)"), ("share_polar", "|lat| ≥ 60°")]
for i, (col, lab) in enumerate(cats):
    for dx, dec, fc in [(-0.18, "2000-2009", "white"), (0.18, "2090-2099", GREY)]:
        s_ = aoa_tab[(aoa_tab.spec == spec) & (aoa_tab.decade == dec)]
        ax.bar(i + dx, s_[col].mean() * 100, width=0.34, color=fc, edgecolor="k", lw=0.5, label=dec[:4] + "s" if i == 0 else None)
        ax.scatter(np.full(4, i + dx), s_[col] * 100, s=5, c=[ESM_C[e] for e in s_.ESM], lw=0.2, edgecolor="k", zorder=3)
ax.set_xticks(range(3), [c[1] for c in cats], fontsize=5.75); ax.set_ylim(0, 118); ax.set_yticks(range(0, 101, 20)); ax.spines["left"].set_bounds(0, 100); ax.set_ylabel("Area-weighted share (%)"); ax.legend(loc="upper right", bbox_to_anchor=(0.99, 0.99), ncol=1, fontsize=max(6.67 - 0.5, 5.2), frameon=False)
ax = fig.add_subplot(gs[1, 0]); letter(ax, "d", x=-0.2); title(ax, "CMIP6 inputs vs Tara observations")
u_ = units.groupby("variable")[["median_ratio_cmip6_over_tara", "spearman_rho"]].median().reindex(["temperature", "salinity", "chlorophyll_a", "PO4", "NO3", "Fe", "Si", "CO3"])
VN = {"temperature": "Temperature", "salinity": "Salinity", "chlorophyll_a": "Chlorophyll a", "PO4": "Phosphate", "NO3": "Nitrate", "Fe": "Iron", "Si": "Silicate", "CO3": "'CO3' (Tara: pigment)"}
ax.scatter(u_.median_ratio_cmip6_over_tara, u_.spearman_rho, s=10, color=[A_C if v in ("CO3", "Si") else "k" for v in u_.index], zorder=3)
OFF = {"temperature": (1.15, 0.03), "salinity": (1.15, -0.04), "chlorophyll_a": (1.15, -0.03), "PO4": (1.15, 0.02)}
for v, r in u_.iterrows():
    fx, dy = OFF.get(v, (1.15, 0))
    ax.text(r.median_ratio_cmip6_over_tara * fx, r.spearman_rho + dy, VN[v], fontsize=6.67, va="center")
ax.set_xscale("log"); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); ax.axvline(1, color=LGREY, lw=0.6, ls="--"); ax.set_xlim(0.5, 60); ax.set_ylim(0, 1)
ax.set_xlabel("Median ratio, CMIP6 / Tara (co-located)"); ax.set_ylabel("Spearman ρ across 83 stations")
ax = fig.add_subplot(gs[1, 1]); letter(ax, "e", x=-0.3); title(ax, "Temperature-only projection: fronts")
tt = pd.read_csv(A29 / "cmip6_front_trajectories_by_projection.csv"); tt = tt[tt.projection == "temperature-only"]
for esm, ge in tt.groupby("ESM"):
    for hem, ls, mk in [("North", "-", "o"), ("South", "--", "s")]:
        g_ = ge[ge.hemisphere == hem]; ax.plot(g_.year, g_.front, ls=ls, marker=mk, ms=1.8, lw=0.8, color=ESM_C[esm])
ax.set_xlabel("Decade"); ax.set_ylabel("Front latitude (° N solid, ° S dashed)")
ax = fig.add_subplot(gs[1, 2]); letter(ax, "f", x=-0.3); title(ax, "Latitudinal profiles by model")
for esm, ge in latyr.groupby("ESM"):
    for yr, ls in [(ge.year.min(), ":"), (ge.year.max(), "-")]:
        g_ = ge[ge.year == yr]; ax.plot(g_.mean_D, g_.latitude_mid, color=ESM_C[esm], lw=0.7, ls=ls)
ax.axvline(0, color=LGREY, lw=0.5, ls="--"); ax.set_xlim(-1, 1); ax.set_ylim(-80, 80); ax.set_xlabel("Zonal-mean D"); ax.set_ylabel("Latitude (°)")
ax.legend(handles=[Line2D([], [], color="k", ls=":", lw=0.7, label="2000s"), Line2D([], [], color="k", lw=0.7, label="2090s")] +
          [Line2D([], [], color=ESM_C[e], lw=1, label=e) for e in ESM_C], loc="center left", fontsize=5.98, labelspacing=0.15)
save(fig, "ExtendedDataFig5_projection_robustness", {"a_b_aoa_share_2deg": pd.DataFrame({"lat": LAT.ravel(), "lon": LON.ravel(), "aoa_2000s": aoa_maps["aoa2000"].ravel(), "aoa_2090s": aoa_maps["aoa2090"].ravel()}).dropna(),
                                                     "c_aoa_shares": aoa_tab, "d_units": units, "e_temperature_only_fronts": tt, "f_profiles": latyr[latyr.year.isin([2004.5, 2094.5])]})

# ================================================================== Extended Data Fig. 4: marker recovery (rebuilt in house style)
# marker-recovery tables come from MR (header)
sd7 = pd.read_csv(MR / "Supplementary_Data_7_GVMAG_traits_marker_recovery.csv")
sd7 = sd7[sd7.order.isin(["Imitervirales", "Algavirales"])]
mh = pd.read_csv(MR / "marker_recovery_adjusted_genome_comparisons.csv")
fig = plt.figure(figsize=(180 * MM, 78 * MM))
gs = fig.add_gridspec(1, 3, left=0.07, right=0.965, bottom=0.2, top=0.86, wspace=0.55, width_ratios=[0.8, 1.2, 1])
ax = fig.add_subplot(gs[0]); letter(ax, "a", x=-0.35); title(ax, "Core-marker recovery")
for xi, (o, col) in enumerate([("Imitervirales", I_C), ("Algavirales", A_C)]):
    vals = sd7.loc[sd7.order == o, "marker_recovery_score"]
    cnt = vals.value_counts().sort_index()
    for sc_, n_ in cnt.items():
        ax.scatter(xi, sc_, s=n_ * 1.2, color=col, lw=0, alpha=1, zorder=3)
        ax.text(xi + 0.22, sc_, str(n_), fontsize=6.44, va="center")
ax.set_xticks([0, 1], ["Imitervirales", "Algavirales"]); ax.set_xlim(-0.5, 1.6); ax.set_ylim(2.5, 6.5); ax.set_yticks([3, 4, 5, 6])
ax.set_ylabel("Six-marker recovery score")
r0 = mh.iloc[0]
ax.text(0.5, 0.02, f"median 5 vs 5; P = {r0.recovery_score_p:.2f}", transform=ax.transAxes, ha="center", fontsize=6.67)
ax = fig.add_subplot(gs[1]); letter(ax, "b", x=-0.3); title(ax, "Contrasts after stratifying by recovery")
mhf = mh[np.isfinite(mh.MH_adjusted_odds_ratio)].reset_index(drop=True)
yy = np.arange(len(mhf))[::-1]
for y_, r in zip(yy, mhf.itertuples()):
    ax.plot([r.MH_95CI_low, r.MH_95CI_high], [y_, y_], color="k", lw=0.8); ax.plot(r.MH_adjusted_odds_ratio, y_, "o", ms=3.5, color="#333333", lw=0)
ax.set_xscale("log"); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); ax.axvline(1, color=GREY, lw=0.5, ls="--"); ax.set_yticks(yy, mhf.trait, fontsize=6.67)
tca = mh[~np.isfinite(mh.MH_adjusted_odds_ratio)]
if len(tca):
    t_ = tca.iloc[0]; ax.text(0.98, 0.02, f"TCA cycle: {int(t_.Imitervirales_positive)}/{int(t_.Imitervirales_total)} vs {int(t_.Algavirales_positive)}/{int(t_.Algavirales_total)}\n(OR not finite; CMH {fmt_p(t_.CMH_p)})",
                              transform=ax.transAxes, ha="right", va="bottom", fontsize=6.325)
ax.set_xlabel("Mantel–Haenszel odds ratio, Imitervirales vs Algavirales (95% CI)")
ax = fig.add_subplot(gs[2]); letter(ax, "c", x=-0.3); title(ax, "Pathway breadth by\nrecovery score")
for o, col, dx in [("Imitervirales", I_C, -0.06), ("Algavirales", A_C, 0.06)]:
    s_ = sd7[sd7.order == o].groupby("marker_recovery_score").n_pathways_present.agg(["mean", "sem"])
    ax.errorbar(s_.index + dx, s_["mean"], yerr=1.96 * s_["sem"].fillna(0), color=col, marker="o", ms=3, lw=0.8, elinewidth=0.6, capsize=1.5, label=o)
ax.set_xticks([3, 4, 5, 6]); ax.set_xlabel("Six-marker recovery score"); ax.set_ylabel("Number of pathway categories")
ax.legend(loc="upper right", bbox_to_anchor=(0.99, 0.99), fontsize=6.67, frameon=False)
ax.text(0.98, 0.40, f"adjusted difference +{r0.adjusted_pathway_difference:.2f}\n(95% CI {r0.adjusted_pathway_95CI_low:.2f}–{r0.adjusted_pathway_95CI_high:.2f})",
        transform=ax.transAxes, va="top", ha="right", fontsize=6.67)
save(fig, "ExtendedDataFig4_marker_recovery", {"a_scores": sd7[["genome_id", "order", "marker_recovery_score"]], "b_mantel_haenszel": mh,
                                               "c_pathways_by_score": sd7[["genome_id", "order", "marker_recovery_score", "n_pathways_present"]]})
