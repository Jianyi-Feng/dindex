# -*- coding: utf-8 -*-
# Redrawn from the 1.15x font-size build; the original file is unchanged.
"""Extended Data Fig. 4 (): projection robustness, with the isotherm comparison added as panels g and h.

Panels a-f reproduce the / figure from make_figures.py; g and h are new and come from
analysis/climate_velocity_comparison.py. Built standalone so the other figures are not regenerated.
Requires cartopy.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
sys.path.insert(0, str(PKG / "tools"))
import paths as P  # noqa: E402  (repository layout is defined in tools/paths.py)

A29 = P.MODEL_STATS
DATA = P.RAW
SRC = P.outdir(P.OUT_SOURCE_DATA)
OUTFIG = P.outdir(P.OUT)
MM = 1 / 25.4
plt.rcParams.update({
    # real Arial for text AND mathtext (was Liberation Sans + DejaVu Sans mixed)
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial", "mathtext.tt": "Arial", "mathtext.cal": "Arial:italic",
    "mathtext.default": "regular",
    "font.size": 7.475, "axes.titlesize": 7.82, "axes.labelsize": 7.475, "xtick.labelsize": 7.13,
    "ytick.labelsize": 7.13, "legend.fontsize": 6.9, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.9, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "axes.titlepad": 3, "axes.labelpad": 2, "savefig.dpi": 600,
})
I_C, A_C, GREY, LGREY = "#0072B2", "#D55E00", "#7F7F7F", "#D9D9D9"
ESM_C = {"CESM2-WACCM": "#E69F00", "IPSL-CM6A-LR": "#56B4E9", "NorESM2-MM": "#009E73", "UKESM1-0-LL": "#CC79A7"}


def letter(ax, s, x=-0.2, y=1.08):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=9.2, fontweight="bold", va="bottom", ha="left")


def title(ax, s):
    ax.set_title(s, loc="left", fontweight="bold")


maps = np.load(P.CMIP6 / "cmip6_maps_2deg.npz")
D0m, D1m = np.nanmean(maps["D2000"], 0), np.nanmean(maps["D2090"], 0)
lon2 = -179 + 2 * np.arange(180); lat2 = -89 + 2 * np.arange(90)
LON, LAT = np.meshgrid(lon2, lat2)
aoa_maps = np.load(P.CMIP6 / "cmip6_aoa_2deg.npz")
units = pd.read_csv(A29 / "cmip6_vs_tara_predictor_units.csv")
aoa_tab = pd.read_csv(A29 / "cmip6_area_of_applicability.csv")
latyr = pd.read_csv(DATA / "CMIP6_Dindex_latitude_year.csv")
pair = pd.read_csv(A29 / "cmip6_front_vs_isotherm_by_model.csv")
iso = pd.read_csv(A29 / "cmip6_isotherm_displacement.csv")
head = pd.read_csv(A29 / "climate_velocity_headline_statistics.csv").set_index("statistic")
proj = ccrs.Robinson(central_longitude=-150)

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


fig = plt.figure(figsize=(180 * MM, 140 * MM))
gs = fig.add_gridspec(3, 3, left=0.09, right=0.96, bottom=0.08, top=0.95, wspace=0.35, hspace=0.45,
                      height_ratios=[0.8, 1, 1])

for k, (key, ttl) in enumerate([("aoa2000", "Environmental support, 2000s"), ("aoa2090", "Environmental support, 2090s")]):
    ax = fig.add_subplot(gs[0, k], projection=proj); ax.set_global()
    pm = ax.pcolormesh(np.r_[lon2 - 1, 181], np.r_[lat2 - 1, 91], aoa_maps[key] * 100, cmap="Greys_r", vmin=0, vmax=100,
                       transform=ccrs.PlateCarree(), rasterized=True)
    ax.contour(LON, LAT, D0m if k == 0 else D1m, levels=[0], colors=[I_C], linewidths=0.6, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor="#F2E6D9", lw=0, zorder=2); ax.coastlines(lw=0.2, color="#8C8C8C", zorder=3)
    ax.spines["geo"].set_linewidth(0.4)
    ax.text(-0.02, 1.04, "ab"[k], transform=ax.transAxes, fontsize=9.2, fontweight="bold")
    ax.set_title(ttl, loc="left", fontweight="bold", x=0.06)
cax = ax.inset_axes([-0.9, -0.12, 0.8, 0.045]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
cb.set_label("Grid cells within environmental-support domain (%, four-model mean); blue, D = 0", fontsize=7.245, labelpad=1)
cb.ax.tick_params(labelsize=6.9, length=1.5, pad=1)

ax = fig.add_subplot(gs[0, 2]); letter(ax, "c", x=-0.22, y=1.02); title(ax, "Share within support domain")
spec = "unweighted, 7 predictors (no CO3)"
cats = [("share_all", "All ocean"), ("share_near_front", "Near fronts\n(|D| < 0.15)"), ("share_polar", "|lat| ≥ 60°")]
for i, (col, lab) in enumerate(cats):
    for dx, dec, fc in [(-0.18, "2000-2009", "white"), (0.18, "2090-2099", GREY)]:
        s_ = aoa_tab[(aoa_tab.spec == spec) & (aoa_tab.decade == dec)]
        ax.bar(i + dx, s_[col].mean() * 100, width=0.34, color=fc, edgecolor="k", lw=0.5, label=dec[:4] + "s" if i == 0 else None)
        ax.scatter(np.full(4, i + dx), s_[col] * 100, s=5, c=[ESM_C[e] for e in s_.ESM], lw=0.2, edgecolor="k", zorder=3)
ax.set_xticks(range(3), [c[1] for c in cats], fontsize=7.13); ax.set_ylim(0, 118); ax.set_yticks(range(0, 101, 20))
ax.spines["left"].set_bounds(0, 100); ax.set_ylabel("Area-weighted share")
ax.legend(loc="upper center", ncol=2, fontsize=6.555, frameon=True, facecolor="white", edgecolor="none", framealpha=1)

ax = fig.add_subplot(gs[1, 0]); letter(ax, "d", x=-0.2); title(ax, "CMIP6 inputs vs Tara observations")
u_ = units.groupby("variable")[["median_ratio_cmip6_over_tara", "spearman_rho"]].median().reindex(
    ["temperature", "salinity", "chlorophyll_a", "PO4", "NO3", "Fe", "Si", "CO3"])
VN = {"temperature": "Temperature", "salinity": "Salinity", "chlorophyll_a": "Chlorophyll a", "PO4": "Phosphate",
      "NO3": "Nitrate", "Fe": "Iron", "Si": "Silicate", "CO3": "'CO3': mislabelled pigment"}
# distinguish the excluded field from a retained predictor that merely agrees weakly across datasets
for v, r in u_.iterrows():
    if v == "CO3":
        ax.scatter(r.median_ratio_cmip6_over_tara, r.spearman_rho, s=24, marker="X", color=A_C,
                   lw=0.3, edgecolor="k", zorder=4)
    else:
        ax.scatter(r.median_ratio_cmip6_over_tara, r.spearman_rho, s=10, marker="o",
                   color=A_C if v == "Si" else "k", zorder=3)
OFF = {"temperature": (1.15, 0.03, "left"), "salinity": (1.15, -0.04, "left"),
       "chlorophyll_a": (1.15, -0.03, "left"), "PO4": (1.15, 0.02, "left"),
       "CO3": (1.0, 0.075, "center"), "Si": (1.18, -0.045, "left")}
for v, r in u_.iterrows():
    fx, dy, ha = OFF.get(v, (1.15, 0, "left"))
    ax.text(r.median_ratio_cmip6_over_tara * fx, r.spearman_rho + dy, VN[v], fontsize=6.44, va="center", ha=ha)
ax.set_xscale("log"); ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}")); ax.axvline(1, color=LGREY, lw=0.6, ls="--"); ax.set_xlim(0.5, 60); ax.set_ylim(0, 1)
ax.set_xlabel("Median ratio, CMIP6 / Tara (co-located)"); ax.set_ylabel("Spearman ρ across 83 stations")
ax.legend(handles=[Line2D([], [], color="k", marker="o", ms=2.6, lw=0, label="Final predictor"),
                   Line2D([], [], color=A_C, marker="o", ms=2.6, lw=0, label="Retained, weak agreement"),
                   Line2D([], [], color=A_C, marker="X", ms=3.4, lw=0, label="Excluded from final models")],
          loc="upper right", bbox_to_anchor=(1.02, 0.62), fontsize=5.98, labelspacing=0.25, handletextpad=0.4)

ax = fig.add_subplot(gs[1, 1]); letter(ax, "e", x=-0.3); title(ax, "Temperature-only\nprojection: fronts")
tt = pd.read_csv(A29 / "cmip6_front_trajectories_by_projection.csv")
tt = tt[tt.projection == "temperature-only"]
for esm, ge in tt.groupby("ESM"):
    for hem, ls, mk in [("North", "-", "o"), ("South", "--", "s")]:
        g_ = ge[ge.hemisphere == hem]; ax.plot(g_.year, g_.front, ls=ls, marker=mk, ms=1.8, lw=0.8, color=ESM_C[esm])
ax.set_xlabel("Decade"); ax.set_ylabel("Front latitude (° N solid, ° S dashed)")

ax = fig.add_subplot(gs[1, 2]); letter(ax, "f", x=-0.3); ax.set_title("Latitudinal profiles by model", loc="left", fontweight="bold", pad=8)
for esm, ge in latyr.groupby("ESM"):
    for yr, ls in [(ge.year.min(), ":"), (ge.year.max(), "-")]:
        g_ = ge[ge.year == yr]; ax.plot(g_.mean_D, g_.latitude_mid, color=ESM_C[esm], lw=0.7, ls=ls)
ax.axvline(0, color=LGREY, lw=0.5, ls="--"); ax.set_xlim(-1, 1); ax.set_ylim(-80, 80)
ax.set_xlabel("Zonal-mean D"); ax.set_ylabel("Latitude (°)")
ax.legend(handles=[Line2D([], [], color="k", ls=":", lw=0.7, label="2000s"), Line2D([], [], color="k", lw=0.7, label="2090s")] +
          [Line2D([], [], color=ESM_C[e], lw=1, label=e) for e in ESM_C], loc="center left", fontsize=6.325, labelspacing=0.15)

# ------------------------------------------------------------------ g: front rate versus isotherm rate
ax = fig.add_subplot(gs[2, 0]); letter(ax, "g", x=-0.2); title(ax, "Front rate versus isotherm rate")
lim = (0, 95)
ax.plot(lim, lim, color=GREY, ls="--", lw=0.6, zorder=1)
for _, r in pair.iterrows():
    ax.scatter(r.isotherm, r.D_front, s=14, marker="o" if r.hemisphere == "North" else "s",
               facecolor=ESM_C[r.ESM], edgecolor="k", lw=0.3, zorder=3)
rr = head.loc["correlation of front rate with isotherm rate"]
rn = head.loc["front-to-isotherm rate ratio, North"]; rs = head.loc["front-to-isotherm rate ratio, South"]
ax.text(0.04, 0.97, f"r = {rr.value:.2f}, P = {rr.P:.3f}\nratio: north {rn.value:.2f} (P = {rn.P:.2f}),\nsouth {rs.value:.2f} (P = {rs.P:.3f})",
        transform=ax.transAxes, fontsize=6.67, va="top", linespacing=1.4)
ax.set_xlim(*lim); ax.set_ylim(*lim); ax.set_xticks(range(0, 91, 20)); ax.set_yticks(range(0, 91, 20))
ax.set_xlabel("Isotherm rate (km per decade)"); ax.set_ylabel("D = 0 front rate (km per decade)")
ax.legend(handles=[Line2D([], [], color="k", marker="o", ms=3, lw=0, label="North"),
                   Line2D([], [], color="k", marker="s", ms=3, lw=0, label="South"),
                   Line2D([], [], color=GREY, ls="--", lw=0.6, label="1:1")],
          loc="lower right", fontsize=6.67, labelspacing=0.2)

# ------------------------------------------------------------------ h: trajectories of front and isotherm
ax = fig.add_subplot(gs[2, 1:]); letter(ax, "h", x=-0.1); title(ax, "Poleward movement of the fronts and of the isotherms they began on")
ens = pd.read_csv(A29 / "cmip6_front_trajectories_by_projection.csv")
ens = ens[ens.projection == "final ensemble (7 predictors)"]
grid = np.load(P.CMIP6 / "cmip6_grid.npz", allow_pickle=True)
cols = {c: i for i, c in enumerate(map(str, grid["columns"]))}
DECS = [f"{y}-{y + 9}" for y in range(2000, 2100, 10)]
YEARS = [int(d[:4]) + 4.5 for d in DECS]


def crossing(profile, sign, lo=15, hi=65):
    p = profile[(sign * profile.index >= lo) & (sign * profile.index <= hi)]
    la, d = p.index.to_numpy() * sign, p.to_numpy()
    o = np.argsort(la); la, d = la[o], d[o]
    for i in range(len(la) - 1):
        if d[i] > 0 >= d[i + 1]:
            return la[i] + (la[i + 1] - la[i]) * d[i] / (d[i] - d[i + 1])
    return np.nan


iso_traj = []
for _, r in iso.iterrows():
    sign = 1 if r.hemisphere == "North" else -1
    for dec, yr in zip(DECS, YEARS):
        a = grid[f"{r.ESM}__{dec}"]
        lat = np.floor(a[:, cols["latitude"]]) + 0.5
        prof = pd.DataFrame({"lat": lat, "v": a[:, cols["Temperature"]]}).groupby("lat").v.mean()
        iso_traj.append(dict(ESM=r.ESM, hemisphere=r.hemisphere, year=yr,
                             lat=crossing(prof - r.isotherm_degC, sign)))
iso_traj = pd.DataFrame(iso_traj)

for hem, off in [("North", 0), ("South", 0)]:
    for esm in ESM_C:
        f_ = ens[(ens.ESM == esm) & (ens.hemisphere == hem)].sort_values("year")
        i_ = iso_traj[(iso_traj.ESM == esm) & (iso_traj.hemisphere == hem)].sort_values("year")
        ax.plot(f_.year, f_.front, color=ESM_C[esm], lw=0.9, marker="o" if hem == "North" else "s", ms=1.6)
        ax.plot(i_.year, i_.lat, color=ESM_C[esm], lw=0.7, ls=":", alpha=0.9)
ax.set_xlabel("Decade (mid-year)"); ax.set_ylabel("Absolute front latitude (°)")
ax.text(0.015, 0.97, "circles, northern front; squares, southern front; dotted, isotherm", transform=ax.transAxes,
        fontsize=6.67, va="top")
ax.legend(handles=[Line2D([], [], color="k", lw=0.9, label="D = 0 front"),
                   Line2D([], [], color="k", lw=0.7, ls=":", label="Isotherm")] +
          [Line2D([], [], color=ESM_C[e], lw=1, label=e) for e in ESM_C],
          loc="upper left", bbox_to_anchor=(0.0, 0.94), ncol=3, fontsize=6.21, labelspacing=0.2, columnspacing=0.8)

_fill_canvas_width(fig)
_tighten_row_gaps(fig)
_crop_canvas_height(fig)
fig.savefig(OUTFIG / "ExtendedDataFig5_projection_robustness.pdf")
fig.savefig(OUTFIG / "ExtendedDataFig5_projection_robustness.svg")
fig.savefig(OUTFIG / "ExtendedDataFig5_projection_robustness.png", dpi=600)
plt.close(fig)
with pd.ExcelWriter(SRC / "ExtendedDataFig5_projection_robustness_source_data.xlsx") as xw:
    pd.DataFrame({"lat": LAT.ravel(), "lon": LON.ravel(), "aoa_2000s": aoa_maps["aoa2000"].ravel(),
                  "aoa_2090s": aoa_maps["aoa2090"].ravel()}).dropna().to_excel(xw, sheet_name="a_b_aoa_share_2deg", index=False)
    aoa_tab.to_excel(xw, sheet_name="c_aoa_shares", index=False)
    units.to_excel(xw, sheet_name="d_units", index=False)
    tt.to_excel(xw, sheet_name="e_temperature_only_fronts", index=False)
    latyr[latyr.year.isin([2004.5, 2094.5])].to_excel(xw, sheet_name="f_profiles", index=False)
    pair.to_excel(xw, sheet_name="g_front_vs_isotherm", index=False)
    iso_traj.to_excel(xw, sheet_name="h_isotherm_trajectories", index=False)
print("wrote ExtendedDataFig5_projection_robustness (panels a-h)")
