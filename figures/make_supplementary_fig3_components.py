"""Supplementary Fig. 12: projected Imitervirales and Algavirales abundances underlying D.

Rebuilt from the final seven-predictor final ensemble (analysis/cmip6_grid.npz), replacing
Fig. S14 of the earlier supplementary package, which used the superseded eight-field  model.
House style follows make_figures.py. Requires cartopy (as for Fig. 5 / Extended Data Fig. 4).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
AN = PKG / "processed_tables" / "model_and_statistics"
SRC = PKG / "outputs" / "source_data"
OUT = PKG / "outputs" / "supplementary"
for _d in (SRC, OUT):
    _d.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({
    # real Arial for text AND mathtext (was Liberation Sans + DejaVu Sans mixed)
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "mathtext.sf": "Arial", "mathtext.tt": "Arial", "mathtext.cal": "Arial:italic",
    "mathtext.default": "regular",
    "font.size": 6, "axes.titlesize": 6.5, "axes.labelsize": 6, "xtick.labelsize": 5.5,
    "ytick.labelsize": 5.5, "legend.fontsize": 5.5, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
    "lines.linewidth": 0.9, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False, "axes.titlepad": 3, "axes.labelpad": 2,
    "savefig.dpi": 600,
})
I_C, A_C, GREY, LGREY = "#0072B2", "#D55E00", "#7F7F7F", "#D9D9D9"
ESM_C = {"CESM2-WACCM": "#E69F00", "IPSL-CM6A-LR": "#56B4E9", "NorESM2-MM": "#009E73", "UKESM1-0-LL": "#CC79A7"}
ESMS = list(ESM_C)
DECS = [f"{y}-{y+9}" for y in range(2000, 2100, 10)]

z = np.load(AN / "cmip6_grid.npz", allow_pickle=True)
cols = {c: i for i, c in enumerate(map(str, z["columns"]))}
lon2 = -179 + 2 * np.arange(180); lat2 = -89 + 2 * np.arange(90)
LON, LAT = np.meshgrid(lon2, lat2)


def order_mean(a, order):
    """Fraction-averaged predicted relative abundance of one order (pico and femto)."""
    return 0.5 * (a[:, cols[f"{order}_pred_Pico"]] + a[:, cols[f"{order}_pred_Femto"]])


def grid2(lon, lat, val):
    """Mean of val on the 2-degree grid used by the main projection figures."""
    xi = np.clip(((lon + 180) / 2).astype(int), 0, 179)
    yi = np.clip(((lat + 90) / 2).astype(int), 0, 89)
    s = np.zeros((90, 180)); n = np.zeros((90, 180))
    np.add.at(s, (yi, xi), val); np.add.at(n, (yi, xi), 1)
    return np.where(n > 0, s / np.maximum(n, 1), np.nan)


maps, lat_prof, traj_rows = {}, [], []
for esm in ESMS:
    for dec in DECS:
        a = z[f"{esm}__{dec}"]
        lon, lat = a[:, cols["longitude"]], a[:, cols["latitude"]]
        for order in ("Imitervirales", "Algavirales"):
            v = order_mean(a, order)
            if dec in ("2000-2009", "2090-2099"):
                maps[(order, dec, esm)] = grid2(lon, lat, v)
            lb = np.floor(lat).astype(int)
            df = pd.DataFrame({"lat1": lb, "v": v}).groupby("lat1", as_index=False).v.mean()
            df["order"], df["ESM"], df["decade"] = order, esm, dec
            lat_prof.append(df)
            band = pd.cut(lat, np.arange(-90, 91, 15))
            g = pd.DataFrame({"band": band.astype(str), "v": v}).groupby("band", as_index=False).v.mean()
            g["order"], g["ESM"], g["decade"] = order, esm, dec
            traj_rows.append(g)
lat_prof = pd.concat(lat_prof, ignore_index=True)
traj = pd.concat(traj_rows, ignore_index=True)
traj["year"] = traj.decade.str.slice(0, 4).astype(int) + 5

M = {(o, d): np.nanmean(np.stack([maps[(o, d, e)] for e in ESMS]), 0) for o in ("Imitervirales", "Algavirales") for d in ("2000-2009", "2090-2099")}
DELTA = {o: M[(o, "2090-2099")] - M[(o, "2000-2009")] for o in ("Imitervirales", "Algavirales")}
stack = {o: np.stack([maps[(o, "2090-2099", e)] - maps[(o, "2000-2009", e)] for e in ESMS]) for o in DELTA}
agree = {o: np.all(np.isfinite(s), 0) & (np.all(s > 0, 0) | np.all(s < 0, 0)) for o, s in stack.items()}

vmax = max(np.nanmax(M[k]) for k in M)
ABU = LinearSegmentedColormap.from_list("abu", ["#FFFFFF", "#7FB7D9", "#08306B"])
DCM = LinearSegmentedColormap.from_list("d", ["#762A83", "#F7F7F7", "#1B7837"])
proj = ccrs.Robinson(central_longitude=-150)

fig = plt.figure(figsize=(180 * MM, 128 * MM))
gs = fig.add_gridspec(3, 3, left=0.07, right=0.985, bottom=0.085, top=0.96, wspace=0.28, hspace=0.55,
                      height_ratios=[0.58, 0.58, 1.0])
letters = iter("abcdefghi")


def lt(ax, s, ttl):
    ax.text(-0.02, 1.05, s, transform=ax.transAxes, fontsize=8, fontweight="bold")
    ax.set_title(ttl, loc="left", fontweight="bold", x=0.06, fontsize=6)


for row, order in enumerate(("Imitervirales", "Algavirales")):
    for k, (dec, lab) in enumerate([("2000-2009", "2000s"), ("2090-2099", "2090s (SSP5-8.5)")]):
        ax = fig.add_subplot(gs[row, k], projection=proj); ax.set_global()
        pm = ax.pcolormesh(np.r_[lon2 - 1, 181], np.r_[lat2 - 1, 91], M[(order, dec)], cmap=ABU, vmin=0, vmax=vmax,
                           transform=ccrs.PlateCarree(), rasterized=True)
        ax.add_feature(cfeature.LAND, facecolor=LGREY, lw=0, zorder=2); ax.coastlines(lw=0.2, color="#8C8C8C", zorder=3)
        ax.spines["geo"].set_linewidth(0.4)
        lt(ax, next(letters), f"{order}, {lab}")
        if row == 1 and k == 0:
            cax = ax.inset_axes([0.1, -0.13, 0.8, 0.05]); cb = fig.colorbar(pm, cax=cax, orientation="horizontal")
            cb.set_label("Predicted relative abundance", fontsize=5, labelpad=1)
            cb.ax.tick_params(labelsize=5.5, length=1.5, pad=1); cb.outline.set_linewidth(0.3)
    ax = fig.add_subplot(gs[row, 2], projection=proj); ax.set_global()
    lim = float(np.nanpercentile(np.abs(np.concatenate([DELTA[o][np.isfinite(DELTA[o])] for o in DELTA])), 99))
    pm2 = ax.pcolormesh(np.r_[lon2 - 1, 181], np.r_[lat2 - 1, 91], DELTA[order], cmap=DCM,
                        norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), transform=ccrs.PlateCarree(), rasterized=True)
    dis = np.isfinite(DELTA[order]) & ~agree[order]
    ax.scatter(LON[dis][::2], LAT[dis][::2], s=0.25, color="#4D4D4D", lw=0, transform=ccrs.PlateCarree(), zorder=4)
    ax.add_feature(cfeature.LAND, facecolor=LGREY, lw=0, zorder=2); ax.coastlines(lw=0.2, color="#8C8C8C", zorder=3)
    ax.spines["geo"].set_linewidth(0.4)
    lt(ax, next(letters), f"{order}, 2090s − 2000s")
    if row == 1:
        cax = ax.inset_axes([0.1, -0.13, 0.8, 0.05]); cb = fig.colorbar(pm2, cax=cax, orientation="horizontal", extend="both")
        cb.set_label("Change in relative abundance; dots: models disagree on sign", fontsize=5, labelpad=1)
        cb.ax.tick_params(labelsize=5.5, length=1.5, pad=1); cb.outline.set_linewidth(0.3)

ax = fig.add_subplot(gs[2, 0])
ax.text(-0.25, 1.05, next(letters), transform=ax.transAxes, fontsize=8, fontweight="bold")
ax.set_title("Latitudinal means", loc="left", fontweight="bold")
for order, col in (("Imitervirales", I_C), ("Algavirales", A_C)):
    for dec, ls in (("2000-2009", "-"), ("2090-2099", "--")):
        g = lat_prof[(lat_prof.order == order) & (lat_prof.decade == dec)].groupby("lat1", as_index=False).v.mean()
        ax.plot(g.v, g.lat1, ls=ls, color=col, lw=0.9)
ax.set_xlabel("Predicted relative abundance"); ax.set_ylabel("Latitude (°)")
ax.set_yticks(range(-90, 91, 30)); ax.axhline(0, color=LGREY, lw=0.4)
ax.legend(handles=[Line2D([], [], color=I_C, label="Imitervirales"), Line2D([], [], color=A_C, label="Algavirales"),
                   Line2D([], [], color="k", ls="-", label="2000s"), Line2D([], [], color="k", ls="--", label="2090s")],
          loc="lower right", fontsize=5.5, handlelength=1.6, labelspacing=0.25)

ax = fig.add_subplot(gs[2, 1])
ax.text(-0.25, 1.05, next(letters), transform=ax.transAxes, fontsize=8, fontweight="bold")
ax.set_title("Change by latitude band", loc="left", fontweight="bold")
bands = [b for b in sorted(traj.band.unique(), key=lambda s: float(s.split(",")[0][1:]))]
def _blab(b):
    lo, hi = (float(x) for x in b.strip("(]").split(", "))
    f = lambda v: f"{abs(int(v))}° {'N' if v > 0 else 'S'}" if v else "0°"
    return f"{f(lo)}–{f(hi)}"
blab = [_blab(b) for b in bands]
w = 0.38
for j, (order, col) in enumerate((("Imitervirales", I_C), ("Algavirales", A_C))):
    d0 = traj[(traj.order == order) & (traj.decade == "2000-2009")].groupby("band").v.mean()
    d1 = traj[(traj.order == order) & (traj.decade == "2090-2099")].groupby("band").v.mean()
    ax.barh(np.arange(len(bands)) + (w / 2 if j == 0 else -w / 2), [d1[b] - d0[b] for b in bands], height=w, color=col, label=order)
ax.axvline(0, color="k", lw=0.4); ax.set_yticks(range(len(bands)), blab, fontsize=5.8)
ax.set_xlabel("Change in relative abundance, 2090s − 2000s"); ax.set_ylabel("")
ax.legend(loc="lower right", fontsize=5.5)

ax = fig.add_subplot(gs[2, 2])
ax.text(-0.25, 1.05, next(letters), transform=ax.transAxes, fontsize=8, fontweight="bold")
ax.set_title("Decadal trajectories, 45–60° N", loc="left", fontweight="bold")
sel = [b for b in bands if b.strip("(]").split(", ")[0] in ("45", "45.0")][0]
for order, ls in (("Imitervirales", "-"), ("Algavirales", "--")):
    for esm in ESMS:
        g = traj[(traj.order == order) & (traj.ESM == esm) & (traj.band == sel)].sort_values("year")
        ax.plot(g.year, g.v, ls=ls, color=ESM_C[esm], lw=0.8, marker="o" if order == "Imitervirales" else "s", ms=1.6)
ax.set_xlabel("Decade (mid-year)"); ax.set_ylabel("Predicted relative abundance")
ax.legend(handles=[Line2D([], [], color=ESM_C[e], label=e) for e in ESMS] +
                  [Line2D([], [], color="k", ls="-", marker="o", ms=1.6, label="Imitervirales"),
                   Line2D([], [], color="k", ls="--", marker="s", ms=1.6, label="Algavirales")],
          loc="center left", bbox_to_anchor=(0.0, 0.42), fontsize=5.3, handlelength=1.6, labelspacing=0.2)

OUT.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT / "SupplementaryFig12_projected_components.pdf")
fig.savefig(OUT / "SupplementaryFig12_projected_components.svg")
fig.savefig(OUT / "SupplementaryFig12_projected_components.png", dpi=600)
plt.close(fig)
with pd.ExcelWriter(SRC / "SupplementaryFig3_source_data.xlsx") as xw:
    pd.DataFrame({"lat": LAT.ravel(), "lon": LON.ravel(),
                  **{f"{o}_{d}": M[(o, dd)].ravel() for o in ("Imitervirales", "Algavirales")
                     for d, dd in (("2000s", "2000-2009"), ("2090s", "2090-2099"))}}
                 ).dropna(subset=["Imitervirales_2000s"]).to_excel(xw, sheet_name="a_f_maps_2deg", index=False)
    lat_prof[lat_prof.decade.isin(["2000-2009", "2090-2099"])].to_excel(xw, sheet_name="g_latitudinal", index=False)
    traj.to_excel(xw, sheet_name="h_bands_decadal", index=False)
print("wrote Supplementary Fig. 12; abundance range 0 -", round(vmax, 3))
