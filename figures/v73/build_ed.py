"""Extended Data Figs 1-9 of the manuscript in the shared house style (see build_main.py). Drawn from the
deposited source data (processed_tables/figure_source_data). Run this script from anywhere:

    python figures/v73/build_ed.py          # all Extended Data figures
    python figures/v73/build_ed.py 7        # only Extended Data Fig. 7

Paths resolve from the repository root; the data folder can be moved with DINDEX_INPUTS and the output folder
with DINDEX_OUT. A few working files of the original author-side job are optional and each has an exact
fallback to the deposited tables (see figures/v73/README.md).
"""
import os, sys, json
import numpy as np, pandas as pd, matplotlib.pyplot as plt, matplotlib as mpl
import cartopy.crs as ccrs, cartopy.feature as cf
import gvstyle as G; G.use()
from build_main import (lab, grid, violin, MON, PKG, INPUTS, OUT, SD, SD59, SUP, RAW,
                        STUDY_JSON, DATA_EXTRA, rd, rd59)


def between_year_sensitivity():
    """Extended Data Fig. 2a input: P of the observed event shift under an assumed between-year s.d.
    With DINDEX_STUDY_JSON set this reads the author-local analysis JSON (.../strengthening_2026-11-23/
    results/02b_between_year_sensitivity.json); otherwise it rebuilds the same values from the deposited
    sheet ExtendedDataFig2_heatwave::a_between_year_sensitivity plus the s.d. of the five pre-event year
    means recorded in Supplementary Table 3."""
    if STUDY_JSON:
        with open(os.path.join(STUDY_JSON, '02b_between_year_sensitivity.json'), encoding='utf-8') as fh:
            return json.load(fh)
    d = rd('ExtendedDataFig2_heatwave', 'a_between_year_sensitivity')
    grid_ = {f'{r.between_year_sd:g}': {'P': float(r.P), 'null_sd': float(r.null_sd_of_shift)} for r in d.itertuples()}
    t3 = pd.read_excel(INPUTS / 'supplementary_tables' / 'Supplementary_Table_3_heatwave_tests.xlsx',
                       sheet_name='model_based_null')
    note = t3['note'].astype(str).str.extract(r'sd of the five pre-event year means ([\d.]+)')[0].dropna()
    if note.empty:
        raise RuntimeError('Supplementary Table 3 does not record the s.d. of the five pre-event year means')
    return {'grid': grid_, 'sd_of_pre_event_year_means': float(note.iloc[0])}

DK, MD, LT = (G.INK, '#6C7A90', '#C3CAD6') if G.THEME == 'vivid' else ('#1A1A1A', '#7A7A7A', '#C4C4C4')
EV = G.EVT if G.THEME == 'vivid' else DK; HO = G.HOST if G.THEME == 'vivid' else DK; HOL = '#9FD3BC' if G.THEME == 'vivid' else LT
def pint(ax, vals, los, his, labels, cols=None):
    for i, (v, l, h) in enumerate(zip(vals, los, his)):
        y = len(vals) - 1 - i; c = (cols[i] if cols else DK); ax.plot([l, h], [y, y], color=c, lw=0.9); ax.plot(v, y, 'o', color=c, ms=3.6)
    ax.set_yticks(range(len(vals))); ax.set_yticklabels(labels[::-1]); ax.set_ylim(-0.6, len(vals) - 0.4)

def ed1():
    F = 'ExtendedDataFig1_index_properties'; a = rd(F, 'a_D_vs_logratio'); b = rd(F, 'b_combined_share'); c = rd(F, 'c_slopes'); d = rd(F, 'd_zero_points')
    fig = G.figure(183, 92); gs = fig.add_gridspec(2, 2, left=0.135, right=0.985, top=0.94, bottom=0.095, wspace=0.62, hspace=0.62, width_ratios=[1.0, 1.0])
    ax = fig.add_subplot(gs[0]); x = np.linspace(-8, 8, 300); ax.plot(x, np.tanh(x / 2), color=LT, lw=2.2, zorder=1)
    for ds, kw in (('SPOT', dict(color=DK, lw=0)), ('Tara', dict(facecolor='white', edgecolor=MD, lw=0.5))):
        dd = a[a.dataset.str.contains(ds, case=False)]; ax.scatter(dd.ln_I_over_A, dd.D, s=5, zorder=3, label=ds, **kw)
    ax.set_xlabel('ln(I/A)'); ax.set_ylabel('D'); ax.set_title('D is a bounded log-ratio'); ax.legend(loc='upper left', handletextpad=0.1); lab(ax, 'a')
    ax = fig.add_subplot(gs[1]); fc = [c_ for c_ in b.columns if 'fraction' in c_][0]; order_ = ['0.02-0.2um', '0.2-1.2um', '1.2-80um']; rj = np.random.default_rng(1)
    for i, f_ in enumerate(order_):
        v = b[b[fc] == f_].I_plus_A.values; G.raincloud(ax, v, i, [DK, MD, '#8E99AB'][i], width=0.36, seed=i)
    ax.set_xticks(range(3)); ax.set_xticklabels(['Viral', 'Pico', 'Micro/\nnano']); ax.set_ylabel('(I + A) share of PolB'); ax.set_ylim(0.35, 1.02); ax.set_title('Share held by the two orders'); lab(ax, 'b')
    from scipy.stats import norm
    ax = fig.add_subplot(gs[2]); nm = ['SPOT, seasonal SST', 'SPOT, SST anomaly', 'Tara, all stations', 'Tara, 13–24 °C', 'SPOT viral, pre-2014']; assert len(c) == 5; xx = np.linspace(-0.03, 0.30, 400)
    for i, r in enumerate(c.itertuples()):                        # half-eye: the sampling distribution implied by the estimate and its 95% interval, with the interval itself underneath
        y0 = len(c) - 1 - i; col = [DK, DK, MD, MD, DK][i]; se = (r.ci_high - r.ci_low) / (2 * 1.96); dn = norm.pdf(xx, r.value, se); dn = dn / dn.max() * 0.72; inside = (xx >= r.ci_low) & (xx <= r.ci_high)
        ax.fill_between(xx, y0, y0 + dn, color=col, alpha=0.13, lw=0); ax.fill_between(xx[inside], y0, y0 + dn[inside], color=col, alpha=0.30, lw=0); ax.plot(xx, y0 + dn, color=col, lw=0.5)
        ax.plot([r.ci_low, r.ci_high], [y0, y0], color=col, lw=1.1, solid_capstyle='round'); ax.plot(r.value, y0, 'o', color=col, ms=3.8, mec='white', mew=0.5, zorder=4); ax.text(0.298, y0 + 0.12, f'n = {int(r.n)}', ha='right', va='bottom', fontsize=5.5, color=MD)
    ax.axvline(0, color=MD, lw=0.5, ls=(0, (3, 2))); ax.set_yticks(range(5)); ax.set_yticklabels(nm[::-1]); ax.set_ylim(-0.35, 4.95); ax.set_xlim(-0.03, 0.30); ax.set_xlabel('Slope of D on temperature (per °C)'); ax.set_title('Temperature sensitivity'); lab(ax, 'c', 23)
    # d: the fitted relationships themselves, and where each crosses D = 0
    ax = fig.add_subplot(gs[3]); sp = rd('Fig1_derive', 'a_samples'); sp = sp[sp.SST.notna()]; ta = pd.read_csv(SUP + 'Supplementary_Data_4_Tara_Dindex_predictors.csv'); ta = ta[ta.D_index_qc == 'ok'].copy()
    ta['femto'] = (ta.size_fraction == 'Femto-size').astype(int); ta['L'] = np.log((ta.I_relative + 1e-3) / (ta.A_relative + 1e-3)); X = np.c_[np.ones(len(ta)), ta.temperature, ta.femto]; bt = np.linalg.lstsq(X, ta.L.values, rcond=None)[0]
    zs = dict(zip(d.statistic, zip(d.value, d.ci_low, d.ci_high))); curves = []; ax.axhline(0, color=MD, lw=0.5, ls=(0, (3, 2)), zorder=1)
    SP = {'micro': ('SST at D=0 (micro)', 'SPOT micro/nano', '-'), 'pico': ('SST at D=0 (pico)', 'SPOT pico', (0, (4, 1.5))), 'viral': ('SST at D=0 (viral)', 'SPOT viral', (0, (1.2, 1.2)))}
    for k, (key, nm_, ls) in SP.items():
        g_ = sp[sp.frac == k]; b_ = np.polyfit(g_.SST, g_.D_index, 1); assert abs(-b_[1] / b_[0] - zs[key][0]) < 1e-3; t_ = np.linspace(g_.SST.min(), g_.SST.max(), 50); ax.plot(t_, np.polyval(b_, t_), color=DK, lw=1.2, ls=ls, zorder=3, label=nm_)
        curves += [dict(series=nm_, temperature=a_, fitted_D=b) for a_, b in zip(t_, np.polyval(b_, t_))]
    for fem, key, nm_, ls in ((0, 'Tara temperature at D=0, pico', 'Tara pico-size', '-'), (1, 'Tara temperature at D=0, femto', 'Tara femto-size', (0, (4, 1.5)))):
        z0 = -(bt[0] + bt[2] * fem) / bt[1]; assert abs(z0 - zs[key][0]) < 1e-3; t_ = np.linspace(-1.8, 30.5, 120); yy = np.tanh((bt[0] + bt[1] * t_ + bt[2] * fem) / 2); ax.plot(t_, yy, color='#98A2B3', lw=1.2, ls=ls, zorder=2, label=nm_)
        curves += [dict(series=nm_, temperature=a_, fitted_D=b) for a_, b in zip(t_, yy)]
    for n_, (key, col) in enumerate((('Tara temperature at D=0, pico', '#98A2B3'), ('Tara temperature at D=0, femto', '#98A2B3'), ('SST at D=0 (micro)', DK), ('SST at D=0 (pico)', DK), ('SST at D=0 (viral)', DK))):
        v, lo, hi = zs[key]; ax.plot([v], [0], 'o', color=col, ms=3.6, mec='white', mew=0.5, zorder=5)
    # the intervals sit in a strip under the curves so that they do not overprint one another
    axs = ax.inset_axes([0, 0.02, 1, 0.20], sharex=ax); axs.patch.set_alpha(0); order_ = ['Tara temperature at D=0, pico', 'Tara temperature at D=0, femto', 'SST at D=0 (micro)', 'SST at D=0 (pico)', 'SST at D=0 (viral)']
    for n_, key in enumerate(order_):
        v, lo, hi = zs[key]; col = '#98A2B3' if 'Tara' in key else DK; axs.plot([lo, min(hi, 30.5)], [n_, n_], color=col, lw=1.0, solid_capstyle='round'); axs.plot(v, n_, 'o', color=col, ms=3.0, mec='white', mew=0.4)
    axs.set_ylim(5.2, -0.8); axs.axis('off')
    ax.set_xlim(-2, 31); ax.set_ylim(-1.95, 1.62); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.spines['left'].set_bounds(-1, 1); ax.set_xlabel('Temperature (°C)'); ax.set_ylabel('Fitted D'); ax.yaxis.set_label_coords(-0.115, 0.55)
    ax.text(-1.5, -1.2, 'Crossing temperature, 95% interval', fontsize=5.5, color=MD, va='bottom'); ax.legend(loc='upper left', ncol=2, columnspacing=1.0, fontsize=5.5, labelspacing=0.15, handlelength=2.2, borderaxespad=0.0); ax.set_title('Where the orders are equal'); lab(ax, 'd', 18)
    EX = {'d_fitted_curves': pd.DataFrame(curves)}
    G.save(fig, OUT + F); return EX

def ed2():
    F = 'ExtendedDataFig2_heatwave'; b = rd(F, 'b_permutations'); c = rd(F, 'c_same_month'); S = between_year_sensitivity()
    fig = G.figure(183, 66); gs = fig.add_gridspec(1, 3, left=0.065, right=0.985, top=0.92, bottom=0.2, wspace=0.5, width_ratios=[0.8, 1.0, 1.25])
    ax = fig.add_subplot(gs[0]); g = sorted((float(k), v['P']) for k, v in S['grid'].items()); ax.plot([x for x, _ in g], [p for _, p in g], '-o', color=DK, ms=3)
    ax.axhline(0.05, color=MD, lw=0.6, ls=(0, (3, 2))); ax.axvline(S['sd_of_pre_event_year_means'], color=LT, lw=2.2, zorder=0); ax.set_yscale('log'); ax.set_xlabel('Assumed between-year s.d. of D'); ax.set_ylabel('«i»P«/i» of the observed shift'.replace('«i»','').replace('«/i»','')); ax.set_title('Between-year sensitivity'); lab(ax, 'a')
    ax = fig.add_subplot(gs[1]); b = b.sort_values('shift'); ev = b.years.astype(str).str.contains('2014') & b.years.astype(str).str.contains('2015')
    ax.barh(range(len(b)), b['shift'], color=np.where(ev, EV, LT), lw=0); ax.set_yticks(range(len(b))); ax.set_yticklabels([' + '.join(__import__('re').findall(r'20\d\d', str(v))) for v in b.years], fontsize=5.2); ax.axvline(0, color=DK, lw=0.4); ax.set_xlabel('Mean residual shift in D'); ax.set_title('Year-label permutations'); lab(ax, 'b', 14)
    ax = fig.add_subplot(gs[2]); c['date'] = pd.to_datetime(c.date); ax.bar(range(len(c)), c['diff'], color=np.where(c['diff'] >= 0, EV if G.THEME == 'vivid' else '#4D4D4D', LT), lw=0); ax.axhline(0, color=DK, lw=0.4)
    ax.set_xticks(range(len(c))); ax.set_xticklabels([f"{MON[t.month - 1]}" for t in c.date], fontsize=5.5); yr = c.date.dt.year
    for y_ in sorted(yr.unique()): ix = np.where(yr == y_)[0]; ax.text(ix.mean(), -0.19, str(y_), transform=ax.get_xaxis_transform(), ha='center', fontsize=6)
    ax.set_ylabel('D minus pre-event mean of same month'); ax.set_title('Same-month differences'); lab(ax, 'c')
    G.save(fig, OUT + F); return {'a_between_year_sensitivity': pd.DataFrame([dict(between_year_sd=float(k), P=v['P'], null_sd_of_shift=v['null_sd']) for k, v in S['grid'].items()]).sort_values('between_year_sd')}

def ed3():
    F = 'ExtendedDataFig3_global_robustness'; a = rd(F, 'a_station_D_chl'); b = rd(F, 'b_lobo'); c = rd(F, 'c_oof_residuals_by_basin'); d = rd(F, 'd_transfer'); e = rd('Fig3_validate', 'e_permutation_importance')
    fig = G.figure(183, 112); g1 = fig.add_gridspec(1, 3, left=0.065, right=0.985, top=0.95, bottom=0.6, wspace=0.55, width_ratios=[0.8, 1.2, 1.2]); g2 = fig.add_gridspec(1, 2, left=0.065, right=0.985, top=0.44, bottom=0.11, wspace=0.55, width_ratios=[0.8, 1.6])
    ax = fig.add_subplot(g1[0]); sc = ax.scatter(a.log10_chl, a.D, c=a['T'], cmap=G.CMAP_T, vmin=-4, vmax=30, s=11, edgecolor=DK, lw=0.3)
    cax = ax.inset_axes([0.5, 0.93, 0.42, 0.03]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal', ticks=[0, 15, 30]); cb.outline.set_linewidth(0.4); cax.tick_params(length=1.5, width=0.4, labelsize=5.5, pad=1); cax.set_title('Temperature (°C)', fontsize=5.5, pad=2, loc='left', fontweight='normal')
    ax.set_xlabel('log$_{10}$ chlorophyll a (mg m$^{-3}$)'); ax.set_ylabel('Station-mean D'); ax.set_ylim(-1.08, 1.5); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_title('D against chlorophyll'); lab(ax, 'a')
    ax = fig.add_subplot(g1[1]); bb = b[b.basin.notna()]; mods = list(dict.fromkeys(bb.model)); bas = sorted(bb.basin.unique())
    for j, (mname, col) in enumerate(zip(mods[:2], (DK, LT))):
        v = [bb[(bb.basin == k) & (bb.model == mname)].RMSE.mean() for k in bas]; ax.barh(np.arange(len(bas)) + (0.19 - 0.38 * j), v, 0.36, color=col, lw=0, label=['Temperature + size fraction', 'Random forest'][j])
    ax.set_yticks(range(len(bas))); ax.set_yticklabels(bas); ax.set_xlabel('RMSE in held-out basin (D units)'); ax.set_xlim(0, 1.5); ax.set_xticks([0, 0.5, 1.0]); ax.set_title('Error in each held-out basin'); ax.legend(loc='lower right', ncol=1, handletextpad=0.3, borderaxespad=0.2); lab(ax, 'b', 19)
    ax = fig.add_subplot(g1[2]); bs = sorted(c.ocean_basin.unique()); rj = np.random.default_rng(2)
    for i, k in enumerate(bs):
        v = c[c.ocean_basin == k].residual_D.values; ax.scatter(i + rj.uniform(-0.18, 0.18, len(v)), v, s=4, color=LT, lw=0); ax.plot([i - 0.26, i + 0.26], [np.median(v)] * 2, color=DK, lw=1.3)   # kept as points + median: seven narrow groups, one with n = 5, are too cluttered as rainclouds
    ax.axhline(0, color=MD, lw=0.5, ls=(0, (3, 2))); ax.set_xticks(range(len(bs))); ax.set_xticklabels([k.replace(' Ocean', '').replace(' Sea', '').replace('Mediterranean', 'Medit.') for k in bs], rotation=55, ha='right', fontsize=5.5); ax.set_ylabel('Observed − predicted D'); ax.set_title('Residuals by basin'); lab(ax, 'c')
    ax = fig.add_subplot(g2[0]); ax.plot([-1, 1], [-1, 1], color=LT, lw=0.7)
    for k in ('viral', 'pico'):
        dd = d[d.frac == k]; st = G.FRAC[k]; ax.scatter(dd.D_tara, dd.D_index, s=8, marker=st['marker'], facecolor='none', edgecolor=st['color'], lw=0.55, label=st['label'])
    ax.set_xlabel('D predicted from Tara calibration'); ax.set_ylabel('Observed D at SPOT'); ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.set_title('Tara calibration applied to SPOT'); ax.legend(loc='upper left', handletextpad=0.1); lab(ax, 'd')
    ax = fig.add_subplot(g2[1]); e = e.sort_values('mean', ascending=False).reset_index(drop=True); rp = (pd.read_csv(os.path.join(DATA_EXTRA, 'permutation_importance_replicates.csv')) if DATA_EXTRA
          else rd('ExtendedDataFig3_global_robustness', 'e_permutation_replicates')); from scipy.stats import gaussian_kde
    NM = {'temperature': 'Temperature', 'log_chl': 'Chlorophyll a', 'femto': 'Size fraction', 'salinity': 'Salinity', 'log_Fe': 'Iron', 'log_PO4': 'Phosphate', 'log_NO3': 'Nitrate', 'log_Si': 'Silicate'}; xx = np.linspace(-0.25, 1.45, 500); rj = np.random.default_rng(5)
    for i, r in enumerate(e.itertuples()):                   # ridgeline of all 1,000 permutation replicates per predictor, with the mean and central 95% range underneath
        v = rp[rp.feature == r.feature].dR2.values; assert len(v) == r.count and abs(v.mean() - r.mean) < 1e-9; y0 = len(e) - 1 - i; col = DK if i == 0 else MD
        dn = gaussian_kde(v)(xx); dn = dn / dn.max() * 0.78; k_ = dn > 0.006; ax.fill_between(xx[k_], y0, y0 + dn[k_], color=col, alpha=0.22, lw=0, zorder=2); ax.plot(xx[k_], y0 + dn[k_], color=col, lw=0.6, zorder=3)
        ax.scatter(v, y0 - 0.10 - rj.uniform(0, 0.16, len(v)), s=0.6, color=col, alpha=0.25, lw=0, zorder=1, rasterized=True)
        ax.plot([r.ci_low, r.ci_high], [y0, y0], color=col, lw=1.1, solid_capstyle='round', zorder=4); ax.plot(r.mean, y0, 'o', color=col, ms=3.8, mec='white', mew=0.5, zorder=5)
    ax.axvline(0, color=MD, lw=0.5, ls=(0, (3, 2)), zorder=0); ax.set_yticks(range(len(e))); ax.set_yticklabels([NM[x] for x in e.feature][::-1]); ax.set_ylim(-0.4, len(e) - 0.05); ax.set_xlim(-0.25, 1.45)
    ax.set_xlabel('Loss of out-of-fold R² when the predictor is permuted'); ax.set_title('Predictor importance'); lab(ax, 'e', 22)
    G.save(fig, OUT + F); return {'e_permutation_importance': e, 'e_permutation_replicates': rp}

def ed4():
    F = 'ExtendedDataFig4_within_order_heatwave'; a = rd(F, 'a_residuals'); b = rd(F, 'b_affinity_shuffle_null'); c = rd(F, 'c_families'); a['date'] = pd.to_datetime(a.date)
    fig = G.figure(183, 70); gs = fig.add_gridspec(2, 3, left=0.07, right=0.985, top=0.93, bottom=0.17, wspace=0.5, hspace=0.18, width_ratios=[1.5, 0.6, 1.1])
    for i, (o, col) in enumerate((('Imitervirales', G.IMI), ('Algavirales', G.ALG))):
        ax = fig.add_subplot(gs[i, 0]); ax.axvspan(pd.Timestamp('2014-01-01'), pd.Timestamp('2016-01-01'), color=G.BAND, lw=0); ax.axhline(0, color=LT, lw=0.5); ax.plot(a.date, a[o], '-o', color=col, ms=2.2, lw=0.7)
        ax.text(0.01, 0.9, o, transform=ax.transAxes, fontsize=6, color=col, va='top'); ax.set_xlim(pd.Timestamp('2009-05-01'), pd.Timestamp('2016-02-01'))
        if i == 0: ax.tick_params(labelbottom=False); ax.set_title('Within-order thermal affinity through time'); lab(ax, 'a')
    fig.text(0.012, 0.55, 'Seasonal residual of mean thermal affinity (°C)', rotation=90, va='center', fontsize=6.5)
    ax = fig.add_subplot(gs[:, 1]); ax.hist(b.null_event_shift, bins=50, orientation='horizontal', color=LT, lw=0); ax.axhline(b.observed.iloc[0], color=EV, lw=1.0); ax.set_xticks([]); ax.spines['bottom'].set_visible(False); ax.set_ylabel('Event shift (°C)'); ax.set_title('Shuffle null'); lab(ax, 'b')
    ax = fig.add_subplot(gs[:, 2]); c = c.sort_values('mean_tara_affinity').reset_index(drop=True)
    for i, r in c.iterrows():
        col = G.IMI if r.order == 'Imitervirales' else G.ALG; ax.plot([100 * r.mean_share_pre, 100 * r.mean_share_event], [i, i], color=col, lw=0.9); ax.plot(100 * r.mean_share_pre, i, 'o', mfc='white', mec=col, ms=4, mew=0.8); ax.plot(100 * r.mean_share_event, i, 'o', color=col, ms=4)
    ax.plot([], [], 'o', mfc='white', mec=DK, ms=4, mew=0.8, label='2009–2013'); ax.plot([], [], 'o', color=DK, ms=4, label='2014–2015'); ax.legend(loc='upper right', handletextpad=0.1)
    ax.set_yticks(range(len(c))); ax.set_yticklabels([f"{r.family.replace('_', ' ')}  {r.mean_tara_affinity:.1f} °C" for r in c.itertuples()]); ax.set_xlabel('Share of giant-virus PolB coverage (%)'); ax.set_ylim(-0.6, len(c) - 0.4); ax.set_title('Family shares before and during'); lab(ax, 'c', 24)
    G.save(fig, OUT + F)

def ed5():
    F = 'ExtendedDataFig5_host_balance'; b = rd59(F, 'a_station'); c = rd59(F, 'b_variance_partition'); d = rd59(F, 'c_pair_rank')
    fig = G.figure(183, 64); gs = fig.add_gridspec(1, 3, left=0.075, right=0.985, top=0.9, bottom=0.17, wspace=0.36)
    ax = fig.add_subplot(gs[0]); x = 100 * b.mam_share
    sc = ax.scatter(x.clip(lower=1e-3), b.D, c=b['T'], cmap=G.CMAP_T, vmin=-4, vmax=30, s=12, edgecolor=DK, lw=0.3); ax.set_xscale('log'); ax.set_xlabel('Mamiellophyceae (% of 18S reads)'); ax.set_ylabel('Station-mean D'); ax.set_title('D against Mamiellophyceae')
    cax = ax.inset_axes([0.05, 0.1, 0.3, 0.03]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal', ticks=[0, 15, 30]); cb.outline.set_linewidth(0.4); cax.tick_params(length=1.5, width=0.4, labelsize=5.5, pad=1); cax.set_title('Temperature (°C)', fontsize=5.5, pad=2, loc='left', fontweight='normal'); lab(ax, 'a')
    ax = fig.add_subplot(gs[1])
    for i, r in c.iterrows():
        bot = 0
        for k, col, nm in (('unique_temperature', G.N_MID, 'Temperature only'), ('shared', G.N_LIGHT, 'Shared'), ('unique_host', HO, 'Mamiellophyceae only')):
            ax.bar(i, r[k], 0.5, bottom=bot, color=col, lw=0, label=nm if i == 0 else None); bot += r[k]
    ax.set_xticks([0, 1]); ax.set_xticklabels(['All', 'Non-Arctic']); ax.set_xlim(-0.7, 1.7); ax.set_ylim(0, 1.12); ax.set_yticks([0, 0.25, 0.5, 0.75, 1]); ax.set_ylabel('Variance in D explained (R²)'); ax.set_title('Temperature or Mamiellophyceae?'); h, l = ax.get_legend_handles_labels(); ax.legend(h[::-1], l[::-1], loc='upper center', ncol=1, labelspacing=0.2, borderaxespad=0.1); lab(ax, 'b')
    ax = fig.add_subplot(gs[2]); mam = d.X.str.contains('Mamiello') | d.Y.str.contains('Mamiello'); bins = np.linspace(0, 0.85, 18)
    ax.hist([d.abs_rho_D[~mam], d.abs_rho_D[mam]], bins=bins, stacked=True, color=[LT, HO], lw=0, label=['Other pairs', 'Pairs with Mamiellophyceae'])
    pre = d[(d.X.str.contains('Prymnes') & d.Y.str.contains('Mamiello')) | (d.Y.str.contains('Prymnes') & d.X.str.contains('Mamiello'))].abs_rho_D.iloc[0]; ax.plot(pre, -1.2, '^', color=HO, ms=4.5, clip_on=False)
    ax.set_xlabel('|Spearman ρ| with station-mean D'); ax.set_ylabel('Number of protist-group pairs'); ax.set_title('All 210 protist-group pairs'); ax.legend(loc='upper left'); lab(ax, 'c')
    G.save(fig, OUT + F)

def ed6():
    F = 'ExtendedDataFig6_SPOT_host_series'; a = rd(F, 'a_mamiellophyceae_series'); b = rd(F, 'b_climatology'); c = rd(F, 'c_single_year_shifts'); d = rd(F, 'd_e_deseasonalized')
    fig = G.figure(183, 122); g1 = fig.add_gridspec(1, 2, left=0.075, right=0.985, top=0.89, bottom=0.58, wspace=0.26, width_ratios=[2.5, 1]); g2 = fig.add_gridspec(1, 3, left=0.075, right=0.985, top=0.43, bottom=0.1, wspace=0.5)
    ax = fig.add_subplot(g1[0]); a['t'] = pd.to_datetime(a.ym + '-15'); ax.axvspan(pd.Timestamp('2014-01-01'), pd.Timestamp('2016-01-01'), color=G.BAND, lw=0)
    for fr, col, nm in (('1.2-80um', HOL, '1.2–80 µm'), ('0.2-1.2um', HO, '0.2–1.2 µm')):
        dd = a[a.fraction == fr].sort_values('t'); ax.plot(dd.t, (100 * dd.s_M).clip(lower=0.03), '-o', color=col, ms=2.2, lw=0.6, label=nm)
    ax.set_yscale('log'); ax.set_ylim(0.02, 400); ax.set_ylabel('Mamiellophyceae (% of 18S reads)'); ax.set_xlabel('Year'); ax.set_title('Mamiellophyceae through time'); h, l = ax.get_legend_handles_labels(); ax.legend(h[::-1], l[::-1], loc='upper left', ncol=2); lab(ax, 'a')
    ax = fig.add_subplot(g1[1]); z = lambda v: (v - np.nanmean(v)) / np.nanstd(v); oc = [k for k in b.columns if k.startswith('s_M_1.2')][0]
    ax.plot(b.month, z(b.s_M), '-s', color=HO, ms=2.6, label='Mamiellophyceae, 0.2–1.2 µm'); ax.plot(b.month, z(b[oc]), '-s', color=HOL, ms=2.6, ls=(0, (3, 1.5)), label='Mamiellophyceae, 1.2–80 µm'); ax.plot(b.month, z(b.AG_vir), '-o', color=G.ALG, ms=2.8, label='Algavirales share')
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON); ax.set_ylabel('Monthly mean (standardized)'); ax.set_title('Annual cycles'); ax.legend(loc='lower center', fontsize=5.5, labelspacing=0.2, borderaxespad=0.2); ax.axhline(0, color=LT, lw=0.5); lab(ax, 'b')
    ax = fig.add_subplot(g2[0]); rj = np.random.default_rng(4); gps = list(dict.fromkeys(c.group))
    for i, g_ in enumerate(gps):
        dd = c[c.group == g_]; o = dd[dd.year != 2015]; e_ = dd[dd.year == 2015]; G.raincloud(ax, np.log10(o.fold_change.values), i, MD, width=0.34, seed=i, ms=4); ax.scatter([i], np.log10(e_.fold_change), s=30, marker='D', color=EV, edgecolor='white', lw=0.5, zorder=7)
    ax.axhline(0, color=MD, lw=0.5, ls=(0, (3, 2))); ax.set_yticks(np.log10([0.3, 0.5, 1, 2, 3])); ax.set_yticklabels(['0.3', '0.5', '1', '2', '3']);   # drawn in log10 units so that the density is estimated on the plotted scale
    ax.set_xticks(range(len(gps))); ax.set_xticklabels([g_.replace('Mamiellophyceae', 'Mamiello-\nphyceae') for g_ in gps]); ax.set_xlim(-0.6, len(gps) - 0.4)
    ax.set_ylabel('Season-adjusted share (fold)'); ax.set_title('Each year vs the rest'); lab(ax, 'c')
    d['yr'] = d.ym.str[:4].astype(int)
    for j, (col_, nm, k) in enumerate((('AG_pico_resid', 'Algavirales share, picoplankton\n(residual, percentage points)', 'd'), ('AG_vir_resid', 'Algavirales share, viral fraction\n(residual, percentage points)', 'e'))):
        ax = fig.add_subplot(g2[j + 1]); dd = d.dropna(subset=['M_resid', col_]); ax.axhline(0, color=LT, lw=0.5); ax.axvline(0, color=LT, lw=0.5)
        o = dd[dd.yr != 2015]; e_ = dd[dd.yr == 2015]; ax.scatter(o.M_resid, 100 * o[col_], s=10, facecolor='white', edgecolor=DK, lw=0.55); ax.scatter(e_.M_resid, 100 * e_[col_], s=12, color=EV, lw=0)
        ax.set_xlabel('Mamiellophyceae share (residual, log$_{10}$)'); ax.set_ylabel(nm); ax.set_title(['Coupling, picoplankton', 'Coupling, viral fraction'][j]); lab(ax, k)
    G.save(fig, OUT + F)

def esm_lines(ax, d, ycol, xcol='year'):
    for e_, st in G.ESM.items():
        for hm, fill in (('North', DK), ('South', 'white')):
            dd = d[(d.ESM == e_) & (d.hemisphere == hm)].sort_values(xcol)
            if len(dd): ax.plot(dd[xcol], dd[ycol].abs(), marker=st['marker'], color=MD, lw=0.6, mfc=fill, mec=DK, ms=2.8, mew=0.5)
def esm_legend(ax, **kw):
    hs = [mpl.lines.Line2D([], [], marker=st['marker'], color='none', mfc='white', mec=DK, ms=3.2, mew=0.6, label=e_) for e_, st in G.ESM.items()]
    hs += [mpl.lines.Line2D([], [], marker='o', color='none', mfc=DK, mec=DK, ms=3.2, label='North'), mpl.lines.Line2D([], [], marker='o', color='none', mfc='white', mec=DK, ms=3.2, mew=0.6, label='South')]
    ax.legend(handles=hs, fontsize=5.5, handletextpad=0.1, **kw)

def ed7():
    F = 'ExtendedDataFig7_projection'; m = rd(F, 'a_b_c_maps_2deg_4model_mean'); d = rd(F, 'd_front_trajectories'); e = rd(F, 'e_rates_by_projection'); f = rd(F, 'f_agreement_by_band')
    fig = G.figure(183, 104); cm = G.cmap_orders()
    for i, (col, vm, k) in enumerate((('D_2000s', 1, 'a'), ('D_2090s', 1, 'b'), ('dD', 0.4, 'c'))):
        ax = fig.add_axes([0.02 + i * 0.33, 0.6, 0.3, 0.34], projection=ccrs.Robinson(central_longitude=-160)); ax.set_global(); ax.spines['geo'].set_linewidth(0.5)
        lon, lat, z = grid(m, col); pm = ax.pcolormesh(lon, lat, z, cmap=cm, vmin=-vm, vmax=vm, transform=ccrs.PlateCarree(), rasterized=True, shading='nearest'); ax.add_feature(cf.LAND, facecolor='#CFCFCF', lw=0, zorder=4)
        if col != 'dD': ax.contour(lon, lat, z, levels=[0], colors=[DK], linewidths=0.8, transform=ccrs.PlateCarree(), zorder=5)
        else: dis = m[(~m.all_models_agree.astype(bool)) & m.dD.notna()].iloc[::3]; ax.scatter(dis.lon, dis.lat, s=0.25, color=DK, lw=0, transform=ccrs.PlateCarree())
        fig.text(0.02 + i * 0.33, 0.945, k, fontsize=8, fontweight='bold'); fig.text(0.05 + i * 0.33, 0.945, ['Mean D, 2000s', 'Mean D, 2090s', 'Change in D, 2000s to 2090s'][i], fontsize=7, fontweight='bold')
        if i in (1, 2):
            cax = fig.add_axes([0.02 + i * 0.33 + (0.0 if i == 2 else -0.165) + 0.07, 0.575, 0.16, 0.011]); cb = fig.colorbar(pm, cax=cax, orientation='horizontal', ticks=[-vm, 0, vm]); cb.outline.set_linewidth(0.4); cax.tick_params(length=2, width=0.5, labelsize=6); cb.set_label('D' if i == 1 else 'ΔD, 2090s − 2000s', labelpad=1)
    ax = fig.add_axes([0.065, 0.1, 0.25, 0.34]); dd = d[d.projection.str.contains('final ensemble')]; esm_lines(ax, dd, 'front'); ax.set_xlabel('Decade (mid-year)'); ax.set_ylabel('Latitude of D = 0 front (°)'); ax.set_title('Front latitude by model'); esm_legend(ax, loc='upper left', ncol=2, columnspacing=0.6); ax.set_ylim(33, 58); lab(ax, 'd')
    ax = fig.add_axes([0.47, 0.1, 0.21, 0.34]); pr = [p_ for p_ in dict.fromkeys(e.projection) if 'archived' not in p_]
    for i, p_ in enumerate(pr):
        for hm, off, fill in (('North', 0.14, DK), ('South', -0.14, 'white')):
            x = e[(e.projection == p_) & (e.hemisphere == hm)]; y = len(pr) - 1 - i + off; ax.plot([x.km_per_decade.min(), x.km_per_decade.max()], [y, y], color=MD, lw=0.6)
            for _, r in x.iterrows(): ax.plot(r.km_per_decade, y, marker=G.ESM[r.ESM]['marker'], mfc=fill, mec=DK, ms=3, mew=0.5)
    short = {'final ensemble (7 predictors)': 'Final ensemble', 'final ensemble, AOA cells only': 'Ensemble,\nsupport domain'}; ax.set_yticks(range(len(pr))); ax.set_yticklabels([short.get(p_, p_.split('(')[0].strip().capitalize()) for p_ in pr][::-1]); ax.set_xlim(0, 125); ax.set_xlabel('Poleward rate (km per decade)'); ax.set_title('Rates across designs'); lab(ax, 'e', 22)
    ax = fig.add_axes([0.78, 0.1, 0.2, 0.34]); fb = f[f.hemisphere == 'both']; ax.bar(range(len(fb)), 100 * fb.share_sign_agreement, 0.62, color=MD, lw=0); ax.set_xticks(range(len(fb))); ax.set_xticklabels(['<30°', '30–45°', '45–60°', '>60°'][:len(fb)]); ax.set_ylim(0, 100)
    ax.set_xlabel('Absolute latitude'); ax.set_ylabel('Models agree on sign (% of area)'); ax.set_title('Model agreement'); lab(ax, 'f', 11)
    G.save(fig, OUT + F)

def ed8():
    F = 'ExtendedDataFig8_projection_robustness'; m = rd(F, 'a_b_aoa_share_2deg'); c = rd(F, 'c_aoa_shares'); u = rd(F, 'd_units'); e = rd(F, 'e_temperature_only_fronts'); f = rd(F, 'f_profiles'); g = rd(F, 'g_front_vs_isotherm'); h = rd(F, 'h_isotherm_trajectories')
    fr = rd('ExtendedDataFig7_projection', 'd_front_trajectories'); fr = fr[fr.projection.str.contains('final ensemble')]
    fig = G.figure(183, 150)
    for i, (col, k) in enumerate((('aoa_2000s', 'a'), ('aoa_2090s', 'b'))):
        ax = fig.add_axes([0.02 + i * 0.33, 0.72, 0.3, 0.24], projection=ccrs.Robinson(central_longitude=-160)); ax.set_global(); ax.spines['geo'].set_linewidth(0.5)
        lon, lat, z = grid(m, col); pm = ax.pcolormesh(lon, lat, 100 * z, cmap='Greys_r', vmin=0, vmax=100, transform=ccrs.PlateCarree(), rasterized=True, shading='nearest'); ax.add_feature(cf.LAND, facecolor='#E9E3D8', lw=0, zorder=4); fig.text(0.02 + i * 0.33, 0.962, k, fontsize=8, fontweight='bold'); fig.text(0.05 + i * 0.33, 0.962, ['Support domain, 2000s', 'Support domain, 2090s'][i], fontsize=7, fontweight='bold')
    cax = fig.add_axes([0.23, 0.705, 0.2, 0.01]); cb = fig.colorbar(pm, cax=cax, orientation='horizontal', ticks=[0, 50, 100]); cb.outline.set_linewidth(0.4); cax.tick_params(length=2, width=0.5, labelsize=6); cb.set_label('Grid cells within support domain (%)', labelpad=1)
    ax = fig.add_axes([0.76, 0.74, 0.22, 0.2]); c0 = c[c.spec == c.spec.iloc[0]]; dec = [c0.decade.min(), c0.decade.max()]
    for j, (dc, col) in enumerate(zip(dec, ('white', MD))):
        for i, k in enumerate(('share_all', 'share_near_front', 'share_polar')):
            v = c0[c0.decade == dc][k] * 100; ax.bar(i + (j - 0.5) * 0.36, v.mean(), 0.34, color=col, edgecolor=DK, lw=0.5, label=['2000s', '2090s'][j] if i == 0 else None); ax.plot(np.full(len(v), i + (j - 0.5) * 0.36), v, '.', color=DK, ms=2.5)
    ax.set_xticks(range(3)); ax.set_xticklabels(['All\nocean', 'Near\nfronts', '≥60°']); ax.set_ylim(0, 118); ax.set_ylabel('Area within support (%)'); ax.set_title('Share within support'); ax.legend(loc='upper center', ncol=2, bbox_to_anchor=(0.5, 1.28)); lab(ax, 'c')
    ax = fig.add_axes([0.065, 0.39, 0.24, 0.22]); uu = u.groupby('variable').agg(ratio=('median_ratio_cmip6_over_tara', 'median'), rho=('spearman_rho', 'median')).reset_index(); ax.axvline(1, color=LT, lw=0.5)
    fin = ~uu.variable.isin(['Si', 'CO3']); ax.scatter(uu.ratio[fin], uu.rho[fin], s=14, color=DK, lw=0); ax.scatter(uu.ratio[uu.variable == 'Si'], uu.rho[uu.variable == 'Si'], s=14, facecolor='white', edgecolor=DK, lw=0.7); ax.scatter(uu.ratio[uu.variable == 'CO3'], uu.rho[uu.variable == 'CO3'], s=16, marker='x', color=DK, lw=0.9); NM = {'temperature': ('Temperature', (4, 3)), 'salinity': ('Salinity', (4, -5)), 'Fe': ('Iron', (4, 0)), 'PO4': ('Phosphate', (4, 3)), 'chlorophyll_a': ('Chlorophyll a', (4, -5)), 'NO3': ('Nitrate', (4, 0)), 'Si': ('Silicate', (4, 0)), 'CO3': ('Excluded field', (-4, 6))}
    [ax.annotate(NM[r.variable][0], (r.ratio, r.rho), xytext=NM[r.variable][1], textcoords='offset points', fontsize=5.5, va='center', ha='right' if r.variable == 'CO3' else 'left') for r in uu.itertuples()]
    ax.set_xscale('log'); ax.set_xlabel('Median ratio, CMIP6 / Tara'); ax.set_ylabel('Spearman ρ across 83 stations'); ax.set_ylim(0, 1.12); ax.set_title('Model inputs vs observations'); lab(ax, 'd')
    ax = fig.add_axes([0.40, 0.39, 0.24, 0.22]); esm_lines(ax, e, 'front'); ax.set_xlabel('Decade (mid-year)'); ax.set_ylabel('Front latitude (°)'); ax.set_title('Temperature-only calibration'); lab(ax, 'e')
    ax = fig.add_axes([0.735, 0.39, 0.245, 0.22]); pers = [f.period.min(), f.period.max()]
    for e_, st in G.ESM.items():
        for pr_, ls in zip(pers, ((0, (1, 1)), '-')):
            dd = f[(f.ESM == e_) & (f.period == pr_)].sort_values('latitude_mid'); ax.plot(dd.mean_D, dd.latitude_mid, color=MD if ls == '-' else LT, lw=0.7, ls=ls)
    ax.axvline(0, color=LT, lw=0.5); ax.plot([], [], color=LT, lw=0.9, ls=(0, (1, 1)), label='2000s'); ax.plot([], [], color=MD, lw=0.9, label='2090s'); ax.legend(loc='lower right'); ax.set_xlabel('Zonal-mean D'); ax.set_ylabel('Latitude (°)'); ax.set_xlim(-1, 1); ax.set_title('Zonal-mean D by model'); lab(ax, 'f')
    ax = fig.add_axes([0.065, 0.075, 0.24, 0.22]); ax.plot([0, 95], [0, 95], color=LT, lw=0.7)
    for _, r in g.iterrows(): ax.plot(r.isotherm, r.D_front, marker=G.ESM[r.ESM]['marker'], mfc=DK if r.hemisphere == 'North' else 'white', mec=DK, ms=4, mew=0.6)
    ax.set_xlabel('Isotherm rate (km per decade)'); ax.set_ylabel('D = 0 front rate (km per decade)'); ax.set_xlim(0, 95); ax.set_ylim(0, 108); ax.set_title('Front rate vs isotherm rate'); esm_legend(ax, loc='upper left'); lab(ax, 'g')
    ax = fig.add_axes([0.40, 0.075, 0.58, 0.22]); esm_lines(ax, fr, 'front')
    for e_ in G.ESM:
        for hm in ('North', 'South'):
            dd = h[(h.ESM == e_) & (h.hemisphere == hm)].sort_values('year'); ax.plot(dd.year, dd.lat.abs(), color=MD, lw=0.7, ls=(0, (1, 1)))
    ax.plot([], [], color=MD, lw=0.7, label='D = 0 front'); ax.plot([], [], color=MD, lw=0.7, ls=(0, (1, 1)), label='Isotherm'); ax.legend(loc='upper left', ncol=2); ax.set_xlabel('Decade (mid-year)'); ax.set_ylabel('Absolute latitude (°)'); ax.set_title('Fronts and their isotherms'); lab(ax, 'h')
    G.save(fig, OUT + F)

def ed9():
    F = 'ExtendedDataFig9_marker_recovery'; a = rd(F, 'a_scores'); b = rd(F, 'b_mantel_haenszel'); c = rd(F, 'c_pathways_by_score'); sc = [k for k in a.columns if 'score' in k][0]; sc2 = [k for k in c.columns if 'score' in k][0]
    fig = G.figure(183, 62); gs = fig.add_gridspec(1, 3, left=0.06, right=0.985, top=0.92, bottom=0.2, wspace=0.34, width_ratios=[0.62, 1.7, 0.9])
    ax = fig.add_subplot(gs[0]); ct = a.groupby(['order', sc]).size()
    for j, (o, col) in enumerate((('Imitervirales', G.IMI), ('Algavirales', G.ALG))):
        for s_, n in ct[o].items(): ax.barh(s_ + (0.19 - 0.38 * j), n, 0.36, color=col, lw=0)
    ax.set_ylabel('Six-marker recovery score'); ax.set_xlabel('Genomes'); ax.set_yticks(sorted(a[sc].unique())); ax.plot([], [], 's', color=G.IMI, ms=4, label='Imitervirales'); ax.plot([], [], 's', color=G.ALG, ms=4, label='Algavirales'); ax.set_title('Marker recovery'); ax.legend(loc='lower right', handletextpad=0.1); lab(ax, 'a')
    # b: prevalence of each trait within each recovery stratum shared by the two orders -- the comparison the stratified test makes, shown directly (and defined even where an odds ratio is not)
    ax = fig.add_subplot(gs[1]); gt = pd.read_csv(SUP + 'Supplementary_Data_7_GVMAG_functional_traits.csv'); gt = gt[gt.order.isin(['Imitervirales', 'Algavirales'])]
    TC = {'Heat-shock proteins': ('pw_heat_shock_proteins', 'Heat-\nshock'), 'DNA repair': ('pw_DNA_repair', 'DNA\nrepair'), 'Cold-shock proteins': ('pw_cold_shock_proteins', 'Cold-\nshock'), 'Type-I rhodopsin': ('type_I_rhodopsin', 'Type-I\nrhodopsin'), 'TCA cycle': ('pw_TCA', 'TCA\ncycle')}; assert list(b.trait) == list(TC)
    strata = [4, 5, 6]; rows = []
    def wil(k, n, z=1.96):
        p_ = k / n; den = 1 + z * z / n; c_ = (p_ + z * z / (2 * n)) / den; h = z * np.sqrt(p_ * (1 - p_) / n + z * z / (4 * n * n)) / den; return c_ - h, c_ + h
    for ti, (tr_, (col_, short)) in enumerate(TC.items()):
        if ti % 2 == 0: ax.axvspan(ti * 4 - 0.75, ti * 4 + 2.75, color=LT, alpha=0.22, lw=0, zorder=0)
        for o, colr, off in (('Imitervirales', G.IMI, -0.13), ('Algavirales', G.ALG, 0.13)):
            xs, ys = [], []
            for si, st in enumerate(strata):
                g_ = gt[(gt.order == o) & (gt.marker_recovery_score == st)]; n = len(g_); k = int((g_[col_] > 0).sum()); lo, hi = wil(k, n); x_ = ti * 4 + si + off; xs.append(x_); ys.append(100 * k / n)
                ax.plot([x_, x_], [100 * lo, 100 * hi], color=colr, lw=0.8, zorder=2); rows.append(dict(trait=tr_, order=o, marker_recovery_score=st, genomes_with_trait=k, genomes=n, percent=100 * k / n, wilson_low=100 * lo, wilson_high=100 * hi))
            ax.plot(xs, ys, '-o', color=colr, lw=0.9, ms=3.2, mec='white', mew=0.4, zorder=3)
    for r in b.itertuples():   # the strata reproduce the deposited totals
        assert r.Imitervirales_positive == int((gt[gt.order == 'Imitervirales'][TC[r.trait][0]] > 0).sum()) and r.Algavirales_positive == int((gt[gt.order == 'Algavirales'][TC[r.trait][0]] > 0).sum())
    ax.set_xticks([ti * 4 + si for ti in range(5) for si in range(3)]); ax.set_xticklabels(['4', '5', '6'] * 5, fontsize=5.5); ax.set_xlim(-0.75, 18.75); ax.set_ylim(-4, 136); ax.set_yticks([0, 25, 50, 75, 100]); ax.spines['left'].set_bounds(0, 100)
    for ti, (tr_, (col_, short)) in enumerate(TC.items()): ax.text(ti * 4 + 1, 135, short, ha='center', va='top', fontsize=5.8, linespacing=1.25)
    ax.set_xlabel('Six-marker recovery score'); ax.set_ylabel('Genomes with trait (%)'); ax.set_title('Traits within recovery strata'); lab(ax, 'b', 10)
    EX = {'b_prevalence_by_stratum': pd.DataFrame(rows)}
    ax = fig.add_subplot(gs[2])
    for o, col, off in (('Imitervirales', G.IMI, -0.06), ('Algavirales', G.ALG, 0.06)):
        g_ = c[c.order == o].groupby(sc2).n_pathways_present.agg(['mean', 'sem', 'count']); g_ = g_[g_['count'] >= 3]; ax.errorbar(g_.index + off, g_['mean'], yerr=1.96 * g_['sem'], fmt='-o', color=col, ms=3, lw=0.8, capsize=1.5)
    ax.set_xticks(sorted(c[sc2].unique())); ax.set_xlabel('Six-marker recovery score'); ax.set_ylabel('Pathway categories per genome'); ax.set_title('Pathway breadth by recovery'); lab(ax, 'c')
    G.save(fig, OUT + F); return EX

if __name__ == '__main__':
    import pickle; os.makedirs(OUT, exist_ok=True); which = sys.argv[1:] or [str(i) for i in range(1, 10)]; extra = {}
    for w in which:
        r = globals()['ed' + w]()
        if r: extra['ED' + w] = r
        print('built ED', w)
    if extra: pickle.dump(extra, open(OUT + f'extra_source_ed_{"".join(which)}.pkl', 'wb'))
