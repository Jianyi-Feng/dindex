"""Main figures 1-5 of the manuscript in the shared house style: colour = the two orders only; no panel titles;
no statistics or sentences inside panels (they are in the legends). Drawn from the deposited source data
(processed_tables/figure_source_data). Run this script from anywhere:

    python figures/v73/build_main.py        # all main figures
    python figures/v73/build_main.py 2      # only Fig. 2

Paths resolve from the repository root; the data folder can be moved with DINDEX_INPUTS and the output folder
with DINDEX_OUT. A few working files of the original author-side job are optional and each has an exact
fallback to the deposited tables (see figures/v73/README.md).
"""
import os, sys, json, itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import stats
import cartopy.crs as ccrs, cartopy.feature as cf

import gvstyle as G; G.use()

# ---- paths: everything resolves from the repository, never from a machine-specific absolute path --------
PKG = Path(__file__).resolve().parents[2]
INPUTS = Path(os.environ.get('DINDEX_INPUTS') or (PKG / 'processed_tables'))
_OUT = os.environ.get('DINDEX_OUT') or str(PKG / 'outputs' / 'v73')
OUT = str(Path(_OUT)) + os.sep
SD = SD59 = str(INPUTS / 'figure_source_data') + os.sep
SUP = str(INPUTS / 'supplementary_data') + os.sep
RAW = str(PKG / 'raw_tables') + os.sep
STUDY_JSON = os.environ.get('DINDEX_STUDY_JSON')      # optional: .../strengthening_2026-11-23/results
DATA_EXTRA = os.environ.get('DINDEX_DATA_EXTRA')      # optional: .../data_extra

# The deposited V72 package renamed three source-data sheets; map the working names used by the drawing code.
SHEET_ALIAS = {
    ('Fig3_validate', 'e_permutation_importance'): ('ExtendedDataFig3_global_robustness', 'e_permutation_importance'),
    ('Fig3_validate', 'f_year_assignments'): ('Fig3_validate', 'e_year_assignments'),
    ('Fig3_validate', 'g_per_phylotype'): ('Fig3_validate', 'f_per_phylotype'),
}


def rd(f, s):
    """Read one deposited figure source-data sheet (with the rename map above applied)."""
    f, s = SHEET_ALIAS.get((f, s), (f, s))
    return pd.read_excel(SD + f + '_source_data.xlsx', s)


rd59 = rd

# The original job imported `polb_cti.timeseries` from a sibling pipeline for the Newey-West slope in Fig. 2d.
# Set DINDEX_POLB_CTI to that folder to use it. Without it the fallback below is used: Fig. 2d draws the fitted
# regression line, and Newey-West changes the standard error, not the point estimate, so the OLS slope is the
# same number.
try:
    _polb = os.environ.get('DINDEX_POLB_CTI')
    if _polb:
        sys.path.insert(0, _polb)
    from polb_cti import timeseries as ts
except Exception:                                       # optional dependency, absent in the shipped package
    class _Timeseries:
        @staticmethod
        def _nw_slope(x, y):
            x = np.asarray(x, float); y = np.asarray(y, float)
            return {'slope': np.polyfit(x, y, 1)[0]}
    ts = _Timeseries()


def within_order_json():
    """Fig. 3g input: within-order Spearman correlations of the seasonal residual of D with niche axes.
    With DINDEX_STUDY_JSON set this reads the author-local analysis JSON (.../strengthening_2026-11-23/
    results/03_within_order.json); otherwise it rebuilds the same coefficients from the deposited sheet
    Fig3_validate::g_niche_axes."""
    if STUDY_JSON:
        with open(os.path.join(STUDY_JSON, '03_within_order.json'), encoding='utf-8') as fh:
            return json.load(fh)
    d = rd('Fig3_validate', 'g_niche_axes')
    c2 = {}
    for o, g in d.groupby('order', sort=False):
        c2[o] = {r.niche_axis: {'rho': float(r.spearman_rho), 'P': float(r.P), 'n': int(r.n)} for r in g.itertuples()}
    return {'C2': c2}

MON = list('JFMAMJJASOND')
def lab(ax, s, x_mm=9, y=0.008):
    G.label(ax, s)          # deferred: placed at save time, left of the y-axis label text

# ============================================================== Fig. 1
def fig1():
    G.use(1.16)
    a = rd('Fig1_derive', 'a_samples'); a['date'] = pd.to_datetime(a.date)
    fig = G.figure(183, 104); gs = fig.add_gridspec(2, 3, left=0.07, right=0.985, top=0.945, bottom=0.085, hspace=0.46, wspace=0.34)
    gs1 = fig.add_gridspec(1, 3, left=0.097, right=0.985, top=0.945, bottom=0.5954, wspace=0.34)
    ax = fig.add_subplot(gs1[0, 0]); ax.axvspan(pd.Timestamp('2014-01-01'), pd.Timestamp('2016-08-31'), color=G.BAND, lw=0)
    for i, k in enumerate(['viral', 'pico', 'micro']):
        d = a[a.frac == k]; ax.scatter(d.date, np.full(len(d), i), s=7, marker=G.FRAC[k]['marker'], color=G.FRAC[k]['color'], lw=0)
    ax.set_yticks([0, 1, 2]); ax.set_yticklabels(['Viral', 'Pico', 'Micro/nano']); ax.set_ylim(-0.7, 2.7); ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    ax.set_xlim(pd.Timestamp('2009-01-01'), pd.Timestamp('2020-01-01')); ax.set_xlabel('Year'); ax.set_title('Sampling, 2009–2019'); lab(ax, 'a', 11)
    for j, (sh, k) in enumerate((('b_viral_monthly', 'b'), ('c_pico_monthly', 'c'), ('d_micro_monthly', 'd'))):
        d = rd('Fig1_derive', sh); ax = fig.add_subplot(gs1[0, 1] if j == 0 else gs1[0, 2] if j == 1 else gs[1, 0])
        for q, c, nm in (('I', G.IMI, 'Imitervirales'), ('A', G.ALG, 'Algavirales')):
            ax.fill_between(d.month, d[q + 'q1'], d[q + 'q3'], color=c, alpha=0.18, lw=0); ax.plot(d.month, d[q + 'm'], '-o', color=c, ms=2.6, label=nm)
        ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON); ax.set_ylim(0, 1.22); ax.set_yticks([0, 0.25, 0.5, 0.75, 1]); ax.set_ylabel('Relative abundance'); ax.set_xlabel('Month')
        ax.set_title(['Viral fraction, 0.02–0.2 µm', 'Picoplankton, 0.2–1.2 µm', 'Micro/nanoplankton, 1.2–80 µm'][j])
        if j == 0: ax.legend(loc='upper center', ncol=2, columnspacing=1.2, borderaxespad=0.1)
        lab(ax, k)
    ax = fig.add_subplot(gs[1, 1]); ax.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2)))
    from scipy import stats as _sst
    dd = a.dropna(subset=['SST']).copy()
    _fit = {}
    for k in ['viral', 'pico', 'micro']:
        d = dd[dd.frac == k]; st = G.FRAC[k]
        ax.scatter(d.SST, d.D_index, s=8, marker=st['marker'], facecolor='none', edgecolor=st['color'], lw=0.55)
        _sl, _ic = np.polyfit(d.SST, d.D_index, 1)
        xs = np.linspace(d.SST.min(), d.SST.max(), 50)
        ax.plot(xs, _ic + _sl * xs, color=st['color'], lw=1.1)
        _rr, _pp = _sst.pearsonr(d.SST, d.D_index)
        _fit[k] = dict(slope=float(_sl), intercept=float(_ic), r=float(_rr), P=float(_pp), n=int(len(d)))
    for _i, _k in enumerate(['viral', 'pico', 'micro']):
        _yy = 0.985 - 0.082 * _i; _st_ = G.FRAC[_k]; _f = _fit[_k]
        ax.plot([0.045], [_yy], marker=_st_['marker'], ms=3.0, mfc='none', mec=_st_['color'], mew=0.6,
                ls='none', transform=ax.transAxes, clip_on=False)
        _pt = '$P$ < 0.001' if _f['P'] < 1e-3 else '$P$ = ' + G.sci(_f['P'])
        ax.text(0.085, _yy, '%s  $r$ = %.2f, %s' % (_st_['label'], _f['r'], _pt), transform=ax.transAxes,
                va='center', ha='left', fontsize=G.fs(4.9), color=G.INK)
    ax.set_xlabel('SST (\u00b0C)'); ax.set_ylabel('D'); ax.set_title('D against temperature')
    ax.set_ylim(-1.05, 1.55); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); lab(ax, 'e')
    f = rd('Fig1_derive', 'f_monthly_D'); ax = fig.add_subplot(gs[1, 2]); ax.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2)))
    sm = rd('Fig1_derive', 'a_samples').copy(); sm['month'] = pd.to_datetime(sm.date).dt.month; fm = sm.groupby(['frac', 'month']).D_index.agg(['mean', 'size']).reset_index()
    FL = {'viral': ('#8D99AE', '-', 'o', 'Viral'), 'pico': ('#AEB8C8', (0, (4, 1.5)), 's', 'Pico'), 'micro': ('#C6CEDA', (0, (1.2, 1.2)), '^', 'Micro/nano')}     # light tints: context behind the pooled mean
    for k, (c_, ls, mk, nm) in FL.items():
        d_ = fm[fm.frac == k]; ax.plot(d_.month, d_['mean'], ls=ls, marker=mk, color=c_, lw=0.9, ms=2.2, mew=0, zorder=2, label=nm)
    ax.fill_between(f.month, f['mean'] - 1.96 * f['sem'], f['mean'] + 1.96 * f['sem'], color=G.INK, alpha=0.10, lw=0, zorder=3); ax.plot(f.month, f['mean'], '-o', color=G.INK, ms=2.8, lw=1.2, zorder=4, label='All fractions')
    ax.set_xticks(range(1, 13)); ax.set_xticklabels(MON); ax.set_ylim(-1.05, 1.05); ax.set_ylabel('D'); ax.set_xlabel('Month'); ax.set_title('Seasonal cycle of D')
    ax.legend(loc='upper left', ncol=2, columnspacing=0.8, handlelength=1.5, handletextpad=0.3, fontsize=G.fs(5.5), borderaxespad=0.1); lab(ax, 'f')
    EX1 = {'f_monthly_D_by_fraction': fm.rename(columns={'frac': 'size_fraction', 'mean': 'mean_D', 'size': 'n_samples'})}
    EX1['e_regression'] = pd.DataFrame([dict(size_fraction=k, slope=v['slope'], intercept=v['intercept'], r=v['r'], P=v['P'], n=v['n']) for k, v in _fit.items()])
    G.save(fig, OUT + 'Fig1_derive'); return EX1

# ============================================================== Fig. 2
def fig2():
    G.use(1.16)
    rng = np.random.default_rng(20261123)
    m = pd.read_csv(RAW + 'SPOT_sample_Dindex.csv'); m['date'] = pd.to_datetime(m.date); m['year'] = m.date.dt.year; m['month'] = m.date.dt.month
    v = m[m.size_fraction_ascii == '0.02-0.2um'].sort_values('date').reset_index(drop=True)
    yrs, mo, y = v.year.values, v.month.values.astype(float), v.D_index.values; pre = yrs < 2014; YE = sorted(set(yrs))
    H = lambda q: np.column_stack([np.ones(len(q)), np.sin(2 * np.pi * q / 12), np.cos(2 * np.pi * q / 12)])
    def shift(yy, ev):
        e = np.isin(yrs, ev); b = np.linalg.lstsq(H(mo[~e]), yy[~e], rcond=None)[0]; r = yy - H(mo) @ b; return r[e].mean() - r[~e].mean()
    b0 = np.linalg.lstsq(H(mo[pre]), y[pre], rcond=None)[0]; res = y - H(mo) @ b0; obs = shift(y, (2014, 2015)); assign = {p: shift(y, p) for p in itertools.combinations(YE, 2)}
    rp = res[pre]; tp = (v.date[pre] - v.date[pre].iloc[0]).dt.days.values / 30.44; gap = np.diff(tp); ok = (gap > .5) & (gap < 1.5)
    phi = np.corrcoef(rp[:-1][ok], rp[1:][ok])[0, 1]; sd = rp.std(ddof=3); tall = (v.date - v.date.iloc[0]).dt.days.values / 30.44; seas = H(mo) @ b0; null = np.empty(20000)
    for k in range(len(null)):
        e = np.empty(len(tall)); e[0] = rng.normal(0, sd)
        for i in range(1, len(tall)):
            p = phi ** max(tall[i] - tall[i - 1], 0.0); e[i] = p * e[i - 1] + rng.normal(0, sd * np.sqrt(max(1 - p * p, 0)))
        null[k] = shift(seas + e, (2014, 2015))
    fig = G.figure(183, 100)
    g1 = fig.add_gridspec(1, 2, left=0.062, right=0.985, top=0.955, bottom=0.60, wspace=0.34); g2 = fig.add_gridspec(1, 3, left=0.062, right=0.985, top=0.46, bottom=0.105, wspace=0.42, width_ratios=[1.12, 0.98, 0.98])
    axa, axb = fig.add_subplot(g1[0]), fig.add_subplot(g1[1]); axc, axd, axe = fig.add_subplot(g2[0]), fig.add_subplot(g2[1]), fig.add_subplot(g2[2])
    span = (pd.Timestamp('2014-01-01'), pd.Timestamp('2016-01-01'))
    mon = v.groupby(v.date.dt.to_period('M')).agg(a=('SST_anomaly', 'mean'), mei=('MEI', 'mean')); t = mon.index.to_timestamp()
    axa.axvspan(*span, color=G.BAND, lw=0); axa.bar(t, mon.a, width=24, color=np.where(mon.a >= 0, G.N_DARK, G.N_MID if G.THEME == 'vivid' else '#BDBDBD'), lw=0); axa.axhline(0, color=G.INK, lw=0.4, ls=(0, (3, 2))); axa.set_ylabel('SST anomaly (°C)')
    a2 = axa.twinx(); a2.plot(t, mon.mei, color=G.MUTED, lw=0.7, ls=(0, (3, 1.5))); a2.set_ylabel('MEI', color=G.MUTED); a2.spines['right'].set_visible(True); a2.spines['right'].set_color(G.MUTED); a2.tick_params(colors=G.MUTED); a2.set_ylim(-3.4, 3.4); axa.set_ylim(-3.6, 3.9); axa.set_title('Heatwave and El Niño at SPOT')
    axb.axvspan(*span, color=G.BAND, lw=0); axb.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2))); tt = pd.date_range('2009-07-01', '2016-01-01', freq='7D')
    axb.plot(tt, H(tt.month + tt.day / 30) @ b0, color=G.MUTED, lw=0.8, ls=(0, (3, 1.5)), label='Seasonal expectation'); axb.plot(v.date, y, color=G.LIGHT, lw=0.6, zorder=1)
    axb.scatter(v.date[pre], y[pre], s=9, color=G.PRE, lw=0, zorder=3, label='2009–2013'); axb.scatter(v.date[~pre], y[~pre], s=11, color=G.EVT, lw=0, zorder=3, label='2014–2015')
    axb.set_ylim(-1.05, 1.45); axb.set_ylabel('D'); axb.legend(loc='upper left', handlelength=1.8); axb.set_title('D against its seasonal cycle')
    for a_ in (axa, axb): a_.set_xlim(pd.Timestamp('2009-05-01'), pd.Timestamp('2016-02-01'))
    axc.hist(null, bins=60, color=G.LIGHT, lw=0); ym = axc.get_ylim()[1]; others = [s for p, s in assign.items() if p != (2014, 2015)]
    axc.plot(others, np.full(len(others), -0.045 * ym), '|', color=G.PRE, ms=5, mew=0.8, clip_on=False); axc.plot([obs], [-0.045 * ym], '|', color=G.EVT, ms=7, mew=1.4, clip_on=False)
    axc.axvline(obs, color=G.EVT, lw=1.0); axc.set_ylim(-0.09 * ym, ym * 1.05); axc.set_yticks([]); axc.spines['left'].set_visible(False); axc.set_xlabel('Shift in seasonal residual of D'); axc.set_title('Observed shift against a simulated null')
    an = v.SST_anomaly.values
    for mask, col in ((pre, G.PRE), (~pre, G.EVT)):
        f = ts._nw_slope(pd.Series(an[mask]), pd.Series(res[mask])); xx = np.array([an[mask].min(), an[mask].max()]); ic = res[mask].mean() - f['slope'] * an[mask].mean()
        axd.scatter(an[mask], res[mask], s=10, color=col, lw=0, zorder=3); axd.plot(xx, ic + f['slope'] * xx, color=col, lw=1.1, zorder=4)
    axd.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2)), zorder=1); axd.set_xlabel('SST anomaly (°C)'); axd.set_ylabel('Seasonal residual of D'); axd.set_title('Response to SST anomaly')
    Xp = np.column_stack([np.ones(pre.sum()), v.SST[pre]]); bb = np.linalg.lstsq(Xp, y[pre], rcond=None)[0]; u = y[pre] - Xp @ bb
    se = np.sqrt((u ** 2).sum() / (pre.sum() - 2) / ((v.SST[pre] - v.SST[pre].mean()) ** 2).sum()); dT = v.SST[~pre].mean() - v.SST[pre].mean()
    ex, lo, hi = bb[1] * dT, (bb[1] - 1.96 * se) * dT, (bb[1] + 1.96 * se) * dT; ob = y[~pre].mean() - y[pre].mean(); yr = [y[yrs == q].mean() - y[pre].mean() for q in (2014, 2015)]
    # e: the comparison drawn on the data it comes from -- D against SST, the pre-event fit carried to the event's mean SST, and the observed event mean above it
    T = v.SST.values; Tp, Te, Dp = T[pre].mean(), T[~pre].mean(), y[pre].mean(); xx = np.linspace(T.min() - 0.2, T.max() + 0.2, 50)
    axe.fill_between(xx, Dp + (bb[1] - 1.96 * se) * (xx - Tp), Dp + (bb[1] + 1.96 * se) * (xx - Tp), color=G.PRE, alpha=0.22, lw=0, zorder=1); axe.plot(xx, Dp + bb[1] * (xx - Tp), color=G.MUTED, lw=1.0, zorder=2)
    axe.scatter(T[pre], y[pre], s=7, color=G.PRE, alpha=0.75, lw=0, zorder=3); axe.scatter(T[~pre], y[~pre], s=8, color=G.EVT, alpha=0.55, lw=0, zorder=3)
    axe.plot([Te, Te], [Dp + lo, Dp + hi], color=G.MUTED, lw=1.0, zorder=5); axe.plot([Te], [Dp + ex], 'o', mfc='white', mec=G.MUTED, ms=5, mew=1.0, zorder=6); axe.plot([Tp], [Dp], 'o', color=G.MUTED, mec='white', mew=0.6, ms=5.5, zorder=6)
    axe.annotate('', xy=(Te, Dp + ob - 0.035), xytext=(Te, Dp + ex + 0.05), arrowprops=dict(arrowstyle='-|>', color=G.EVT, lw=0.9, shrinkA=0, shrinkB=0, mutation_scale=6), zorder=5); axe.plot([Te], [Dp + ob], 'o', color=G.EVT, mec='white', mew=0.6, ms=6, zorder=7)
    for q, d_ in zip((2014, 2015), yr): axe.plot([T[yrs == q].mean()], [Dp + d_], 'o', mfc='white', mec=G.EVT, ms=3.6, mew=0.8, zorder=6)
    lk = dict(arrowstyle='-', lw=0.45, shrinkA=1, shrinkB=3.5)
    axe.annotate('Observed mean,\n2014–2015', xy=(Te, Dp + ob), xytext=(19.5, 1.2), color=G.EVT, fontsize=G.fs(6), va='center', ha='left', linespacing=1.2, arrowprops=dict(color=G.EVT, **lk), zorder=8)
    axe.annotate('Expected from\npre-event slope', xy=(Te, Dp + ex), xytext=(19.6, -1.22), color=G.MUTED, fontsize=G.fs(6), va='center', ha='left', linespacing=1.2, arrowprops=dict(color=G.MUTED, **lk), zorder=8)
    axe.set_xlim(T.min() - 0.4, T.max() + 0.4); axe.set_ylim(-1.45, 1.42); axe.set_yticks([-1, -0.5, 0, 0.5, 1]); axe.spines['left'].set_bounds(-1, 1); axe.set_xlabel('SST (°C)'); axe.set_ylabel('D'); axe.set_title('Observed vs expected')
    for ax, s in zip((axa, axb, axc, axd, axe), 'abcde'): lab(ax, s, 9 if s != 'c' else 3)
    G.save(fig, OUT + 'Fig2_test')
    src = {'c_null': pd.DataFrame({'null_event_shift': null}), 'c_assignments': pd.DataFrame([{'event_years': f'{p[0]}+{p[1]}', 'shift': s} for p, s in assign.items()]),
           'd_residual_vs_anomaly': pd.DataFrame({'sample_id': v.sample_id, 'date': v.date.dt.date, 'period': np.where(pre, '2009-2013', '2014-2015'), 'SST_anomaly': an, 'D_seasonal_residual': res}),
           'e_expected_vs_observed': pd.DataFrame([{'quantity': 'expected from pre-event D-SST slope', 'dD': ex, 'lo': lo, 'hi': hi}, {'quantity': 'observed', 'dD': ob}, {'quantity': 'observed, 2014 only', 'dD': yr[0]}, {'quantity': 'observed, 2015 only', 'dD': yr[1]}]),
           'e_points_and_fit': pd.concat([pd.DataFrame({'series': np.where(pre, 'sample 2009-2013', 'sample 2014-2015'), 'date': v.date.dt.date, 'SST': T, 'D': y}),
                                          pd.DataFrame([dict(series='mean 2009-2013', SST=Tp, D=Dp), dict(series='expected at event mean SST', SST=Te, D=Dp + ex, D_low=Dp + lo, D_high=Dp + hi), dict(series='observed mean 2014-2015', SST=Te, D=Dp + ob),
                                                        dict(series='observed mean 2014', SST=T[yrs == 2014].mean(), D=Dp + yr[0]), dict(series='observed mean 2015', SST=T[yrs == 2015].mean(), D=Dp + yr[1]), dict(series='pre-event fit', slope=bb[1], slope_se=se, intercept=bb[0])])], ignore_index=True)}
    return src

# ============================================================== Fig. 3
def density(ax, x, y, color, xr, yr, levels=(0.25, 0.5, 0.75, 0.95)):
    """Filled kernel-density background: nested regions holding 25, 50, 75 and 95 % of the points, darker towards the mode."""
    from scipy.stats import gaussian_kde
    import matplotlib.colors as mc
    m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]; k = gaussian_kde(np.vstack([x, y])); gx, gy = np.meshgrid(np.linspace(*xr, 140), np.linspace(*yr, 140)); z = k(np.vstack([gx.ravel(), gy.ravel()])).reshape(gx.shape)
    zs = np.sort(z.ravel())[::-1]; cs = np.cumsum(zs) / zs.sum(); th = sorted(zs[np.searchsorted(cs, l)] for l in levels); rgb = mc.to_rgb(color)
    ax.contourf(gx, gy, z, levels=th + [z.max() * 1.01], colors=[rgb + (a,) for a in (0.07, 0.12, 0.18, 0.25)], zorder=1); ax.contour(gx, gy, z, levels=th[:1], colors=[rgb + (0.35,)], linewidths=0.4, zorder=1)
def fig3():
    G.use(1.30)
    st = rd('Fig3_validate', 'a_b_station_means'); oof = rd('Fig3_validate', 'c_oof_predictions'); tr = rd('Fig3_validate', 'd_model_transfer'); ya = rd('Fig3_validate', 'f_year_assignments'); ph = rd('Fig3_validate', 'g_per_phylotype')
    W = within_order_json(); st['arctic'] = st.basin.str.contains('Arctic'); oof['arctic'] = oof.ocean_basin.str.contains('Arctic')
    fig = G.figure(183, 165)
    axa = fig.add_axes([0.025, 0.628, 0.555, 0.355], projection=ccrs.Robinson()); axa.set_global(); axa.add_feature(cf.LAND, facecolor='#EAEAEA', lw=0); axa.coastlines(lw=0.25, color='#BDBDBD'); axa.spines['geo'].set_linewidth(0.5)
    sc = axa.scatter(st.lon, st.lat, c=st.D, cmap=G.cmap_orders(), vmin=-1, vmax=1, s=36, edgecolor='white', lw=0.5, transform=ccrs.PlateCarree(), zorder=3)
    for _lat in (-35, 35):
        axa.plot(np.linspace(-180, 180, 100), np.full(100, _lat), transform=ccrs.PlateCarree(), color='#9A9A9A', lw=0.4, ls=':', zorder=2)
    cax = fig.add_axes([0.583, 0.700, 0.012, 0.205]); cb = fig.colorbar(sc, cax=cax, orientation='vertical', ticks=[-1, -0.5, 0, 0.5, 1]); cb.outline.set_linewidth(0.4); cax.tick_params(length=2, width=0.5, labelsize=G.fs(6)); cax.set_title('D-index', fontsize=G.fs(6), pad=3)
    fig.text(0.025, 0.983, 'a', fontsize=G.fs(8), fontweight='bold'); fig.text(0.055, 0.983, 'Tara Oceans: station-mean D', fontsize=G.fs(6.6), fontweight='bold')
    axb = fig.add_axes([0.708, 0.665, 0.272, 0.30]); na, ar = st[~st.arctic], st[st.arctic]; axb.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2)))
    axb.scatter(na['T'], na.D, s=11, color=G.INK, lw=0, zorder=3, label='Other oceans'); axb.scatter(ar['T'], ar.D, s=11, facecolor='white', edgecolor=G.INK, lw=0.6, zorder=3, label='Arctic')
    for d_, ls in ((st, '-'), (na, (0, (4, 2)))):
        b = np.polyfit(d_['T'], d_.D, 1); xx = np.array([d_['T'].min(), d_['T'].max()]); axb.plot(xx, np.polyval(b, xx), color=G.MUTED, lw=1.0, ls=ls, zorder=2)
    r_all, p_all = stats.pearsonr(st['T'], st.D); r_na, p_na = stats.pearsonr(na['T'], na.D)
    pa = '$P$ < 0.001' if p_all < 1e-3 else '$P$ = ' + G.sci(p_all); pn = '$P$ < 0.001' if p_na < 1e-3 else '$P$ = ' + G.sci(p_na)
    axb.text(0.03, 0.98, 'All stations  $r$ = %.2f, %s\nNon-Arctic  $r$ = %.2f, %s' % (r_all, pa, r_na, pn),
             transform=axb.transAxes, va='top', ha='left', fontsize=G.fs(5.4), color=G.INK, linespacing=1.5)
    axb.set_xlabel('Temperature (°C)'); axb.set_ylabel('Station-mean D'); axb.set_ylim(-1.12, 1.42); axb.set_yticks([-1, -0.5, 0, 0.5, 1]); axb.legend(loc='lower right', handletextpad=0.2); axb.set_title('D against temperature')
    axc = fig.add_axes([0.075, 0.375, 0.20, 0.19]); axc.plot([-1, 1], [-1, 1], color=G.LIGHT, lw=0.7, ls=(0, (3, 2)))
    for d_, kw in ((oof[~oof.arctic], dict(color=G.INK, lw=0)), (oof[oof.arctic], dict(facecolor='white', edgecolor=G.INK, lw=0.5))): axc.scatter(d_.observed_D, d_.predicted_D, s=6, zorder=3, **kw)
    axc.set_xlabel('Observed D'); axc.set_ylabel('Predicted D'); axc.set_xticks([-1, 0, 1]); axc.set_yticks([-1, 0, 1]); axc.set_title('Out-of-fold prediction')
    xo, yp = oof.observed_D.values, oof.predicted_D.values; r2 = 1 - ((xo - yp) ** 2).sum() / ((xo - xo.mean()) ** 2).sum()
    axc.text(0.97, 0.04, '$R^2$ = %.2f, $n$ = %d' % (r2, len(oof)), transform=axc.transAxes, va='bottom', ha='right', fontsize=G.fs(5.6), color=G.INK)
    axd = fig.add_axes([0.40, 0.375, 0.33, 0.19])
    NAME = {'Temperature +\nfraction (D)': 'Temperature only', 'Ridge, 7 env\n+ fraction': 'Ridge', 'Random forest, 7 env\n+ fraction (log-ratio)': 'Random forest (log-ratio)', 'Random forest, 7 env\n+ fraction (D)': 'Random forest (D)', 'Ensemble (ridge + random forest), 7 env + fraction': 'Ensemble'}
    assert set(tr.model) == set(NAME); HI = {'Ensemble': G.INK}
    for r in tr.itertuples():
        nm = NAME[r.model]; c_ = HI.get(nm, G.MUTED); hi = nm in HI; a_, b_ = r.station_blocked_R2_mean, r.leave_one_basin_out_R2
        if pd.notna(r.station_blocked_R2_min): axd.plot([0, 0], [r.station_blocked_R2_min, r.station_blocked_R2_max], color=c_, lw=2.6 if hi else 2.0, alpha=0.25, solid_capstyle='round', zorder=2)
        axd.plot([0, 1], [a_, b_], color=c_, lw=1.5 if hi else 0.9, zorder=4 if hi else 3); axd.plot([0, 1], [a_, b_], 'o', color=c_, ms=4.2 if hi else 3.2, mec='white', mew=0.5, zorder=5)
        axd.text(1.06, b_, f'{nm}  {b_:.2f}', va='center', ha='left', fontsize=G.fs(5.8), color=c_, fontweight='bold' if hi else 'normal')
    axd.text(-0.06, tr.station_blocked_R2_mean.max(), f'{tr.station_blocked_R2_mean.max():.2f}', va='center', ha='right', fontsize=G.fs(5.8), color=G.INK, fontweight='bold')
    axd.set_xlim(-0.22, 1.05); axd.set_ylim(0.1, 0.7); axd.set_xticks([0, 1]); axd.set_xticklabels(['Stations withheld', 'Whole basin withheld']); axd.set_ylabel('Out-of-sample $R^2$ for D'); axd.set_yticks([0.2, 0.4, 0.6]); axd.spines['bottom'].set_bounds(0, 1); axd.set_title('Transfer across basins')
    axe = fig.add_axes([0.075, 0.065, 0.235, 0.225]); grp = ['within', 'Imitervirales', 'Algavirales']; cols = {'within': G.INK, 'Imitervirales': G.IMI, 'Algavirales': G.ALG}; rj = np.random.default_rng(3)
    for i, g in enumerate(grp):
        d_ = ya[ya['index'] == g]; ev = d_[d_.event_years.astype(str).str.contains('2014') & d_.event_years.astype(str).str.contains('2015')]; ot = d_.drop(ev.index)
        G.raincloud(axe, ot['shift'].values, i, G.MUTED, width=0.36, seed=i, ms=3); axe.scatter([i], ev['shift'], s=30, marker='D', color=cols[g], edgecolor='white', lw=0.5, zorder=7)
    axe.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2))); axe.set_xticks(range(3)); axe.set_xticklabels(['Both orders', 'Imitervirales', 'Algavirales'], rotation=30, ha='right', rotation_mode='anchor'); axe.tick_params(axis='x', labelsize=G.fs(5.4)); axe.set_xlim(-0.5, 2.5); axe.set_ylim(-2.5, 2.3); axe.set_ylabel('Shift in mean thermal affinity (°C)'); axe.set_title('Within-order affinity shift')
    axf = fig.add_axes([0.415, 0.065, 0.270, 0.225])
    for o, c in (('Imitervirales', G.IMI), ('Algavirales', G.ALG)):
        d_ = ph[ph.order == o]; density(axf, d_.tara_affinity.values, d_.event_change_log10.values, c, (-6, 36), (-2.6, 3.0)); axf.scatter(d_.tara_affinity, d_.event_change_log10, s=3 + 60 * np.sqrt(d_.mean_share.clip(0, 0.2)), color=c, alpha=0.35, lw=0, zorder=2)
        q = pd.qcut(d_.tara_affinity, 6, duplicates='drop'); md = d_.groupby(q, observed=True).agg(x=('tara_affinity', 'median'), y=('event_change_log10', 'median')); axf.plot(md.x, md.y, '-s', color=c, ms=3, lw=1.1, zorder=4, label=o)
    axf.axhline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2))); axf.set_xlabel('Thermal affinity (°C)'); axf.set_ylabel('Change in within-order share (log$_{10}$)'); axf.set_ylim(-1.9, 2.3); axf.set_xlim(-3, 31.5); axf.legend(loc='upper left', handlelength=1.6, ncol=2, columnspacing=1.0); axf.set_title('Phylotype response vs affinity')
    axg = fig.add_axes([0.775, 0.065, 0.205, 0.225]); keys = ['temperature', 'log10_chlorophyll', 'abs_latitude', 'nitrate', 'salinity']; lb = ['Temperature', 'Chlorophyll a', 'Latitude', 'Nitrate', 'Salinity']; out = []
    for o, c, off in (('Imitervirales', G.IMI, 0.16), ('Algavirales', G.ALG, -0.16)):
        vv = [abs(W['C2'][o][k]['rho']) for k in keys]; yy = np.arange(len(keys))[::-1] + off; axg.hlines(yy, 0, vv, color=c, lw=0.9); axg.plot(vv, yy, 'o', color=c, ms=3.4)
        out += [dict(order=o, niche_axis=k, spearman_rho=W['C2'][o][k]['rho'], P=W['C2'][o][k]['P'], n=W['C2'][o][k]['n']) for k in keys]
    axg.set_yticks(range(len(keys))); axg.set_yticklabels(lb[::-1]); axg.set_xlim(0, 0.34); axg.set_xlabel('|Spearman $\\rho$|'); axg.set_title('Which niche axis?')
    for ax, s, xm in ((axb, 'b', 9), (axc, 'c', 9), (axd, 'd', 9), (axe, 'e', 9), (axf, 'f', 10), (axg, 'g', 15)): lab(ax, s, xm)
    G.save(fig, OUT + 'Fig3_validate'); return {'g_niche_axes': pd.DataFrame(out)}

# ============================================================== Fig. 4
def grid(d, col):
    p = d.pivot(index='lat', columns='lon', values=col); return p.columns.values, p.index.values, p.values
def fig4():
    G.use(1.30)
    a = rd('Fig4_hosts_and_fronts', 'a_station'); b = rd59('Fig4_hosts_and_fronts', 'b_cross_match'); c = rd('Fig4_hosts_and_fronts', 'c_variance_partition'); d = rd('Fig4_hosts_and_fronts', 'd_map_2deg_support_masked'); e = rd('Fig4_hosts_and_fronts', 'e_rates_support_domain')
    fig = G.figure(183, 132)
    ax = fig.add_axes([0.07, 0.62, 0.24, 0.34]); sc = ax.scatter(a.D_host, a.D, c=a['T'], cmap=G.CMAP_T, vmin=-4, vmax=30, s=13, edgecolor=G.INK, lw=0.3)
    cax = ax.inset_axes([0.06, 0.845, 0.34, 0.036]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal', ticks=[0, 15, 30]); cb.outline.set_linewidth(0.4); cax.tick_params(length=1.5, width=0.4, labelsize=G.fs(5.8), pad=1); cax.set_title('Temperature (°C)', fontsize=G.fs(6), pad=2, loc='left', fontweight='normal')
    ax.set_xlabel('Host balance, D$_{\mathrm{host}}$'); ax.set_ylabel('Station-mean D'); ax.set_xlim(-1.08, 1.08); ax.set_ylim(-1.08, 1.45); ax.set_yticks([-1, -0.5, 0, 0.5, 1]); ax.set_title('D against the host balance'); lab(ax, 'a')
    ax = fig.add_axes([0.470, 0.625, 0.168, 0.325]); ax.axvline(0, color=G.MUTED, lw=0.5, ls=(0, (3, 2)))
    for i, (fam, col, nm) in enumerate((('IM_01', G.IMI, 'Mesomimi-\nviridae'), ('AG_01', G.ALG, 'Prasino-\nviridae'))):
        for hg, kw, off in (('Prymnesiophyceae', dict(mfc='white'), 0.13), ('Mamiellophyceae', dict(), -0.13)):
            r = b[(b.virus_family == fam) & (b.host_group == hg)].rho.iloc[0]; ax.plot([0, r], [1 - i + off] * 2, color=col, lw=0.9); ax.plot(r, 1 - i + off, 'o', color=col, mec=col, ms=4.2, mew=0.8, **kw)
    ax.plot([], [], 'o', mfc='white', mec=G.INK, ms=4, mew=0.8, label='Prymnesiophytes'); ax.plot([], [], 'o', color=G.INK, ms=4, label='Mamiellophyceae'); ax.legend(loc='lower left', bbox_to_anchor=(0.02, 0.02), handletextpad=0.2, borderaxespad=0.1, frameon=True, facecolor='white', edgecolor='none', framealpha=1.0, fontsize=G.fs(5.0))
    ax.set_yticks([1, 0]); ax.set_yticklabels(['Mesomimiviridae', 'Prasinoviridae']); ax.tick_params(axis='y', labelsize=G.fs(6.2)); ax.set_ylim(-1.0, 1.6); ax.set_xlim(-0.8, 0.9); ax.set_xlabel('Spearman $\\rho$, virus vs host share'); ax.set_title('Virus families and their hosts'); lab(ax, 'b', 15)
    ax = fig.add_axes([0.765, 0.62, 0.215, 0.34])
    for i, r in c.iterrows():
        bot = 0
        for k, col, nm in (('unique_temperature', G.N_MID, 'Temperature only'), ('shared', G.N_LIGHT, 'Shared'), ('unique_host', G.HOST, 'Host balance only')):
            ax.bar(i, r[k], 0.55, bottom=bot, color=col, lw=0, label=nm if i == 0 else None); bot += r[k]
    ax.set_xticks([0, 1]); ax.set_xticklabels(['All', 'Non-Arctic']); ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 1.0); ax.set_ylabel('Variance in D explained ($R^2$)'); ax.set_title('Temperature or hosts?'); h, l = ax.get_legend_handles_labels(); ax.legend(h[::-1], l[::-1], loc='upper right'); lab(ax, 'c')
    ax = fig.add_axes([0.03, 0.09, 0.62, 0.43], projection=ccrs.Robinson(central_longitude=-160)); ax.set_global(); ax.spines['geo'].set_linewidth(0.5)
    lon, lat, z = grid(d, 'dD'); pm = ax.pcolormesh(lon, lat, z, cmap=G.cmap_orders(), vmin=-0.4, vmax=0.4, transform=ccrs.PlateCarree(), rasterized=True, shading='nearest')
    _, _, sup = grid(d.assign(s=(~d.support.astype(bool)).astype(float)), 's'); ax.contourf(lon, lat, sup, levels=[0.5, 1.5], colors=['#F2F2F2'], transform=ccrs.PlateCarree())
    dis = d[(~d.all_models_agree.astype(bool)) & d.dD.notna()].iloc[::3]; ax.scatter(dis.lon, dis.lat, s=0.25, color=G.INK, lw=0, transform=ccrs.PlateCarree())
    ax.add_feature(cf.LAND, facecolor='#CFCFCF', lw=0, zorder=4)
    for col_, ls in (('D_2000s', (0, (3, 1.5))), ('D_2090s', '-')):
        _, _, zz = grid(d, col_); ax.contour(lon, lat, zz, levels=[0], colors=[G.INK], linewidths=0.9, linestyles=[ls], transform=ccrs.PlateCarree(), zorder=5)
    ax.plot([], [], color=G.INK, lw=0.9, ls=(0, (3, 1.5)), label='D = 0, 2000s'); ax.plot([], [], color=G.INK, lw=0.9, label='D = 0, 2090s'); fig.legend(*ax.get_legend_handles_labels(), loc='lower left', bbox_to_anchor=(0.03, 0.035), ncol=2)
    cax = fig.add_axes([0.36, 0.062, 0.22, 0.011]); cb = fig.colorbar(pm, cax=cax, orientation='horizontal', ticks=[-0.4, 0, 0.4], extend='both'); cb.outline.set_linewidth(0.4); cax.tick_params(length=2, width=0.5, labelsize=G.fs(6)); cb.set_label('ΔD, 2090s − 2000s', labelpad=1)
    fig.text(0.03, 0.525, 'd', fontsize=G.fs(8), fontweight='bold'); fig.text(0.062, 0.525, 'Projected change in D by the 2090s', fontsize=G.fs(7.2), fontweight='bold')
    ax = fig.add_axes([0.772, 0.115, 0.198, 0.40])
    for j, hm in enumerate(['North', 'South']):
        dd = e[e.hemisphere == hm]; yv = 1 - j; ax.plot([dd.km_per_decade.min(), dd.km_per_decade.max()], [yv, yv], color=G.MUTED, lw=0.7); ax.plot([dd.km_per_decade.mean()] * 2, [yv - 0.16, yv + 0.16], color=G.INK, lw=1.6)
        for _, r in dd.iterrows(): ax.plot(r.km_per_decade, yv, marker=G.ESM[r.ESM]['marker'], mfc='white', mec=G.INK, ms=3.6, mew=0.7, ls='none', label=r.ESM if j == 0 else None)
    ax.set_yticks([1, 0]); ax.set_yticklabels(['Northern front', 'Southern front']); ax.set_ylim(-0.75, 1.5); ax.set_xlim(0, 100); ax.set_xlabel('Poleward rate (km per decade)'); ax.set_title('Front movement'); ax.legend(loc='lower right', handletextpad=0.1, fontsize=G.fs(4.8), borderaxespad=0); lab(ax, 'e', 11)
    G.save(fig, OUT + 'Fig4_hosts_and_fronts')

# ============================================================== Fig. 5
def wilson(k, n, z=1.96):
    p = k / n; den = 1 + z * z / n; c = (p + z * z / (2 * n)) / den; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den; return c - h, c + h
def violin(ax, data, pos, col, w=0.34):
    vp = ax.violinplot(data, positions=[pos], widths=2 * w, showextrema=False); 
    for b in vp['bodies']: b.set_facecolor(col); b.set_alpha(0.22); b.set_edgecolor(col); b.set_linewidth(0.6)
    q1, md, q3 = np.percentile(data, [25, 50, 75]); ax.plot([pos, pos], [q1, q3], color=col, lw=1.6, solid_capstyle='butt'); ax.plot(pos, md, 'o', mfc='white', mec=col, ms=3, mew=0.8)
def fig5():
    G.use(1.16)
    a = rd('Fig5_interpret', 'a_trait_prevalence'); b = rd('Fig5_interpret', 'b_pathway_richness'); c = rd('Fig5_interpret', 'c_microdiversity')
    nI, nA = (b.order == 'Imitervirales').sum(), (b.order == 'Algavirales').sum()
    fig = G.figure(183, 66); gs = fig.add_gridspec(1, 4, left=0.152, right=0.985, top=0.875, bottom=0.17, wspace=0.40, width_ratios=[1.3, 0.5, 0.62, 1.3])
    gt = pd.read_csv(SUP + 'Supplementary_Data_7_GVMAG_functional_traits.csv'); gt = gt[gt.order.isin(['Imitervirales', 'Algavirales'])].copy()
    TC = {'Heat-shock proteins': 'pw_heat_shock_proteins', 'DNA repair': 'pw_DNA_repair', 'Cold-shock proteins': 'pw_cold_shock_proteins', 'Type-I rhodopsin': 'type_I_rhodopsin', 'TCA cycle': 'pw_TCA'}; assert list(a.trait) == list(TC)
    M = (gt[list(TC.values())] > 0).astype(int); M.columns = list(TC); gt['n5'] = M.sum(axis=1).values
    for r in a.itertuples(): assert M.loc[(gt.order == 'Imitervirales').values, r.trait].sum() == r.kI and M.loc[(gt.order == 'Algavirales').values, r.trait].sum() == r.kA     # the matrix reproduces the deposited counts
    ax = fig.add_subplot(gs[0]); x0 = 0; GAP = 9; ticks = []
    import matplotlib.colors as mc
    for o, col, n in (('Imitervirales', G.IMI, nI), ('Algavirales', G.ALG, nA)):
        sel = (gt.order == o).values; od = np.lexsort([-M[sel][t].values for t in list(TC)[::-1]] + [-gt.n5[sel].values]); mm = M[sel].values[od].T; assert mm.shape[1] == n
        ax.imshow(mm, extent=(x0, x0 + n, len(TC), 0), aspect='auto', interpolation='nearest', cmap=mc.ListedColormap([mc.to_rgba(col, 0.10), mc.to_rgba(col, 1.0)]), vmin=0, vmax=1)
        ax.text(x0 + n / 2, -0.18, f'{o} ({n})', ha='center', va='bottom', fontsize=G.fs(6), color=col); ticks.append((x0, n)); x0 += n + GAP
    for yv in range(1, len(TC)): ax.axhline(yv, color='white', lw=1.6)
    ax.set_xlim(0, x0 - GAP); ax.set_ylim(len(TC), -0.02); ax.set_yticks(np.arange(len(TC)) + 0.5); ax.set_yticklabels(list(TC)); ax.set_xticks([]); ax.tick_params(length=0)
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_xlabel('Genomes, one column each'); ax.set_title('Trait repertoire', pad=12); lab(ax, 'a', 21, 0.05)
    axp = fig.add_subplot(gs[1]); bm, bp = ax.get_position(), axp.get_position(); axp.set_position([bm.x1 + 0.012, bp.y0, bp.x1 - bm.x1 - 0.03, bp.height])
    for i, r in a.iterrows():
        yv = i + 0.5; axp.plot([100 * r.kA / nA, 100 * r.kI / nI], [yv, yv], color=G.LIGHT, lw=1.6, zorder=1, solid_capstyle='round')
        for k, n, col in ((r.kI, nI, G.IMI), (r.kA, nA, G.ALG)):
            lo, hi = wilson(k, n); axp.plot([100 * lo, 100 * hi], [yv] * 2, color=col, lw=0.9, zorder=2); axp.plot(100 * k / n, yv, 'o', color=col, ms=3.8, mec='white', mew=0.4, zorder=3)
    axp.set_ylim(len(TC), -0.02); axp.set_yticks([]); axp.spines['left'].set_visible(False); axp.set_xlim(-6, 106); axp.set_xticks([0, 50, 100]); axp.set_xlabel('Genomes (%)'); axp.set_title('Prevalence', pad=12)
    EX5 = {'a_trait_matrix': pd.concat([gt[['genome_id', 'order']].reset_index(drop=True), M.reset_index(drop=True)], axis=1)}
    ax = fig.add_subplot(gs[2])
    for j, (o, col) in enumerate((('Imitervirales', G.IMI), ('Algavirales', G.ALG))): G.raincloud(ax, b[b.order == o].n_pathways_present.values, j, col, width=0.4, seed=j)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Imitervirales', 'Algavirales'], rotation=40, ha='right', rotation_mode='anchor'); ax.tick_params(axis='x', labelsize=G.fs(4.9)); ax.set_ylabel('Pathway categories per genome'); ax.set_xlim(-0.55, 1.55); ax.set_title('Pathway breadth', pad=12); lab(ax, 'b', 9, 0.05)
    ax = fig.add_subplot(gs[3]); fr = list(dict.fromkeys(c.size_fraction))
    for i, f_ in enumerate(fr):
        for o, col, off in (('Imitervirales', G.IMI, -0.25), ('Algavirales', G.ALG, 0.25)):
            dd = c[(c.size_fraction == f_) & (c.order == o)].median_pi.dropna().values
            if len(dd) > 2: G.raincloud(ax, dd, i + off, col, width=0.2, seed=i, ms=1.6)
    ax.set_xticks(range(len(fr))); ax.set_xticklabels(['Viral', 'Pico', 'Micro/nano'][:len(fr)]); ax.set_ylabel('Nucleotide diversity ($\pi$) per genome'); ax.set_ylim(0, None); ax.set_xlim(-0.55, len(fr) - 0.45); ax.set_title('Population microdiversity', pad=12); lab(ax, 'c', 9, 0.05)
    G.save(fig, OUT + 'Fig5_interpret'); return EX5

if __name__ == '__main__':
    import os; os.makedirs(OUT, exist_ok=True); import pickle
    which = sys.argv[1:] or ['1', '2', '3', '4', '5']; extra = {}
    for w in which:
        r = globals()['fig' + w]()
        if r: extra['Fig' + w] = r
    if extra: pickle.dump(extra, open(OUT + f'extra_source_main_{"".join(which)}.pkl', 'wb'))
    print('built', which)
