# V73 figure build scripts (main Figs. 1–5 and Extended Data Figs. 1–9)

These are the production scripts for the submitted figures: the V73 build that renders every main figure
and every Extended Data figure to the Nature production spec (183 or 89 mm wide, Arial 5–7 pt, panel labels
8 pt bold, 600 dpi PDF). They are the scripts that were used for the V73 submission package; the same four
files also exist in the author-side submission tree, unchanged apart from the path handling and optional-input
fallbacks listed under [Changes from the author-side scripts](#changes-from-the-author-side-scripts).

## Scripts

| File | Purpose |
|---|---|
| `build_main.py` | Main Figs. 1–5 (`fig1()` … `fig5()`); writes `<figure>.pdf` / `<figure>.png` and a `extra_source_main_*.pkl` with the regenerated per-figure source data |
| `build_ed.py` | Extended Data Figs. 1–9 (`ed1()` … `ed9()`); writes the PDF/PNG pair and `extra_source_ed_*.pkl` |
| `gvstyle.py` | The house style: colours (Imitervirales blue `#0072B2`, Algavirales vermilion `#D55E00`, event purple, host green), figure sizes, rcParams, `save()` |
| `qa_figures.py` | Objective check of every PDF in a folder: width 89/183 mm, height ≤ 170 mm, text 5–7 pt, panel labels 8 pt, Arial |

## Running

```bash
# from the repository root
python figures/v73/build_main.py            # all main figures
python figures/v73/build_main.py 2 4        # only Figs. 2 and 4
python figures/v73/build_ed.py              # all Extended Data figures
python figures/v73/build_ed.py 7            # only Extended Data Fig. 7
python figures/v73/qa_figures.py            # check the PDFs written to outputs/v73
```

Python ≥ 3.9 with `numpy`, `pandas`, `scipy`, `matplotlib` and `cartopy` (map panels); `openpyxl` is needed to
read the deposited `.xlsx` source data, and `pymupdf` for `qa_figures.py`. Arial must be installed for the
intended typeface (see `figures/README.md` for font installation); without it the scripts fall back to
DejaVu Sans and `qa_figures.py` reports a non-Arial font.

The scripts write to `outputs/v73/` at the repository root (`outputs/` is git-ignored). Input paths and output
paths are resolved from the script location, so they can be run from any working directory.

## Inputs

All figure inputs are the deposited source-data workbooks in `processed_tables/figure_source_data/`
(`<workbook>_source_data.xlsx`, one sheet per panel or panel group), the deposited Supplementary Data CSVs in
`processed_tables/supplementary_data/`, and two original observation tables in `raw_tables/`. The worksheet
name of every panel is visible in the scripts, e.g. `rd('Fig3_validate', 'a_b_station_means')`.

| Environment variable | Default | Meaning |
|---|---|---|
| `DINDEX_INPUTS` | `<repo>/processed_tables` | Root of the deposited tables (`figure_source_data/`, `supplementary_data/`, `supplementary_tables/`) |
| `DINDEX_OUT` | `<repo>/outputs/v73` | Output folder for the PDF/PNG figures and the `.pkl` source-data files |
| `DINDEX_STUDY_JSON` | unset | Optional author-local analysis folder (see below); when set, Fig. 3g and Extended Data Fig. 2a read their JSON files instead of the deposited fallbacks |
| `DINDEX_DATA_EXTRA` | unset | Optional author-local `data_extra` folder with `permutation_importance_replicates.csv`; when unset, Extended Data Fig. 3e reads the identical deposited sheet |
| `DINDEX_POLB_CTI` | unset | Optional author-local `polb_cti` pipeline (Newey–West slope used to draw the Fig. 2d regression line) |

The scripts require `raw_tables/` next to `processed_tables/` because Fig. 2 is drawn from the original
sample-level table (`raw_tables/SPOT_sample_Dindex.csv`, 167 samples, identical to Supplementary Data 2 except
for date-string formatting and a few ≤1e-9 rounding differences). If `DINDEX_INPUTS` points elsewhere, keep
`raw_tables/` in the repository or copy the sample table next to the alternative inputs.

## Figure → script → inputs

### Main figures (`build_main.py`)

| Figure | Function | Deposited source data | Additional inputs |
|---|---|---|---|
| Fig. 1 | `fig1()` | `Fig1_derive` (`a_samples`, `b_viral_monthly`, `c_pico_monthly`, `d_micro_monthly`, `f_monthly_D`) | — |
| Fig. 2 | `fig2()` | none (figures computed from the sample table) | `raw_tables/SPOT_sample_Dindex.csv` |
| Fig. 3 | `fig3()` | `Fig3_validate` (`a_b_station_means`, `c_oof_predictions`, `d_model_transfer`, `f_year_assignments` → `e_year_assignments`, `g_per_phylotype` → `f_per_phylotype`); panel g also reads `Fig3_validate::g_niche_axes` when `DINDEX_STUDY_JSON` is unset | optional `03_within_order.json` |
| Fig. 4 | `fig4()` | `Fig4_hosts_and_fronts` (`a_station`, `b_cross_match`, `c_variance_partition`, `d_map_2deg_support_masked`, `e_rates_support_domain`) | — |
| Fig. 5 | `fig5()` | `Fig5_interpret` (`a_trait_prevalence`, `b_pathway_richness`, `c_microdiversity`) | `supplementary_data/Supplementary_Data_7_GVMAG_functional_traits.csv` |

### Extended Data figures (`build_ed.py`)

| Figure | Function | Deposited source data | Additional inputs |
|---|---|---|---|
| ED Fig. 1 | `ed1()` | `ExtendedDataFig1_index_properties` (`a_D_vs_logratio`, `b_combined_share`, `c_slopes`, `d_zero_points`) | `Fig1_derive::a_samples`, `supplementary_data/Supplementary_Data_4_Tara_Dindex_predictors.csv` |
| ED Fig. 2 | `ed2()` | `ExtendedDataFig2_heatwave` (`b_permutations`, `c_same_month`; panel a via `a_between_year_sensitivity` when `DINDEX_STUDY_JSON` is unset) | fallback s.d. from `supplementary_tables/Supplementary_Table_3_heatwave_tests.xlsx` (`model_based_null`); optional `02b_between_year_sensitivity.json` |
| ED Fig. 3 | `ed3()` | `ExtendedDataFig3_global_robustness` (`a_station_D_chl`, `b_lobo`, `c_oof_residuals_by_basin`, `d_transfer`, plus `e_permutation_importance` from `Fig3_validate` and `e_permutation_replicates` when `DINDEX_DATA_EXTRA` is unset) | optional `data_extra/permutation_importance_replicates.csv` |
| ED Fig. 4 | `ed4()` | `ExtendedDataFig4_within_order_heatwave` (`a_residuals`, `b_affinity_shuffle_null`, `c_families`) | — |
| ED Fig. 5 | `ed5()` | `ExtendedDataFig5_host_balance` (`a_station`, `b_variance_partition`, `c_pair_rank`) | — |
| ED Fig. 6 | `ed6()` | `ExtendedDataFig6_SPOT_host_series` (`a_mamiellophyceae_series`, `b_climatology`, `c_single_year_shifts`, `d_e_deseasonalized`) | — |
| ED Fig. 7 | `ed7()` | `ExtendedDataFig7_projection` (`a_b_c_maps_2deg_4model_mean`, `d_front_trajectories`, `e_rates_by_projection`, `f_agreement_by_band`) | — |
| ED Fig. 8 | `ed8()` | `ExtendedDataFig8_projection_robustness` (`a_b_aoa_share_2deg`, `c_aoa_shares`, `d_units`, `e_temperature_only_fronts`, `f_profiles`, `g_front_vs_isotherm`, `h_isotherm_trajectories`); one panel also reads `ExtendedDataFig7_projection::d_front_trajectories` | — |
| ED Fig. 9 | `ed9()` | `ExtendedDataFig9_marker_recovery` (`a_scores`, `b_mantel_haenszel`, `c_pathways_by_score`) | `supplementary_data/Supplementary_Data_7_GVMAG_functional_traits.csv` |

The deposited `processed_tables/figure_source_data/` workbooks are the authoritative per-panel source data
shipped with the manuscript; these scripts regenerate the figures from them, and their `.pkl` outputs hold the
additional per-figure tables that the original author-side build used to assemble the deposited workbooks.

## Changes from the author-side scripts

The numerical plotting code (panel layout, colours, statistics, labels, sizes) is unchanged. Only the
following path/initialisation edits were made so that the scripts run from the repository:

1. **Path constants.** The author-side header
   `SD = '../submission/SPOT_manuscript_V57/03_source_data/'; PK = '../submission/SPOT_manuscript_V72/'; SD59 = PK + '03_source_data/'; OUT = 'out_current/'`
   is replaced by the repository-relative block above (`PKG`, `INPUTS`, `OUT`, `SD`/`SD59`, `SUP`, `RAW`).
   The original job read some sheets from the V57 working tree and the ridge + random-forest ensemble sheets
   from the V72 tree; the deposited V72 package contains **every** sheet, so both reads now resolve to
   `processed_tables/figure_source_data/` and `rd` / `rd59` are the same function.
2. **Three sheet renames** introduced when the V72 deposit was assembled are mapped in `SHEET_ALIAS`:
   `Fig3_validate::e_permutation_importance` → `ExtendedDataFig3_global_robustness::e_permutation_importance`;
   `Fig3_validate::f_year_assignments` → `Fig3_validate::e_year_assignments`;
   `Fig3_validate::g_per_phylotype` → `Fig3_validate::f_per_phylotype`. The values were verified to be
   identical (same rows and numbers).
3. **`polb_cti` import (Fig. 2d).** The original job imported `polb_cti.timeseries` from a sibling pipeline.
   The pipeline is not part of the deposit, so the import is optional (`DINDEX_POLB_CTI`); without it an
   ordinary least-squares `_nw_slope` fallback is used. Fig. 2d draws the fitted regression line, and the
   Newey–West estimator changes the standard error, not the point estimate, so the plotted line is identical.
4. **`strengthening_2026-11-23` JSONs.** Fig. 3g (`03_within_order.json`) and ED Fig. 2a
   (`02b_between_year_sensitivity.json`) read those files when `DINDEX_STUDY_JSON` is set and otherwise rebuild
   the same numbers from the deposited sheets (`Fig3_validate::g_niche_axes`; `ExtendedDataFig2_heatwave::
   a_between_year_sensitivity` plus the s.d. of the five pre-event year means recorded in Supplementary
   Table 3). The rebuilds were verified against the deposit (identical values; the ED2 vertical line is 0.028).
5. **`data_extra/permutation_importance_replicates.csv`.** ED Fig. 3e reads it when `DINDEX_DATA_EXTRA` is set;
   otherwise it uses the deposited `ExtendedDataFig3_global_robustness::e_permutation_replicates` sheet, whose
   per-feature counts (1000) and means reproduce the deposited `e_permutation_importance` values to <1e-9,
   satisfying the script's internal assertion.
6. **Output folder creation.** Both scripts now call `os.makedirs(OUT, exist_ok=True)` before writing (the
   author-side `build_ed.py` assumed the folder already existed).
7. **Arial discovery (`gvstyle.py`).** Font files are now also searched in the macOS and Windows font
   directories and in `~/.fonts`, in addition to the Linux `msttcorefonts` path.
8. **`qa_figures.py`** takes the folder argument optionally and defaults to `outputs/v73`.

## Verification status in the package

`build_main.py`, `build_ed.py`, `gvstyle.py` and `qa_figures.py` compile with the repository's packaging
interpreter, and the three fallback paths (sheet aliases, niche-axis coefficients, between-year sensitivity
grid + pre-event s.d., permutation replicates) were checked numerically against the deposited tables. A full
re-render of all fourteen figures needs `cartopy` and the Arial font on the render host and is the
responsibility of the person reproducing the figures; the deposited PDFs are part of the submission and are
not duplicated in this repository.
