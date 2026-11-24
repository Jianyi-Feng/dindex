# Manuscript analysis scripts

Code for the analyses reported in *Recurrent warming-associated shifts in marine giant-virus dominance across scales*.
Scripts are grouped by the part of the paper they produce and are listed in the
order they were run.

Every script resolves its input and output paths relative to its own location:

```python
PKG = Path(__file__).resolve().parent.parent     # repository root
RAW = PKG / "raw_tables"                         # raw input tables
PROC = PKG / "processed_tables"                  # analysis and figure outputs
```

The tables these scripts read and write are part of this repository, so nothing has to be downloaded
first. `RUN.md` lists which scripts additionally need the two large gridded CMIP6 fields or external
raw model output.

## Index derivation and the SPOT time series

| Script | What it does |
|---|---|
| `inherited/dindex_analyses.py` | Core SPOT analyses: D by size fraction, temperature sensitivity at SPOT, seasonal climatology, zero points, heatwave residuals, Tara calibration and transfer |
| `inherited/dindex_analyses_extra.py` | Additional SPOT summaries (seasonal statistics, event tests) |
| `pre2014_pair_derivation.py` | Screens all 36 pairs of the nine giant-virus orders among the 89 pre-2014 metagenomes, and refits the D–SST slope on that subset |

## Marine-heatwave inference

Heatwave residuals, the exact year-label permutation over the 21 two-of-seven year assignments,
the calendar-year bootstrap and the same-month comparisons are produced by
`inherited/dindex_analyses.py`. Event tests are tabulated in
`processed_tables/supplementary_data/Supplementary_Data_3_*` and
`processed_tables/supplementary_tables/Supplementary_Table_3_*`.

## Tara Oceans calibration, validation and benchmarks

| Script | What it does |
|---|---|
| `tara_benchmarks.py` | Station-blocked cross-validation and leave-one-basin-out validation for the temperature, ridge and random-forest models, plus out-of-fold permutation importance (seven environmental predictors) |
| `lobo_predictions.py` | Pooled leave-one-basin-out predictions used for the basin-level figure |
| `nested_cv_without_pigment.py` | Nested station-blocked validation of the 14-regressor ensemble after excluding the mislabelled pigment field |
| `tara_carbonate_table_audit.py` | Shows that the field labelled carbonate in the source table is identical to the HPLC median of 19'-butanoyloxyfucoxanthin for all 24,241 overlapping non-missing rows |
| `carbonate_sensitivity.py` | Quantifies the effect of that mislabelled field (audit only; not used for the final maps) |
| `genome_and_validation_tables.py` | Assembles the genome-level and model-validation tables |
| `cmip6_fronts_table.py` | Builds the gene/pathway glossary |

## CMIP6 projections and dominance fronts

| Script | What it does |
|---|---|
| `cmip6_extract.py` | Regrids the source model output onto the common grid (`archived__cmip6_grid.npz`, eight-field version) |
| `cmip6_units_check.py` | Checks that the CMIP6 predictor units match those of the Tara training table |
| `Predict_D-index_noCO3.py` | The final ecological model: 14 candidate regressors and seven environmental predictors, with the mislabelled pigment field excluded |
| `projection_rerun.py` → `cmip6_analyses.py` | Applies that model to historical and SSP5-8.5 fields; produces the current `cmip6_grid.npz` and the front trajectories and rates |
| `cmip6_gridded_analyses.py` | Front detection, front rates, sign agreement and the eight-field sensitivity |
| `cmip6_aoa_fronts.py` | Re-estimates fronts restricted to the environmental-support domain |
| `climate_velocity_comparison.py` | Tracks the sea-surface isotherm on which each model's 2000s front lay, and computes the velocity of climate change on the 2° grid |

## Genomes, microdiversity and association network

| Script | What it does |
|---|---|
| `pnps_and_network.py` | pN/pS contrasts under coverage thresholds and the descriptive 18S–giant-virus association network |
| `assemble_supplementary_data.py` | Assembles Supplementary Data 2–10 from the analysis outputs |
| `historical_new_data_analyses.py` | **Historical — do not rerun.** Provenance source for the earlier-draft supplementary tables; it writes earlier-draft tables and figures |

## Figures and checks

`figures/make_figures.py` (main Figs. 1–5 and Extended Data Figs. 1, 2, 3, 5),
`figures/make_extended_data_fig4.py`, `figures/make_supplementary_figs4_5_pnps_network.py`,
`figures/make_supplementary_fig3_components.py` and `figures/build_supplementary_figures.py`
regenerate the published panels and their per-figure source data into `outputs/`.

`figures/v73/` holds the V73 production scripts used for the submitted figures:
`build_main.py` (main Figs. 1–5), `build_ed.py` (Extended Data Figs. 1–9), `gvstyle.py` (house style)
and `qa_figures.py` (Nature-spec PDF check). `figures/v73/README.md` maps every figure to its
source-data sheets and lists the optional author-local inputs and their deposited fallbacks.

Pre-submission consistency checks were run in the author-side submission tree; they are not part of
this repository.

Manuscript text and the Word-file builders are not part of this repository.

## Randomness

Random seeds are set inside each script (`random_state=1` for forests, `seed=42` for the inner
randomised search, `seed=0` for bootstraps). Software and database versions and the critical
parameters are listed in Supplementary Table 11 of the manuscript.
