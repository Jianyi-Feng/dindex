# Reproducing the manuscript figures, tables and data

This document follows the chain from the original observations to the published figures:

```
public archive (SRA / JGI / figshare / ENA / Zenodo / ESGF)
        │  raw reads, genomes, model fields, amplicon tables
        ▼
repository: raw_tables/                      original observation and analysis tables
        │  analysis/*.py                      D-index, calibration, heatwave, CMIP6, genome, network scripts
        ▼
repository: processed_tables/                deposited outputs (authoritative)
        ├── supplementary_data/              Supplementary Data 1–14
        ├── supplementary_tables/            Supplementary Tables 1–20
        ├── figure_source_data/              one workbook per main / Extended Data figure
        ├── model_and_statistics/            model fits, robustness and audit tables
        ├── genomes_and_microdiversity/      genome, marker-recovery and inStrain tables
        ├── phylogenies/                     tree files and alignment metadata
        └── cmip6_gridded_fields/            2° CMIP6 gridded fields (three .npz arrays)
        │  figures/v73/*.py                  production figure scripts (main Figs. 1–5, Extended Data Figs. 1–9)
        ▼
outputs/v73/  (PDF + PNG figures, regenerated per-figure source data; git-ignored)
```

Nothing in this chain needs to be downloaded to inspect the tables or the figure code: the repository ships
`raw_tables/` and `processed_tables/`. Downloads are only needed to rebuild the derived tables from the
original reads/models, or to run the map and projection scripts that read the three large gridded CMIP6 arrays.

## Data sources and what each provides

| Archive | Identifier | What it provides | Used for |
|---|---|---|---|
| NCBI Sequence Read Archive | [PRJNA1536449](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1536449) | Newly generated SPOT metagenomic reads (2014–2019 samples) | PolB abundances, sample-level D, order composition, microdiversity |
| JGI Genome Portal | proposal ID `2799` | Previously released SPOT viral metagenomes (the 2009–2013 samples) | the pre-event part of the SPOT time series |
| NCBI BioProject | [PRJNA814250](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA814250) | Previously released 0.2–1.2 µm size-fraction samples | cross-checking the size-fraction response |
| figshare | [10.6084/m9.figshare.34419060](https://doi.org/10.6084/m9.figshare.34419060) | NCLDV and *Mirusviricota* MAGs (nucleotide/protein), PolB abundance tables, sample-level D, genome trait and marker-recovery tables, microdiversity, out-of-fold predictions, gridded CMIP6 fields and predictions, phylogenies | the complete deposit behind `processed_tables/` (including the three `.npz` CMIP6 arrays) |
| ENA | `PRJEB1787`, `PRJEB1788`, `PRJEB9740` | Tara Oceans metagenomes | Tara calibration and validation |
| Zenodo | [10.5281/zenodo.3768510](https://doi.org/10.5281/zenodo.3768510) | Tara Oceans 18S V9 OTU table | 18S host series (Extended Data Fig. 6, Supplementary Data 13/13b) |
| ESGF | [esgf-node.llnl.gov](https://esgf-node.llnl.gov/) | CMIP6 model output for CESM2-WACCM, IPSL-CM6A-LR, NorESM2-MM, UKESM1-0-LL (historical + SSP5-8.5) | projection and dominance-front analyses (Extended Data Figs. 7–8, Supplementary Data 9/9b) |
| GitHub release assets | see `README.md` "Data" | `cmip6_grid.npz`, `cmip6_maps_2deg.npz`, `cmip6_aoa_2deg.npz` (too large for git) | `analysis/cmip6_*.py`, map panels of the figure scripts |

The three `.npz` arrays are also available in the figshare deposit under `07_CMIP6_gridded_fields/`; put them
in `processed_tables/cmip6_gridded_fields/` and `tools/check_data_layout.py` will confirm the layout.

## From raw tables to analysis scripts

`raw_tables/` holds the original observation tables and the intermediate tables the analyses start from.
The main reads are:

| Raw table | Read by |
|---|---|
| `SPOT_sample_Dindex.csv` | `analysis/inherited/dindex_analyses.py`, `analysis/inherited/dindex_analyses_extra.py`, `analysis/historical_new_data_analyses.py`, `analysis/pre2014_pair_derivation.py`, `analysis/assemble_supplementary_data.py`, `figures/make_figures.py`, `figures/v73/build_main.py` (Fig. 2) |
| `SPOT_PolB_order_relative_abundance_long.csv` | `analysis/historical_new_data_analyses.py`, `analysis/pre2014_pair_derivation.py` |
| `SPOT_heatwave_Dindex_stats.csv` | `analysis/inherited/dindex_analyses.py`, `analysis/assemble_supplementary_data.py`, `figures/make_figures.py` |
| `SPOT_seasonal_Dindex_stats.csv` | seasonal baselines and event tests |
| `Tara_observed_Dindex.csv` | `analysis/inherited/dindex_analyses.py`, `analysis/tara_benchmarks.py`, `analysis/lobo_predictions.py`, `analysis/nested_cv_without_pigment.py`, `analysis/cmip6_gridded_analyses.py`, `analysis/cmip6_units_check.py`, `analysis/carbonate_sensitivity.py`, `analysis/tara_carbonate_table_audit.py`, `figures/make_figures.py` |
| `CMIP6_Dindex_latitude_year.csv` | `analysis/cmip6_gridded_analyses.py`, `analysis/inherited/dindex_analyses.py`, `figures/make_figures.py`, `figures/make_extended_data_fig4.py` |
| `Dindex_front_trajectory.csv`, `Dindex_front_regression.csv` | `analysis/cmip6_gridded_analyses.py`, `figures/make_figures.py` |
| `Dindex_model_validation.csv`, `Dindex_permutation_importance.csv` | model-validation and permutation-importance tables/panels (`analysis/inherited/dindex_analyses.py`, `figures/make_figures.py`) |
| `SPOT_18S_GV_network_nodes.csv`, `SPOT_18S_GV_network_edges.csv` | `analysis/pnps_and_network.py`, `figures/make_supplementary_figs4_5_pnps_network.py` |
| `SPOT_GV_and_ref_rhodopsin_metadata.xlsx` | `analysis/historical_new_data_analyses.py` (rhodopsin metadata table) |

Every analysis script resolves its paths from its own location (`PKG = Path(__file__).resolve().parent.parent`;
see `analysis/README.md`), so the scripts run from a fresh clone. `RUN.md` lists which scripts additionally
need the gridded CMIP6 fields.

## From analysis scripts to processed tables

The deposited `processed_tables/` are the authoritative outputs; the scripts that regenerate them are:

| Processed output | Script(s) |
|---|---|
| `processed_tables/supplementary_data/` (Supplementary Data 1–14) | `analysis/assemble_supplementary_data.py` plus the analysis scripts that create the underlying tables (`tara_benchmarks.py`, `lobo_predictions.py`, `cmip6_analyses.py`, `genome_and_validation_tables.py`, the host-series and microdiversity outputs). Supplementary Data 14 (curated gene/pathway list) is a curated input table, not a computed output |
| `processed_tables/supplementary_tables/` (Supplementary Tables 1–20) | produced by the analysis scripts listed in `analysis/README.md`; the gene/pathway glossary is built by `analysis/cmip6_fronts_table.py`, the association network and pN/pS tables by `analysis/pnps_and_network.py`, and the genome/model-validation tables by `analysis/genome_and_validation_tables.py` |
| `processed_tables/figure_source_data/` (one workbook per figure) | the figure scripts (`figures/v73/build_main.py`, `figures/v73/build_ed.py`) write the per-panel tables they draw as `.pkl` files into `outputs/v73/`; the deposited workbooks are the authoritative copies of those tables |
| `processed_tables/model_and_statistics/` | `analysis/inherited/dindex_analyses.py`, `analysis/tara_benchmarks.py`, `analysis/lobo_predictions.py`, `analysis/nested_cv_without_pigment.py`, the CMIP6 scripts and the audit scripts (`carbonate_sensitivity.py`, `tara_carbonate_table_audit.py`, `cmip6_units_check.py`) |
| `processed_tables/genomes_and_microdiversity/` | `analysis/genome_and_validation_tables.py`, `analysis/pnps_and_network.py` and the inStrain-derived genome tables (from the MAGs in the figshare deposit) |
| `processed_tables/phylogenies/` | tree files and metadata imported from the phylogenetic analysis of the NCLDV/Mirusviricota MAGs (raster trees are shown in Supplementary Figs. 1–2) |
| `processed_tables/cmip6_gridded_fields/` | `analysis/cmip6_extract.py` (regridding), `analysis/projection_rerun.py` → `analysis/cmip6_analyses.py` (historical + SSP5-8.5 application of the calibrated model) |

## Figure generation (main Figs. 1–5, Extended Data Figs. 1–9)

The submitted figures are built by the scripts in `figures/v73/`:

```bash
python figures/v73/build_main.py     # main Figs. 1–5  -> outputs/v73/*.pdf, *.png
python figures/v73/build_ed.py       # Extended Data Figs. 1–9
python figures/v73/qa_figures.py     # Nature-spec check of the PDFs
```

The per-figure mapping from drawing function to source-data workbook, and the optional author-local inputs with
their exact deposited fallbacks, are listed in `figures/v73/README.md`. In short:

| Figure | Drawing code | Source-data workbook |
|---|---|---|
| Fig. 1 | `build_main.py:fig1()` | `figure_source_data/Fig1_derive_source_data.xlsx` |
| Fig. 2 | `build_main.py:fig2()` | `raw_tables/SPOT_sample_Dindex.csv` (the deposited `Fig2_test` workbook is the rendered result table) |
| Fig. 3 | `build_main.py:fig3()` | `figure_source_data/Fig3_validate_source_data.xlsx` |
| Fig. 4 | `build_main.py:fig4()` | `figure_source_data/Fig4_hosts_and_fronts_source_data.xlsx` |
| Fig. 5 | `build_main.py:fig5()` | `figure_source_data/Fig5_interpret_source_data.xlsx` + `supplementary_data/Supplementary_Data_7_*` |
| Extended Data Figs. 1–9 | `build_ed.py:ed1()` … `ed9()` | `figure_source_data/ExtendedDataFig1…9_*_source_data.xlsx` |

The earlier scripts in `figures/` (`make_figures.py`, `make_extended_data_fig4.py`,
`make_supplementary_fig3_components.py`, `make_supplementary_figs4_5_pnps_network.py`,
`build_supplementary_figures.py`) render the same panels through the `dindex` package API and assemble the
Supplementary Figures and their per-figure source data. `figures/README.md` documents their render order, the
matplotlib version used for the published PDFs, and the Arial requirement.

Supplementary Figs. 1 and 2 are raster exports of the iTOL phylogeny and the gene heatmap; they are generated
outside these scripts and carried through unchanged.

## Environment and randomness

- Python ≥ 3.9; `requirements.txt` pins the version used for the published PDFs (matplotlib 3.11.2 on Linux;
  the figure scripts add `cartopy` for map panels and `pymupdf` for `qa_figures.py`).
- R (base R only) for `R/dindex.R` and `tests/test_dindex.R`.
- Arial must be installed on the render host (see `figures/README.md`); it is a licensed font and is not
  shipped in this repository.
- Seeds are fixed inside the scripts: `random_state=1` for the random forests, `seed=42` for the inner
  randomised search, `seed=0` for the bootstraps, `20261123` for the Fig. 2 null simulation.
- Software and database versions are listed in Supplementary Table 11.

## Checks

```bat
python tools\check_data_layout.py     :: every table referenced by analysis/ and figures/ is present
set PYTHONPATH=src && python -m pytest -q
Rscript tests/test_dindex.R
```

`tools/check_data_layout.py` resolves the legacy working-tree paths through `tools/paths.py` and is the fastest
way to tell whether a copied tree is complete. `MANIFEST_sha256.csv` at the repository root records the size and
SHA-256 of every tracked file, so a download can be checked row by row; `release/SHA256SUMS.txt` does the same
for the reviewer release archive.
