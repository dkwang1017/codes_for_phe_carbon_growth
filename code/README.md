# Code for: The carbon assimilation–tree growth relationship is stronger in anomalously early springs but weaker where springs advance faster

This repository contains the analysis code used for the tree-ring data
processing, the model-based phenology extraction, and the baseline GPP–RWI
mixed-effects model reported in the manuscript. Methodological details for each
step are given in the **Methods** section of the paper.

| Script | Language | What it does |
|---|---|---|
| `codes/tree_ring_data_processing.R` | R | Detrends raw ITRDB ring-width series with a 30-year cubic smoothing spline (50% frequency response), builds site chronologies with a biweight robust mean, computes the expressed population signal (EPS), and exports chronologies for 1982–2018 with sample depth ≥ 5. |
| `codes/cmip6_trendy_gpp_phenology_extraction.R` | R | Derives annual SOS and EOS from monthly CMIP6/TRENDY GPP fields: gap filling by cubic spline interpolation, seasonal curve fitting with the Elmore double-logistic function in `phenofit`, and extraction at 50% of the seasonal amplitude. |
| `codes/quantify_gpp_rwi_relationship.py` | Python (calls R) | Fits the baseline mixed-effects model of ring-width index on current- and previous-year GPP, with categorical year effects, a GPP-grid random intercept, and chronology-specific random intercept and slopes. Writes fixed effects, variance components, and convergence diagnostics. |

## 1. System requirements

- **Operating system:** Windows 10/11 (developed on Windows 10, 64-bit). The code uses only cross-platform R and Python packages and should run on Linux or macOS without modification.
- **R:** version 4.2.3, with the packages below (versions on the analysis workstation):

  | Package | Version |
  |---|---|
  | dplR | 1.7.8 |
  | phenofit | 0.3.9 |
  | R.matlab | 3.7.0 |
  | doParallel | 1.0.17 |
  | foreach | 1.5.2 |
  | imputeTS | 3.3 |
  | terra | 1.7-18 |
  | lubridate | 1.9.3 |
  | lme4 | 1.1-32 |

- **Python:** version 3.9, with:

  | Package | Version |
  |---|---|
  | numpy | 1.23.5 |
  | pandas | 1.4.4 |

- **Rscript** must be on the system `PATH`: `quantify_gpp_rwi_relationship.py` runs the mixed model through R and `lme4`.
- No non-standard hardware is required. Gridded phenology extraction is the most memory-intensive step; a machine with ≥ 16 GB RAM is recommended for 0.5° global arrays.

## 2. Installation

R packages:

```r
install.packages(c("dplR", "phenofit", "R.matlab", "doParallel",
                   "foreach", "imputeTS", "terra", "lubridate", "lme4"))
```

Python packages:

```bash
pip install numpy pandas
```

## 3. Quick start (demo)

No demo dataset is bundled. To check that the environment is set up correctly,
run the Python script in audit-only mode on any prepared input table; it checks
the input columns and reports the sample-selection counts without fitting:

```bash
python codes/quantify_gpp_rwi_relationship.py \
    --input annual_analysis_table.csv \
    --product GLASS \
    --output-dir outputs/glass \
    --prepare-only
```

Expected output: `outputs/glass/M0_sample_membership.csv` and
`outputs/glass/M0_sample_audit.json`, the latter reporting the number of input
rows, complete-case rows, retained chronologies, and GPP grid cells.

## 4. Instructions for use

### 4.1 Tree-ring processing

Edit the two paths at the top of `tree_ring_data_processing.R`:

```r
input_dir  <- "path/to/ITRDB/treering/measurements"   # .rwl files, searched recursively
output_dir <- "path/to/output/tree_ring_chronologies"
```

then run:

```bash
Rscript codes/tree_ring_data_processing.R
```

Outputs, written to `output_dir`:

- `<site>.crn` — one chronology file per site
- `EPS_global.csv` — EPS for each site
- `Global_RWI_full.csv` — all chronologies, all years
- `Global_RWI_1982_2018_sampdepth_ge_5.csv` — filtered chronologies used downstream
- `processing_log_global.csv`, `session_info_tree_ring_processing.txt`

### 4.2 CMIP6 and TRENDY phenology extraction

Edit `input_dir` and `output_dir` at the top of
`cmip6_trendy_gpp_phenology_extraction.R`. The script expects one `.mat` file
per model containing a monthly GPP array (latitude × longitude × month) covering
1980–2020. Edit `model_variable_map` to match the variable names inside the
files; any CMIP6 model added there must also be listed in `cmip6_models` so that
the shorter 1982–2014 window is applied.

The analysis windows follow the manuscript: **TRENDY 1982–2018, CMIP6
1982–2014**. The provided example extracts the **50% amplitude threshold**
(the `stop()` guard rejects other thresholds).

```bash
Rscript codes/cmip6_trendy_gpp_phenology_extraction.R
```

Outputs, written to `output_dir`: `CMIP6_TRENDY_GPP_<model>_SOS_<y1>_<y2>_TRS50.tif`,
the matching `_EOS_` file, and `session_info_CMIP6_TRENDY_GPP_phenology.txt`.

### 4.3 Baseline GPP–RWI mixed model

The script starts from an annual, preprocessed chronology–year table. Required
columns: `site`, `Year`, `Lat`, `Lon`, `grid_key`, `RWI_detrended_z`, plus the
two product-specific GPP columns (`GPP_current_iav_z` and `GPP_previous_iav_z`
for GLASS; `BESS_GPP_current_iav_z` and `BESS_GPP_previous_iav_z` for BESS).
`GPP_previous_iav_z` must refer to the actual calendar year *t* − 1.

```bash
python codes/quantify_gpp_rwi_relationship.py \
    --input annual_analysis_table.csv \
    --product GLASS \
    --output-dir outputs/glass
```

Add `--rscript /path/to/Rscript` if R is not on the `PATH`, and
`--prepare-only` to export the sample audit without fitting.

The fitted model is

```r
RWI_detrended_z ~ GPP_current_M0_z + GPP_previous_M0_z + Year +
    (1 + GPP_current_M0_z + GPP_previous_M0_z | site) + (1 | grid_key)
```

fitted by maximum likelihood, with sum-to-zero coding of the year factor and
95% confidence intervals and two-sided P values based on an asymptotic normal
(Wald) reference.

Outputs, written to `output_dir`:

- `M0_key_fixed_effects.csv` — current- and previous-year GPP coefficients
- `M0_fixed_effects_all.csv` — all fixed effects, including year terms
- `M0_chronology_random_effects.csv`, `M0_chronology_random_effects_covariance.csv`, `M0_grid_random_intercepts.csv`
- `M0_model_information.json` — likelihood, R², singular-fit flag, optimizer code, R/lme4 versions
- `M0_sample_membership.csv`, `M0_sample_audit.json` — rows retained and counts

## 5. Notes

- `tree_ring_data_processing.R` exports EPS values but does not itself apply the
  EPS ≥ 0.85 screen; the chronology inclusion criteria used in the manuscript
  are applied when the analysis sample is assembled.
- Raster extents in `cmip6_trendy_gpp_phenology_extraction.R` are assigned as
  −180°–180° E, 30°–90° N, i.e. the input arrays are assumed to be pre-cropped
  to the study domain north of 30° N.

## 6. License and repository

The code is released under the **MIT License**.

Source repository: <https://github.com/dkwang1017/codes_for_phe_carbon_growth>

A complete description of each step is provided in the **Methods** section of
the manuscript (subsections *Tree-ring data*, *Relationship between carbon
assimilation and tree growth*, and *Earth system model simulations*).
