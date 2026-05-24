# Phenology–Growth Coupling Analysis Code

This folder contains analysis scripts supporting the manuscript:

**Advancing spring phenology weakens the coupling between carbon assimilation and tree growth across the Northern Hemisphere**

The scripts are provided to help editors and reviewers inspect the main analytical workflow and review selected processing steps when the required input datasets are available. A small set of processed data files is provided in the adjacent `../source_data/` folder to illustrate representative input formats.

## Folder structure

Recommended package structure:

```text
codes_for_phe_carbon_growth/
├── LICENSE
├── codes/
│   ├── README.md
│   ├── tree_ring_data_processing.R
│   ├── rwi_gpp_lag_correlation.m
│   ├── gpp_growth_msd_lag_analysis.m
│   ├── cmip6_trendy_gpp_phenology_extraction.R
│   ├── sem_analysis.R
│   └── xgboost_shap_analysis.py
└── source_data/
    ├── README_source_data.md
    ├── ITRDB_v8.0_1982_2018_10years_30N.csv
    └── sos_trend_cmip6_trendy_datasets.zip
```

## Main analysis scripts

| Script | Main purpose | Main outputs |
|---|---|---|
| `tree_ring_data_processing.R` | Processes raw ITRDB tree-ring width files, detrends individual series, constructs site-level ring-width index chronologies, calculates EPS, and exports chronology tables. | RWI chronology tables, EPS table, processing log, and session information. |
| `rwi_gpp_lag_correlation_and_mapping.m` | Calculates site-level lagged Pearson correlations between RWI and 0–12-month lagged accumulated GPP, estimates SOS trends at tree-ring sites, and optionally maps site-level results to a 0.05° grid. | Site-level RWI-GPP lag-correlation table, filtered table, and optional GeoTIFF maps of maximum correlation, optimal lag, and SOS trend. |
| `gpp_growth_msd_lag_analysis.m` | Calculates lagged mean squared deviation between GPP and tree-growth proxies for site-level RWI and gridded aboveground biomass. It also performs the first-vs-last 10-year split-period analysis. | Site-level MSD/lag CSV files, gridded MSD/lag GeoTIFFs, and split-period change table. |
| `cmip6_trendy_gpp_phenology_extraction.R` | Extracts modeled GPP-based start-of-season and end-of-season metrics from CMIP6 and TRENDY monthly GPP fields using `phenofit`. | Annual SOS and EOS GeoTIFF stacks for each model. |
| `structural_equation_modeling_analysis.R` | Performs piecewise structural equation modeling for biome-specific pathways linking SOS trend changes to changes in GPP-growth coupling. | Standardized path coefficients, R² values, d-separation tests, Fisher's C statistics, and diagnostic plots. |
| `xgboost_shap_analysis.py` | Fits XGBoost regression models and calculates SHAP values to assess the relative importance of phenological, climatic, and edaphic predictors. | Model performance, grid-search results, selected hyperparameters, SHAP feature importance, and diagnostic figures. |

## Selected processed data

The adjacent `../source_data/` folder contains selected processed files for reviewer inspection and to illustrate representative input formats. These files are not a complete raw-data archive. Complete data sources are described in the manuscript Data Availability statement.

Please see `../source_data/README_source_data.md` for file-level descriptions and variable information.

## Software requirements

### R

Recommended: R >= 4.0

Required packages:

```r
R.matlab
phenofit
doParallel
foreach
imputeTS
terra
lubridate
dplR
piecewiseSEM
dplyr
readr
```

### Python

Recommended: Python >= 3.9

Required packages:

```text
pandas
numpy
matplotlib
scikit-learn
xgboost
shap
```

### MATLAB

The MATLAB workflow requires standard MATLAB functions plus:

- Statistics and Machine Learning Toolbox
- Mapping Toolbox

The MATLAB scripts use functions including `readgeoraster`, `georasterref`, `geotiffwrite`, `regress`, `detrend`, `corr`, `corrcoef`, and `polyfit`.

## Input data and path settings

The scripts require raw or processed inputs derived from the following data streams:

- International Tree-Ring Data Bank tree-ring width files (`*.rwl`)
- Monthly GPP products used for site-level and gridded assimilation-growth analyses
- Annual SOS raster stacks
- Chloris aboveground biomass products
- CMIP6 and TRENDY model GPP outputs
- Site-level phenological, climatic, soil, and vegetation covariates used for SEM and XGBoost-SHAP analyses

Users should edit the path variables at the beginning of each script before running the analyses. Path placeholders are intentionally retained to avoid hard-coded local file paths.

## Recommended running order

A typical workflow is:

```text
1. tree_ring_data_processing.R
2. rwi_gpp_lag_correlation_and_mapping.m
3. gpp_growth_msd_lag_analysis.m
4. cmip6_trendy_gpp_phenology_extraction.R
5. structural_equation_modeling_analysis.R
6. xgboost_shap_analysis.py
```

The exact order may vary depending on which processed intermediate datasets are already available.

## Notes on computational cost

The gridded AGB-GPP MSD analysis in `gpp_growth_msd_lag_analysis.m` can be computationally intensive because it loops through raster cells and evaluates multiple lag windows. For rapid inspection, the script allows this component to be disabled by setting:

```matlab
run_agb_grid_msd = false;
```

The site-level RWI-GPP analyses are substantially lighter and are generally more suitable for first-pass inspection.

## Workflow notes

- The scripts export processing logs, selected model parameters, model diagnostics, and session information where applicable.
- The Python XGBoost-SHAP workflow uses a fixed `random_state`.
- The R and MATLAB workflows retain explicit user-defined parameters near the top of each script.
- Large raw data files, intermediate raster stacks, and manuscript-scale outputs are not stored in this package and should be obtained from the original data providers or regenerated from the full input datasets.

## License

The code and documentation are distributed under the MIT License. Processed data files may include derivatives of third-party public datasets; reuse of the underlying datasets should follow the licenses, citation requirements, and terms of use specified by the original data providers.

## Citation

If this package is used, please cite the associated manuscript once available.

## Contact

For questions about the analysis workflow, please contact the manuscript authors.
