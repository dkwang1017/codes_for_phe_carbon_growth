# Source data for the manuscript

**Manuscript title:** Advancing spring phenology weakens the coupling between carbon assimilation and tree growth across the Northern Hemisphere

This folder contains a small set of processed data files prepared for editorial and reviewer inspection. These files are provided to illustrate representative data formats and selected intermediate products used by the analysis scripts. They are not intended to serve as a complete raw-data archive. The complete data sources are described in the manuscript Data Availability statement.

## Folder contents

| File | Description | Spatial/temporal coverage | Main variables or contents | Notes |
|---|---|---|---|---|
| `ITRDB_v8.0_1982_2018_10years_30N.csv` | Site-year tree-ring chronology data derived from ITRDB v8.0 for northern mid- and high-latitude sites. | Northern Hemisphere sites >30°N; 1982–2018. | `site`, `year`, `std`, `samp_depth`, `Lat`, `Lon`, `Species`, `Elevation`, `EPS`. | `std` is the standardized ring-width index chronology; `samp_depth` is chronology sample depth; `EPS` is expressed population signal. |
| `sos_trend_cmip6_trendy_datasets.zip` | Model-derived start-of-season trend rasters from CMIP6 and TRENDY GPP phenology extraction. | CMIP6 files: 1980–2014; TRENDY files: 1980–2020. Northern mid- and high-latitude domain where applicable. | GeoTIFF files named as `sos_trend_slope_<model>_<period>[_CMIP6].tif`. | Values represent the linear trend/slope of GPP-derived SOS timing for each model. Positive values indicate delayed SOS; negative values indicate advanced SOS. |

## Variable dictionary

### `ITRDB_v8.0_1982_2018_10years_30N.csv`

| Variable | Meaning | Unit / format |
|---|---|---|
| `site` | ITRDB site identifier | Character |
| `year` | Calendar year | Year |
| `std` | Standardized ring-width index chronology | Dimensionless |
| `samp_depth` | Number of tree-ring series contributing to the chronology in a given year | Count |
| `Lat` | Site latitude | Decimal degrees |
| `Lon` | Site longitude | Decimal degrees |
| `Species` | Tree species code | Character |
| `Elevation` | Site elevation | m a.s.l. |
| `EPS` | Expressed population signal for chronology quality assessment | Dimensionless |

### `sos_trend_cmip6_trendy_datasets.zip`

The ZIP archive contains GeoTIFF rasters for model-specific SOS trends. The file naming convention is:

```text
sos_trend_slope_<model_name>_<start_year>_<end_year>[_CMIP6].tif
```

where `<model_name>` identifies the CMIP6 or TRENDY model, and `<start_year>_<end_year>` indicates the temporal period used for trend estimation.

## Relationship to analysis code

These processed files are used together with the scripts provided in the accompanying `../codes/` folder. In particular:

- `ITRDB_v8.0_1982_2018_10years_30N.csv` provides site-level tree-ring chronology data used for tree-growth analyses and GPP-growth coupling calculations.
- `sos_trend_cmip6_trendy_datasets.zip` provides model-level SOS trend rasters used to compare observed and model-derived phenological trends.

The full workflow also requires additional public datasets, including satellite GPP/phenology products, reanalysis climate variables, biomass products, and Earth system model outputs. These should be downloaded from their original repositories as described in the manuscript Data Availability statement.

## Recommended use

1. Unzip `sos_trend_cmip6_trendy_datasets.zip` before running scripts that require model-specific SOS trend rasters.
2. Keep file names unchanged, because the analysis scripts use file names to identify model names and analysis periods.
3. Update local paths in the script user-settings sections before running the analyses.
4. Do not treat these files as complete raw-data archives; they are selected processed data and intermediate products provided for reviewer inspection.

## Data reuse and licensing

The code and documentation are distributed under the license provided in the package root. The data files in this folder may include processed derivatives of third-party public datasets. Reuse of the underlying datasets should follow the licenses, citation requirements, and terms of use specified by the original data providers.
