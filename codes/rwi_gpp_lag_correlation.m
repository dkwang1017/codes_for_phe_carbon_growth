%% rwi_gpp_lag_correlation_and_mapping.m
% -------------------------------------------------------------------------
% Purpose
% -------------------------------------------------------------------------
% This script calculates the lagged relationship between tree-ring width
% index (RWI) and monthly gross primary productivity (GPP) at tree-ring
% sites. For each site, annual GPP is accumulated over 12-month windows
% shifted from 0 to 12 months, and the maximum Pearson correlation between
% RWI and lagged accumulated GPP is retained. The script also estimates the
% site-level trend in the start of the growing season (SOS) and optionally
% maps site-level results to a 0.05-degree grid.
%
% Main outputs
% -------------------------------------------------------------------------
% 1. site_RWI_GPP_lag_correlation_results_1982_2018.csv
%    Site-level maximum RWI-GPP correlation, optimal lag, correlation P value,
%    SOS trend and SOS trend P value.
%
% 2. site_RWI_GPP_lag_correlation_results_1982_2018_filtered.csv
%    Filtered site-level table after removing outliers and invalid values.
%
% 3. Optional GeoTIFF files for site-level gridded visualization:
%    - site_max_correlation_005deg.tif
%    - site_lag_months_005deg.tif
%    - site_sos_trend_005deg.tif
%
% Required input data
% -------------------------------------------------------------------------
% 1. A site-level RWI table with at least the following columns:
%    site, year, Lat, Lon, std
%
% 2. A monthly GPP .mat file containing a 3-D array:
%    latitude x longitude x month
%
% 3. A multi-layer SOS GeoTIFF, where the first layer corresponds to
%    sos_first_year.
%
% Notes for users
% -------------------------------------------------------------------------
% - Replace the generic paths below with local file paths.
% - The function latlon_to_rowcol() assumes a regular global 0.05-degree grid
%   with longitude from -180 to 180 and latitude from 90 to -90.
% - If your raster uses a different grid origin or orientation, adjust
%   latlon_to_rowcol() accordingly.
% - The script is intended to document the analysis workflow used for
%   site-level RWI-GPP lagged correlation analysis.
% -------------------------------------------------------------------------

clear; clc;

%% 1. User settings

% Input files
rwi_table_path = "path/to/ITRDB_RWI_site_table_1982_2018_30N.csv";
gpp_mat_path   = "path/to/monthly_GPP_1981_2018_30N.mat";
sos_tif_path   = "path/to/SOS_1982_2022.tif";

% Output folder
output_dir = "path/to/output/RWI_GPP_lag_correlation";
if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

% Variable name in the GPP .mat file
% Change this if your .mat file uses a different variable name.
gpp_variable_name = "GLASS_gpp_gC_sum_monthly_30N";

% Grid settings
n_rows_global = 3600;
n_cols_global = 7200;

% The output mapping in this example only covers 30-90 degrees N.
% A 0.05-degree grid over 30-90 degrees N has 1200 rows and 7200 columns.
n_rows_nh30 = 1200;
n_cols_nh30 = 7200;

% Temporal settings
gpp_start_year = 1981;
analysis_start_year = 1982;
analysis_end_year   = 2018;
sos_first_year = 1982;

% Lag settings
lag_months = 0:12;

% Filtering settings for the optional gridded visualization step
apply_iqr_filter = true;


%% 2. Load data

fprintf("Loading RWI site table...\n");
data = readtable(rwi_table_path, "Delimiter", ",");

fprintf("Loading monthly GPP data...\n");
gpp_struct = load(gpp_mat_path);
if isfield(gpp_struct, gpp_variable_name)
    GPP = gpp_struct.(gpp_variable_name);
else
    error("The specified GPP variable '%s' was not found in the .mat file.", gpp_variable_name);
end
clear gpp_struct;

fprintf("Loading SOS raster stack...\n");
SOS = imread(sos_tif_path);


%% 3. Site-level lagged RWI-GPP correlation analysis

required_columns = ["site", "year", "Lat", "Lon", "std"];
missing_columns = setdiff(required_columns, string(data.Properties.VariableNames));
if ~isempty(missing_columns)
    error("The RWI table is missing required columns: %s", strjoin(missing_columns, ", "));
end

% Restrict to the target analysis period.
data = data(data.year >= analysis_start_year & data.year <= analysis_end_year, :);

sites = unique(data.site);
results = cell(0, 10);

fprintf("Starting site-level lagged correlation analysis for %d sites...\n", length(sites));

for i = 1:length(sites)

    fprintf("Processing site %d/%d: %s\n", i, length(sites), string(sites{i}));

    site_data = data(strcmp(data.site, sites{i}), :);

    lat = site_data.Lat(1);
    lon = site_data.Lon(1);
    [r, c] = latlon_to_rowcol(lat, lon, n_rows_global, n_cols_global);

    if isnan(r) || isnan(c)
        warning("Site %s has invalid coordinates and was skipped.", string(sites{i}));
        continue;
    end

    RWI_values = site_data.std;
    years = site_data.year;

    start_year_site = min(years);
    end_year_site = max(years);
    num_years = length(years);

    best_corr = -Inf;
    best_lag = NaN;
    best_p_value = NaN;

    for lag = lag_months

        gpp_annual = nan(num_years, 1);

        for j = 1:num_years

            current_year = years(j);

            % 12-month accumulated GPP window shifted backward by lag months.
            start_idx = (current_year - gpp_start_year) * 12 + 1 - lag;
            end_idx = start_idx + 11;

            if start_idx > 0 && end_idx <= size(GPP, 3)
                gpp_slice = squeeze(GPP(r, c, start_idx:end_idx));
                gpp_annual(j) = sum(gpp_slice, "omitnan");
            end
        end

        valid_idx = ~isnan(gpp_annual) & ~isnan(RWI_values);

        if sum(valid_idx) > 2
            gpp_det = detrend(gpp_annual(valid_idx));
            rwi_valid = RWI_values(valid_idx);

            [corr_value, p_value] = corr(rwi_valid, gpp_det, "Type", "Pearson", "Rows", "complete");

            if ~isnan(corr_value) && corr_value > best_corr
                best_corr = corr_value;
                best_lag = lag;
                best_p_value = p_value;
            end
        end
    end

    % Estimate SOS trend for the site.
    sos_series = double(squeeze(SOS(r, c, years - sos_first_year + 1)));
    [sos_trend, sos_p_value] = compute_linear_trend(years, sos_series);

    results = [results; {sites{i}, start_year_site, end_year_site, lat, lon, ...
        best_corr, best_lag, best_p_value, sos_trend, sos_p_value}];

end

results_table = cell2table(results, "VariableNames", ...
    {'site', 'start_year', 'end_year', 'Lat', 'Lon', ...
     'max_correlation', 'lag_months', 'p_value_corr', ...
     'sos_trend', 'p_value_sos'});

output_csv = fullfile(output_dir, "site_RWI_GPP_lag_correlation_results_1982_2018.csv");
writetable(results_table, output_csv);
fprintf("Saved site-level lagged correlation results to:\n%s\n", output_csv);


%% 4. Optional filtering and gridded visualization

results_table_filt = results_table;

if apply_iqr_filter
    results_table_filt.sos_trend(results_table_filt.sos_trend == 0) = NaN;
    results_table_filt.sos_trend = remove_outliers_iqr(results_table_filt.sos_trend);
    results_table_filt.max_correlation = remove_outliers_iqr(results_table_filt.max_correlation);
end

% Remove invalid rows for gridded visualization.
rows_to_remove = isnan(results_table_filt.max_correlation) | ...
                 isinf(results_table_filt.max_correlation) | ...
                 isnan(results_table_filt.lag_months) | ...
                 isnan(results_table_filt.sos_trend);

results_table_filt(rows_to_remove, :) = [];

filtered_csv = fullfile(output_dir, "site_RWI_GPP_lag_correlation_results_1982_2018_filtered.csv");
writetable(results_table_filt, filtered_csv);
fprintf("Saved filtered site-level results to:\n%s\n", filtered_csv);

% Initialize 30-90 degrees N grids.
max_correlation_grid = NaN(n_rows_nh30, n_cols_nh30);
lag_months_grid      = NaN(n_rows_nh30, n_cols_nh30);
sos_trend_grid       = NaN(n_rows_nh30, n_cols_nh30);

for i = 1:height(results_table_filt)

    lat = results_table_filt.Lat(i);
    lon = results_table_filt.Lon(i);

    [row_global, col_global] = latlon_to_rowcol(lat, lon, n_rows_global, n_cols_global);

    % Convert global row index to a 30-90N grid row index.
    % For a global grid ordered from north to south, 30-90N corresponds to
    % the first 1200 rows.
    row_nh30 = row_global;
    col_nh30 = col_global;

    if row_nh30 >= 1 && row_nh30 <= n_rows_nh30 && ...
       col_nh30 >= 1 && col_nh30 <= n_cols_nh30

        max_correlation_grid(row_nh30, col_nh30) = results_table_filt.max_correlation(i);
        lag_months_grid(row_nh30, col_nh30)      = results_table_filt.lag_months(i);
        sos_trend_grid(row_nh30, col_nh30)       = results_table_filt.sos_trend(i);
    end
end

% Write GeoTIFF files.
geo_ref = georasterref("RasterSize", [n_rows_nh30, n_cols_nh30], ...
                       "Latlim", [30, 90], ...
                       "Lonlim", [-180, 180]);

geotiffwrite(fullfile(output_dir, "site_max_correlation_005deg.tif"), ...
             flipud(max_correlation_grid), geo_ref);

geotiffwrite(fullfile(output_dir, "site_lag_months_005deg.tif"), ...
             flipud(lag_months_grid), geo_ref);

geotiffwrite(fullfile(output_dir, "site_sos_trend_005deg.tif"), ...
             flipud(sos_trend_grid), geo_ref);

fprintf("Saved optional GeoTIFF maps to:\n%s\n", output_dir);
fprintf("All processing completed.\n");


%% Helper functions

function [row_idx, col_idx] = latlon_to_rowcol(lat, lon, n_rows, n_cols)
%LATLON_TO_ROWCOL Convert latitude/longitude to row/column indices.
%
% This function assumes a regular geographic grid:
% - longitude range: -180 to 180 degrees
% - latitude range: 90 to -90 degrees
% - rows ordered from north to south
% - columns ordered from west to east
%
% Inputs
% - lat, lon: site coordinates in decimal degrees
% - n_rows, n_cols: raster dimensions
%
% Outputs
% - row_idx, col_idx: integer row and column indices

    if isnan(lat) || isnan(lon) || lat < -90 || lat > 90 || lon < -180 || lon > 180
        row_idx = NaN;
        col_idx = NaN;
        return;
    end

    lat_res = 180 / n_rows;
    lon_res = 360 / n_cols;

    row_idx = floor((90 - lat) / lat_res) + 1;
    col_idx = floor((lon + 180) / lon_res) + 1;

    row_idx = max(1, min(n_rows, row_idx));
    col_idx = max(1, min(n_cols, col_idx));
end


function [slope, p_value] = compute_linear_trend(years, values)
%COMPUTE_LINEAR_TREND Estimate linear trend and P value using ordinary least squares.

    years = years(:);
    values = values(:);

    valid_idx = ~isnan(years) & ~isnan(values);
    valid_years = years(valid_idx);
    valid_values = values(valid_idx);

    if length(valid_years) > 1
        p = polyfit(valid_years, valid_values, 1);
        slope = p(1);

        X = [ones(length(valid_years), 1), valid_years(:)];
        [~, ~, ~, ~, stats] = regress(valid_values(:), X);
        p_value = stats(3);
    else
        slope = NaN;
        p_value = NaN;
    end
end


function x_filtered = remove_outliers_iqr(x)
%REMOVE_OUTLIERS_IQR Replace outliers outside 1.5 x IQR with NaN.

    x_filtered = x;
    valid = ~isnan(x) & ~isinf(x);

    if sum(valid) < 4
        return;
    end

    q1 = prctile(x(valid), 25);
    q3 = prctile(x(valid), 75);
    iqr_value = q3 - q1;

    lower_bound = q1 - 1.5 * iqr_value;
    upper_bound = q3 + 1.5 * iqr_value;

    outliers = x < lower_bound | x > upper_bound;
    x_filtered(outliers) = NaN;
end
