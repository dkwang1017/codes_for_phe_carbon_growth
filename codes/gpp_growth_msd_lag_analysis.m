% gpp_growth_msd_lag_analysis.m
%
% Purpose
% -------
% This script calculates the temporal association between carbon assimilation
% and tree growth proxies using mean squared deviation (MSD) and lagged GPP
% accumulation windows. It supports two analyses used in the manuscript:
%
%   1. Site-level tree-ring width index (RWI) versus monthly GPP.
%   2. Gridded aboveground biomass (AGB) versus monthly GPP.
%
% The script also includes an optional split-period analysis for RWI sites,
% where MSD and optimal lag are calculated separately for the first and last
% 10-year periods of each site chronology.
%
% Main outputs
% ------------
%   1. CSV file containing site-level minimum MSD, optimal lag, Pearson
%      correlation at optimal lag, MSD components, and SOS trend.
%   2. GeoTIFF files containing gridded AGB-GPP minimum MSD and optimal lag.
%   3. CSV file containing split-period changes in MSD, lag, correlation
%      and MSD components for RWI sites.
%
% Notes for users
% ---------------
% Replace the paths in the "User settings" section with local paths to your
% processed input datasets. The helper functions are included at the end of
% this file, so no external custom MATLAB functions are required.
%
% Required MATLAB toolboxes/functions
% -----------------------------------
% - Statistics and Machine Learning Toolbox: regress
% - Mapping Toolbox: readgeoraster, georasterref, geotiffwrite
% - Standard MATLAB functions: readtable, load, detrend, corrcoef, polyfit
%
% -------------------------------------------------------------------------

clear; clc;

%% ========================= User settings ===============================

% Analysis switches
run_site_level_msd      = true;
run_agb_grid_msd        = false;  % set to true for full gridded AGB-GPP analysis
run_split_period_msd    = true;

% Generic input paths. Replace with local file paths.
rwi_table_path          = "path/to/ITRDB_RWI_site_table_1982_2018.csv";
gpp_mat_path            = "path/to/monthly_GPP_1981_2018_30N.mat";
sos_tif_path            = "path/to/SOS_1982_2022.tif";
agb_mat_path            = "path/to/Chloris_AGB_2003_2019_30N.mat";

% Output directory
output_dir              = "path/to/output/MSD_lag_analysis";
if ~exist(output_dir, "dir")
    mkdir(output_dir);
end

% Variable names in input .mat files
gpp_variable_name       = "GLASS_gpp_gC_sum_monthly_30N";
agb_variable_name       = "Chlorisbiomass_1200720017_03_19_30N";

% Spatial settings for 0.05-degree grids covering 30–90 degrees N
n_rows                  = 3600;
n_cols                  = 7200;
lat_limits              = [30 90];
lon_limits              = [-180 180];

% Temporal settings
gpp_start_year          = 1981;       % first year represented in monthly GPP stack
site_analysis_years     = 1982:2018;  % years used for RWI-GPP analyses
agb_analysis_years      = 2003:2018;  % years used for AGB-GPP analyses

% Lag settings
lag_months              = 0:12;
min_valid_years         = 3;
split_window_years      = 10;         % compare first and last 10-year periods

%% ========================= Load common data ============================

fprintf("Loading monthly GPP data...\n");
gpp_struct = load(gpp_mat_path, gpp_variable_name);
GPP = gpp_struct.(gpp_variable_name);

fprintf("Loading SOS data...\n");
% Use readgeoraster rather than imread so that multi-band or multi-layer
% GeoTIFF stacks are read as annual raster layers rather than as ordinary
% image channels.
[SOS, SOS_ref] = readgeoraster(sos_tif_path); %#ok<NASGU>
SOS = double(SOS);

% The SOS stack is assumed to begin in 1982, with one annual layer per year.
% This check prevents silent errors when only one layer is read or when the
% raster stack does not cover the required analysis period.
required_sos_layers = max(site_analysis_years) - 1982 + 1;
assert(size(SOS, 3) >= required_sos_layers, ...
    "SOS raster stack must contain at least %d annual layers starting in 1982.", ...
    required_sos_layers);

%% ========================= 1. Site-level RWI-GPP MSD ===================

if run_site_level_msd

    fprintf("\n=== Running site-level RWI-GPP MSD analysis ===\n");

    data = readtable(rwi_table_path, "Delimiter", ",");

    % Expected columns in the RWI table:
    %   site, Lat, Lon, year, std
    required_columns = ["site", "Lat", "Lon", "year", "std"];
    assert(all(ismember(required_columns, string(data.Properties.VariableNames))), ...
        "RWI table must contain columns: site, Lat, Lon, year, std.");

    sites = unique(data.site);
    results = {};

    for i = 1:length(sites)

        fprintf("Processing site %d/%d: %s\n", i, length(sites), string(sites(i)));

        site_data = data(strcmp(string(data.site), string(sites(i))), :);

        % Restrict to the target analysis period.
        site_data = site_data(ismember(site_data.year, site_analysis_years), :);
        if height(site_data) < min_valid_years
            continue;
        end

        lat = site_data.Lat(1);
        lon = site_data.Lon(1);
        [r, c] = latlon_to_rowcol(lat, lon, n_rows, n_cols, lat_limits, lon_limits);

        RWI   = site_data.std(:);
        years = site_data.year(:);

        best = compute_best_msd_lag(RWI, years, GPP, r, c, gpp_start_year, lag_months, min_valid_years);

        [sos_trend, sos_p] = compute_sos_trend(SOS, r, c, years);

        results = [results; {string(sites(i)), min(years), max(years), lat, lon, ...
            best.MSD, best.lag, best.r, best.d.bias, best.d.var, best.d.phase, ...
            sos_trend, sos_p}]; %#ok<AGROW>
    end

    varNames = {'site','start_year','end_year','Lat','Lon', ...
        'min_MSD','best_lag','r_at_best', ...
        'd_bias','d_var','d_phase', ...
        'sos_trend','sos_p_value'};

    results_table = cell2table(results, "VariableNames", varNames);
    out_csv = fullfile(output_dir, "site_RWI_GPP_MSD_lag_results_1982_2018.csv");
    writetable(results_table, out_csv);

    fprintf("Site-level MSD results saved to:\n%s\n", out_csv);
end

%% ========================= 2. Gridded AGB-GPP MSD ======================

if run_agb_grid_msd

    fprintf("\n=== Running gridded AGB-GPP MSD analysis ===\n");

    agb_struct = load(agb_mat_path, agb_variable_name);
    AGB_all = agb_struct.(agb_variable_name);

    % The example manuscript analysis used 2003–2018, corresponding to the
    % first 16 annual layers in the processed Chloris AGB stack.
    n_agb_years = length(agb_analysis_years);
    AGB = AGB_all(:, :, 1:n_agb_years);

    [height, width, ~] = size(AGB);
    best_lag_all = nan(height, width);
    min_msd_all  = nan(height, width);

    fprintf("Processing %d rows x %d columns...\n", height, width);

    for r = 1:height

        fprintf("Processing row %d/%d...\n", r, height);

        for c = 1:width

            bio_series = squeeze(AGB(r, c, :));

            if all(isnan(bio_series) | bio_series == 0)
                continue;
            end

            best = compute_best_msd_lag(bio_series, agb_analysis_years(:), GPP, ...
                r, c, gpp_start_year, lag_months, min_valid_years);

            min_msd_all(r, c)  = best.MSD;
            best_lag_all(r, c) = best.lag;
        end
    end

    GeoRef = georasterref("RasterSize", [height width], ...
        "Latlim", lat_limits, "Lonlim", lon_limits);

    lag_tif = fullfile(output_dir, "AGB_GPP_best_lag_2003_2018_MSD.tif");
    msd_tif = fullfile(output_dir, "AGB_GPP_min_MSD_2003_2018.tif");

    geotiffwrite(lag_tif, flipud(best_lag_all), GeoRef);
    geotiffwrite(msd_tif, flipud(min_msd_all), GeoRef);

    fprintf("AGB-GPP optimal lag saved to:\n%s\n", lag_tif);
    fprintf("AGB-GPP minimum MSD saved to:\n%s\n", msd_tif);
end

%% ========================= 3. Split-period RWI-GPP MSD =================

if run_split_period_msd

    fprintf("\n=== Running split-period RWI-GPP MSD analysis ===\n");

    data = readtable(rwi_table_path, "Delimiter", ",");
    required_columns = ["site", "Lat", "Lon", "year", "std"];
    assert(all(ismember(required_columns, string(data.Properties.VariableNames))), ...
        "RWI table must contain columns: site, Lat, Lon, year, std.");

    sites = unique(data.site);
    results = {};

    for i = 1:length(sites)

        fprintf("Processing site %d/%d: %s\n", i, length(sites), string(sites(i)));

        site_data = data(strcmp(string(data.site), string(sites(i))), :);
        site_data = site_data(ismember(site_data.year, site_analysis_years), :);

        if height(site_data) < 2 * split_window_years
            fprintf("  Skipped: fewer than %d records required for first/last %d-year comparison.\n", ...
                2 * split_window_years, split_window_years);
            continue;
        end

        lat = site_data.Lat(1);
        lon = site_data.Lon(1);
        [r, c] = latlon_to_rowcol(lat, lon, n_rows, n_cols, lat_limits, lon_limits);

        years = site_data.year(:);
        RWI   = site_data.std(:);

        % Compare the first and last 10-year periods of each site chronology,
        % matching the manuscript description. If the chronology contains gaps,
        % the first and last 10 available records within site_analysis_years are used.
        [sorted_years, sort_idx] = sort(years);
        sorted_RWI = RWI(sort_idx);

        years1 = sorted_years(1:split_window_years);
        years2 = sorted_years(end-split_window_years+1:end);
        RWI1   = sorted_RWI(1:split_window_years);
        RWI2   = sorted_RWI(end-split_window_years+1:end);

        if length(years1) < split_window_years || length(years2) < split_window_years
            fprintf("  Skipped: one split period has fewer than %d records.\n", split_window_years);
            continue;
        end

        res1 = compute_best_msd_lag(RWI1, years1, GPP, r, c, gpp_start_year, lag_months, min_valid_years);
        res2 = compute_best_msd_lag(RWI2, years2, GPP, r, c, gpp_start_year, lag_months, min_valid_years);

        [sos_trend, sos_p] = compute_sos_trend(SOS, r, c, years);

        dMSD    = res2.MSD     - res1.MSD;
        dlag    = res2.lag     - res1.lag;
        dr      = res2.r       - res1.r;
        ddbias  = res2.d.bias  - res1.d.bias;
        ddvar   = res2.d.var   - res1.d.var;
        ddphase = res2.d.phase - res1.d.phase;

        results = [results; {string(sites(i)), min(years), max(years), lat, lon, ...
            res1.MSD, res1.lag, res1.r, res1.d.bias, res1.d.var, res1.d.phase, ...
            res2.MSD, res2.lag, res2.r, res2.d.bias, res2.d.var, res2.d.phase, ...
            dMSD, dlag, dr, ddbias, ddvar, ddphase, ...
            sos_trend, sos_p}]; %#ok<AGROW>
    end

    varNames = {'site','start_year','end_year','Lat','Lon', ...
        'MSD1','lag1','r1','d_bias1','d_var1','d_phase1', ...
        'MSD2','lag2','r2','d_bias2','d_var2','d_phase2', ...
        'dMSD','dlag','dr','dd_bias','dd_var','dd_phase', ...
        'sos_trend','sos_p_value'};

    results_table = cell2table(results, "VariableNames", varNames);
    out_csv = fullfile(output_dir, "site_RWI_GPP_MSD_lag_first_last_10yr_1982_2018.csv");
    writetable(results_table, out_csv);

    fprintf("Split-period MSD results saved to:\n%s\n", out_csv);
end

fprintf("\nAll requested MSD analyses completed.\n");

%% ========================= Helper functions ============================

function best = compute_best_msd_lag(growth_series, years, GPP, row_idx, col_idx, gpp_start_year, lag_months, min_valid_years)
% compute_best_msd_lag
%
% Calculates the minimum MSD between a tree growth proxy and lagged annual
% GPP accumulation across a set of lagged 12-month windows.
%
% Inputs
% ------
% growth_series  : annual growth proxy, such as RWI or AGB
% years          : years corresponding to growth_series
% GPP            : monthly GPP array with dimensions row x col x month
% row_idx        : row index for the target grid cell
% col_idx        : column index for the target grid cell
% gpp_start_year : first year represented in the monthly GPP stack
% lag_months     : vector of lag windows, usually 0:12
% min_valid_years: minimum number of paired observations
%
% Outputs
% -------
% best.MSD       : minimum MSD across lag windows
% best.lag       : lag month corresponding to minimum MSD
% best.r         : Pearson correlation at the selected lag
% best.d         : MSD decomposition terms: bias, variance and phase

    growth_series = growth_series(:);
    years = years(:);

    best.MSD = NaN;
    best.lag = NaN;
    best.r = NaN;
    best.d = struct("bias", NaN, "var", NaN, "phase", NaN);

    current_min = inf;

    for lag = lag_months

        gpp_annual = nan(length(years), 1);

        for j = 1:length(years)
            yr = years(j);
            idx1 = (yr - gpp_start_year) * 12 + 1 - lag;
            idx2 = idx1 + 11;

            if idx1 >= 1 && idx2 <= size(GPP, 3)
                gpp_annual(j) = sum(GPP(row_idx, col_idx, idx1:idx2), "omitnan");
            end
        end

        valid = ~isnan(gpp_annual) & ~isnan(growth_series) & gpp_annual > 0;
        if sum(valid) < min_valid_years
            continue;
        end

        x = growth_series(valid);
        y = gpp_annual(valid);

        % Detrend and amplitude-normalize both series before MSD
        % decomposition to reduce the influence of units and long-term
        % trends.
        x = detrend(x);
        y = detrend(y);

        x_abs = max(abs(x));
        y_abs = max(abs(y));

        if x_abs > 0
            x = x / x_abs;
        end

        if y_abs > 0
            y = y / y_abs;
        end

        mux = mean(x, "omitnan");
        muy = mean(y, "omitnan");
        sx  = std(x, "omitnan");
        sy  = std(y, "omitnan");

        R = corrcoef(x, y, "Rows", "complete");
        if numel(R) >= 4
            rho = R(1, 2);
        else
            rho = NaN;
        end
        if isnan(rho)
            rho = 0;
        end

        d_bias  = (mux - muy)^2;
        d_var   = (sx - sy)^2;
        d_phase = 2 * sx * sy * (1 - rho);
        MSD     = d_bias + d_var + d_phase;

        if MSD < current_min
            current_min = MSD;
            best.MSD = MSD;
            best.lag = lag;
            best.r = rho;
            best.d = struct("bias", d_bias, "var", d_var, "phase", d_phase);
        end
    end
end

function [sos_trend, sos_p] = compute_sos_trend(SOS, row_idx, col_idx, years)
% compute_sos_trend
%
% Calculates the linear trend in SOS for the years corresponding to the
% target site. The SOS raster stack is assumed to start in 1982, with one
% annual layer per year.

    years = years(:);
    sos_index = years - 1982 + 1;

    valid_layer = sos_index >= 1 & sos_index <= size(SOS, 3);
    years = years(valid_layer);
    sos_index = sos_index(valid_layer);

    if isempty(years)
        sos_trend = NaN;
        sos_p = NaN;
        return;
    end

    sos_series = double(squeeze(SOS(row_idx, col_idx, sos_index)));
    valid = ~isnan(sos_series);

    if sum(valid) > 1
        vy = years(valid);
        vs = sos_series(valid);

        p = polyfit(vy, vs, 1);
        sos_trend = p(1);

        X = [ones(length(vy), 1), vy(:)];
        [~, ~, ~, ~, stats] = regress(vs(:), X);
        sos_p = stats(3);
    else
        sos_trend = NaN;
        sos_p = NaN;
    end
end

function [row_idx, col_idx] = latlon_to_rowcol(lat, lon, n_rows, n_cols, lat_limits, lon_limits)
% latlon_to_rowcol
%
% Converts latitude and longitude coordinates to row and column indices for
% a regular geographic grid. This function assumes that row 1 corresponds to
% the northern edge of the grid and that columns increase from west to east.
%
% Inputs
% ------
% lat, lon   : geographic coordinates
% n_rows     : number of raster rows
% n_cols     : number of raster columns
% lat_limits : [min_lat max_lat]
% lon_limits : [min_lon max_lon]

    min_lat = lat_limits(1);
    max_lat = lat_limits(2);
    min_lon = lon_limits(1);
    max_lon = lon_limits(2);

    lat_res = (max_lat - min_lat) / n_rows;
    lon_res = (max_lon - min_lon) / n_cols;

    row_idx = floor((max_lat - lat) / lat_res) + 1;
    col_idx = floor((lon - min_lon) / lon_res) + 1;

    row_idx = max(1, min(n_rows, row_idx));
    col_idx = max(1, min(n_cols, col_idx));
end
