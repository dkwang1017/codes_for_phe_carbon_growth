################################################################################
# CMIP6 and TRENDY GPP phenology extraction
#
#
# Main workflow:
#   1. Read monthly GPP arrays from MATLAB .mat files for each CMIP6 or TRENDY model.
#   2. Select the target analysis period.
#   3. Reshape gridded GPP data into pixel-level time series.
#   4. Remove invalid values and retain pixels with sufficient temporal coverage.
#   5. Fill short gaps by spline interpolation.
#   6. Fit seasonal curves using the Elmore double-logistic method in phenofit.
#   7. Extract SOS and EOS using a 50% relative threshold of seasonal amplitude.
#   8. Export annual SOS and EOS layers as GeoTIFF files.
#
#
# Required R packages:
#   R.matlab, phenofit, doParallel, foreach, imputeTS, terra, lubridate
################################################################################

# ---- 1. Load packages ---------------------------------------------------------
required_packages <- c(
  "R.matlab", "phenofit", "doParallel", "foreach", "imputeTS", "terra", "lubridate"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg))
  }
}

library(R.matlab)
library(phenofit)
library(doParallel)
library(foreach)
library(imputeTS)
library(terra)
library(lubridate)

# ---- 2. User-defined parameters ---------------------------------------------
# Input folder containing model-specific .mat files.
input_dir <- "path/to/CMIP6_and_TRENDY/GPP/mat_files"

# Output folder for annual SOS and EOS GeoTIFF files.
output_dir <- "path/to/output/CMIP6_TRENDY_GPP_phenology"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Temporal settings.
# The analysis period follows the manuscript: TRENDY simulations were analysed
# for 1982-2018, and CMIP6 simulations for 1982-2014, because the CMIP6
# historical period ends in 2014.
nptperyear <- 12              # monthly data
threshold <- 0.5              # this example implements the 50% threshold only
if (threshold != 0.5) stop("This example extracts the 50% threshold only.")
years <- 1982:2018            # target analysis period for TRENDY simulations
cmip6_years <- 1982:2014      # target analysis period for CMIP6 simulations

# CMIP6 models in 'model_variable_map' that should use the shorter 1982-2014
# window. Leave empty when every model listed there is a TRENDY simulation.
cmip6_models <- character(0)

# Data quality-control threshold.
# Pixels are retained if at least 70% of monthly observations are valid.
min_valid_fraction <- 0.70
min_valid_points_year <- nptperyear * min_valid_fraction

# Spatial metadata for output rasters.
output_extent <- c(-180, 180, 30, 90)
output_crs <- "EPSG:4326"

# Number of CPU cores for parallel processing.
n_cores <- max(1, parallel::detectCores() - 1)

# Map model names to variable names inside the .mat files.
# Add or remove CMIP6 and TRENDY models as needed. Any CMIP6 model added here
# must also be listed in 'cmip6_models' above so that the 1982-2014 window is
# applied to it.
model_variable_map <- list(
  "CLM5_0" = "CLM5.0.S3.gpp",
  "IBIS"   = "IBIS.S3.gpp",
  "JULES"  = "JULES.S3.gpp"
)

# ---- 3. Phenology extraction function ---------------------------------------
extract_gpp_phenology <- function(gpp_ts, observation_dates, years,
                                  nptperyear = 12,
                                  min_valid_points_total,
                                  min_valid_points_year,
                                  threshold = 0.5) {
  # Replace non-finite values with NA.
  gpp_ts[is.infinite(gpp_ts) | is.nan(gpp_ts)] <- NA

  # Exclude time series with insufficient valid observations.
  if (sum(!is.na(gpp_ts)) < min_valid_points_total) {
    return(rep(NA_real_, length(years) * 2))
  }

  # Fill missing values using spline interpolation.
  gpp_ts <- imputeTS::na_interpolation(gpp_ts, option = "spline")

  ts_df <- data.frame(
    date = observation_dates,
    y = gpp_ts,
    w = rep(1, length(gpp_ts))
  )

  result_vec <- rep(NA_real_, length(years) * 2)

  for (i in seq_along(years)) {
    year_i <- years[i]
    year_data <- ts_df[lubridate::year(ts_df$date) == year_i, ]

    if (sum(!is.na(year_data$y)) < min_valid_points_year) next

    pheno_i <- tryCatch({
      input <- phenofit::check_input(
        t = year_data$date,
        y = year_data$y,
        w = year_data$w,
        nptperyear = nptperyear
      )

      brks <- phenofit::season_mov(
        input,
        list(rFUN = "smooth_wHANTS", wFUN = "wTSM", iters = 3)
      )

      fit <- phenofit::curvefits(
        input,
        brks,
        options = list(methods = "Elmore", wFUN = "wTSM", iters = 3)
      )

      pheno <- phenofit::get_pheno(fit, TRS = threshold)$doy$Elmore
      c(
        mean(pheno$TRS5.sos, na.rm = TRUE),
        mean(pheno$TRS5.eos, na.rm = TRUE)
      )
    }, error = function(e) {
      c(NA_real_, NA_real_)
    })

    result_vec[(i - 1) * 2 + 1] <- pheno_i[1]
    result_vec[(i - 1) * 2 + 2] <- pheno_i[2]
  }

  return(result_vec)
}

# ---- 4. Main model loop ------------------------------------------------------
for (model_name in names(model_variable_map)) {
  message("Processing model: ", model_name)

  # Apply the ensemble-specific analysis window.
  model_years <- if (model_name %in% cmip6_models) cmip6_years else years
  n_years_model <- length(model_years)
  min_valid_points_total <- nptperyear * min_valid_fraction * n_years_model

  model_var <- model_variable_map[[model_name]]
  input_file <- file.path(input_dir, paste0(model_name, "_S3_gpp_1980_2020.mat"))

  if (!file.exists(input_file)) {
    warning("Input file not found: ", input_file, ". Skipping this model.")
    next
  }

  # Read MATLAB file and extract the model-specific GPP array.
  mat_data <- R.matlab::readMat(input_file)
  if (!model_var %in% names(mat_data)) {
    warning("Variable '", model_var, "' not found in ", basename(input_file), ". Skipping this model.")
    next
  }

  gpp_array <- mat_data[[model_var]]

  # Select target years using monthly mid-month timestamps.
  total_layers <- dim(gpp_array)[3]
  observation_dates_all <- seq(as.Date("1980-01-15"), by = "1 month", length.out = total_layers)
  selected_index <- which(lubridate::year(observation_dates_all) %in% model_years)

  observation_dates <- observation_dates_all[selected_index]
  gpp_array <- gpp_array[, , selected_index, drop = FALSE]
  total_layers <- length(selected_index)

  # Reshape raster cube into a pixel x time matrix.
  n_lat <- dim(gpp_array)[1]
  n_lon <- dim(gpp_array)[2]
  n_pixels <- n_lat * n_lon

  gpp_matrix <- matrix(NA_real_, nrow = n_pixels, ncol = total_layers)
  for (t in seq_len(total_layers)) {
    gpp_matrix[, t] <- as.vector(gpp_array[, , t])
  }

  # Replace invalid values and apply unit conversion if required by input data.
  # The multiplier below follows the original processing script. Modify it if
  # the input GPP data are already in the desired units.
  gpp_matrix[gpp_matrix == -9999] <- NA_real_
  gpp_matrix <- gpp_matrix * 1e9

  # Retain pixels with sufficient temporal coverage.
  valid_rows <- apply(gpp_matrix, 1, function(x) sum(!is.na(x)) >= min_valid_points_total)
  gpp_matrix_valid <- gpp_matrix[valid_rows, , drop = FALSE]

  if (nrow(gpp_matrix_valid) == 0) {
    warning("No valid pixels found for model: ", model_name)
    next
  }

  gpp_list <- split(gpp_matrix_valid, row(gpp_matrix_valid))

  # Run phenology extraction in parallel.
  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  pheno_results <- foreach::foreach(
    i = seq_along(gpp_list),
    .combine = rbind,
    .packages = c("phenofit", "imputeTS", "lubridate")
  ) %dopar% {
    extract_gpp_phenology(
      gpp_ts = gpp_list[[i]],
      observation_dates = observation_dates,
      years = model_years,
      nptperyear = nptperyear,
      min_valid_points_total = min_valid_points_total,
      min_valid_points_year = min_valid_points_year,
      threshold = threshold
    )
  }

  parallel::stopCluster(cl)

  # Reconstruct full pixel-level result matrix.
  full_results <- matrix(NA_real_, nrow = n_pixels, ncol = n_years_model * 2)
  full_results[valid_rows, ] <- pheno_results

  # Convert annual SOS and EOS vectors back to raster layers.
  sos_layers <- vector("list", n_years_model)
  eos_layers <- vector("list", n_years_model)

  for (i in seq_along(model_years)) {
    sos_matrix <- matrix(full_results[, (i - 1) * 2 + 1], nrow = n_lat, ncol = n_lon)
    eos_matrix <- matrix(full_results[, (i - 1) * 2 + 2], nrow = n_lat, ncol = n_lon)

    sos_layers[[i]] <- terra::rast(sos_matrix)
    eos_layers[[i]] <- terra::rast(eos_matrix)
  }

  sos_raster <- terra::rast(sos_layers)
  eos_raster <- terra::rast(eos_layers)

  terra::ext(sos_raster) <- output_extent
  terra::ext(eos_raster) <- output_extent
  terra::crs(sos_raster) <- output_crs
  terra::crs(eos_raster) <- output_crs

  names(sos_raster) <- paste0("SOS_", model_years)
  names(eos_raster) <- paste0("EOS_", model_years)

  # Write output files.
  sos_file <- file.path(output_dir, paste0("CMIP6_TRENDY_GPP_", model_name, "_SOS_", min(model_years), "_", max(model_years), "_TRS", threshold * 100, ".tif"))
  eos_file <- file.path(output_dir, paste0("CMIP6_TRENDY_GPP_", model_name, "_EOS_", min(model_years), "_", max(model_years), "_TRS", threshold * 100, ".tif"))

  terra::writeRaster(sos_raster, sos_file, overwrite = TRUE)
  terra::writeRaster(eos_raster, eos_file, overwrite = TRUE)

  message("Finished model: ", model_name)
  message("SOS output: ", sos_file)
  message("EOS output: ", eos_file)
}

# ---- 5. Save session information -------------------------------------------
session_file <- file.path(output_dir, "session_info_CMIP6_TRENDY_GPP_phenology.txt")
writeLines(capture.output(sessionInfo()), con = session_file)
message("All TRENDY GPP phenology extraction completed.")
