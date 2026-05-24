################################################################################
# Tree-ring data processing
# Manuscript: Advancing spring phenology weakens the coupling between carbon
# assimilation and tree growth across the Northern Hemisphere
#
# Purpose
# -------
# This script processes raw ITRDB tree-ring width files (*.rwl), detrends each
# tree-ring series, constructs site-level ring-width index (RWI) chronologies,
# calculates expressed population signal (EPS), and exports the resulting
# chronology table used for subsequent carbon assimilation-growth analyses.
#
# Main workflow
# 1. Recursively search all *.rwl files under the ITRDB measurements directory.
# 2. Exclude files containing "noaa" in the file name, following the original
#    data-cleaning workflow.
# 3. Read raw tree-ring width data in Tucson format using dplR::read.rwl().
# 4. Detrend raw tree-ring series using a cubic smoothing spline with a 30-year
#    50% frequency-response cutoff.
# 5. Build site-level standard chronologies using a biweight robust mean.
# 6. Calculate EPS for each site.
# 7. Export individual *.crn files, EPS statistics, the full RWI chronology, and
#    the filtered RWI chronology for 1982-2022 with sample depth > 5.
#
# Notes
# -----
# This script is a streamlined global processing workflow. It does not process
# continents or regions separately. Instead, it treats all ITRDB measurement files
# under the user-defined input directory as one combined collection and produces
# one global chronology table.
################################################################################

# ==============================================================================
# 0. Required packages
# ==============================================================================

required_packages <- c("dplR")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package '%s' is required but not installed.", pkg), call. = FALSE)
  }
}

library(dplR)

# ==============================================================================
# 1. User settings
# ==============================================================================

# Root directory containing raw ITRDB *.rwl measurement files.
# Edit this path to the local directory where the ITRDB measurement files are stored.
# The script searches this directory recursively; therefore, files can be organized
# in regional subfolders or stored directly under this folder.
input_dir <- "path/to/ITRDB/treering/measurements"

# Output directory for derived chronologies and summary tables.
# Edit this path if a different output location is preferred.
output_dir <- "path/to/output/tree_ring_chronologies"

# Analysis period and chronology quality threshold.
start_year <- 1982
end_year <- 2022
min_sample_depth <- 5

# Detrending settings. These match the manuscript Methods:
# cubic smoothing spline with a 50% frequency-response cutoff at 30 years.
spline_nyrs <- 30
spline_f <- 0.5

# Whether to prewhiten chronology series. The manuscript workflow uses raw
# standardized chronology without prewhitening.
use_prewhiten <- FALSE

# ==============================================================================
# 2. Helper functions
# ==============================================================================

create_dir_if_missing <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

safe_read_rwl <- function(file_path) {
  tryCatch(
    {
      read.rwl(file_path, header = TRUE, format = "tucson")
    },
    error = function(e) {
      message(sprintf("[SKIP] Failed to read RWL file: %s", basename(file_path)))
      message(sprintf("       Reason: %s", e$message))
      return(NULL)
    }
  )
}

safe_detrend <- function(rwl_data, site_id) {
  tryCatch(
    {
      detrend(rwl_data, method = "Spline", nyrs = spline_nyrs, f = spline_f)
    },
    error = function(e) {
      message(sprintf("[SKIP] Detrending failed for site: %s", site_id))
      message(sprintf("       Reason: %s", e$message))
      return(NULL)
    }
  )
}

safe_chron <- function(rwi_data, site_id) {
  tryCatch(
    {
      chron(rwi_data, prefix = "RWI", biweight = TRUE, prewhiten = use_prewhiten)
    },
    error = function(e) {
      message(sprintf("[SKIP] Chronology construction failed for site: %s", site_id))
      message(sprintf("       Reason: %s", e$message))
      return(NULL)
    }
  )
}

safe_rwi_stats <- function(rwi_data, site_id) {
  tryCatch(
    {
      stats <- rwi.stats(rwi_data, ids = NULL, period = "max")
      eps_value <- stats$eps
      if (length(eps_value) == 0) eps_value <- NA_real_
      data.frame(site = site_id, eps = as.numeric(eps_value), stringsAsFactors = FALSE)
    },
    error = function(e) {
      message(sprintf("[WARN] EPS calculation failed for site: %s", site_id))
      message(sprintf("       Reason: %s", e$message))
      data.frame(site = site_id, eps = NA_real_, stringsAsFactors = FALSE)
    }
  )
}

standardize_chronology_table <- function(crn_data, site_id) {
  # dplR::chron() returns a data frame with years as row names and typically
  # two columns: standardized chronology and sample depth. This helper function
  # makes the output explicit and robust for downstream analysis.
  out <- as.data.frame(crn_data)
  out$year <- suppressWarnings(as.integer(rownames(out)))
  out$site <- site_id

  # Reorder columns for readability: site, year, chronology columns.
  out <- out[, c("site", "year", setdiff(names(out), c("site", "year"))), drop = FALSE]
  rownames(out) <- NULL
  out
}

# ==============================================================================
# 3. Main processing
# ==============================================================================

start_time <- Sys.time()
create_dir_if_missing(output_dir)

message("Searching for raw ITRDB RWL files...")
rwl_files <- dir(input_dir, pattern = "\\.rwl$", full.names = TRUE, recursive = TRUE)
rwl_files <- rwl_files[!grepl("noaa", basename(rwl_files), ignore.case = TRUE)]

if (length(rwl_files) == 0) {
  stop(sprintf("No .rwl files were found under: %s", input_dir), call. = FALSE)
}

message(sprintf("Number of RWL files found: %d", length(rwl_files)))
message(sprintf("Output directory: %s", output_dir))

all_eps <- list()
all_crn <- list()
processing_log <- data.frame(
  site = character(),
  file = character(),
  status = character(),
  note = character(),
  stringsAsFactors = FALSE
)

for (i in seq_along(rwl_files)) {
  current_file <- rwl_files[i]
  site_id <- sub("\\.rwl$", "", basename(current_file), ignore.case = TRUE)

  message(sprintf("[%d/%d] Processing site: %s", i, length(rwl_files), site_id))

  rwl_data <- safe_read_rwl(current_file)
  if (is.null(rwl_data) || !any(!is.na(rwl_data))) {
    processing_log <- rbind(
      processing_log,
      data.frame(site = site_id, file = current_file, status = "skipped", note = "empty_or_read_failed")
    )
    next
  }

  rwi_data <- safe_detrend(rwl_data, site_id)
  if (is.null(rwi_data)) {
    processing_log <- rbind(
      processing_log,
      data.frame(site = site_id, file = current_file, status = "skipped", note = "detrending_failed")
    )
    next
  }

  crn_data <- safe_chron(rwi_data, site_id)
  if (is.null(crn_data)) {
    processing_log <- rbind(
      processing_log,
      data.frame(site = site_id, file = current_file, status = "skipped", note = "chronology_failed")
    )
    next
  }

  eps_data <- safe_rwi_stats(rwi_data, site_id)
  all_eps[[site_id]] <- eps_data

  crn_table <- standardize_chronology_table(crn_data, site_id)
  all_crn[[site_id]] <- crn_table

  crn_file <- file.path(output_dir, paste0(site_id, ".crn"))
  tryCatch(
    {
      write.crn(crn_data, fname = crn_file)
    },
    error = function(e) {
      message(sprintf("[WARN] Failed to write CRN file for site: %s", site_id))
      message(sprintf("       Reason: %s", e$message))
    }
  )

  processing_log <- rbind(
    processing_log,
    data.frame(site = site_id, file = current_file, status = "processed", note = "ok")
  )
}

# ==============================================================================
# 4. Export summary tables
# ==============================================================================

if (length(all_crn) == 0) {
  stop("No valid chronologies were produced. Please check input files and paths.", call. = FALSE)
}

combined_crn <- do.call(rbind, all_crn)
combined_eps <- do.call(rbind, all_eps)

# Identify sample-depth column. dplR chronology tables usually contain
# 'samp.depth'. This check makes the script safer if column names differ.
sample_depth_col <- intersect(names(combined_crn), c("samp.depth", "sample.depth", "n.series"))
if (length(sample_depth_col) == 0) {
  warning("No sample-depth column was detected. The filtered RWI table will only use the year filter.")
  filtered_rwi <- subset(combined_crn, year >= start_year & year <= end_year)
} else {
  sample_depth_col <- sample_depth_col[1]
  filtered_rwi <- combined_crn[
    combined_crn$year >= start_year &
      combined_crn$year <= end_year &
      combined_crn[[sample_depth_col]] > min_sample_depth,
    ,
    drop = FALSE
  ]
}

# Output files.
eps_file <- file.path(output_dir, "EPS_global.csv")
crn_full_file <- file.path(output_dir, "Global_RWI_full.csv")
rwi_filtered_file <- file.path(output_dir, sprintf("Global_RWI_%d_%d_sampdepth_gt_%d.csv", start_year, end_year, min_sample_depth))
log_file <- file.path(output_dir, "processing_log_global.csv")
session_file <- file.path(output_dir, "session_info_tree_ring_processing.txt")

write.csv(combined_eps, eps_file, row.names = FALSE)
write.csv(combined_crn, crn_full_file, row.names = FALSE)
write.csv(filtered_rwi, rwi_filtered_file, row.names = FALSE)
write.csv(processing_log, log_file, row.names = FALSE)

sink(session_file)
cat("Tree-ring data processing session information\n")
cat("====================================================\n\n")
cat("Input directory:\n", input_dir, "\n\n", sep = "")
cat("Output directory:\n", output_dir, "\n\n", sep = "")
cat("Number of input RWL files:\n", length(rwl_files), "\n", sep = "")
cat("Number of processed sites:\n", length(all_crn), "\n", sep = "")
cat("Analysis period:\n", start_year, "-", end_year, "\n", sep = "")
cat("Minimum sample depth:\n", min_sample_depth, "\n\n", sep = "")
cat("Detrending method: cubic smoothing spline\n")
cat("Spline nyrs:\n", spline_nyrs, "\n", sep = "")
cat("Spline f:\n", spline_f, "\n", sep = "")
cat("Chronology construction: biweight robust mean\n")
cat("Prewhitening:\n", use_prewhiten, "\n\n", sep = "")
print(sessionInfo())
sink()

end_time <- Sys.time()
message("Processing completed.")
message(sprintf("EPS table saved to: %s", eps_file))
message(sprintf("Full RWI table saved to: %s", crn_full_file))
message(sprintf("Filtered RWI table saved to: %s", rwi_filtered_file))
message(sprintf("Processing log saved to: %s", log_file))
message(sprintf("Session info saved to: %s", session_file))
message(sprintf("Total runtime: %.2f minutes", as.numeric(difftime(end_time, start_time, units = "mins"))))
