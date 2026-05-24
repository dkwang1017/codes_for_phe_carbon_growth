################################################################################
# Structural equation modeling for phenology-driven changes in GPP-growth coupling
#
# Purpose:
#   This script performs piecewise structural equation modeling (SEM) to evaluate
#   potential environmental and physiological pathways linking changes in the
#   start-of-season (SOS) trend to changes in the coupling between gross primary
#   productivity and tree growth (dR). Separate SEMs are fitted for evergreen
#   needleleaf forests (ENF) and deciduous broadleaf forests (DBF).
#
# Main workflow:
#   1. Load the prepared site-level data table.
#   2. Filter records by forest type.
#   3. Remove incomplete records and extreme values.
#   4. Standardize variables before SEM fitting.
#   5. Fit piecewise SEMs for ENF and DBF using linear models.
#   6. Export standardized path coefficients, R-squared values, d-separation
#      tests and Fisher's C statistics.
#
# Input:
#   A processed CSV file containing site-level variables, including:
#     - site
#     - IGBP
#     - dR
#     - delta_sos_trend
#     - delta_sm_spring
#     - delta_GDD_spring
#     - delta_lai_spring
#     - delta_gpp_spring
#
# Output:
#   CSV files containing SEM coefficients, R-squared values, d-separation results
#   and Fisher's C statistics for ENF and DBF.
#
# Notes for users:
#   Replace input_file and output_dir with local paths before running this script.
#   This script assumes that the input data have already been harmonized across
#   tree-ring sites, phenology products, GPP datasets and environmental covariates.
################################################################################

# =========================
# 0. Load required packages
# =========================
required_packages <- c("piecewiseSEM", "dplyr", "readr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(piecewiseSEM)
library(dplyr)
library(readr)

# =========================
# 1. User-defined paths
# =========================
input_file <- "path/to/sem_environmental_drivers.csv"
output_dir <- "path/to/output/SEM_analysis"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =========================
# 2. Analysis parameters
# =========================
site_id_pattern <- "^[A-Za-z]+[0-9]+$"
outlier_sd_threshold <- 3

# Variables used in ENF and DBF SEMs
vars_enf <- c(
  "delta_sm_spring",
  "delta_GDD_spring",
  "delta_lai_spring",
  "dR",
  "delta_sos_trend"
)

vars_dbf <- c(
  "delta_sm_spring",
  "delta_lai_spring",
  "delta_gpp_spring",
  "dR",
  "delta_sos_trend"
)

# =========================
# 3. Helper functions
# =========================
preprocess_sem_data <- function(df, selected_vars, outlier_sd = 3) {
  # Select variables, remove incomplete cases, remove extreme values based on
  # a +/- outlier_sd standard deviation criterion, and standardize all variables.
  df_clean <- df %>%
    dplyr::select(dplyr::all_of(selected_vars)) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::filter(
      dplyr::if_all(
        dplyr::everything(),
        ~ abs((.x - mean(.x, na.rm = TRUE)) / sd(.x, na.rm = TRUE)) <= outlier_sd
      )
    ) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ as.numeric(scale(.x))))

  return(df_clean)
}

safe_write_csv <- function(x, file_name) {
  readr::write_csv(as.data.frame(x), file.path(output_dir, file_name))
}

export_sem_outputs <- function(sem_model, model_name) {
  # Export standardized path coefficients, R-squared values, d-separation tests,
  # and Fisher's C goodness-of-fit statistics.
  coefficients <- piecewiseSEM::coefs(sem_model, standardize = "scale")
  r2_values <- piecewiseSEM::rsquared(sem_model)

  dsep_results <- try(piecewiseSEM::dSep(sem_model), silent = TRUE)
  fisher_c <- try(piecewiseSEM::fisherC(sem_model), silent = TRUE)

  safe_write_csv(coefficients, paste0(model_name, "_standardized_path_coefficients.csv"))
  safe_write_csv(r2_values, paste0(model_name, "_r_squared.csv"))

  if (!inherits(dsep_results, "try-error")) {
    safe_write_csv(dsep_results, paste0(model_name, "_d_separation_tests.csv"))
  }

  if (!inherits(fisher_c, "try-error")) {
    safe_write_csv(fisher_c, paste0(model_name, "_fishers_C.csv"))
  }

  # Save printed SEM summary for reviewer inspection.
  sink(file.path(output_dir, paste0(model_name, "_SEM_summary.txt")))
  cat("========== ", model_name, " SEM ==========" , "\n\n", sep = "")
  cat("Sample size:", nrow(sem_model$data), "\n\n")
  cat("--- Standardized path coefficients ---\n")
  print(coefficients)
  cat("\n--- R-squared values ---\n")
  print(r2_values)
  cat("\n--- d-separation tests ---\n")
  print(dsep_results)
  cat("\n--- Fisher's C statistic ---\n")
  print(fisher_c)
  sink()
}

# =========================
# 4. Load and filter data
# =========================
dat <- readr::read_csv(input_file, show_col_types = FALSE) %>%
  dplyr::filter(grepl(site_id_pattern, site))

dat_enf <- dat %>% dplyr::filter(IGBP == "ENF")
dat_dbf <- dat %>% dplyr::filter(IGBP == "DBF")

# =========================
# 5. Preprocess data
# =========================
dat_enf_sem <- preprocess_sem_data(
  dat_enf,
  selected_vars = vars_enf,
  outlier_sd = outlier_sd_threshold
)

dat_dbf_sem <- preprocess_sem_data(
  dat_dbf,
  selected_vars = vars_dbf,
  outlier_sd = outlier_sd_threshold
)

readr::write_csv(dat_enf_sem, file.path(output_dir, "ENF_SEM_input_standardized.csv"))
readr::write_csv(dat_dbf_sem, file.path(output_dir, "DBF_SEM_input_standardized.csv"))

# =========================
# 6. Fit SEM for ENF
# =========================
# Ecological interpretation:
#   The ENF model evaluates whether changes in SOS trend affect dR through changes
#   in spring soil moisture, spring growing-degree days and spring LAI.
sem_enf <- piecewiseSEM::psem(
  lm(delta_sm_spring ~ delta_sos_trend, data = dat_enf_sem),
  lm(delta_GDD_spring ~ delta_sm_spring + delta_sos_trend, data = dat_enf_sem),
  lm(delta_lai_spring ~ delta_sm_spring + delta_GDD_spring + delta_sos_trend, data = dat_enf_sem),
  lm(dR ~ delta_sm_spring + delta_GDD_spring + delta_lai_spring, data = dat_enf_sem),
  data = dat_enf_sem
)

# =========================
# 7. Fit SEM for DBF
# =========================
# Ecological interpretation:
#   The DBF model evaluates whether changes in SOS trend affect dR through changes
#   in spring soil moisture, spring LAI and spring GPP.
sem_dbf <- piecewiseSEM::psem(
  lm(delta_sm_spring ~ delta_sos_trend, data = dat_dbf_sem),
  lm(delta_lai_spring ~ delta_sm_spring + delta_sos_trend, data = dat_dbf_sem),
  lm(delta_gpp_spring ~ delta_sm_spring + delta_lai_spring + delta_sos_trend, data = dat_dbf_sem),
  lm(dR ~ delta_sm_spring + delta_lai_spring + delta_gpp_spring, data = dat_dbf_sem),
  data = dat_dbf_sem
)

# =========================
# 8. Export outputs
# =========================
export_sem_outputs(sem_enf, "ENF")
export_sem_outputs(sem_dbf, "DBF")

# Optional quick visualization. These plots are intended for diagnostic inspection
# only; manuscript figures may require custom plotting based on exported tables.
pdf(file.path(output_dir, "SEM_diagnostic_plots.pdf"), width = 8, height = 6)
plot(sem_enf)
plot(sem_dbf)
dev.off()

# Save R session information for reproducibility.
sink(file.path(output_dir, "session_info_SEM_analysis.txt"))
print(sessionInfo())
sink()

message("SEM analysis completed. Outputs were saved to: ", output_dir)
