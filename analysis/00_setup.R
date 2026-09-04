# ============================================================
# 00_setup.R
# Project setup, package checks, constants, and helper functions
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "janitor",
  "lme4",
  "emmeans",
  "broom.mixed",
  "ggplot2",
  "patchwork",
  "openxlsx",
  "DescTools",
  "scales",
  "tibble",
  "here"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the analysis."
  )
}

# Project-relative paths
data_path <- Sys.getenv("HSI_DATA_PATH", unset = here::here("data", "Silent Trial Data.xlsx"))
derived_dir <- here::here("derived")
table_dir <- here::here("output", "tables")
figure_dir <- here::here("output", "figures")

dir.create(derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_path)) {
  stop(
    "Source workbook not found at: ", data_path,
    "\nPlace the workbook in data/Silent Trial Data.xlsx ",
    "or set the HSI_DATA_PATH environment variable."
  )
}

# Study constants
N_CLINICIANS_EXPECTED <- 23L
N_PRESENTATIONS_EXPECTED <- 293L
N_ELIGIBLE_ROWS_EXPECTED <- 6739L
N_PAIRED_ROWS_EXPECTED <- 6733L

HSI_GREEN_CUTOFF <- 0.069
HSI_RED_CUTOFF <- 0.34
HSI_GREEN_DISPLAY <- 0.0689

MANAGEMENT_LABELS <- c(
  "1" = "Discharge",
  "2" = "Repeat ultrasound",
  "3" = "Diuretic renogram",
  "4" = "Surgical referral"
)

RISK_LABELS <- c(
  "1" = "Green",
  "2" = "Yellow",
  "3" = "Red"
)

glmer_control <- lme4::glmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 200000)
)

format_p <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

format_or_ci <- function(or, low, high) {
  sprintf("%.2f (%.2f–%.2f)", or, low, high)
}

format_n_percent <- function(x, denominator = length(x)) {
  n_value <- sum(x, na.rm = TRUE)
  sprintf("%d (%.1f%%)", n_value, 100 * n_value / denominator)
}

format_median_iqr <- function(x, digits = 1) {
  fmt <- paste0("%.", digits, "f")
  sprintf(
    "%s (%s–%s)",
    sprintf(fmt, stats::median(x, na.rm = TRUE)),
    sprintf(fmt, stats::quantile(x, 0.25, na.rm = TRUE)),
    sprintf(fmt, stats::quantile(x, 0.75, na.rm = TRUE))
  )
}

assert_equal <- function(actual, expected, label) {
  if (!identical(as.integer(actual), as.integer(expected))) {
    stop(label, ": expected ", expected, ", observed ", actual)
  }
}
