# ============================================================
# 06_exploratory_review_order.R
# Exploratory review-order analysis
# ============================================================

# IMPORTANT:
# Exact session boundaries are not encoded in the source workbook used here.
# This script therefore divides the common case sequence into seven
# approximately equal ordered blocks as a proxy for the seven review sessions.
# If the exact session mapping becomes available, replace this proxy before
# describing the variable as an exact session number.

if (!exists("table_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")

obj <- readRDS(objects_path)
analysis_paired <- obj$analysis_paired
clinician_lookup <- obj$clinician_lookup

review_block_roster <- analysis_paired |>
  dplyr::distinct(presentation_order) |>
  dplyr::arrange(presentation_order) |>
  dplyr::mutate(
    review_block = dplyr::ntile(presentation_order, 7)
  )

concordance_data <- analysis_paired |>
  dplyr::left_join(
    review_block_roster,
    by = "presentation_order"
  ) |>
  dplyr::left_join(
    clinician_lookup |>
      dplyr::select(
        clinician_id,
        specialty,
        experience_group
      ),
    by = "clinician_id"
  ) |>
  dplyr::mutate(
    review_block = as.integer(review_block),

    pre_hsi_concordant = dplyr::case_when(
      risk_group == "Green" & pre_decision %in% c(1, 2) ~ 1L,
      risk_group == "Yellow" & pre_decision %in% c(2, 3) ~ 1L,
      risk_group == "Red" & pre_decision %in% c(3, 4) ~ 1L,
      TRUE ~ 0L
    ),

    post_hsi_concordant = dplyr::case_when(
      risk_group == "Green" & post_decision %in% c(1, 2) ~ 1L,
      risk_group == "Yellow" & post_decision %in% c(2, 3) ~ 1L,
      risk_group == "Red" & post_decision %in% c(3, 4) ~ 1L,
      TRUE ~ 0L
    ),

    became_concordant = as.integer(
      pre_hsi_concordant == 0 &
        post_hsi_concordant == 1
    )
  )

session_summary <- concordance_data |>
  dplyr::group_by(review_block) |>
  dplyr::summarise(
    observations = dplyr::n(),
    pre_hsi_discrepant_percent =
      100 * mean(pre_hsi_concordant == 0),
    post_hsi_discrepant_percent =
      100 * mean(post_hsi_concordant == 0),
    became_concordant_percent_among_initially_discordant =
      100 * mean(
        became_concordant[pre_hsi_concordant == 0]
      ),
    .groups = "drop"
  )

fit_review_order_model <- function(data, outcome, label) {
  dat <- data |>
    dplyr::filter(!is.na(.data[[outcome]]))

  model <- lme4::glmer(
    stats::as.formula(
      paste0(
        outcome,
        " ~ review_block + ",
        "(1 | clinician_id) + (1 | presentation_id)"
      )
    ),
    data = dat,
    family = stats::binomial,
    control = glmer_control
  )

  coefs <- summary(model)$coefficients
  beta <- coefs["review_block", "Estimate"]
  se <- coefs["review_block", "Std. Error"]
  p <- coefs["review_block", "Pr(>|z|)"]

  tibble::tibble(
    Outcome = label,
    `OR per review block` = exp(beta),
    `95% CI` = format_or_ci(
      exp(beta),
      exp(beta - 1.96 * se),
      exp(beta + 1.96 * se)
    ),
    `p-value` = format_p(p)
  )
}

review_order_models <- dplyr::bind_rows(
  fit_review_order_model(
    concordance_data,
    "pre_hsi_concordant",
    "Initial recommendation HSI-concordant"
  ),
  fit_review_order_model(
    concordance_data,
    "post_hsi_concordant",
    "Post-HSI recommendation HSI-concordant"
  ),
  fit_review_order_model(
    dplyr::filter(concordance_data, pre_hsi_concordant == 0),
    "became_concordant",
    "Became HSI-concordant among initially discordant"
  )
)

openxlsx::write.xlsx(
  list(
    "Review Block Summary" = session_summary,
    "Trend Models" = review_order_models
  ),
  file = file.path(table_dir, "exploratory_review_order.xlsx"),
  overwrite = TRUE
)

message(
  "Exploratory review-order analysis completed using seven equal ordered blocks."
)
