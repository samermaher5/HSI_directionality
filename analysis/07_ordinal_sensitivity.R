# ============================================================
# 07_ordinal_sensitivity.R
# Ordinal mixed-model sensitivity analysis
# ============================================================

if (!requireNamespace("ordinal", quietly = TRUE)) {
  stop("Package 'ordinal' is required for this sensitivity analysis.")
}

if (!exists("table_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")

obj <- readRDS(objects_path)
analysis_paired <- obj$analysis_paired

ordinal_data <- analysis_paired |>
  dplyr::mutate(
    post_plan_ordered = ordered(
      post_plan,
      levels = unname(MANAGEMENT_LABELS)
    ),
    risk_group = factor(
      risk_group,
      levels = c("Green", "Yellow", "Red")
    )
  ) |>
  droplevels()

ordinal_continuous_model <- ordinal::clmm(
  post_plan_ordered ~
    hsi_prob_10 +
    initial_plan +
    (1 | clinician_id) +
    (1 | presentation_id),
  data = ordinal_data,
  link = "logit",
  Hess = TRUE,
  control = ordinal::clmm.control(
    maxIter = 1000,
    maxLineIter = 1000,
    gradTol = 1e-4
  )
)

ordinal_risk_model <- ordinal::clmm(
  post_plan_ordered ~
    risk_group +
    initial_plan +
    (1 | clinician_id) +
    (1 | presentation_id),
  data = ordinal_data,
  link = "logit",
  Hess = TRUE,
  control = ordinal::clmm.control(
    maxIter = 1000,
    maxLineIter = 1000,
    gradTol = 1e-4
  )
)

ordinal_risk_null <- ordinal::clmm(
  post_plan_ordered ~
    initial_plan +
    (1 | clinician_id) +
    (1 | presentation_id),
  data = ordinal_data,
  link = "logit",
  Hess = TRUE,
  control = ordinal::clmm.control(
    maxIter = 1000,
    maxLineIter = 1000,
    gradTol = 1e-4
  )
)

extract_clmm_results <- function(model, model_label) {
  coef_table <- as.data.frame(coef(summary(model))) |>
    tibble::rownames_to_column("term")

  coef_table |>
    dplyr::filter(!grepl("\\|", term)) |>
    dplyr::mutate(
      OR = exp(Estimate),
      CI_low = exp(Estimate - 1.96 * `Std. Error`),
      CI_high = exp(Estimate + 1.96 * `Std. Error`),
      `OR (95% CI)` = format_or_ci(OR, CI_low, CI_high),
      `p-value` = format_p(`Pr(>|z|)`),
      Model = model_label
    )
}

ordinal_continuous_results <- extract_clmm_results(
  ordinal_continuous_model,
  "Continuous displayed HSI probability"
)

ordinal_risk_results <- extract_clmm_results(
  ordinal_risk_model,
  "HSI risk category"
)

ordinal_risk_lrt <- stats::anova(
  ordinal_risk_null,
  ordinal_risk_model
)

openxlsx::write.xlsx(
  list(
    "Continuous HSI" = ordinal_continuous_results,
    "Risk Category" = ordinal_risk_results,
    "Risk Category LRT" = as.data.frame(ordinal_risk_lrt)
  ),
  file = file.path(table_dir, "ordinal_sensitivity.xlsx"),
  overwrite = TRUE
)

message("Ordinal sensitivity analysis completed.")
