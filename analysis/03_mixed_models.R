# ============================================================
# 03_mixed_models.R
# Primary mixed-effects analyses for HSI directionality
# ============================================================

if (!exists("table_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")

obj <- readRDS(objects_path)
analysis_paired <- obj$analysis_paired

model_data <- analysis_paired |>
  dplyr::mutate(
    escalation = as.integer(post_decision > pre_decision),
    deescalation = as.integer(post_decision < pre_decision),
    transition_surgery = as.integer(post_decision == 4 & pre_decision != 4),
    transition_discharge = as.integer(post_decision == 1 & pre_decision != 1)
  )

# ----------------------------
# Eligible analysis populations
# ----------------------------

data_escalation <- model_data |>
  dplyr::filter(pre_decision < 4) |>
  droplevels()

data_deescalation <- model_data |>
  dplyr::filter(pre_decision > 1) |>
  droplevels()

data_transition_surgery <- model_data |>
  dplyr::filter(pre_decision %in% c(2, 3)) |>
  droplevels()

# A transition to discharge is possible from any non-discharge starting plan.
# The continuous discharge model is intentionally univariable because all
# observed discharge transitions arose from repeat ultrasound; adding
# initial_plan produces separation and an unstable Hessian.
data_transition_discharge <- model_data |>
  dplyr::filter(pre_decision != 1) |>
  droplevels()

analysis_sets <- list(
  escalation = data_escalation,
  deescalation = data_deescalation,
  transition_surgery = data_transition_surgery,
  transition_discharge = data_transition_discharge
)

# ----------------------------
# Continuous HSI models: univariable + adjusted
# ----------------------------

fit_continuous_models <- function(dat, outcome, outcome_label, adjust_initial = TRUE) {

  univ_formula <- stats::as.formula(
    paste0(
      outcome,
      " ~ hsi_prob_10 + (1 | clinician_id) + (1 | presentation_id)"
    )
  )

  univ_model <- lme4::glmer(
    univ_formula,
    data = dat,
    family = stats::binomial,
    control = glmer_control
  )

  adjusted_model <- NULL
  if (adjust_initial) {
    adjusted_formula <- stats::as.formula(
      paste0(
        outcome,
        " ~ hsi_prob_10 + initial_plan + ",
        "(1 | clinician_id) + (1 | presentation_id)"
      )
    )

    adjusted_model <- lme4::glmer(
      adjusted_formula,
      data = dat,
      family = stats::binomial,
      control = glmer_control
    )
  }

  extract_hsi <- function(model, model_type) {
    coef_table <- summary(model)$coefficients
    beta <- coef_table["hsi_prob_10", "Estimate"]
    se <- coef_table["hsi_prob_10", "Std. Error"]
    p <- coef_table["hsi_prob_10", "Pr(>|z|)"]

    tibble::tibble(
      Outcome = outcome_label,
      Model = model_type,
      Observations = nrow(dat),
      Events = sum(dat[[outcome]], na.rm = TRUE),
      OR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      `OR (95% CI)` = format_or_ci(
        exp(beta),
        exp(beta - 1.96 * se),
        exp(beta + 1.96 * se)
      ),
      `p-value` = format_p(p),
      singular_fit = lme4::isSingular(model, tol = 1e-4)
    )
  }

  results <- extract_hsi(
    univ_model,
    "Univariable mixed-effects"
  )

  if (!is.null(adjusted_model)) {
    results <- dplyr::bind_rows(
      results,
      extract_hsi(adjusted_model, "Adjusted mixed-effects")
    )
  }

  list(
    univariable = univ_model,
    adjusted = adjusted_model,
    results = results
  )
}

continuous_fits <- list(
  escalation = fit_continuous_models(
    data_escalation,
    "escalation",
    "Escalation",
    adjust_initial = TRUE
  ),
  deescalation = fit_continuous_models(
    data_deescalation,
    "deescalation",
    "De-escalation",
    adjust_initial = TRUE
  ),
  transition_surgery = fit_continuous_models(
    data_transition_surgery,
    "transition_surgery",
    "Transition to surgical referral",
    adjust_initial = TRUE
  ),
  transition_discharge = fit_continuous_models(
    data_transition_discharge,
    "transition_discharge",
    "Transition to discharge",
    adjust_initial = FALSE
  )
)

table5_long <- dplyr::bind_rows(
  lapply(continuous_fits, `[[`, "results")
)

table5 <- table5_long |>
  dplyr::select(Outcome, Model, `OR (95% CI)`, `p-value`) |>
  tidyr::pivot_wider(
    names_from = Model,
    values_from = c(`OR (95% CI)`, `p-value`)
  ) |>
  dplyr::rename(
    `Univariable mixed-effects OR (95% CI)` =
      `OR (95% CI)_Univariable mixed-effects`,
    `Univariable p-value` =
      `p-value_Univariable mixed-effects`,
    `Adjusted mixed-effects OR (95% CI)` =
      `OR (95% CI)_Adjusted mixed-effects`,
    `Adjusted p-value` =
      `p-value_Adjusted mixed-effects`
  )

# ----------------------------
# Categorical HSI models
# ----------------------------

fit_categorical_model <- function(dat, outcome, adjust_initial = TRUE) {

  full_formula <- if (adjust_initial) {
    stats::as.formula(
      paste0(
        outcome,
        " ~ risk_group + initial_plan + ",
        "(1 | clinician_id) + (1 | presentation_id)"
      )
    )
  } else {
    stats::as.formula(
      paste0(
        outcome,
        " ~ risk_group + ",
        "(1 | clinician_id) + (1 | presentation_id)"
      )
    )
  }

  reduced_formula <- if (adjust_initial) {
    stats::as.formula(
      paste0(
        outcome,
        " ~ initial_plan + ",
        "(1 | clinician_id) + (1 | presentation_id)"
      )
    )
  } else {
    stats::as.formula(
      paste0(
        outcome,
        " ~ (1 | clinician_id) + (1 | presentation_id)"
      )
    )
  }

  full_model <- lme4::glmer(
    full_formula,
    data = dat,
    family = stats::binomial,
    control = glmer_control
  )

  reduced_model <- lme4::glmer(
    reduced_formula,
    data = dat,
    family = stats::binomial,
    control = glmer_control
  )

  lrt <- stats::anova(reduced_model, full_model, test = "Chisq")
  overall_p <- lrt$`Pr(>Chisq)`[2]

  em <- emmeans::emmeans(
    full_model,
    ~ risk_group,
    type = "response",
    weights = "proportional"
  ) |>
    as.data.frame()

  # emmeans uses asymp.LCL/asymp.UCL for these models.
  em <- em |>
    dplyr::rename(
      probability = prob,
      ci_lower = asymp.LCL,
      ci_upper = asymp.UCL
    )

  list(
    full = full_model,
    reduced = reduced_model,
    overall_p = overall_p,
    probabilities = em
  )
}

categorical_fits <- list(
  escalation = fit_categorical_model(
    data_escalation,
    "escalation",
    adjust_initial = TRUE
  ),
  deescalation = fit_categorical_model(
    data_deescalation,
    "deescalation",
    adjust_initial = TRUE
  )
)

categorical_labels <- c(
  escalation = "Escalation",
  deescalation = "De-escalation"
)

adjusted_probs_long <- dplyr::bind_rows(
  lapply(names(categorical_fits), function(nm) {
    fit <- categorical_fits[[nm]]
    fit$probabilities |>
      dplyr::mutate(
        Outcome = categorical_labels[[nm]],
        `Overall p-value` = fit$overall_p
      )
  })
)

table4 <- adjusted_probs_long |>
  dplyr::mutate(
    formatted_probability = sprintf(
      "%.1f%% (%.1f–%.1f%%)",
      100 * probability,
      100 * ci_lower,
      100 * ci_upper
    )
  ) |>
  dplyr::select(
    Outcome,
    risk_group,
    formatted_probability,
    `Overall p-value`
  ) |>
  tidyr::pivot_wider(
    names_from = risk_group,
    values_from = formatted_probability
  ) |>
  dplyr::mutate(
    `Overall p-value` = format_p(`Overall p-value`)
  ) |>
  dplyr::rename(
    `Green risk, adjusted probability (95% CI)` = Green,
    `Yellow risk, adjusted probability (95% CI)` = Yellow,
    `Red risk, adjusted probability (95% CI)` = Red
  )

# ----------------------------
# Any-revision model by risk category
# ----------------------------

model_any_change <- lme4::glmer(
  changed ~
    risk_group +
    (1 | clinician_id) +
    (1 | presentation_id),
  data = model_data,
  family = stats::binomial,
  control = glmer_control
)

model_any_change_reduced <- lme4::glmer(
  changed ~
    (1 | clinician_id) +
    (1 | presentation_id),
  data = model_data,
  family = stats::binomial,
  control = glmer_control
)

any_change_lrt <- stats::anova(
  model_any_change_reduced,
  model_any_change,
  test = "Chisq"
)

any_change_fixed_effects <- broom.mixed::tidy(
  model_any_change,
  effects = "fixed",
  conf.int = TRUE,
  exponentiate = TRUE
)

# ----------------------------
# Export
# ----------------------------

model_diagnostics <- tibble::tibble(
  Outcome = names(continuous_fits),
  univariable_singular_fit = vapply(
    continuous_fits,
    function(x) lme4::isSingular(x$univariable, tol = 1e-4),
    logical(1)
  ),
  adjusted_singular_fit = vapply(
    continuous_fits,
    function(x) {
      if (is.null(x$adjusted)) return(NA)
      lme4::isSingular(x$adjusted, tol = 1e-4)
    },
    logical(1)
  )
)

openxlsx::write.xlsx(
  list(
    "Table 4" = table4,
    "Table 5" = table5,
    "Table 5 Long" = table5_long,
    "Adjusted Probabilities Long" = adjusted_probs_long,
    "Any Change Fixed Effects" = any_change_fixed_effects,
    "Any Change LRT" = as.data.frame(any_change_lrt),
    "Model Diagnostics" = model_diagnostics
  ),
  file = file.path(table_dir, "mixed_models.xlsx"),
  overwrite = TRUE
)

saveRDS(
  list(
    model_data = model_data,
    analysis_sets = analysis_sets,
    continuous_fits = continuous_fits,
    categorical_fits = categorical_fits,
    adjusted_probs_long = adjusted_probs_long,
    table4 = table4,
    table5 = table5,
    any_change_model = model_any_change,
    any_change_lrt = any_change_lrt
  ),
  file = file.path(derived_dir, "model_results.rds")
)

message("Primary mixed-effects models completed.")
