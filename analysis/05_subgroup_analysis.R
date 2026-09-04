# ============================================================
# 05_subgroup_analysis.R
# Exploratory subgroup analyses by specialty and experience
# ============================================================

if (!exists("table_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
models_path <- file.path(derived_dir, "model_results.rds")

if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")
if (!file.exists(models_path)) source("analysis/03_mixed_models.R")

obj <- readRDS(objects_path)
mod <- readRDS(models_path)

analysis_paired <- obj$analysis_paired
clinician_lookup <- obj$clinician_lookup

df_subgroup <- analysis_paired |>
  dplyr::mutate(
    escalation = as.integer(post_decision > pre_decision),
    deescalation = as.integer(post_decision < pre_decision),
    transition_surgery = as.integer(post_decision == 4)
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
  dplyr::filter(
    !is.na(specialty),
    !is.na(experience_group)
  )

analysis_sets <- list(
  escalation = df_subgroup |>
    dplyr::filter(pre_decision < 4) |>
    droplevels(),
  deescalation = df_subgroup |>
    dplyr::filter(pre_decision > 1) |>
    droplevels(),
  transition_surgery = df_subgroup |>
    dplyr::filter(pre_decision %in% c(2, 3)) |>
    droplevels()
)

run_subgroup_model <- function(
    data,
    outcome,
    outcome_label,
    subgroup_var,
    subgroup_label,
    adjust_for_specialty = FALSE
) {

  dat <- data |>
    dplyr::filter(
      !is.na(.data[[outcome]]),
      !is.na(.data[[subgroup_var]])
    ) |>
    dplyr::mutate(
      subgroup = factor(.data[[subgroup_var]])
    ) |>
    droplevels()

  specialty_term <- if (adjust_for_specialty) " + specialty" else ""

  full_formula <- stats::as.formula(
    paste0(
      outcome,
      " ~ hsi_prob_10 * subgroup + initial_plan",
      specialty_term,
      " + (1 | clinician_id) + (1 | presentation_id)"
    )
  )

  reduced_formula <- stats::as.formula(
    paste0(
      outcome,
      " ~ hsi_prob_10 + subgroup + initial_plan",
      specialty_term,
      " + (1 | clinician_id) + (1 | presentation_id)"
    )
  )

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

  interaction_test <- stats::anova(
    reduced_model,
    full_model,
    test = "Chisq"
  )

  p_interaction <- interaction_test$`Pr(>Chisq)`[2]

  trends_raw <- emmeans::emtrends(
    full_model,
    specs = ~ subgroup,
    var = "hsi_prob_10"
  ) |>
    summary(infer = TRUE) |>
    as.data.frame()

  trend_col <- grep("\\.trend$", names(trends_raw), value = TRUE)[1]

  subgroup_trends <- trends_raw |>
    dplyr::mutate(
      OR_per_10pp = exp(.data[[trend_col]]),
      CI_low = exp(asymp.LCL),
      CI_high = exp(asymp.UCL)
    ) |>
    dplyr::select(
      subgroup,
      OR_per_10pp,
      CI_low,
      CI_high,
      p.value
    )

  subgroup_counts <- dat |>
    dplyr::group_by(subgroup) |>
    dplyr::summarise(
      clinicians_n = dplyr::n_distinct(clinician_id),
      observations_n = dplyr::n(),
      events_n = sum(.data[[outcome]], na.rm = TRUE),
      event_rate = events_n / observations_n,
      .groups = "drop"
    )

  subgroup_trends |>
    dplyr::left_join(subgroup_counts, by = "subgroup") |>
    dplyr::transmute(
      Outcome = outcome_label,
      `Subgroup type` = subgroup_label,
      Subgroup = subgroup,
      `Clinicians, n` = clinicians_n,
      `Observations, n` = observations_n,
      `Events, n/N (%)` = sprintf(
        "%d/%d (%.1f%%)",
        events_n,
        observations_n,
        100 * event_rate
      ),
      `OR per 10-percentage-point increase in displayed HSI risk (95% CI)` =
        format_or_ci(OR_per_10pp, CI_low, CI_high),
      `Within-subgroup p-value` = format_p(p.value),
      `Interaction p-value` = format_p(p_interaction),
      `Experience model adjusted for specialty` = adjust_for_specialty,
      `Singular fit` = lme4::isSingular(full_model, tol = 1e-4)
    )
}

subgroup_results <- dplyr::bind_rows(
  run_subgroup_model(
    analysis_sets$escalation,
    "escalation",
    "Escalation",
    "specialty",
    "Clinician specialty"
  ),
  run_subgroup_model(
    analysis_sets$deescalation,
    "deescalation",
    "De-escalation",
    "specialty",
    "Clinician specialty"
  ),
  run_subgroup_model(
    analysis_sets$transition_surgery,
    "transition_surgery",
    "Transition to surgical referral",
    "specialty",
    "Clinician specialty"
  ),
  run_subgroup_model(
    analysis_sets$escalation,
    "escalation",
    "Escalation",
    "experience_group",
    "Clinician experience level",
    adjust_for_specialty = TRUE
  ),
  run_subgroup_model(
    analysis_sets$deescalation,
    "deescalation",
    "De-escalation",
    "experience_group",
    "Clinician experience level",
    adjust_for_specialty = TRUE
  ),
  run_subgroup_model(
    analysis_sets$transition_surgery,
    "transition_surgery",
    "Transition to surgical referral",
    "experience_group",
    "Clinician experience level",
    adjust_for_specialty = TRUE
  )
)

clinician_counts <- clinician_lookup |>
  dplyr::count(
    specialty,
    experience_group,
    name = "clinicians_n"
  )

openxlsx::write.xlsx(
  list(
    "Supplementary Table S2" = subgroup_results,
    "Clinician Counts" = clinician_counts
  ),
  file = file.path(table_dir, "supplementary_subgroups.xlsx"),
  overwrite = TRUE
)

saveRDS(
  subgroup_results,
  file = file.path(derived_dir, "subgroup_results.rds")
)

message(
  "Exploratory subgroup analysis completed. ",
  "Experience-level models were adjusted for clinician specialty."
)
