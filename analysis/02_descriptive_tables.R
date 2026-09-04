# ============================================================
# 02_descriptive_tables.R
# Descriptive results and manuscript Tables 1-3
# ============================================================

if (!exists("table_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")

obj <- readRDS(objects_path)
analysis_all <- obj$analysis_all
analysis_paired <- obj$analysis_paired

presentation_level <- analysis_all |>
  dplyr::distinct(presentation_id, .keep_all = TRUE) |>
  dplyr::mutate(
    risk_group = factor(
      risk_colour,
      levels = 1:3,
      labels = unname(RISK_LABELS)
    )
  )

# ----------------------------
# Table 1: presentation characteristics by HSI risk group
# ----------------------------

make_table1_column <- function(dat) {
  n <- nrow(dat)

  c(
    format_median_iqr(dat$age, 1),
    format_median_iqr(dat$apd_us1, 1),
    format_n_percent(dat$sfu_us1 == 1, n),
    format_n_percent(dat$sfu_us1 == 2, n),
    format_n_percent(dat$sfu_us1 == 3, n),
    format_n_percent(dat$sfu_us1 == 4, n),
    format_n_percent(dat$renal_scan == 1, n),
    format_n_percent(dat$surgery == 1, n),
    format_median_iqr(dat$prediction_number, 3)
  )
}

table1 <- tibble::tibble(
  Characteristic = c(
    "Age, months",
    "APD on ultrasound, mm",
    "SFU grade 1",
    "SFU grade 2",
    "SFU grade 3",
    "SFU grade 4",
    "Renal scan performed",
    "Surgery during follow-up",
    "Displayed HSI-predicted probability of surgery"
  ),
  Overall = make_table1_column(presentation_level),
  `Green risk` = make_table1_column(
    dplyr::filter(presentation_level, risk_group == "Green")
  ),
  `Yellow risk` = make_table1_column(
    dplyr::filter(presentation_level, risk_group == "Yellow")
  ),
  `Red risk` = make_table1_column(
    dplyr::filter(presentation_level, risk_group == "Red")
  )
)

# ----------------------------
# Overall change summary
# ----------------------------

overall_change_summary <- analysis_paired |>
  dplyr::summarise(
    complete_pairs = dplyr::n(),
    changed_n = sum(changed),
    changed_percent = 100 * mean(changed),
    escalation_n = sum(change_direction == "Escalation"),
    escalation_percent = 100 * mean(change_direction == "Escalation"),
    deescalation_n = sum(change_direction == "De-escalation"),
    deescalation_percent = 100 * mean(change_direction == "De-escalation"),
    unchanged_n = sum(change_direction == "No change"),
    unchanged_percent = 100 * mean(change_direction == "No change")
  )

change_rate_ci <- stats::binom.test(
  x = sum(analysis_paired$changed),
  n = nrow(analysis_paired)
)$conf.int

# ----------------------------
# Table 2: 4 x 4 pre/post transition matrix
# ----------------------------

transition_counts <- table(
  analysis_paired$initial_plan,
  analysis_paired$post_plan
)

transition_row_percent <- 100 * prop.table(transition_counts, margin = 1)

table2_long <- as.data.frame(transition_counts) |>
  dplyr::rename(
    `Initial management plan` = Var1,
    `Post-HSI management plan` = Var2,
    n = Freq
  ) |>
  dplyr::mutate(
    `Initial management plan` = as.character(`Initial management plan`),
    `Post-HSI management plan` = as.character(`Post-HSI management plan`)
  ) |>
  dplyr::group_by(`Initial management plan`) |>
  dplyr::mutate(
    row_total = sum(n),
    row_percent = 100 * n / row_total,
    `n (% of row)` = sprintf("%d (%.1f%%)", n, row_percent)
  ) |>
  dplyr::ungroup()

table2 <- table2_long |>
  dplyr::select(
    `Initial management plan`,
    `Post-HSI management plan`,
    `n (% of row)`
  ) |>
  tidyr::pivot_wider(
    names_from = `Post-HSI management plan`,
    values_from = `n (% of row)`
  )

# ----------------------------
# Stuart-Maxwell test
# ----------------------------

stuart_maxwell <- DescTools::StuartMaxwellTest(transition_counts)

stuart_maxwell_result <- tibble::tibble(
  statistic = unname(stuart_maxwell$statistic),
  df = unname(stuart_maxwell$parameter),
  p_value = unname(stuart_maxwell$p.value)
)

# ----------------------------
# Overall recommendation distributions
# ----------------------------

overall_distribution <- analysis_paired |>
  dplyr::select(initial_plan, post_plan) |>
  tidyr::pivot_longer(
    cols = c(initial_plan, post_plan),
    names_to = "time",
    values_to = "management_plan"
  ) |>
  dplyr::mutate(
    time = dplyr::recode(
      time,
      initial_plan = "Pre-HSI",
      post_plan = "Post-HSI"
    )
  ) |>
  dplyr::count(time, management_plan, name = "n") |>
  dplyr::group_by(time) |>
  dplyr::mutate(percent = 100 * n / sum(n)) |>
  dplyr::ungroup()

# ----------------------------
# Table 3: risk-stratified directional shifts
# ----------------------------

risk_direction <- analysis_paired |>
  dplyr::group_by(risk_group) |>
  dplyr::summarise(
    n_pairs = dplyr::n(),
    changed_n = sum(changed),
    changed_percent = 100 * mean(changed),
    escalation_n = sum(change_direction == "Escalation"),
    escalation_percent = 100 * mean(change_direction == "Escalation"),
    deescalation_n = sum(change_direction == "De-escalation"),
    deescalation_percent = 100 * mean(change_direction == "De-escalation"),
    unchanged_n = sum(change_direction == "No change"),
    unchanged_percent = 100 * mean(change_direction == "No change"),
    .groups = "drop"
  )

risk_distribution <- analysis_paired |>
  dplyr::select(risk_group, initial_plan, post_plan) |>
  tidyr::pivot_longer(
    cols = c(initial_plan, post_plan),
    names_to = "time",
    values_to = "management_plan"
  ) |>
  dplyr::mutate(
    time = dplyr::recode(
      time,
      initial_plan = "Pre-HSI",
      post_plan = "Post-HSI"
    )
  ) |>
  dplyr::count(risk_group, time, management_plan, name = "n") |>
  dplyr::group_by(risk_group, time) |>
  dplyr::mutate(percent = 100 * n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::select(risk_group, time, management_plan, percent) |>
  tidyr::pivot_wider(
    names_from = time,
    values_from = percent
  ) |>
  dplyr::mutate(delta_pp = `Post-HSI` - `Pre-HSI`) |>
  dplyr::select(risk_group, management_plan, delta_pp) |>
  tidyr::pivot_wider(
    names_from = management_plan,
    values_from = delta_pp
  )

table3 <- risk_direction |>
  dplyr::left_join(risk_distribution, by = "risk_group") |>
  dplyr::transmute(
    `HSI risk group` = risk_group,
    `Complete pairs` = n_pairs,
    `Any change, n (%)` = sprintf("%d (%.1f%%)", changed_n, changed_percent),
    `Discharge Δ, pp` = sprintf("%+.1f", Discharge),
    `Repeat US Δ, pp` = sprintf("%+.1f", `Repeat ultrasound`),
    `Renogram Δ, pp` = sprintf("%+.1f", `Diuretic renogram`),
    `Surgical referral Δ, pp` = sprintf("%+.1f", `Surgical referral`),
    `Escalation, n (%)` = sprintf("%d (%.1f%%)", escalation_n, escalation_percent),
    `De-escalation, n (%)` = sprintf("%d (%.1f%%)", deescalation_n, deescalation_percent)
  )

# ----------------------------
# Export
# ----------------------------

openxlsx::write.xlsx(
  list(
    "Table 1" = table1,
    "Table 2" = table2,
    "Table 2 Long" = table2_long,
    "Table 3" = table3,
    "Overall Change Summary" = overall_change_summary,
    "Overall Change 95CI" = data.frame(
      lower_percent = 100 * change_rate_ci[1],
      upper_percent = 100 * change_rate_ci[2]
    ),
    "Overall Distribution" = overall_distribution,
    "Stuart-Maxwell" = stuart_maxwell_result
  ),
  file = file.path(table_dir, "descriptive_tables.xlsx"),
  overwrite = TRUE
)

saveRDS(
  list(
    table1 = table1,
    table2 = table2,
    table3 = table3,
    overall_change_summary = overall_change_summary,
    overall_distribution = overall_distribution,
    transition_counts = transition_counts,
    stuart_maxwell_result = stuart_maxwell_result
  ),
  file = file.path(derived_dir, "descriptive_results.rds")
)

message("Descriptive tables exported.")
