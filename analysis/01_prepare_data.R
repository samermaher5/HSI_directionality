# ============================================================
# 01_prepare_data.R
# Build the analysis dataset from the authorized source workbook
# ============================================================

if (!exists("data_path")) source("analysis/00_setup.R")

# ----------------------------
# Import
# ----------------------------

# Import this sheet as text so known annotation cells are handled explicitly
# rather than being silently coerced during readxl type guessing.
decision_raw <- readxl::read_excel(
  data_path,
  sheet = "Decision Making Data",
  col_types = "text"
) |>
  janitor::clean_names() |>
  dplyr::mutate(
    raw_decision_row = dplyr::row_number(),
    clinician_id = as.numeric(clinician_id),
    patient_case_id = as.numeric(dplyr::na_if(patient_case_id, "??")),
    initial_decision = as.numeric(dplyr::na_if(initial_decision, "?")),
    post_model_decision = as.numeric(
      dplyr::na_if(post_model_decision, "(Not on ppt)")
    )
  )

patient_raw <- readxl::read_excel(
  data_path,
  sheet = "Patient Variables"
) |>
  janitor::clean_names() |>
  dplyr::mutate(raw_patient_row = dplyr::row_number())

clinician_raw <- readxl::read_excel(
  data_path,
  sheet = "Demographic"
) |>
  janitor::clean_names()

required_decision_cols <- c(
  "clinician_id",
  "patient_case_id",
  "initial_decision",
  "post_model_decision"
)

required_patient_cols <- c(
  "patient_case_id",
  "prediction_number",
  "colour_1_green_2_yellow_3_red"
)

missing_decision_cols <- setdiff(required_decision_cols, names(decision_raw))
missing_patient_cols <- setdiff(required_patient_cols, names(patient_raw))

if (length(missing_decision_cols) > 0) {
  stop("Missing Decision Making Data columns: ",
       paste(missing_decision_cols, collapse = ", "))
}

if (length(missing_patient_cols) > 0) {
  stop("Missing Patient Variables columns: ",
       paste(missing_patient_cols, collapse = ", "))
}

# ----------------------------
# Exclude case 267
# ----------------------------

# Case 267 is excluded because no HSI output was shown.
decision_cleaning <- decision_raw |>
  dplyr::filter(
    !is.na(patient_case_id),
    patient_case_id != 267
  )

patient_cleaning <- patient_raw |>
  dplyr::filter(
    !is.na(patient_case_id),
    patient_case_id != 267
  )

# ----------------------------
# Verify review count and common sequence
# ----------------------------

rows_per_clinician <- decision_cleaning |>
  dplyr::count(clinician_id, name = "n_presentations")

assert_equal(
  dplyr::n_distinct(decision_cleaning$clinician_id),
  N_CLINICIANS_EXPECTED,
  "Number of clinicians"
)

if (!all(rows_per_clinician$n_presentations == N_PRESENTATIONS_EXPECTED)) {
  stop("Not every clinician has exactly 293 eligible presentations.")
}

decision_ordered <- decision_cleaning |>
  dplyr::arrange(clinician_id, raw_decision_row) |>
  dplyr::group_by(clinician_id) |>
  dplyr::mutate(presentation_order = dplyr::row_number()) |>
  dplyr::ungroup()

presentation_alignment <- decision_ordered |>
  dplyr::group_by(presentation_order) |>
  dplyr::summarise(
    n_clinicians = dplyr::n(),
    n_unique_case_ids = dplyr::n_distinct(patient_case_id),
    .groups = "drop"
  )

alignment_problems <- presentation_alignment |>
  dplyr::filter(
    n_clinicians != N_CLINICIANS_EXPECTED |
      n_unique_case_ids != 1
  )

if (nrow(alignment_problems) > 0) {
  stop("Clinicians did not review the same case sequence.")
}

# ----------------------------
# Create stable presentation IDs
# ----------------------------

# Original case IDs 2900 and 2901 each occur more than once but represent
# distinct case presentations. Preserve them using presentation-specific IDs.

presentation_roster <- decision_ordered |>
  dplyr::group_by(presentation_order) |>
  dplyr::summarise(
    patient_case_id = dplyr::first(patient_case_id),
    .groups = "drop"
  ) |>
  dplyr::group_by(patient_case_id) |>
  dplyr::mutate(occurrence_within_case = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    presentation_id = dplyr::case_when(
      patient_case_id == 2900 & occurrence_within_case == 1 ~ "2900_A",
      patient_case_id == 2900 & occurrence_within_case == 2 ~ "2900_B",
      patient_case_id == 2901 & occurrence_within_case == 1 ~ "2901_A",
      patient_case_id == 2901 & occurrence_within_case == 2 ~ "2901_B",
      TRUE ~ as.character(patient_case_id)
    )
  )

decision_clean <- decision_ordered |>
  dplyr::left_join(
    presentation_roster |>
      dplyr::select(
        presentation_order,
        patient_case_id,
        presentation_id
      ),
    by = c("presentation_order", "patient_case_id")
  )

patient_clean <- patient_cleaning |>
  dplyr::arrange(raw_patient_row) |>
  dplyr::group_by(patient_case_id) |>
  dplyr::mutate(occurrence_within_case = dplyr::row_number()) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    presentation_id = dplyr::case_when(
      patient_case_id == 2900 & occurrence_within_case == 1 ~ "2900_A",
      patient_case_id == 2900 & occurrence_within_case == 2 ~ "2900_B",
      patient_case_id == 2901 & occurrence_within_case == 1 ~ "2901_A",
      patient_case_id == 2901 & occurrence_within_case == 2 ~ "2901_B",
      TRUE ~ as.character(patient_case_id)
    )
  )

patient_duplicates <- patient_clean |>
  dplyr::count(presentation_id) |>
  dplyr::filter(n != 1)

if (nrow(patient_duplicates) > 0) {
  stop("At least one presentation_id has multiple Patient Variables rows.")
}

decision_duplicates <- decision_clean |>
  dplyr::count(clinician_id, presentation_id) |>
  dplyr::filter(n != 1)

if (nrow(decision_duplicates) > 0) {
  stop("At least one clinician-presentation combination is not unique.")
}

decision_without_patient <- decision_clean |>
  dplyr::distinct(presentation_id) |>
  dplyr::anti_join(
    patient_clean |> dplyr::distinct(presentation_id),
    by = "presentation_id"
  )

patient_without_decision <- patient_clean |>
  dplyr::distinct(presentation_id) |>
  dplyr::anti_join(
    decision_clean |> dplyr::distinct(presentation_id),
    by = "presentation_id"
  )

if (nrow(decision_without_patient) > 0 || nrow(patient_without_decision) > 0) {
  stop("Decision and patient-variable sheets do not fully align.")
}

# ----------------------------
# Join patient/model data
# ----------------------------

patient_variables <- patient_clean |>
  dplyr::rename(
    risk_colour = colour_1_green_2_yellow_3_red
  ) |>
  dplyr::select(
    dplyr::any_of(c(
      "presentation_id",
      "apd_us1",
      "sfu_us1",
      "surgery",
      "renal_scan",
      "prediction_number",
      "risk_colour",
      "duplicate_id",
      "age"
    ))
  )

analysis_all <- decision_clean |>
  dplyr::left_join(
    patient_variables,
    by = "presentation_id"
  )

assert_equal(
  nrow(analysis_all),
  N_ELIGIBLE_ROWS_EXPECTED,
  "Eligible clinician-presentation observations"
)

assert_equal(
  dplyr::n_distinct(analysis_all$presentation_id),
  N_PRESENTATIONS_EXPECTED,
  "Eligible case presentations"
)

# ----------------------------
# Complete paired analysis dataset
# ----------------------------

analysis_paired <- analysis_all |>
  dplyr::filter(
    !is.na(initial_decision),
    !is.na(post_model_decision),
    !is.na(prediction_number),
    !is.na(risk_colour)
  ) |>
  dplyr::mutate(
    clinician_id = factor(clinician_id),
    presentation_id = factor(presentation_id),
    pre_decision = as.integer(initial_decision),
    post_decision = as.integer(post_model_decision),
    initial_plan = factor(
      pre_decision,
      levels = 1:4,
      labels = unname(MANAGEMENT_LABELS)
    ),
    post_plan = factor(
      post_decision,
      levels = 1:4,
      labels = unname(MANAGEMENT_LABELS)
    ),
    changed = as.integer(post_decision != pre_decision),
    decision_change = post_decision - pre_decision,
    change_direction = factor(
      dplyr::case_when(
        decision_change > 0 ~ "Escalation",
        decision_change < 0 ~ "De-escalation",
        TRUE ~ "No change"
      ),
      levels = c("Escalation", "De-escalation", "No change")
    ),
    step_size = abs(decision_change),
    risk_group = factor(
      risk_colour,
      levels = 1:3,
      labels = unname(RISK_LABELS)
    ),
    hsi_prob_10 = prediction_number / 0.10
  )

assert_equal(
  nrow(analysis_paired),
  N_PAIRED_ROWS_EXPECTED,
  "Complete paired observations"
)

# ----------------------------
# Clinician lookup
# ----------------------------

clinician_lookup <- clinician_raw |>
  dplyr::mutate(
    clinician_id = factor(clinician_id),
    specialty = dplyr::case_when(
      specialty_uro_1_neph_2 == 1 ~ "Pediatric urology",
      specialty_uro_1_neph_2 == 2 ~ "Pediatric nephrology",
      TRUE ~ NA_character_
    ),
    designation_code = as.character(designation),
    experience_group = dplyr::case_when(
      designation_code %in% c("1", "3") ~ "Resident/novice APP",
      designation_code == "2" ~ "Fellow",
      designation_code %in% c("4", "5") ~ "Staff/expert APP",
      TRUE ~ NA_character_
    ),
    experience_group = factor(
      experience_group,
      levels = c(
        "Resident/novice APP",
        "Fellow",
        "Staff/expert APP"
      )
    )
  ) |>
  dplyr::select(
    dplyr::any_of(c(
      "clinician_id",
      "specialty",
      "experience_group",
      "designation",
      "years_of_experience"
    ))
  )

# ----------------------------
# Save local derived objects
# ----------------------------

saveRDS(
  list(
    analysis_all = analysis_all,
    analysis_paired = analysis_paired,
    presentation_roster = presentation_roster,
    patient_clean = patient_clean,
    clinician_raw = clinician_raw,
    clinician_lookup = clinician_lookup
  ),
  file = file.path(derived_dir, "analysis_objects.rds")
)

message(
  "Data preparation complete: ",
  nrow(analysis_all), " eligible rows; ",
  nrow(analysis_paired), " complete paired rows."
)
