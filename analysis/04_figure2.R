# ============================================================
# 04_figure2.R
# Figure 2: continuous dose-response + categorical risk estimates
# ============================================================

if (!exists("figure_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
models_path <- file.path(derived_dir, "model_results.rds")

if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")
if (!file.exists(models_path)) source("analysis/03_mixed_models.R")

obj <- readRDS(objects_path)
mod <- readRDS(models_path)

analysis_all <- obj$analysis_all
analysis_paired <- obj$analysis_paired
adjusted_probs_long <- mod$adjusted_probs_long

figure2_base <- analysis_paired |>
  dplyr::mutate(
    escalation = as.integer(post_decision > pre_decision),
    deescalation = as.integer(post_decision < pre_decision)
  )

data_escalation_fig <- figure2_base |>
  dplyr::filter(pre_decision < 4) |>
  droplevels()

data_deescalation_fig <- figure2_base |>
  dplyr::filter(pre_decision > 1) |>
  droplevels()

presentation_scores <- analysis_all |>
  dplyr::distinct(
    presentation_id,
    prediction_number,
    risk_colour
  ) |>
  dplyr::mutate(
    risk_group = factor(
      risk_colour,
      levels = 1:3,
      labels = unname(RISK_LABELS)
    )
  )

# ----------------------------
# Spline mixed-effects models
# ----------------------------

model_escalation_spline <- lme4::glmer(
  escalation ~
    splines::ns(prediction_number, df = 3) +
    initial_plan +
    (1 | clinician_id) +
    (1 | presentation_id),
  data = data_escalation_fig,
  family = stats::binomial,
  control = glmer_control
)

model_deescalation_spline <- lme4::glmer(
  deescalation ~
    splines::ns(prediction_number, df = 3) +
    initial_plan +
    (1 | clinician_id) +
    (1 | presentation_id),
  data = data_deescalation_fig,
  family = stats::binomial,
  control = glmer_control
)

make_adjusted_curve <- function(
    model,
    model_data,
    outcome_label,
    n_grid = 200,
    n_draws = 2000,
    seed = 2026
) {

  risk_grid <- seq(
    min(model_data$prediction_number, na.rm = TRUE),
    max(model_data$prediction_number, na.rm = TRUE),
    length.out = n_grid
  )

  plan_weights <- model_data |>
    dplyr::count(initial_plan, name = "n") |>
    dplyr::mutate(weight = n / sum(n)) |>
    dplyr::select(initial_plan, weight)

  prediction_grid <- tidyr::crossing(
    prediction_number = risk_grid,
    initial_plan = plan_weights$initial_plan
  ) |>
    dplyr::left_join(plan_weights, by = "initial_plan") |>
    dplyr::mutate(
      initial_plan = factor(
        initial_plan,
        levels = levels(model_data$initial_plan)
      )
    )

  fixed_terms <- stats::delete.response(
    stats::terms(model, fixed.only = TRUE)
  )

  X <- stats::model.matrix(fixed_terms, data = prediction_grid)

  beta_hat <- lme4::fixef(model)
  X <- X[, names(beta_hat), drop = FALSE]

  set.seed(seed)
  beta_draws <- MASS::mvrnorm(
    n = n_draws,
    mu = beta_hat,
    Sigma = as.matrix(stats::vcov(model))
  )

  point_prob <- stats::plogis(X %*% beta_hat)
  sim_prob <- stats::plogis(X %*% t(beta_draws))

  prediction_grid |>
    dplyr::mutate(
      row_id = dplyr::row_number(),
      point_probability = as.numeric(point_prob)
    ) |>
    dplyr::group_by(prediction_number) |>
    dplyr::group_modify(
      ~ {
        rows <- .x$row_id
        weighted_sim <- sweep(
          sim_prob[rows, , drop = FALSE],
          1,
          .x$weight,
          "*"
        )
        adjusted_sim <- colSums(weighted_sim)

        tibble::tibble(
          adjusted_probability = sum(.x$point_probability * .x$weight),
          ci_lower = stats::quantile(adjusted_sim, 0.025),
          ci_upper = stats::quantile(adjusted_sim, 0.975)
        )
      }
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(Outcome = outcome_label)
}

figure2_curve_data <- dplyr::bind_rows(
  make_adjusted_curve(
    model_escalation_spline,
    data_escalation_fig,
    "Escalation"
  ),
  make_adjusted_curve(
    model_deescalation_spline,
    data_deescalation_fig,
    "De-escalation"
  )
) |>
  dplyr::mutate(
    Outcome = factor(
      Outcome,
      levels = c("Escalation", "De-escalation")
    )
  )

# ----------------------------
# Panel A: continuous HSI probability
# ----------------------------

outcome_colours <- c(
  "Escalation" = "#D55E00",
  "De-escalation" = "#0072B2"
)

panel_a <- ggplot2::ggplot(
  figure2_curve_data,
  ggplot2::aes(
    x = prediction_number,
    y = adjusted_probability,
    colour = Outcome,
    linetype = Outcome
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = ci_lower,
      ymax = ci_upper,
      fill = Outcome
    ),
    alpha = 0.12,
    colour = NA,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(linewidth = 1.2) +
  ggplot2::geom_rug(
    data = presentation_scores,
    ggplot2::aes(x = prediction_number),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.14,
    colour = "grey35",
    linewidth = 0.38
  ) +
  ggplot2::scale_colour_manual(values = outcome_colours) +
  ggplot2::scale_fill_manual(values = outcome_colours) +
  ggplot2::scale_linetype_manual(
    values = c("Escalation" = "solid", "De-escalation" = "dashed")
  ) +
  ggplot2::scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = ggplot2::expansion(mult = c(0.02, 0.03))
  ) +
  ggplot2::labs(
    x = "Displayed HSI-predicted probability of surgery",
    y = "Adjusted probability of management revision",
    colour = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "top"
  )

# ----------------------------
# Panel B: categorical HSI risk group
# ----------------------------

panel_b_data <- adjusted_probs_long |>
  dplyr::filter(
    Outcome %in% c("Escalation", "De-escalation")
  ) |>
  dplyr::mutate(
    risk_group = factor(
      risk_group,
      levels = c("Green", "Yellow", "Red")
    ),
    Outcome = factor(
      Outcome,
      levels = c("Escalation", "De-escalation")
    ),
    x_num = as.numeric(risk_group),
    x_plot = dplyr::if_else(
      Outcome == "Escalation",
      x_num - 0.06,
      x_num + 0.06
    )
  )

# Light risk-category backgrounds are supplementary visual cues only;
# category names and outcome line types ensure the figure does not rely
# exclusively on red/green colour discrimination.
panel_b <- ggplot2::ggplot() +
  ggplot2::annotate(
    "rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
    fill = "#D9F2D9", alpha = 0.30
  ) +
  ggplot2::annotate(
    "rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
    fill = "#FFF4CC", alpha = 0.30
  ) +
  ggplot2::annotate(
    "rect", xmin = 2.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
    fill = "#F9D6D5", alpha = 0.30
  ) +
  ggplot2::geom_line(
    data = panel_b_data,
    ggplot2::aes(
      x = x_plot,
      y = probability,
      colour = Outcome,
      linetype = Outcome,
      group = Outcome
    ),
    linewidth = 1.1,
    show.legend = FALSE
  ) +
  ggplot2::geom_errorbar(
    data = panel_b_data,
    ggplot2::aes(
      x = x_plot,
      ymin = ci_lower,
      ymax = ci_upper,
      colour = Outcome
    ),
    width = 0.05,
    linewidth = 0.85,
    show.legend = FALSE
  ) +
  ggplot2::geom_point(
    data = panel_b_data,
    ggplot2::aes(
      x = x_plot,
      y = probability,
      colour = Outcome,
      shape = Outcome
    ),
    size = 3.1,
    show.legend = FALSE
  ) +
  ggplot2::scale_colour_manual(values = outcome_colours) +
  ggplot2::scale_linetype_manual(
    values = c("Escalation" = "solid", "De-escalation" = "dashed")
  ) +
  ggplot2::scale_shape_manual(
    values = c("Escalation" = 16, "De-escalation" = 17)
  ) +
  ggplot2::scale_x_continuous(
    breaks = 1:3,
    labels = c("Green", "Yellow", "Red"),
    limits = c(0.5, 3.5)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = ggplot2::expansion(mult = c(0.02, 0.03))
  ) +
  ggplot2::labs(
    x = "HSI risk category",
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 12.5) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "none"
  )

figure2 <- panel_a + panel_b +
  patchwork::plot_layout(
    widths = c(1.15, 0.95),
    guides = "collect"
  ) +
  patchwork::plot_annotation(
    title = "Relationship between HSI output and directional management revision",
    tag_levels = "A",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(face = "bold", size = 14)
    )
  ) &
  ggplot2::theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(figure_dir, "Figure_2.png"),
  plot = figure2,
  width = 13.5,
  height = 7.2,
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(figure_dir, "Figure_2.pdf"),
  plot = figure2,
  width = 13.5,
  height = 7.2,
  bg = "white"
)

openxlsx::write.xlsx(
  list(
    "Panel A Continuous Curves" = figure2_curve_data,
    "Panel B Categorical Data" = panel_b_data,
    "Presentation Scores" = presentation_scores
  ),
  file = file.path(table_dir, "Figure_2_source_data.xlsx"),
  overwrite = TRUE
)

message("Figure 2 exported.")
