# ============================================================
# optional_sankey.R
# Optional alluvial/Sankey display of the Table 2 transition matrix
# Not part of the core manuscript analysis unless explicitly used.
# ============================================================

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  stop("Package 'ggalluvial' is required for this optional figure.")
}

if (!exists("figure_dir")) source("analysis/00_setup.R")

objects_path <- file.path(derived_dir, "analysis_objects.rds")
if (!file.exists(objects_path)) source("analysis/01_prepare_data.R")

obj <- readRDS(objects_path)
analysis_paired <- obj$analysis_paired

transitions <- analysis_paired |>
  dplyr::count(
    Pre = initial_plan,
    Post = post_plan,
    name = "n"
  ) |>
  dplyr::mutate(
    Pre = factor(Pre, levels = unname(MANAGEMENT_LABELS)),
    Post = factor(Post, levels = unname(MANAGEMENT_LABELS))
  )

sankey <- ggplot2::ggplot(
  transitions,
  ggplot2::aes(
    axis1 = Pre,
    axis2 = Post,
    y = n
  )
) +
  ggalluvial::geom_alluvium(
    ggplot2::aes(fill = Pre),
    width = 0.12,
    alpha = 0.75
  ) +
  ggalluvial::geom_stratum(
    width = 0.15,
    fill = "white",
    colour = "black"
  ) +
  ggplot2::geom_text(
    stat = "stratum",
    ggplot2::aes(label = after_stat(stratum)),
    size = 4
  ) +
  ggplot2::scale_x_discrete(
    limits = c("Pre-HSI", "Post-HSI"),
    expand = c(0.08, 0.08)
  ) +
  ggplot2::labs(
    title = "Management transitions following HSI exposure",
    x = NULL,
    y = "Clinician-presentation observations",
    fill = "Pre-HSI recommendation"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(figure_dir, "Optional_Sankey.png"),
  sankey,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  file.path(figure_dir, "Optional_Sankey.pdf"),
  sankey,
  width = 11,
  height = 7,
  bg = "white"
)

message("Optional Sankey/alluvial figure exported.")
