#!/usr/bin/env Rscript
# =============================================================================
# Script:  05_pca_statistics.R
# Author:  Rocio Rodriguez
# Purpose: Fit one pooled PCA and describe PC1–PC4 and axis-retention diagnostics.
# Input:   outputs/data_var_clusters.csv
# Outputs: Calanoida_data.csv, pca_results.rds, eigenvalue/trait tables and plots.
# =============================================================================

# 1. Settings and packages -----------------------------------------------------
source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["PCA_statistics"]]
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(FactoMineR)
library(ggplot2)
library(ggtext)
library(patchwork)

# 2. Read clustered observations ----------------------------------------------
dat <- read_csv(file.path(OUT_DIR, "data_var_clusters.csv"), show_col_types = FALSE) |>
  filter(!(Region == "Subarctic" & Cluster == "ST")) |>
  mutate(Region = factor(Region, levels = REGION_LEVELS),
         Cluster = factor(Cluster, levels = CLUSTER_LEVELS))

# 3. Fit the pooled PCA once --------------------------------------------------
# Preserve the original method: multiply each morphology row by sqrt(abundance),
# then center and scale columns in FactoMineR. This is a transformation of the
# response matrix, distinct from supplying observation weights to PCA.
weighted_morphology <- sweep(as.matrix(dat[MORPH_VARIABLES]), 1, sqrt(dat$abundance), FUN = "*")
pca <- PCA(as.data.frame(weighted_morphology), scale.unit = TRUE,
           ncp = length(MORPH_VARIABLES), graph = FALSE)
eigenvalues <- as.data.frame(pca$eig)
names(eigenvalues) <- c("eigenvalue", "percent_var", "cumulative_percent")
eigenvalues$PC <- seq_len(nrow(eigenvalues))

# Retain Kaiser axes for the exported dataset and always include PC1–PC4 for
# the figures. Eigenvalue > 1 is a retention rule, not a significance test.
retained_axes <- sort(unique(c(1:4, which(eigenvalues$eigenvalue > 1))))
scores <- as.data.frame(pca$ind$coord[, retained_axes, drop = FALSE])
names(scores) <- paste0("PC", retained_axes)
dat <- bind_cols(dat, scores)
write_csv(dat, file.path(OUT_DIR, "Calanoida_data.csv"))
write_csv(eigenvalues, file.path(OUT_DIR, "pca_eigenvalues.csv"))
saveRDS(list(fit = pca, data = dat, eigenvalues = eigenvalues,
             retained_axes = retained_axes, variables = MORPH_VARIABLES),
        file.path(OUT_DIR, "pca_results.rds"))

# 4. Original axis diagnostics and trait heatmaps -----------------------------
res.pca <- pca
eig_df <- eigenvalues
kaiser_percent <- 100 / length(MORPH_VARIABLES)
PCS_TO_PLOT <- 1:4
# ---- 2a) Eigenvalue plot ----
p_eigen <- ggplot(eig_df, aes(x = PC, y = eigenvalue)) +
  geom_point() +
  geom_line() +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "Principal Component",
    y = "Eigenvalue"
  ) +
  theme_classic()

save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_eigen.png"),
  plot     = p_eigen,
  width    = 5,
  height   = 4,
  dpi      = 300
)

# ---- 2b) Percent variance plot ----
p_var <- ggplot(eig_df, aes(x = PC, y = percent_var)) +
  geom_point() +
  geom_line() +
  geom_hline(
    yintercept = kaiser_percent,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "Principal Component",
    y = "% Variance"
  ) +
  theme_classic()

save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_variance_kaiser.png"),
  plot     = p_var,
  width    = 5,
  height   = 4,
  dpi      = 300
)

# ---- 2c) Cumulative variance plot ----
p_cumvar <- ggplot(eig_df, aes(x = PC, y = cumulative_percent)) +
  geom_point() +
  geom_line() +
  geom_hline(
    yintercept = 90,
    linetype = "dashed",
    color = "red"
  ) +
  labs(
    x = "Principal Component",
    y = "% Cumulative Variance"
  ) +
  theme_classic()

save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_p_cumvar.png"),
  plot     = p_cumvar,
  width    = 5,
  height   = 4,
  dpi      = 300
)

# ---- 2d) Combined significant-axis plot, horizontal with tags ----
significant_axis_plot <- 
  p_eigen + p_var + p_cumvar +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(
    tag_levels = "A",
    theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain"))
  )

save_plot_formats(
  filename = file.path(PLOT_DIR, "significant_axis.png"),
  plot     = significant_axis_plot,
  width    = 15,
  height   = 4.5,
  dpi      = 300
)

# ---- 2e) Combined significant-axis plot, one column with tags ----
significant_axis_plot_column <-
  p_eigen + p_var + p_cumvar +
  patchwork::plot_layout(ncol = 1) +
  patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain")))

save_plot_formats(
  filename = file.path(PLOT_DIR, "significant_axis_column.png"),
  plot     = significant_axis_plot_column,
  width    = 5,
  height   = 12,
  dpi      = 300,
  bg       = "white"
)

# -----------------------------------------------------------------------------
# 3) Coordinates / cos2 / contributions: ordered, colored labels; PC1-PC4
# -----------------------------------------------------------------------------

# Use the shared category assignments, in the original displayed order.
morph_plot_groups <- lapply(MORPH_GROUPS, function(x) paste0(x, "_YJ"))
morph_plot_groups$Transparency <- c(
  "mean_YJ", "stddev_YJ", "mode_YJ", "min_YJ", "max_YJ", "intden_YJ",
  "median_YJ", "skew_YJ", "kurt_YJ", "cv_YJ", "sr_YJ", "histcum1_YJ",
  "histcum2_YJ", "histcum3_YJ", "nb1_YJ", "meanpos_YJ", "nb2_YJ"
)
dim_keep <- paste0("Dim.", PCS_TO_PLOT)
variable_key <- data.frame(
  variable = unlist(morph_plot_groups, use.names = FALSE),
  morph_group = rep(names(morph_plot_groups), lengths(morph_plot_groups)),
  stringsAsFactors = FALSE
)

# Keep the same variable order in every panel, even when values are small.
variable_key <- variable_key |>
  dplyr::filter(variable %in% rownames(res.pca$var$coord)) |>
  dplyr::mutate(
    label = sub("_YJ$", "", variable),
    label_color = unname(MORPH_COLORS[morph_group])
  )

plot_variables <- variable_key$variable

# ggtext renders colored axis labels; the _YJ suffix is hidden only in the plots.
colored_labels <- stats::setNames(
  sprintf(
    "<span style='color:%s;'>%s</span>",
    variable_key$label_color,
    variable_key$label
  ),
  variable_key$variable
)

group_legend <- paste(
  sprintf(
    "<span style='color:%s;'>%s</span>",
    unname(MORPH_COLORS[names(morph_plot_groups)]),
    names(morph_plot_groups)
  ),
  collapse = " &nbsp;|&nbsp; "
)

pc_labels <- stats::setNames(
  paste0("PC", PCS_TO_PLOT),
  dim_keep
)

# Horizontal lines separate categories, accounting for ggplot's bottom-up y-axis.
group_sizes <- vapply(
  names(morph_plot_groups),
  function(group) sum(variable_key$morph_group == group),
  integer(1)
)
group_sizes <- group_sizes[group_sizes > 0]
group_boundaries <- head(cumsum(group_sizes), -1)
separator_y <- length(plot_variables) - group_boundaries + 0.5

make_variable_heatmap <- function(values, title, fill_title,
                                  palette, fill_limits, direction = 1) {
  # Subset before pivoting: only the requested PCs enter the plot data.
  plot_data <- as.data.frame(values[, dim_keep, drop = FALSE])
  plot_data$variable <- rownames(plot_data)
  
  plot_data <- plot_data |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(dim_keep),
      names_to = "PC",
      values_to = "value"
    ) |>
    dplyr::mutate(
      PC = factor(PC, levels = dim_keep),
      variable = factor(variable, levels = rev(plot_variables))
    )
  
  # Keep numbers readable on the darker ends of each palette.
  if (palette == "RdBu") {
    plot_data$text_color <- ifelse(abs(plot_data$value) > 0.7, "white", "grey10")
  } else {
    plot_data$text_color <- ifelse(
      plot_data$value > fill_limits[1] + 0.65 * diff(fill_limits),
      "white", "grey10"
    )
  }
  
  p <- ggplot(plot_data, aes(x = PC, y = variable, fill = value)) +
    geom_tile(color = "white", linewidth = 0.25) +
    geom_text(aes(label = sprintf("%.2f", value), color = text_color), size = 4) +
    scale_color_identity(guide = "none") +
    geom_hline(yintercept = separator_y, color = "grey55", linewidth = 0.5) +
    scale_x_discrete(
      limits = dim_keep,
      labels = pc_labels,
      position = "bottom",
      expand = expansion(mult = 0, add = 0)
    ) +
    scale_y_discrete(
      limits = rev(plot_variables),
      labels = function(x) unname(colored_labels[x]),
      expand = expansion(mult = 0, add = 0)
    ) +
    labs(
      title = title,
      subtitle = group_legend,
      x = NULL,
      y = NULL,
      fill = fill_title
    ) +
    guides(
      fill = guide_colorbar(
        direction = "horizontal", title.position = "top"
      )
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(size = 11, face = "bold"),
      axis.text.y = ggtext::element_markdown(size = 10, margin = margin(r = 8)),
      # Set the rendered left-axis element explicitly. A more specific
      # theme element can otherwise display the HTML tags as plain text.
      axis.text.y.left = ggtext::element_markdown(
        size = 10, hjust = 1, margin = margin(r = 8)
      ),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = ggtext::element_markdown(size = 10, margin = margin(b = 10)),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 12, 10, 10)
    )
  
  # Each metric has its own continuous palette and appropriate numeric range.
  p + scale_fill_distiller(
    palette = palette,
    direction = direction,
    limits = fill_limits
  )
}

# ---- 3a-c) Three panels with identical axes, labels, and row order ----
p_coord <- make_variable_heatmap(
  res.pca$var$coord, "Variable coordinates", "Correlation",
  palette = "RdBu", fill_limits = c(-1, 1), direction = -1
)
p_cos2 <- make_variable_heatmap(
  res.pca$var$cos2, "Variable cos²", "cos²",
  palette = "Greens", fill_limits = c(0, 1)
)

# Plot original percentages; set the upper limit from the PCs actually displayed.
contrib_max <- max(res.pca$var$contrib[, dim_keep, drop = FALSE], na.rm = TRUE)

p_contrib <- make_variable_heatmap(
  res.pca$var$contrib, "Variable contributions", "Contribution (%)",
  palette = "Purples", fill_limits = c(0, contrib_max)
)

heatmap_height <- max(10, 0.22 * length(plot_variables) + 2)

save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_coord_heatmap.png"),
  plot = p_coord, width = 8, height = heatmap_height, dpi = 300, bg = "white"
)
save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_cos2_heatmap.png"),
  plot = p_cos2, width = 8, height = heatmap_height, dpi = 300, bg = "white"
)
save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_contrib_heatmap.png"),
  plot = p_contrib, width = 8, height = heatmap_height, dpi = 300, bg = "white"
)

# ---- 3d) Combined panels: coordinates, cos2, then contributions ----
morph_variable_pca_plot <-
  p_coord + p_cos2 + p_contrib +
  patchwork::plot_layout(ncol = 3) +
  patchwork::plot_annotation(tag_levels = "A", theme = ggplot2::theme(plot.tag = ggplot2::element_text(face = "plain")))

save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_morphological_variable_contributions_combined.png"),
  plot = morph_variable_pca_plot,
  width = 20, height = 12, dpi = 300, bg = "white"
)


# Save the metric tables alongside the common fitted PCA.
for (metric in c("coord", "cos2", "contrib")) {
  values <- as.data.frame(pca$var[[metric]][, 1:4, drop = FALSE]) |>
    rownames_to_column("variable") |>
    pivot_longer(-variable, names_to = "PC", values_to = "value") |>
    mutate(PC = sub("Dim.", "PC", PC, fixed = TRUE)) |>
    left_join(variable_key, by = "variable") |>
    rename(category = morph_group, color = label_color)
  suffix <- c(coord = "coordinates", cos2 = "cos2", contrib = "contributions")[[metric]]
  write_csv(values, file.path(OUT_DIR, paste0("pca_trait_", suffix, ".csv")))
}
write_csv(dat |>
            group_by(Region, Cluster) |>
            summarise(
              n_rows = n(),
              n_D = sum(D_N == "D", na.rm = TRUE),
              n_N = sum(D_N == "N", na.rm = TRUE),
              n_D_minus_N = n_D - n_N,
              sum_abundance = sum(abundance, na.rm = TRUE),
              abundance_D = sum(abundance[D_N == "D"], na.rm = TRUE),
              abundance_N = sum(abundance[D_N == "N"], na.rm = TRUE),
              abundance_D_minus_N = abundance_D - abundance_N,
              mean_abundance = mean(abundance, na.rm = TRUE),
              .groups = "drop"
            ),
          file.path(OUT_DIR, "cluster_rows_per_region.csv"))
