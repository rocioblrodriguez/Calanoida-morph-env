#!/usr/bin/env Rscript
# =============================================================================
# Script:  07_pca_two_clusters.R
# Author:  Rocio Rodriguez
# Purpose: Show the saved PCA with the independent Small/Large assignments.
# Style:   Original supplied plot construction, dimensions, and combinations.
#          Shared paths/colors, PNG/PDF/SVG export, and plain panel tags.
# =============================================================================

source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["PCA_2clusters"]]
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

POINT_SIZE <- 3.5
POINT_ALPHA <- 1
CENTROID_SIZE <- 7
TRAIT_LABEL_SCALE <- 14
TRAIT_LABEL_SIZE <- 4.5
PLOT_BASE_SIZE <- 16

COMBINED_BASE_SIZE <- 22
COMBINED_TRAIT_LABEL_SIZE <- 5.5
CLUSTER_LEGEND_POINT_SIZE <- 6
REGION_LEGEND_POINT_SIZE <- 5

cluster_levels <- c("S", "L")
cluster_cols <- c(S = "#FED976", L = "#800026")
cluster_labels <- c(S = "S = Small", L = "L = Large")
region_levels <- c("Subarctic", "Transition", "Subtropical")
region_shapes <- c(Subarctic = 21, Transition = 22, Subtropical = 23)

save_figure <- function(plot, name, width, height) {
  save_plot_formats(file.path(PLOT_DIR, paste0(name, ".png")), plot,
                    width = width, height = height, units = "in", dpi = FIGURE_DPI,
                    bg = "transparent", limitsize = FALSE)
}

results <- readRDS(file.path(OUT_DIR, "pca_results.rds"))
Cal_data <- results$data |>
  filter(!(Region == "Subarctic" & Cluster == "ST")) |>
  rename(Cluster_original = Cluster) |>
  mutate(Cluster = factor(Cluster2, levels = cluster_levels))
Cal_data_ordered <- arrange(Cal_data, abundance)
res.pca <- results$fit
pc_names <- paste0("PC", 1:4)
var_perc <- res.pca$eig[, 2]
names(var_perc) <- paste0("PC", seq_along(var_perc))
pc_label <- function(pc_name) paste0(pc_name, " (", round(var_perc[pc_name], 1), "%)")
PAIR_WIDTH <- 12
PAIR_HEIGHT <- 8
COMBINED_WIDTH <- 25
COMBINED_HEIGHT <- 16
table_directory <- file.path(OUT_DIR, "PCA_2clusters_SL_allRegions")
dir.create(table_directory, recursive = TRUE, showWarnings = FALSE)
cluster_summary <- Cal_data |>
  group_by(Cluster) |>
  summarise(n_observations = n(), total_abundance = sum(abundance),
            mean_esd_mm = mean(esd * process_particle_pixel_size_mm), .groups = "drop")
centroids <- Cal_data |>
  group_by(Cluster, Region) |>
  summarise(across(all_of(pc_names), ~ weighted.mean(.x, abundance)),
            total_abundance = sum(abundance), .groups = "drop")
loadings_all <- as.data.frame(res.pca$var$coord[, 1:4, drop = FALSE])
names(loadings_all) <- pc_names
trait_coordinates <- loadings_all[unname(TRAIT_PROXIES), , drop = FALSE]
rownames(trait_coordinates) <- names(TRAIT_PROXIES)
trait_coordinates["SMALL", ] <- -trait_coordinates["LARGE", ]
trait_coordinates$label <- rownames(trait_coordinates)

make_pooled_pca_pair <- function(x_pc, y_pc, combined = FALSE) {
  pair_traits <- trait_coordinates
  pair_traits$x <- pair_traits[[x_pc]] * TRAIT_LABEL_SCALE
  pair_traits$y <- pair_traits[[y_pc]] * TRAIT_LABEL_SCALE
  base_size <- if (combined) COMBINED_BASE_SIZE else PLOT_BASE_SIZE
  trait_size <- if (combined) COMBINED_TRAIT_LABEL_SIZE else TRAIT_LABEL_SIZE

  ggplot2::ggplot(Cal_data_ordered, ggplot2::aes(x = .data[[x_pc]], y = .data[[y_pc]])) +
    ggplot2::geom_point(
      ggplot2::aes(fill = Cluster, colour = Cluster, shape = Region),
      size = POINT_SIZE, alpha = POINT_ALPHA, stroke = 0.3
    ) +
    ggplot2::geom_point(
      data = centroids,
      ggplot2::aes(fill = Cluster, shape = Region),
      colour = "white", size = CENTROID_SIZE, stroke = 1.4,
      show.legend = FALSE
    ) +
    ggplot2::geom_label(
      data = pair_traits,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE, fill = "white", colour = "black",
      size = trait_size, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::scale_shape_manual(values = region_shapes, limits = region_levels, drop = FALSE) +
    ggplot2::scale_fill_manual(
      values = cluster_cols, limits = cluster_levels, labels = cluster_labels, drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = cluster_cols, limits = cluster_levels, labels = cluster_labels, drop = FALSE
    ) +
    ggplot2::labs(
      x = pc_label(x_pc), y = pc_label(y_pc),
      fill = "Cluster", colour = "Cluster", shape = "Region"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1, ncol = 2, byrow = TRUE,
        override.aes = list(shape = 21, size = CLUSTER_LEGEND_POINT_SIZE, colour = "black", alpha = 1)
      ),
      colour = "none",
      shape = ggplot2::guide_legend(
        order = 2, ncol = 3, byrow = TRUE,
        override.aes = list(size = REGION_LEGEND_POINT_SIZE, fill = "grey80", colour = "black", stroke = 1, alpha = 1)
      )
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (combined) "bottom" else "right",
      legend.box = "vertical",
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      legend.key = ggplot2::element_rect(fill = "transparent", colour = NA)
    )
}

# ---- SAVE SIX POOLED AXIS-PAIR PLOTS ------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
pc_pairs <- list(
  PC1_PC2 = c("PC1", "PC2"), PC1_PC3 = c("PC1", "PC3"),
  PC1_PC4 = c("PC1", "PC4"), PC2_PC3 = c("PC2", "PC3"),
  PC2_PC4 = c("PC2", "PC4"), PC3_PC4 = c("PC3", "PC4")
)
all_region_plots <- vector("list", length(pc_pairs))
names(all_region_plots) <- names(pc_pairs)
for (pair_name in names(pc_pairs)) {
  axes <- pc_pairs[[pair_name]]
  p <- make_pooled_pca_pair(axes[1], axes[2])
  all_region_plots[[pair_name]] <- p
  save_figure(
    p, paste0("pca_", pair_name, "_2clusters_SL_Calanoida_allRegions"),
    PAIR_WIDTH, PAIR_HEIGHT
  )
}

# ---- COMBINED SIX-PANEL FIGURE: ONE CLUSTER AND ONE REGION LEGEND -------------
combined_panels <- lapply(pc_pairs, function(axes) {
  make_pooled_pca_pair(axes[1], axes[2], combined = TRUE)
})
combined_all_regions <- (
  patchwork::wrap_plots(combined_panels, ncol = 3, nrow = 2) +
    patchwork::plot_layout(guides = "collect")
) & ggplot2::theme(
  legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical"
)
save_figure(
  combined_all_regions, "pca_2clusters_SL_Calanoida_allRegions_allPairs_combined",
  COMBINED_WIDTH, COMBINED_HEIGHT
)


write_csv(Cal_data, file.path(table_directory, "Calanoida_PCA_scores_2clusters_SL_allRegions.csv"))
write_csv(cluster_summary, file.path(table_directory, "cluster_summary_2clusters_SL_allRegions.csv"))
write_csv(centroids, file.path(table_directory, "abundance_weighted_centroids_Cluster_Region.csv"))
write_csv(results$eigenvalues, file.path(table_directory, "PCA_variance_2clusters_SL_allRegions.csv"))
