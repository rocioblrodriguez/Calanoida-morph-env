#!/usr/bin/env Rscript
# =============================================================================
# Script:  06_pca_plots.R
# Author:  Rocio Rodriguez
# Purpose: Plot the saved PCA with six morphotypes, pooled and regional.
# Style:   Original supplied plot construction, dimensions, and combinations.
#          Shared paths/colors, PNG/PDF/SVG export, and plain panel tags.
# =============================================================================

source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["PCA_6clusters"]]
library(dplyr)
library(readr)
library(ggplot2)
library(ggimage)
library(patchwork)

results <- readRDS(file.path(OUT_DIR, "pca_results.rds"))
Cal_data <- results$data |>
  dplyr::filter(!(Region == "Subarctic" & Cluster == "ST"))
res.pca <- results$fit
IN_CSV <- file.path(OUT_DIR, "Calanoida_data.csv")
cluster_cols <- CLUSTER_COLORS
cluster_icons <- CLUSTER_ICONS
cluster_labels <- CLUSTER_LABELS[c("ST", "MT", "LT", "SD", "MD", "LD")]
region_shapes <- REGION_SHAPES
regions <- reg_order <- REGION_LEVELS
icon_size_min <- 0.08
icon_size_max <- 0.20
var_perc <- res.pca$eig[, 2]
names(var_perc) <- paste0("PC", seq_along(var_perc))
pc_label <- function(pc_name) paste0(pc_name, " (", round(var_perc[pc_name], 1), "%)")
loadings_all <- as.data.frame(res.pca$var$coord[, 1:4, drop = FALSE])
composite_map <- as.list(TRAIT_PROXIES)
pca.vars1 <- loadings_all[unlist(composite_map), , drop = FALSE]
rownames(pca.vars1) <- names(composite_map)
pca.vars1 <- rbind(pca.vars1, SMALL = -pca.vars1["LARGE", ])
pca.vars1$label <- rownames(pca.vars1)
label_scale <- 18
Cal_data_ordered <- arrange(Cal_data, abundance)

# Original icons are scaled using mean raw ESD across all cluster/region nodes.
icon_metadata <- Cal_data |>
  group_by(Cluster, Region) |>
  summarise(mean_esd = weighted.mean(esd, abundance), .groups = "drop") |>
  mutate(icon = unname(cluster_icons[as.character(Cluster)]),
         icon_size = icon_sizes(mean_esd, icon_size_min, icon_size_max))
centroids <- Cal_data |>
  group_by(Region, Cluster) |>
  summarise(across(all_of(paste0("PC", 1:4)), ~ weighted.mean(.x, abundance)),
            .groups = "drop") |>
  left_join(icon_metadata, by = c("Region", "Cluster"))
write_csv(centroids, file.path(OUT_DIR, "pca_centroids_cluster_region.csv"))
all_region_panels <- list()

save_region_pca_pair <- function(x_dim, y_dim) {
  x_pc <- paste0("PC", x_dim)
  y_pc <- paste0("PC", y_dim)
  pair_name <- paste(x_pc, y_pc, sep = "_")
  
  # Use the loadings for the selected axes to position trait labels.
  pair_traits <- pca.vars1
  pair_traits$x <- pair_traits[[paste0("Dim.", x_dim)]] * label_scale
  pair_traits$y <- pair_traits[[paste0("Dim.", y_dim)]] * label_scale
  
  # Recalculate centroids on the selected axes; reuse icon sizes based on ESD.
  pair_means <- Cal_data %>%
    dplyr::group_by(Cluster, Region) %>%
    dplyr::summarise(
      wmean_x = stats::weighted.mean(.data[[x_pc]], w = abundance, na.rm = TRUE),
      wmean_y = stats::weighted.mean(.data[[y_pc]], w = abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      icon_metadata %>% dplyr::select(Cluster, Region, icon, icon_size),
      by = c("Cluster", "Region")
    )
  
  # All regions together: color by Cluster and shape by Region.
  p_all_pair <- ggplot(Cal_data_ordered, aes(x = .data[[x_pc]], y = .data[[y_pc]])) +
    geom_point(aes(fill = Cluster, shape = Region, color = Cluster), size = 2, stroke = 0.3) +
    geom_point(
      data = pair_means,
      aes(x = wmean_x, y = wmean_y, fill = Cluster, shape = Region),
      color = "white", size = 5, stroke = 1.4, show.legend = FALSE
    ) +
    geom_label(
      data = pair_traits,
      aes(x = x, y = y, label = label),
      fill = "white", color = "black", size = 4.5, fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_shape_manual(values = region_shapes) +
    scale_fill_manual(values = cluster_cols, breaks = names(cluster_labels), labels = cluster_labels) +
    scale_color_manual(values = cluster_cols, breaks = names(cluster_labels), labels = cluster_labels) +
    theme_classic(base_size = 12) +
    labs(
      title = NULL, subtitle = NULL,
      x = pc_label(x_pc), y = pc_label(y_pc),
      fill = "Cluster", color = "Cluster", shape = "Region"
    ) +
    guides(
      fill = guide_legend(order = 1, ncol = 3, byrow = TRUE, override.aes = list(size = 5)),
      color = "none",
      shape = guide_legend(order = 2, override.aes = list(size = 4, stroke = 1))
    ) +
    theme(
      legend.position   = "right",
      plot.background   = element_rect(fill = "transparent", color = NA),
      panel.background  = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key        = element_rect(fill = "transparent", color = NA),
      plot.title        = element_blank(),
      plot.subtitle     = element_blank()
    )
  
  save_plot_formats(
    filename = file.path(PLOT_DIR, paste0("pca_", pair_name, "_clusters_Calanoida_allRegions.png")),
    plot = p_all_pair, width = 12, height = 8, dpi = 300
  )
  all_region_panels[[pair_name]] <<- p_all_pair
  
  plots_by_region <- list()
  
  for (reg in regions) {
    cal_reg <- dplyr::filter(Cal_data, Region == reg)
    
    reg_means <- dplyr::filter(pair_means, Region == reg)
    
    p_reg_pair <- ggplot(cal_reg, aes(x = .data[[x_pc]], y = .data[[y_pc]])) +
      geom_point(
        aes(fill = Cluster), shape = 21, size = 2.5, stroke = 0.1,
        color = "black", show.legend = TRUE
      ) +
      ggimage::geom_image(
        data = reg_means,
        aes(x = wmean_x, y = wmean_y, image = icon, size = icon_size),
        asp = 1, inherit.aes = FALSE, show.legend = FALSE
      ) +
      scale_size_identity() +
      scale_fill_manual(
        values = cluster_cols, limits = names(cluster_labels),
        breaks = names(cluster_labels), labels = cluster_labels, drop = FALSE
      ) +
      geom_label(
        data = pair_traits,
        aes(x = x, y = y, label = label),
        fill = "white", color = "black", size = 4.5, fontface = "bold",
        inherit.aes = FALSE
      ) +
      scale_x_continuous(limits = c(-25, 20)) +
      scale_y_continuous(limits = c(-25, 25)) +
      labs(
        title = NULL, subtitle = NULL,
        x = pc_label(x_pc), y = pc_label(y_pc), fill = "Cluster"
      ) +
      theme_classic(base_size = 20) +
      guides(fill = guide_legend(override.aes = list(size = 5), ncol = 3, byrow = TRUE)) +
      theme(
        panel.border      = element_rect(colour = "black", size = 1.5, fill = "transparent"),
        legend.position   = "bottom",
        axis.line         = element_blank(),
        plot.background   = element_rect(fill = "transparent", color = NA),
        panel.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key        = element_rect(fill = "transparent", color = NA),
        plot.title        = element_blank(),
        plot.subtitle     = element_blank()
      )
    
    plots_by_region[[as.character(reg)]] <- p_reg_pair
    
    reg_clean <- gsub(" ", "", as.character(reg))
    save_plot_formats(
      filename = file.path(PLOT_DIR, paste0("pca_", pair_name, "_clusters_Calanoida_", reg_clean, ".png")),
      plot = p_reg_pair, width = 10, height = 8, dpi = 300
    )
  }
  
  reg_present_pair <- reg_order[reg_order %in% names(plots_by_region)]
  if (length(reg_present_pair) > 0) {
    plots_to_combine_pair <- lapply(reg_present_pair, function(r) plots_by_region[[r]])
    
    # Vertical panels: one shared legend on the right, with six stacked keys.
    plots_column <- lapply(plots_to_combine_pair, function(p) {
      p + ggplot2::guides(
        fill = ggplot2::guide_legend(
          ncol = 1, byrow = TRUE, override.aes = list(size = 5)
        )
      )
    })
    
    combined_pair <- (
      patchwork::wrap_plots(plots_column, ncol = 1) +
        patchwork::plot_layout(guides = "collect")
    ) & ggplot2::theme(
      legend.position = "right",
      legend.direction = "vertical"
    )
    
    save_plot_formats(
      filename = file.path(PLOT_DIR, paste0("pca_", pair_name, "_clusters_Calanoida_3regions_combined.png")),
      plot = combined_pair, width = 12, height = 20, dpi = 300
    )
    
    # Horizontal panels: one shared bottom legend, three columns by two rows.
    plots_row <- lapply(plots_to_combine_pair, function(p) {
      p + ggplot2::guides(
        fill = ggplot2::guide_legend(
          ncol = 3, nrow = 2, byrow = TRUE, override.aes = list(size = 5)
        )
      )
    })
    
    combined_pair_row <- (
      patchwork::wrap_plots(plots_row, nrow = 1) +
        patchwork::plot_layout(guides = "collect")
    ) & ggplot2::theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
    
    save_plot_formats(
      filename = file.path(PLOT_DIR, paste0("pca_", pair_name, "_clusters_Calanoida_3regions_combined_row.png")),
      plot = combined_pair_row, width = 25, height = 8, dpi = 300
    )
  }
  
  return(plots_by_region)
}

pca_plot_list_by_region_PC1_PC2 <- save_region_pca_pair(1, 2)
pca_plot_list_by_region_PC3_PC4 <- save_region_pca_pair(3, 4)
pca_plot_list_by_region_PC1_PC4 <- save_region_pca_pair(1, 4)
pca_plot_list_by_region_PC1_PC3 <- save_region_pca_pair(1, 3)
pca_plot_list_by_region_PC2_PC3 <- save_region_pca_pair(2, 3)
pca_plot_list_by_region_PC2_PC4 <- save_region_pca_pair(2, 4)

# -----------------------------------------------------------------------------
# FINAL COMBINED FIGURE: all-region plots for all six PC pairs
# -----------------------------------------------------------------------------
all_region_pair_order <- c(
  "PC1_PC2", "PC1_PC3", "PC1_PC4",
  "PC2_PC3", "PC2_PC4", "PC3_PC4"
)

all_region_plots <- lapply(all_region_pair_order, function(pair_name) {
  p <- all_region_panels[[pair_name]]
  
  p +
    ggplot2::labs(title = NULL, subtitle = NULL, tag = NULL) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        order = 1, ncol = 3, nrow = 2, byrow = TRUE,
        override.aes = list(shape = 21, size = 5, colour = "black")
      ),
      color = "none",
      shape = ggplot2::guide_legend(
        order = 2, ncol = 3, nrow = 1, byrow = TRUE,
        override.aes = list(size = 4, stroke = 1, fill = "grey80", colour = "black")
      )
    )
})

combined_all_regions <- (
  patchwork::wrap_plots(all_region_plots, ncol = 3, nrow = 2) +
    patchwork::plot_layout(guides = "collect")
) & ggplot2::theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.box = "vertical",
  plot.title = ggplot2::element_blank(),
  plot.subtitle = ggplot2::element_blank(),
  plot.tag = ggplot2::element_blank()
)

save_plot_formats(
  filename = file.path(PLOT_DIR, "pca_clusters_Calanoida_allRegions_allPairs_combined.png"),
  plot = combined_all_regions, width = 16, height = 10, dpi = 300
)

message("Done!")
message("  - Input data: ", IN_CSV)
message("  - Plots saved in: ", PLOT_DIR)
