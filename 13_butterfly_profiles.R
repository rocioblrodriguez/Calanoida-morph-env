#!/usr/bin/env Rscript
# =============================================================================
# Script:  13_butterfly_profiles.R
# Author:  Rocio Rodriguez
# Purpose: Draw original butterfly profiles from script 09 abundance profiles.
# Style:   Original supplied plot construction, dimensions, and combinations.
#          Shared paths/colors, PNG/PDF/SVG export, and plain panel tags.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Reads the resident/migrator depth profiles saved by script 09 rather than recalculating behavior classifications.
# - Converts each Region x Cluster profile to relative abundance across depth for the mirrored day/night butterfly display.
# - Day and night are plotted on opposite sides of zero, while pattern distinguishes residents from migrators.
# - Produces one profile figure per region and a combined figure with a single shared legend.
# -----------------------------------------------------------------------------


source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["Butterfly"]]
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpattern)
library(patchwork)
library(cowplot)

profiles <- readRDS(file.path(OUT_DIR, "behavior_results.rds"))$profiles
segments <- profiles |>
  left_join(select(NET_BINS, depth_mid, depth_m = depth_max), by = "depth_mid") |>
  transmute(Region, Cluster, depth_mid, depth_m, Residents_D = Residents, Residents_N = Residents,
            Migrators_D, Migrators_N) |>
  pivot_longer(c(Residents_D, Residents_N, Migrators_D, Migrators_N),
               names_to = c("mig_status", "D_N"), names_sep = "_", values_to = "value") |>
  filter(value > 0) |>
  group_by(Region, Cluster, mig_status, D_N) |>
  mutate(value_plot = value / sum(value)) |>
  ungroup() |>
  mutate(value_signed = if_else(D_N == "D", value_plot, -value_plot),
         depth_layer = factor(depth_m, levels = rev(NET_BINS$depth_max)),
         Region = factor(Region, levels = REGION_LEVELS),
         Cluster = factor(Cluster, levels = NODE_CLUSTER_LEVELS),
         D_N = factor(D_N, levels = c("N", "D")),
         mig_status = factor(mig_status, levels = c("Residents", "Migrators")))
cluster_cols <- CLUSTER_COLORS
region_tags <- REGION_TAGS
regions <- REGION_LEVELS

make_butterfly_plot <- function(
    plot_df,
    xlab_text,
    x_limits = NULL,
    tag_text = NULL
) {
  
  p_base <- ggplot(
    plot_df,
    aes(
      x       = value_signed,
      y       = depth_layer,
      fill    = Cluster,
      pattern = mig_status,
      alpha   = D_N
    )
  ) +
    ggpattern::geom_col_pattern(
      width           = 0.9,
      position        = position_stack(reverse = TRUE),
      colour          = "black",
      linewidth       = 0.25,
      pattern_angle   = 45,
      pattern_fill    = "white",
      pattern_colour  = "white",
      pattern_density = 0.1,
      pattern_spacing = 0.06,
      pattern_alpha   = 1
    ) +
    geom_vline(
      xintercept = 0,
      linewidth  = 0.8,
      linetype   = 2
    ) +
    scale_fill_manual(
      values = cluster_cols,
      drop   = FALSE,
      name   = "Cluster"
    ) +
    scale_alpha_manual(
      values = c("N" = 0.9, "D" = 0.6)
    ) +
    scale_pattern_manual(
      name   = "DVM status",
      values = c("Residents" = "stripe", "Migrators" = "none")
    ) +
    scale_x_continuous(
      limits = x_limits,
      labels = abs,
      name   = xlab_text
    ) +
    guides(
      fill = guide_legend(
        order = 1,
        ncol  = 1
      ),
      pattern = guide_legend(
        order = 2,
        ncol  = 1,
        override.aes = list(
          fill            = "darkgrey",
          colour          = "black",
          alpha           = 1,
          pattern_fill    = "white",
          pattern_colour  = "white",
          pattern_density = 0.06,
          pattern_spacing = 0.01
        )
      ),
      alpha = "none"
    ) +
    labs(
      y   = "Depth (m)",
      x   = xlab_text,
      tag = tag_text
    ) +
    facet_grid(
      Cluster ~ .,
      switch = "y",
      drop = FALSE
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.title      = element_text(size = 18, face = "bold"),
      legend.text       = element_text(size = 16),
      legend.key.size   = grid::unit(0.8, "cm"),
      strip.text.y.left = element_blank(),
      strip.background  = element_blank(),
      plot.tag          = element_text(size = 20, face = "plain"),
      plot.tag.position = c(0.02, 0.98),
      plot.margin       = margin(t = 25, r = 10, b = 10, l = 10)
    )
  
  p_with_legend <- p_base +
    theme(
      legend.position = "right",
      legend.box      = "vertical"
    )
  
  p_no_legend <- p_base +
    theme(
      legend.position = "none"
    )
  
  plot_with_legend <- cowplot::ggdraw(p_with_legend) +
    cowplot::draw_label(
      "Night",
      x        = 0.28,
      y        = 0.985,
      fontface = "bold",
      size     = 14
    ) +
    cowplot::draw_label(
      "Day",
      x        = 0.72,
      y        = 0.985,
      fontface = "bold",
      size     = 14
    )
  
  plot_no_legend <- cowplot::ggdraw(p_no_legend) +
    cowplot::draw_label(
      "Night",
      x        = 0.28,
      y        = 0.985,
      fontface = "bold",
      size     = 14
    ) +
    cowplot::draw_label(
      "Day",
      x        = 0.72,
      y        = 0.985,
      fontface = "bold",
      size     = 14
    )
  
  return(
    list(
      plot_with_legend = plot_with_legend,
      plot_no_legend   = plot_no_legend,
      legend_source    = p_with_legend
    )
  )
}


plots_relative_with_legend <- setNames(
  vector("list", length(regions)),
  as.character(regions)
)

plots_relative_no_legend <- setNames(
  vector("list", length(regions)),
  as.character(regions)
)

legend_source_relative <- NULL

for (reg in regions) {
  
  reg_char <- as.character(reg)
  reg_tag  <- region_tags[reg_char]
  
  segs_relative <- filter(segments, Region == reg)

  p_relative <- make_butterfly_plot(
    plot_df   = segs_relative,
    xlab_text = "Relative abundance",
    x_limits  = c(-2, 2),
    tag_text  = reg_tag
  )
  
  plots_relative_with_legend[[reg_char]] <- p_relative$plot_with_legend
  plots_relative_no_legend[[reg_char]]   <- p_relative$plot_no_legend
  
  if (is.null(legend_source_relative)) {
    legend_source_relative <- p_relative$legend_source
  }
  
  save_plot_formats(
    filename = file.path(PLOT_DIR, paste0("butterfly_relative_", reg_char, ".png")),
    plot     = p_relative$plot_with_legend,
    device   = ragg::agg_png,
    width    = 8,
    height   = 12,
    dpi      = 300,
    bg       = "transparent"
  )
}


# -----------------------------------------------------------------------------
# 10. COMBINED FIGURE WITH ONE SHARED LEGEND
# -----------------------------------------------------------------------------

shared_legend <- cowplot::get_legend(
  legend_source_relative +
    theme(
      legend.position = "right",
      legend.box      = "vertical"
    )
)

combined_panels <- patchwork::wrap_plots(
  plots_relative_no_legend,
  nrow = 1,
  ncol = 3
)

combined_relative <- cowplot::plot_grid(
  combined_panels,
  shared_legend,
  nrow       = 1,
  rel_widths = c(1, 0.13)
)

save_plot_formats(
  filename = file.path(PLOT_DIR, "butterfly_relative_ALLregions.png"),
  plot     = combined_relative,
  device   = ragg::agg_png,
  width    = 22,
  height   = 12,
  dpi      = 300,
  bg       = "transparent"
)



write_csv(segments, file.path(OUT_DIR, "butterfly_relative_profile_values.csv"))
