#!/usr/bin/env Rscript
# =============================================================================
# Script:  09_behavior.R
# Author:  Rocio Rodriguez
# Purpose: Restore the supplied behavior calculations and plotting dimensions.
# Input:   outputs/Calanoida_data.csv; cluster icons from 00_setup.R.
# Outputs: Original behavior panels and summary; PNG/PDF/SVG/EPS; inputs for script 13.
# Source:  08.Behavior_mig&res_plots(1).R
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Converts day/night vertical abundance profiles into resident and migrator components for each Region x Cluster group.
# - Calculates weighted mean depth (WMD) and the original SE-based uncertainty intervals for the behavior summaries.
# - Summarizes resident/migrator proportions and migrator DVM amplitude from the day-night WMD difference.
# - Produces the resident/migrator proportion, WMD, and DVM panels and saves intermediate behavior results for the butterfly profiles in script 13.
# -----------------------------------------------------------------------------


source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["Behavior"]]

# -----------------------------------------------------------------------------
# Title:     Residents vs Migrators summary panels
# Author:    Rocio Rodriguez
# Date:      2026-10-20
#
# Description:
#   This script generates publication-ready panels from Calanoida net-depth profiles:
#     (1) Proportion of Residents vs Migrators by Cluster and Region
#     (2) All-data depth distribution violin with Migrator WMD points and Residents mean marker
#     (3) WMD CI-only panel with Migrators and Residents offset
#     (4) Migrator DVM amplitude with propagated uncertainty
#   Cluster colors are shown in individual plots and once in each combined figure,
#   alongside the separate Time and DVM status legends.
#   Text sizes and export dimensions are controlled in the USER INPUT section.
#
#   It also saves a simple CSV summary by Region x Cluster with:
#     - proportion Residents
#     - proportion Migrators
#     - Migrator Day/Night WMD
#     - Resident Day/Night WMD
#     - Resident mean WMD
#     - Migrator DVM amplitude
#
# Outputs saved to OUT_DIR:
#   - cluster_region_migrator_resident_WMD_DVM_summary.csv
#
# Outputs saved to PLOT_DIR:
#   - proportion_resident_vs_migrator.png
#   - WMD_violin_ALLDATA_plus_RESIDENT_meanPoint.png
#   - WMD_points_TRUE_CI_offset_migrators_residents.png
#   - DVM_migrators.png
#   - COMBINED_proportion_WMDviolin_DVMmigrators.png
#   - COMBINED_proportion_WMDci_DVMmigrators.png
# -----------------------------------------------------------------------------

# ---- USER INPUT -------------------------------------------------------------


if (!dir.exists(OUT_DIR))  dir.create(OUT_DIR, recursive = TRUE)
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

IN_CSV <- file.path(OUT_DIR, "Calanoida_data.csv")

CRUISE   <- "nh1208"
CATEGORY <- "Calanoida"
STATUS   <- "validated"

DEPTH_MAX    <- 1200
DEPTH_BREAKS <- seq(0, DEPTH_MAX, by = 100)

# Text sizes in points, shared by individual plots and combined figures.
# Increase these values here to enlarge all plot text in future exports.
plot_base_size    <- 35
axis_title_size   <- 35
axis_text_size    <- 35
strip_text_size   <- 35
legend_title_size <- 35
legend_text_size  <- 35
caption_text_size <- 35
panel_tag_size    <- 35

legend_key_size_cm <- 1.5

# Export dimensions in inches. The wider individual plots accommodate larger text.
individual_plot_width  <- 18
individual_plot_height <- 6
combined_plot_width    <- 28
combined_plot_height   <- 30
plot_dpi               <- 300

# Original net-depth codes (deepest boundary of each sampled layer).
# These values restore the supplied behavior script's WMD and DVM calculations.
net_depth_map <- c(
  n1 = 1000,
  n2 = 800,
  n3 = 600,
  n4 = 400,
  n5 = 200,
  n6 = 100,
  n7 = 50,
  n8 = 25
)

# Cluster order and colors
cluster_order <- c("ST", "SD", "MT", "MD", "LT", "LD")

cluster_colors <- c(
  "ST" = "#FED976",
  "MT" = "#FEB24C",
  "LT" = "#FD8D3C",
  "SD" = "#FC4E2A",
  "MD" = "#BD0026",
  "LD" = "#800026"
)

# Cluster icon PNGs for DVM plot
cluster_icons <- c(
  ST = file.path(ICON_DIR, "Asset 1.png"),
  MT = file.path(ICON_DIR, "Asset 2.png"),
  LT = file.path(ICON_DIR, "Asset 3.png"),
  SD = file.path(ICON_DIR, "Asset 4.png"),
  MD = file.path(ICON_DIR, "Asset 5.png"),
  LD = file.path(ICON_DIR, "Asset 6.png")
)

icon_size_min <- 0.09
icon_size_max <- 0.15

# ---- PACKAGES ---------------------------------------------------------------

required_pkgs <- c(
  "readr", "dplyr", "tidyr", "tibble",
  "ggplot2", "ggpubr", "ggimage",
  "ggpattern", "ragg", "grid"
)

new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))

if (length(new_pkgs) > 0) {
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

invisible(lapply(required_pkgs, function(p) library(p, character.only = TRUE)))

# ---- READ + FILTER ----------------------------------------------------------

Cal_data <- readr::read_csv(IN_CSV, show_col_types = FALSE) %>%
  dplyr::filter(
    Cruise == CRUISE,
    annotation_category %in% CATEGORY,
    annotation_status == STATUS,
    !(Region == "Subarctic" & Cluster == "ST")
  ) %>%
  dplyr::mutate(
    Region = factor(
      Region,
      levels = c("Subarctic", "Transition", "Subtropical")
    ),
    Cluster = factor(Cluster, levels = cluster_order),
    D_N = factor(D_N, levels = c("D", "N")),
    net = factor(net),
    net_depth = as.numeric(unname(net_depth_map[as.character(net)]))
  )

# ---- HELPERS ----------------------------------------------------------------

# Export the original plot at its original width, height, dpi, and background.
# Only the output device changes for the additional PDF, SVG, and EPS copies.
save_behavior_figure <- function(filename, plot, ...) {
  original_args <- list(...)
  stem <- tools::file_path_sans_ext(filename)
  for (extension in FIGURE_FORMATS) {
    args <- c(list(filename = paste0(stem, ".", extension), plot = plot), original_args)
    if (extension == "pdf") args$device <- grDevices::cairo_pdf
    if (extension == "svg") args$device <- svglite::svglite
    if (extension == "eps") args$device <- cairo_eps
    do.call(ggplot2::ggsave, args)
  }
}

# Use this theme for every panel, before extracting the shared legends.
behavior_plot_theme <- function() {
  ggplot2::theme_classic(base_size = plot_base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(size = plot_base_size),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      strip.text = ggplot2::element_text(size = strip_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.key.size = grid::unit(legend_key_size_cm, "cm"),
      legend.position = "right",
      plot.caption = ggplot2::element_text(size = caption_text_size, hjust = 0)
    )
}

scale_icon_size <- function(df, size_col = "size_mean") {
  if (size_col %in% names(df) && any(is.finite(df[[size_col]]))) {
    rng <- range(df[[size_col]], na.rm = TRUE)
    
    if (is.finite(rng[1]) && is.finite(rng[2]) && diff(rng) > 0) {
      df <- df %>%
        dplyr::mutate(
          icon_size = icon_size_min +
            (.data[[size_col]] - rng[1]) *
            (icon_size_max - icon_size_min) / diff(rng)
        )
    } else {
      df$icon_size <- (icon_size_min + icon_size_max) / 2
    }
  } else {
    df$icon_size <- (icon_size_min + icon_size_max) / 2
  }
  
  df
}

make_circle_polys <- function(df_pts, rx = 0.14, ry = 14, n = 60) {
  theta <- seq(0, 2 * pi, length.out = n)
  
  out <- vector("list", nrow(df_pts))
  
  for (i in seq_len(nrow(df_pts))) {
    x_center <- if ("x_pos" %in% names(df_pts)) {
      df_pts$x_pos[i]
    } else {
      df_pts$Cluster_num[i]
    }
    
    out[[i]] <- data.frame(
      Region = df_pts$Region[i],
      Cluster = df_pts$Cluster[i],
      Cluster_num = df_pts$Cluster_num[i],
      id = paste0("resmean_", i),
      x = x_center + rx * cos(theta),
      y = df_pts$WMD_res_mean[i] + ry * sin(theta)
    )
  }
  
  dplyr::bind_rows(out)
}

# -----------------------------------------------------------------------------
# 1) Depth profiles per Region x Cluster x net_depth
# -----------------------------------------------------------------------------

wide_depth <- Cal_data %>%
  dplyr::group_by(Region, Cluster, net_depth, D_N) %>%
  dplyr::summarise(
    ab = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = D_N,
    values_from = ab,
    values_fill = 0
  ) %>%
  dplyr::mutate(
    migr_D = pmax(D - N, 0),
    migr_N = pmax(N - D, 0),
    res = pmin(D, N)
  )

migr_long <- wide_depth %>%
  dplyr::select(Region, Cluster, net_depth, migr_D, migr_N) %>%
  tidyr::pivot_longer(
    cols = c(migr_D, migr_N),
    names_to = "D_N",
    values_to = "abundance"
  ) %>%
  dplyr::mutate(
    D_N = ifelse(D_N == "migr_D", "D", "N"),
    D_N = factor(D_N, levels = c("D", "N")),
    Component = "Migrators"
  ) %>%
  dplyr::filter(abundance > 0)

res_long <- wide_depth %>%
  dplyr::transmute(
    Region,
    Cluster,
    net_depth,
    res_D = res,
    res_N = res
  ) %>%
  tidyr::pivot_longer(
    cols = c(res_D, res_N),
    names_to = "D_N",
    values_to = "abundance"
  ) %>%
  dplyr::mutate(
    D_N = ifelse(D_N == "res_D", "D", "N"),
    D_N = factor(D_N, levels = c("D", "N")),
    Component = "Residents"
  ) %>%
  dplyr::filter(abundance > 0)

violin_df <- dplyr::bind_rows(res_long, migr_long) %>%
  dplyr::mutate(
    Component = factor(Component, levels = c("Residents", "Migrators"))
  )

# -----------------------------------------------------------------------------
# 2) WMD and original SE-based intervals
# -----------------------------------------------------------------------------

# Preserve the supplied approximation: k is the number of occupied bins;
# SE_w = SD_w / sqrt(k), and interval limits are WMD +/- 1.96 * SE_w.
# These are the original SE-derived intervals, rather than depth percentiles.
wmd_stats <- violin_df %>%
  dplyr::group_by(Region, Cluster, Component, D_N) %>%
  dplyr::summarise(
    k = sum(abundance > 0, na.rm = TRUE),
    total_abundance = sum(abundance, na.rm = TRUE),
    
    WMD = ifelse(
      total_abundance > 0,
      sum(abundance * net_depth, na.rm = TRUE) / total_abundance,
      NA_real_
    ),
    
    var_w = ifelse(
      total_abundance > 1,
      sum(abundance * (net_depth - WMD)^2, na.rm = TRUE) /
        (total_abundance - 1),
      NA_real_
    ),
    
    SD_w = sqrt(var_w),
    SE_w = SD_w / sqrt(pmax(1, k)),
    
    WMD_lower = WMD - 1.96 * SE_w,
    WMD_upper = WMD + 1.96 * SE_w,
    
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(WMD)) %>%
  dplyr::mutate(
    WMD_lower = pmax(WMD_lower, 0),
    WMD_upper = pmin(WMD_upper, DEPTH_MAX)
  )

# -----------------------------------------------------------------------------
# 3) Proportion Residents vs Migrators
# -----------------------------------------------------------------------------

wide_depth_norm <- wide_depth %>%
  dplyr::group_by(Region, Cluster) %>%
  dplyr::mutate(
    sumD = sum(D, na.rm = TRUE),
    sumN = sum(N, na.rm = TRUE),
    Dp = dplyr::if_else(sumD > 0, D / sumD, NA_real_),
    Np = dplyr::if_else(sumN > 0, N / sumN, NA_real_)
  ) %>%
  dplyr::ungroup()

mig_by_cluster_region <- wide_depth_norm %>%
  dplyr::group_by(Region, Cluster) %>%
  dplyr::summarise(
    mig_prop = 0.5 * sum(abs(Dp - Np), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    mig_prop = pmin(pmax(mig_prop, 0), 1),
    res_prop = dplyr::if_else(!is.na(mig_prop), 1 - mig_prop, NA_real_),
    Region = factor(Region, levels = c("Subarctic", "Transition", "Subtropical")),
    Cluster = factor(Cluster, levels = cluster_order)
  )

plot_df_prop <- mig_by_cluster_region %>%
  tidyr::pivot_longer(
    cols = c(res_prop, mig_prop),
    names_to = "DVM_status",
    values_to = "prop"
  ) %>%
  dplyr::mutate(
    DVM_status = dplyr::recode(
      DVM_status,
      mig_prop = "Migrators",
      res_prop = "Residents"
    ),
    DVM_status = factor(DVM_status, levels = c("Residents", "Migrators"))
  )

proportion_migvsres_cluster <- ggplot2::ggplot(
  plot_df_prop,
  ggplot2::aes(
    x = Cluster,
    y = prop,
    fill = Cluster,
    pattern = DVM_status
  )
) +
  ggpattern::geom_col_pattern(
    colour = "black",
    linewidth = 0.25,
    width = 0.85,
    position = "stack",
    pattern_angle = 45,
    pattern_fill = "white",
    pattern_colour = "white",
    pattern_density = 0.1,
    pattern_spacing = 0.04,
    pattern_alpha = 1
  ) +
  ggplot2::facet_wrap(~ Region, nrow = 1) +
  ggplot2::scale_fill_manual(
    name = "Cluster",
    values = cluster_colors,
    breaks = cluster_order,
    limits = cluster_order,
    drop = FALSE
  ) +
  ggpattern::scale_pattern_manual(
    name = "DVM status",
    values = c(
      "Residents" = "stripe",
      "Migrators" = "none"
    )
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      title = "Cluster",
      title.position = "top",
      order = 1,
      override.aes = list(pattern = "none", alpha = 1)
    ),
    pattern = ggplot2::guide_legend(
      title = "DVM status",
      title.position = "top",
      order = 2,
      override.aes = list(
        fill = "darkgrey",
        colour = "black",
        alpha = 1,
        pattern_fill = "white",
        pattern_colour = "white",
        pattern_density = 0.08,
        pattern_spacing = 0.04
      )
    )
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Relative abundance"
  ) +
  behavior_plot_theme()

# -----------------------------------------------------------------------------
# 4) Data prep for WMD panels
# -----------------------------------------------------------------------------

wmd_mig_pts <- wmd_stats %>%
  dplyr::filter(Component == "Migrators") %>%
  dplyr::mutate(
    Cluster_num = as.numeric(factor(Cluster, levels = cluster_order)),
    x_pos = Cluster_num - 0.15
  )

res_daynight <- wmd_stats %>%
  dplyr::filter(Component == "Residents") %>%
  dplyr::mutate(
    Cluster_num = as.numeric(factor(Cluster, levels = cluster_order)),
    x_pos = Cluster_num + 0.15
  )

res_mean <- wmd_stats %>%
  dplyr::filter(Component == "Residents") %>%
  dplyr::group_by(Region, Cluster) %>%
  dplyr::summarise(
    WMD_res_mean = mean(WMD, na.rm = TRUE),
    
    SE_res = sqrt(sum(SE_w^2, na.rm = TRUE)) / dplyr::n(),
    
    WMD_res_lower = WMD_res_mean - 1.96 * SE_res,
    WMD_res_upper = WMD_res_mean + 1.96 * SE_res,
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    WMD_res_lower = pmax(WMD_res_lower, 0),
    WMD_res_upper = pmin(WMD_res_upper, DEPTH_MAX),
    Cluster_num = as.numeric(factor(Cluster, levels = cluster_order)),
    x_pos = Cluster_num + 0.18
  ) %>%
  dplyr::filter(is.finite(WMD_res_mean))

res_mean_polys <- make_circle_polys(
  res_mean,
  rx = 0.12,
  ry = 24,
  n = 60
) %>%
  dplyr::mutate(
    poly_grp = interaction(id, drop = TRUE)
  )

pd2 <- ggplot2::position_dodge(width = 0.80)

# -----------------------------------------------------------------------------
# 5) WMD violin panel
# -----------------------------------------------------------------------------

p_all_violin_plus_resmean <- ggplot2::ggplot(
  Cal_data %>%
    dplyr::mutate(
      Cluster_num = as.numeric(factor(Cluster, levels = cluster_order))
    ),
  ggplot2::aes(
    x = Cluster_num,
    y = net_depth
  )
) +
  ggplot2::geom_violin(
    ggplot2::aes(
      weight = abundance,
      fill = Cluster,
      colour = Cluster,
      alpha = D_N,
      group = interaction(Cluster_num, D_N)
    ),
    trim = TRUE,
    scale = "width",
    position = pd2
  ) +
  
  ggplot2::geom_point(
    data = dplyr::filter(wmd_mig_pts, D_N == "D"),
    ggplot2::aes(
      x = Cluster_num,
      y = WMD,
      colour = Cluster,
      shape = "Day",
      group = interaction(Cluster_num, D_N)
    ),
    inherit.aes = FALSE,
    fill = "white",
    size = 6,
    stroke = 1.0,
    position = pd2,
    show.legend = c(shape = TRUE, colour = FALSE, fill = FALSE, alpha = FALSE)
  ) +
  
  ggplot2::geom_point(
    data = dplyr::filter(wmd_mig_pts, D_N == "N"),
    ggplot2::aes(
      x = Cluster_num,
      y = WMD,
      colour = Cluster,
      fill = Cluster,
      shape = "Night",
      group = interaction(Cluster_num, D_N)
    ),
    inherit.aes = FALSE,
    size = 6,
    stroke = 1.0,
    position = pd2,
    show.legend = c(shape = TRUE, colour = FALSE, fill = FALSE, alpha = FALSE)
  ) +
  
  ggpattern::geom_polygon_pattern(
    data = res_mean_polys,
    ggplot2::aes(
      x = x,
      y = y,
      group = poly_grp,
      fill = Cluster,
      colour = Cluster
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.03,
    pattern_spacing = 0.03,
    pattern_fill = "white",
    pattern_colour = "white",
    show.legend = FALSE
  ) +
  
  ggplot2::facet_wrap(~ Region, nrow = 1) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(cluster_order),
    labels = cluster_order
  ) +
  ggplot2::scale_fill_manual(
    values = cluster_colors,
    breaks = cluster_order,
    limits = cluster_order
  ) +
  ggplot2::scale_colour_manual(
    values = cluster_colors,
    breaks = cluster_order,
    limits = cluster_order
  ) +
  ggplot2::scale_alpha_manual(
    values = c("D" = 0.25, "N" = 0.75),
    guide = "none"
  ) +
  ggplot2::scale_shape_manual(
    name = "Time",
    values = c("Day" = 21, "Night" = 21)
  ) +
  ggplot2::scale_y_reverse(
    name = "Depth (m)",
    limits = c(DEPTH_MAX, 0)
  ) +
  ggplot2::labs(x = NULL) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      title = "Cluster",
      title.position = "top",
      order = 1,
      override.aes = list(alpha = 1, colour = "black")
    ),
    colour = "none",
    alpha = "none",
    shape = ggplot2::guide_legend(
      title = "Time",
      title.position = "top",
      order = 2,
      override.aes = list(
        colour = "black",
        fill = c("white", "black")
      )
    )
  ) +
  behavior_plot_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# -----------------------------------------------------------------------------
# 6) WMD CI-only panel with offset Migrators and Residents
# -----------------------------------------------------------------------------

p_WMD_only <- ggplot2::ggplot() +
  
  ggplot2::geom_errorbar(
    data = dplyr::filter(wmd_mig_pts, D_N == "D"),
    ggplot2::aes(
      x = x_pos,
      ymin = WMD,
      ymax = WMD_upper,
      colour = Cluster
    ),
    width = 0.15,
    linewidth = 1.5,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_errorbar(
    data = dplyr::filter(wmd_mig_pts, D_N == "N"),
    ggplot2::aes(
      x = x_pos,
      ymin = WMD_lower,
      ymax = WMD,
      colour = Cluster
    ),
    width = 0.15,
    linewidth = 1.5,
    show.legend = FALSE
  ) +
  
  ggplot2::geom_point(
    data = dplyr::filter(wmd_mig_pts, D_N == "D"),
    ggplot2::aes(
      x = x_pos,
      y = WMD,
      colour = Cluster,
      shape = "Day"
    ),
    fill = "white",
    size = 11,
    stroke = 1.0,
    show.legend = c(shape = TRUE, colour = FALSE, fill = FALSE)
  ) +
  
  ggplot2::geom_point(
    data = dplyr::filter(wmd_mig_pts, D_N == "N"),
    ggplot2::aes(
      x = x_pos,
      y = WMD,
      colour = Cluster,
      fill = Cluster,
      shape = "Night"
    ),
    size = 11,
    stroke = 1.0,
    # Enable fill so the cluster legend is drawn for this panel and its composite.
    show.legend = c(shape = TRUE, colour = FALSE, fill = TRUE)
  ) +
  
  ggplot2::geom_errorbar(
    data = res_mean,
    ggplot2::aes(
      x = x_pos,
      ymin = WMD_res_lower,
      ymax = WMD_res_upper,
      colour = Cluster
    ),
    width = 0.15,
    linewidth = 1.5,
    show.legend = FALSE
  ) +
  
  ggpattern::geom_polygon_pattern(
    data = res_mean_polys,
    ggplot2::aes(
      x = x,
      y = y,
      group = poly_grp,
      fill = Cluster,
      colour = Cluster
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    pattern = "stripe",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.012,
    pattern_fill = "white",
    pattern_colour = "white",
    show.legend = FALSE
  ) +
  
  ggplot2::facet_wrap(~ Region, nrow = 1) +
  ggplot2::scale_x_continuous(
    breaks = seq_along(cluster_order),
    labels = cluster_order
  ) +
  ggplot2::scale_fill_manual(
    values = cluster_colors,
    breaks = cluster_order,
    limits = cluster_order
  ) +
  ggplot2::scale_colour_manual(
    values = cluster_colors,
    breaks = cluster_order,
    limits = cluster_order
  ) +
  ggplot2::scale_shape_manual(
    name = "Time",
    values = c("Day" = 21, "Night" = 21)
  ) +
  ggplot2::scale_y_reverse(
    name = "Depth (m)",
    limits = c(DEPTH_MAX, 0)
  ) +
  ggplot2::labs(
    x = NULL,
    caption = "Migrators are offset left; Residents are offset right."
  ) +
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      title = "Cluster",
      title.position = "top",
      order = 1,
      override.aes = list(shape = 21, colour = "black", size = 4.5, alpha = 1)
    ),
    colour = "none",
    shape = ggplot2::guide_legend(
      title = "Time",
      title.position = "top",
      order = 2,
      override.aes = list(
        colour = "black",
        fill = c("white", "black")
      )
    )
  ) +
  behavior_plot_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# -----------------------------------------------------------------------------
# 7) Migrators DVM with propagated CI
# -----------------------------------------------------------------------------

WMD_migr_data <- migr_long %>%
  dplyr::group_by(Region, Cluster, D_N) %>%
  dplyr::summarise(
    k = sum(abundance > 0, na.rm = TRUE),
    total_abundance = sum(abundance, na.rm = TRUE),
    
    WMD = ifelse(
      total_abundance > 0,
      sum(abundance * net_depth, na.rm = TRUE) / total_abundance,
      NA_real_
    ),
    
    var_w = ifelse(
      total_abundance > 1,
      sum(abundance * (net_depth - WMD)^2, na.rm = TRUE) /
        (total_abundance - 1),
      NA_real_
    ),
    
    SD_w = sqrt(var_w),
    SE_w = SD_w / sqrt(pmax(1, k)),
    
    .groups = "drop"
  )

DVM_migr_data <- WMD_migr_data %>%
  dplyr::group_by(Region, Cluster) %>%
  dplyr::filter(all(c("D", "N") %in% D_N)) %>%
  dplyr::summarise(
    WMD_Day = WMD[D_N == "D"],
    WMD_Night = WMD[D_N == "N"],
    
    DVM = WMD_Day - WMD_Night,
    
    SE_D = SE_w[D_N == "D"],
    SE_N = SE_w[D_N == "N"],
    
    CI_Lower = DVM - 1.96 * sqrt(SE_D^2 + SE_N^2),
    CI_Upper = DVM + 1.96 * sqrt(SE_D^2 + SE_N^2),
    
    .groups = "drop"
  ) %>%
  dplyr::filter(is.finite(DVM))

# Scale icons using esd if present
if ("esd" %in% names(Cal_data)) {
  size_by_cluster <- Cal_data %>%
    dplyr::group_by(Region, Cluster) %>%
    dplyr::summarise(
      size_mean = mean(esd, na.rm = TRUE),
      .groups = "drop"
    )
  
  DVM_migr_data <- DVM_migr_data %>%
    dplyr::left_join(size_by_cluster, by = c("Region", "Cluster"))
} else {
  DVM_migr_data$size_mean <- NA_real_
}

DVM_migr_data <- scale_icon_size(DVM_migr_data, "size_mean")

icons_exist <- file.exists(cluster_icons)
cluster_icons_ok <- cluster_icons[icons_exist]

DVM_migr_icons <- DVM_migr_data %>%
  dplyr::mutate(
    icon = unname(cluster_icons_ok[as.character(Cluster)])
  ) %>%
  dplyr::filter(!is.na(icon))

dvm_plot_mig <- ggplot2::ggplot(
  DVM_migr_data,
  ggplot2::aes(
    x = Cluster,
    y = DVM,
    color = Cluster
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = CI_Lower,
      ymax = CI_Upper
    ),
    width = 0.25
  ) +
  ggimage::geom_image(
    data = DVM_migr_icons,
    ggplot2::aes(
      x = Cluster,
      y = DVM,
      image = icon,
      size = icon_size
    ),
    inherit.aes = FALSE,
    asp = 1,
    alpha = 0.9
  ) +
  ggplot2::scale_size_identity() +
  ggplot2::scale_color_manual(
    values = cluster_colors,
    breaks = cluster_order,
    limits = cluster_order
  ) +
  ggplot2::labs(
    x = "Cluster",
    y = "DVM (m) - Migrators",
    color = "Cluster"
  ) +
  ggplot2::facet_wrap(~ Region, nrow = 1) +
  behavior_plot_theme() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

# -----------------------------------------------------------------------------
# 8) EXPORT SIMPLE SUMMARY CSV: proportions, WMD, and DVM by Region x Cluster
# -----------------------------------------------------------------------------

prop_summary <- mig_by_cluster_region %>%
  dplyr::select(
    Region,
    Cluster,
    prop_residents = res_prop,
    prop_migrators = mig_prop
  )

wmd_simple <- wmd_stats %>%
  dplyr::select(
    Region,
    Cluster,
    Component,
    D_N,
    WMD
  ) %>%
  tidyr::pivot_wider(
    names_from = c(Component, D_N),
    values_from = WMD,
    names_glue = "{Component}_{D_N}_WMD"
  )

res_wmd_simple <- res_mean %>%
  dplyr::select(
    Region,
    Cluster,
    Residents_mean_WMD = WMD_res_mean
  )

dvm_simple <- DVM_migr_data %>%
  dplyr::select(
    Region,
    Cluster,
    Migrators_DVM = DVM
  )

cluster_region_summary <- prop_summary %>%
  dplyr::left_join(wmd_simple, by = c("Region", "Cluster")) %>%
  dplyr::left_join(res_wmd_simple, by = c("Region", "Cluster")) %>%
  dplyr::left_join(dvm_simple, by = c("Region", "Cluster")) %>%
  dplyr::arrange(Region, Cluster)

readr::write_csv(
  cluster_region_summary,
  file.path(OUT_DIR, "cluster_region_migrator_resident_WMD_DVM_summary.csv")
)

message(
  "Saved simple cluster-region summary CSV: ",
  file.path(OUT_DIR, "cluster_region_migrator_resident_WMD_DVM_summary.csv")
)

# -----------------------------------------------------------------------------
# 9) Combined plot 1: Proportion -> WMD violin -> Migrator DVM
# -----------------------------------------------------------------------------

p_violin_leg_cluster <- p_all_violin_plus_resmean +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::guides(
    shape = "none",
    fill = ggplot2::guide_legend(
      title = "Cluster",
      title.position = "top",
      order = 1,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(alpha = 1, colour = "black")
    )
  )

p_violin_leg_time <- p_all_violin_plus_resmean +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::guides(fill = "none")

p_prop_leg <- proportion_migvsres_cluster +
  ggplot2::theme(legend.position = "bottom") +
  # The cluster key is extracted separately; keep only status in this block.
  ggplot2::guides(fill = "none")

legend_cluster_violin <- ggpubr::get_legend(p_violin_leg_cluster)
legend_time_violin <- ggpubr::get_legend(p_violin_leg_time)
legend_status_violin <- ggpubr::get_legend(p_prop_leg)

legend_row_violin <- ggpubr::ggarrange(
  ggpubr::as_ggplot(legend_cluster_violin),
  ggpubr::as_ggplot(legend_time_violin),
  ggpubr::as_ggplot(legend_status_violin),
  ncol = 3,
  widths = c(1.5, 1, 1)
)

p_prop_comb1 <- proportion_migvsres_cluster +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(x = NULL)

p_violin_comb <- p_all_violin_plus_resmean +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank()
  ) +
  ggplot2::labs(x = NULL)

p_dvm_comb1 <- dvm_plot_mig +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank()
  )

combined_main_violin <- ggpubr::ggarrange(
  p_prop_comb1,
  p_violin_comb,
  p_dvm_comb1,
  ncol = 1,
  nrow = 3,
  labels = c("A", "B", "C"),
  font.label = list(size = panel_tag_size, face = "plain"),
  heights = c(1.0, 1.2, 1.0)
)

combined_prop_WMDviolin_dvmmig <- ggpubr::ggarrange(
  combined_main_violin,
  legend_row_violin,
  ncol = 1,
  heights = c(1, 0.18)
)

# -----------------------------------------------------------------------------
# 10) Combined plot 2: Proportion -> WMD CI-only -> Migrator DVM
# -----------------------------------------------------------------------------

p_wmd_leg_cluster <- p_WMD_only +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::guides(
    shape = "none",
    fill = ggplot2::guide_legend(
      title = "Cluster",
      title.position = "top",
      order = 1,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(shape = 21, colour = "black", size = 7, alpha = 1)
    )
  )

p_wmd_leg_time <- p_WMD_only +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::guides(fill = "none")

legend_cluster_ci <- ggpubr::get_legend(p_wmd_leg_cluster)
legend_time_ci <- ggpubr::get_legend(p_wmd_leg_time)
legend_status_ci <- ggpubr::get_legend(p_prop_leg)

legend_row_ci <- ggpubr::ggarrange(
  ggpubr::as_ggplot(legend_cluster_ci),
  ggpubr::as_ggplot(legend_time_ci),
  ggpubr::as_ggplot(legend_status_ci),
  ncol = 3,
  widths = c(1.5, 1, 1)
)

p_prop_comb2 <- proportion_migvsres_cluster +
  ggplot2::theme(legend.position = "none") +
  ggplot2::labs(x = NULL)

p_wmd_comb <- p_WMD_only +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    plot.caption = ggplot2::element_blank()
  ) +
  ggplot2::labs(x = NULL)

p_dvm_comb2 <- dvm_plot_mig +
  ggplot2::theme(
    legend.position = "none",
    strip.text = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank()
  )

combined_main_ci <- ggpubr::ggarrange(
  p_prop_comb2,
  p_wmd_comb,
  p_dvm_comb2,
  ncol = 1,
  nrow = 3,
  labels = c("A", "B", "C"),
  font.label = list(size = panel_tag_size, face = "plain"),
  heights = c(1.0, 1.2, 1.0)
)

combined_prop_WMDci_dvmmig <- ggpubr::ggarrange(
  combined_main_ci,
  legend_row_ci,
  ncol = 1,
  heights = c(1, 0.18)
)

# -----------------------------------------------------------------------------
# 11) SAVE PLOTS
# -----------------------------------------------------------------------------

save_behavior_figure(
  filename = file.path(PLOT_DIR, "proportion_resident_vs_migrator.png"),
  plot = proportion_migvsres_cluster,
  device = ragg::agg_png,
  width = individual_plot_width,
  height = individual_plot_height,
  dpi = plot_dpi,
  bg = "transparent"
)

save_behavior_figure(
  filename = file.path(PLOT_DIR, "WMD_violin_ALLDATA_plus_RESIDENT_meanPoint.png"),
  plot = p_all_violin_plus_resmean,
  width = individual_plot_width,
  height = individual_plot_height,
  dpi = plot_dpi
)

save_behavior_figure(
  filename = file.path(PLOT_DIR, "WMD_points_TRUE_CI_offset_migrators_residents.png"),
  plot = p_WMD_only,
  width = individual_plot_width,
  height = individual_plot_height,
  dpi = plot_dpi
)

save_behavior_figure(
  filename = file.path(PLOT_DIR, "DVM_migrators.png"),
  plot = dvm_plot_mig,
  width = individual_plot_width,
  height = individual_plot_height,
  dpi = plot_dpi
)

save_behavior_figure(
  filename = file.path(PLOT_DIR, "COMBINED_proportion_WMDviolin_DVMmigrators.png"),
  plot = combined_prop_WMDviolin_dvmmig,
  width = combined_plot_width,
  height = combined_plot_height,
  dpi = plot_dpi
)

save_behavior_figure(
  filename = file.path(PLOT_DIR, "COMBINED_proportion_WMDci_DVMmigrators.png"),
  plot = combined_prop_WMDci_dvmmig,
  width = combined_plot_width,
  height = combined_plot_height,
  dpi = plot_dpi
)

# -----------------------------------------------------------------------------
# 12) EXPORT OBJECTS
# -----------------------------------------------------------------------------

assign("Cal_data", Cal_data, envir = .GlobalEnv)
assign("wide_depth", wide_depth, envir = .GlobalEnv)
assign("wide_depth_norm", wide_depth_norm, envir = .GlobalEnv)

assign("violin_df", violin_df, envir = .GlobalEnv)
assign("wmd_stats", wmd_stats, envir = .GlobalEnv)
assign("wmd_mig_pts", wmd_mig_pts, envir = .GlobalEnv)
assign("res_daynight", res_daynight, envir = .GlobalEnv)
assign("res_mean", res_mean, envir = .GlobalEnv)

assign("WMD_migr_data", WMD_migr_data, envir = .GlobalEnv)
assign("DVM_migr_data", DVM_migr_data, envir = .GlobalEnv)

assign("cluster_region_summary", cluster_region_summary, envir = .GlobalEnv)

assign("proportion_migvsres_cluster", proportion_migvsres_cluster, envir = .GlobalEnv)
assign("p_all_violin_plus_resmean", p_all_violin_plus_resmean, envir = .GlobalEnv)
assign("p_WMD_only", p_WMD_only, envir = .GlobalEnv)
assign("dvm_plot_mig", dvm_plot_mig, envir = .GlobalEnv)
assign("combined_prop_WMDviolin_dvmmig", combined_prop_WMDviolin_dvmmig, envir = .GlobalEnv)
assign("combined_prop_WMDci_dvmmig", combined_prop_WMDci_dvmmig, envir = .GlobalEnv)

message("Done! Plots saved in: ", PLOT_DIR)
message("Simple summary CSV saved in: ", OUT_DIR)
message("WMD CIs are still calculated for plotting, but they are not exported in the CSV.")
message("Generated both WMD violin and WMD CI-only offset panels.")

# 13) Save profile inputs for the existing butterfly script --------------------
# Script 13 maps these bins back to its original depth labels. Attach that coordinate to the same
# net-bin abundances; this does not change any behavior statistic or plot above.
profiles_for_butterfly <- wide_depth %>%
  dplyr::left_join(
    dplyr::select(NET_BINS, depth_max, depth_mid),
    by = c("net_depth" = "depth_max")
  ) %>%
  dplyr::transmute(
    Region, Cluster, depth_mid, D, N,
    Residents = res, Migrators_D = migr_D, Migrators_N = migr_N
  )
readr::write_csv(profiles_for_butterfly, file.path(OUT_DIR, "behavior_depth_profiles.csv"))
readr::write_csv(wmd_stats, file.path(OUT_DIR, "behavior_SE_depth_intervals.csv"))
saveRDS(
  list(
    profiles = profiles_for_butterfly,
    proportions = prop_summary,
    migrator_depth = dplyr::filter(wmd_stats, Component == "Migrators"),
    resident_depth = res_mean,
    dvm = DVM_migr_data,
    depth_convention = "Original net depths: 25, 50, 100, 200, 400, 600, 800, 1000 m",
    interval_method = "Original weighted-variance/occupied-bin SE approximation"
  ),
  file.path(OUT_DIR, "behavior_results.rds")
)
