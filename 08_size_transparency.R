#!/usr/bin/env Rscript
# =============================================================================
# Script:  08_size_transparency.R
# Author:  Rocio Rodriguez
# Purpose: Draw original size/transparency panels with abundance marginals.
# Style:   Original supplied plot construction, dimensions, and combinations.
#          Shared paths/colors, PNG/PDF/SVG export, and plain panel tags.
# =============================================================================

source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["Size_transparency"]]
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggimage)

IN_CSV <- file.path(OUT_DIR, "Calanoida_data.csv")
Cop_data <- read_csv(IN_CSV, show_col_types = FALSE) |>
  dplyr::filter(!(Region == "Subarctic" & Cluster == "ST"))
X_LIMS_ESD <- c(0, 5)
Y_LIMS_TRANS <- c(100, 230)
N_BINS <- 30
icon_size_min <- 0.06
icon_size_max <- 0.15
cluster_order <- CLUSTER_LEVELS
cluster_cols <- CLUSTER_COLORS
cluster_icons <- CLUSTER_ICONS
region_letter <- REGION_TAGS

# ---- preprocess ----
Cop_data <- Cop_data %>%
  dplyr::mutate(
    major_mm_ind = as.numeric(major * process_particle_pixel_size_mm),
    esd_mm_ind   = as.numeric(esd * process_particle_pixel_size_mm),
    Transparency = .data[["mean"]],
    Region       = factor(
      Region,
      levels = c("Subarctic", "Transition", "Subtropical")
    ),
    Cluster      = trimws(as.character(Cluster)),
    Cluster      = factor(
      Cluster,
      levels = cluster_order
    )
  )

# Keep the same geographic order even if input rows are ordered differently.
regions <- levels(Cop_data$Region)
regions <- regions[regions %in% as.character(Cop_data$Region)]

# ---- lists to store plots and summaries ----
transp_esd_marginal_by_region <- list()
transp_esd_main_by_region     <- list()
summary_by_region_cluster     <- list()

# ---- build plots per Region ----
for (reg in regions) {
  
  reg_char <- as.character(reg)
  
  region_data <- Cop_data %>%
    dplyr::filter(Region == reg)
  
  
  # ---------------------------------------------------------------------------
  # Summary:
  # - mean ESD and mean Transparency
  # - 5–95% quantiles used only for plot error bars
  # - raw trait means used for CSV output
  # ---------------------------------------------------------------------------
  
  df_summary <- region_data %>%
    dplyr::group_by(Cluster, Region) %>%
    dplyr::summarise(
      mean_Transparency = mean(Transparency, na.rm = TRUE),
      
      mean_esd_mm       = mean(esd_mm_ind, na.rm = TRUE),
      mean_esd_px       = mean(.data[["esd"]], na.rm = TRUE),
      
      q05_Transparency  = stats::quantile(Transparency, 0.05, na.rm = TRUE),
      q95_Transparency  = stats::quantile(Transparency, 0.95, na.rm = TRUE),
      q05_esd_mm        = stats::quantile(esd_mm_ind, 0.05, na.rm = TRUE),
      q95_esd_mm        = stats::quantile(esd_mm_ind, 0.95, na.rm = TRUE),
      
      mean_major_mm     = mean(major_mm_ind, na.rm = TRUE),
      mean_major_px     = mean(.data[["major"]], na.rm = TRUE),
      
      mean_esd          = mean(.data[["esd"]], na.rm = TRUE),
      mean_mean         = mean(.data[["mean"]], na.rm = TRUE),
      mean_mode         = mean(.data[["mode"]], na.rm = TRUE),
      mean_circ         = mean(.data[["circ."]], na.rm = TRUE),
      mean_symetrieh    = mean(.data[["symetrieh"]], na.rm = TRUE),
      mean_perimareaexc = mean(.data[["perimareaexc"]], na.rm = TRUE),
      mean_slope        = mean(.data[["slope"]], na.rm = TRUE),
      
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      icon = cluster_icons[as.character(Cluster)]
    )
  
  # ---------------------------------------------------------------------------
  # Store summary for CSV.
  # No CI / quantile columns are included in this table.
  # ---------------------------------------------------------------------------
  
  summary_by_region_cluster[[reg_char]] <- df_summary %>%
    dplyr::select(
      Region,
      Cluster,
      
      mean_esd_mm,
      mean_esd_px,
      mean_Transparency,
      
      mean_major_mm,
      mean_major_px,
      mean_esd,
      mean_mean,
      mean_mode,
      mean_circ,
      mean_symetrieh,
      mean_perimareaexc,
      mean_slope
    )
  
  # ---- scale mean_esd_mm to icon size ----
  df_summary$icon_size <- icon_sizes(df_summary$mean_esd_mm, icon_size_min, icon_size_max)

  # ---------------------------------------------------------------------------
  # MAIN plot: error bars + icons at means
  # ---------------------------------------------------------------------------
  
  p_center <- ggplot(
    df_summary,
    aes(x = mean_esd_mm, y = mean_Transparency)
  ) +
    geom_errorbar(
      aes(
        ymin = q05_Transparency,
        ymax = q95_Transparency,
        color = Cluster
      ),
      width = 0.08,
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    geom_errorbarh(
      aes(
        xmin = q05_esd_mm,
        xmax = q95_esd_mm,
        color = Cluster
      ),
      height = 2.2,
      linewidth = 1.2,
      show.legend = TRUE
    ) +
    ggimage::geom_image(
      aes(image = icon, size = icon_size),
      asp = 1,
      inherit.aes = TRUE,
      show.legend = FALSE
    ) +
    scale_size_identity() +
    scale_color_manual(
      values = cluster_cols,
      limits = cluster_order,
      breaks = cluster_order,
      drop = FALSE,
      labels = c(
        "ST = Small Transparent",
        "SD = Small Dark",
        "MT = Medium Transparent",
        "MD = Medium Dark",
        "LT = Large Transparent",
        "LD = Large Dark"
      )
    ) +
    scale_x_continuous(
      limits = X_LIMS_ESD,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = Y_LIMS_TRANS,
      expand = c(0, 0)
    ) +
    labs(
      x = "ESD (mm)",
      y = "Transparency",
      color = "Cluster"
    ) +
    theme_classic(base_size = 24) +
    guides(
      color = guide_legend(
        ncol = 3,
        override.aes = list(linewidth = 2)
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.box.just = "left",
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5)
    )
  
  # Use a labeled copy for main-only figures. The marginal layout keeps the
  # original center plot, with its region title above the top distribution.
  p_main <- p_center +
    labs(title = reg_char, tag = unname(region_letter[reg_char])) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.tag = element_text(face = "plain", size = 26)
    )
  
  transp_esd_main_by_region[[reg_char]] <- p_main
  
  # ---------------------------------------------------------------------------
  # Top marginal: abundance-weighted histogram of ESD
  # ---------------------------------------------------------------------------
  
  p_top <- ggplot(
    region_data,
    aes(x = esd_mm_ind, weight = abundance, fill = Cluster)
  ) +
    geom_area(
      stat = "bin",
      bins = N_BINS,
      position = "stack",
      show.legend = FALSE
    ) +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_continuous(
      limits = X_LIMS_ESD,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 125000),
      expand = c(0, 0)
    ) +
    labs(
      title = reg_char,
      x = NULL,
      y = expression("Abundance " ~ m^{-2})
    ) +
    theme_classic(base_size = 24) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      axis.title.x    = element_blank(),
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.line.x     = element_blank(),
      legend.position = "none",
      plot.margin     = margin(b = 2, l = 8, r = 5, t = 25)
    )
  
  # ---------------------------------------------------------------------------
  # Right marginal: abundance-weighted histogram of Transparency
  # ---------------------------------------------------------------------------
  
  p_right <- ggplot(
    region_data,
    aes(x = Transparency, weight = abundance, fill = Cluster)
  ) +
    geom_area(
      stat = "bin",
      bins = N_BINS,
      position = "stack",
      show.legend = FALSE
    ) +
    coord_flip() +
    scale_fill_manual(values = cluster_cols, guide = "none") +
    scale_x_continuous(
      limits = Y_LIMS_TRANS,
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, 25000),
      breaks = c(5000, 15000),
      labels = c("5000", "15000"),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = expression("Abundance " ~ m^{-2})
    ) +
    theme_classic(base_size = 24) +
    theme(
      axis.title.y    = element_blank(),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.line.y     = element_blank(),
      legend.position = "none",
      plot.margin     = margin(t = 5, l = 2, r = 25, b = 5)
    )
  
  # ---------------------------------------------------------------------------
  # Combine marginal plot
  # ---------------------------------------------------------------------------
  
  design <- c(
    patchwork::area(t = 1, l = 1, b = 1, r = 3), # top
    patchwork::area(t = 1, l = 4, b = 1, r = 4), # spacer/tag
    patchwork::area(t = 2, l = 1, b = 4, r = 3), # center
    patchwork::area(t = 2, l = 4, b = 4, r = 4)  # right
  )
  
  spacer_tag <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::xlim(0, 1) +
    ggplot2::ylim(0, 1) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1,
      label = region_letter[reg_char],
      fontface = "plain",
      hjust = 1,
      vjust = 1,
      size = 10
    ) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(t = 25, r = 5, b = 0, l = 0)
    )
  
  combined <- p_top + spacer_tag + p_center + p_right +
    patchwork::plot_layout(design = design)
  
  transp_esd_marginal_by_region[[reg_char]] <- combined
  
  # ---------------------------------------------------------------------------
  # Save per-region plots
  # ---------------------------------------------------------------------------
  
  save_plot_formats(
    filename = file.path(
      PLOT_DIR,
      paste0("transp_esd_summary_", reg_char, ".png")
    ),
    plot   = p_main,
    width  = 10,
    height = 10,
    dpi    = 300,
    bg     = "transparent"
  )
  
  save_plot_formats(
    filename = file.path(
      PLOT_DIR,
      paste0("transp_esd_summary_marginal_", reg_char, ".png")
    ),
    plot   = combined,
    width  = 12,
    height = 12,
    dpi    = 300,
    bg     = "transparent"
  )
}

# -----------------------------------------------------------------------------
# Save CSV: mean raw traits per Region x Cluster
# -----------------------------------------------------------------------------

if (length(summary_by_region_cluster) > 0) {
  
  summary_out <- dplyr::bind_rows(summary_by_region_cluster) %>%
    dplyr::arrange(Region, Cluster)
  
  summary_csv_path <- file.path(
    OUT_DIR,
    "transp_esd_means_by_region_cluster.csv"
  )
  
  readr::write_csv(summary_out, summary_csv_path)
  
  message("Saved summary CSV: ", summary_csv_path)

  # -------------------------------------------------------------------------
  # Publication-style comparison table
  # -------------------------------------------------------------------------
  # Each gray row reports the overall cluster mean. The three white rows below
  # it report regional means and their percentage difference from the value in
  # each comparison region. Missing combinations remain NA (for example,
  # Subarctic ST after the requested filter).

  cluster_overall <- Cop_data %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(
      mean_esd_mm = mean(esd_mm_ind, na.rm = TRUE),
      mean_Transparency = mean(Transparency, na.rm = TRUE),
      .groups = "drop"
    )

  table_cluster_order <- rev(cluster_order)
  table_region_order <- c("Subarctic", "Transition", "Subtropical")

  fmt_number <- function(x) {
    ifelse(is.finite(x), sprintf("%.2f", x), "NA")
  }

  fmt_difference <- function(value, reference, same_region = FALSE) {
    if (same_region || !is.finite(value) || !is.finite(reference) || reference == 0) {
      return("")
    }
    sprintf("%d%%", round(100 * (value / reference - 1)))
  }

  comparison_rows <- list()
  row_index <- 1L

  for (cluster_name in table_cluster_order) {
    overall_row <- cluster_overall %>%
      dplyr::filter(as.character(Cluster) == cluster_name)

    comparison_rows[[row_index]] <- data.frame(
      row_type = "cluster",
      label = cluster_name,
      esd = if (nrow(overall_row)) fmt_number(overall_row$mean_esd_mm[1]) else "NA",
      esd_subarctic = "", esd_transition = "", esd_subtropical = "",
      transparency = if (nrow(overall_row)) fmt_number(overall_row$mean_Transparency[1]) else "NA",
      transparency_subarctic = "", transparency_transition = "",
      transparency_subtropical = "",
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L

    cluster_regional <- summary_out %>%
      dplyr::filter(as.character(Cluster) == cluster_name)

    esd_reference <- stats::setNames(
      cluster_regional$mean_esd_mm,
      as.character(cluster_regional$Region)
    )
    transparency_reference <- stats::setNames(
      cluster_regional$mean_Transparency,
      as.character(cluster_regional$Region)
    )

    for (region_name in table_region_order) {
      regional_row <- cluster_regional %>%
        dplyr::filter(as.character(Region) == region_name)

      esd_value <- if (nrow(regional_row)) regional_row$mean_esd_mm[1] else NA_real_
      transparency_value <- if (nrow(regional_row)) regional_row$mean_Transparency[1] else NA_real_

      get_reference <- function(values, name) {
        if (name %in% names(values)) unname(values[[name]]) else NA_real_
      }

      comparison_rows[[row_index]] <- data.frame(
        row_type = "region",
        label = region_name,
        esd = fmt_number(esd_value),
        esd_subarctic = fmt_difference(
          esd_value, get_reference(esd_reference, "Subarctic"), region_name == "Subarctic"
        ),
        esd_transition = fmt_difference(
          esd_value, get_reference(esd_reference, "Transition"), region_name == "Transition"
        ),
        esd_subtropical = fmt_difference(
          esd_value, get_reference(esd_reference, "Subtropical"), region_name == "Subtropical"
        ),
        transparency = fmt_number(transparency_value),
        transparency_subarctic = fmt_difference(
          transparency_value, get_reference(transparency_reference, "Subarctic"), region_name == "Subarctic"
        ),
        transparency_transition = fmt_difference(
          transparency_value, get_reference(transparency_reference, "Transition"), region_name == "Transition"
        ),
        transparency_subtropical = fmt_difference(
          transparency_value, get_reference(transparency_reference, "Subtropical"), region_name == "Subtropical"
        ),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }

  comparison_table <- dplyr::bind_rows(comparison_rows)
  readr::write_csv(
    comparison_table,
    file.path(OUT_DIR, "morphological_cluster_regional_comparison_table.csv")
  )

  table_headers <- c(
    "Morphological\nclusters",
    "Average esd\n(mm)",
    "Difference esd\n(mm) to\nSubarctic",
    "Difference esd\n(mm) to\nTransition",
    "Difference esd\n(mm) to\nSubtropical",
    "Average\nTransparency",
    "Difference\ntransparency to\nSubarctic",
    "Difference\ntransparency to\nTransition",
    "Difference\ntransparency to\nSubtropical"
  )
  table_columns <- c(
    "label", "esd", "esd_subarctic", "esd_transition", "esd_subtropical",
    "transparency", "transparency_subarctic", "transparency_transition",
    "transparency_subtropical"
  )
  column_widths <- c(1.55, 1.45, 1.55, 1.55, 1.65, 1.55, 1.65, 1.65, 1.75)
  x_edges <- c(0, cumsum(column_widths))
  x_centers <- head(x_edges, -1) + column_widths / 2
  header_height <- 1.45
  body_height <- 0.49
  n_body <- nrow(comparison_table)
  total_height <- header_height + n_body * body_height

  body_rectangles <- data.frame(
    ymin = header_height + (seq_len(n_body) - 1) * body_height,
    ymax = header_height + seq_len(n_body) * body_height,
    fill = ifelse(comparison_table$row_type == "cluster", "#D0D0D0", "white")
  )

  cell_labels <- do.call(
    rbind,
    lapply(seq_len(n_body), function(i) {
      data.frame(
        column = seq_along(table_columns),
        x = x_centers,
        y = header_height + (i - 0.5) * body_height,
        label = as.character(unlist(comparison_table[i, table_columns], use.names = FALSE)),
        row_type = comparison_table$row_type[i],
        stringsAsFactors = FALSE
      )
    })
  )

  table_plot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = body_rectangles,
      ggplot2::aes(
        xmin = 0, xmax = max(x_edges), ymin = ymin, ymax = ymax, fill = fill
      ),
      colour = NA
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_rect(
      ggplot2::aes(xmin = 0, xmax = max(x_edges), ymin = 0, ymax = header_height),
      fill = "white", colour = NA
    ) +
    ggplot2::geom_text(
      data = data.frame(x = x_centers, label = table_headers),
      ggplot2::aes(x = x, y = header_height / 2, label = label),
      size = 3.45, fontface = "bold", lineheight = 0.95
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cell_labels, column == 1, row_type == "region"),
      ggplot2::aes(x = x_edges[1] + 0.07, y = y, label = label),
      hjust = 0, size = 3.35, fontface = "plain"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cell_labels, column == 1, row_type == "cluster"),
      ggplot2::aes(x = x_edges[1] + 0.07, y = y, label = label),
      hjust = 0, size = 3.35, fontface = "bold"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cell_labels, column != 1, row_type == "region"),
      ggplot2::aes(x = x, y = y, label = label),
      size = 3.25, fontface = "plain"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(cell_labels, column != 1, row_type == "cluster"),
      ggplot2::aes(x = x, y = y, label = label),
      size = 3.25, fontface = "bold"
    ) +
    ggplot2::geom_segment(
      data = data.frame(x = x_edges),
      ggplot2::aes(x = x, xend = x, y = 0, yend = total_height),
      linewidth = 0.25, colour = "#777777"
    ) +
    ggplot2::geom_segment(
      data = data.frame(
        y = c(0, header_height, header_height + seq(4, n_body, by = 4) * body_height)
      ),
      ggplot2::aes(x = 0, xend = max(x_edges), y = y, yend = y),
      linewidth = 0.25, colour = "#777777"
    ) +
    ggplot2::scale_y_reverse(limits = c(total_height, 0), expand = c(0, 0)) +
    ggplot2::scale_x_continuous(limits = range(x_edges), expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 8, 8, 8))

  save_plot_formats(
    filename = file.path(PLOT_DIR, "morphological_cluster_regional_comparison_table.png"),
    plot = table_plot,
    width = 16,
    height = 10.5,
    dpi = 300,
    bg = "white"
  )

  message(
    "Saved formatted comparison table: ",
    file.path(PLOT_DIR, "morphological_cluster_regional_comparison_table.[png/pdf/svg/eps]")
  )
}

# -----------------------------------------------------------------------------
# Combined figures: ONE ROW per figure, one shared legend at the bottom
# Main:     transp_esd_byregion_main.png
# Marginal: transp_esd_byregion_marginal.png
# Each region retains its own top + right marginal distributions.
# -----------------------------------------------------------------------------

# Only the first regional panel supplies the shared cluster legend. Its color
# scale includes all six clusters, even when ST is absent from the Subarctic.
# Suppress guides explicitly: a bottom legend theme propagates into nested
# patchworks and would otherwise reactivate marginal fill legends.
# The original per-region plots keep their own legends for separate exports.
combine_with_one_legend <- function(plots) {
  panels <- lapply(seq_along(plots), function(i) {
    panel <- plots[[i]]
    if (i > 1L) {
      panel <- panel & ggplot2::guides(colour = "none", fill = "none")
    }
    panel
  })
  
  patchwork::wrap_plots(
    panels,
    nrow = 1,
    guides = "collect"
  ) &
    ggplot2::theme(legend.position = "bottom")
}

if (length(transp_esd_main_by_region) > 0) {
  
  combined_all_main <- combine_with_one_legend(transp_esd_main_by_region)
  
  save_plot_formats(
    filename = file.path(PLOT_DIR, "transp_esd_byregion_main.png"),
    plot     = combined_all_main,
    width    = 10 * length(transp_esd_main_by_region),
    height   = 10,
    dpi      = 300,
    bg       = "transparent"
  )
}

if (length(transp_esd_marginal_by_region) > 0) {
  
  combined_all_marg <- combine_with_one_legend(transp_esd_marginal_by_region)
  
  save_plot_formats(
    filename = file.path(PLOT_DIR, "transp_esd_byregion_marginal.png"),
    plot     = combined_all_marg,
    width    = 12 * length(transp_esd_marginal_by_region),
    height   = 12,
    dpi      = 300,
    bg       = "transparent"
  )
}
