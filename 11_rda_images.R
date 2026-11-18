#!/usr/bin/env Rscript
# =============================================================================
# Script:  11_rda_images.R
# Author:  Rocio Rodriguez
# Purpose: Draw original RDA image panels from the saved models in script 10.
# Style:   Original supplied plot construction, dimensions, and combinations.
#          Shared paths/colors, PNG/PDF/SVG export, and plain panel tags.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Uses only the RDA models already fitted and saved by script 10.
# - Builds morphospace panels with representative organism images, abundance-weighted node centroids, morphological loading labels, and environmental vectors.
# - Generates regional and Region x Cluster daytime panels without changing the underlying model fits.
# - Saves the individual panels and the original combined figure layouts in multiple graphics formats.
# -----------------------------------------------------------------------------


source("R_Scripts/00_setup.R", local = TRUE)
PLOT_DIR <- PLOT_DIRS[["RDA"]]
library(readr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggimage)
library(vegan)
library(morphr)
library(png)
library(patchwork)

results <- readRDS(file.path(OUT_DIR, "rda_results.rds"))
dat <- read_csv(file.path(OUT_DIR, "Calanoida_data.csv"), show_col_types = FALSE) |>
  filter(!(Region == "Subarctic" & Cluster == "ST"))
PNG_ROOT_DIR <- file.path(OUT_DIR, "morphr_png_rda")
cluster_order <- NODE_CLUSTER_LEVELS
cluster_icons <- CLUSTER_ICONS
cluster_cols <- CLUSTER_COLORS
# Preserve the RDA figure's original environmental-vector palette.
env_cols <- c(Temp = "red", Sal = "#fc8d62", O2 = "#8da0cb", Fluor = "#66c2a5")
env_labels <- c(Temp = "Temperature (°C)", Sal = "Salinity (PSU)",
                O2 = "Dissolved O₂ (µmol kg⁻¹)", Fluor = "Fluorescence (mg m⁻³)")
USE_FULL_ENV_TEXT <- FALSE
# Original RDA trait-label proxies; these do not affect the fitted model.
composite_map <- c(TRANSPARENT = "mean_YJ", LARGE = "major_YJ", DARK = "stddev_YJ",
                   CIRCULAR = "circ._YJ", SYMMETRICAL = "symetriev_YJ")
crop_bottom_px <- 20

# Use one organism nearest the center of each occupied cell. A 9 × 9 grid
# increases morphospace coverage while matching the 0.09 relative image width.
steps          <- 10
size_scale     <- 0.10
gamma          <- 0.30
# Keep narrow organisms visible without changing their relative width ranking.
min_floor      <- 0.65

arrow_scale    <- 0.16
morph_scale    <- 0.08

env_label_mult  <- 1.2
env_label_nudge <- 0.014

icon_size_min <- 0.06
icon_size_max <- 0.15

full_model_label_size <- 5

# Larger text used only when panels are assembled into combined figures.
# Individual plots retain their original sizes.
combined_full_model_label_size <- 6
combined_morph_label_size      <- 6
combined_env_label_size        <- 6
combined_plot_base_size        <- 26
combined_tag_size              <- 26


grid_reps <- function(df_scores, steps, axes = c(1, 2)) {
  xr <- range(df_scores$Axis1, na.rm = TRUE)
  yr <- range(df_scores$Axis2, na.rm = TRUE)
  
  xb <- seq(xr[1], xr[2], length.out = steps + 1)
  yb <- seq(yr[1], yr[2], length.out = steps + 1)
  
  cx <- (xb[-1] + xb[-length(xb)]) / 2
  cy <- (yb[-1] + yb[-length(yb)]) / 2
  
  df_bins <- df_scores %>%
    dplyr::mutate(
      bin_x = cut(Axis1, xb, include.lowest = TRUE, labels = FALSE),
      bin_y = cut(Axis2, yb, include.lowest = TRUE, labels = FALSE),
      cen_x = cx[bin_x],
      cen_y = cy[bin_y],
      dist2 = (Axis1 - cen_x)^2 + (Axis2 - cen_y)^2
    )
  
  df_bins %>%
    dplyr::filter(!is.na(bin_x), !is.na(bin_y)) %>%
    dplyr::group_by(bin_x, bin_y) %>%
    dplyr::slice_min(dist2, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}


process_pngs <- function(df_sel, out_dir, bottom = 20) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  df_sel <- df_sel %>%
    dplyr::mutate(png_path = file.path(out_dir, paste0(id, ".png")))
  
  for (i in seq_len(nrow(df_sel))) {
    in_jpg  <- df_sel$img_path[i]
    out_png <- df_sel$png_path[i]
    
    if (file.exists(out_png)) next
    
    x <- morphr::img_read(in_jpg)
    
    dims <- dim(x)
    w <- dims[1]
    h <- dims[2]
    
    x <- x[1:w, 1:max(1, h - bottom), , , drop = FALSE] %>%
      morphr::img_extract_largest(quiet = TRUE) %>%
      morphr::img_invert() %>%
      morphr::img_make_horizontal() %>%
      morphr::img_invert() %>%
      morphr::img_make_transparent()
    
    morphr::img_write(x, out_png)
    
    rm(x)
    
    if (i %% 10 == 0) gc()
  }
  
  df_sel
}

png_bbox_width <- function(path) {
  arr <- png::readPNG(path)
  
  if (length(dim(arr)) == 3 && dim(arr)[3] == 4) {
    alpha <- arr[, , 4]
    cols  <- which(colSums(alpha > 0) > 0)
    
    if (length(cols)) {
      return(max(cols) - min(cols) + 1)
    }
  }
  
  dim(arr)[2]
}

compute_sizes <- function(df_sel, size_scale, gamma = 1, min_floor = 0) {
  w <- vapply(df_sel$png_path, png_bbox_width, numeric(1))
  r <- w / median(w)
  r_adj <- pmax(r^gamma, min_floor)
  
  dplyr::mutate(df_sel, size_plot = size_scale * r_adj)
}

env_vectors <- function(fit, df_env, arrow_scale, axes = c(1, 2)) {
  env_df <- as.data.frame(df_env[, names(env_cols), drop = FALSE])
  rownames(env_df) <- df_env$id
  
  ef <- vegan::envfit(
    fit,
    env_df,
    choices = axes,
    permutations = 0
  )
  out <- as.data.frame(ef$vectors$arrows)
  
  out$var <- rownames(out)
  out <- out[, c(1, 2, 3), drop = FALSE]
  colnames(out)[1:2] <- c("Axis1", "Axis2")
  
  out$Axis1 <- out$Axis1 * arrow_scale
  out$Axis2 <- out$Axis2 * arrow_scale
  out$var <- factor(out$var, levels = names(env_cols))
  
  out$lab_x <- out$Axis1 * env_label_mult + env_label_nudge * sign(out$Axis1)
  out$lab_y <- out$Axis2 * env_label_mult + env_label_nudge * sign(out$Axis2)
  
  out
}


morph_loadings <- function(fit, composite_map, morph_scale, axes = c(1, 2)) {
  L <- as.data.frame(
    vegan::scores(fit, display = "species", choices = axes, scaling = 1)
  )
  
  colnames(L) <- c("Axis1", "Axis2")
  
  P <- L[unname(composite_map), , drop = FALSE]
  rownames(P) <- names(composite_map)
  P <- rbind(P, SMALL = -P["LARGE", , drop = FALSE])
  
  P$label <- rownames(P)
  P$Axis1 <- P$Axis1 * morph_scale
  P$Axis2 <- P$Axis2 * morph_scale
  
  P
}


# ---- PLOT FUNCTION ----------------------------------------------------------

plot_region <- function(sel,
                        fit,
                        env_fit,
                        pca_vars,
                        cluster_centers,
                        axes = c(1, 2),
                        base_size = 20,
                        full_model_size = full_model_label_size,
                        morph_label_size = 4,
                        env_label_size = 4) {
  
  ev_constrained <- vegan::eigenvals(fit, model = "constrained")
  if (length(ev_constrained) < max(axes)) {
    stop("Requested axes ", paste(axes, collapse = " and "), " are not available.")
  }
  
  pct_constrained <- round(
    100 * ev_constrained[axes] / sum(ev_constrained),
    1
  )
  
  full_model_r2_pct <- 100 * fit$CCA$tot.chi / fit$tot.chi
  
  x_lab <- sprintf(
    "RDA%d (%.1f%% of constrained variance)",
    axes[1], pct_constrained[1]
  )
  y_lab <- sprintf(
    "RDA%d (%.1f%% of constrained variance)",
    axes[2], pct_constrained[2]
  )
  
  cluster_legend_df <- tibble::tibble(
    Cluster = factor(cluster_order, levels = cluster_order),
    x = 999,
    y = 999
  )
  
  env_legend_df <- tibble::tibble(
    var = factor(names(env_cols), levels = names(env_cols)),
    x = 999,
    y = 999
  )
  
  env_text <- if (USE_FULL_ENV_TEXT) {
    env_labels[as.character(env_fit$var)]
  } else {
    as.character(env_fit$var)
  }
  
  ggplot(sel, aes(Axis1, Axis2)) +
    
    ggimage::geom_image(
      aes(image = png_path, size = size_plot),
      by = "width",
      asp = 1,
      show.legend = FALSE
    ) +
    
    scale_size_identity() +
    coord_equal() +
    labs(
      x = x_lab,
      y = y_lab
    ) +
    annotate(
      "label",
      x = -0.6, y = 0.6,
      label = sprintf("Full model R² = %.2f%%", full_model_r2_pct),
      hjust = 0, vjust = 1,
      size = full_model_size,
      fontface = "bold",
      fill = "black",
      color = "white",
      label.size = 0.25
    ) +
    theme_classic(base_size = base_size) +
    theme(legend.position = "bottom") +
    
    geom_point(
      data = cluster_legend_df,
      aes(x = x, y = y, fill = Cluster),
      shape = 21,
      size = 5,
      color = "black",
      alpha = 0,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    
    scale_fill_manual(
      values = cluster_cols,
      drop = FALSE,
      name = "Cluster"
    ) +
    
    ggimage::geom_image(
      data = cluster_centers,
      aes(x = Axis1, y = Axis2, image = icon_path, size = icon_size),
      inherit.aes = FALSE,
      asp = 1,
      show.legend = FALSE
    ) +
    
    geom_segment(
      data = env_fit,
      aes(x = 0, y = 0, xend = Axis1, yend = Axis2, color = var),
      arrow = grid::arrow(length = grid::unit(0.15, "inches")),
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    geom_text(
      data = env_fit,
      aes(x = lab_x, y = lab_y, color = var),
      label = env_text,
      inherit.aes = FALSE,
      size = env_label_size,
      fontface = "bold",
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    
    geom_point(
      data = env_legend_df,
      aes(x = x, y = y, color = var),
      shape = 16,
      size = 5,
      alpha = 0,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    
    scale_color_manual(
      values = env_cols,
      labels = env_labels,
      drop = FALSE,
      name = "Environment"
    ) +
    
    guides(
      color = guide_legend(
        order = 1,
        nrow = 2,
        byrow = TRUE,
        override.aes = list(alpha = 1, shape = 16, size = 5)
      ),
      fill = guide_legend(
        order = 2,
        override.aes = list(alpha = 1, size = 5, shape = 21)
      )
    ) +
    
    geom_label(
      data = pca_vars,
      aes(x = Axis1, y = Axis2, label = label),
      inherit.aes = FALSE,
      fill = "white",
      color = "black",
      size = morph_label_size,
      label.size = 0.25,
      fontface = "bold",
      check_overlap = TRUE
    ) +
    
    scale_x_continuous(limits = c(-0.6, 0.6)) +
    scale_y_continuous(limits = c(-0.6, 0.6))
}


plot_region_cluster <- function(sel, fit, env_fit, pca_vars,
                                cluster_name, cluster_center,
                                base_size = 20) {
  
  ev <- vegan::eigenvals(fit, model = "constrained")
  
  if (length(ev) < 2) {
    stop("This RDA has fewer than two constrained axes.")
  }
  
  pct <- round(100 * ev[1:2] / sum(ev), 1)
  
  full_model_r2_pct <- 100 * fit$CCA$tot.chi / fit$tot.chi
  
  x_lab <- sprintf("RDA1 (%.1f%% constrained variance)", pct[1])
  y_lab <- sprintf("RDA2 (%.1f%% constrained variance)", pct[2])
  
  env_legend_df <- tibble::tibble(
    var = factor(names(env_cols), levels = names(env_cols)),
    x = 999,
    y = 999
  )
  
  env_text <- if (USE_FULL_ENV_TEXT) {
    env_labels[as.character(env_fit$var)]
  } else {
    as.character(env_fit$var)
  }
  
  ggplot2::ggplot(sel, ggplot2::aes(Axis1, Axis2)) +
    
    ggimage::geom_image(
      ggplot2::aes(image = png_path, size = size_plot),
      by = "width",
      asp = 1,
      show.legend = FALSE
    ) +
    
    ggplot2::scale_size_identity() +
    ggplot2::coord_equal() +
    ggplot2::labs(x = x_lab, y = y_lab) +
    ggplot2::annotate(
      "label",
      x = -0.6, y = 0.6,
      label = sprintf("Full model R² = %.2f%%", full_model_r2_pct),
      hjust = 0, vjust = 1,
      size = full_model_label_size,
      fontface = "bold",
      fill = "black",
      color = "white",
      label.size = 0.25
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(legend.position = "bottom") +
    
    # Cluster icon at the abundance-weighted centroid of this node
    ggimage::geom_image(
      data = cluster_center,
      ggplot2::aes(
        x = Axis1,
        y = Axis2,
        image = icon_path,
        size = icon_size
      ),
      inherit.aes = FALSE,
      asp = 1,
      show.legend = FALSE
    ) +
    
    # Environmental vectors
    ggplot2::geom_segment(
      data = env_fit,
      ggplot2::aes(
        x = 0, y = 0,
        xend = Axis1, yend = Axis2,
        color = var
      ),
      arrow = grid::arrow(length = grid::unit(0.15, "inches")),
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    
    ggplot2::geom_text(
      data = env_fit,
      ggplot2::aes(x = lab_x, y = lab_y, color = var),
      label = env_text,
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold",
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    
    # Environment legend
    ggplot2::geom_point(
      data = env_legend_df,
      ggplot2::aes(x = x, y = y, color = var),
      shape = 16,
      size = 5,
      alpha = 0,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    
    ggplot2::scale_color_manual(
      values = env_cols,
      labels = env_labels,
      drop = FALSE,
      name = "Environment"
    ) +
    
    ggplot2::guides(
      color = ggplot2::guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(alpha = 1, shape = 16, size = 5)
      )
    ) +
    
    # Morphological loading labels
    ggplot2::geom_label(
      data = pca_vars,
      ggplot2::aes(x = Axis1, y = Axis2, label = label),
      inherit.aes = FALSE,
      fill = "white",
      color = "black",
      size = 4.0,
      label.size = 0.25,
      fontface = "bold",
      check_overlap = TRUE
    ) +
    
    ggplot2::scale_x_continuous(limits = c(-0.6, 0.6)) +
    ggplot2::scale_y_continuous(limits = c(-0.6, 0.6))
}



# Prepare only graphical inputs from the fitted model; no refitting here.
prepare_panel <- function(result, axes, png_directory, node = FALSE) {
  fit <- result$fit
  scores <- as.data.frame(vegan::scores(fit, display = "sites", choices = axes, scaling = 1))
  names(scores) <- c("Axis1", "Axis2")
  scores$id <- result$ids
  metadata <- dat[match(result$ids, dat$id), ]
  scores <- left_join(scores, select(metadata, id, Cluster, abundance, esd), by = "id")
  scores$Cluster <- factor(scores$Cluster, levels = cluster_order)
  scores$img_path <- file.path(IMAGES_DIR, paste0(scores$id, ".jpg"))
  selected <- grid_reps(scores, steps, axes)
  selected <- process_pngs(selected, png_directory, crop_bottom_px)
  selected <- compute_sizes(selected, size_scale, gamma, min_floor)
  vectors <- env_vectors(fit, select(metadata, id, all_of(ENV_RAW)), arrow_scale, axes)
  traits <- morph_loadings(fit, composite_map, morph_scale, axes)
  centers <- scores |>
    group_by(Cluster) |>
    summarise(Axis1 = weighted.mean(Axis1, abundance),
              Axis2 = weighted.mean(Axis2, abundance),
              mean_esd = weighted.mean(esd, abundance), .groups = "drop") |>
    mutate(icon_path = unname(cluster_icons[as.character(Cluster)]))
  # A node has one icon: the original explicitly uses the middle icon size.
  centers$icon_size <- if (node) rep((icon_size_min + icon_size_max) / 2, nrow(centers)) else
    icon_sizes(centers$mean_esd, icon_size_min, icon_size_max)
  list(sel = selected, fit = fit, env_fit = vectors, pca_vars = traits, cluster_centers = centers)
}

combined_theme <- theme(
  legend.position = "bottom",
  axis.title = element_text(size = combined_plot_base_size),
  axis.text = element_text(size = combined_plot_base_size - 4),
  legend.title = element_text(size = combined_plot_base_size - 2),
  legend.text = element_text(size = combined_plot_base_size - 4),
  plot.tag = element_text(face = "plain", size = combined_tag_size)
)
combine_rda <- function(panels, ncol, nrow) {
  (wrap_plots(panels, ncol = ncol, nrow = nrow, byrow = TRUE) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "A", theme = theme(
      plot.tag = element_text(face = "plain", size = combined_tag_size),
      legend.position = "bottom"))) & combined_theme
}

# Regional individual panels, regional rows, and the original 2 × 3 combinations.
all_panels <- list()
axis_pairs <- list(RDA12 = c(1, 2), RDA34 = c(3, 4))
for (diel in names(DIEL_RUNS)) {
  all_panels[[diel]] <- list()
  for (pair_name in names(axis_pairs)) {
    axes <- axis_pairs[[pair_name]]
    directory <- file.path(PLOT_DIR, diel, pair_name)
    combined_panels <- setNames(vector("list", length(REGION_LEVELS)), REGION_LEVELS)
    for (region in REGION_LEVELS) {
      args <- prepare_panel(results$regional[[diel]][[region]], axes,
                            file.path(PNG_ROOT_DIR, diel, pair_name, region))
      individual <- do.call(plot_region, c(args, list(axes = axes, base_size = 20,
        full_model_size = full_model_label_size, morph_label_size = 4, env_label_size = 4)))
      combined_panels[[region]] <- do.call(plot_region, c(args, list(axes = axes,
        base_size = combined_plot_base_size, full_model_size = combined_full_model_label_size,
        morph_label_size = combined_morph_label_size, env_label_size = combined_env_label_size)))
      save_plot_formats(file.path(directory, paste0("RDA_", diel, "_", region, "_", pair_name, ".png")),
                        individual, width = 12, height = 8, dpi = 300)
    }
    save_plot_formats(file.path(directory, paste0("RDA_", diel, "_ALLregions_combined_", pair_name, ".png")),
                      combine_rda(combined_panels, 3, 1), width = 36, height = 9, dpi = 300)
    all_panels[[diel]][[pair_name]] <- combined_panels
  }
  six <- c(all_panels[[diel]]$RDA12, all_panels[[diel]]$RDA34)
  save_plot_formats(file.path(PLOT_DIR, "combined", paste0("RDA12_RDA34_", diel, "_2x3_A-F.png")),
                    combine_rda(six, 3, 2), width = 30, height = 18, dpi = 300)
}

# Original nine-panel combination. The source used different PNG/PDF sizes;
# preserve each; SVG follows the original PNG dimensions and EPS follows PDF.
nine <- c(all_panels$Day$RDA12, all_panels$Night$RDA12, all_panels$DayNight$RDA12)
combined_nine <- combine_rda(nine, 3, 3)
combined_path <- file.path(PLOT_DIR, "combined", "RDA12_Day_Night_DayNight_3x3_A-I.png")
save_plot_formats(combined_path, combined_nine, width = 40, height = 35, dpi = 300,
                  formats = c("png", "svg"))
save_plot_formats(combined_path, combined_nine, width = 30, height = 27, dpi = 300,
                  formats = c("pdf", "eps"))

# Daytime node panels reuse their own saved model; no extra night/combined fits.
for (node in names(results$nodes)) {
  result <- results$nodes[[node]]
  args <- prepare_panel(result, c(1, 2), file.path(PNG_ROOT_DIR, "Day", "Region_x_Cluster", node), node = TRUE)
  panel <- plot_region_cluster(args$sel, args$fit, args$env_fit, args$pca_vars,
                               result$model$Cluster, args$cluster_centers)
  save_plot_formats(file.path(PLOT_DIR, "Day", "Region_x_Cluster", paste0("RDA_Day_", node, ".png")),
                    panel, width = 12, height = 8, dpi = 300)
}
