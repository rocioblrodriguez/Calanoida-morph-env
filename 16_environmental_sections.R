#!/usr/bin/env Rscript
# Script:  16_environmental_sections.R
# =============================================================================
# Latitude-depth environmental sections, with sampling dots colored by Region
# Author: Rocio Rodriguez
#
# Input: outputs/data_gradients.csv, created by
#        01.Download&Process_Gradients.R
#
# Layout: Temperature | Fluorescence
#         Oxygen      | Salinity
#
# Creates a separate four-panel figure for each cruise, using ALL three regions.
# Background colors show environmental values; dot colors show Region.
# Environmental colorbars stay beside their panels; the Region legend is shared.
#
# IMPORTANT ABOUT THE BACKGROUND:
# This script uses the CTD means already calculated for the eight MOCNESS depth
# bins, rather than rereading raw CTD casts. Each station/bin contributes once,
# regardless of the number of EcoTaxa images. By default, shape-preserving cubic
# interpolation (monoH.FC) is applied between bin midpoints, then between station
# latitudes, giving curved contours through the input values. Values at the
# shallowest and deepest midpoints are held constant to their bin edges. There
# is no extrapolation beyond sampled latitude/depth coverage. Cruises are kept
# separate because they were sampled in different years.
# Environmental values are not averaged across day/night or MOCNESS tows here:
# those tows already share station/net CTD means from the processing script.
# =============================================================================

# 1. USER SETTINGS ------------------------------------------------------------

# -----------------------------------------------------------------------------
# Workflow notes:
# - Builds latitude-depth sections for temperature, fluorescence, dissolved oxygen, and salinity separately for each cruise.
# - Uses the station/net CTD means already generated in script 01, after removing organism-level repetition so image abundance does not weight the environment.
# - Interpolates first through depth and then across station latitude only within sampled coverage.
# - Overlays sampling locations colored by region and uses fixed environmental scales to support comparison between cruises.
# -----------------------------------------------------------------------------

source("R_Scripts/00_setup.R", local = TRUE)
INPUT_CSV <- file.path(OUT_DIR, "data_gradients.csv")
PLOT_DIR <- PLOT_DIRS[["Environmental_sections"]]
OUT_DIR <- PLOT_DIR

# NULL = every cruise in the input; alternatively c("nh1208") or c("oc473").
CRUISES_TO_PLOT <- NULL

# Dot colors follow the Region column in your processing script.
REGION_COLORS <- c(
  "Subarctic"   = "#0000FF",
  "Transition"  = "#FFD700",
  "Subtropical" = "#FF0000"
)

# Low environmental values = purple; high values = red, as in the example.
SECTION_COLORS <- c(
  "#AD4FFF", "#763CFF", "#3E80FF", "#35C8DB",
  "#3DC56B", "#E2DF48", "#FFB541", "#F23E32"
)

BASE_SIZE         <- 14
TITLE_SIZE        <- 15
AXIS_TITLE_SIZE   <- 13
AXIS_TEXT_SIZE    <- 11
DOT_SIZE          <- 2.2
DOT_HALO          <- 0.4  # thin white outline for visibility; 0 disables it
CONTOUR_TEXT_SIZE <- 2.8  # ggplot text size, in mm
SHOW_CONTOURS     <- TRUE
SHOW_CONTOUR_TEXT <- TRUE
ADD_PANEL_TAGS    <- FALSE
SAVE_INDIVIDUAL   <- FALSE

FIG_WIDTH_IN  <- 12
FIG_HEIGHT_IN <- 9
PNG_DPI       <- 400

# "smooth" = curved, shape-preserving interpolation; "linear" = original method.
INTERPOLATION <- "smooth"
N_LATITUDE    <- 600
DEPTH_STEP    <- 1  # interpolation grid, in metres

# Extra smoothing of raster colors during display only.
# FALSE removes this extra display blur, keeping the interpolated values and
# smooth contour positions unchanged.
RASTER_INTERPOLATE <- TRUE

# Fixed environmental color scales, identical for both cruises.
# Colorbars show the maximum at the top and minimum at the bottom.
# Out-of-range values use the nearest end color; input values are unchanged.
COLOR_LIMITS <- list(
  Temp  = c(0, 20),
  Fluor = c(0, 2.5),
  O2    = c(0, 300),
  Sal   = c(25, 35)
)
COLOR_BREAKS <- list(
  Temp  = seq(0, 20, by = 5),
  Fluor = seq(0, 2.5, by = 0.5),
  O2    = seq(0, 300, by = 50),
  Sal   = seq(25, 35, by = 2.5)
)

# Edit these to control which contour lines/labels appear.
CONTOUR_BREAKS <- list(
  Temp  = seq(0, 20, by = 5),
  Fluor = c(0.1, 0.2, 0.5, 1, 1.5, 2, 2.5),
  O2    = c(0, 20, 50, 100, 150, 200, 250, 300),
  Sal   = seq(25, 35, by = 0.5)
)

# 2. PACKAGES -----------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "ggplot2", "cowplot")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs)) {
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}
invisible(lapply(required_pkgs, library, character.only = TRUE))

if (length(INTERPOLATION) != 1L || !INTERPOLATION %in% c("smooth", "linear")) {
  stop('INTERPOLATION must be "smooth" or "linear".')
}

if (!file.exists(INPUT_CSV)) {
  stop("Input not found: ", INPUT_CSV,
       "\nRun your download/processing script first, or edit INPUT_CSV.")
}
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# 3. READ ONLY THE COLUMNS NEEDED FOR THESE FIGURES -----------------------------
env_vars <- c("Temp", "Fluor", "O2", "Sal")
required_cols <- c("Cruise", "Station", "Region", "net", "lat", env_vars)
header <- names(readr::read_csv(INPUT_CSV, n_max = 0, show_col_types = FALSE))
missing_cols <- setdiff(required_cols, header)
if (length(missing_cols)) {
  stop("Missing columns in INPUT_CSV: ", paste(missing_cols, collapse = ", "))
}

dat <- readr::read_csv(
  INPUT_CSV,
  col_select = dplyr::all_of(required_cols),
  col_types = readr::cols(
    Cruise = readr::col_character(), Station = readr::col_character(),
    Region = readr::col_character(), net = readr::col_character(),
    lat = readr::col_double(), Temp = readr::col_double(),
    Fluor = readr::col_double(), O2 = readr::col_double(),
    Sal = readr::col_double()
  ),
  show_col_types = FALSE
)

# Use the SAME depth-bin definitions as your CTD processing script.
net_bins <- data.frame(
  net = c("n8", "n7", "n6", "n5", "n4", "n3", "n2", "n1"),
  depth_min = c(0, 25, 50, 100, 200, 400, 600, 800),
  depth_max = c(25, 50, 100, 200, 400, 600, 800, 1000)
)
net_bins$depth_mid <- (net_bins$depth_min + net_bins$depth_max) / 2

n_before <- nrow(dat)
dat <- dat |>
  dplyr::filter(
    !is.na(Cruise), !is.na(Station), is.finite(lat),
    Region %in% names(REGION_COLORS), net %in% net_bins$net
  )
if (nrow(dat) < n_before) {
  warning(n_before - nrow(dat), " rows excluded because sampling metadata were missing or unrecognized.")
}
if (!nrow(dat)) stop("No usable sampling records were found.")

# One latitude per cruise/station: mean of distinct recorded tow latitudes.
# Deduplicating first prevents EcoTaxa image counts from weighting this mean.
station_meta <- dat |>
  dplyr::distinct(Cruise, Station, Region, lat) |>
  dplyr::group_by(Cruise, Station, Region) |>
  dplyr::summarise(latitude = mean(lat), .groups = "drop")
if (anyDuplicated(station_meta[c("Cruise", "Station")])) {
  stop("A cruise/station has more than one Region assignment. Check the input.")
}

# The original processing script assigns one environmental record to each
# cruise/station/net. Remove organism-level repetitions BEFORE interpolation.
env_records <- dat |>
  dplyr::distinct(Cruise, Station, net, Temp, Fluor, O2, Sal)
if (anyDuplicated(env_records[c("Cruise", "Station", "net")])) {
  stop("Multiple environmental records found for the same cruise/station/net.\n",
       "Use the original data_gradients.csv from your processing script.")
}
samples <- env_records |>
  dplyr::left_join(station_meta, by = c("Cruise", "Station")) |>
  dplyr::left_join(net_bins, by = "net") |>
  dplyr::mutate(Region = factor(Region, levels = names(REGION_COLORS))) |>
  dplyr::arrange(Cruise, latitude, depth_mid)

available_cruises <- sort(unique(samples$Cruise))
if (is.null(CRUISES_TO_PLOT)) CRUISES_TO_PLOT <- available_cruises
if (!length(CRUISES_TO_PLOT) || any(!CRUISES_TO_PLOT %in% available_cruises)) {
  stop("Check CRUISES_TO_PLOT. Available cruises: ",
       paste(available_cruises, collapse = ", "))
}

# 4. PANEL SETTINGS -----------------------------------------------------------
panels <- list(
  list(key = "Temp", title = expression("Temperature (" * degree * "C)"),
       max_depth = 1000, depth_breaks = seq(0, 1000, 250)),
  list(key = "Fluor", title = expression("Fluorescence (mg m"^{-3} * ")"),
       max_depth = 200, depth_breaks = seq(0, 200, 50)),
  list(key = "O2", title = expression("Oxygen (" * mu * "mol kg"^{-1} * ")"),
       max_depth = 1000, depth_breaks = seq(0, 1000, 250)),
  list(key = "Sal", title = "Salinity (PSU)",
       max_depth = 1000, depth_breaks = seq(0, 1000, 250))
)

# Check the requested scale limits before building the plots.
for (variable in env_vars) {
  lim <- COLOR_LIMITS[[variable]]
  if (length(lim) != 2L || any(!is.finite(lim)) || lim[1] >= lim[2]) {
    stop("COLOR_LIMITS must contain two increasing finite values for ", variable)
  }
}

# 5. INTERPOLATE VERTICALLY, THEN ACROSS LATITUDES ------------------------------
# Use the same interpolation for the background and its contour lines.
# Restrict evaluation to observed coordinates: splinefun would otherwise
# extrapolate by default. With only two distinct coordinates, use a straight line.
interpolate_axis <- function(x, values, xout) {
  keep <- is.finite(x) & is.finite(values)
  x <- x[keep]
  values <- values[keep]
  result <- rep(NA_real_, length(xout))
  n_unique <- length(unique(x))
  if (n_unique < 2L) return(result)
  inside <- is.finite(xout) & xout >= min(x) & xout <= max(x)
  if (!any(inside)) return(result)

  if (INTERPOLATION == "smooth" && n_unique >= 3L) {
    curve <- stats::splinefun(x, values, method = "monoH.FC", ties = mean)
    result[inside] <- curve(xout[inside])
  } else {
    result[inside] <- stats::approx(
      x, values, xout = xout[inside], rule = 1, ties = mean
    )$y
  }
  result
}

make_section <- function(cruise_samples, variable, max_depth) {
  x <- seq(min(cruise_samples$latitude), max(cruise_samples$latitude),
           length.out = N_LATITUDE)
  y <- sort(unique(c(seq(0, max_depth, by = DEPTH_STEP), max_depth)))
  stations <- split(cruise_samples, cruise_samples$Station)
  station_lat <- vapply(stations, function(s) s$latitude[1], numeric(1))

  profiles <- vapply(stations, function(s) {
    s <- s[is.finite(s[[variable]]), , drop = FALSE]
    s <- s[order(s$depth_mid), , drop = FALSE]
    if (!nrow(s)) return(rep(NA_real_, length(y)))
    # Constant extension is restricted to the end bins' own boundaries.
    yy <- c(s$depth_min[1], s$depth_mid, tail(s$depth_max, 1))
    vv <- c(s[[variable]][1], s[[variable]], tail(s[[variable]], 1))
    interpolate_axis(yy, vv, xout = y)
  }, numeric(length(y)))

  if (length(unique(station_lat)) < 2) {
    stop("At least two different station latitudes are needed for a section.")
  }
  z <- matrix(NA_real_, nrow = length(x), ncol = length(y))
  for (j in seq_along(y)) {
    ok <- is.finite(profiles[j, ])
    if (length(unique(station_lat[ok])) >= 2) {
      z[, j] <- interpolate_axis(station_lat[ok], profiles[j, ok], xout = x)
    }
  }
  if (!any(is.finite(z))) stop("Insufficient profile coverage for ", variable)
  raster <- expand.grid(latitude = x, depth = y)
  raster$value <- as.vector(z)
  list(x = x, y = y, z = z, raster = raster)
}

# Contours and one label per level, using base R (no extra contour packages).
make_contours <- function(section, levels) {
  empty <- data.frame(latitude = numeric(), depth = numeric(),
                      level = numeric(), piece = integer())
  zrange <- range(section$z, na.rm = TRUE)
  levels <- levels[levels > zrange[1] & levels < zrange[2]]
  if (!length(levels)) return(list(lines = empty, labels = empty))
  paths <- grDevices::contourLines(
    x = section$x, y = section$y, z = section$z, levels = levels
  )
  if (!length(paths)) return(list(lines = empty, labels = empty))
  lines <- dplyr::bind_rows(lapply(seq_along(paths), function(i) {
    data.frame(latitude = paths[[i]]$x, depth = paths[[i]]$y,
               level = paths[[i]]$level, piece = i)
  }))
  xr <- diff(range(section$x))
  yr <- diff(range(section$y))
  lengths <- vapply(paths, function(p) {
    sum(sqrt((diff(p$x) / xr)^2 + (diff(p$y) / yr)^2))
  }, numeric(1))
  best <- vapply(split(seq_along(paths), vapply(paths, function(p) p$level, numeric(1))),
                 function(ids) ids[which.max(lengths[ids])], integer(1))
  labels <- dplyr::bind_rows(lapply(best, function(i) {
    p <- paths[[i]]
    arc <- c(0, cumsum(sqrt((diff(p$x) / xr)^2 + (diff(p$y) / yr)^2)))
    k <- which.min(abs(arc - max(arc) / 2))
    data.frame(latitude = p$x[k], depth = p$y[k],
               level = p$level, piece = i)
  }))
  # Keep text off the panel borders; retain all contour lines.
  labels <- labels |>
    dplyr::filter(
      latitude > min(section$x) + 0.03 * xr,
      latitude < max(section$x) - 0.03 * xr,
      depth > min(section$y) + 0.03 * yr,
      depth < max(section$y) - 0.03 * yr
    )
  list(lines = lines, labels = labels)
}

# 6. DRAW ONE PANEL -----------------------------------------------------------
make_panel <- function(cruise_samples, spec, tag) {
  section <- make_section(cruise_samples, spec$key, spec$max_depth)
  contours <- make_contours(section, CONTOUR_BREAKS[[spec$key]])
  dots <- cruise_samples |>
    dplyr::filter(depth_mid <= spec$max_depth, is.finite(.data[[spec$key]]))

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = section$raster,
      ggplot2::aes(x = latitude, y = depth, fill = value),
      interpolate = RASTER_INTERPOLATE
    ) +
    ggplot2::scale_fill_gradientn(
      colors = SECTION_COLORS, limits = COLOR_LIMITS[[spec$key]],
      breaks = COLOR_BREAKS[[spec$key]],
      oob = scales::squish, na.value = "white", name = NULL,
      guide = ggplot2::guide_colorbar(
        barheight = grid::unit(36, "mm"), barwidth = grid::unit(4, "mm")
      )
    )
  if (SHOW_CONTOURS && nrow(contours$lines)) {
    p <- p + ggplot2::geom_path(
      data = contours$lines,
      ggplot2::aes(x = latitude, y = depth, group = piece),
      color = "black", linewidth = 0.3, alpha = 0.75,
      linejoin = "round", lineend = "round"
    )
    if (SHOW_CONTOUR_TEXT && nrow(contours$labels)) {
      p <- p + ggplot2::geom_text(
        data = contours$labels,
        ggplot2::aes(x = latitude, y = depth, label = level),
        size = CONTOUR_TEXT_SIZE, check_overlap = TRUE
      )
    }
  }
  if (DOT_HALO > 0) {
    p <- p + ggplot2::geom_point(
      data = dots, ggplot2::aes(x = latitude, y = depth_mid),
      color = "white", size = DOT_SIZE + DOT_HALO
    )
  }
  p +
    ggplot2::geom_point(
      data = dots,
      ggplot2::aes(x = latitude, y = depth_mid, color = Region),
      size = DOT_SIZE
    ) +
    ggplot2::scale_color_manual(values = REGION_COLORS, drop = FALSE, guide = "none") +
    ggplot2::scale_x_continuous(
      limits = range(section$x), breaks = pretty(range(section$x), n = 6),
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::scale_y_reverse(
      limits = c(spec$max_depth, 0), breaks = spec$depth_breaks,
      expand = ggplot2::expansion(mult = 0)
    ) +
    ggplot2::labs(
      title = spec$title, x = expression("Latitude (" * degree * "N)"),
      y = "Depth (m)", tag = if (ADD_PANEL_TAGS) tag else NULL
    ) +
    ggplot2::theme_bw(base_size = BASE_SIZE) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 0.5),
      plot.title = ggplot2::element_text(size = TITLE_SIZE, hjust = 0, margin = ggplot2::margin(b = 8)),
      plot.tag = ggplot2::element_text(size = TITLE_SIZE + 1, face = "plain"),
      axis.title = ggplot2::element_text(size = AXIS_TITLE_SIZE),
      axis.text = ggplot2::element_text(size = AXIS_TEXT_SIZE, color = "black"),
      legend.position = "right", legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(10, 10, 5, 10)
    )
}

# Build a single Region legend independently, leaving panel colorbars in place.
legend_data <- data.frame(
  x = seq_along(REGION_COLORS), y = 1,
  Region = factor(names(REGION_COLORS), levels = names(REGION_COLORS))
)
region_legend <- cowplot::get_legend(
  ggplot2::ggplot(legend_data, ggplot2::aes(x, y, color = Region)) +
    ggplot2::geom_point(size = 3.5) +
    ggplot2::scale_color_manual(values = REGION_COLORS, drop = FALSE, name = "Region") +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::theme_void(base_size = BASE_SIZE) +
    ggplot2::theme(legend.position = "bottom")
)

# 7. SAVE FIGURES -------------------------------------------------------------
save_figure <- function(plot, basename, width, height) {
  png_path <- file.path(OUT_DIR, paste0(basename, ".png"))
  save_plot_formats(png_path, plot, width = width, height = height,
                    units = "in", dpi = PNG_DPI, bg = "white")
}

for (cruise_id in CRUISES_TO_PLOT) {
  message("Creating environmental sections for ", cruise_id, " ...")
  cruise_samples <- samples[samples$Cruise == cruise_id, , drop = FALSE]
  plots <- lapply(seq_along(panels), function(i) {
    make_panel(cruise_samples, panels[[i]], LETTERS[i])
  })
  four_panels <- cowplot::plot_grid(
    plotlist = plots, ncol = 2, align = "hv", axis = "tblr"
  )
  combined <- cowplot::plot_grid(
    four_panels, region_legend, ncol = 1, rel_heights = c(1, 0.065)
  )
  safe_cruise <- gsub("[^A-Za-z0-9_-]", "_", cruise_id)
  save_figure(combined, paste0("Environmental_sections_", safe_cruise, "_4panels"),
              FIG_WIDTH_IN, FIG_HEIGHT_IN)
  if (SAVE_INDIVIDUAL) {
    for (i in seq_along(plots)) {
      individual <- cowplot::plot_grid(
        plots[[i]], region_legend, ncol = 1, rel_heights = c(1, 0.12)
      )
      save_figure(individual,
                  paste0("Environmental_", safe_cruise, "_", panels[[i]]$key),
                  width = 7, height = 5.5)
    }
  }
}

message("Done. Figures are in: ", OUT_DIR)

