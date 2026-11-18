#!/usr/bin/env Rscript
# Script:  15_sampling_map.R
# -----------------------------------------------------------------------------
# Title:     Northeast Pacific ocean currents, oceanographic regions,
#            and NH1208 MOCNESS sampling stations
# Author:    Rocio Rodriguez
# Date:      2026-10-19
#
# PURPOSE:
#   Create a map of the Northeast Pacific showing the oceanographic setting
#   and sampling locations used for the NH1208 cruise.
#
# DATA:
#   A) SAMPLING STATIONS:
#        - Reads station coordinates from data_gradients.csv.
#        - Uses NH1208 MOCNESS stations only.
#        - Stations are classified into three oceanographic regions:
#            Subarctic:    Stations 7 and 11
#            Transition:   Stations 15 and 21
#            Subtropical:  Stations 27 and 34
#
#   B) OCEAN SURFACE CURRENTS:
#        - Major North Pacific surface currents are obtained from the
#          ESRI Major Ocean Surface Current FeatureServer.
#        - Current direction is represented using arrows.
#        - Major currents are labeled directly on the map.
#
#   C) OCEANOGRAPHIC REGIONS:
#        - Eastern Pacific Subarctic Gyre
#        - North Pacific Polar Front
#        - North Pacific Subtropical Gyre
#        - Labels are manually positioned to provide geographic context
#          for the sampling transect.
#
# MAP EXTENT:
#   - Longitude: 180°W to 105°W
#   - Latitude:    5°N to 65°N
#
# MAP LAYERS:
#   1) Land boundaries
#   2) Major North Pacific surface currents and current direction
#   3) Oceanographic region labels
#   4) NH1208 MOCNESS sampling stations colored by region
#   5) Longitude and latitude coordinates
#
# OUTPUT:
#   Saves a high-resolution PNG map to:
#     plots/map_ESRI_currents_NH1208_regions_oceanographic_labels.png
#
# NOTE:
#   Oceanographic region labels are manually positioned for visualization
#   and represent broad geographic features rather than formal boundaries.
# -----------------------------------------------------------------------------

# =============================================================================
# USER INPUT
# =============================================================================


# -----------------------------------------------------------------------------
# Workflow notes:
# - Extracts NH1208 station coordinates and region assignments from the processed Chapter 1 dataset.
# - Adds major North Pacific surface-current geometry and land boundaries to provide oceanographic context for the transect.
# - Uses manually positioned broad-region/current labels for visualization; these labels are not analytical boundaries.
# - Saves the publication map showing sampling stations, oceanographic regions, and major circulation features.
# -----------------------------------------------------------------------------

source("R_Scripts/00_setup.R", local = TRUE)
DATA_CSV <- file.path(OUT_DIR, "data_gradients.csv")
PLOT_DIR <- PLOT_DIRS[["Map"]]
OUT_MAP <- file.path(
  PLOT_DIR,
  "map_ESRI_currents_NH1208_regions_oceanographic_labels.png"
)

if (!dir.exists(dirname(OUT_MAP))) {
  dir.create(dirname(OUT_MAP), recursive = TRUE)
}


# =============================================================================
# 0) PACKAGES
# =============================================================================

required_pkgs <- c(
  "sf",
  "ggplot2",
  "dplyr",
  "readr",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggrepel"
)

new_pkgs <- setdiff(
  required_pkgs,
  rownames(installed.packages())
)

if (length(new_pkgs) > 0) {
  install.packages(
    new_pkgs,
    repos = "https://cloud.r-project.org"
  )
}

invisible(
  lapply(
    required_pkgs,
    function(p) library(p, character.only = TRUE)
  )
)


# =============================================================================
# 1) READ NH1208 MOCNESS STATION DATA
# =============================================================================

dat <- readr::read_csv(
  DATA_CSV,
  show_col_types = FALSE
)

stations_tbl <- dat |>
  
  dplyr::mutate(
    Cruise = tolower(Cruise),
    Station = as.numeric(Station),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  ) |>
  
  dplyr::filter(
    Cruise == "nh1208",
    !is.na(Station),
    is.finite(lon),
    is.finite(lat)
  ) |>
  
  dplyr::mutate(
    Region = dplyr::case_when(
      Station %in% c(7, 11)  ~ "Subarctic",
      Station %in% c(15, 21) ~ "Transition",
      Station %in% c(27, 34) ~ "Subtropical",
      TRUE ~ NA_character_
    )
  ) |>
  
  dplyr::filter(
    !is.na(Region)
  ) |>
  
  dplyr::distinct(
    Station,
    Region,
    lon,
    lat
  ) |>
  
  dplyr::mutate(
    Region = factor(
      Region,
      levels = c(
        "Subarctic",
        "Transition",
        "Subtropical"
      )
    )
  )

message("NH1208 stations found:")
print(stations_tbl)

if (nrow(stations_tbl) == 0) {
  stop(
    "No NH1208 stations found. ",
    "Check Cruise, Station, lon, and lat columns."
  )
}


# Convert station table to sf object

stations_sf <- sf::st_as_sf(
  stations_tbl,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)


# =============================================================================
# 2) READ ESRI MAJOR OCEAN SURFACE CURRENT LAYER
# =============================================================================

currents_url <- paste0(
  "https://services4.arcgis.com/InH2ctI0mXwZbjMf/arcgis/rest/services/",
  "Major%20Ocean%20Surface%20Current/FeatureServer/0/query?",
  "where=Ocean_Basin%3D%27NP%27",
  "&outFields=*",
  "&returnGeometry=true",
  "&outSR=4326",
  "&f=geojson"
)

currents <- sf::st_read(
  currents_url,
  quiet = FALSE
)

if (nrow(currents) == 0) {
  stop(
    "No current features were returned from the ESRI service."
  )
}

message(
  "Number of current features read: ",
  nrow(currents)
)

message("Current layer columns:")
print(names(currents))


# =============================================================================
# 3) READ WORLD LAND POLYGONS
# =============================================================================

world <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
)

currents <- sf::st_transform(
  currents,
  4326
)

world <- sf::st_transform(
  world,
  4326
)

stations_sf <- sf::st_transform(
  stations_sf,
  4326
)


# =============================================================================
# 4) DEFINE NORTHEAST PACIFIC MAP EXTENT
# =============================================================================

lon_min <- -180
lon_max <- -105

lat_min <- 5
lat_max <- 65


bbox_nep <- sf::st_bbox(
  c(
    xmin = lon_min,
    xmax = lon_max,
    ymin = lat_min,
    ymax = lat_max
  ),
  crs = sf::st_crs(4326)
)


# Crop currents

currents_nep <- suppressWarnings(
  sf::st_crop(
    currents,
    bbox_nep
  )
)


# Crop land

world_nep <- suppressWarnings(
  sf::st_crop(
    world,
    bbox_nep
  )
)


# Crop stations

stations_nep <- suppressWarnings(
  sf::st_crop(
    stations_sf,
    bbox_nep
  )
)


if (nrow(currents_nep) == 0) {
  stop(
    "No current features fall inside the Northeast Pacific bounding box."
  )
}

if (nrow(stations_nep) == 0) {
  stop(
    "No NH1208 stations fall inside the Northeast Pacific bounding box."
  )
}


message("Currents inside map:")
print(
  unique(currents_nep$Name)
)

message("Stations inside map:")
print(stations_nep)


# =============================================================================
# 5) CREATE CURRENT LABEL POSITIONS
# =============================================================================

currents_for_labels_proj <- currents_nep |>
  
  dplyr::filter(
    !is.na(Name)
  ) |>
  
  dplyr::mutate(
    Current_full_name = as.character(Name)
  ) |>
  
  sf::st_transform(
    3857
  )


# Get representative point along each current

label_pts_proj <- sf::st_point_on_surface(
  currents_for_labels_proj
)

label_pts <- sf::st_transform(
  label_pts_proj,
  4326
)


# Build table containing label coordinates

label_tbl <- currents_for_labels_proj |>
  
  sf::st_transform(
    4326
  ) |>
  
  sf::st_drop_geometry() |>
  
  dplyr::select(
    Current_full_name,
    Temperature
  ) |>
  
  dplyr::bind_cols(
    
    sf::st_coordinates(
      label_pts
    ) |>
      
      as.data.frame() |>
      
      dplyr::rename(
        lon = X,
        lat = Y
      )
  )


# Average positions if current occurs in multiple segments

label_tbl <- label_tbl |>
  
  dplyr::group_by(
    Current_full_name,
    Temperature
  ) |>
  
  dplyr::summarise(
    lon = mean(
      lon,
      na.rm = TRUE
    ),
    
    lat = mean(
      lat,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# =============================================================================
# 6) OCEANOGRAPHIC REGION LABELS
# =============================================================================
#
# Coordinates can be adjusted manually here.
#
# IMPORTANT:
# Eastern Pacific Subarctic Gyre has been moved NORTH:
#
#       OLD = 52 degrees N
#       NEW = 57 degrees N
#
# =============================================================================

region_labels <- data.frame(
  
  label = c(
    "Eastern Pacific\nSubarctic Gyre",
    "North Pacific\nPolar Front",
    "North Pacific\nSubtropical Gyre"
  ),
  
  lon = c(
    -148,
    -145,
    -138
  ),
  
  lat = c(
    57,     # moved higher
    42,
    24
  )
)


# =============================================================================
# 7) LONGITUDE AND LATITUDE AXIS BREAKS
# =============================================================================

lon_breaks <- seq(
  -180,
  -110,
  by = 10
)

lat_breaks <- seq(
  10,
  60,
  by = 10
)


# Longitude labels

lon_labels <- function(x) {
  paste0(
    abs(x),
    "°W"
  )
}


# Latitude labels

lat_labels <- function(x) {
  paste0(
    x,
    "°N"
  )
}


# =============================================================================
# 8) PLOT
# =============================================================================

p <- ggplot2::ggplot() +
  
  
  # ---------------------------------------------------------------------------
# LAND
# ---------------------------------------------------------------------------

ggplot2::geom_sf(
  data = world_nep,
  fill = "black",
  colour = "black",
  linewidth = 0.2
) +
  
  
  # ---------------------------------------------------------------------------
# ESRI OCEAN CURRENTS
# ---------------------------------------------------------------------------

ggplot2::geom_sf(
  data = currents_nep,
  
  ggplot2::aes(
    colour = Temperature
  ),
  
  linewidth = 1.15,
  alpha = 0.75,
  
  arrow = grid::arrow(
    length = grid::unit(
      0.20,
      "cm"
    ),
    type = "closed"
  ),
  
  lineend = "round"
) +
  
  
  # ---------------------------------------------------------------------------
# CURRENT NAME LABELS - WHITE HALO
# ---------------------------------------------------------------------------

ggrepel::geom_text_repel(
  data = label_tbl,
  
  ggplot2::aes(
    x = lon,
    y = lat,
    label = Current_full_name
  ),
  
  colour = "white",
  
  # Increased from 4.0
  size = 5.0,
  
  fontface = "bold",
  segment.colour = NA,
  
  box.padding = 0.50,
  point.padding = 0.30,
  
  max.overlaps = Inf
) +
  
  
  # ---------------------------------------------------------------------------
# CURRENT NAME LABELS - BLACK TEXT
# ---------------------------------------------------------------------------

ggrepel::geom_text_repel(
  data = label_tbl,
  
  ggplot2::aes(
    x = lon,
    y = lat,
    label = Current_full_name
  ),
  
  colour = "black",
  
  # Increased from 3.3
  size = 4.2,
  
  fontface = "bold",
  segment.colour = NA,
  
  box.padding = 0.50,
  point.padding = 0.30,
  
  max.overlaps = Inf
) +
  
  
  # ---------------------------------------------------------------------------
# OCEANOGRAPHIC REGION LABELS
# ---------------------------------------------------------------------------

ggplot2::geom_label(
  data = region_labels,
  
  ggplot2::aes(
    x = lon,
    y = lat,
    label = label
  ),
  
  fill = "white",
  colour = "black",
  
  label.size = 0.30,
  
  # Increased from 3.4
  size = 4.5,
  
  fontface = "bold",
  alpha = 0.88,
  lineheight = 0.95
) +
  
  
  # ---------------------------------------------------------------------------
# NH1208 MOCNESS STATIONS
# ---------------------------------------------------------------------------

ggplot2::geom_sf(
  data = stations_nep,
  
  ggplot2::aes(
    fill = Region
  ),
  
  shape = 21,
  size = 6.0,
  colour = "white",
  stroke = 1.3
) +
  
  
  # ---------------------------------------------------------------------------
# MAP EXTENT
# ---------------------------------------------------------------------------

ggplot2::coord_sf(
  xlim = c(
    lon_min,
    lon_max
  ),
  
  ylim = c(
    lat_min,
    lat_max
  ),
  
  expand = FALSE
) +
  
  
  # ---------------------------------------------------------------------------
# LONGITUDE AXIS
# ---------------------------------------------------------------------------

ggplot2::scale_x_continuous(
  name = "Longitude (°)",
  breaks = lon_breaks,
  labels = lon_labels,
  expand = c(0, 0)
) +
  
  
  # ---------------------------------------------------------------------------
# LATITUDE AXIS
# ---------------------------------------------------------------------------

ggplot2::scale_y_continuous(
  name = "Latitude (°)",
  breaks = lat_breaks,
  labels = lat_labels,
  expand = c(0, 0)
) +
  
  
  # ---------------------------------------------------------------------------
# CURRENT COLORS
# ---------------------------------------------------------------------------

ggplot2::scale_colour_manual(
  
  values = c(
    "Warm" = "red",
    "Cold" = "blue"
  ),
  
  na.value = "grey45",
  
  name = "Current type",
  
  guide = "none"
) +
  
  
  # ---------------------------------------------------------------------------
# REGION COLORS
# ---------------------------------------------------------------------------

ggplot2::scale_fill_manual(
  
  values = REGION_COLORS,
  
  breaks = c(
    "Subarctic",
    "Transition",
    "Subtropical"
  ),
  
  name = "Region"
) +
  
  
  # ---------------------------------------------------------------------------
# THEME
# ---------------------------------------------------------------------------

ggplot2::theme_classic() +
  
  ggplot2::theme(
    
    
    # Background --------------------------------------------------------------
    
    plot.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    
    panel.background = ggplot2::element_rect(
      fill = "white",
      colour = NA
    ),
    
    
    # Border ------------------------------------------------------------------
    
    panel.border = ggplot2::element_rect(
      fill = NA,
      colour = "black",
      linewidth = 0.6
    ),
    
    
    # Axis titles -------------------------------------------------------------
    
    axis.title = ggplot2::element_text(
      colour = "black",
      size = 16
    ),
    
    
    # Coordinate numbers ------------------------------------------------------
    
    axis.text = ggplot2::element_text(
      colour = "black",
      size = 14
    ),
    
    
    # Axis ticks --------------------------------------------------------------
    
    axis.ticks = ggplot2::element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    
    # Legend ------------------------------------------------------------------
    
    legend.position = "left",
    
    legend.title = ggplot2::element_text(
      face = "bold",
      colour = "black",
      size = 15
    ),
    
    legend.text = ggplot2::element_text(
      colour = "black",
      size = 13
    ),
    
    
    # Margins -----------------------------------------------------------------
    
    plot.margin = ggplot2::margin(
      t = 12,
      r = 12,
      b = 12,
      l = 12
    )
  )


# =============================================================================
# 9) SAVE MAP
# =============================================================================

save_plot_formats(
  OUT_MAP,
  p,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)

message(
  "Wrote map: ",
  OUT_MAP
)
