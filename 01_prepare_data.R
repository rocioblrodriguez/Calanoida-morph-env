#!/usr/bin/env Rscript
# Script:  01_prepare_data.R
# -----------------------------------------------------------------------------
# Title:     Preprocess EcoTaxa with CTD-derived Environmental Data (simple)
# Author:    Rocio Rodriguez
# Date:      2025-10-24
# Purpose:   1) Load EcoTaxa detections
#            2) Parse IDs and compute density/abundance
#            3) Add Day/Night and Station/Region metadata
#            4) Read Sea-Bird CTD casts & pick nearest cast to each Station
#            5) Derive MOC-like binned CTD environmental variables
#            6) Merge CTD env back into EcoTaxa table and write CSV
# -----------------------------------------------------------------------------

# ---- USER INPUT ----
# Paths come from 00_setup.R; the original processing logic below is retained.

# -----------------------------------------------------------------------------
# Workflow notes:
# - Reads the raw EcoTaxa export and parses object identifiers into cruise, MOCNESS tow, net, fraction, and scan fields.
# - Applies the project-specific sampled-volume corrections before calculating organism density and vertically integrated abundance.
# - Assigns day/night status, station, and oceanographic region metadata to each observation.
# - Reads CTD casts, associates stations with the appropriate environmental profiles, and summarizes CTD variables over the MOCNESS depth intervals.
# - Merges the resulting environmental information back onto the organism-level EcoTaxa table used by all later scripts.
# -----------------------------------------------------------------------------

source("R_Scripts/00_setup.R", local = TRUE)
OUT_CSV <- file.path(OUT_DIR, "data_gradients.csv")
# ---- packages ----
required_pkgs <- c("readr","dplyr","tidyr","lubridate","stringr","purrr","sf","oce")
new_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(new_pkgs) > 0) install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(required_pkgs, function(p) library(p, character.only = TRUE)))

# 1) Load EcoTaxa export -------------------------------------------------------
data_raw <- readr::read_tsv(ECOTAXA_TSV, col_names = TRUE, show_col_types = FALSE)

# 2) Process detections --------------------------------------------------------
data <- data_raw |>
  tidyr::separate(
    col    = object_id,
    into   = c("cruise","moc","net","fraction","scan"),
    extra  = "drop",
    remove = FALSE
  ) |>
  dplyr::mutate(
    sample_tot_vol = dplyr::case_when(
      cruise == "nh1208" & moc == "m7"  & net %in% c( "n5") ~ 772.4,
      cruise == "nh1208" & moc == "m19" & net %in% c( "n6") ~ 282.8,
      cruise == "nh1208" & moc == "m19" & net %in% c("n7") ~ 278.7,
      cruise == "nh1208" & moc == "m19" & net %in% c("n8") ~ 436.9,
      cruise == "oc473"  & moc == "m4"  & net %in% c("n3") ~ 1186.2,
      TRUE ~ sample_tot_vol
    ),
    density   = acq_sub_part / sample_tot_vol,
    abundance = density * (object_depth_max - object_depth_min),
    process_particle_pixel_size_mm = dplyr::if_else(
      is.nan(process_particle_pixel_size_mm), 0.0053, process_particle_pixel_size_mm
    ),
    object_lon = dplyr::if_else(
      cruise == "nh1208" & moc == "m9", -141.4467, object_lon
    ),
    Cruise = cruise
  ) |>
  dplyr::rename_with(~ gsub("^object_", "", .x))

# 3) Day/Night flag from time --------------------------------------------------
data <- data |>
  dplyr::mutate(
    time_period = suppressWarnings(lubridate::hms(time)),
    D_N = dplyr::case_when(
      !is.na(time_period) & lubridate::hour(time_period) >= 6 & lubridate::hour(time_period) < 18 ~ "D",
      TRUE ~ "N"
    )
  ) |>
  dplyr::select(-time_period)

# 3b) Map (Cruise, moc) -> Station --------------------------------------------
data <- data |>
  dplyr::mutate(
    Station = dplyr::case_when(
      # Oceanus = oc473
      Cruise == "oc473"  & moc %in% c("m4","m5")    ~ 5,
      Cruise == "oc473"  & moc %in% c("m6","m7")    ~ 8,
      Cruise == "oc473"  & moc %in% c("m8","m9")    ~ 13,
      Cruise == "oc473"  & moc %in% c("m12","m13")  ~ 21,
      Cruise == "oc473"  & moc %in% c("m14","m15")  ~ 26,
      Cruise == "oc473"  & moc %in% c("m16","m17")  ~ 31,
      # New Horizon = nh1208
      Cruise == "nh1208" & moc %in% c("m7","m8")    ~ 7,
      Cruise == "nh1208" & moc %in% c("m9","m10")   ~ 11,
      Cruise == "nh1208" & moc %in% c("m11","m12")  ~ 15,
      Cruise == "nh1208" & moc %in% c("m15","m16")  ~ 21,
      Cruise == "nh1208" & moc %in% c("m19","m20")  ~ 27,
      Cruise == "nh1208" & moc %in% c("m25","m26")  ~ 34,
      TRUE ~ NA_real_
    )
  )

# 4) Region assignment ---------------------------------------------------------
data <- data |>
  dplyr::mutate(
    Region = dplyr::case_when(
      Cruise == "nh1208" & Station %in% c(7, 11)  ~ "Subarctic",
      Cruise == "nh1208" & Station %in% c(15, 21) ~ "Transition",
      Cruise == "nh1208" & Station %in% c(27, 34) ~ "Subtropical",
      Cruise == "oc473"  & Station %in% c(5, 8)   ~ "Subtropical",
      Cruise == "oc473"  & Station %in% c(13, 21) ~ "Transition",
      Cruise == "oc473"  & Station %in% c(26, 31) ~ "Subarctic",
      TRUE ~ NA_character_
    )
  )

# 5) Read CTD & pick nearest cast per station ---------------------------------
files <- list.files(CTD_DIR, pattern = "\\.[Cc][Nn][Vv]$", full.names = TRUE, recursive = TRUE)
read_one <- function(path) {
  ctd <- oce::read.ctd.sbe(path) |> oce::ctdTrim(method = "downcast")
  cruise_id <- basename(dirname(path))
  as.data.frame(ctd[["data"]]) |>
    dplyr::mutate(
      path      = path,
      file      = basename(path),
      cast_id   = stringr::str_remove(basename(path), "\\.[Cc][Nn][Vv]$"),
      cruise    = cruise_id,
      longitude = ctd@metadata$longitude,
      latitude  = ctd@metadata$latitude,
      station   = if (!is.null(ctd@metadata$station)) ctd@metadata$station else NA_character_,
      startTime = ctd@metadata$startTime
    )
}
dfs <- purrr::map(files, ~ tryCatch(read_one(.x),
                                    error = function(e) { warning(sprintf("Skipping %s: %s", .x, e$message)); NULL }))
all_ctd <- dplyr::bind_rows(Filter(Negate(is.null), dfs))

casts <- all_ctd |>
  dplyr::distinct(cast_id, cruise, longitude, latitude) |>
  dplyr::filter(is.finite(longitude), is.finite(latitude))

stations_tbl <- data |>
  dplyr::mutate(lat = as.numeric(lat), lon = as.numeric(lon)) |>
  dplyr::filter(is.finite(lat), is.finite(lon)) |>
  dplyr::distinct(Cruise, Station, lon, lat)

stations_sf <- sf::st_as_sf(stations_tbl, coords = c("lon","lat"), crs = 4326, remove = FALSE)
casts_sf    <- sf::st_as_sf(casts,       coords = c("longitude","latitude"), crs = 4326, remove = FALSE)

nearest_tbl <- unique(stations_sf$Cruise) |>
  purrr::map_dfr(function(cru){
    stas <- dplyr::filter(stations_sf, Cruise == cru)
    cts  <- dplyr::filter(casts_sf, cruise == cru)
    if (nrow(cts) == 0 || nrow(stas) == 0) return(dplyr::tibble())
    idx  <- sf::st_nearest_feature(stas, cts)
    dist <- as.numeric(sf::st_distance(stas, cts[idx,], by_element = TRUE)) / 1000
    sf::st_drop_geometry(stas) |>
      dplyr::bind_cols(
        sf::st_drop_geometry(cts[idx, c("cast_id","longitude","latitude")]) |>
          dplyr::rename(cast_lon = longitude, cast_lat = latitude)
      ) |>
      dplyr::mutate(distance_km = dist)
  })

nearest_ranked <- nearest_tbl |>
  dplyr::group_by(Cruise, Station) |>
  dplyr::arrange(distance_km, .by_group = TRUE) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::ungroup()

nearest_one <- nearest_ranked |>
  dplyr::filter(dplyr::case_when(Station == 11 ~ rank == 2L, TRUE ~ rank == 1L)) |>
  dplyr::select(-rank)

sel_casts <- nearest_one |>
  dplyr::filter(!is.na(cast_id)) |>
  dplyr::pull(cast_id) |>
  unique()

meta_cast_station <- nearest_one |>
  dplyr::filter(!is.na(cast_id)) |>
  dplyr::distinct(Cruise, cast_id, Station)

# 6) Derive CTD env per Station × net (MOC-like bins) -------------------------
vars <- intersect(c("temperature","salinity","oxygen","fluorescence"), names(all_ctd))

profiles_df <- all_ctd |>
  dplyr::filter(cast_id %in% sel_casts, depth >= 0, depth <= 1000) |>
  dplyr::select(cruise, cast_id, depth, dplyr::all_of(vars)) |>
  tidyr::pivot_longer(cols = dplyr::all_of(vars), names_to = "variable", values_to = "value") |>
  dplyr::inner_join(meta_cast_station, by = c("cast_id", "cruise" = "Cruise"))

net_breaks <- c(0, 25, 50, 100, 200, 400, 600, 800, 1000)
net_labels <- c("n8","n7","n6","n5","n4","n3","n2","n1")

profiles_net_means <- profiles_df |>
  dplyr::mutate(
    net = cut(depth, breaks = net_breaks, labels = net_labels, include.lowest = TRUE, right = TRUE),
    net = factor(net, levels = net_labels)
  ) |>
  dplyr::group_by(cruise, Station, net, variable) |>
  dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

env_ctd <- profiles_net_means |>
  tidyr::pivot_wider(names_from = variable, values_from = value) |>
  dplyr::rename(
    Cruise = cruise,
    Temp   = temperature,
    Sal    = salinity,
    O2     = oxygen,
    Fluor  = fluorescence
  )

# 7) Merge CTD env into EcoTaxa table -----------------------------------------
data_merged <- data |>
  dplyr::left_join(env_ctd, by = c("Cruise","Station","net"))

# 8) Write out -----------------------------------------------------------------
readr::write_csv(data_merged, OUT_CSV)
message("Wrote: ", OUT_CSV)


