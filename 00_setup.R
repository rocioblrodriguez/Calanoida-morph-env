#!/usr/bin/env Rscript
# =============================================================================
# Script:  00_setup.R
# Author:  Rocio Rodriguez
# Purpose: Define the paths, analysis settings, variables, and figure style
#          shared by the CHP1 scripts. Run from the revised analysis folder.
# =============================================================================

# -----------------------------------------------------------------------------
# Workflow notes:
# - Centralizes project paths, analysis constants, descriptor groups, cluster/region ordering, colors, and icon locations.
# - All downstream scripts source this file so shared definitions remain consistent across analyses and figures.
# - Morphological descriptor categories defined here are descriptive groupings used throughout the workflow.
# - This script defines settings only; it does not fit models or modify the input observations.
# -----------------------------------------------------------------------------


# 1. Project paths ------------------------------------------------------------
# Set the working directory to:
# /Users/rociorodriguez/Git_work/CHP1/Final_R/CHP1_revised_scripts
# Read the existing inputs from the shared CHP1 data directory below. The TSV
# and JPGs come from the same October 2026 export. Preserve the CTD cruise subfolders
# and original JPG and icon filenames. New results stay in the revised folder.
PROJECT_DIR <- normalizePath(".", winslash = "/")
DATA_DIR <- "/Users/rociorodriguez/Git_work/CHP1/data"
EXPORT_DIR <- file.path(DATA_DIR, "export_TSV_5421_job260837")
OUT_DIR <- file.path(PROJECT_DIR, "outputs")
PLOT_DIR <- file.path(PROJECT_DIR, "plots")
CTD_DIR <- file.path(DATA_DIR, "CTD")
IMAGES_DIR <- file.path(EXPORT_DIR, "export__TSV_5421_20261031_2221")
ICON_DIR <- file.path(DATA_DIR, "1x")
ECOTAXA_TSV <- file.path(EXPORT_DIR, "ecotaxa_export__TSV_5421_20261031_2221.tsv")
DIAGNOSTICS_CSV <- file.path(DATA_DIR, "morphology_diagnostics_long.csv")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Keep all figures for each section together. No existing files are moved.
PLOT_DIRS <- list(
  PCA_statistics = file.path(PLOT_DIR, "PCA", "statistics"),
  PCA_6clusters = file.path(PLOT_DIR, "PCA", "6clusters"),
  PCA_2clusters = file.path(PLOT_DIR, "PCA", "2clusters"),
  Size_transparency = file.path(PLOT_DIR, "Size_transparency"),
  Behavior = file.path(PLOT_DIR, "Behavior"),
  RDA = file.path(PLOT_DIR, "RDA"),
  Centroid_environment_links = file.path(PLOT_DIR, "Centroid_environment_links"),
  Butterfly = file.path(PLOT_DIR, "Butterfly"),
  Silhouette = file.path(PLOT_DIR, "Silhouette"),
  Map = file.path(PLOT_DIR, "Map"),
  Environmental_sections = file.path(PLOT_DIR, "Environmental_sections")
)
invisible(lapply(PLOT_DIRS, dir.create, recursive = TRUE, showWarnings = FALSE))

# 2. Analysis settings --------------------------------------------------------
# Processing retains both cruises. The morphology analyses use NH1208 only.
CRUISE <- "nh1208"
CATEGORY <- "Calanoida"
STATUS <- "validated"
SEED <- 123
N_PERM <- 999
KMEANS_NSTART <- 100
REGION_LEVELS <- c("Subarctic", "Transition", "Subtropical")
CLUSTER_LEVELS <- c("ST", "SD", "MT", "MD", "LT", "LD")
NODE_CLUSTER_LEVELS <- c("ST", "MT", "LT", "SD", "MD", "LD")
ENV_RAW <- c("Temp", "Sal", "O2", "Fluor")
ENV_STD <- paste0(ENV_RAW, "_std")
DIEL_RUNS <- list(Day = "D", Night = "N", DayNight = c("D", "N"))

# 3. Morphological descriptors ------------------------------------------------
# These four groups are the single category definition used throughout.
# In particular, fcons is Complexity; histogram/intensity descriptors are
# Transparency. This classification affects labels and category summaries.
MORPH_GROUPS <- list(
  Transparency = c(
    "mean", "stddev", "mode", "min", "max", "intden", "median", "skew",
    "kurt", "cv", "sr", "histcum1", "histcum2", "histcum3", "nb1", "nb2",
    "meanpos"
  ),
  Size = c("area", "major", "minor", "feret", "area_exc", "skelarea", "convarea", "esd"),
  Shape = c(
    "perim.", "circ.", "fractal", "symetrieh", "symetriev", "symetriehc",
    "symetrievc", "convperim", "thickr", "elongation", "circex", "cdexc"
  ),
  Complexity = c("slope", "fcons", "perimareaexc", "feretareaexc", "perimferet", "perimmajor")
)


# Preserve the original 43-column analysis order while using the shared groups
# above to assign categories. The suffix identifies transformed descriptors.
MORPH_VARIABLES <- paste0(c(
  "area", "major", "minor", "feret", "area_exc", "skelarea", "convarea", "esd",
  "mean", "stddev", "mode", "min", "max", "intden", "median", "skew", "kurt",
  "cv", "sr", "perim.", "circ.", "fractal", "symetriehc", "symetrievc",
  "symetrieh", "symetriev", "convperim", "fcons", "thickr", "elongation",
  "circex", "cdexc", "slope", "histcum1", "histcum2", "histcum3", "nb1",
  "nb2", "meanpos", "perimareaexc", "feretareaexc", "perimferet", "perimmajor"
), "_YJ")

# Transform the same full descriptor set used by the original correlation
# script. Explicit names replace the previous positional selection (24:90).
ALL_MORPH_RAW <- c(
  "area", "mean", "stddev", "mode", "min", "max", "x", "y", "xm", "ym",
  "perim.", "bx", "by", "width", "height", "major", "minor", "angle", "circ.",
  "feret", "intden", "median", "skew", "kurt", "%area", "xstart", "ystart",
  "area_exc", "fractal", "skelarea", "slope", "histcum1", "histcum2", "histcum3",
  "xmg5", "ymg5", "nb1", "nb2", "nb3", "compentropy", "compmean", "compslope",
  "compm1", "compm2", "compm3", "symetrieh", "symetriev", "symetriehc",
  "symetrievc", "convperim", "convarea", "fcons", "thickr", "tag", "esd",
  "elongation", "range", "meanpos", "centroids", "cv", "sr", "perimareaexc",
  "feretareaexc", "perimferet", "perimmajor", "circex", "cdexc"
)

# Final six-cluster solution: fixed 20-descriptor candidate selected from the
# silhouette analysis. The plotting scripts never select variables.
CLUSTERING_GROUPS <- list(
  Size = c(
    "area_exc", "esd", "convarea", "feret", "skelarea"
  ),
  Transparency = c(
    "histcum1", "meanpos", "intden", "skew", "nb2"
  ),
  Shape = c(
    "perim.", "cdexc", "circex", "symetriev", "convperim"
  ),
  Complexity = c(
    "slope", "perimferet", "fractal", "fcons", "feretareaexc"
  )
)

CLUSTERING_VARIABLES <- paste0(
  unlist(CLUSTERING_GROUPS, use.names = FALSE),
  "_YJ"
)
# Net intervals are shared by processing, behavior, and environmental sections.
# Behavior and butterfly figures retain their original net-depth codes.
# Midpoints are used by the environmental sections.
NET_BINS <- data.frame(
  net = c("n8", "n7", "n6", "n5", "n4", "n3", "n2", "n1"),
  depth_min = c(0, 25, 50, 100, 200, 400, 600, 800),
  depth_max = c(25, 50, 100, 200, 400, 600, 800, 1000)
)
NET_BINS$depth_mid <- (NET_BINS$depth_min + NET_BINS$depth_max) / 2

# 4. Figure style -------------------------------------------------------------
# Use these named mappings in every plot, including combined figures.
CLUSTER_COLORS <- c(
  ST = "#FED976", SD = "#FC4E2A", MT = "#FEB24C",
  MD = "#BD0026", LT = "#FD8D3C", LD = "#800026"
)
CLUSTER_LABELS <- c(
  ST = "ST = Small Transparent", SD = "SD = Small Dark",
  MT = "MT = Medium Transparent", MD = "MD = Medium Dark",
  LT = "LT = Large Transparent", LD = "LD = Large Dark"
)
REGION_COLORS <- c(Subarctic = "#0000FF", Transition = "#FFA500", Subtropical = "#FF0000")
REGION_SHAPES <- c(Subarctic = 21, Transition = 22, Subtropical = 23)
REGION_TAGS <- c(Subarctic = "A", Transition = "B", Subtropical = "C")
# Preserve the original icon filenames; no asset copies or renaming are needed.
CLUSTER_ICONS <- c(
  ST = file.path(ICON_DIR, "Asset 1.png"),
  SD = file.path(ICON_DIR, "Asset 4.png"),
  MT = file.path(ICON_DIR, "Asset 2.png"),
  MD = file.path(ICON_DIR, "Asset 5.png"),
  LT = file.path(ICON_DIR, "Asset 3.png"),
  LD = file.path(ICON_DIR, "Asset 6.png")
)
MORPH_COLORS <- c(Transparency = "#2C7FB8", Size = "#E66101", Shape = "#238B45", Complexity = "#AE2D9F")
ENV_COLORS <- c(Temp = "#E64B35", Sal = "#FC8D62", O2 = "#4DBBD5", Fluor = "#00A087")
ENV_LABELS <- c(
  Temp = "Temperature (°C)", Sal = "Salinity (PSU)",
  O2 = "Oxygen (µmol kg⁻¹)", Fluor = "Fluorescence (mg m⁻³)"
)
TRAIT_PROXIES <- c(
  TRANSPARENT = "mean_YJ", LARGE = "major_YJ", DARK = "stddev_YJ",
  CIRCULAR = "circ._YJ", SYMMETRICAL = "symetrieh_YJ"
)
FIGURE_DPI <- 300
FIGURE_FORMATS <- c("png", "pdf", "svg", "eps")
# Font sizes, legend positions, and export sizes remain in each plotting script.

source(file.path(PROJECT_DIR, "R_Scripts", "helpers.R"), local = TRUE)
