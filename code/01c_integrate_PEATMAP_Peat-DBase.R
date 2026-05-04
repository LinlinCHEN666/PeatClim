# ============================================================
# Script: 01c_integrate_PEATMAP_Peat-DBase.R
# Purpose:
#   Integrate prepared PEATMAP and Peat-DBase occurrence datasets
#   into a single modern peatland occurrence dataset on the
#   WorldClim 5 arc-min grid for downstream SDM use.
#
# Inputs:
#   - Prepared PEATMAP occurrence CSV from 01a_prepare_PEATMAP_occurrence.R
#   - Prepared Peat-DBase occurrence CSV from 01b_prepare_Peat-DBase_occurrence.R
#   - WorldClim v2.1 5 arc-min bioclimatic raster
#     (bio1 used only as grid template)
#
# Outputs:
#   - Integrated PEATMAP and Peat-DBase occurrence CSV
#   - CSV summary of counts by processing step
#
# Notes:
#   - Both input datasets are expected to contain columns:
#       lon, lat, fid, VALUE
#   - The final output keeps one record per 5 arc-min grid cell.
#   - Duplicate cells are removed after combining the two datasets.
#   - New fid values are assigned after integration.
#
# External data requirements:
#   PEATMAP dataset
#   - Citation:
#       Xu, J., Morris, P. J., Liu, J., and Holden, J. (2017).
#       PEATMAP: Refining estimates of global peatland distribution
#       based on a meta-analysis. University of Leeds. [Dataset]
#       https://doi.org/10.5518/252
#   - Download from:
#       https://archive.researchdata.leeds.ac.uk/251/
#
#   Peat-DBase dataset
#   - Citation:
#       Skye, J. et al. (2025).
#       Peat-DBase v.1: A Compiled Database of Global Peat Depth Measurements
#       (1.0.0b). Zenodo.
#   - Download from:
#       https://doi.org/10.5281/zenodo.15530644
#
#   WorldClim v2.1 bioclimatic variables
#   - Citation:
#       Fick, S. E. and Hijmans, R. J. (2017).
#       WorldClim 2: new 1 km spatial resolution climate surfaces
#       for global land areas. International Journal of Climatology,
#       37(12), 4302-4315.
#   - Data used here:
#       WorldClim v2.1 bioclimatic variables at 5 arc-min resolution
#       (bio1 used only as grid template in this script)
#   - Download from:
#       https://www.worldclim.org/data/worldclim21.html
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

# ---------------------------- #
# User settings
# ---------------------------- #

# Rasterization rule:
# FALSE = centre-of-cell rule
# TRUE  = any touched cell counts as presence
touches_rule <- FALSE

# ---------------------------- #
# Helper function
# ---------------------------- #

create_directory <- function(path) {
  if (length(path) != 1) {
    stop("Error: 'path' must be a single string. Received: ", toString(path))
  }
  
  if (!dir.exists(path)) {
    tryCatch({
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      if (!dir.exists(path)) {
        stop("Failed to create directory: ", path)
      }
    }, error = function(e) {
      stop("Error while creating directory ", path, ": ", e$message)
    })
  }
}

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_occurrence <- here("data_processed", "occurrence")
dir_wc         <- here("data_external", "WorldClim", "wc2.1_5m_bio")

create_directory(dir_occurrence)

# Final outputs
if (touches_rule == 'FALSE'){
  out_global_csv  <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min.csv")
  out_summary_csv <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min_summary.csv")
} else if (touches_rule == 'TRUE'){
  out_global_csv  <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min_alt.csv")
  out_summary_csv <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min_summary_alt.csv")
}

# ---------------------------- #
# Input files
# ---------------------------- #
if (touches_rule == "FALSE") {
  path_peatmap <- file.path(dir_occurrence, "PEATMAP_global_5min.csv")
} else if (touches_rule == 'TRUE'){
  path_peatmap <- file.path(dir_occurrence, "PEATMAP_global_5min_alt.csv")
}

path_peatdb  <- file.path(dir_occurrence, "Peat-DBase_global_5min.csv")
template_tif <- file.path(dir_wc, "wc2.1_5m_bio_1.tif")

if (!file.exists(path_peatmap)) {
  stop(
    "Prepared PEATMAP occurrence file not found.\n\n",
    "Expected location:\n",
    path_peatmap, "\n\n",
    "Run 01a_prepare_PEATMAP_occurrence.R first."
  )
}

if (!file.exists(path_peatdb)) {
  stop(
    "Prepared Peat-DBase occurrence file not found.\n\n",
    "Expected location:\n",
    path_peatdb, "\n\n",
    "Run 01b_prepare_Peat-DBase_occurrence.R first."
  )
}

if (!file.exists(template_tif)) {
  stop(
    "WorldClim bioclimatic data not found.\n\n",
    "This project requires WorldClim version 2.1 bioclimatic variables ",
    "at 5 arc-min resolution (1970-2000 average).\n",
    "Download from: https://www.worldclim.org/data/worldclim21.html\n\n",
    "Expected location:\n",
    here("data_external", "WorldClim", "wc2.1_5m_bio")
  )
}

PEATMAP <- read_csv(path_peatmap, show_col_types = FALSE)
PeatDB  <- read_csv(path_peatdb, show_col_types = FALSE)
template <- rast(template_tif)[[1]]

message("PEATMAP occurrence file loaded: ", path_peatmap)
message("Peat-DBase occurrence file loaded: ", path_peatdb)
message("Template raster loaded: ", template_tif)
message("Resolution: ", paste(res(template), collapse = ", "))
message("CRS: ", crs(template))

if (nlyr(template) != 1) {
  stop("Template raster must have exactly one layer.")
}

if (!isTRUE(all.equal(res(template), c(5/60, 5/60), tolerance = 1e-8))) {
  stop("Template raster does not have 5 arc-min resolution.")
}

required_cols <- c("lon", "lat", "fid", "VALUE")

missing_peatmap <- setdiff(required_cols, names(PEATMAP))
missing_peatdb  <- setdiff(required_cols, names(PeatDB))

if (length(missing_peatmap) > 0) {
  stop(
    "PEATMAP input is missing required columns: ",
    paste(missing_peatmap, collapse = ", ")
  )
}

if (length(missing_peatdb) > 0) {
  stop(
    "Peat-DBase input is missing required columns: ",
    paste(missing_peatdb, collapse = ", ")
  )
}

# ---------------------------- #
# Processing
# ---------------------------- #

n_peatmap <- nrow(PEATMAP)
n_peatdb  <- nrow(PeatDB)

# Add source labels for traceability before combining
PEATMAP <- PEATMAP %>%
  mutate(source = "PEATMAP")

PeatDB <- PeatDB %>%
  mutate(source = "Peat-DBase")

# 1. Combine the two prepared occurrence datasets
PMPD <- bind_rows(PEATMAP, PeatDB)

if (nrow(PMPD) == 0) {
  stop("No records found after combining PEATMAP and Peat-DBase.")
}

# 2. Assign each point to a 5 arc-min WorldClim grid cell
PMPD <- PMPD %>%
  mutate(grid5m_cell = cellFromXY(template, cbind(lon, lat)))

# 3. Remove records outside template coverage
PMPD_valid <- PMPD %>%
  filter(!is.na(grid5m_cell))

if (nrow(PMPD_valid) == 0) {
  stop("No integrated records fall within the WorldClim template extent.")
}

# 4. Keep one record per 5 arc-min grid cell
#    Keep first by source order: PEATMAP rows come first, then Peat-DBase rows.
PMPD_dedup <- PMPD_valid %>%
  group_by(grid5m_cell) %>%
  slice(1) %>%
  ungroup()

if (nrow(PMPD_dedup) == 0) {
  stop("No integrated records remain after grid-cell deduplication.")
}

# 5. Final SDM-style output
global <- PMPD_dedup %>%
  transmute(
    lon = lon,
    lat = lat,
    fid = seq_len(n()),
    VALUE = 1L
  )

# ---------------------------- #
# Write final outputs
# ---------------------------- #

write_csv(global, out_global_csv)
message("Integrated occurrence CSV saved: ", out_global_csv)

processing_summary <- tibble(
  step = c(
    "PEATMAP_INPUT",
    "PEATDB_INPUT",
    "COMBINED",
    "WITH_TEMPLATE_COVERAGE",
    "GRID_CELL_DEDUPLICATED",
    "FINAL_OUTPUT"
  ),
  n_records = c(
    n_peatmap,
    n_peatdb,
    nrow(PMPD),
    nrow(PMPD_valid),
    nrow(PMPD_dedup),
    nrow(global)
  )
)

write_csv(processing_summary, out_summary_csv)
message("Processing summary saved: ", out_summary_csv)
