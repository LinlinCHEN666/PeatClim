# ============================================================
# Script: 01b_prepare_Peat-DBase_occurrence.R
# Purpose:
#   Prepare modern peatland occurrence points from Peat-DBase
#   on the WorldClim 5 arc-min grid for downstream SDM use.
#
# Inputs:
#   - Peat-DBase CSV file
#   - WorldClim v2.1 5 arc-min bioclimatic raster
#     (bio1 used only as grid template)
#
# Outputs:
#   - Global Peat-DBase occurrence CSV
#   - CSV summary of counts by processing step
#
# Notes:
#   - This script keeps only peat-presence records with depth_cm >= 30.
#   - Records flagged with error_found == FALSE are retained.
#   - Records are restricted to sample_duplication_flag in {1, 2, 3, 4, 5}.
#   - The final output keeps one record per 5 arc-min grid cell.
#   - When multiple records fall in the same grid cell, the record with the
#     smallest (location_id, sample_id) is retained.
#
# External data requirements:
#   Peat-DBase dataset
#   - Citation:
#       Skye, J., Melton, J. R., Goldblatt, C., Saumier, L., Gallego-Sala, A.,
#       Garneau, M., Winton, R. S., Bahati, E. B., Benavides, J. C.,
#       Fedorchuk, L., Imani, G., Kagaba, C., Kansiime, F., Lamentowicz, M.,
#       Mbasi, M., Wochal, D., Czerwiński, S., Landowski, J., Landowska, J.,
#       Maire, V., Väliranta, M. M., Cole, L. E. S., Davies, M. A., Sun, J.,
#       and Wang, Y. (2025).
#       Peat-DBase v.1: A Compiled Database of Global Peat Depth Measurements
#       (1.0.0b). Zenodo.
#   - Download from:
#       https://doi.org/10.5281/zenodo.15530644
#   - Associated publication:
#       Skye, J. et al. (2025).
#       Peat-DBase v.1: a compiled database of global peat depth measurements.
#       Earth System Science Data, 17, 7313-7330.
#       https://doi.org/10.5194/essd-17-7313-2025
#   - Place extracted files under:
#       data_external/Peat-DBase/
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
#   - Place extracted files under:
#       data_external/WorldClim/
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
})

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

dir_peatdb    <- here("data_external", "Peat-DBase")
dir_wc        <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_processed <- here("data_processed", "occurrence")

create_directory(dir_processed)

# Final outputs
out_global_csv  <- file.path(dir_processed, "Peat-DBase_global_5min.csv")
out_summary_csv <- file.path(dir_processed, "Peat-DBase_global_5min_summary.csv")

# ---------------------------- #
# Input files
# ---------------------------- #

peatdb_file  <- file.path(dir_peatdb, "Peat_DBase_version_1_0_0_b.csv")

if (!file.exists(peatdb_file)) {
  stop(
    "Peat-DBase CSV file not found.\n\n",
    "This project requires:\n",
    "Peat_DBase_version_1_0_0_b.csv\n\n",
    "Expected location:\n",
    peatdb_file, "\n\n",
    "Download the exact dataset version from:\n",
    "https://doi.org/10.5281/zenodo.15530644"
  )
}
df <- read_csv(peatdb_file, show_col_types = FALSE)
message("Peat-DBase CSV loaded: ", peatdb_file)

template_tif <- file.path(dir_wc, "wc2.1_5m_bio_1.tif")
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


template <- rast(template_tif)[[1]]
message("Template raster loaded: ", template_tif)
message("Resolution: ", paste(res(template), collapse = ", "))
message("CRS: ", crs(template))

if (nlyr(template) != 1) {
  stop("Template raster must have exactly one layer.")
}

if (!isTRUE(all.equal(res(template), c(5/60, 5/60), tolerance = 1e-8))) {
  stop("Template raster does not have 5 arc-min resolution.")
}

# ---------------------------- #
# Processing
# ---------------------------- #
# See Skye et al. (2025) for detailed data-cleaning steps and rationale.

n_raw <- nrow(df)

# 1. Remove records flagged as erroneous
df1 <- df %>%
  filter(error_found == FALSE)

# 2. Keep accepted duplication-flag categories
df2 <- df1 %>%
  filter(sample_duplication_flag %in% c(1, 2, 3, 4, 5))

# 3. Keep peat-presence records only
df3 <- df2 %>%
  filter(depth_cm >= 30)

# 4. Assign records to the 5 arc-min WorldClim grid
df4 <- df3 %>%
  mutate(grid5m_cell = cellFromXY(template, cbind(lon, lat)))

# 5. Remove records outside template coverage
pts_vect <- terra::vect(df4, geom = c("lon", "lat"), crs = "EPSG:4326")
df4$bio1 <- terra::extract(template, pts_vect)[, 2]

df5 <- df4 %>%
  filter(!is.na(grid5m_cell), !is.na(bio1))

# 6. Keep one record per 5 arc-min grid cell:
#    the one with the smallest (location_id, sample_id)
df6 <- df5 %>%
  group_by(grid5m_cell) %>%
  slice_min(
    order_by = tibble(location_id, sample_id),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

if (nrow(df6) == 0) {
  stop("No Peat-DBase occurrence records remain after filtering.")
}

# 7. Final SDM-style output
global <- df6 %>%
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
message("Global Peat-DBase CSV saved: ", out_global_csv)

processing_summary <- tibble(
  step = c(
    "RAW",
    "ERROR_FILTERED",
    "DUPLICATION_FILTERED",
    "DEPTH_GE_30",
    "WITH_TEMPLATE_COVERAGE",
    "GRID_CELL_DEDUPLICATED",
    "FINAL_OUTPUT"
  ),
  n_records = c(
    n_raw,
    nrow(df1),
    nrow(df2),
    nrow(df3),
    nrow(df5),
    nrow(df6),
    nrow(global)
  )
)

write_csv(processing_summary, out_summary_csv)
message("Processing summary saved: ", out_summary_csv)

