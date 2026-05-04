# ============================================================
# Script: 01a_prepare_PEATMAP_occurrence.R
# Purpose:
#   Prepare modern peatland occurrence points from PEATMAP polygons
#   on the WorldClim 5 arc-min grid for downstream SDM use.
#
# Inputs:
#   - PEATMAP regional polygon shapefiles
#   - WorldClim v2.1 5 arc-min bioclimatic raster 
#     (bio1 used only as grid template)
#
# Outputs:
#   - Global merged PEATMAP occurrence CSV
#   - CSV summary of counts by region and processing step
#
# Notes:
#   - Presence is defined by polygon-to-grid rasterization.
#   - With touches = FALSE, a cell is present only if the cell centre falls within a peat polygon.
#   - With touches = TRUE, any touched cell counts as presence. (An alternative data-conversion rule)
#   - The final output keeps one record per 5 arc-min grid cell.
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
#   - Associated publication:
#       Xu, J., Morris, P. J., Liu, J., and Holden, J. (2018). 
#       PEATMAP: Refining estimates of global peatland distribution 
#       based on a meta-analysis. CATENA, 160, 134-140.
#   - Place extracted files under:
#       data_external/PEATMAP/
#
#   WorldClim v2.1 bioclimatic variables
#   - Citation:
#       Fick, S. E. and Hijmans, R. J. (2017).
#       WorldClim 2: new 1 km spatial resolution climate surfaces
#       for global land areas. International Journal of Climatology,
#       37(12), 4302-4315.
#   - Data used:
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

dir_peatmap   <- here("data_external", "PEATMAP")
dir_wc   <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_processed     <- here("data_processed", "occurrence")

create_directory(dir_processed)

# Final outputs
if (touches_rule == 'FALSE'){
  out_global_csv   <- file.path(dir_processed, paste0("PEATMAP_global_5min.csv"))
  out_summary_csv  <- file.path(dir_processed, paste0("PEATMAP_global_5min_summary.csv"))
} else if (touches_rules == 'TRUE'){
  out_global_csv   <- file.path(dir_processed, paste0("PEATMAP_global_5min_alt.csv"))
  out_summary_csv  <- file.path(dir_processed, paste0("PEATMAP_global_5min_summary_alt.csv"))
}

# ---------------------------- #
# Input files
# ---------------------------- #
# PEATMAP is provided as a set of regional shapefiles. Define paths to each region's shapefile.
peat_files <- c(
  Africa        = file.path(dir_peatmap, "Africa", "AF_Peatland.shp"),
  EA            = file.path(dir_peatmap, "Asia", "EA_Peatland.shp"),
  NNA           = file.path(dir_peatmap, "Asia", "Histosols_Hokkaido_Mongolia_North_Korea.shp"),
  NEA           = file.path(dir_peatmap, "Asia", "NEA_Peatland.shp"),
  SEA           = file.path(dir_peatmap, "Asia", "SEA_Peatland.shp"),
  SIB           = file.path(dir_peatmap, "Asia", "SIB_Peatland.shp"),
  BIE           = file.path(dir_peatmap, "Europe", "British Isles Peatland.shp"),
  FE            = file.path(dir_peatmap, "Europe", "Finland_Peatland.shp"),
  NE            = file.path(dir_peatmap, "Europe", "Norway_Peatland.shp"),
  OE            = file.path(dir_peatmap, "Europe", "Other_European_Peatland.shp"),
  SE            = file.path(dir_peatmap, "Europe", "Sweden_Peatland.shp"),
  CNA           = file.path(dir_peatmap, "North_America", "Canada_Peatland.shp"),
  ONA           = file.path(dir_peatmap, "North_America", "Other_North_American_Peatlands.shp"),
  USA           = file.path(dir_peatmap, "North_America", "USA_Peatland.shp"),
  Oceania       = file.path(dir_peatmap, "Oceania", "Oceania_Peatland.shp"),
  South_America = file.path(dir_peatmap, "South_America", "SA_Peatland.shp")
)

missing_shps <- peat_files[!file.exists(peat_files)]

if (length(missing_shps) > 0) {
  stop(
    "Required PEATMAP shapefiles were not found.\n\n",
    "Download the PEATMAP dataset from:\n",
    "https://archive.researchdata.leeds.ac.uk/251/\n",
    "DOI: 10.5518/252\n\n",
    "After downloading, extract the zip files into:\n",
    here("data_external", "PEATMAP"), "\n\n",
    "Missing files:\n",
    paste(names(missing_shps), missing_shps, sep = " -> ", collapse = "\n")
  )
}


# WorldClim template raster (bio1 used only as grid template)
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

# sanity checks
if (nlyr(template) != 1) {
  stop("Template raster must have exactly one layer.")
}

if (!isTRUE(all.equal(res(template), c(5/60, 5/60), tolerance = 1e-8))) {
  stop("Template raster does not have 5 arc-min resolution.")
}

# ---------------------------- #
# Helper function
# ---------------------------- #

process_region <- function(region_name,
                           shp_path,
                           template_raster,
                           touches = FALSE) {
  message("\n==============================")
  message("Processing region: ", region_name)
  message("==============================")
  
  # Read shapefile
  peat <- vect(shp_path)
  
  n_features_raw <- nrow(peat)
  
  if (n_features_raw == 0) {
    warning("No features found in shapefile for region: ", region_name)
    peat_df <- tibble(
      lon = numeric(0),
      lat = numeric(0),
      region = character(0)
    )
    
    summary_row <- tibble(
      region = region_name,
      shapefile = shp_path,
      n_features_raw = 0L,
      n_presence_cells = 0L
    )
    
    return(list(points = peat_df, summary = summary_row))
  }
  
  # Reproject to template CRS if needed
  if (!same.crs(peat, template_raster)) {
    peat <- project(peat, crs(template_raster))
  }
  
  # Crop template to region extent to speed up rasterization
  template_crop <- crop(template_raster, ext(peat), snap = "out")
  
  # Rasterize polygons to the template grid
  peat_raster <- rasterize(
    x = peat,
    y = template_crop,
    field = 1,
    background = 0,
    touches = touches
  )
  
  names(peat_raster) <- "peat_presence"
  
  # Extract occupied cells
  vals <- values(peat_raster, mat = FALSE)
  peat_cells <- which(vals == 1)
  
  if (length(peat_cells) == 0) {
    warning("No peat cells found after rasterization for region: ", region_name)
    peat_df <- tibble(
      lon = numeric(0),
      lat = numeric(0),
      region = character(0)
    )
  } else {
    xy <- xyFromCell(peat_raster, peat_cells)
    peat_df <- tibble(
      lon = xy[, 1],
      lat = xy[, 2],
      region = region_name
    )
  }
  
  message("Presence cells: ", nrow(peat_df))
  
  summary_row <- tibble(
    region = region_name,
    shapefile = shp_path,
    n_features_raw = n_features_raw,
    n_presence_cells = nrow(peat_df)
  )
  
  return(list(points = peat_df, summary = summary_row))
}

# ---------------------------- #
# Process all PEATMAP regions
# ---------------------------- #

all_results <- lapply(seq_along(peat_files), function(i) {
  process_region(
    region_name = names(peat_files)[i],
    shp_path = peat_files[i],
    template_raster = template,
    touches = touches_rule
  )
})

all_points_list <- lapply(all_results, `[[`, "points")
all_summaries   <- lapply(all_results, `[[`, "summary")

all_points_raw <- bind_rows(all_points_list)
regional_summary <- bind_rows(all_summaries)

message("\nTotal rows before removing duplicates: ", nrow(all_points_raw))

if (nrow(all_points_raw) == 0) {
  stop("No PEATMAP occurrence cells were generated from any regional shapefile.")
}

# ---------------------------- #
# Remove duplicate grid cells
# ---------------------------- #
# Regional files may overlap; final product should keep one row per 5 arc-min cell.

all_points_unique <- all_points_raw %>%
  distinct(lon, lat)

message("Total unique grid cells after duplicate removal: ", nrow(all_points_unique))

if (nrow(all_points_unique) == 0) {
  stop("No unique PEATMAP occurrence cells remain after duplicate removal.")
}

# ---------------------------- #
# Drop records outside template data coverage
# ---------------------------- #

pts_vect <- terra::vect(all_points_unique, geom = c("lon", "lat"), crs = "EPSG:4326")
all_points_unique$bio1 <- terra::extract(template, pts_vect)[, 2]


all_points_filtered <- all_points_unique %>%
  filter(!is.na(bio1)) %>%
  mutate(
    fid = seq_len(n()),
    VALUE = 1L
  ) %>%
  select(lon, lat, fid, VALUE)

message("Total grid cells after removing cells outside template coverage: ", nrow(all_points_filtered))

# ---------------------------- #
# Write final outputs
# ---------------------------- #

write_csv(all_points_filtered, out_global_csv)
message("Global PEATMAP CSV saved: ", out_global_csv)

processing_summary <- regional_summary %>%
  mutate(touches = touches_rule) %>%
  bind_rows(
    tibble(
      region = "GLOBAL",
      shapefile = NA_character_,
      n_features_raw = NA_integer_,
      n_presence_cells = nrow(all_points_raw),
      touches = touches_rule
    ),
    tibble(
      region = "GLOBAL_UNIQUE",
      shapefile = NA_character_,
      n_features_raw = NA_integer_,
      n_presence_cells = nrow(all_points_unique),
      touches = touches_rule
    ),
    tibble(
      region = "GLOBAL_FILTERED",
      shapefile = NA_character_,
      n_features_raw = NA_integer_,
      n_presence_cells = nrow(all_points_filtered),
      touches = touches_rule
    )
  )

write_csv(processing_summary, out_summary_csv)
message("Processing summary saved: ", out_summary_csv)
