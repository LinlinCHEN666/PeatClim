# ============================================================
# Script: 02a_upscale_WorldClim_to_N48.R
# Purpose:
#   Upscale WorldClim v2.1 bioclimatic variables from 5 arc-min
#   resolution to N48 resolution and mask them using a 
#   HadCM3BL Pre-industrial land-sea mask .
#
# Inputs:
#   - WorldClim v2.1 bioclimatic rasters at 5 arc-min resolution
#     (bio1 to bio19)
#   - Target N48 land-sea mask raster
#
# Outputs:
#   - WorldClim bioclimatic rasters aggregated to N48 resolution
#     and masked to land cells only
#
# Notes:
#   - The target grid is defined by the supplied land-sea mask.
#   - Ocean cells in the output are set to NA using the mask.
#   - Output rasters are rotated to longitude range [-180, 180].
#
# External data requirements:
#
#   WorldClim v2.1 bioclimatic variables
#   - Citation:
#       Fick, S. E. and Hijmans, R. J. (2017).
#       WorldClim 2: new 1 km spatial resolution climate surfaces
#       for global land areas. International Journal of Climatology,
#       37(12), 4302-4315.
#   - Data used here:
#       WorldClim v2.1 bioclimatic variables at 5 arc-min resolution
#   - Download from:
#       https://www.worldclim.org/data/worldclim21.html
#   - Place extracted files under:
#       data_external/WorldClim/
#
#   Target land-sea mask
#   - Place the mask file under:
#       data_external/land_sea_mask/
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

here::i_am("code/02a_upscale_WorldClim_to_N48.R")

# ---------------------------- #
# Helper function
# ---------------------------- #

create_directory <- function(path) {
  if (length(path) != 1) {
    stop("Error: 'path' must be a single string.")
  }
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(path)) {
    stop("Failed to create directory: ", path)
  }
}

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_wc_input   <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_mask       <- here("data_external", "HadCM3BL")
dir_wc_output  <- here("data_processed", "WorldClim_N48", "wc2.1_N48_bio")

create_directory(dir_wc_output)

mask_file <- file.path(dir_mask, "tdezc1.qrparm.mask_lsm.nc")

# ---------------------------- #
# Input checks
# ---------------------------- #

bio_files <- file.path(dir_wc_input, paste0("wc2.1_5m_bio_", 1:19, ".tif"))

missing_bio <- bio_files[!file.exists(bio_files)]
if (length(missing_bio) > 0) {
  stop(
    "Missing WorldClim input files:\n",
    paste(missing_bio, collapse = "\n")
  )
}

if (!file.exists(mask_file)) {
  stop("Land-sea mask file not found: ", mask_file)
}

# ---------------------------- #
# Load target land-sea mask
# ---------------------------- #

land_mask <- rast(mask_file)

# Convert ocean cells to NA
land_mask[land_mask == 0] <- NA

# Ensure CRS is defined
crs(land_mask) <- "EPSG:4326"

target_ext <- ext(land_mask)
target_res <- res(land_mask)
target_crs <- crs(land_mask)

message("Loaded target land-sea mask: ", mask_file)
message("Target extent: ", paste(as.vector(target_ext), collapse = ", "))
message("Target resolution: ", paste(target_res, collapse = ", "))
message("Target CRS: ", target_crs)

# Create target grid aligned to the land-sea mask
target_grid <- rast(
  ext = target_ext,
  resolution = target_res,
  crs = target_crs
)

# ---------------------------- #
# Aggregate and mask all 19 variables
# ---------------------------- #

for (b in 1:19) {
  message("Processing bio", b, "...")
  
  input_file <- file.path(dir_wc_input, paste0("wc2.1_5m_bio_", b, ".tif"))
  input_raster <- rast(input_file)
  
  # Aggregate 5 arc-min data to target N48 resolution using mean
  agg_factor <- round(target_res / res(input_raster))
  
  agg_raster <- aggregate(
    input_raster,
    fact = agg_factor,
    fun = mean,
    na.rm = TRUE
  )
  
  # Align aggregated raster to exact target grid
  agg_projected <- project(
    agg_raster,
    target_grid,
    method = "bilinear"
  )
  
  # Mask ocean cells
  masked_raster <- mask(agg_projected, land_mask)
  
  # Rotate longitude from [0, 360] to [-180, 180] if needed
  masked_raster_rot <- terra::rotate(masked_raster)
  
  output_file <- file.path(dir_wc_output, paste0("wc2.1_N48_bio_", b, ".tif"))
  writeRaster(masked_raster_rot, output_file, overwrite = TRUE)
  
  message("Saved: ", output_file)
}

message("Finished preparing WorldClim N48 rasters.")
