# ============================================================
# Script: 02b_prepare_PI_19_biovars.R
# Purpose:
#   Create the 19 bioclimatic variables for the HadCM3BL Pre-industrial climate
#   simulation from monthly temperature and precipitation fields,
#   and mask the results using the corresponding land-sea mask.
#
# Inputs:
#   - Monthly maximum temperature NetCDF file
#   - Monthly minimum temperature NetCDF file
#   - Monthly precipitation NetCDF file
#   - land-sea mask NetCDF file
#
# Outputs:
#   - NetCDF file containing 19 masked bioclimatic variables
#
# Notes:
#   - This script is restricted to the PI (tdezc1) job.
#   - Monthly precipitation is converted from mm/day to mm/month
#     by multiplying by 30 before calculating bioclimatic variables.
#   - The 19 bioclimatic variables are generated using
#     dismo::biovars(prec, tmin, tmax).
#
# External data requirements:
#   - Monthly tdezc1 temperature and precipitation fields
#   - tdezc1 land-sea mask
# ============================================================

suppressPackageStartupMessages({
  library(dismo)
  library(raster)
  library(ncdf4)
  library(terra)
  library(here)
})

here::i_am("code/02b_prepare_PI_19_biovars.R")

# ---------------------------- #
# User settings
# ---------------------------- #

jobid <- "tdezc1"

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

dir_input <- here("data_external", "HadCM3BL")
dir_mask  <- here("data_external", "HadCM3BL")
dir_out   <- here("data_processed", "HadCM3BL")

create_directory(dir_out)

path_tmax <- file.path(dir_input, paste0(jobid, "_tempmax_av.nc"))
path_tmin <- file.path(dir_input, paste0(jobid, "_tempmin_av.nc"))
path_prec <- file.path(dir_input, paste0(jobid, "_precipmon_av.nc"))
path_lsm  <- file.path(dir_mask,  paste0(jobid, ".qrparm.mask_lsm.nc"))

path_out  <- file.path(dir_out,   paste0(jobid, "_bioclim.nc"))

# ---------------------------- #
# Input checks
# ---------------------------- #

required_files <- c(path_tmax, path_tmin, path_prec, path_lsm)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input files:\n",
    paste(missing_files, collapse = "\n")
  )
}

message("Processing jobid: ", jobid)

# ---------------------------- #
# Load monthly climate data
# ---------------------------- #

# Each file is expected to have 12 layers corresponding to months.
tmax <- stack(path_tmax)
tmin <- stack(path_tmin)
prec <- stack(path_prec) * 30  # convert from mm/day to mm/month

if (nlayers(tmax) != 12) stop("tmax does not contain 12 monthly layers.")
if (nlayers(tmin) != 12) stop("tmin does not contain 12 monthly layers.")
if (nlayers(prec) != 12) stop("prec does not contain 12 monthly layers.")

# ---------------------------- #
# Create bioclimatic variables
# ---------------------------- #

# biovars order: precipitation, tmin, tmax
bio <- biovars(prec, tmin, tmax)

message("Created ", nlayers(bio), " bioclimatic variables.")
message("Variable names: ", paste(names(bio), collapse = ", "))

# ---------------------------- #
# Apply land-sea mask
# ---------------------------- #

lsm <- stack(path_lsm)

bioclim_land <- terra::mask(bio, lsm, maskvalue = 0)

message("Applied land-sea mask.")

# ---------------------------- #
# Write output NetCDF
# ---------------------------- #

# Extract coordinate information from one layer
r <- bio[[1]]
lon_vals <- xFromCol(r, 1:ncol(r))
lat_vals <- yFromRow(r, 1:nrow(r))
lat_vals <- sort(lat_vals) # Ensure latitude values are in ascending order (south to north)
lat_dim <- ncdim_def("lat", "degrees_north", lat_vals)
lon_dim <- ncdim_def("lon", "degrees_east", lon_vals)

# Create variable definitions for each bioclimatic layer
bio_vars <- list()
for(i in 1:nlayers(bioclim_land)) {
  varname <- names(bioclim_land)[i]  # e.g., "bio1", "bio2", etc.
  # Define each variable with dimensions ordered as (lat, lon)
  bio_vars[[i]] <- ncvar_def(varname, units = "", dim = list(lon_dim, lat_dim),
                             missval = -9999, longname = varname, prec = "float")
}

if(file.exists(path_out)) file.remove(path_out)
nc_out <- nc_create(paste0(path_out), bio_vars)

for(i in 1:nlayers(bioclim_land)) {
  varname <- names(bioclim_land)[i]
  data_matrix <- as.matrix(bioclim_land[[i]])
  # Flip the rows so that the lowest latitude is first
  data_matrix <- data_matrix[nrow(data_matrix):1, ]
  # Transpose the matrix so that its dimensions become: 
  # first dimension = lat (73), second dimension = lon (96)
  data_matrix <- t(data_matrix)
  ncvar_put(nc_out, varname, data_matrix)
}

# Add axis attributes 
ncatt_put(nc_out, "lat", "axis", "Y")
ncatt_put(nc_out, "lon", "axis", "X")

nc_close(nc_out)
