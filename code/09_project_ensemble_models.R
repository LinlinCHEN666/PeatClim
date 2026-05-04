# ============================================================
# Script: 09_project_ensemble_models.R
# Purpose:
#   Project biomod2 ensemble models to a selected climate scenario 
#   using prepared predictor rasters.
#
# Inputs:
#   - Config CSV file
#   - Ensemble-model output from 08_build_ensemble_models.R
#   - Climate predictor rasters for the selected climate_id
#
# Outputs:
#   - biomod2 ensemble projection object
#   - projection plot of ensemble mean (PNG and PDF)
#   - projection plot of ensemble cv (PNG and PDF)
#
# Notes:
#   - For this script, climate_id could choose "WorldClim_5min", "WorldClim_N48" or "PI". 
#     Modify script to include other climate data if needed.
#   - This script is intended to be run in HPC by command line:
#       Rscript code/09_project_ensemble_models.R config/config_project_ensemble_models.csv <jobid>
#   - Production runs may be submitted through SLURM using:
#       hpc/submit_r_job.sh
#   - If required on HPC, PROJ/GDAL environment variables should be set in the shell launcher,
#     e.g.:Sys.setenv(PROJ_LIB = ".../share/proj")
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(biomod2)
  library(terra)
  library(tibble)
  library(here)
  library(viridis)
  library(sf)
  library(ggplot2)
  library(scales)
})

here::i_am("code/09_project_ensemble_models.R")

# ---------------------------- #
# Helper functions
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
# Read command-line arguments
# ---------------------------- #

args <- commandArgs(trailingOnly = TRUE)
config_file <- args[1]
jobid <- args[2]

# ---------------------------- #
# Read config
# ---------------------------- #
config_path <- here(config_file)

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

config <- read.csv(
  config_path,
  comment.char = "#",
  stringsAsFactors = FALSE
)

params <- config[config$jobid == jobid, ]

if (nrow(params) == 0) {
  stop("No settings found for jobid: ", jobid)
}

if (nrow(params) != 1) {
  stop("Expected exactly 1 row for jobid ", jobid, ", found ", nrow(params))
}

# ---------------------------- #
# Extract parameters
# ---------------------------- #

myRespName   <- params$myRespName[[1]]
env_variables <- eval(parse(text = params$env_variables[[1]]))
jobid_ensemble <- params$jobid_ensemble[[1]]
nb_cpu       <- params$nb_cpu[[1]]
algo <- params$algo[[1]] #"all", c('EMmean','EMcv')
climate_id   <- params$climate_id[[1]]
age          <- params$age[[1]]

algo <- if (algo == "all") {
  "all"
}else {
  eval(parse(text = algo)) 
}


# ---------------------------- #
# Project paths
# ---------------------------- #

dir_proj    <- here("outputs", "09_projections")
job_proj_dir      <- file.path(dir_proj, myRespName, jobid)

create_directory(dir_proj)
create_directory(job_proj_dir)

# ---------------------------- #
# Load environmental predictors
# ---------------------------- #

message("Load environmental variables")
if (climate_id == "WorldClim_5min") {
  dir_clim <- here("data_external", "WorldClim", "wc2.1_5m_bio")
  env_files <- file.path(dir_clim, paste0("wc2.1_5m_bio_", env_variables, ".tif"))
  
  missing_env <- env_files[!file.exists(env_files)]
  if (length(missing_env) > 0) {
    stop(
      "Missing WorldClim predictor files:\n",
      paste(missing_env, collapse = "\n")
    )
  }
  
  myExpl <- rast(env_files)
  names(myExpl) <- paste0("bio", env_variables)
  if (ext(myExpl)[2] > 181) {
    myExpl <- terra::rotate(myExpl)
  }
  crs(myExpl) <- "EPSG:4326"
  
} else if (climate_id == "WorldClim_N48"){
  dir_clim <- here("data_processed","WorldClim_N48","wc2.1_N48_bio")
  env_files <- file.path(dir_clim, paste0("wc2.1_N48_bio_", env_variables, ".tif"))
  
  missing_env <- env_files[!file.exists(env_files)]
  if (length(missing_env) > 0) {
    stop(
      "Missing WorldClim predictor files:\n",
      paste(missing_env, collapse = "\n")
    )
  }
  
  myExpl <- rast(env_files)
  names(myExpl) <- paste0("bio", env_variables)
  if (ext(myExpl)[2] > 181) {
    myExpl <- terra::rotate(myExpl)
  }
  crs(myExpl) <- "EPSG:4326"
  
} else if (climate_id == "PI") {
  dir_clim <- here("data_processed","HadCM3BL")
  env_files <- file.path(dir_clim,"tdezc1_bioclim.nc")
  
  missing_env <- env_files[!file.exists(env_files)]
  if (length(missing_env) > 0) {
    stop(
      "Missing WorldClim predictor files:\n",
      paste(missing_env, collapse = "\n")
    )
  }
  
  myExpl_all <- rast(env_files)
  myExpl <- myExpl_all[[env_variables]]
  if (ext(myExpl)[2] > 181) {
    myExpl <- terra::rotate(myExpl)
  }
  crs(myExpl) <- "EPSG:4326"
} else {
  stop("Unsupported climate_id: ", climate_id,
       "\nPlease prepare predictors and add support for this climate_id in the script.")
}

message("Environmental predictors loaded: ", paste(names(myExpl), collapse = ", "))

# ---------------------------- #
# Load ensemble model object
# ---------------------------- #

ensemble_model_path <- file.path(
  here(),
  myRespName,
  paste0(myRespName,".", jobid_ensemble, ".ensemble.models.out"))

if (!file.exists(ensemble_model_path)) {
  stop(
    "Ensemble model object not found: ", ensemble_model_path, "\n",
    "Run 08_build_ensemble_models.R first for jobid_ensemble = ", jobid_ensemble
  )
}

myBiomodEM <- get(load(ensemble_model_path))

message("Loaded ensemble model object: ", ensemble_model_path)
message("Built ensemble models available:")
print(get_built_models(myBiomodEM))

# ---------------------------- #
# Select ensemble members
# ---------------------------- #

if (length(algo) == 1 && algo == "all") {
  models_to_project <- "all"
} else {
  models_to_project <- get_built_models(
    obj = myBiomodEM,
    algo = algo #c('EMmean', 'EMcv')
  )
}

cat("models_to_project:\n")
print(models_to_project)

# ---------------------------- #
# Project ensemble models
# ---------------------------- #

message("Project ensemble models")

myBiomodEMProj <- BIOMOD_EnsembleForecasting(
  bm.em = myBiomodEM,
  proj.name = jobid,
  new.env = myExpl,
  models.chosen = models_to_project,
  #metric.binary = "all",
  #metric.filter = "all",
  nb.cpu = nb_cpu,
  on_0_1000 = TRUE,
  do.stack = FALSE,
  keep.in.memory = FALSE
)

message("Ensemble-model project completed. Now plot projection figures...")
# ---------------------------- #
# plot ensemble mean and cv projection, in Robinson projection
# ---------------------------- #

## Prepare raster (0–1) and project to Robinson #####
r <- rast(myBiomodEMProj@proj.out@link[1:2]) # first 2 layers are EMmean and EMcv
r1 <- r/1000

crs(r1) <- "EPSG:4326"
r1 <- crop(r1, ext(-180, 180, -90, 90), snap="in")

robin_crs <- st_crs("ESRI:54030")
r_proj <- project(r1, "ESRI:54030", method="bilinear")

print(r_proj[[2]])

dx <- res(r_proj)[1]
dy <- res(r_proj)[2]


## Graticules #####
grat_ll <- st_graticule(
  lat = seq(-60, 60, by=30),
  lon = seq(-180, 180, by=60)
)
grat_r <- st_transform(st_as_sf(grat_ll), robin_crs)

## Oval Robinson frame polygon #####
make_robinson_frame_sf <- function(crs_out, n = 720, eps = 1e-6) {
  lon <- seq(-180, 180, length.out=n)
  lat <- seq(-90 + eps, 90 - eps, length.out=n)
  
  top    <- cbind(lon, rep( 90 - eps, n))
  right  <- cbind(rep(180, n), rev(lat))
  bottom <- cbind(rev(lon), rep(-90 + eps, n))
  left   <- cbind(rep(-180, n), lat)
  
  xy <- rbind(top, right, bottom, left, top[1,])
  
  poly_ll <- st_sfc(st_polygon(list(xy)), crs = 4326)
  st_transform(st_as_sf(poly_ll), crs_out)
}

frame_r <- make_robinson_frame_sf(robin_crs)
bb <- st_bbox(frame_r)

frame_v <- vect(frame_r)          # convert sf -> terra vect
r_proj_m <- mask(r_proj, frame_v) # sets cells outside oval to NA


## Plot #####
df_r <- list()
p <- list()
for (i in 1:length(names(r_proj))){
  
  df_r[[i]] <- as.data.frame(r_proj_m[[i]], xy = TRUE, na.rm = TRUE)
  names(df_r[[i]])[3] <- "prob"
  
  # conditional fill scale
  fill_scale <- if (i == 1) {
    scale_fill_viridis_c(
      option = "viridis",
      limits = c(0, 1),
      oob = squish,
      na.value = "white",
      name = "Probability",
      breaks = seq(0, 1, 0.2),
      labels = sprintf("%.2f", seq(0, 1, 0.2))
    )
  } else if (i == 2) {
    scale_fill_viridis_c(
      option = "viridis",
      na.value = "white",
      name = "CV"
    )
  }
  
  
  p[[i]] <- ggplot() +
    
    # graticules
    geom_sf(data = grat_r, color = "grey80", linewidth = 0.35) +
    
    #  use geom_tile
    geom_tile(data = df_r[[i]], aes(x = x, y = y, fill = prob),width = dx, height = dy) +
    
    # oval frame last
    geom_sf(data = frame_r, fill = NA, color = "black", linewidth = 0.9) +
    
    # colors + legend
    fill_scale +
    
    coord_sf(
      crs = robin_crs,
      datum = NA,
      xlim = c(bb["xmin"], bb["xmax"]),
      ylim = c(bb["ymin"], bb["ymax"]),
      expand = FALSE
    ) +
    
    guides(
      fill = guide_colorbar(
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(18, "cm"),
        barheight = unit(0.6, "cm"),
        ticks = TRUE
      )
    ) +
    
    labs(title = names(r[[i]])) +
    
    # WHITE everywhere
    theme_void(base_size = 12) +
    theme(
      # backgrounds
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      
      #legend.box.spacing = unit(1.2, "pt"),      # reduce gap between map and legend
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      # title formatting
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 14,margin = ggplot2::margin(b = 2)),
      #plot.subtitle = element_text(hjust = 0.5, size = 12, margin = ggplot2::margin(b = 2)),   # margin=increases gap between subtitle and map
      
      # legend formatting
      legend.position = "bottom",
      legend.direction = "horizontal",
      
      plot.margin = ggplot2::margin(8,0,0,0)
    )
  
  ggsave(
    filename = file.path(job_proj_dir,paste0("proj_",jobid,"_",names(r[[i]]),".png")),
    plot = p[[i]],
    width = 10, height = 6, dpi = 300, bg = "white"
  )
  ggsave(
    filename = file.path(job_proj_dir,paste0("proj_",jobid,"_",names(r[[i]]),".pdf")),
    plot = p[[i]],
    width = 10, height = 6, dpi = 300, bg = "white"
  )
}

message("Ensemble mean projection completed successfully.")
