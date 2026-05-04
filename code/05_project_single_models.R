# ============================================================
# Script: 05_project_single_models.R
# Purpose:
#   Project trained single biomod2 models to a selected climate scenario 
#   using WorldClim or other prepared climate inputs.
#
# Inputs:
#   - Config CSV file
#   - Single-model output from 04_build_single_models.R
#   - Climate predictor rasters for the selected climate_id
#
# Outputs:
#   - biomod2 single-model projection object
#   - projection plot of the clamping mask (PNG and PDF)
#   - projection plot of individual runs (PNG and PDF)
#
# Notes:
#   - This script is intended to be run in HPC by command line:
#       Rscript code/05_project_single_models.R config/config_project_single_models.csv <jobid>
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

here::i_am("code/05_project_single_models.R")

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

myRespName   <- params$myRespName[[1]]
nb_cpu       <- params$nb.cpu[[1]]
env_variables <- eval(parse(text = params$env.variables[[1]]))
jobid_build  <- params$jobid_build[[1]]
climate_id   <- params$climate_id[[1]]

message("jobid: ", jobid)
message("myRespName: ", myRespName)
message("jobid_build: ", jobid_build)
message("climate_id: ", climate_id)


# ---------------------------- #
# Project paths
# ---------------------------- #

dir_wc          <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_proj        <- here("outputs", "05_projections")
job_proj_dir  <- file.path(dir_proj, myRespName, jobid)

create_directory(dir_proj)
create_directory(job_proj_dir)

# ---------------------------- #
# Load climate predictors
# ---------------------------- #

message("Load environmental variables")

if (climate_id == "WorldClim_5min") {
  env_files <- file.path(dir_wc, paste0("wc2.1_5m_bio_", env_variables, ".tif"))
  
  missing_env <- env_files[!file.exists(env_files)]
  if (length(missing_env) > 0) {
    stop(
      "Missing WorldClim predictor files:\n",
      paste(missing_env, collapse = "\n")
    )
  }
  
  myExpl <- rast(env_files)
  names(myExpl) <- paste0("bio", env_variables)
  crs(myExpl) <- "EPSG:4326"
  
} else {
  stop("Unsupported climate_id in current script: ", climate_id,
       "\nPlease prepare predictors and add support for this climate_id in the script.")
}

message("Environmental predictors loaded: ", paste(names(myExpl), collapse = ", "))

# ---------------------------- #
# Load trained model output
# ---------------------------- #
model_out_path <- file.path(
  here(),
  myRespName,
  paste0(myRespName, ".", jobid_build, ".models.out")
)

if (!file.exists(model_out_path)) {
  stop(
    "Trained model object not found: ", model_out_path, "\n",
    "Run 04_build_single_models.R first for jobid_build = ", jobid_build
  )
}

myBiomodModelOut <- get(load(model_out_path))

message("Loaded model object: ", model_out_path)
message("Built models available:")
print(get_built_models(myBiomodModelOut))

# ---------------------------- #
# Project single models
# ---------------------------- #

message("Project single models")

selected_models <- get_built_models(myBiomodModelOut)[
  !grepl("allData|allRun", get_built_models(myBiomodModelOut))
]

myBiomodProj <- BIOMOD_Projection(
  bm.mod = myBiomodModelOut,
  proj.name = jobid,
  new.env = myExpl,
  models.chosen = selected_models,
  build.clamping.mask = TRUE,
  seed.val = 42,
  nb.cpu = nb_cpu
)

message("Single-model project completed. Now plot clamping and projection figures...")


# ---------------------------- #
# Helper: plot Robinson function
# ---------------------------- #

plot_robinson_rasters <- function(
    r,
    out_dir,
    filename_prefix,
    type = "probability",
    width = 10,
    height = 6,
    dpi = 300,
    save_pdf = TRUE
) {
  # Required packages:
  # terra, sf, ggplot2, scales, viridis
  
  if (!inherits(r, "SpatRaster")) {
    stop("'r' must be a terra SpatRaster.")
  }
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ------------------------------------------------------------
  # Prepare raster
  # ------------------------------------------------------------
  if (type == "probability") {
    r1 <- r / 1000
  } else {
    r1 <- r
  }
  
  crs(r1) <- "EPSG:4326"
  r1 <- crop(r1, ext(-180, 180, -90, 90), snap = "in")
  
  robin_crs <- st_crs("ESRI:54030")
  
  if (type == "clamp") {
    r_proj <- project(r1, "ESRI:54030", method = "near")
  } else {
    r_proj <- project(r1, "ESRI:54030", method = "bilinear")
  }
  
  dx <- res(r_proj)[1]
  dy <- res(r_proj)[2]
  
  # ------------------------------------------------------------
  # Graticules
  # ------------------------------------------------------------
  grat_ll <- st_graticule(
    lat = seq(-60, 60, by = 30),
    lon = seq(-180, 180, by = 60)
  )
  grat_r <- st_transform(st_as_sf(grat_ll), robin_crs)
  
  # ------------------------------------------------------------
  # Robinson frame
  # ------------------------------------------------------------
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
  
  frame_v <- vect(frame_r)
  r_proj_m <- mask(r_proj, frame_v)
  
  # ------------------------------------------------------------
  # Plot each layer
  # ------------------------------------------------------------
  for (i in seq_len(nlyr(r_proj_m))) {
    
    df_r <- as.data.frame(r_proj_m[[i]], xy = TRUE, na.rm = TRUE)
    names(df_r)[3] <- "value"
    
    if (type == "clamp") {
      df_r$value <- factor(df_r$value, levels = sort(unique(df_r$value)))
    }
    
    if (type == "probability") {
      plot_title <- names(r[[i]])
      
      fill_scale <- scale_fill_viridis_c(
        option = "viridis",
        limits = c(0, 1),
        oob = squish,
        na.value = "white",
        name = "Probability",
        breaks = seq(0, 1, 0.2),
        labels = sprintf("%.2f", seq(0, 1, 0.2))
      )
      
      fill_guide <- guides(
        fill = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(18, "cm"),
          barheight = unit(0.6, "cm"),
          ticks = TRUE
        )
      )
      
    } else if (type == "clamp") {
      plot_title <- "Clamping mask"
      
      fill_scale <- scale_fill_viridis_d(
        option = "viridis",
        na.value = "white",
        name = "Value"
      )
      
      fill_guide <- guides(
        fill = guide_legend(
          title.position = "top",
          title.hjust = 0.5
        )
      )
      
    } else {
      plot_title <- names(r[[i]])
      
      fill_scale <- scale_fill_viridis_c(
        option = "viridis",
        na.value = "white"
      )
      
      fill_guide <- guides(
        fill = guide_colorbar(
          title.position = "top",
          title.hjust = 0.5,
          barwidth = unit(18, "cm"),
          barheight = unit(0.6, "cm"),
          ticks = TRUE
        )
      )
    }
    
    p <- ggplot() +
      geom_sf(data = grat_r, color = "grey80", linewidth = 0.35) +
      geom_tile(
        data = df_r,
        aes(x = x, y = y, fill = value),
        width = dx,
        height = dy
      ) +
      geom_sf(data = frame_r, fill = NA, color = "black", linewidth = 0.9) +
      fill_scale +
      fill_guide +
      coord_sf(
        crs = robin_crs,
        datum = NA,
        xlim = c(bb["xmin"], bb["xmax"]),
        ylim = c(bb["ymin"], bb["ymax"]),
        expand = FALSE
      ) +
      labs(title = plot_title) +
      theme_void(base_size = 12) +
      theme(
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        legend.background = element_rect(fill = "white", color = NA),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        legend.box.margin = ggplot2::margin(0, 0, 0, 0),
        plot.title = element_text(
          hjust = 0.5,
          face = "bold",
          size = 14,
          margin = ggplot2::margin(b = 2)
        ),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = ggplot2::margin(8, 0, 0, 0)
      )
    
    ggsave(
      filename = file.path(
        out_dir,
        
        paste0(filename_prefix, "_", names(r[[i]]), ".png")
      ),
      plot = p,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
    
    if (save_pdf) {
      ggsave(
        filename = file.path(
          out_dir,
          paste0(filename_prefix, "_", names(r[[i]]), ".pdf")
        ),
        plot = p,
        width = width,
        height = height,
        bg = "white"
      )
    }
  }
}

# ---------------------------- #
# plot each projected layer separately, in Robinson projection
# ---------------------------- #

r_prob <- rast(myBiomodProj@proj.out@link[1])

plot_robinson_rasters(
  r = r_prob,
  out_dir = job_proj_dir,
  filename_prefix = paste0("proj_", jobid),
  type = "probability",
  save_pdf = TRUE
)

# ---------------------------- #
# plot ClampingMask
# ---------------------------- #

# #This mask values will correspond to the number of variables in each pixel that are out of their calibration / validation range, identifying locations where predictions are uncertain.
r_clamp_path <- file.path(here(),myRespName,paste0("proj_",jobid),paste0("proj_",jobid,"_ClampingMask.tif"))
r_clamp <- rast(r_clamp_path)

plot_robinson_rasters(
  r = r_clamp,
  out_dir = job_proj_dir,
  filename_prefix = paste0("proj_", jobid, "_ClampingMask"),
  type = "clamp",
  save_pdf = TRUE
)
