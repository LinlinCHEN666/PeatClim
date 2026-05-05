# ============================================================
# Script: 10c_evaluate_HTLTPeatland.R
# Purpose:
#   Evaluate the weighted HTLTPeatland probability map for the
#   WorldClim_5min case by calculating TP, FN, and sensitivity
#   against peatland occurrence data at 5 arc-min resolution.
#
# Inputs:
#   - Weighted HTLTPeatland projection raster from
#     10b_integrate_and_project_HTLTPeatland.R
#   - Global peatland occurrence data at 5 arc-min resolution from 
#     01c_integrate_PEATMAP_Peat-DBase.R
#
# Outputs:
#   - CSV file of TP, FN, and sensitivity for user-defined thresholds
#   - Maps showing thresholded predictions (PNG and PDF)
#
# Notes:
#   - This script is restricted to climate_id == "WordlClim_5min".
#   - Evaluation is performed against global peatland 5 arc-min occurrence data.
#   - Sensitivity is defined as TP / (TP + FN).
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(sf)
  library(scales)
  library(here)
})

here::i_am("code/10c_evaluate_HTLTPeatland.R")

# ---------------------------- #
# User settings
# ---------------------------- #

jobht <- "EM01_GBM013_WC5min"
joblt <- "EM01_GBM012_WC5min"
climate_id <- "WorldClim_5min"
myRespName <- "HTLTPeatland"

thresholds_prob <- c(0.5, 0.6, 0.7, 0.8)

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_occurrence <- here("data_processed", "occurrence")
dir_proj       <- here("data_processed", "projections")
dir_figures    <- here("outputs")

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

proj_path <- file.path(dir_proj, paste0(jobht, "_", joblt, "_", climate_id, ".tif"))
out_csv   <- file.path(
  dir_figures,
  paste0("10c_HTLTPeatland_sensitivity_", jobht, "_", joblt, ".csv")
)

# ---------------------------- #
# read peatland occurrence data
# ---------------------------- #

peatland_path <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min.csv")

if (!file.exists(peatland_path)) {
  stop("Peatland occurrence file not found: ", peatland_path)
}

peatland_data <- read_csv(peatland_path, show_col_types = FALSE)

if (!all(c("lon", "lat", "VALUE") %in% names(peatland_data))) {
  stop("peatland_data must contain columns: lon, lat, VALUE")
}

message("Loaded global peatland 5 arc-min occurrence data: ", nrow(peatland_data), " records")

# ---------------------------- #
# read HTLTPeatland projection raster
# ---------------------------- #

if (!file.exists(proj_path)) {
  stop("Projection tif file not found: ", proj_path)
}

p_global <- rast(proj_path)

if (ext(p_global)[2] > 181) {
  p_global <- terra::rotate(p_global)
}

prob <- terra::extract(p_global, peatland_data[, c("lon", "lat")])
colnames(prob) <- c("ID", "prob")

# ---------------------------- #
# Helper function:
# calculate TP, FN, sensitivity
# ---------------------------- #

calculate_sensitivity_htlt <- function(threshold_prob) {
  threshold <- threshold_prob * 1000
  eval_df <- bind_cols(peatland_data, prob["prob"]) %>%
    mutate(
      pred_presence = ifelse(prob >= threshold, 1L, 0L)
    )
  
  TP <- sum(eval_df$VALUE == 1 & eval_df$pred_presence == 1, na.rm = TRUE)
  FN <- sum(eval_df$VALUE == 1 & eval_df$pred_presence == 0, na.rm = TRUE)
  
  sensitivity <- if ((TP + FN) == 0) NA_real_ else TP / (TP + FN)
  
  tibble(
    threshold_prob,
    TP = TP,
    FN = FN,
    sensitivity = sensitivity
  )
}

# ---------------------------- #
# Optional helper function:
# thresholded Robinson map
# ---------------------------- #

plot_threshold_map <- function(threshold_prob, out_png, out_pdf) {
  
  threshold <- threshold_prob * 1000
  r <- ifel(p_global >= threshold, p_global, 0) / 1000
  crs(r) <- "EPSG:4326"
  r <- crop(r, ext(-180, 180, -90, 90), snap = "in")
  
  robin_crs <- st_crs("ESRI:54030")
  r_proj <- project(r, "ESRI:54030", method = "bilinear")
  
  dx <- res(r_proj)[1]
  dy <- res(r_proj)[2]
  
  df_r <- as.data.frame(r_proj, xy = TRUE, na.rm = TRUE)
  names(df_r)[3] <- "prob"
  
  grat_ll <- st_graticule(
    lat = seq(-60, 60, by = 30),
    lon = seq(-180, 180, by = 60)
  )
  grat_r <- st_transform(st_as_sf(grat_ll), robin_crs)
  
  make_robinson_frame_sf <- function(crs_out, n = 720, eps = 1e-6) {
    lon <- seq(-180, 180, length.out = n)
    lat <- seq(-90 + eps, 90 - eps, length.out = n)
    
    top    <- cbind(lon, rep( 90 - eps, n))
    right  <- cbind(rep(180, n), rev(lat))
    bottom <- cbind(rev(lon), rep(-90 + eps, n))
    left   <- cbind(rep(-180, n), lat)
    
    xy <- rbind(top, right, bottom, left, top[1, ])
    poly_ll <- st_sfc(st_polygon(list(xy)), crs = 4326)
    st_transform(st_as_sf(poly_ll), crs_out)
  }
  
  frame_r <- make_robinson_frame_sf(robin_crs)
  bb <- st_bbox(frame_r)
  
  # add peatland occurrence points for reference
  
  df_point <- peatland_data
  
  v <- vect(df_point, geom = c("lon", "lat"), crs = "EPSG:4326")
  res_5m <- 5 / 60
  
  # global raster template
  r_temp <- rast(
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    resolution = res_5m,
    crs = "EPSG:4326"
  )
  
  # rasterize points to presence grid
  r_point <- rasterize(v, r_temp, field = "VALUE", fun = "max", background = NA)
  
  # project to Robinson
  r_point_robin <- project(r_point, "ESRI:54030", method = "near")
  
  # convert to dataframe for ggplot
  df_point_robin <- as.data.frame(r_point_robin, xy = TRUE, na.rm = TRUE)
  names(df_point_robin)[3] <- "presence"
  
  p <- ggplot() +
    geom_sf(data = grat_r, color = "grey80", linewidth = 0.35) +
    geom_tile(data = df_r, aes(x = x, y = y, fill = prob), width = dx, height = dy) +
    geom_tile(data = df_point_robin, aes(x = x, y = y),
              fill = "#E31A1C", alpha = 0.6,  width = res(r_point_robin)[1], height = res(r_point_robin)[2]) +
    
    geom_sf(data = frame_r, fill = NA, color = "black", linewidth = 0.9) +
    scale_fill_viridis_c(
      limits = c(0, 1),
      oob = squish,
      na.value = "white",
      name = "Probability",
      breaks = seq(0, 1, 0.2),
      labels = sprintf("%.2f", seq(0, 1, 0.2))
    ) +
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
    labs(title = paste0("    threshold = ", threshold_prob)) +
    theme_void(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      legend.background = element_rect(fill = "white", color = NA),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 2)),
      legend.position = "bottom",
      legend.direction = "horizontal",
      plot.margin = ggplot2::margin(8, 0, 0, 0)
    )
  
  ggsave(out_png, plot = p, width = 12, height = 6.5, dpi = 300, bg = "white")
  ggsave(out_pdf, plot = p, width = 12, height = 6.5, dpi = 300, bg = "white")
  
}

# ---------------------------- #
# Run evaluation
# ---------------------------- #

results <- bind_rows(lapply(thresholds_prob, function(thr) {
  calculate_sensitivity_htlt(
    threshold = thr
  )
}))

write_csv(results, out_csv)
message("Saved sensitivity summary: ", out_csv)

# threshold maps
for (thr in thresholds_prob) {
  message("Generating threshold map for threshold = ", thr)
  out_png <- file.path(
    dir_figures,
    paste0("10c_HTLTPeatland_threshold_", thr, ".png")
  )
  out_pdf <- file.path(
    dir_figures,
    paste0("10c_HTLTPeatland_threshold_", thr, ".pdf")
  )
  plot_threshold_map(
    threshold_prob = thr,
    out_png = out_png,
    out_pdf = out_pdf
  )
}
