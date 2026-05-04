# ============================================================
# Script: 10b_integrate_and_project_HTLTPeatland.R
# Purpose:
#   Combine HTPeatland and LTPeatland ensemble probabilities
#   into HTLTPeatland probability using a temperature-dependent
#   logistic weighting, and plot the resulting global map.
#
# Inputs:
#   - HTPeatland ensemble projection raster
#   - LTPeatland ensemble projection raster
#   - Climate data (bio1) of a selected climate scenario based on ensemble projection
#
# Outputs:
#   - Weighted HTLTPeatland raster (GeoTIFF)
#   - Global Robinson-projection map (PNG and PDF)
#
# Notes:
#   - For this script, climate_id could choose "WorldClim_5min", "WorldClim_N48" or "PI".
#     Modify script to include other climate data if needed.
#   - The blended probability is:
#       P_HTLTPeatland = w_HT * P_HT + w_LT * P_LT
#   - Logistic weights are defined from bio1 using split temperature t0.
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(viridis)
  library(ggplot2)
  library(scales)
  library(sf)
  library(here)
})

here::i_am("code/10b_integrate_and_project_HTLTPeatland.R")

# ---------------------------- #
# User settings
# ---------------------------- #

jobht <- "EM01_GBM013_WC5min"
joblt <- "EM01_GBM012_WC5min"
climate_id <- "WorldClim_5min" # "WorldClim_5min", "WorldClim_N48", or "PI"
t0 <- 11.37   # split temperature (°C)
w  <- 2       # transition width (°C)

# ---------------------------- #
# Project paths
# ---------------------------- #
dir_proj_ht   <- here("HTPeatland", paste0("proj_", jobht))
dir_proj_lt   <- here("LTPeatland", paste0("proj_", joblt))

if (climate_id == "WorldClim_5min") {

  dir_clim <- here("data_external", "WorldClim", "wc2.1_5m_bio")
  bio1_path <- file.path(dir_clim, "wc2.1_5m_bio_1.tif")
  if (!file.exists(bio1_path)) {
    stop("bio1 raster not found: ", bio1_path)
  }
  bio1 <- rast(bio1_path)
  crs(bio1) <- "EPSG:4326"
  
} else if (climate_id == "WorldClim_N48"){

  dir_clim <- here("data_processed","WorldClim_N48","wc2.1_N48_bio")
  bio1_path <- file.path(dir_clim, "wc2.1_N48_bio_1.tif")
  if (!file.exists(bio1_path)) {
    stop("bio1 raster not found: ", bio1_path)
  }
  bio1 <- rast(bio1_path)
  crs(bio1) <- "EPSG:4326"
  
} else if (climate_id == "PI") {

  dir_clim <- here("data_processed","HadCM3BL")
  bio_path <- file.path(dir_clim,"tdezc1_bioclim.nc")
  if (!file.exists(bio_path)) {
    stop("bio1 raster not found: ", bio_path)
  }
  bio <- rast(bio_path)
  bio1 <- bio[['bio1']]
  crs(bio1) <- "EPSG:4326"
}

dir_proj_htlt     <- here("data_processed","projections","HTLTPeatland")
dir_figures   <- here("outputs", "10b_projections")

dir.create(dir_proj_htlt, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

p_ht_path <- file.path(
  dir_proj_ht,
  "individual_projections",
  "HTPeatland_EMmeanByTSS_mergedData_mergedRun_mergedAlgo.tif"
)

p_lt_path <- file.path(
  dir_proj_lt,
  "individual_projections",
  "LTPeatland_EMmeanByTSS_mergedData_mergedRun_mergedAlgo.tif"
)

if (!file.exists(p_ht_path)) {
  stop("HT raster not found: ", p_ht_path)
}
if (!file.exists(p_lt_path)) {
  stop("LT raster not found: ", p_lt_path)
}

out_tif <- file.path(dir_proj_htlt, paste0("HTLTPeatland_",jobht, "_", joblt, ".tif"))
out_png <- file.path(dir_figures, paste0("HTLTPeatland_",jobht, "_", joblt, ".png"))
out_pdf <- file.path(dir_figures, paste0("HTLTPeatland_",jobht, "_", joblt, ".pdf"))

# ---------------------------- #
# Load input rasters
# ---------------------------- #

p_ht <- rast(p_ht_path)
p_lt <- rast(p_lt_path)


message("Loaded HT raster: ", p_ht_path)
message("Loaded LT raster: ", p_lt_path)
message("Loaded bio1 raster: ", bio_path)

# Rotate if needed
if (ext(p_ht)[2] > 181) {
  p_ht <- terra::rotate(p_ht)
}
if (ext(p_lt)[2] > 181) {
  p_lt <- terra::rotate(p_lt)
}
if (ext(bio1)[2] > 181) {
  bio1 <- terra::rotate(bio1)
}

# ---------------------------- #
# Temperature-gated soft blend
# ---------------------------- #

w_ht <- 1 / (1 + exp(-(bio1 - t0) / w))
w_lt <- 1 - w_ht

p_global <- (p_ht * w_ht) + (p_lt * w_lt)

writeRaster(p_global, filename = out_tif, overwrite = TRUE)
message("Saved blended raster: ", out_tif)

# ---------------------------- #
# Prepare raster for Robinson map
# ---------------------------- #

r <- p_global / 1000
crs(r) <- "EPSG:4326"
r <- crop(r, ext(-180, 180, -90, 90), snap = "in")

robin_crs <- st_crs("ESRI:54030")
r_proj <- project(r, "ESRI:54030", method = "bilinear")

dx <- res(r_proj)[1]
dy <- res(r_proj)[2]

df_r <- as.data.frame(r_proj, xy = TRUE, na.rm = TRUE)
names(df_r)[3] <- "prob"

# ---------------------------- #
# Graticules and frame
# ---------------------------- #

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

# ---------------------------- #
# Plot Robinson map
# ---------------------------- #

title_txt <- paste0("HTLTPeatland probability (",climate_id,")")

p <- ggplot() +
  geom_sf(data = grat_r, color = "grey80", linewidth = 0.35) +
  geom_tile(data = df_r, aes(x = x, y = y, fill = prob), width = dx, height = dy) +
  geom_sf(data = frame_r, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_viridis_c(
    option = "viridis",
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
  labs(title = title_txt) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14,margin = ggplot2::margin(b = 2)),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = ggplot2::margin(8, 0, 0, 0)
  )

ggsave(
  filename = out_png,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 10,
  height = 6,
  bg = "white"
)

message("Saved PNG: ", out_png)
message("Saved PDF: ", out_pdf)