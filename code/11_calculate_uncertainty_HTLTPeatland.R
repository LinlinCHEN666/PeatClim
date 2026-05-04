# ============================================================
# Script: 11_calculate_uncertainty_HTLTPeatland.R
# Purpose:
#   Calculate propagated uncertainty (coefficient of variation, CV)
#   for the integrated HTLTPeatland projection by propagating
#   uncertainty from the HTPeatland and LTPeatland ensemble projections.
#
# Inputs:
#   - HTPeatland ensemble projection output
#   - LTPeatland ensemble projection output
#   - WorldClim v2.1 bio1 raster at 5 arc-min resolution
# Outputs:
#   - Propagated CV raster for HTLTPeatland
#   - CSV summary of CV for all land cells and suitable cells
#   - Robinson-projection figure of propagated CV (PNG and PDF)
#
# Notes:
#   - Propagated uncertainty is calculated as:
#       sd_HTLT = sqrt((w_HT^2 * sd_HT^2) + (w_LT^2 * sd_LT^2))
#       CV_HTLT = sd_HTLT / P_HTLT
#   - Here:
#       sd_HT = CV_HT * P_HT
#       sd_LT = CV_LT * P_LT
#   - Logistic weights are derived from bio1 using split temperature t0.
#   - This script design for WorldClim 5 arc-minute projections. 
#     Change input file if use other climate scenarios.
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(biomod2)
  library(sf)
  library(ggplot2)
  library(scales)
  library(readr)
  library(tibble)
  library(here)
})

here::i_am("code/11_calculate_uncertainty_HTLTPeatland.R")

# ---------------------------- #
# User settings
# ---------------------------- #

t0 <- 11.37
w  <- 2
thres <- 0.6

myRespName    <- "HTLTPeatland"
myRespName_ht <- "HTPeatland"
myRespName_lt <- "LTPeatland"

jobht   <- "EM01_GBM013_WC5min"
joblt   <- "EM01_GBM012_WC5min"
jobhtlt <- paste0(jobht,"_",joblt)

# ---------------------------- #
# Helper functions
# ---------------------------- #

create_directory <- function(path) {
  if (length(path) != 1) stop("'path' must be a single string.")
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Failed to create directory: ", path)
}

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

get_projection_paths <- function(bm_proj_obj) {
  links <- bm_proj_obj@proj.out@link
  
  cv_path <- links[grepl("EMcv", links)]
  mean_path <- links[
    grepl("EMmean", links) &
      !grepl("bin", links) &
      !grepl("filt", links)
  ]
  
  if (length(cv_path) != 1) stop("Expected 1 EMcv file, found: ", length(cv_path))
  if (length(mean_path) != 1) stop("Expected 1 EMmean file, found: ", length(mean_path))
  
  list(cv = cv_path, mean = mean_path)
}

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_wc        <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_proj      <- here()
dir_figures   <- here("outputs")

job_dir_htlt <- file.path(dir_figures, "11_uncertainty", paste0("proj_", jobhtlt))

create_directory(job_dir_htlt)
create_directory(dir_figures)

bio1_path <- file.path(dir_wc, "wc2.1_5m_bio_1.tif")

myBiomodEMProj_ht_path <- file.path(
  dir_proj, myRespName_ht, paste0("proj_", jobht),
  paste0(myRespName_ht, ".", jobht, ".ensemble.projection.out")
)

myBiomodEMProj_lt_path <- file.path(
  dir_proj, myRespName_lt, paste0("proj_", joblt),
  paste0(myRespName_lt, ".", joblt, ".ensemble.projection.out")
)

out_tif <- file.path(job_dir_htlt, paste0(myRespName, "_propagated_cv.tif"))
out_png <- file.path(job_dir_htlt, paste0(myRespName, "_propagated_cv.png"))
out_pdf <- file.path(job_dir_htlt, paste0(myRespName, "_propagated_cv.pdf"))
out_csv <- file.path(job_dir_htlt, paste0(myRespName, "_propagated_cv_summary.csv"))

# ---------------------------- #
# Input checks
# ---------------------------- #

required_files <- c(myBiomodEMProj_ht_path, myBiomodEMProj_lt_path, bio1_path)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

# ---------------------------- #
# Load HT and LT ensemble projections
# ---------------------------- #

myBiomodEMProj_ht <- get(load(myBiomodEMProj_ht_path))
myBiomodEMProj_lt <- get(load(myBiomodEMProj_lt_path))

paths_ht <- get_projection_paths(myBiomodEMProj_ht)
paths_lt <- get_projection_paths(myBiomodEMProj_lt)

cv_ht <- rast(paths_ht$cv) / 1000
p_ht  <- rast(paths_ht$mean) / 1000

cv_lt <- rast(paths_lt$cv) / 1000
p_lt  <- rast(paths_lt$mean) / 1000

compareGeom(p_ht, p_lt, cv_ht, cv_lt, stopOnError = TRUE)

# ---------------------------- #
# Load and align bio1
# ---------------------------- #

bio1 <- rast(bio1_path)

bio1_range <- global(bio1, range, na.rm = TRUE)
if (bio1_range[1, 2] > 100) {
  message("BIO1 appears to be scaled by 10. Dividing by 10.")
  bio1 <- bio1 / 10
}

if (!compareGeom(p_ht, bio1, stopOnError = FALSE)) {
  message("BIO1 geometry does not match prediction rasters. Resampling BIO1 to prediction grid.")
  bio1 <- project(bio1, p_ht, method = "bilinear")
  bio1 <- resample(bio1, p_ht, method = "bilinear")
}

# ---------------------------- #
# Propagate uncertainty
# ---------------------------- #

# Convert CV to SD
sd_ht <- cv_ht * p_ht
sd_lt <- cv_lt * p_lt

# Logistic weights from temperature
w_ht <- 1 / (1 + exp(-(bio1 - t0) / w))
w_lt <- 1 - w_ht

# Integrated probability
p_htlt <- (p_ht * w_ht) + (p_lt * w_lt)

# Propagated SD and CV
sd_htlt <- sqrt((w_ht^2 * sd_ht^2) + (w_lt^2 * sd_lt^2))

eps <- 1e-6
cv_htlt <- ifel(p_htlt > eps, sd_htlt / p_htlt, NA)

writeRaster(cv_htlt, filename = out_tif, overwrite = TRUE)
message("Saved propagated CV raster: ", out_tif)

# ---------------------------- #
# Summarize propagated CV
# ---------------------------- #

q_global <- terra::global(
  cv_htlt,
  fun = function(x, ...) quantile(x, probs = c(0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
)

global_summary <- data.frame(
  area = "global",
  mean_cv   = as.numeric(global(cv_htlt, fun = "mean",   na.rm = TRUE)[1, 1]),
  median_cv = as.numeric(global(cv_htlt, fun = median,   na.rm = TRUE)[1, 1]),
  q25_cv    = as.numeric(q_global[1]),
  q75_cv    = as.numeric(q_global[3]),
  iqr_cv    = as.numeric(q_global[3] - q_global[1]),
  p95_cv    = as.numeric(q_global[4]),
  max_cv    = as.numeric(global(cv_htlt, fun = "max",    na.rm = TRUE)[1, 1])
)

# summarize CV for suitable area, mask using probability threshold
cv_suitable <- ifel(p_htlt >= thres, cv_htlt, NA)

q_suitable <- terra::global(
  cv_suitable,
  fun = function(x, ...) quantile(x, probs = c(0.25, 0.5, 0.75, 0.95), na.rm = TRUE)
)

suitable_summary <- data.frame(
  area = "suitable",
  mean_cv   = as.numeric(global(cv_suitable, fun = "mean",   na.rm = TRUE)[1, 1]),
  median_cv = as.numeric(global(cv_suitable, fun = median,   na.rm = TRUE)[1, 1]),
  q25_cv    = as.numeric(q_suitable[1]),
  q75_cv    = as.numeric(q_suitable[3]),
  iqr_cv    = as.numeric(q_suitable[3] - q_suitable[1]),
  p95_cv    = as.numeric(q_suitable[4]),
  max_cv    = as.numeric(global(cv_suitable, fun = "max",    na.rm = TRUE)[1, 1])
)

cv_summary <- rbind(global_summary, suitable_summary)
write_csv(cv_summary, out_csv)
message("Saved CV summary: ", out_csv)

# ---------------------------- #
# Plot Robinson map
# ---------------------------- #

r <- cv_htlt
crs(r) <- "EPSG:4326"
r <- crop(r, ext(-180, 180, -90, 90), snap = "in")

robin_crs <- st_crs("ESRI:54030")
r_proj <- project(r, "ESRI:54030", method = "bilinear")

dx <- res(r_proj)[1]
dy <- res(r_proj)[2]

df_r <- as.data.frame(r_proj, xy = TRUE, na.rm = TRUE)
names(df_r)[3] <- "cv"

grat_ll <- st_graticule(
  lat = seq(-60, 60, by = 30),
  lon = seq(-180, 180, by = 60)
)
grat_r <- st_transform(st_as_sf(grat_ll), robin_crs)

frame_r <- make_robinson_frame_sf(robin_crs)
bb <- st_bbox(frame_r)

p <- ggplot() +
  geom_sf(data = grat_r, color = "grey80", linewidth = 0.35) +
  geom_tile(data = df_r, aes(x = x, y = y, fill = cv), width = dx, height = dy) +
  geom_sf(data = frame_r, fill = NA, color = "black", linewidth = 0.9) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(0, 0.03),
    oob = squish,
    na.value = "white",
    name = "CV",
    breaks = seq(0, 0.03, 0.01),
    labels = sprintf("%.2f", seq(0, 0.03, 0.01))
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
  labs(title = "HTLTPeatland propagated CV") +
  theme_void(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14,
      margin = ggplot2::margin(b = 8)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = ggplot2::margin(15, 0, 0, 0)
  )

ggsave(out_png, plot = p, width = 10, height = 6, dpi = 300, bg = "white")
ggsave(out_pdf, plot = p, width = 10, height = 6, bg = "white")

message("Saved PNG: ", out_png)
message("Saved PDF: ", out_pdf)