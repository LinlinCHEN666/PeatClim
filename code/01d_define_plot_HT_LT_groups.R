# ============================================================
# Script: 01d_define_plot_HT_LT_groups.R
# Purpose:
#   Define HT- and LT- peatland groups from the integrated
#   PEATMAP + Peat-DBase occurrence dataset using WorldClim
#   climate variables on the 5 arc-min grid, 
#   and create figures for the resulting groups.
#
# Inputs:
#   - Integrated occurrence CSV from 01c_integrate_PEATMAP_Peat-DBase.R
#   - WorldClim v2.1 bioclimatic rasters:
#       bio1  = annual mean temperature
#       bio12 = annual precipitation
#
# Outputs:
#   - HTPeatland occurrence CSV
#   - LTPeatland occurrence CSV
#   - HTLTPeatland occurrence CSV annotated with extracted climate values
#     and group labels
#   - CSV summary of counts and classification boundary
#   - Scatter plot of bio1 vs bio12 with HT/LT groups (PNG and PDF)
#   - Global Robinson-projection map of HT/LT peatland groups (PNG and PDF)
#
# Notes:
#   - HT/LT groups are defined using k-means clustering (k = 2) on bio1.
#   - The temperature boundary is defined as the mean of the two cluster centres.
#   - Records with missing bio1 or bio12 are removed before classification.
#
# External data requirements:
#   WorldClim v2.1 bioclimatic variables
#   - Citation:
#       Fick, S. E. and Hijmans, R. J. (2017).
#       WorldClim 2: new 1 km spatial resolution climate surfaces
#       for global land areas. International Journal of Climatology,
#       37(12), 4302-4315.
#   - Data used here:
#       WorldClim v2.1 bioclimatic variables at 5 arc-min resolution
#       (bio1 and bio12 used in this script)
#   - Download from:
#       https://www.worldclim.org/data/worldclim21.html
#   - Place extracted files under:
#       data_external/WorldClim/
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(tibble)
  library(readr)
  library(here)
  library(ggplot2)
  library(scales)
})

here::i_am("code/01d_define_plot_HT_LT_groups.R")

# ---------------------------- #
# User settings
# ---------------------------- #

# Rasterization rule:
# FALSE = centre-of-cell rule
# TRUE  = any touched cell counts as presence
touches_rule <- FALSE

kmeans_seed <- 123

x_limits <- c(0, 7000)
y_limits <- c(-20, 30)

group_levels <- c("HT", "LT")
group_labels <- c("HTPeatland", "LTPeatland")
group_colors <- c(
  HTPeatland = "#E31A1C",
  LTPeatland = "#1F78B4"
)

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

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_occurrence <- here("data_processed", "occurrence")
dir_wc         <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_figures    <- here("outputs")

create_directory(dir_occurrence)
create_directory(dir_figures)

if (touches_rule == FALSE) {
  path_input        <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min.csv")
  out_annotated_csv <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated.csv")
  out_summary_csv   <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_summary.csv")
  
  out_scatter_png <- file.path(dir_figures, "01d_PMPD_HT_LT_scatter.png")
  out_map_png     <- file.path(dir_figures, "01d_PMPD_HT_LT_map.png")
  out_scatter_pdf <- file.path(dir_figures, "01d_PMPD_HT_LT_scatter.pdf")
  out_map_pdf     <- file.path(dir_figures, "01d_PMPD_HT_LT_map.pdf")
} else {
  path_input        <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_global_5min_alt.csv")
  out_annotated_csv <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated_alt.csv")
  out_summary_csv   <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_summary_alt.csv")
  
  out_scatter_png <- file.path(dir_figures, "01d_PMPD_HT_LT_scatter_alt.png")
  out_map_png     <- file.path(dir_figures, "01d_PMPD_HT_LT_map_alt.png")
  out_scatter_pdf <- file.path(dir_figures, "01d_PMPD_HT_LT_scatter_alt.pdf")
  out_map_pdf     <- file.path(dir_figures, "01d_PMPD_HT_LT_map_alt.pdf")
}

bio1_file  <- file.path(dir_wc, "wc2.1_5m_bio_1.tif")
bio12_file <- file.path(dir_wc, "wc2.1_5m_bio_12.tif")

# ---------------------------- #
# Input files
# ---------------------------- #

if (!file.exists(path_input)) {
  stop(
    "Integrated occurrence file not found.\n\n",
    "Expected location:\n",
    path_input, "\n\n",
    "Run 01c_integrate_PEATMAP_Peat-DBase.R first."
  )
}

if (!file.exists(bio1_file)) {
  stop(
    "WorldClim bio1 raster not found.\n\n",
    "Expected location:\n",
    bio1_file
  )
}

if (!file.exists(bio12_file)) {
  stop(
    "WorldClim bio12 raster not found.\n\n",
    "Expected location:\n",
    bio12_file
  )
}

occ <- read_csv(path_input, show_col_types = FALSE)
bio1 <- rast(bio1_file)[[1]]
bio12 <- rast(bio12_file)[[1]]

message("Integrated occurrence file loaded: ", path_input)
message("bio1 raster loaded: ", bio1_file)
message("bio12 raster loaded: ", bio12_file)
message("Resolution: ", paste(res(bio1), collapse = ", "))
message("CRS: ", crs(bio1))

required_cols <- c("lon", "lat", "fid", "VALUE")
missing_cols <- setdiff(required_cols, names(occ))

if (length(missing_cols) > 0) {
  stop(
    "Input occurrence file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (nlyr(bio1) != 1 || nlyr(bio12) != 1) {
  stop("Each climate raster must have exactly one layer.")
}

if (!isTRUE(all.equal(res(bio1), c(5/60, 5/60), tolerance = 1e-8))) {
  stop("bio1 raster does not have 5 arc-min resolution.")
}

if (!isTRUE(all.equal(res(bio12), c(5/60, 5/60), tolerance = 1e-8))) {
  stop("bio12 raster does not have 5 arc-min resolution.")
}

# ---------------------------- #
# Extract climate variables
# ---------------------------- #

n_input <- nrow(occ)

pts_vect <- terra::vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326")

occ$bio1 <- terra::extract(bio1, pts_vect)[, 2]
occ$bio12 <- terra::extract(bio12, pts_vect)[, 2]

occ_clim <- occ %>%
  filter(!is.na(bio1), !is.na(bio12))

if (nrow(occ_clim) == 0) {
  stop("No occurrence records remain after extracting bio1 and bio12.")
}

message("Records after removing missing bio1/bio12: ", nrow(occ_clim))

# ---------------------------- #
# Define HT/LT groups
# ---------------------------- #

set.seed(kmeans_seed)
k2 <- kmeans(occ_clim$bio1, centers = 2)

centres <- sort(as.numeric(k2$centers))
boundary <- mean(centres)

occ_clim$group <- ifelse(occ_clim$bio1 < boundary, "LT", "HT")

message("K-means centres (sorted): ", paste(round(centres, 4), collapse = ", "))
message("HT/LT boundary (bio1): ", round(boundary, 4))

# ---------------------------- #
# Prepare and write data outputs
# ---------------------------- #

annotated <- occ_clim %>%
  mutate(
    kmeans_seed = kmeans_seed,
    boundary_bio1 = boundary
  ) %>%
  select(lon, lat, fid, VALUE, bio1, bio12, group, kmeans_seed, boundary_bio1)

write_csv(annotated, out_annotated_csv)
message("Annotated HT/LT CSV saved: ", out_annotated_csv)

lt <- occ_clim %>%
  filter(group == "LT")

ht <- occ_clim %>%
  filter(group == "HT")

processing_summary <- tibble(
  metric = c(
    "input_records",
    "records_with_bio1_bio12",
    "ht_records",
    "lt_records",
    "kmeans_seed",
    "bio1_cluster_center_low",
    "bio1_cluster_center_high",
    "bio1_boundary"
  ),
  value = c(
    n_input,
    nrow(occ_clim),
    nrow(ht),
    nrow(lt),
    kmeans_seed,
    centres[1],
    centres[2],
    boundary
  )
)

write_csv(processing_summary, out_summary_csv)
message("HT/LT summary CSV saved: ", out_summary_csv)

# ---------------------------- #
# Prepare plotting data
# ---------------------------- #

dat <- annotated %>%
  mutate(
    group = factor(group, levels = group_levels, labels = group_labels)
  )

# ---------------------------- #
# 1) Scatter plot: bio12 vs bio1
# ---------------------------- #

p_scatter <- ggplot(dat) +
  geom_point(
    aes(x = bio12, y = bio1, color = group),
    shape = 20,
    alpha = 0.4,
    size = 1
  ) +
  scale_color_manual(
    values = group_colors,
    breaks = group_labels
  ) +
  guides(
    color = guide_legend(override.aes = list(alpha = 1, size = 3))
  ) +
  geom_hline(
    yintercept = boundary,
    linetype = "dashed",
    linewidth = 1,
    color = "black"
  ) +
  scale_y_continuous(
    breaks = sort(unique(c(pretty(dat$bio1), boundary))),
    labels = scales::number_format(accuracy = 0.01),
    minor_breaks = NULL
  ) +
  coord_cartesian(
    xlim = x_limits,
    ylim = y_limits
  ) +
  labs(
    x = "BIO12: Annual Precipitation (mm)",
    y = "BIO1: Annual Mean Temperature (°C)",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    legend.position = c(0.98, 0.02),
    legend.justification = c("right", "bottom"),
    legend.background = element_rect(fill = NA, colour = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = out_scatter_png,
  plot = p_scatter,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = out_scatter_pdf,
  plot = p_scatter,
  width = 8,
  height = 5,
  units = "in"
)

message("Scatter plot saved: ", out_scatter_png)
message("Scatter plot saved: ", out_scatter_pdf)

# ---------------------------- #
# 2) Robinson map
# ---------------------------- #

r <- bio1
crs(r) <- "EPSG:4326"
r <- crop(r, ext(-180, 180, -90, 90), snap = "in")

robin_crs <- st_crs("ESRI:54030")
r_proj <- project(r, "ESRI:54030", method = "bilinear")

df_r <- as.data.frame(r_proj, xy = TRUE, na.rm = TRUE)
names(df_r)[3] <- "bio1"

grat_ll <- st_graticule(
  lat = seq(-60, 60, by = 30),
  lon = seq(-180, 180, by = 60)
)
grat_r <- st_transform(st_as_sf(grat_ll), robin_crs)

frame_r <- make_robinson_frame_sf(robin_crs)
bb <- st_bbox(frame_r)

df_cat <- dat %>%
  mutate(group_id = ifelse(group == "HTPeatland", 1, 2)) %>%
  select(lon, lat, group_id)

v <- vect(df_cat, geom = c("lon", "lat"), crs = "EPSG:4326")

r_temp <- rast(
  xmin = -180,
  xmax = 180,
  ymin = -90,
  ymax = 90,
  resolution = res(bio1)[1],
  crs = "EPSG:4326"
)

r_cat <- rasterize(v, r_temp, field = "group_id")
r_robin <- project(r_cat, "ESRI:54030", method = "near")

df_robin <- as.data.frame(r_robin, xy = TRUE, na.rm = TRUE)
names(df_robin)[3] <- "group_id"

df_robin$group <- factor(
  df_robin$group_id,
  levels = c(1, 2),
  labels = c("HTPeatland", "LTPeatland")
)

p_base <- ggplot() +
  geom_sf(data = grat_r, color = "grey80", linewidth = 0.35) +
  geom_raster(data = df_r, aes(x = x, y = y), fill = "grey85") +
  geom_sf(data = frame_r, fill = NA, color = "black", linewidth = 0.9) +
  coord_sf(
    crs = robin_crs,
    datum = NA,
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"]),
    expand = FALSE
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1))
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.margin = ggplot2::margin(0, 0, 0, 0),
    legend.box.margin = ggplot2::margin(0, 0, 0, 0),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.margin = ggplot2::margin(8, 0, 0, 0)
  )

p_map <- p_base +
  geom_raster(data = df_robin, aes(x = x, y = y, fill = group)) +
  scale_fill_manual(values = group_colors, name = NULL)

ggsave(
  filename = out_map_png,
  plot = p_map,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = out_map_pdf,
  plot = p_map,
  width = 10,
  height = 6,
  units = "in"
)

message("HT/LT map saved: ", out_map_png)
message("HT/LT map saved: ", out_map_pdf)