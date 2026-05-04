# ============================================================
# Script: 04_calculate_bioclim_correlations.R
# Purpose:
#   Calculate Pearson correlations among the 19 WorldClim v2.1
#   bioclimatic variables at 5 arc-min resolution and create
#   a correlation heatmap.
#
# Inputs:
#   - WorldClim v2.1 bioclimatic rasters (bio1 to bio19)
#
# Outputs:
#   - 19 × 19 Pearson correlation matrix (CSV)
#   - heatmaps of the correlation matrix (PNG and PDF)
#
# Notes:
#   - Correlations are calculated using non-missing raster cells only,
#     and therefore represent land areas covered by WorldClim.
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
#   - Download from:
#       https://www.worldclim.org/data/worldclim21.html
#   - Place extracted files under:
#       data_external/WorldClim/
# ============================================================

suppressPackageStartupMessages({
  library(terra)
  library(reshape2)
  library(ggplot2)
  library(readr)
  library(here)
})

here::i_am("code/03_calculate_bioclim_correlations.R")

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_wc        <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_figures   <- here("outputs")

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

bio_files <- file.path(dir_wc, paste0("wc2.1_5m_bio_", 1:19, ".tif"))

missing_files <- bio_files[!file.exists(bio_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing WorldClim files:\n",
    paste(missing_files, collapse = "\n")
  )
}

out_csv <- file.path(dir_figures, "03_bioclim_correlation_matrix.csv")
out_png <- file.path(dir_figures, "03_bioclim_correlation_heatmap.png")
out_pdf <- file.path(dir_figures, "03_bioclim_correlation_heatmap.pdf")

# ---------------------------- #
# Load bioclimatic variables
# ---------------------------- #

climate_stack <- rast(bio_files)
names(climate_stack) <- paste0("bio", 1:19)

message("Loaded bioclimatic variables: ", paste(names(climate_stack), collapse = ", "))

# ---------------------------- #
# Extract raster values
# ---------------------------- #
# Convert all raster layers to a data frame

climate_data <- as.data.frame(climate_stack, xy = FALSE, na.rm = TRUE)

if (nrow(climate_data) == 0) {
  stop("No complete raster cells available for correlation calculation.")
}

message("Number of raster cells used: ", nrow(climate_data))

# ---------------------------- #
# Compute Pearson correlation matrix
# ---------------------------- #

correlation_matrix <- cor(climate_data, method = "pearson")

write_csv(
  tibble::rownames_to_column(
    as.data.frame(correlation_matrix),
    var = "variable"
  ),
  out_csv
)
message("Saved correlation matrix: ", out_csv)

# ---------------------------- #
# plot correlation heatmap
# ---------------------------- #

cor_long <- melt(correlation_matrix)

p <- ggplot(cor_long, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson\ncorrelation"
  ) +
  labs(
    title = "Correlation heatmap of 19 bioclimatic variables",
    x = "Bioclimatic variable",
    y = "Bioclimatic variable"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

ggsave(out_png, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(out_pdf, plot = p, width = 10, height = 8, bg = "white")

message("Saved heatmap PNG: ", out_png)
message("Saved heatmap PDF: ", out_pdf)
