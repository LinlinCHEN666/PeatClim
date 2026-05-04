# ============================================================
# Script: 07_plot_TSS_across_predictor_sets.R
# Purpose:
#   Plot TSS (mean ± SD) for GBM single-model runs across
#   peatland groups and predictor sets.
#
# Inputs:
#   - Single-model biomod2 output files from 04_build_single_models.R
#   - Model job IDs and PA labels based on config/config_build_single_models.csv
#
# Outputs:
#   - Figure of TSS across models and predictor sets (PNG and PDF)
#   - CSV file containing the plotted summary statistics
#
# Notes:
#   - It compares GBM models trained with:
#       1. All 19 BIO predictors
#       2. Selected BIO predictors
#   - TSS is summarized across RUN1-RUN3 and PA1-PA10.
# ============================================================

suppressPackageStartupMessages({
  library(biomod2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(here)
  library(tibble)
})

here::i_am("code/07_plot_TSS_across_predictor_sets.R")

# ---------------------------- #
# User settings
# ---------------------------- #

algo <- "GBM"
runs <- paste0("RUN", 1:3)
PAs  <- paste0("PA", 1:3)

# Model metadata 
model_table <- tibble(
  Category = c(
    "GlobalPeatland", "GlobalPeatland",
    "LTPeatland",     "LTPeatland",
    "HTPeatland",     "HTPeatland"
  ),
  Series = c(
    "19bios", "select_bios",
    "19bios", "select_bios",
    "19bios", "select_bios"
  ),
  jobid = c(
    "GBM008", "GBM011",
    "GBM009", "GBM012",
    "GBM010", "GBM013"
  )
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

# ---------------------------- #
# Project paths
# ---------------------------- #
dir_models   <- here()
dir_figures  <- here("outputs")

create_directory(dir_figures)

out_png <- file.path(dir_figures, "07_TSS_across_predictors.png")
out_pdf <- file.path(dir_figures, "07_TSS_across_predictors.pdf")
out_csv <- file.path(dir_figures, "07_TSS_across_predictors_data.csv")

# ---------------------------- #
# Helper functions
# ---------------------------- #

get_single_model_out <- function(base_dir, myRespName, jobid) { 
  file_out <- file.path(
    base_dir,
    myRespName,
    paste0(myRespName, ".", jobid, ".models.out")
  )
  
  if (!file.exists(file_out)) {
    stop("Could not find model output file: ", file_out)
  }
  
  obj <- get(load(file_out))
  message("Loaded: ", file_out)
  message("Failed models: ", paste(obj@models.failed, collapse = ", "))
  
  obj
}

get_tss_summary <- function(model_out, algo = "GBM",
                            runs = paste0("RUN", 1:3),
                            PAs = paste0("PA", 1:10)) {
  eval_df <- get_evaluations(
    obj = model_out,
    algo = algo,
    metric.eval = "TSS",
    run = runs,
    PA = PAs
  )
  
  tibble(
    calibration_mean = mean(eval_df$calibration, na.rm = TRUE),
    calibration_sd   = sd(eval_df$calibration, na.rm = TRUE),
    validation_mean  = mean(eval_df$validation, na.rm = TRUE),
    validation_sd    = sd(eval_df$validation, na.rm = TRUE)
  )
}

# ---------------------------- #
# Load model outputs and summarize TSS
# ---------------------------- #

summary_list <- vector("list", nrow(model_table))

for (i in seq_len(nrow(model_table))) {
  row_i <- model_table[i, ]
  
  model_out <- get_single_model_out(
    base_dir = dir_models,
    myRespName = row_i$Category,
    jobid = row_i$jobid
  )
  
  summary_list[[i]] <- get_tss_summary(
    model_out = model_out,
    algo = algo,
    runs = runs,
    PAs = PAs
  ) %>%
    mutate(
      Category = row_i$Category,
      Series = row_i$Series,
      jobid = row_i$jobid
    )
}

sum_df <- bind_rows(summary_list)

# ---------------------------- #
# Reshape for plotting
# ---------------------------- #

plot_df <- sum_df %>%
  pivot_longer(
    cols = c(
      calibration_mean, validation_mean,
      calibration_sd, validation_sd
    ),
    names_to = c("Dataset", ".value"),
    names_pattern = "(calibration|validation)_(mean|sd)"
  ) %>%
  mutate(
    Category = factor(
      Category,
      levels = c("GlobalPeatland", "LTPeatland", "HTPeatland")
    ),
    Series = factor(
      Series,
      levels = c("19bios", "select_bios")
    ),
    Dataset = factor(
      Dataset,
      levels = c("calibration", "validation")
    )
  )

write_csv(plot_df, out_csv)

# ---------------------------- #
# Plot
# ---------------------------- #

pd <- position_dodge(width = 0.35)

p <- ggplot(plot_df, aes(x = Category, y = mean)) +
  geom_point(
    aes(
      shape = Series,
      color = Category,
      alpha = Dataset,
      group = interaction(Series, Dataset)
    ),
    position = pd,
    size = 4
  ) +
  geom_errorbar(
    aes(
      ymin = mean - sd,
      ymax = mean + sd,
      group = interaction(Series, Dataset)
    ),
    position = pd,
    width = 0.12,
    linewidth = 0.65,
    color = "grey30"
  ) +
  scale_shape_manual(
    values = c("19bios" = 16, "select_bios" = 17),
    labels = c(
      "19bios" = "All 19 BIOs",
      "select_bios"  = "Selected BIOs"
    )
  ) +
  scale_color_manual(
    values = c(
      "GlobalPeatland" = "#FF7F00",
      "LTPeatland"     = "#1F78B4",
      "HTPeatland"     = "#E31A1C"
    ),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c(
      calibration = 1,
      validation  = 0.6
    ),
    labels = c(
      calibration = "Calibration",
      validation  = "Validation"
    )
  ) +
  labs(
    x = NULL,
    y = "TSS (mean ± SD)",
    alpha = "Evaluation set",
    shape = "Predictor set"
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.4),
    panel.grid.minor.y = element_blank(),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    axis.title   = element_text(size = 13),
    axis.text    = element_text(size = 12),
    legend.position = c(0.98, 0.01),
    legend.justification = c(1, 0),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = NA, color = NA)
  )

ggsave(
  filename = out_png,
  plot = p,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 6,
  height = 6,
  units = "in"
)

message("Saved PNG: ", out_png)
message("Saved PDF: ", out_pdf)
message("Saved source data: ", out_csv)