# ============================================================
# Script: 06_plot_TSS_across_algos_and_PAs.R 
# Purpose:
#   Plot TSS (mean ± SD) across a series of single-model runs
#   for different algorithms and pseudo-absence settings.
#
# Inputs:
#   - biomod2 single-model output files from 04_build_single_models.R
#   - Model job IDs and PA labels based on config/config_build_single_models.csv
#
# Outputs:
#   - Figure of TSS across algorithms and PA settings (PNG and PDF)
#   - CSV file containing the plotted summary statistics
#
# Notes:
#   - This script extracts TSS values for RUN1, RUN2, and RUN3.
#   - It plots calibration and validation TSS separately.
#   - Panel selection is controlled by the `panel` setting below.
# ============================================================

suppressPackageStartupMessages({
  library(biomod2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(here)
  library(tibble)
})

here::i_am("code/06_plot_TSS_across_algos_and_PAs.R")

# ---------------------------- #
# User settings
# ---------------------------- #

myRespName <- "GlobalPeatland"
panel <- 2

if (panel == 1) {
  models <- c("ANN", "CTA", "FDA", "GAM")
  figure_base <- "06a_TSS_across_algos_PAs_panel1"
} else if (panel == 2) {
  models <- c("GBM", "GLM", "SRE", "XGBOOST")
  figure_base <- "06b_TSS_across_algos_PAs_panel2"
} else {
  stop("panel must be 1 or 2.")
}

# Job series and PA labels based on config/config_build_single_models.csv
job_suffix <- c("001", "002", "003", "004", "005", "006", "007")

pa_labels <- c(
  "10×100",
  "20×100",
  "10×1000",
  "10×5000",
  "10×10,000",
  "10×140,000",
  "3×450,000"
)

if (length(job_suffix) != length(pa_labels)) {
  stop("job_suffix and pa_labels must have the same length.")
}

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

dir_models   <- here(myRespName)
dir_figures  <- here("outputs")

create_directory(dir_figures)

out_png      <- file.path(dir_figures, paste0(figure_base, ".png"))
out_pdf      <- file.path(dir_figures, paste0(figure_base, ".pdf"))
out_csv      <- file.path(dir_figures, paste0(figure_base, "_data.csv"))

# ---------------------------- #
# Helper function
# ---------------------------- #

extract_tss_summary <- function(model_out_path) {
  obj <- get(load(model_out_path))
  
  tss <- get_evaluations(
    obj,
    metric.eval = "TSS",
    run = c("RUN1", "RUN2", "RUN3")
  )
  
  tss %>%
    summarise(
      calibration_mean = mean(calibration, na.rm = TRUE),
      calibration_sd   = sd(calibration, na.rm = TRUE),
      validation_mean  = mean(validation, na.rm = TRUE),
      validation_sd    = sd(validation, na.rm = TRUE)
    )
}

# ---------------------------- #
# Collect summary statistics
# ---------------------------- #

all_metrics <- list()

for (model in models) {
  job_ids <- paste0(model, job_suffix)
  
  model_results <- vector("list", length(job_ids))
  
  for (i in seq_along(job_ids)) {
    model_out_path <- file.path(
      dir_models,
      paste0(myRespName, ".", job_ids[i], ".models.out")
    )
    
    message("Reading: ", model_out_path)
    
    if (!file.exists(model_out_path)) {
      warning("Missing file: ", model_out_path)
      model_results[[i]] <- tibble(
        calibration_mean = NA_real_,
        calibration_sd   = NA_real_,
        validation_mean  = NA_real_,
        validation_sd    = NA_real_
      )
      next
    }
    
    model_results[[i]] <- tryCatch(
      extract_tss_summary(model_out_path),
      error = function(e) {
        warning("Failed to read TSS from ", model_out_path, ": ", e$message)
        tibble(
          calibration_mean = NA_real_,
          calibration_sd   = NA_real_,
          validation_mean  = NA_real_,
          validation_sd    = NA_real_
        )
      }
    )
  }
  
  model_df <- bind_rows(model_results) %>%
    mutate(
      Model = model,
      jobid = job_ids,
      PA = pa_labels
    )
  
  all_metrics[[model]] <- model_df
}

all_metrics_df <- bind_rows(all_metrics)

# ---------------------------- #
# Reshape for plotting
# ---------------------------- #

plot_df <- all_metrics_df %>%
  pivot_longer(
    cols = c(
      calibration_mean, validation_mean,
      calibration_sd, validation_sd
    ),
    names_to = c("metric", ".value"),
    names_pattern = "(calibration|validation)_(mean|sd)"
  ) %>%
  mutate(
    PA = factor(PA, levels = pa_labels),
    Model = factor(Model, levels = models),
    metric = factor(metric, levels = c("calibration", "validation"))
  )

write_csv(plot_df, out_csv)

# ---------------------------- #
# Plot
# ---------------------------- #

pd <- position_dodge(width = 0.35)

p_base <- ggplot(
  plot_df,
  aes(x = PA, y = mean, color = metric, group = metric)
) +
  geom_point(position = pd, size = 3) +
  geom_line(position = pd) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.2,
    position = pd,
    color = "grey30"
  ) +
  geom_text(
    data = filter(plot_df, metric == "calibration"),
    aes(label = round(mean, 2)),
    position = pd,
    vjust = -0.8,
    size = 3,
    color = "black"
  ) +
  geom_text(
    data = filter(plot_df, metric == "validation"),
    aes(label = round(mean, 2)),
    position = pd,
    vjust = 1.8,
    size = 3,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      calibration = "#FF7F00",
      validation  = "#FDBF6F"
    ),
    labels = c(
      calibration = "Calibration",
      validation  = "Validation"
    )
  ) +
  scale_x_discrete(limits = pa_labels) +
  coord_cartesian(ylim = c(0.40, 1.02)) +
  labs(
    x = "",
    y = "TSS (mean ± SD)",
    color = "Evaluation set"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.placement = "outside",
    panel.spacing = unit(0, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA)
  )

if (panel == 1){
  p <- p_base +
    facet_wrap(~ Model, ncol = 1) 
} else if (panel == 2){
  p <- p_base +
    facet_wrap(~ Model,ncol = 1,
               labeller = as_labeller(c(GBM = 'GBM', GLM = 'GLM', SRE = 'SRE', XGBOOST = "XGBoost")))
}


ggsave(
  filename = out_png,
  plot = p,
  width = 4.6,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename = out_pdf,
  plot = p,
  width = 4.6,
  height = 8,
  units = "in"
)

message("Saved PNG: ", out_png)
message("Saved PDF: ", out_pdf)
message("Saved plotted data: ", out_csv)
