# ============================================================
# Script: 08_build_ensemble_models.R
# Purpose:
#   Build biomod2 ensemble models from previously trained single-model outputs.
#
# Inputs:
#   - Config CSV file
#   - Single-model biomod2 output from 04_build_single_models.R
#
# Outputs:
#   - biomod2 ensemble-model output object
#   - Ensemble model evaluation plots (PNG and PDF)
#
# Notes:
#   - This script is intended to be run in HPCby command line:
#       Rscript code/08_build_ensemble_models.R config/config_build_ensemble_models.csv <jobid>
#   - Production runs may be submitted through SLURM using:
#       hpc/submit_r_job.sh
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(biomod2)
  library(tibble)
  library(here)
})

here::i_am("code/08_build_ensemble_models.R")

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

# ---------------------------- #
# Extract parameters
# ---------------------------- #

myRespName <- params$myRespName[[1]]
jobid_single <- params$jobid_single[[1]]
em_by <- params$em_by[[1]]
metric_select <- params$metric_select[[1]]
metric_select_thresh <- params$metric_select_thresh[[1]]
nb_cpu <- params$nb_cpu[[1]]

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_models    <- here()
single_model_path <- file.path(dir_models, myRespName, paste0(myRespName,".",jobid_single, ".models.out"))

dir_figures   <- here("outputs", "08_diagnostics")
job_fig_dir      <- file.path(dir_figures, myRespName, jobid)

create_directory(dir_figures)
create_directory(job_fig_dir)

# ---------------------------- #
# Read single-model output
# ---------------------------- #

if (!file.exists(single_model_path)) {
  stop(
    "Single-model output not found: ", single_model_path, "\n",
    "Run 04_build_single_models.R first for jobid_single = ", jobid_single
  )
}

myBiomodModelOut <- get(load(single_model_path))

# ---------------------------- #
# Build ensemble model
# ---------------------------- #

selected_models <- get_built_models(myBiomodModelOut)[
  !grepl("allData|allRun", get_built_models(myBiomodModelOut))
]

message("Selected models for ensemble: ", paste(selected_models, collapse = ", "))

message("Build ensemble models:...")

myEnsembleModelOut <- BIOMOD_EnsembleModeling(
  bm.mod = myBiomodModelOut,
  models.chosen = selected_models ,
  em.by = em_by,
  em.algo = c("EMmean", "EMcv"),
  metric.select = c(metric_select),
  metric.select.thresh = c(metric_select_thresh),
  metric.eval = c("TSS", "ROC"),
  var.import = 10,
  EMci.alpha = 0.05,
  seed.val = 42,
  nb.cpu = nb_cpu
)

print("myEnsembleModelOut")
print(myEnsembleModelOut)

# ---------------------------- #
# Evaluation plots
# ---------------------------- #

p_eval <- bm_PlotEvalBoxplot(bm.out = myEnsembleModelOut, group.by = c("algo", "algo"))

png(
  file.path(job_fig_dir, paste0("EvalBoxplot_EM_", myRespName, "_", jobid, ".png")),
  width = 8, height = 4, units = "in", res = 300
)
print(p_eval)
dev.off()

pdf(
  file.path(job_fig_dir, paste0("EvalBoxplot_EM_", myRespName, "_", jobid, ".pdf")),
  width = 8, height = 4
)
print(p_eval)
dev.off()

message("Ensemble-model building completed successfully.")