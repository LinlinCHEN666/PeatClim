# ============================================================
# Script: 04_build_single_models.R
# Purpose:
#   Format biomod2 input data and build single SDM models
#   from prepared peatland occurrence datasets 
#   using WorldClim v2.1 bioclimatic predictors.
#
# Inputs:
#   - Config CSV file
#   - Prepared occurrence CSV selected by occur_id
#   - WorldClim v2.1 5 arc-min bioclimatic rasters
#
# Outputs:
#   - biomod2 single-model output object
#   - Model evaluation plots (PNG and PDF)
#   - Model diagnostic plots: Variable importance and response-curve (PNG and PDF)
#
# Notes:
#   - This script is intended to be run in HPC by command line:
#     Rscript code/04_build_single_models.R config/config_build_single_models.csv <jobid>
#     <jobid> should match a row in the config CSV to specify model settings and input data, such as: -j GBM001.
#   - Production runs may be submitted through SLURM using: 
#     hpc/submit_r_job.sh   
#   - Pseudo-absence generation and model fitting are combined in this script for a single job configuration.
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
#
#   Occurrence datasets
#   - Prepared by:
#       01a_prepare_PEATMAP_occurrence.R
#       01b_prepare_Peat-DBase_occurrence.R
#       01c_integrate_PEATMAP_Peat-DBase.R
#       02_define_plot_HT_LT_groups.R
# ============================================================

suppressPackageStartupMessages({
  library(biomod2)
  library(terra)
  library(dplyr)
  library(readr)
  library(tibble)
  library(here)
  library(ggplot2)
})

here::i_am("code/04_build_single_models.R")

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
# Project paths
# ---------------------------- #

dir_occurrence  <- here("data_processed", "occurrence")
dir_wc          <- here("data_external", "WorldClim", "wc2.1_5m_bio")
dir_figures     <- here("outputs","04_diagnostics")

create_directory(dir_figures)

config_path <- here(config_file)

if (!file.exists(config_path)) {
  stop("Config file not found: ", config_path)
}

# ---------------------------- #
# Read config
# ---------------------------- #

config <- read_csv(config_path, comment = "#", show_col_types = FALSE)

params <- config %>%
  filter(jobid == !!jobid)

if (nrow(params) == 0) {
  stop("No settings found for jobid: ", jobid)
}

if (nrow(params) != 1) {
  stop("Expected exactly 1 row for jobid ", jobid, ", found ", nrow(params))
}

myRespName      <- params$myRespName[[1]]
occur_id        <- params$occur_id[[1]]
nb_cpu          <- params$`nb.cpu`[[1]]
env_variables <- eval(parse(text = params$env.variables)) # e.g. c(1:19), c(1,12,15) # No quote
PA_nb_rep       <- params$`PA.nb.rep`[[1]]
PA_nb_absences  <- params$`PA.nb.absences`[[1]]
PA_strategy     <- params$`PA.strategy`[[1]]
CV_strategy     <- params$`CV.strategy`[[1]]
CV_nb_rep       <- params$`CV.nb.rep`[[1]]
CV_perc         <- params$`CV.perc`[[1]]
OPT_strategy    <- params$`OPT.strategy`[[1]]
model_value       <- params$model[[1]] # e.g. "all", c('GBM'),"c('ANN', 'GBM')" # Be care of the quote
models <- if (model_value == "all") {
  "all"
}else {
  eval(parse(text = model_value)) 
}

message("jobid: ", jobid)
message("myRespName: ", myRespName)
message("occur_id: ", occur_id)
message("env_variables: ", paste(env_variables, collapse = ", "))
message("models: ", paste(models, collapse = ", "))

# ---------------------------- #
# Define and load occurrence inputs
# ---------------------------- #
if (occur_id == 'PEATMAP_global'){
  occ_file <- file.path(dir_occurrence, "PEATMAP_global_5min.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
} else if (occur_id == 'PMPD_global') {
  occ_file <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
  occ <- occ %>% select('lon', 'lat', 'fid', 'VALUE')
} else if (occur_id == 'PMPD_HT') {
  occ_file <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
  occ <- occ %>% filter(group == 'HT') %>% select('lon', 'lat', 'fid', 'VALUE')
} else if (occur_id == 'PMPD_LT') {
  occ_file <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
  occ <- occ %>% filter(group == 'LT') %>% select('lon', 'lat', 'fid', 'VALUE')
} else if (occur_id == 'PMPD_global_alt') {
  occ_file <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated_alt.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
  occ <- occ %>% select('lon', 'lat', 'fid', 'VALUE')
} else if (occur_id == 'PMPD_HT_alt') {
  occ_file <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated_alt.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
  occ <- occ %>% filter(group == 'HT') %>% select('lon', 'lat', 'fid', 'VALUE')
} else if (occur_id == 'PMPD_LT_alt') {
  occ_file <- file.path(dir_occurrence, "PEATMAP_Peat-DBase_HTLT_5min_annotated_alt.csv")
  occ <- read_csv(occ_file, show_col_types = FALSE)
  occ <- occ %>% filter(group == 'LT') %>% select('lon', 'lat', 'fid', 'VALUE')
} else {
  stop("Invalid occur_id: ", occur_id)
}

if (!file.exists(occ_file)) {
  stop("Occurrence file not found: ", occ_file)
}

required_occ_cols <- c("lon", "lat", "fid", "VALUE")
missing_occ_cols <- setdiff(required_occ_cols, names(occ))

if (length(missing_occ_cols) > 0) {
  stop(
    "Occurrence file is missing required columns: ",
    paste(missing_occ_cols, collapse = ", ")
  )
}

myRespXY <- occ %>% select(lon, lat)
myResp   <- as.numeric(occ$VALUE)

message("Occurrence file loaded: ", occ_file)
message("Number of occurrence records: ", nrow(occ))


# ---------------------------- #
# Define output paths
# ---------------------------- #

job_fig_dir    <- file.path(dir_figures, myRespName, jobid)

create_directory(job_fig_dir)

# ---------------------------- #
# Load environmental predictors
# ---------------------------- #

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

message("Environmental predictors loaded: ", paste(names(myExpl), collapse = ", "))

# ---------------------------- #
# Format biomod2 input
# ---------------------------- #

message("Formatting biomod2 data with pseudo-absences")

myResp.PA <- ifelse(myResp == 1, 1, NA)

myBiomodData.r <- BIOMOD_FormatingData(
  resp.var       = myResp.PA,
  expl.var       = myExpl,
  resp.xy        = myRespXY,
  resp.name      = myRespName,
  PA.nb.rep      = PA_nb_rep,
  PA.nb.absences = PA_nb_absences,
  PA.strategy    = PA_strategy,
  filter.raster  = TRUE
)


# ---------------------------- #
# Build single models
# ---------------------------- #

message("Building single biomod2 models")

myBiomodModelOut <- BIOMOD_Modeling(
  bm.format    = myBiomodData.r,
  modeling.id  = jobid,
  models       = models,
  CV.strategy  = CV_strategy,
  CV.nb.rep    = CV_nb_rep,
  CV.perc      = CV_perc,
  OPT.strategy = OPT_strategy,
  var.import   = 10,
  metric.eval  = c("TSS", "ROC"),
  seed.val     = 42,
  nb.cpu       = nb_cpu
)

message("Single-model training completed. Now plot evaluation and diagnostic figures...")

# ---------------------------- #
# plots
# ---------------------------- #

# ---------------------------- #
message("Saving evaluation plots")
# ---------------------------- #

p_eval <- bm_PlotEvalBoxplot(
  bm.out = myBiomodModelOut,
  group.by = c("algo", "algo"),
  main = paste0(myRespName, " ", jobid)
)

png(file.path(job_fig_dir, paste0("EvalBoxplot_algo_algo_", myRespName, "_", jobid, ".png")),
    width = 8, height = 4, units = "in", res = 300)
print(p_eval)
dev.off()

pdf(file.path(job_fig_dir, paste0("EvalBoxplot_algo_algo_", myRespName, "_", jobid, ".pdf")),
    width = 8, height = 4)
print(p_eval)
dev.off()

# ---------------------------- #
message("Saving variable importance plots")
# ---------------------------- #

p0 <- bm_PlotVarImpBoxplot(
  bm.out = myBiomodModelOut, 
  group.by = c('expl.var', 'algo', 'algo'),
  
  do.plot  = FALSE)

df <- p0$plot$data
colnames(df)

df <- df %>% filter(run %in% c("RUN1", "RUN2", "RUN3"))

# calculate mean and sd for each variable
df_summary <- df %>%
  group_by(expl.var) %>%
  summarise(mean_var_imp = mean(.data[["var.imp"]], na.rm = TRUE),
            sd_var_imp = sd(.data[["var.imp"]], na.rm = TRUE))

write.csv(df_summary, file.path(job_fig_dir,paste0("VarImpSummary_",myRespName,"_",jobid,".csv")), row.names = FALSE)

# Reorder bios on the x axis from bio1 … bio19
bio_levels <- paste0("bio", 1:19)
df$expl.var <- factor(df$expl.var, levels = bio_levels)

# Decide what counts as “too far away”
#    Use a robust upper limit based on the 99th percentile of all values.
#    (Change 0.99 to 0.95 if you want to be stricter.)
upper_cap <- quantile(df[["var.imp"]], 0.99, na.rm = TRUE) 

# Drop the extreme outliers for plotting only, and count them
df_in <- df %>% filter(.data[["var.imp"]] <= upper_cap)
n_dropped <- nrow(df) - nrow(df_in)

# get range of ylim
# Find the variable with the highest kept value
bio_max <- df_in %>% slice_max(.data[["var.imp"]], n = 1) %>% pull(expl.var)

# Get its standard deviation
sd_max <- df_in %>%
  filter(expl.var == bio_max) %>%
  summarise(sd_val = sd(.data[["var.imp"]], na.rm = TRUE)) %>%
  pull(sd_val)

ymax <- max(df_in[["var.imp"]], na.rm = TRUE)

# define color
if (myRespName == 'GlobalPeatland'){
  fill_color <- "#FF7F00"
} else if(myRespName == 'LTPeatland'){
  fill_color <- "#1F78B4"
} else if(myRespName == 'HTPeatland'){
  fill_color <- "#E31A1C"
} else {
  fill_color <- "#FF7F00" # default color
}

# Redraw the plot with:
p <- ggplot(df, aes(x = expl.var, y = .data[["var.imp"]], group = interaction(expl.var, algo)#, 
                    #fill = algo
)) +
  geom_boxplot(width = 0.6,
               linewidth = 0.4, # for margin
               fill = fill_color, # for box
               color = "black" # for margin
  ) + #outlier.shape = NA,
  # geom_jitter(width = 0.12, size = 1.2, alpha = 0.5) + # show raw data points
  scale_x_discrete(drop = FALSE,
                   labels = function(x) sub("^bio", "BIO", x)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  coord_cartesian(ylim = c(0, ymax + sd_max + 0.02)) + #ylim = c(0, upper_cap) #ylim = c(0, ymax + sd_max + 0.02)
  labs(x = NULL, y = NULL, title = myRespName) +
  
  theme_bw(base_size = 12) +                      # white background
  theme(
    text = element_text(family = "Arial"), 
    legend.position = "none",#"right",
    
    #axis.title = element_text(size = 13), # face = "bold"
    axis.text.y = element_text(size =12),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    
    panel.grid = element_blank(),                 # cleaner look
    
    plot.title = element_text(size = 13, face ='bold', hjust = 0.5), # set to middle
    
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = NA,color = NA)
    
  )

# save plot
ggsave(
  filename = file.path(job_fig_dir,paste0("VarImpBoxplot_",myRespName,"_",jobid,".png")),
  plot = p,                                         
  width = 6.2,                                                 
  height = 3,
  units = "in",                                                
  dpi = 300                                                   
)
ggsave(filename = file.path(job_fig_dir,paste0("VarImpBoxplot_",myRespName,"_",jobid,".pdf")), 
       plot=p, width = 6.2, height = 3, device = cairo_pdf)


# ---------------------------- #
message("Saving response curves plots: max")
# ---------------------------- #

pred <- get_built_models(myBiomodModelOut)
for (i in 1:length(models)){
  
  pred_subset <- pred[ grepl(models[i], pred) & !grepl("allRun", pred) ]
  print(pred_subset)
  
  get_formal_data(myBiomodModelOut, "expl.var.names")

  vars_max <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                                    models.chosen = pred_subset,
                                    show.variables = get_formal_data(myBiomodModelOut, "expl.var.names"),
                                    # do.bivariate = TRUE,
                                    do.plot = FALSE,
                                    do.progress = TRUE,
                                    main = paste0('Response curves for ',myRespName,' (max)\n',jobid),
                                    fixed.var = 'max')  
  png(file.path(job_fig_dir, paste0("ResponseCurves_max_", myRespName, "_", jobid, ".png")),
      width = 12, height = 8, units = "in", res = 300)
  print(vars_max)
  dev.off()
  
  pdf(file.path(job_fig_dir, paste0("ResponseCurves_max_", myRespName, "_", jobid, ".pdf")),
      width = 12, height = 8)
  print(vars_max)
  dev.off()
  
}

# ---------------------------- #
message("Saving response curves plots: mean")
# ---------------------------- #

pred <- get_built_models(myBiomodModelOut)
for (i in 1:length(models)){
  
  pred_subset <- pred[ grepl(models[i], pred) & !grepl("allRun", pred) ]
  print(pred_subset)
  
  get_formal_data(myBiomodModelOut, "expl.var.names")

  vars_mean <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                                     models.chosen = pred_subset, 
                                     show.variables = get_formal_data(myBiomodModelOut, "expl.var.names"),
                                     # do.bivariate = TRUE,
                                     do.plot = FALSE,
                                     do.progress = TRUE,
                                     main = paste0('Response curves for ',myRespName,' (mean)\n',jobid),
                                     fixed.var = 'mean')  
  png(file.path(job_fig_dir, paste0("ResponseCurves_mean_", myRespName, "_", jobid, ".png")),
      width = 12, height = 8, units = "in", res = 300)
  print(vars_mean)
  dev.off()
  
  pdf(file.path(job_fig_dir, paste0("ResponseCurves_mean_", myRespName, "_", jobid, ".pdf")),
      width = 12, height = 8)
  print(vars_mean)
  dev.off()
}

# ---------------------------- #
message("Saving response curves plots: min")
# ---------------------------- #
pred <- get_built_models(myBiomodModelOut)
for (i in 1:length(models)){
  
  pred_subset <- pred[ grepl(models[i], pred) & !grepl("allRun", pred) ]
  print(pred_subset)
  
  get_formal_data(myBiomodModelOut, "expl.var.names")

  vars_min <- bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                                    models.chosen = pred_subset, 
                                    show.variables = get_formal_data(myBiomodModelOut, "expl.var.names"),
                                    # do.bivariate = TRUE,
                                    do.plot = FALSE,
                                    do.progress = TRUE,
                                    main = paste0('Response curves for ',myRespName,' (min)\n',jobid),
                                    fixed.var = 'min')  
  png(file.path(job_fig_dir, paste0("ResponseCurves_min_", myRespName, "_", jobid, ".png")),
      width = 12, height = 8, units = "in", res = 300)
  print(vars_min)
  dev.off()
  
  pdf(file.path(job_fig_dir, paste0("ResponseCurves_min_", myRespName, "_", jobid, ".pdf")),
      width = 12, height = 8)
  print(vars_min)
  dev.off()
}

# ---------------------------- #
message("Saving response curves plots: maxmeanmin")
# ---------------------------- #
algo <- models # e.g. "GBM", "ANN", "RF", "GLM", "GAM"

# convert to dataframe vars --------------------------------------------- #
vars_max_df <- as.data.frame(vars_max[["tab"]]) # obs = nb.pts * nvars * nruns. default nb.pts = 100
vars_mean_df <- as.data.frame(vars_mean[["tab"]]) # nb.pts is the number of x-points per curve. 
vars_min_df <- as.data.frame(vars_min[["tab"]]) # 100 points means 100 evenly spaced values between the variable min and max

# Variable names
vars <- levels(vars_max_df[["expl.name"]])
if (is.null(vars) || length(vars) == 0) {
  vars <- unique(as.character(vars_max_df[["expl.name"]]))
}

n_vars <- length(vars)

if (n_vars == 0) {
  stop("No variables found in `vars`.")
}

# ------------------------------------------------------------
# Basic checks
# ------------------------------------------------------------
required_cols <- c("expl.name", "expl.val", "pred.val")

for (obj_name in c("vars_max_df", "vars_mean_df", "vars_min_df")) {
  obj <- get(obj_name)
  miss <- setdiff(required_cols, names(obj))
  if (length(miss) > 0) {
    stop(sprintf(
      "%s is missing required columns: %s",
      obj_name,
      paste(miss, collapse = ", ")
    ))
  }
}

# If this script errors, device still gets closed
on.exit(dev.off(), add = TRUE)


# ------------------------------------------------------------
# Helper: summarise one variable from one data frame
# ------------------------------------------------------------
summarise_curve <- function(df, var_name) {
  d <- df[df$expl.name == var_name, c("expl.val", "pred.val")]
  
  if (nrow(d) == 0) {
    return(NULL)
  }
  
  mean_df <- aggregate(pred.val ~ expl.val, d, mean)
  min_df  <- aggregate(pred.val ~ expl.val, d, min)
  max_df  <- aggregate(pred.val ~ expl.val, d, max)
  
  o <- order(mean_df$expl.val)
  
  data.frame(
    x    = mean_df$expl.val[o],
    mean = mean_df$pred.val[o],
    ymin = min_df$pred.val[match(mean_df$expl.val[o], min_df$expl.val)],
    ymax = max_df$pred.val[match(mean_df$expl.val[o], max_df$expl.val)]
  )
}

# ------------------------------------------------------------
# Helper: plot one panel
# ------------------------------------------------------------
plot_response_curve_panel <- function(
    var_name,
    vars_max_df,
    vars_mean_df,
    vars_min_df,
    y_lim = c(0, 1),
    alpha_val = 0.4
) {
  max_s  <- summarise_curve(vars_max_df,  var_name)
  mean_s <- summarise_curve(vars_mean_df, var_name)
  min_s  <- summarise_curve(vars_min_df,  var_name)
  
  # If no data at all for this variable, draw an empty panel
  if (is.null(max_s) && is.null(mean_s) && is.null(min_s)) {
    plot.new()
    title(main = toupper(var_name), line = 0.5)
    text(0.5, 0.5, "No data", cex = 1)
    return(invisible(NULL))
  }
  
  # Collect x-range from available summaries
  x_vals <- c(
    if (!is.null(max_s))  max_s$x  else numeric(0),
    if (!is.null(mean_s)) mean_s$x else numeric(0),
    if (!is.null(min_s))  min_s$x  else numeric(0)
  )
  
  x_lim <- range(x_vals, na.rm = TRUE)
  
  plot(
    NA,
    xlim = x_lim,
    ylim = y_lim,
    xlab = "",
    ylab = "",
    main = "",
    yaxt = "s",
    xaxt = "s"
  )
  
  title(main = toupper(var_name), line = 0.5, cex.main = 1.2)
  
  # max: red
  if (!is.null(max_s)) {
    polygon(
      c(max_s$x, rev(max_s$x)),
      c(max_s$ymin, rev(max_s$ymax)),
      border = NA,
      col = adjustcolor("red", alpha.f = alpha_val)
    )
    lines(max_s$x, max_s$mean, col = "red", lwd = 0.6)
  }
  
  # mean: black
  if (!is.null(mean_s)) {
    polygon(
      c(mean_s$x, rev(mean_s$x)),
      c(mean_s$ymin, rev(mean_s$ymax)),
      border = NA,
      col = adjustcolor("black", alpha.f = alpha_val * 0.8)
    )
    lines(mean_s$x, mean_s$mean, col = "black", lwd = 0.6)
  }
  
  # min: blue
  if (!is.null(min_s)) {
    polygon(
      c(min_s$x, rev(min_s$x)),
      c(min_s$ymin, rev(min_s$ymax)),
      border = NA,
      col = adjustcolor("blue", alpha.f = alpha_val)
    )
    lines(min_s$x, min_s$mean, col = "blue", lwd = 0.6)
  }
}

# ------------------------------------------------------------
# Helper: choose automatic panel layout for n_vars != 19
# Returns nrow and ncol
# ------------------------------------------------------------
get_auto_layout <- function(n_vars) {
  n_col <- ceiling(sqrt(n_vars))
  n_row <- ceiling(n_vars / n_col)
  list(nrow = n_row, ncol = n_col)
}

# ------------------------------------------------------------
# Helper: draw one-line legend at the bottom
# ------------------------------------------------------------
plot_horizontal_legend <- function(alpha_val = 0.4, cex_text = 1.1) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  x_centers <- c(0.18, 0.50, 0.82)
  labels    <- c("max\u00B1sd", "mean\u00B1sd", "min\u00B1sd")
  fills     <- c(
    adjustcolor("red",   alpha.f = alpha_val),
    adjustcolor("black", alpha.f = alpha_val * 0.8),
    adjustcolor("blue",  alpha.f = alpha_val)
  )
  line_cols <- c("red", "black", "blue")
  
  box_w  <- 0.10
  box_h  <- 0.22
  y_box  <- 0.62
  y_text <- 0.23
  
  for (i in seq_along(labels)) {
    x0 <- x_centers[i] - box_w / 2
    x1 <- x_centers[i] + box_w / 2
    
    rect(
      x0, y_box - box_h / 2,
      x1, y_box + box_h / 2,
      col = fills[i], border = NA
    )
    segments(x0, y_box, x1, y_box, col = line_cols[i], lwd = 0.8)
    text(x_centers[i], y_text, labels[i], cex = cex_text)
  }
}

# ------------------------------------------------------------
# Helper: plot response curves
# ------------------------------------------------------------
plot_response_curves_meanSD <- function(
    vars,
    vars_max_df,
    vars_mean_df,
    vars_min_df,
    n_vars,
    alpha_val = 0.4
) {
  # ------------------------------------------------------------
  # Global graphics settings
  # ------------------------------------------------------------
  par(
    cex      = 1.0,
    cex.axis = 1.0,
    cex.lab  = 1.0,
    cex.main = 1.2
  )
  
  # ------------------------------------------------------------
  # Layout
  # ------------------------------------------------------------
  if (n_vars == 19) {
    panel_mat <- matrix(c(seq_len(19), 0), nrow = 4, byrow = TRUE)
    legend_id <- 20
    
    layout_mat <- rbind(
      panel_mat,
      rep(legend_id, 5)
    )
    
    layout(
      mat = layout_mat,
      heights = c(1, 1, 1, 1, 0.40)
    )
    
  } else {
    lay <- get_auto_layout(n_vars)
    
    n_slots <- lay$nrow * lay$ncol
    panel_mat <- matrix(
      c(seq_len(n_vars), rep(0, n_slots - n_vars)),
      nrow = lay$nrow,
      byrow = TRUE
    )
    
    legend_id <- n_vars + 1
    
    layout_mat <- rbind(
      panel_mat,
      rep(legend_id, lay$ncol)
    )
    
    layout(
      mat = layout_mat,
      heights = c(rep(1, lay$nrow), 0.40)
    )
  }
  
  # ------------------------------------------------------------
  # Draw panels
  # ------------------------------------------------------------
  par(mar = c(2.2, 2.4, 2.0, 0.6))
  
  for (i in seq_along(vars)) {
    plot_response_curve_panel(
      var_name     = vars[i],
      vars_max_df  = vars_max_df,
      vars_mean_df = vars_mean_df,
      vars_min_df  = vars_min_df,
      y_lim        = c(0, 1),
      alpha_val    = alpha_val
    )
  }
  
  # ------------------------------------------------------------
  # Draw legend
  # ------------------------------------------------------------
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot_horizontal_legend(alpha_val = alpha_val, cex_text = 1.1)
}

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
png(
  filename = file.path(
    job_fig_dir,
    paste0("ResponseCurves_meanSD_", myRespName, "_", jobid, ".png")
  ),
  width = 8, height = 5, units = "in", res = 300,
  type = "cairo-png",
  pointsize = 10
)
plot_response_curves_meanSD(
  vars = vars,
  vars_max_df = vars_max_df,
  vars_mean_df = vars_mean_df,
  vars_min_df = vars_min_df,
  n_vars = n_vars,
  alpha_val = 0.4
)
dev.off()

pdf(
  file = file.path(
    job_fig_dir,
    paste0("ResponseCurves_meanSD_", myRespName, "_", jobid, ".pdf")
  ),
  width = 8,
  height = 5,
  pointsize = 10
)
plot_response_curves_meanSD(
  vars = vars,
  vars_max_df = vars_max_df,
  vars_mean_df = vars_mean_df,
  vars_min_df = vars_min_df,
  n_vars = n_vars,
  alpha_val = 0.4
)
dev.off()
