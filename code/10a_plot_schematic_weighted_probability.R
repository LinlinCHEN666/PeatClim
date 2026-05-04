# ============================================================
# Script: 10a_plot_schematic_weighted_probability.R
# Purpose:
#   Create schematic figures illustrating the temperature-dependent
#   logistic weighting used to blend HTPeatland and LTPeatland
#   model probabilities.
#
# Inputs:
#   - No external data files required
#
# Outputs:
#   - Logistic weight schematic figure (PNG and PDF)
#   - Blended-probability heatmap schematic figure (PNG and PDF)
#
# Notes:
#   - These figures are conceptual illustrations of the blending method.
#   - The logistic weighting is defined by:
#       w_HT(T) = 1 / (1 + exp(-(T - t0) / w))
#       w_LT(T) = 1 - w_HT(T)
#   - The blended probability is:
#       P_HTLTPeatland = w_HT * P_HT + w_LT * P_LT
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(here)
  library(tibble)
})

here::i_am("code/10a_plot_schematic_weighted_probability.R")

# ---------------------------- #
# User settings
# ---------------------------- #

t0 <- 11.37   # split temperature (°C)
w  <- 2       # transition width (°C)

Tseq <- seq(-20, 35, by = 0.1)
Tvals <- c(t0 - 10, t0, t0 + 10)

# ---------------------------- #
# Project paths
# ---------------------------- #

dir_figures <- here("outputs")

dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

out_weights_png <- file.path(dir_figures, "10a_weighted_probability_logistic_weights.png")
out_weights_pdf <- file.path(dir_figures, "10a_weighted_probability_logistic_weights.pdf")
out_heatmap_png <- file.path(dir_figures, "10a_weighted_probability_heatmap.png")
out_heatmap_pdf <- file.path(dir_figures, "10a_weighted_probability_heatmap.pdf")

# ---------------------------- #
# Logistic weighting curves
# ---------------------------- #

w_HT <- 1 / (1 + exp(-(Tseq - t0) / w))
w_LT <- 1 - w_HT

df_w <- tibble(
  T = Tseq,
  w_HT = w_HT,
  w_LT = w_LT
) %>%
  pivot_longer(
    cols = c(w_HT, w_LT),
    names_to = "curve",
    values_to = "weight"
  )

p_weights <- ggplot(df_w, aes(T, weight, color = curve)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = t0, linetype = 2) +
  scale_y_continuous(
    limits = c(0, 1.01),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_x_continuous(
    breaks = sort(c(-20, -10, 0, t0, 20, 30)),
    labels = c("-20", "-10", "0", number_format(accuracy = 0.01)(t0), "20", "30")
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c("w_HT", "w_LT"),
    labels = c(
      expression(italic(w)[HTPeatland]),
      expression(italic(w)[LTPeatland])
    ),
    values = c(
      w_HT = "#E31A1C",
      w_LT = "#1F78B4"
    )
  ) +
  labs(
    x = expression("Temperature (" * italic(T) * ", " * degree * "C)"),
    y = expression("Weight (" * italic(w) * ")")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = c(0.90, 0.16),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    plot.background = element_rect(fill = "white", colour = "white")
  )

ggsave(out_weights_png, plot = p_weights, width = 6, height = 4, dpi = 300, bg = "white")
ggsave(out_weights_pdf, plot = p_weights, width = 6, height = 4, bg = "white")

message("Saved logistic-weight schematic: ", out_weights_png)
# ---------------------------- #
# Blended-probability heatmap
# ---------------------------- #

grid_f <- expand_grid(
  T = Tvals,
  P_LT = seq(0, 1, by = 0.01),
  P_HT = seq(0, 1, by = 0.01)
) %>%
  mutate(
    w_HT = 1 / (1 + exp(-(T - t0) / w)),
    P_HTLTPeatland = w_HT * P_HT + (1 - w_HT) * P_LT
  )

p_heatmap <- ggplot(grid_f, aes(P_LT, P_HT, fill = P_HTLTPeatland)) +
  geom_tile(width = 0.01, height = 0.01, na.rm = TRUE) +
  coord_fixed(expand = FALSE) +
  facet_wrap(
    ~ T,
    nrow = 1,
    labeller = label_bquote(italic(T) == .(T) * degree * C)
  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    limits = c(0, 1)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    limits = c(0, 1)
  ) +
  scale_fill_viridis_c(
    limits = c(0, 1),
    name = expression(italic(P)[HTLTPeatland])
  ) +
  labs(
    x = expression(italic(P)[LTPeatland]),
    y = expression(italic(P)[HTPeatland])
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    strip.background = element_rect(fill = "white", colour = NA),
    axis.title.x = element_text(margin = ggplot2::margin(t = 6), size = 12),
    axis.title.y = element_text(margin = ggplot2::margin(r = 6), size = 12),
    legend.key.height = unit(12, "pt"),
    legend.background = element_rect(fill = "white", colour = NA),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(out_heatmap_png, plot = p_heatmap, width = 10, height = 4, dpi = 300, bg = "white")
ggsave(out_heatmap_pdf, plot = p_heatmap, width = 10, height = 4, bg = "white")

message("Saved blended-probability heatmap: ", out_heatmap_png)
