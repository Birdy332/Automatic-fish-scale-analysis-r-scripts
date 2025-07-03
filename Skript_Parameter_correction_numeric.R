# =============================================================================
# Script:      Skript_allometric_relationships_scale_radius.R
# Author:      Christian Vogelmann
# Affiliation: LMU München
# Contact:     c.vogelmann@lmu.de
# Purpose:     Allometric analysis: fish scale radius (r_max) vs. fish length and weight
# Date:        2025-07-02
# Description: Standalone script for the published dataset 'Parameter_correction_numeric.csv'.
#              Produces all model fits and figures for allometric relationships.
# =============================================================================

rm(list=ls())

# ----------------------------
# 1. Libraries
# ----------------------------
library(dplyr)
library(ggplot2)

# ----------------------------
# 2. Read data
# ----------------------------
# Path to CSV (adjust if needed)
input_file <- "Parameter_correction_numeric.csv"
output_folder <- "./"   # Output folder for figures/statistics

# Create output directory if it doesn't exist
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

# Read in the data
dat <- read.csv(input_file, stringsAsFactors = FALSE)

# ----------------------------
# 3. Clean numeric columns
# ----------------------------
dat <- dat %>%
  mutate(
    r_max   = as.numeric(r_max),      # Scale radius in mm
    length  = as.numeric(length),     # Fish length in cm
    weight  = as.numeric(weight)      # Fish weight in g
  ) %>%
  filter(!is.na(r_max), !is.na(length), !is.na(weight), r_max > 0, length > 0, weight > 0)

# ----------------------------
# 4A. Power model: r_max vs. weight
# ----------------------------
model_weight <- lm(log(weight) ~ log(r_max), data = dat)
summary_weight <- summary(model_weight)
a_w <- signif(exp(coef(model_weight)[1]), 3)
b_w <- signif(coef(model_weight)[2], 3)
r2_w <- signif(summary_weight$r.squared, 3)
formel_w <- paste0(
  "y = ", a_w, " x^", b_w,
  "\nlog(y) = ", signif(coef(model_weight)[1], 3), " + ", signif(coef(model_weight)[2], 3), " log(x)",
  "\nR² = ", r2_w
)

# Log-log plot with linear fit
p_w <- ggplot(dat, aes(x = r_max, y = weight)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "firebrick",
    linewidth = 1,
    linetype = "dotted"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  annotate("text", 
           x = min(dat$r_max, na.rm = TRUE) * 1.2, 
           y = max(dat$weight, na.rm = TRUE) * 0.85, 
           label = formel_w, hjust = 0, vjust = 1, size = 4.2, fontface = "bold") +
  labs(
    x = "Scale radius r_max [mm]",
    y = "Fish weight [g]",
    title = "Allometric relationship: Scale radius vs. Fish weight"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(output_folder, "rmax_vs_weight_allometric.png"), plot = p_w, width = 8, height = 5, dpi = 600)

# Original scale (power law curve)
dat$pred_weight <- a_w * dat$r_max^b_w
p_w2 <- ggplot(dat, aes(x = r_max, y = weight)) +
  geom_point(color = "steelblue", size = 2, alpha = 0.7) +
  stat_function(fun = function(x) a_w * x^b_w, color = "firebrick", size = 1) +
  labs(
    x = "Scale radius r_max [mm]",
    y = "Fish weight [g]",
    title = "Allometric relationship: Scale radius vs. Fish weight (original scale)"
  ) +
  annotate("text", 
           x = min(dat$r_max, na.rm = TRUE) * 1.2, 
           y = max(dat$weight, na.rm = TRUE) * 0.85, 
           label = formel_w, hjust = 0, vjust = 1, size = 4.2, fontface = "bold") +
  theme_minimal(base_size = 13)

ggsave(file.path(output_folder, "rmax_vs_weight_originalscale.png"), plot = p_w2, width = 8, height = 5, dpi = 600)

# ----------------------------
# 4B. Power model: r_max vs. length
# ----------------------------
model_length <- lm(log(length) ~ log(r_max), data = dat)
summary_length <- summary(model_length)
a_l <- signif(exp(coef(model_length)[1]), 3)
b_l <- signif(coef(model_length)[2], 3)
r2_l <- signif(summary_length$r.squared, 3)
formel_l <- paste0(
  "y = ", a_l, " x^", b_l,
  "\nlog(y) = ", signif(coef(model_length)[1], 3), " + ", signif(coef(model_length)[2], 3), " log(x)",
  "\nR² = ", r2_l
)

# Log-log plot with linear fit
p_l <- ggplot(dat, aes(x = r_max, y = length)) +
  geom_point(color = "darkorange", size = 2, alpha = 0.7) +
  stat_smooth(
    method = "lm",
    formula = y ~ x,
    se = FALSE,
    color = "firebrick",
    linewidth = 1,
    linetype = "dotted"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  annotate("text", 
           x = min(dat$r_max, na.rm = TRUE) * 1.2, 
           y = max(dat$length, na.rm = TRUE) * 0.85, 
           label = formel_l, hjust = 0, vjust = 1, size = 4.2, fontface = "bold") +
  labs(
    x = "Scale radius r_max [mm]",
    y = "Fish length [cm]",
    title = "Allometric relationship: Scale radius vs. Fish length"
  ) +
  theme_minimal(base_size = 13)

ggsave(file.path(output_folder, "rmax_vs_length_allometric.png"), plot = p_l, width = 8, height = 5, dpi = 600)

# Original scale (power law curve)
dat$pred_length <- a_l * dat$r_max^b_l
p_l2 <- ggplot(dat, aes(x = r_max, y = length)) +
  geom_point(color = "darkorange", size = 2, alpha = 0.7) +
  stat_function(fun = function(x) a_l * x^b_l, color = "firebrick", size = 1) +
  labs(
    x = "Scale radius r_max [mm]",
    y = "Fish length [cm]",
    title = "Allometric relationship: Scale radius vs. Fish length (original scale)"
  ) +
  annotate("text", 
           x = min(dat$r_max, na.rm = TRUE) * 1.2, 
           y = max(dat$length, na.rm = TRUE) * 0.85, 
           label = formel_l, hjust = 0, vjust = 1, size = 4.2, fontface = "bold") +
  theme_minimal(base_size = 13)

ggsave(file.path(output_folder, "rmax_vs_length_originalscale.png"), plot = p_l2, width = 8, height = 5, dpi = 600)

# ----------------------------
# 5. Export numeric (cleaned) table for reproducibility
# ----------------------------
write.csv(dat, file = file.path(output_folder, "Parameter_correction_numeric.csv"), row.names = FALSE)

# ----------------------------
# 6. Print model parameters (for supplement)
# ----------------------------
cat("Weight model: y =", a_w, "* x^", b_w, " | R² =", r2_w, "\n")
cat("Length model: y =", a_l, "* x^", b_l, " | R² =", r2_l, "\n")

# ----------------------------
# 7. Save R session info
# ----------------------------
writeLines(capture.output(sessionInfo()), file.path(output_folder, "session_info.txt"))

# ----------------------------
# 8. (Optional) Show plots in RStudio Viewer
# ----------------------------
print(p_w)
print(p_l)