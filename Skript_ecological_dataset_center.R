# =============================================================================
# Script:      Skript_ecological_dataset_center.R
# Author:      Christian Vogelmann
# Affiliation: LMU München
# Contact:     c.vogelmann@lmu.de
# Purpose:     Comparison of manual verification vs. automatic fish scale measurements.
# Date:        2025-07-02
# Description: Standalone analysis and visualization script for the published
#              dataset 'comparison_all_scales.csv', including center correction.
# =============================================================================

rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(viridis)
library(patchwork)
library(BlandAltmanLeh)

# --- 1. Read merged dataset (provided as supplement) ---
# Path to CSV (adjust if necessary)
data_path <- "comparison_all_scales.csv"
save_path <- "./"   # Output folder (adjust if needed)

merged <- read.csv(data_path, stringsAsFactors = FALSE)

# --- 2. Core statistics table (Bland-Altman, errors, bias etc.) ---
comparison_pairs <- list(
  circuli_r1              = c("circuli r1_val",        "circuli r1_auto"),
  circuli_r2              = c("circuli r2_val",        "circuli r2_auto"),
  circuli_r3              = c("circuli r3_val",        "circuli r3_auto"),
  r1                      = c("r1_val",                "r1_auto"),
  r2                      = c("r2_val",                "r2_auto"),
  r3                      = c("r3_val",                "r3_auto"),
  annuli                  = c("anuli_val",             "anuli_auto"),
  circuli_dist_r1         = c("circuli dist r1_val",   "circuli dist r1_auto"),
  circuli_dist_r2         = c("circuli dist r2_val",   "circuli dist r2_auto"),
  circuli_dist_r3         = c("circuli dist r3_val",   "circuli dist r3_auto"),
  r_max                   = c("r_max [mm]_val",        "r_max [mm]_auto"),
  circuli                 = c("circuli_val",           "circuli_auto"),
  center_x                = c("center_x [pixel]_val",  "center_x [pixel]_auto"),
  center_y                = c("center_y [pixel]_val",  "center_y [pixel]_auto"),
  length_r1               = c("length r1_val",         "length r1_auto"),
  length_r2               = c("length r2_val",         "length r2_auto"),
  length_r3               = c("length r3_val",         "length r3_auto"),
  length_alternative_r1   = c("length alternative r1_val", "length alternative r1_auto"),
  length_alternative_r2   = c("length alternative r2_val", "length alternative r2_auto"),
  length_alternative_r3   = c("length alternative r3_val", "length alternative r3_auto")
)

stat_table <- data.frame()

for (name in names(comparison_pairs)) {
  var_val <- comparison_pairs[[name]][1]
  var_auto <- comparison_pairs[[name]][2]
  
  if (!(var_val %in% names(merged)) | !(var_auto %in% names(merged))) {
    cat("Missing column:", name, "\n")
    next
  }
  
  v <- as.numeric(merged[[var_val]])
  a <- as.numeric(merged[[var_auto]])
  idx <- complete.cases(v, a)
  v <- v[idx]
  a <- a[idx]
  
  cat("Comparison:", name, "| Pairs:", length(v), "\n")
  
  if (length(v) < 2) {
    cat("Too few valid pairs for", name, "\n")
    next
  }
  
  diff <- a - v
  p_norm <- tryCatch(shapiro.test(diff)$p.value, error=function(e) NA)
  if (!is.na(p_norm) && p_norm > 0.05) {
    test <- t.test(a, v, paired = TRUE)
    testtype <- "paired t-test"
  } else {
    test <- wilcox.test(a, v, paired = TRUE)
    testtype <- "Wilcoxon"
  }
  
  rel_err <- mean(abs(diff/v), na.rm = TRUE) * 100
  abs_err <- mean(abs(diff), na.rm = TRUE)
  bias <- mean(diff, na.rm=TRUE)
  loa  <- 1.96 * sd(diff, na.rm=TRUE)
  lo_lower <- bias - loa
  lo_upper <- bias + loa
  
  stat_row <- data.frame(
    Parameter = name,
    N = length(v),
    Mean_valid = mean(v, na.rm=TRUE),
    Mean_auto  = mean(a, na.rm=TRUE),
    SD_valid   = sd(v, na.rm=TRUE),
    SD_auto    = sd(a, na.rm=TRUE),
    AbsError   = abs_err,
    RelErrorPct= rel_err,
    Bias = bias,
    LoA_lower = lo_lower,
    LoA_upper = lo_upper,
    Test      = testtype,
    Test_p    = signif(test$p.value, 3)
  )
  
  stat_table <- rbind(stat_table, stat_row)
}

write.csv(stat_table, file.path(save_path, "comparison_stats_core_variables.csv"), row.names=FALSE)
print(stat_table)

# --- 3. Boxplots und Bland-Altman Plots ---
for (name in names(comparison_pairs)) {
  var_val <- comparison_pairs[[name]][1]
  var_auto <- comparison_pairs[[name]][2]
  
  if (!(var_val %in% names(merged)) | !(var_auto %in% names(merged))) next
  v <- as.numeric(merged[[var_val]])
  a <- as.numeric(merged[[var_auto]])
  idx <- complete.cases(v, a)
  v <- v[idx]; a <- a[idx]
  
  if (length(v) < 2) next
  
  # Boxplot
  png(file.path(save_path, paste0("boxplot_", name, ".png")), width=900, height=600, res=150)
  boxplot(v, a, names = c("Manual", "Automatic"),
          main = paste("Boxplot:", name),
          ylab = name, col = c("skyblue", "salmon"))
  dev.off()
  
  # Bland-Altman-Plot
  png(file.path(save_path, paste0("bland_altman_", name, ".png")), width=900, height=600, res=150)
  means <- (v + a) / 2
  diffs <- a - v
  bias <- mean(diffs, na.rm=TRUE)
  loa  <- 1.96 * sd(diffs, na.rm=TRUE)
  lo_lower <- bias - loa
  lo_upper <- bias + loa
  
  plot(means, diffs,
       main = paste("Bland-Altman:", name),
       xlab = "Mean (Manual/Auto)", ylab = "Difference (Auto - Manual)")
  abline(h = bias, col = "blue", lwd = 2)
  abline(h = lo_upper, col = "red", lty = 2)
  abline(h = lo_lower, col = "red", lty = 2)
  legend("topright", legend = c(
    paste("Bias =", signif(bias, 3)),
    paste("LoA =", signif(lo_lower, 3), "-", signif(lo_upper, 3))
  ), bty = "n")
  dev.off()
}

# --- 4. Center Correction Visualizations and Statistics ---

# Ensure numeric center columns
merged$`center_x [pixel]_auto` <- as.numeric(merged$`center_x [pixel]_auto`)
merged$`center_y [pixel]_auto` <- as.numeric(merged$`center_y [pixel]_auto`)
merged$`center_x [pixel]_val`  <- as.numeric(merged$`center_x [pixel]_val`)
merged$`center_y [pixel]_val`  <- as.numeric(merged$`center_y [pixel]_val`)

# Scatterplot: Manual vs. Automatic centers, with lines
plot_long <- merged %>%
  select(filename, index,
         x_auto = `center_x [pixel]_auto`, y_auto = `center_y [pixel]_auto`,
         x_val  = `center_x [pixel]_val`,  y_val  = `center_y [pixel]_val`) %>%
  pivot_longer(
    cols = c(x_auto, y_auto, x_val, y_val),
    names_to = c("coord", "type"),
    names_pattern = "(x|y)_(auto|val)",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = coord, values_from = value) %>%
  mutate(type = recode(type, auto = "Automatic", val = "Verified"))

seg_df <- merged %>%
  select(filename, index,
         x_auto = `center_x [pixel]_auto`, y_auto = `center_y [pixel]_auto`,
         x_val  = `center_x [pixel]_val`,  y_val  = `center_y [pixel]_val`)

p_scatter <- ggplot() +
  geom_segment(
    data = seg_df,
    aes(x = x_auto, y = y_auto, xend = x_val, yend = y_val),
    color = "gray60", linewidth = 0.35, alpha = 0.6
  ) +
  geom_point(
    data = plot_long,
    aes(x = x, y = y, color = type),
    size = 2.2, alpha = 0.9
  ) +
  scale_color_manual(values = c("Automatic" = "deepskyblue3", "Verified" = "firebrick3")) +
  labs(
    x = "X coordinate (pixel)",
    y = "Y coordinate (pixel)",
    title = "Comparison of Centers: Automatic vs. Verified",
    color = "Center"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")
print(p_scatter)
ggsave(file.path(save_path, "scatter_center_comparison.png"), plot = p_scatter, width = 8, height = 7, dpi = 600)

# Error calculation (in pixel and mm)
dpi <- 3600
mm_per_pixel <- 25.4 / dpi
merged$error_px <- sqrt(
  (merged$`center_x [pixel]_auto` - merged$`center_x [pixel]_val`)^2 +
    (merged$`center_y [pixel]_auto` - merged$`center_y [pixel]_val`)^2
)
merged$error_mm <- merged$error_px * mm_per_pixel

# Histogram: Error in mm (full range)
p_hist_mm_full <- ggplot(merged, aes(x = error_mm)) +
  geom_histogram(binwidth = 0.08, fill = "mediumorchid3", color = "white") +
  labs(
    title = "Distribution of Center Deviations (mm) – All Cases",
    x = "Error (mm)",
    y = "Count"
  ) +
  theme_minimal(base_size = 15)
print(p_hist_mm_full)
ggsave(file.path(save_path, "histogram_center_error_mm_full.png"), plot = p_hist_mm_full, width = 7, height = 5, dpi = 600)

# Histogram: Error in mm (zoom: up to 0.5 mm)
p_hist_mm_zoom <- ggplot(merged[merged$error_mm <= 0.5, ], aes(x = error_mm)) +
  geom_histogram(binwidth = 0.02, fill = "deepskyblue4", color = "white", boundary = 0) +
  scale_x_continuous(
    limits = c(0, 0.5),
    breaks = seq(0, 0.5, 0.1)
  ) +
  labs(
    title = "Distribution of Center Deviations [mm] (up to 0.5 mm)",
    x = "Error [mm]",
    y = "Count"
  ) +
  theme_minimal(base_size = 15)
print(p_hist_mm_zoom)
ggsave(file.path(save_path, "histogram_center_error_mm_zoom.png"), plot = p_hist_mm_zoom, width = 7, height = 5, dpi = 600)

# Barplot: Number of cases with/without correction
merged$error_class <- ifelse(merged$error_px == 0, "No correction", "Correction")

p_bar <- ggplot(merged, aes(x = error_class, fill = error_class)) +
  geom_bar() +
  scale_fill_manual(values = c("No correction" = "forestgreen", "Correction" = "darkred")) +
  labs(
    title = "Cases With and Without Center Correction",
    x = "",
    y = "Number of cases"
  ) +
  theme_minimal(base_size = 15) +
  theme(legend.position = "none")
print(p_bar)
ggsave(file.path(save_path, "center_correction_bar.png"), plot = p_bar, width = 6, height = 5, dpi = 600)

# Statistical Summary (printed and saved)
n_total    <- nrow(merged)
n_nocorr   <- sum(merged$error_px == 0)
n_corr     <- n_total - n_nocorr
n_outliers <- sum(merged$error_mm > 0.5)
mean_error <- mean(merged$error_mm)
median_error <- median(merged$error_mm)
max_error  <- max(merged$error_mm)
percent_nocorr <- round(100 * n_nocorr / n_total, 2)
percent_outlier <- round(100 * n_outliers / n_total, 2)

summary_df <- data.frame(
  Total_cases = n_total,
  Cases_without_correction = n_nocorr,
  Cases_with_correction = n_corr,
  Percent_without_correction = percent_nocorr,
  Mean_error_mm = mean_error,
  Median_error_mm = median_error,
  Max_error_mm = max_error,
  Outlier_cases_gt_0_5mm = n_outliers,
  Percent_outliers_gt_0_5mm = percent_outlier
)

print(summary_df)
write.csv(summary_df, file.path(save_path, "center_error_summary.csv"), row.names = FALSE)

cat("Summary:\n")
cat("Total cases:        ", n_total, "\n")
cat("Without correction: ", n_nocorr, "(", percent_nocorr, "%)\n")
cat("With correction:    ", n_corr, "\n")
cat("Outliers (>0.5 mm): ", n_outliers, "(", percent_outlier, "%)\n")
cat("Mean error (mm):    ", round(mean_error,4), "\n")
cat("Median error (mm):  ", round(median_error,4), "\n")
cat("Max error (mm):     ", round(max_error,4), "\n")

writeLines(capture.output(sessionInfo()), file.path(save_path, "session_info.txt"))