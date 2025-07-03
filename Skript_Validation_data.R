# =============================================================================
# Script:      Skript_validation_manual_vs_automated_circ_area.R
# Author:      Christian Vogelmann
# Affiliation: LMU München
# Contact:     c.vogelmann@lmu.de
# Purpose:     Comparison of manual vs. automated fish scale measurements
# Date:        2025-07-02
# Description: Standalone analysis and visualization script for the
#              dataset 'Validation_data.csv'. Produces all summary statistics
#              and regression/validation plots for scale area and circuli count.
# =============================================================================

rm(list=ls())

# ----------------------------
# 1. Libraries
# ----------------------------
library(ggplot2)

# ----------------------------
# 2. Read data
# ----------------------------

# Path to CSV (adjust if needed)
data_file <- "Validation_data.csv"
output_dir <- "./"   # Output folder for figures/statistics

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read in the data
dat <- read.csv(data_file, stringsAsFactors = FALSE)

# ----------------------------
# 3. Variable assignment
# ----------------------------
# Rename columns if needed for clarity
dat$Manual_area <- as.numeric(dat$Manual_area)             # or adapt if different column names
dat$Automatic_area <- as.numeric(dat$Automatic_area)
dat$Manual_circuli <- as.numeric(dat$Manual_circuli)
dat$Automatic_circuli <- as.numeric(dat$Automatic_circuli)

# ----------------------------
# 4. Regression and statistics for area
# ----------------------------
lm_area <- lm(Automatic_area ~ Manual_area, data = dat)
summary_area <- summary(lm_area)
cor_area <- cor(dat$Manual_area, dat$Automatic_area, use = "complete.obs")
r2_area <- summary_area$r.squared
diff_area <- dat$Automatic_area - dat$Manual_area
bias_area <- mean(diff_area, na.rm=TRUE)
sd_area <- sd(diff_area, na.rm=TRUE)
loa_area_lower <- bias_area - 1.96*sd_area
loa_area_upper <- bias_area + 1.96*sd_area

# ----------------------------
# 5. Plot area: scatter/regression
# ----------------------------
plot_area <- ggplot(dat, aes(x = Manual_area, y = Automatic_area)) +
  geom_point(alpha=0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Manual area [mm²]", y = "Automatic area [mm²]", title = "Regression of scale area: manual vs. automatic") +
  annotate("text", x = max(dat$Manual_area, na.rm=TRUE)*0.75, y = min(dat$Automatic_area, na.rm=TRUE)*1.05,
           label = paste0("r = ", round(cor_area, 2), "\nR² = ", round(r2_area, 2),
                          "\nBias = ", round(bias_area,2), " mm²\nLoA: [", round(loa_area_lower,2), ", ", round(loa_area_upper,2), "]"),
           color="black", size=4) +
  theme_minimal()
ggsave(filename = file.path(output_dir, "Scatterplot_Area.png"), plot = plot_area, width = 7, height = 6, dpi = 600)
print(plot_area)

# ----------------------------
# 6. Regression/statistics for circuli
# ----------------------------
lm_circ <- lm(Automatic_circuli ~ Manual_circuli, data=dat)
summary_circ <- summary(lm_circ)
cor_circ <- cor(dat$Manual_circuli, dat$Automatic_circuli, use = "complete.obs")
r2_circ <- summary_circ$r.squared
diff_circ <- dat$Automatic_circuli - dat$Manual_circuli
bias_circ <- mean(diff_circ, na.rm=TRUE)
sd_circ <- sd(diff_circ, na.rm=TRUE)
loa_circ_lower <- bias_circ - 1.96*sd_circ
loa_circ_upper <- bias_circ + 1.96*sd_circ

# ----------------------------
# 7. Plot circuli: scatter/regression
# ----------------------------
plot_circ <- ggplot(dat, aes(x = Manual_circuli, y = Automatic_circuli)) +
  geom_point(alpha=0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  labs(x = "Manual circuli count", y = "Automatic circuli count", title = "Regression of circuli count: manual vs. automatic") +
  annotate("text", x = max(dat$Manual_circuli, na.rm=TRUE)*0.75, y = min(dat$Automatic_circuli, na.rm=TRUE)*1.05,
           label = paste0("r = ", round(cor_circ, 2), "\nR² = ", round(r2_circ, 2),
                          "\nBias = ", round(bias_circ,2), "\nLoA: [", round(loa_circ_lower,2), ", ", round(loa_circ_upper,2), "]"),
           color="black", size=4) +
  theme_minimal()
ggsave(filename = file.path(output_dir, "Scatterplot_Circuli.png"), plot = plot_circ, width = 7, height = 6, dpi = 600)
print(plot_circ)

# ----------------------------
# 8. Statistical tests
# ----------------------------
shapiro_area <- shapiro.test(diff_area)
shapiro_circ <- shapiro.test(diff_circ)

t_area <- t.test(dat$Automatic_area, dat$Manual_area, paired = TRUE)
w_area <- wilcox.test(dat$Automatic_area, dat$Manual_area, paired = TRUE)
t_circ <- t.test(dat$Automatic_circuli, dat$Manual_circuli, paired = TRUE)
w_circ <- wilcox.test(dat$Automatic_circuli, dat$Manual_circuli, paired = TRUE)

# ----------------------------
# 9. Save summary statistics
# ----------------------------
results <- data.frame(
  Parameter = c("Area", "Circuli"),
  N = c(sum(complete.cases(dat$Automatic_area, dat$Manual_area)),
        sum(complete.cases(dat$Automatic_circuli, dat$Manual_circuli))),
  Correlation = c(cor_area, cor_circ),
  R2 = c(r2_area, r2_circ),
  Bias = c(bias_area, bias_circ),
  LoA_lower = c(loa_area_lower, loa_circ_lower),
  LoA_upper = c(loa_area_upper, loa_circ_upper),
  Shapiro_p = c(shapiro_area$p.value, shapiro_circ$p.value),
  ttest_p = c(t_area$p.value, t_circ$p.value),
  wilcox_p = c(w_area$p.value, w_circ$p.value)
)
write.csv(results, file.path(output_dir, "Validation_statistics.csv"), row.names = FALSE)

# ----------------------------
# 10. Save processed data for reproducibility
# ----------------------------
write.csv(dat, file.path(output_dir, "Validation_data.csv"), row.names = FALSE)

# ----------------------------
# 11. Save R session info
# ----------------------------
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))