# Automatic-fish-scale-analysis – R Script Repository

This repository contains all R scripts used for data analysis, visualization, and statistical validation in the manuscript:

**"Automatic fish scale analysis: age determination, annuli and circuli detection, length and weight back-calculation of coregonid scales"**  
Vogelmann et al., 2025 (in review)

---

##  Included Scripts

###  `Skript_ecological_dataset_center.R`
- Input: `comparison_all_scales.csv`
- Purpose: Comparison of automatically vs. manually verified measurements
- Generates:
  - Core statistics table: bias, LoA, paired tests
  - Boxplots and Bland–Altman plots for all biometric variables
  - Center correction analysis with scatterplots, histograms, and summary table
- Output:  
  - `comparison_stats_core_variables.csv`  
  - `scatter_center_comparison.png`, `histogram_center_error_mm_full.png`, `center_correction_bar.png`, ...  
  - ~20 plots and 2 summary files  

###  `Skript_validation_manual_vs_automated_circ_area.R`
- Input: `Validation_data.csv`
- Purpose: Validation of circuli count and scale area based on microscope measurements
- Generates:
  - Regression statistics (r, R², bias, LoA)
  - Scatterplots for area and circuli
  - Statistical tests (Shapiro, t-test, Wilcoxon)
- Output:  
  - `Validation_statistics.csv`  
  - `Scatterplot_Area.png`, `Scatterplot_Circuli.png`

###  `Skript_allometric_relationships_scale_radius.R`
- Input: `Parameter_correction_numeric.csv`
- Purpose: Calibration of allometric models (scale radius vs. length and weight)
- Generates:
  - Power law regressions (log-log and original scale)
  - Exported parameters for model fitting
- Output:  
  - `rmax_vs_length_allometric.png`, `rmax_vs_weight_allometric.png`, etc.  
  - `Parameter_correction_numeric.csv` (cleaned version)

---

##  Datasets

All required `.csv` input files can be downloaded from the related Figshare dataset:  
https://doi.org/10.6084/m9.figshare.29467970

---

##  Reproducibility

Each script is standalone and can be run in R or RStudio. All generated plots and result tables are saved automatically in the working directory. A complete `session_info.txt` is included in each script run for reproducibility.

---

##  Contact & Software

The full Coregon Analyzer software (including GUI and image-processing routines) is available from the authors upon request and provided under a research-use license.

For questions:  
**Christian Vogelmann**  
LMU München  
 c.vogelmann@lmu.de

---

##  Citation

Vogelmann et al. (2025). _Automatic fish scale analysis: age determination, annuli and circuli detection, length and weight back-calculation of coregonid scales_. *Ecological Informatics*, [in review].

##  License

- All **R scripts** in this repository are licensed under the [MIT License](LICENSE).
- All associated **data files** are archived on Figshare under the [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/):  
  https://doi.org/10.6084/m9.figshare.29467970
- The full executable for automated scale analysis (`Algorithm.py`) is provided via Figshare.  
  The complete program (with interface) is available upon request from the authors.
