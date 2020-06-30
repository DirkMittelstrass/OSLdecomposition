library(Luminescence)
library(OSLdecomposition)
report <- FALSE # Disable auto-reporting to save a lot of computing time

# 120 ka old sibirian sample from Maggi Fuchs
# late background analysis shows age under-estimation

data_path <- system.file("examples", "Oy7-01-14_63-100_1mm_Qz1-1_ed.BIN", package = "OSLdecomposition")
Oy7 <- read_BIN2R(data_path, fastForward = TRUE)

# Step 0 (There are pre-bleaching steps in aliquot 1,2,3,4. They will be renamed to "OSL2")
Oy7 <- RLum.OSL_correction(Oy7)

# Step 1
Oy7 <- RLum.OSL_global_fitting(Oy7, report = report)

# Step 2
Oy7 <- RLum.OSL_decomposition(Oy7, report = report)

# Step 3
Oy7_fast <- analyse_SAR.CWOSL_decomposed(Oy7, OSL.component = "fast")
Oy7_medium <- analyse_SAR.CWOSL_decomposed(Oy7, OSL.component = 2)

# Classic analysis
Oy7_late_background <- analyse_SAR.CWOSL_decomposed(Oy7, 1, 3, 80, 100)

# Compare results
plot_RadialPlot(list(Oy7_late_background, Oy7_fast, Oy7_medium))
