library(Luminescence)
library(OSLdecomposition)

data_path <- system.file("examples", "Oy7-01-14_63-100_1mm_Qz1-1_ed.BIN", package = "OSLdecomposition")
Oy7 <- read_BIN2R(data_path, fastForward = TRUE)

Oy7 <- RLum.OSL_correction(Oy7)

Oy7 <- RLum.OSL_global_fitting(Oy7)

Oy7 <- RLum.OSL_decomposition(Oy7)

Oy7_fast <- analyse_SAR.CWOSL(Oy7, OSL.component = 1)
Oy7_medium <- analyse_SAR.CWOSL(Oy7, OSL.component = 2)

Oy7_late_background <- analyse_SAR.CWOSL(BT594, 1, 3, 80, 100, plot = FALSE)

plot_RadialPlot(list(Oy7_late_background, Oy7_fast, Oy7_medium))