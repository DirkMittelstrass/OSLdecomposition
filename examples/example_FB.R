library(Luminescence)
library(OSLdecomposition)
report <- TRUE # To save a lot of computing time, disable auto-reporting

# FB_10Gy is a dose recovery test with the La Fontainebleau quartz in the Bayreuth lexsyg system with green stimulation
data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
FB <- read_BIN2R(data_path, fastForward = TRUE)


# Adequate background correction is recommended to prevent slow component overestimation.
# Here, the OSL measurements of an empty aliquot (data set entry 11) is subtracted from the other OSL
# measurements
FB_corrected <- RLum.OSL_correction(FB, background = 11)

FB_fitted <- RLum.OSL_global_fitting(FB_corrected,
                                     stimulation_intensity = 50,
                                     stimulation_wavelength = 530,
                                     report = report)

FB_decomposed <-RLum.OSL_decomposition(FB_fitted,
                                       report = report)

FB_fast_De <- analyse_SAR.CWOSL_decomposed(FB_decomposed,
                                           OSL.component = 1)


# be aware that late background subtraction results in overestimated background.error
# for background corrected data
FB_late_background <- analyse_SAR.CWOSL_decomposed(FB, 1, 10, 800, 999)

# black: late background, red: fast component
plot_KDE(list(FB_late_background,
              FB_fast_De))


