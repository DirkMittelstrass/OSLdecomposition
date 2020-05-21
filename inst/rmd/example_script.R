library(Luminescence)
library(OSLdecomposition)

FB <- read_BIN2R(file.choose(), fastForward = TRUE)

FB_corrected <- RLum.OSL_correction(FB, background = 11)

FB_fitted <- RLum.OSL_global_fitting(FB_corrected,
                                     stimulation_intensity = 50,
                                     stimulation_wavelength = 530)

FB_decomposed <-RLum.OSL_decomposition(FB_fitted)

FB_fast_De <- analyse_SAR.CWOSL(FB_decomposed, OSL.component = 1)

# IMPORTANT: Remove "Luminescence:::" before building
# results <- analyse_SAR.CWOSL(lum_data, OSL.component = "fast")



