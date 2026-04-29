# Settings ----------------------------------------------------------------
library(Luminescence)
library(OSLdecomposition)
report <- FALSE #  Enable auto-reporting to get detailed results about information in step 1 and step 2 results

# Load data ---------------------------------------------------------------
# BIN or BINX file of a quartz OSL data set measured on a Riso OSL/TL reader with
# the Single-Aliquot dose Regeneration (SAR) protocol
data_path <- file.choose()
data_raw <- read_BIN2R(data_path, fastForward = TRUE)


# DATA PREPARATION  ----------------------------------------------------

# While optional, it is recommended to check the input data consistency and
# correct it for know issues
data_corrected <- RLum.OSL_correction(
  data_raw,
  limit_duration = 20, # [s], reduces the recording time. This prevents the fitting from being overly weighted towards slow components
  PMT_pulse_pair_resolution = 18, # [ns], PMT dead-time correction, see Kreutzer and Mittelstrass (2026)
  background_sequence = NULL # optional: List index of empty aliquot(s) for background correction
)

# STEP 1: Determine global OSL decay parameter ---------------------------

data_fitted <- RLum.OSL_global_fitting(
  data_corrected,
  stimulation_intensity = 80, # [mW/cm²], optional but necessary to estimate photo-ionisation cross-sections
  stimulation_wavelength = 470, # [nm], optional but necessary to estimate photo-ionisation cross-sections
  report = report
)

# Show result of the global curve fitting
plot_OSLcurve(curve = data_fitted$FITTING$curve,
              components = data_fitted$FITTING$components)

# STEP 2: Determine global OSL decay parameter ---------------------------
data_decomposed <- RLum.OSL_decomposition(
  data_fitted,
  report = report)

# Plot the very first OSL record of the data set and its decomposition result
data_OSL_only <- get_RLum(data_decomposed, recordType = "OSL")
plot_OSLcurve(data_OSL_only[[1]][[1]])

# STEP 3: Equivalent dose (De) calculation and dose statistics ------------------------------------------------
# Use the R package Luminescence function for dose evalution according to the
# protocol of Murray and Wintle (2000)
fast_De <- analyse_SAR.CWOSL(
  data_decomposed,
  OSL.component = data_decomposed$DECOMPOSITION$dominating.component, # Usually the Fast component
  plot = FALSE
)

# For comparison
late_background_De <-
  analyse_SAR.CWOSL(
    data_corrected,
    signal_integral = c(1, 10), # Channel indices
    background_integral = c(800, 999), # Channel indices
    plot = FALSE)

# Show the equivalent dose distribution in accordance to Galbraith & Green (1990)
# black: late background De's, red: fast component De's
plot_AbanicoPlot(
  data = list(late_background_De, fast_De))
