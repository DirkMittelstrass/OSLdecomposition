#' Identifies CW-OSL signal components in RLum.Analysis data sets
#'
#' This function identifies CW-OSL signal components in [RLum.Analysis-class] data sets
#' from the [Luminescence-package]. First, all CW-OSL records are combined to one
#' global average CW-OSL curve, then the multi-exponential fitting approach of
#' Bluszcz and Adamiec (2006) is applied.
#'
#'
#' The workflow is as following:
#'
#' \enumerate{
#'   \item{[sum_OSLcurve]: Combine all measurement of type `record_type` to one global average curve}
#'   \item{[fit_OSLcurves]: Identify  component by a multi-exponential fitting}
#'   \item{Create a `html` report to summarize the results (*optional*)}
#' }
#'
#' The function processes just data sets created within the [Luminescence-package] (Kreutzer et al. 2012).
#' Data sets must be formatted as [RLum.Analysis-class] objects. Output objects will also be
#' [RLum.Analysis-class] objects. The data set should be processed with [R.Lum_correction] and is
#' meant for further analysis with [RLum.OSL_decomposition] afterwards.
#'
#' The `html` report is rendered by the [rmarkdown-package] and saved in the working directory,
#' which is usually the directory of the data set file. The correct display of the report needs
#' a internet connection to allow formula encoding by *MathJax*.
#' The *Rmarkdown* source code can be found with the following command:
#'
#' `system.file("rmd", "report_Step1.Rmd", package = "OSLdecomposition")`
#'
#'
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' Data set of one or multiple CW-OSL measured aliquots
#'
#' @param record_type [character] (*with default*):
#' Type of records, selected by the [RLum.Analysis-class] attribute `@recordType`.
#' Common are: `"OSL"`,`"SGOSL"` or `"IRSL"`
#'
#' @param max_components [numeric] (*with default*):
#' Maximum number of components *K*, see [fit_OSLcurve]
#'
#' @param F_threshold [numeric] (*with default*):
#' Fitting stop criterion, see [fit_OSLcurve]
#'
#' @param stimulation_intensity [numeric] (*with default*):
#' Intensity of optical stimulation in *mW / cm²*. Used to calculate photo-ionisation cross-sections, see [fit_OSLcurve]
#'
#' @param stimulation_wavelength [numeric] (*with default*):
#' Wavelength of optical stimulation in *nm*. Used to calculate photo-ionisation cross-sections, see [fit_OSLcurve]
#'
#' @param report [logical] (*with default*):
#' Creates a `html` report, saves it in the working directory and opens it in your
#' standard browser. The report contains the results and further information
#' on the data processing
#'
#' @param image_format [character] (*with default*):
#' Image format of the automatically saved graphs if `report = TRUE`.
#' Allowed are `.pdf`, `.eps`, `.svg` (vector graphics), `.jpg`, `.png`, `.bmp` (pixel graphics)
#' and more, see [ggplot2::ggsave]. The images are saved in the working directory subfolder `/report_figures`
#'
#' @param verbose [logical] (*with default*):
#' Enables console text output
#'
#'
#' @section Last updates:
#'
#' 2020-11-06, DM: Added roxygen documentation
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite the package the following way:
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' Automated identification and separation of quartz CW-OSL signal components with R, *in preparation*.
#'
#' @seealso [RLum.OSL_correction], [RLum.OSL_decomposition], [sum_OSLcurves], [fit_OSLcurve]
#'
#' @references
#'
#' Bluszcz, A., Adamiec, G., 2006. Application of differential evolution to fitting OSL
#' decay curves. Radiation Measurements 41, 886–891.
#'
#' Kreutzer, S., Schmidt, C., Fuchs, M.C., Dietze, M., Fischer, M., Fuchs, M., 2012.
#' Introducing an R package for luminescence dating analysis. Ancient TL, 30 (1), 1-8.
#'
#' @return
#'
#' The input `object`, a [list] of [RLum.Analysis-class] objects is given back.
#' The returned data set contains a new list element `object[["OSL_COMPONENTS"]]` which provides
#' a [list] of the input parameters, the global average OSL curve from [sum_OSLcurves] and
#' the `output.complex = TRUE` results from [fit_OSLcurve]
#'
#' @examples
#'
#' # 'FB_10Gy' is a dose recovery test with the La Fontainebleau quartz
#' # measured in a lexsyg research with green LED stimulation
#' data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
#' data_set <- Luminescence::read_BIN2R(data_path, fastForward = TRUE)
#'
#' # Check data set and perform background correction
#' data_set_corrected <- RLum.OSL_correction(data_set, background = 11)
#'
#' # Identify components and create report
#' data_set_fitted <- RLum.OSL_global_fitting(data_set_corrected,
#'                                            max_components = 3,
#'                                            stimulation_intensity = 50,
#'                                            stimulation_wavelength = 530)
#'
#' @md
#' @export

RLum.OSL_global_fitting <- function(object,
                                    record_type = "OSL",
                                    max_components = 5,
                                    F_threshold = 150,
                                    stimulation_intensity = 35,
                                    stimulation_wavelength = 470,
                                    report = TRUE,
                                    image_format = "pdf",
                                    verbose = TRUE){

  # Changelog:
  # * 2020-May  , DM: First reasonable version
  # * 2020-11-06, DM: Added roxygen documentation
  #
  # ToDo:
  # * Get stimulation.intensity from @info[["LPOWER"]]
  # * add 'autoname' and other file handling parameters
  # * add background fitting functionality

  # Hidden parameters
  report_format <- "html"

  # Get name of the input object
  object_name <- deparse(substitute(object))

  # define new list object to safely ignore incompatible list elements
  data_set <- list()
  data_set_overhang <- list()

  # Test if object is a list. If not, create a list
  if (is.list(object)) {

    for (i in 1:length(object)) {

      if (class(object[[i]]) == "RLum.Analysis") {

        data_set[[length(data_set) + 1]] <- object[[i]]
      } else {

        element_name <- names(object)[i]

        if (element_name == "OSL_COMPONENTS") {

          warning("Input object already contains Step 1 results. Old results were overwritten")
        }else{

          data_set_overhang[[element_name]] <- object[[i]]
          if (!((element_name == "DECOMPOSITION")  || (element_name=="CORRECTION"))) {
            warning("Input object list element ", i, " is not of type 'RLum.Analysis' and was included in the fitting procedure, but was appended to the result list")}}}}

  } else {

    data_set <- list(object)
    warning("Input is not of type list, but output is of type list")}

  if (length(data_set) == 0) stop("Input object contains no RLum.Analysis data")

  # calc arithmetic mean curve
  if(verbose) cat("STEP 1.1 ----- Build global average curve from all CW-OSL curves -----\n")

  # measure computing time
  time.start <- Sys.time()

  global_curve <- sum_OSLcurves(data_set,
                                record_type = record_type,
                                output.plot = FALSE,
                                verbose = verbose)

  if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")


  # find components via fitting and F-statistics
  if(verbose) cat("STEP 1.2 ----- Perform multi-exponential curve fitting -----\n")

  time.start <- Sys.time()

  fit_data <- fit_OSLcurve(global_curve,
                         K.max = max_components,
                         F.threshold = F_threshold,
                         stimulation.intensity = stimulation_intensity,
                         stimulation.wavelength = stimulation_wavelength,
                         verbose = verbose,
                         output.complex = TRUE)

  # Add 'record_type' to the argument list
  fit_data$parameters$record_type <- record_type

  if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")

  if (report) {
    if(("rmarkdown" %in% rownames(installed.packages())) && ("kableExtra" %in% rownames(installed.packages()))) {

      if(verbose) cat("STEP 1.3 ----- Create report -----\n")

      time.start <- Sys.time()

      # the RMD script has to be located in the "/inst" folder of the project
      # then, it will be installed with the package
      try({

        # for test purposes only:
        #rmd_path <- "C:\\Users\\mitte\\Desktop\\R\\OSLdecomposition\\inst\\rmd\\report_Step1.Rmd"
        rmd_path <- system.file("rmd", "report_Step1.Rmd", package = "OSLdecomposition")

        output_path <- getwd()
        output_file <- paste0(output_path, "/", "report_Step1.", report_format)
        image_path <- paste0(output_path, "/report_figures/")

        if (!(dir.exists(image_path))) dir.create(image_path)
        cat("Save", toupper(image_format), "images to:", image_path, "\n")

        rmarkdown::render(rmd_path,
                          params = list(fit_data = fit_data,
                                        data_set = data_set,
                                        object_name = object_name,
                                        image_format = image_format,
                                        image_path = image_path),
                          output_file = output_file,
                          output_format = paste0(report_format,"_document"),
                          quiet = TRUE)

        cat("Save", toupper(report_format), "report to:", output_file, "\n")

        # ToDo: Replace the following try() outside the big try
        try({

          browseURL(output_file)
          cat("Open", toupper(report_format), "report in your systems standard browser\n")})

        if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")})

    } else {

      warning("Packages 'rmarkdown' and 'kableExtra' are needed to create reports. One or both are missing.")}}

  # Return fitted data
  object <- c(data_set, data_set_overhang, OSL_COMPONENTS = list(fit_data))
  invisible(object)
}
