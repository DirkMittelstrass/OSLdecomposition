#' @title Separate CW-OSL components in RLum.Analysis data sets
#'
#' @description Calculates the CW-OSL signal component intensities for each CW-OSL measurement
#' under the requirement that the decay rates are already given. The signal decomposition
#' process uses an analytical approach described in detail in Mittelstrass (2019) and
#' Mittelstrass et al. (in preparation). This function processes [Luminescence::RLum.Analysis-class] data sets created within the [Luminescence::Luminescence-package] (Kreutzer et al. 2012).
#'
#' The workflow of this function is as follows:
#'
#' \enumerate{
#'   \item [optimise_OSLintervals]: Approximates the optimal integration intervals. Uses the global
#'   average curve as time axis template. If none global average curve is given, one is automatically created using [sum_OSLcurves].
#'   \item [decompose_OSLcurve]: Calculates component intensities for **all** `record_type` measurements.
#'   Uses the `"det"` algorithm if a background correction was performed with [RLum.OSL_correction] or the
#'   `"det+nls"` algorithm if no background correction was performed. For error estimation, the `"empiric"` approach is used.
#'   \item Creates a `html` report to summarize the results (*optional*).
#'}
#'
#' Data sets must be formatted as [Luminescence::RLum.Analysis-class] objects and
#' should have been processed with [RLum.OSL_correction] and [RLum.OSL_global_fitting] beforehand.
#' Output objects are also [Luminescence::RLum.Analysis-class] objects and are
#' meant for equivalent dose determination with [Luminescence::analyse_SAR.CWOSL].
#'
#' If `report = TRUE`, a `html` report of the results is rendered by the [rmarkdown::rmarkdown-package]
#' and saved in the working directory, which is usually the directory of the data file.
#' This report can be displayed, shared and published online without any requirements regarding
#' the operation system or installed software. However, an internet connection is needed to display
#' the *MathJax* encoded equations and special characters.
#' The *Rmarkdown* source code of the report can be found with the following command:
#'
#' `system.file("rmd", "report_Step2.Rmd", package = "OSLdecomposition")`
#'
#'
#' @param object [Luminescence::RLum.Analysis-class] or [list] of [Luminescence::RLum.Analysis-class]
#' (**required**):
#' Data set of one or multiple CW-OSL measured aliquots. The data set must either
#' contain a list element `$FITTING` or the parameter `decay_rates` must
#' be defined.
#'
#' @param record_type [character] (*with default*):
#' Type of records, selected by the [Luminescence::RLum.Analysis-class] attribute `@recordType`.
#' Common are: `"OSL"`,`"SGOSL"` or `"IRSL"`.
#'
#' @param K [numeric] (*with default*):
#' Number of components. Selects the according result table in the `$FITTING` list item of the data set `object`.
#'
#' @param decay_rates [numeric] vector or [data.frame] (*optional*):
#' User-defined component decay rates. If this parameter is defined, the parameter `K` will ignored.
#' If the input object is a [data.frame], then the decay rates must be stored in the column `$lambda`.
#'
#' @param report [logical] (*with default*):
#' Creates a `html` report, saves it in the `report_dir` directory.
#' The report contains the results and detailed information on the data processing.
#'
#' @param report_dir [character] (*optional*):
#' Path of output directory if `report = TRUE`. If `report_dir = NULL` (default),
#' a temporary folder is used which is deleted when the R session is closed.
#' File paths are also allowed as parameter, then a new directory named after the OSL data file
#' will be created.
#'
#' @param image_format [character] (*with default*):
#' Image format of the automatically saved graphs if `report = TRUE` and `report_dir` is set.
#' Allowed are `.pdf`, `.eps`, `.svg` (vector graphics), `.jpg`, `.png`, `.bmp` (pixel graphics)
#' and more, see [ggplot2::ggsave]. The images are saved in the `report_dir` subfolder `/report_figures`.
#' Set `image_format = NULL` if no images shall be saved.
#'
#' @param open_report [logical] (*with default*):
#' If set to `TRUE` a browser window displaying the report will be opened automatically.
#'
#' @param rmd_path [character] (*with default*):
#' **For advanced users:** File path to the [rmarkdown::rmarkdown-package] source code file of the report.
#' This allows to execute a manipulated version of the report.
#'
#' @param verbose [logical] (*with default*):
#' Enables console text output.
#'
#'
#' @section Last updates:
#'
#' 2026-03-02, DM:
#' * Function no longer crashes if record data contains no '@info$IRR_TIME' parameters.
#' * Default number of components 'K' if 'K' is not set is no longer 3. Instead 'K = length(decay_rates)' if 'decay_rates' are set, else 'K = $FITTING$K.selected'.
#' * Made pattern matching of 'record_type' with '@recordType' slot ready for Luminescence package 1.2
#' * Improved input data checks
#' * Existing RLum.OSL decomposition results are now removed with each execution of this function
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite the package the following way:
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' R package OSLdecomposition: Automated identification and separation of quartz CW-OSL signal components, *in preparation*.
#'
#' @seealso [RLum.OSL_global_fitting], [decompose_OSLcurve], [optimise_OSLintervals], [Luminescence::analyse_SAR.CWOSL]
#'
#' @references
#'
#' Kreutzer, S., Schmidt, C., Fuchs, M.C., Dietze, M., Fischer, M., Fuchs, M., 2012.
#' Introducing an R package for luminescence dating analysis. Ancient TL, 30 (1), 1-8.
#'
#' Mittelstraß, D., 2019. Decomposition of weak optically stimulated luminescence signals and
#' its application in retrospective dosimetry at quartz (Master thesis). TU Dresden, Dresden.
#'
#' @return
#'
#' The input `object`, a [list] of [Luminescence::RLum.Analysis-class] objects is returned but with
#' a new list element `object[["DECOMPOSITION"]]`, containing:
#'
#' \itemize{
#'   \item `$decompositon.input` [data.frame]: Set of input components. Relevant is just the column `$lambda`
#'   \item `$results` [data.frame]: Overview table of decomposition
#'   \item `$parameters` [list]: Input and algorithm parameters
#' }
#'
#' The [Luminescence::RLum.Data.Curve-class] attribute `@info` of each CW-OSL record contains the
#' new entry `$COMPONENTS` with the curve-individual signal component parameters.
#' It can be read for example by:
#'
#'  `object[[i]]@records[[j]]@info[["COMPONENTS"]]`
#'
#' @examples
#'
#' #'FB_10Gy' is a dose recovery test with the Fontainebleau quartz
#' # measured in a lexsyg research with green LED stimulation
#' data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
#' data_set <- Luminescence::read_BIN2R(data_path, fastForward = TRUE)
#'
#' # Separate components
#' data_set_decomposed <- RLum.OSL_decomposition(
#' data_set, decay_rates = c(0.8, 0.05))
#'
#' @md
#' @export

RLum.OSL_decomposition <- function(
  object,
  record_type = "OSL",
  K = NA,
  decay_rates = NULL,
  report = FALSE,
  report_dir = NULL,
  image_format = "pdf",
  open_report = TRUE,
  rmd_path = NULL,
  verbose = TRUE
){
  ### Changelog
  # * 2020-May,   DM: First reasonable version
  # * 2020-11-07, DM: Added roxygen documentation; Auto-switch between "det" and "det+nls" depending on background correction
  # * 2020-11-23, SK: Moved report call into utils.R
  # * 2021-02-15, DM: Added new parameter `rmd_path`
  # * 2022-05-02, DM: Added new parameter `open_report` to give control over automatic browser opening
  # * 2023-09-01, DM: Improved input checks to return more helpful messages
  # * 2026-02-26, DM: Function no longer crashes if record data contains no '@info$IRR_TIME' parameters.
  # * 2026-02-26, DM: Default number of components 'K' if 'K' is not set is no longer 3. Instead 'K = length(decay_rates)' if 'decay_rates' are set, else 'K = $FITTING$K.selected'.
  # * 2026-02-26, DM: Made pattern matching of 'record_type' with '@recordType' slot ready for Luminescence package 1.2
  # * 2026-03-02, DM: Improved input data checks
  # * 2026-03-02, DM: Existing RLum.OSL decomposition results are now removed with each execution of this function
  #
  #
  ### ToDo's
  # * read 'lambda.error' if available and transfer it to decompose_OSLcurve for better error calculation
  # * collect warnings from [decompose_OSLcurve()] and others and display them bundled with record-index at the end


  # hidden parameters
  background_fitting <- FALSE
  error_calculation <- "empiric" # "poisson", "empiric", "nls", numeric value
  verbose_performance <- FALSE

  # get name of the input object
  object_name <- deparse(substitute(object))

  # define new list object to safely ignore incompatible list elements
  data_set <- list()
  data_set_overhang <- list()

  # Test if object is a list of RLum.Analysis objects
  if (is.list(object)) {
    for (i in 1:length(object)) {

      if (inherits(object[[i]], "RLum.Analysis")) {

        data_set[[length(data_set) + 1]] <- object[[i]]
      } else {

        element_name <- names(object)[i]
        allowed <- c("Sequence.Header", "FITTING", "CORRECTION")
        not_allowed <- c("DECOMPOSITION")

        if (is.null(element_name)){
          cat("List element no.", i, "is not of type 'RLum.Analysis' and is removed from the data set.\n")

        } else if (element_name %in% not_allowed) {
          cat("Removed old list element", element_name, "to circumvent misleading results.\n")

        } else if (element_name %in% allowed) {
          data_set_overhang[[element_name]] <- object[[i]]

        } else{
          cat("List element", paste0("\"", element_name, "\""), "is not of type 'RLum.Analysis' and removed from the data set.\n")
        }
      }
    }

  } else if (inherits(object, "Risoe.BINfileData")) {
    stop(paste("Data is of type 'Risoe.BINfileData' instead of type 'RLum.Analysis'.",
               "Please apply the Luminescence package function Risoe.BINfileData2RLum.Analysis()",
               "to the data or ensure that read_BIN2R() has 'fastForward = TRUE' set."))

  } else if (inherits(object, "RLum.Analysis")) {
    data_set <- list(object)
    warning("Input was not of type list, but output is of type list.")

  } else {
    stop(paste("Invalid data type: Input object need to be a list of RLum.Analysis objects.",
               "Instead it is of type", class(object)[1]))
  }

  if (length(data_set) == 0) stop("Input data contains no RLum.Analysis objects. Please check if the data import was done correctly.")

  ################### Find out, which algorithm to use ###################

  # default
  algorithm <- "det+nls" # "det", "nls", "det+nls"

  # Was there a background correction?
  if ("CORRECTION" %in% names(data_set_overhang)){
    if ("background_curve" %in% names(data_set_overhang[["CORRECTION"]])
        | (data_set_overhang[["CORRECTION"]]$parameters$subtract_offset > 0)) {

      algorithm <- "det" }}


  ################### Get the decay rates ###################

  global_curve <- NULL
  component_table <- NULL

  if (is.null(decay_rates)) {

    if ("FITTING" %in% names(data_set_overhang)) {

      if (!is.numeric(K)) K <- data_set_overhang$FITTING$K.selected

      component_table <- data_set_overhang$FITTING$component.tables[[K]]
      global_curve <- data_set_overhang$FITTING$curve

    } else {
       stop("Neither contains the input object an element $FITTING (Step 1 results), nor is the argument 'decay_rates' defined")
    }


  } else if ((inherits(decay_rates, "data.frame")) && ("lambda" %in% colnames(decay_rates))) {
    component_table <- decay_rates

  } else if (is.numeric(decay_rates) && (length(decay_rates < 8))) {
    component_table <- data.frame(lambda = decay_rates)

  } else {
    stop("Neither is argument 'decay_rates' of class [data.frame] containing a column $lambda, nor is it a numeric vector with max. 7 elements")
  }

  # If K is still 'NA'
  if (!is.numeric(K)) K <- nrow(component_table)


  # Check if the components are named
  if (!("name" %in% colnames(component_table))) {
    component_table$name <- paste0("Component ", 1:nrow(component_table))
  }

  # Witt NLS, just the nls() errors are available
  if ((algorithm == "nls") &! (error_calculation == "nls")) {
    if (verbose) warning("When algorithm 'nls' is chosen, error.calculation must be also 'nls'. Argument changed to error.calculation='nls'")
    error_calculation <- "nls"}


  ################################ STEP 2.1: Integration intervals ################################
  if (verbose) cat("STEP 2.1 ----- Define signal bin intervals -----\n")

  time.start <- Sys.time()

  if (is.null(global_curve)) {
    cat(record_type, "curve template necessary but not given. Obtain and use global average curve:\n[sum_OSLcurves()]: ")
    global_curve <- sum_OSLcurves(data_set,
                                  record_type = record_type,
                                  output.plot = FALSE,
                                  verbose = verbose)

    if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")
    time.start <- Sys.time()}

  if (verbose) cat("Find intervals with lowest component cross correlation by maximising the denominator determinant in Cramers rule:\n")
  component_table <- optimise_OSLintervals(component_table,
                                           global_curve,
                                           background.component = background_fitting,
                                           verbose = verbose)

  if (verbose_performance) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n")
  if (verbose) cat("\n")

  ################################ STEP 2.2: Decomposition  ################################

  if (verbose) cat("STEP 2.2 ----- Decompose each ", record_type," curve -----\n")
  if (verbose) cat("Calculate signal intensity n in each", record_type, " by '", algorithm,"' algorithm with", error_calculation, "error estimation\n")
  if (verbose) cat("Table of input decay constants and signal bin intervals for [decompose_OSLcurve()]:\n")
  if (verbose) print(subset(component_table, select = c("name", "lambda", "t.start", "t.end", "ch.start", "ch.end")))

  time.start <- Sys.time()

  # Build one big table containing all results
  # So we can easily filter out any statistical aspect we want later
  results <- data.frame(NULL)

  N_records <- 0
  if(verbose) cat("\nProcess aliquots: ")

  for (j in 1:length(data_set)) {

    N_in_aliquot <- 1
    if(verbose) cat("o")

    for (i in 1:length(data_set[[j]]@records)) {

      current_record <- data_set[[j]]@records[[i]]

      if (check_RLum.Data(current_record, record_type, verbose = FALSE)) {

        decomp_table <- decompose_OSLcurve(current_record@data,
                                           component_table,
                                           algorithm = algorithm,
                                           background.fitting = background_fitting,
                                           verbose = FALSE)

        # Add the resulting data.frame to the info section of the records RLum object
        current_record@info[["COMPONENTS"]] <- decomp_table

        # BIN and BINX should contain the irradiation time, necessary to build
        # a dose-response curve. However, XSYG files do not.
        irr_time <- NA
        if ("IRR_TIME" %in% names(current_record@info))
          irr_time <- current_record@info[["IRR_TIME"]]

        # A new line for the big table
        results <- rbind(results, data.frame(list.index = j,
                                             record.index = i,
                                             n = t(decomp_table$n),
                                             n.error = t(decomp_table$n.error),
                                             n.residual = t(decomp_table$n.residual),
                                             initial.signal = t(decomp_table$initial.signal),
                                             IRR_TIME = irr_time))

        data_set[[j]]@records[[i]] <- current_record
        N_in_aliquot <- N_in_aliquot + 1
        N_records <- N_records + 1}}}

  if (verbose) cat("\nSuccessfully decomposed", N_records,"records\n")
  if (verbose_performance) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n")
  if (verbose) cat("\n")

  # Build overview list
  dec_data <- list(parameters = list(record_type = record_type,
                                     decay_rates = decay_rates,
                                     algorithm = algorithm,
                                     error_calculation = error_calculation,
                                     background_fitting = background_fitting),
                   decompositon.input = component_table,
                   results = results)

  ################################ STEP 2.3: Report  ################################
  if (report) {

    if(verbose) cat("STEP 2.3 ----- Create report -----\n")

    # Rebuild needed to call .render_report (won't work from the global environment)
     .render_report(
        nature = "decomposition",
        dec_data = dec_data,
        data_set = data_set,
        object_name = object_name,
        image_format = image_format,
        report_dir = report_dir,
        open_report = open_report,
        rmd_path = rmd_path,
        verbose = verbose)}

# Return output -----------------------------------------------------------
  object <- c(data_set, data_set_overhang, DECOMPOSITION = list(dec_data))
  invisible(object)
}
