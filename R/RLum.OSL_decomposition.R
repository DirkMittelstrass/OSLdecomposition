#' Separate CW-OSL components in RLum.Analysis data sets
#'
#' Calculate the CW-OSL signal component intensites for each CW-OSL measurement
#' under the requirement that the decay rates are already given. The signal decompositon
#' process uses an analytical approach described in details in Mittelstrass (2019) and
#' Mittelstrass et al. (2021). This function processes just [RLum.Analysis-class] data sets created within
#' the [Luminescence-package] (Kreutzer et al. 2012).
#'
#' The workflow of this function is as following:
#'
#' \enumerate{
#'   \item {[optimise_OSLintervals]: Approximate the best suiting integration intervals. Use the global
#'   average curve as time axis template. If none is given, create one using [sum_OSLcurves]}
#'   \item [decompose_OSLcurve]: {Calculate component intensities for **all** `record_type` measurements.
#'   Use the `"det"` algorithm if a background correction was performed with [RLum.OSL_correction] or the
#'   `"det+nls"` algorithm if no background correction was performed. Use the `"empiric"` error estimation approach.}
#'   \item Create a `html` report to summarize the results (*optional*)
#'}
#'
#' Data sets must be formatted as [RLum.Analysis-class] objects and
#' should have been processed with [RLum.OSL_correction] and [RLum.OSL_global_fitting] beforehand.
#' Output objects are also [RLum.Analysis-class] objects and are
#' meant for equivalent dose determination with [Luminescence::analyse_SAR.CWOSL].
#'
#' If `report = TRUE`, a `html` report of the results is rendered by the [rmarkdown-package]
#' and saved in the working directory, which is usually the directory of the data file.
#' This report can be displayed, shared and published online without any requirements to
#' the OS or installed software. But an internet connection is needed to display
#' the *MathJax* encoded equations and special characters.
#' The *Rmarkdown* source code of the report can be found with the following command:
#'
#' `system.file("rmd", "report_Step2.Rmd", package = "OSLdecomposition")`
#'
#'
#'
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' Data set of one or multiple CW-OSL measured aliquots. The data set must either
#' contain a list element `"OSL_COMPONENTS"` or the parameter `decay_rates` must
#' be defined
#'
#' @param record_type [character] (*with default*):
#' Type of records, selected by the [RLum.Analysis-class] attribute `@recordType`.
#' Common are: `"OSL"`,`"SGOSL"` or `"IRSL"`
#'
#' @param K [numeric] (*with default*):
#' Number of components. Should be chosen after OSL component evaluation with
#'
#' @param decay_rates [numeric] vector or [data.frame] (*optional*):
#' User-defined component decay rates. If the
#'
#' @param report [logical] (*with default*):
#' Creates a `html` report, saves it in the working directory and opens it in your
#' standard browser. The report contains the results and further information
#' on the data processing
#'
#' @param image_format [character] (*with default*):
#' Image format of the automatically saved graphs if `report = TRUE`.
#' Allowed are `.pdf`, `.eps`, `.svg` (vector graphics), `.jpg`, `.png`, `.bmp` (pixel graphics)
#' and more, see [ggplot2::ggsave]. The images are saved in the working directory subfolder `/report_figures`.
#' Set `image_format = NULL` if no images shall be saved
#'
#' @param verbose [logical] (*with default*):
#' Enables console text output
#'
#'
#' @section Last updates:
#'
#' 2020-11-07, DM: Added roxygen documentation; Auto-switch between "det" and "det+nls" depending on background correction
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite the package the following way:
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' Automated identification and separation of quartz CW-OSL signal components with R, *in preparation*.
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
#' The input `object`, a [list] of [RLum.Analysis-class] objects is returned but with
#' a new list element `object[["DECOMPOSITION"]]`, containing:
#'
#' \itemize{
#'   \item `$decompositon.input` [data.frame]: Set of input components. Relevant is just the column `$lambda`
#'   \item `$results` [data.frame]: Overview table of decomposition
#'   \item `$parameters` [list]: Input and algorithm parameters
#' }
#'
#' The [RLum.Data.Curve-class] attribute `@info` of each CW-OSL record contains the
#' new entry `$COMPONENTS` with the curve-individual signal component parameters.
#' It can be read for example by:
#'
#'  `object[[i]]@records[[j]]@info[["COMPONENTS"]]`
#'
#' @examples
#'
#' #'FB_10Gy' is a dose recovery test with the La Fontainebleau quartz
#' # measured in a lexsyg research with green LED stimulation
#' data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
#' data_set <- Luminescence::read_BIN2R(data_path, fastForward = TRUE)
#'
#' \dontrun
#' # Separate components and create report
#' data_set_decomposed <- RLum.OSL_decomposition(
#' data_set, decay_rates = c(0.8, 0.05))
#' }
#'
#' @md
#' @export

RLum.OSL_decomposition <- function(
  object,
  record_type = "OSL",
  K = 3,
  decay_rates = NULL,
  report = TRUE,
  image_format = "pdf",
  verbose = TRUE
){
  ### Changelog
  # * 2020-May,   DM: First reasonable version
  # * 2020-11-07, DM: Added roxygen documentation; Auto-switch between "det" and "det+nls" depending on background correction
  #
  ### ToDo's
  # * read 'lambda.error' if available and transfer it to decompose_OSLcurve for better error calculation
  # * collect warnings from [decompose_OSLcurve()] and others and display them bundled with record-index at the end
  # * switch between "det" and "det+nls" automatically, depending on whether the data set was background corrected or not


  # hidden parameters
  background_fitting <- FALSE
  error_calculation <- "empiric" # "poisson", "empiric", "nls", numeric value


  # get name of the input object
  object_name <- deparse(substitute(object))

  # define new list object to safely ignore incompatible list elements
  data_set <- list()
  data_set_overhang <- list()

  # Test if input object is a list
  if (is.list(object)) {

    for (i in 1:length(object)) {

      if (class(object[[i]]) == "RLum.Analysis") {

        data_set[[length(data_set) + 1]] <- object[[i]]

      } else {

        element_name <- names(object)[i]
        if (element_name == "DECOMPOSITION") {
          warning("Input object contained already Step 2 results (list element '$DECOMPOSITION'). Old results overwritten!")
        } else {

          data_set_overhang[[element_name]] <- object[[i]]
          if (!((element_name == "OSL_COMPONENTS")  || (element_name=="CORRECTION"))) {
            warning("Input object list element ", i, " is not of type 'RLum.Analysis' and was included in the decomposition procedure, but was appended to the result list")}}}}

  } else {

    if (class(object) == "RLum.Analysis") {

      data_set <- list(object)
    } else {
      stop("Input object is not a RLum.Analysis object nor a list of RLum.Analysis objects ")
    }
  }

  if (length(data_set) == 0) stop("Input object contains no RLum.Analysis data")

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

    if ("OSL_COMPONENTS" %in% names(data_set_overhang)) {

      if (!is.numeric(K)) K <- data_set_overhang$OSL_COMPONENTS$K.selected

      component_table <- data_set_overhang$OSL_COMPONENTS$component.tables[[K]]
      global_curve <- data_set_overhang$OSL_COMPONENTS$curve

    } else {
       stop("Neither contains the input object an element $OSL_COMPONENTS (Step 1 results), nor is the argument 'decay_rates' defined")
    }


  } else if ((class(decay_rates) == "data.frame") && ("lambda" %in% colnames(decay_rates))) {
    component_table <- decay_rates

  } else if (is.numeric(decay_rates) && (length(decay_rates < 8))) {
    component_table <- data.frame(lambda = decay_rates)

  } else {
    stop("Neither is argument 'decay_rates' of class [data.frame] containing a column $lambda, nor is it a numeric vector with max. 7 elements")
  }

  # Check if the components are named
  if (!("name" %in% colnames(component_table))) {
    component_table$name <- paste0("Component ", 1:nrow(component_table))
  }

  # Whit NLS, just the nls() errors are available
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

  if(verbose)  if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")



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
  if(verbose) cat("\n\n")
  #if(verbose) cat("\       \t|\t Ln ", rep("     \t", times = nrow(test)), "    |\t Tn\n")
  #if(verbose) cat("Aliquot\t| ", paste0("n.", 1:nrow(component_table), "    \t"), "| ", paste0("n.", 1:nrow(component_table), "    \t"))

  for (j in 1:length(data_set)) {

    N_in_aliquot <- 1
    #if(verbose) cat(paste0("\n  #",j,"  \t"))
    if(verbose) cat(".")

    for (i in 1:length(data_set[[j]]@records)) {

      current_record <- data_set[[j]]@records[[i]]

      if (current_record@recordType == record_type) {

        decomp_table <- decompose_OSLcurve(current_record@data,
                                           component_table,
                                           algorithm = algorithm,
                                           background.fitting = background_fitting,
                                           verbose = FALSE)

        # Add the resulting data.frame to the info section of the records RLum object
        current_record@info[["COMPONENTS"]] <- decomp_table

        # A new line for the big table
        results <- rbind(results, data.frame(list.index = j,
                                             record.index = i,
                                             n = t(decomp_table$n),
                                             n.error = t(decomp_table$n.error),
                                             n.residual = t(decomp_table$n.residual),
                                             initial.signal = t(decomp_table$initial.signal),
                                             IRR_TIME = current_record@info[["IRR_TIME"]]))

        data_set[[j]]@records[[i]] <- current_record

        #if(verbose && (N_in_aliquot < 3)) cat("|", paste0(formatC(decomp_table$n, format = "e", digits = 2),"\t"))
        N_in_aliquot <- N_in_aliquot + 1
        N_records <- N_records + 1}}}

  if(verbose) cat("\nSuccessfully decomposed", N_records,"records\n")
  if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")

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
    if(("rmarkdown" %in% rownames(utils::installed.packages())) && ("kableExtra" %in% rownames(utils::installed.packages()))) {

      if(verbose) cat("STEP 2.3 ----- Create report -----\n")
      if(verbose) cat("This process can take a few minutes...\n")

      time.start <- Sys.time()

      # the RMD script has to be located in the "/inst" folder of the project
      # then, it will be installed with the package
      try({

        report_format <- "html"
        # for test purposes:
        #rmd_path <- "C:\\Users\\mitte\\Desktop\\R\\OSLdecomposition\\inst\\rmd\\report_Step2.Rmd"
        rmd_path <- system.file("rmd", "report_Step2.Rmd", package = "OSLdecomposition")

        output_path <- getwd()
        output_file <- paste0(output_path, "/", "report_Step2.", report_format)

        image_path <- NULL
        if (!is.null(image_format)) {

          image_path <- paste0(output_path, "/report_figures/")
          if (!(dir.exists(image_path))) dir.create(image_path)}

        rmarkdown::render(rmd_path,
                          params = list(dec_data = dec_data,
                                        data_set = data_set,
                                        object_name = object_name,
                                        image_format = image_format,
                                        image_path = image_path),
                          output_file = output_file,
                          output_format = paste0(report_format,"_document"),
                          quiet = TRUE)

        cat("Saved", toupper(report_format), "report to:", output_file, "\n")
        if (!is.null(image_format)) cat("Saved", toupper(image_format), "images to:", image_path, "\n")

        # ToDo: Replace the following try() outside the big try
        try({
          utils::browseURL(output_file)
          cat("Opened", toupper(report_format), "report in your systems standard browser\n")})

        if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")})

    } else {

      warning("Packages 'rmarkdown' and 'kableExtra' are needed to create reports. One or both are missing.")}}

  # Return decomposed data
  object <- c(data_set, data_set_overhang, DECOMPOSITION = list(dec_data))
  invisible(object)
}
