#' Title
#'
#' @param object
#' @param record_type
#' @param decay_rates
#' If not defined, the decay rates are taken from object$OSL_decomponents$case.tables[[object$OSL_decomponents$K.selected]].
#' Therefore, the number of components can be changed by altering object$OSL_decomponents$K.selected
#'
#' @param algorithm
#' @param error_calculation
#' @param background_fitting
#' @param report
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
RLum.OSL_decomposition <- function(
  object,
  record_type = "OSL",
  decay_rates = NULL,
  algorithm = "det+nls", # "det", "nls", "det+nls"
  error_calculation = "empiric",   # "poisson", "empiric", "nls", numeric value
  background_fitting = FALSE,
  report = TRUE,
  verbose = TRUE
){
  ### ToDo's
  # - read 'lambda.error' if available and transfer it to decompose_OSLcurve for better error calculation
  # - collect warnings from [decompose_OSLcurve()] and others and display them bundled with record-index at the end


  library(OSLdecomposition)
  library(Luminescence)

  # define new list object to safely ignore incompatible list elements
  data_set <- list()
  data_set_overhang <- list()
  data_is_lone_aliquot <- FALSE

  # Test if input object is a list
  if (is.list(object)) {

    for (i in 1:length(object)) {

      if (class(object[[i]]) == "RLum.Analysis") {

        data_set[[length(data_set) + 1]] <- object[[i]]

      } else {

        element_name <- names(object)[i]
        data_set_overhang[[element_name]] <- object[[i]]

        if (element_name != "OSL_COMPONENTS") {
          warning("List element ", i, " is not of type 'RLum.Analysis' and will not be included in fitting procedure, but will be appended to the result list")
        }
      }
    }

  } else {

    if (class(object) == "RLum.Analysis") {

      data_set <- list(object)
      data_is_lone_aliquot <- TRUE

    } else {
      stop("Input object is not a RLum.Analysis object nor a list of RLum.Analysis objects ")
    }
  }

  if (length(data_set) == 0) stop("Input object contains no RLum.Analysis data")

  ##### Get the decay rates ######

  global_curve <- NULL
  component_table <- NULL

  if (is.null(decay_rates)) {

    if ("OSL_COMPONENTS" %in% names(data_set_overhang)) {
      component_table <- data_set_overhang$OSL_COMPONENTS$case.tables[[data_set_overhang$OSL_COMPONENTS$K.selected]]
      global_curve <- data_set_overhang$OSL_COMPONENTS$curve

    } else {
       stop("Neither contains the input object an element $OSL_COMPONENTS (Step 1 results), nor is the argument 'decay_rates' defined")
    }


  } else if (is(decay_rates, "data.frame") && ("lambda" %in% colnames(decay_rates))) {
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
    error_calculation <- "nls"
  }

  ################################ STEP 2.1: Integration intervals ################################
  if (verbose) cat("STEP 2.1 ----- Define signal bin intervals -----\n")

  # Check if the integration intervals are given
  if (!("t.start" %in% colnames(component_table)) ||
      !("t.end" %in% colnames(component_table)) ||
      !("ch.start" %in% colnames(component_table)) ||
      !("ch.end" %in% colnames(component_table))) {

    time.start <- Sys.time()

    if (is.null(global_curve)) {
      cat(record_type, "curve template necessary but not given. Obtain and use global average curve:\n[sum_OSLcurves()]: ")
      global_curve <- sum_OSLcurves(data_set,
                                    record_type = record_type,
                                    output.plot = FALSE,
                                    verbose = verbose)

      if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")
      time.start <- Sys.time()
    }

    if (verbose) cat("Iterate minimum denominator determinant:\n[calc_OSLintervals()]: ")
    component_table <- calc_OSLintervals(component_table,
                                         global_curve,
                                         background.fitting = background_fitting,
                                         verbose = verbose)

    if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")

  } else {

    if(verbose) cat("Integration intervals are already given. Step skipped\n\n")
  }

  ################################ STEP 2.2: Decomposition  ################################

  if (verbose) cat("STEP 2.2 ----- Decompose each ", record_type," curve -----\n")
  if (verbose) cat("Calculate signal intensity n in each", record_type, " by '", algorithm,"' algorithm with ", error_calculation, " error estimation\n")
  if (verbose) cat("Table of input decay constants and signal bin intervals for [decompose_OSLcurve()]:\n")
  if (verbose) print(subset(component_table, select = c(name, lambda, t.start, t.end, ch.start, ch.end)))

  time.start <- Sys.time()

  # Build summarizing result table for each aliquot


  N_records <- 0
  for (j in 1:length(data_set)) {

    N_in_aliquot <- 1
    for (i in 1:length(data_set[[j]]@records)) {

      if (data_set[[j]]@records[[i]]@recordType == record_type) {

        decomp_table <- decompose_OSLcurve(data_set[[j]]@records[[i]]@data,
                                           component_table,
                                           algorithm = algorithm,
                                           background.fitting = background_fitting,
                                           verbose = FALSE)

        # Add the resulting data.frame to the info section of the records RLum object
        data_set[[j]]@records[[i]]@info[["COMPONENTS"]] <- decomp_table

        N_in_aliquot <- N_in_aliquot + 1
      }
    }
  }

  if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")


  ################################ STEP 2.3: Report  ################################



  invisible(data_set)

}
