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

        data_set_overhang[[length(data_set_overhang) + 1]] <- object[[i]]

        if (names(object)[i] != "COMPONENTS") {
          warning("List element ", i, " is not of type 'RLum.Analysis' and will not be included in fitting procedure")
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

  ### Get the decay rates ###
  if (is.null(decay_rates)) {
 ############## CONTINUE HERE ############################
  }

  # calc arithmetic mean curve
  if(verbose) cat("STEP 2.1 ----- Build arithmetic mean curve from all CW-OSL curves -----\n")

  # measure computing time
  time.start <- Sys.time()

  ##########################

  # prove if object is a list of aliquots or just one single aliquot
  data_is_lone_aliquot <- FALSE
  if (is.list(data_set)) {

  } else {

    # Change Object to List, it will be transformed back before return
    data_is_lone_aliquot <- TRUE
    data_set <- list(data_set)
  }

  # Whit NLS, just the nls() errors are available
  if ((algorithm == "nls") &! (error.calculation == "nls")) {
    if (verbose) warning("When algorithm 'nls' is chosen, error.calculation must be also 'nls'. Argument changed to error.calculation='nls'")
    error.calculation <- "nls"
  }


  # are the integration intervals given?
  if (!("t.start" %in% colnames(components)) ||
      !("t.end" %in% colnames(components)) ||
      !("ch.start" %in% colnames(components)) ||
      !("ch.end" %in% colnames(components))) {
    if (verbose) warning("Integration intervals not provided. calc_OSLintervals() executed")
    components <- calc_OSLintervals(components,
                                    curve,
                                    background.fitting = background.fitting,
                                    verbose = verbose)

  }

  ########## Main loop #############

  for (j in aliquot_selection) {
    if (j < 1 || j > length(data_set)) {
      message("\nWarning: Item ", j," is not a part of the data set. Item skipped")
    } else {
      if (class(data_set[[j]]) != "RLum.Analysis") {
        message("\nWarning: Item ", j," is not of class RLum.Analysis. Item skipped")
      } else {


        if (verbose) writeLines(" ", sep = "\n")
        if (verbose) writeLines(c("Decompose aliquot: ", j), sep = "\t")
        #if (verbose) writeLines(" ", sep = "\n")

        records <- data_set[[j]]@records


        ########## Perform decomposition ###########

        for (i in c(1:length(records))) {
          if (records[[i]]@recordType == record_type) {


            decomp_table <- decompose_OSLcurve(records[[i]]@data,
                                                components,
                                                algorithm = algorithm,
                                                background.fitting = background.fitting,
                                                verbose = FALSE)

            data_set[[j]]@records[[i]]@info[["COMPONENTS"]] <- decomp_table

          }
        }

      }
    }
  }

  invisible(data_set)

}
