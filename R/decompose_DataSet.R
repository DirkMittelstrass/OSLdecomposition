decompose_DataSet <- function(
  data_set,
  components,
  aliquot_selection = NULL,
  record_type = "OSL",
  algorithm = "det+nls", # "det", "nls", "det+nls"
  error.calculation = "empiric",   # "poisson", "empiric", "nls", numeric value
  background.fitting = FALSE,
  verbose = TRUE
){
  ########## Input checks ###########

  # prove if object is a list of aliquots or just one single aliquot
  data_is_lone_aliquot <- FALSE
  if (is.list(data_set)) {

    # prove if a aliquot selection is given. If not, take all aliquots of the data set
    if (is.null(aliquot_selection)
        || is.na(aliquot_selection)
        || length(aliquot_selection) > length(data_set)) {
      aliquot_selection <- c(1:length(data_set))
    }
  } else {

    # Change Object to List, it will be transformed back before return
    data_is_lone_aliquot <- TRUE
    data_set <- list(data_set)
    aliquot_selection <- c(1)
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
      message("Warning: Item ", j," is not a part of the data set. Item skipped")
    } else {
      if (class(data_set[[j]]) != "RLum.Analysis") {
        message("Warning: Item ", j," is not of class RLum.Analysis. Item skipped")
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
