#' Decomposes RLum.Data.Curves into its CW-OSL components
#'
#' @last_changed 2020-04-21
#'
#' @param object
#' @param record_type
#' @param decay_rates
#' If not defined, the decay rates are taken from object$OSL_decomponents$case.tables[[object$OSL_decomponents$K.selected]].
#' Therefore, the number of components can be changed by altering object$OSL_decomponents$K.selected
#'
#' @param error_calculation
#' @param report
#' @param verbose
#'
#' @return
#' @examples

RLum.OSL_decomposition <- function(
  object,
  record_type = "OSL",
  decay_rates = NULL,
  error_calculation = "empiric",   # "poisson", "empiric", "nls", numeric value
  report = TRUE,
  verbose = TRUE
){
  ### ToDo's
  # - read 'lambda.error' if available and transfer it to decompose_OSLcurve for better error calculation
  # - collect warnings from [decompose_OSLcurve()] and others and display them bundled with record-index at the end


  ### currently HIDDEN PARAMETERS ###
  algorithm <- "det+nls" # "det", "nls", "det+nls"
  background_fitting <- FALSE


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
        if (element_name == "DECOMPOSITION") {
          warning("Input object contained already Step 2 results (list element '$DECOMPOSITION'). Old results overwritten!")
        } else {

          data_set_overhang[[element_name]] <- object[[i]]
          if (element_name != "OSL_COMPONENTS") {
            warning("Input object list element ", i, " is not of type 'RLum.Analysis' and was included in the decomposition procedure, but was appended to the result list")
          }
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

  ## Little function to measure elapsed time ##

  time_duration <- function(start, end){
  v1 <- strptime(durations, format='%H:%M:%S')
  v1$hour * 3600 + v1$min * 60 + v1$sec
  return()
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
  if (verbose) cat("Calculate signal intensity n in each", record_type, " by '", algorithm,"' algorithm with", error_calculation, "error estimation\n")
  if (verbose) cat("Table of input decay constants and signal bin intervals for [decompose_OSLcurve()]:\n")
  if (verbose) print(subset(component_table, select = c(name, lambda, t.start, t.end, ch.start, ch.end)))

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
        N_records <- N_records + 1
      }
    }
  }

  if(verbose) cat("\nSuccessfully decomposed", N_records,"records\n")
  if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")

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
    if("rmarkdown" %in% rownames(installed.packages()) == TRUE) {

      if(verbose) cat("STEP 2.3 ----- Create report -----\n")
      library(rmarkdown)

      time.start <- Sys.time()

      # the RMD script has to be located in the "/inst" folder of the project
      # then, it will be installed with the package
      try({

        report_format <- "html"
        # for test purposes:
        rmd_path <- "C:\\Users\\mitte\\Desktop\\R\\OSLdecomposition\\inst\\rmd\\report_Step2.Rmd"
        #rmd_path <- system.file("rmd", "report_Step1.Rmd", package = "OSLdecomposition")
        output_file <- paste0(getwd(), "/", "report_Step2.", report_format)

        rmarkdown::render(rmd_path,
                          params = list(dec_data = dec_data, data_set = data_set),
                          output_file = output_file,
                          output_format = paste0(report_format,"_document"),
                          quiet = TRUE)

        cat("Save", toupper(report_format), "report to:", output_file, "\n")

        # ToDo: Replace the following try() outside the big try
        try({
          browseURL(output_file)
          cat("Open", toupper(report_format), "report in the systems standard browser\n")
        })


        if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
      })

    } else {

      warning("Package 'rmarkdown' is needed to create reports.")
    }
  }

  # Return decomposed data
  object <- c(data_set, data_set_overhang, DECOMPOSITION = list(dec_data))
  invisible(object)
}
