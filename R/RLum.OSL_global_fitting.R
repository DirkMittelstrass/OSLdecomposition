#' Identifies CW-OSL signal components in RLum.Analysis data sets
#'
#' @last_changed 2020-04-19
#'
#' @param object
#' @param max_components
#' @param record_type
#' @param F_threshold
#' @param stimulation_intensity
#' @param stimulation_wavelength
#' @param report
#' @param verbose
#'
#' @return
#' @examples

RLum.OSL_global_fitting <- function(object,
                     max_components = 3,
                     record_type = "OSL",
                     F_threshold = 50,
                     stimulation_intensity = 30,
                     stimulation_wavelength = 470,
                     report = TRUE,
                     verbose = TRUE){

  ### ToDo:
  # - add 'autoname' argument
  # - add file name argument
  # - add file directory argument
  # - add background fitting functionality

  library(OSLdecomposition)
  library(Luminescence)

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

        if (names(object)[i]=="OSL_COMPONENTS") {

          warning("Input object already contains Step 1 results. Old results were overwritten")
        }else{

          data_set_overhang[[length(data_set_overhang) + 1]] <- object[[i]]
          warning("List element ", i, " is not of type 'RLum.Analysis' and was not included in fitting procedure")
        }
      }
    }

  } else {

    data_set <- list(object)
    warning("Input is not of type list, but output is of type list")
  }

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
                         background.fitting = FALSE,
                         verbose = verbose,
                         output.plot = FALSE,
                         output.complex = TRUE)

  # Add 'record_type' to the argument list
  fit_data$parameters$record.type <- record_type

  if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")

  if (report) {
    if(("rmarkdown" %in% rownames(installed.packages())) && ("kableExtra" %in% rownames(installed.packages()))) {

      if(verbose) cat("STEP 1.3 ----- Create report -----\n")
      library(rmarkdown)

      time.start <- Sys.time()

      # the RMD script has to be located in the "/inst" folder of the project
      # then, it will be installed with the package
      try({

        report_format <- "html"
        # for test purposes:
        rmd_path <- "C:\\Users\\mitte\\Desktop\\R\\OSLdecomposition\\inst\\rmd\\report_Step1.Rmd"
        #rmd_path <- system.file("rmd", "report_Step1.Rmd", package = "OSLdecomposition")
        output_file <- paste0(getwd(), "/", "report_Step1.", report_format)

        rmarkdown::render(rmd_path,
                          params = list(fit_data = fit_data,
                                        data_set = data_set,
                                        object_name = object_name),
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

      warning("Packages 'rmarkdown' and 'kableExtra' are needed to create reports. One or both are missing.")
    }
  }

  # Return fitted data
  object <- c(data_set, data_set_overhang, OSL_COMPONENTS = list(fit_data))
  invisible(object)
}
