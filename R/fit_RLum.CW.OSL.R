#' Determines CW-OSL signal components from RLum.Analysis data sets
#'
#'
#'
#'
#'

fit_RLum.CW.OSL <- function(object,
                     max_components = 3,
                     record_type = "OSL",
                     F_threshold = 50,
                     stimulation_intensity = 30,
                     stimulation_wavelength = 470,
                     report = TRUE,
                     report_format = "html",
                     plot = FALSE,
                     verbose = TRUE){

  # ToDo:
  # - add 'autoname' argument
  # - add file name argument
  # - add file directory argument
  # - export not-html reports



  # define new object list to rule out incompatible list elements
  data_set <- list()
  data_set_overhang <- list()

  # Test if object is a list. If not, create a list
  if (is.list(object)) {

    for (i in 1:length(object)) {

      if (class(object[[i]]) == "RLum.Analysis") {

        data_set[[length(data_set) + 1]] <- object[[i]]
      } else {

        if (names(lum_data)[i]=="COMPONENTS") {

          warning("Data set already contains step 1 results. They will be ovewritten")
        }else{

          data_set_overhang[[length(data_set) + 1]] <- object[[i]]
          warning("List element ", i, " is not of type 'RLum.Analysis'")
        }
      }
    }

  } else {

    data_set <- list(object)
  }



  # calc arithmetic mean curve
  if(verbose) cat("STEP 1.1 ----- Build arithmetic mean curve from all CW-OSL curves -----\n")

  # measure computing time
  time.start <- Sys.time()

  global_curve <- sum_OSLcurves(data_set,
                                record_type = record_type,
                                output.plot = plot,
                                plot.first = FALSE,
                                plot.global = plot,
                                title = NULL,
                                verbose = verbose)

  if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")


  # find components via fitting and F-statistics
  if(verbose) cat("STEP 1.2 ----- Perform multi-exponential curve fitting -----\n")

  time.start <- Sys.time()

  C.list <- fit_OSLcurve(global_curve,
                         K.max = max_components,
                         F.threshold = F_threshold,
                         stimulation.intensity = stimulation_intensity,
                         stimulation.wavelength = stimulation_wavelength,
                         applied.time.cut = TRUE,
                         background.fitting = FALSE,
                         verbose = verbose,
                         output.plot = plot)

  if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")

  if (report) {
    if("rmarkdown" %in% rownames(installed.packages()) == TRUE) {

      if(verbose) cat("STEP 1.3 ----- Create analysis report -----\n")
      library(rmarkdown)

      time.start <- Sys.time()

      # the RMD script has to be located in the "/inst" folder of the project
      # then, it will be installed with the package
      rmd_path <- system.file("rmd", "report_Step1.Rmd", package = "OSLdecomposition")
      output_file <- paste0(getwd(), "/", "report_Step1.", report_format)

      rmarkdown::render(rmd_path,
                        params = list(data_set = data_set),
                        output_file = output_file,
                        output_format = paste0(report_format,"_document"),
                        quiet = TRUE)


      # Add Try() and a warning, if it fails
      browseURL(file_path)

      if(verbose) cat("(time needed:", round(as.numeric(Sys.time() - time.start), digits = 2),"s)\n\n")

    } else {

      warning("Package 'rmarkdown' is needed to create reports.")
    }
  }

  # Print results

  object <- c(data_set, data_set_overhang, COMPONENTS = list(C.list))

  invisible(object)

}
