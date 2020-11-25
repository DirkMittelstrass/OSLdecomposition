#' @title Render Report
#'
#' @description Centralised function to render the reports
#'
#' @param nature [character] (**required**): report nature
#'
#' @param dec_data
#'
#' @param fit_data
#'
#' @param data_set
#'
#' @param object_name
#'
#' @param output_dir [character] (*with default*): output folder of the report, it is
#' a temporary folder by default, passed to [rmarkdown::render]
#'
#' @param image_format [character] (*with default*): additional format for output image rendering
#'
#' @param image_path [character] (*width default*): optional image output folder
#'
#' @param report_format [character] (*with default*): report format
#'
#' @param verbose [logical] (*with default*): enables/disables terminal output.
#' If set to `FALSE` no browser window will be opened.
#'
#' @return returns a report in the chosen format and opens a browser window
#'
#' @md
#' @noRd
################################ REPORT  ################################
.render_report <- function(
  nature = "",
  cor_data = NULL,
  dec_data = NULL,
  fit_data = NULL,
  data_set,
  object_name,
  report_dir,
  image_format = NULL,
  report_format = "html",
  verbose = TRUE
){

# Pre-check ---------------------------------------------------------------
  pkgs <- c("rmarkdown", "knitr", "kableExtra")
  if(!all(pkg_test <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))){
    stop(paste0("Creating reports requires the package(s) '",
                  paste(pkgs[!pkg_test], collapse = "', '"),"'.\n",
                "To install this package run install.packages(c('",
                paste(pkgs[!pkg_test], collapse = "', '"),"')) in your R console."),
         call. = FALSE)
  }

  image_path <- NULL

  # set file path for the report
  if(is.null(report_dir)){

    output_dir <- tempdir()
    image_format <- NULL

  } else {

    # If the input was a file path, we need to remove the extension
    report_dir <- tools::file_path_sans_ext(report_dir)

    # Now add a "/" and normalize the path.
    # then it doesn't matter if the user put an "/" at the end of the directory
    output_dir <- normalizePath(paste0(report_dir, "/"))

    # Create Directory if it not already exists
    dir.create(output_dir, showWarnings = FALSE)

    # We may need also a path for the images
    if (!is.null(image_format)) {

      image_path <- paste0(output_dir,"report_figures\\")
      dir.create(image_path, showWarnings = FALSE)}

  }

# Set parameters ----------------------------------------------------------
 if(verbose) cat("This process can take a few minutes...\n")

  # the RMD script has to be located in the "/inst" folder of the project
  # then, it will be installed with the package

  ##select RMD-file
  rmd_ht <- c(
    correction = "report_Step1.Rmd",
    global_fitting = "report_Step1.Rmd",
    decomposition = "report_Step2.Rmd")

  if(!(nature[1] %in% names(rmd_ht))){
    stop("[.render_report()] report nature unknown, supported are: \n",
         paste0(" ",rmd_ht, " -> '", names(rmd_ht), "'\n"))

  }

  ##preset parameters remove NULL
  input_params <- list(
      dec_data = dec_data,
      fit_data = fit_data,
      data_set = data_set,
      object_name = object_name,
      image_format = image_format,
      image_path = image_path)

  input_params <- input_params[!sapply(input_params, is.null)]

  ##set timer
  time.start <- Sys.time()

  try({
    output <- rmarkdown::render(
      input =  system.file("rmd", rmd_ht[nature[1]], package = "OSLdecomposition", mustWork = TRUE),
      params = input_params,
      output_dir = output_dir,
      output_format = paste0(report_format, "_document"),
      quiet = TRUE,
      clean = TRUE
    )


    if(verbose) {

      if(!is.null(report_dir)) {
        cat("Saved", toupper(report_format), "report to:", output_dir, "\n")
        if(!is.null(image_format)) cat("Saved figures as", toupper(image_format),
                                       "files to:", normalizePath(paste0(output_dir,"/report_figures/")), "\n")}

      cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")}

  })

  try({
    if(verbose) {
      utils::browseURL(output)
      cat("Open", toupper(report_format), "report in the systems standard browser\n")
    }
  })

}
