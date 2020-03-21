render_Report <- function(object,
                          output_type = "html"){

  # ToDo:
  # - add file name argument
  # - add file directory argument

  ### Input checks ###
  output_type <- tolower(output_type)

  #require(shiny)
  rmd_path <- system.file("rmd", "report_HTML.Rmd", package = "OSLdecomposition")

  file_path <- paste0(getwd(), "/", "report.", output_type)

  rmarkdown::render(rmd_path,
                    params = list(lum_data = object),
                    output_file = file_path,
                    output_format = paste0(output_type, "_document"),
                    quiet = TRUE)

  # Add Try() and maybe an argument
  browseURL(file_path)

}
