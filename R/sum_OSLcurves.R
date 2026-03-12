#' @title Combine lists of records to one global average curve
#'
#' @description This function adds up all records of the same type and
#' calculates the arithmetic mean signals.
#' This is useful to create global average curve with sufficient signal-to-noise ratio
#' for OSL components identification with [fit_OSLcurve] or to create a signal background
#' reference curve.
#'
#'
#' @param object [Luminescence::RLum.Analysis-class], or [list] of [Luminescence::RLum.Analysis-class], [Luminescence::RLum.Data.Curve-class] or [data.frame] (**required**):
#' Data set of CW-OSL records.
#'
#' @param record_type [character] (*with default*):
#' Type of records which are selected from the input `object` if it is a RLum object,
#' for example: `"OSL"`,`"SGOSL"` or `"IRSL"`. Does not apply for lists of [data.frame].
#'
#' @param selection [numeric] vector (*optional*):
#' Vector specifying the indices of elements (aliquots) of a list of [Luminescence::RLum.Analysis-class] objects
#' which shall be included.
#'
#' @param Y_offset [numeric] (*with default*):
#' Signal offset (background) which will be subtracted from each record.
#'
#' @param verbose [logical] (*with default*):
#' Enables console text output.
#'
#' @param output.plot [logical] (*with default*):
#' Returns a linear plot as well as a pseudoLM-OSL plot with all data points of all records and the average curve.
#'
#' @param theme.set [ggplot2::ggplot2-package] object (*with default*):
#' sets the graphical theme of the output plot.
#' See [ggplot2::theme_bw] for available themes
#'
#' @param plot.first [logical] (*with default*):
#' Plot includes additional drawing of first `record_type` record of first `object` list element.
#'
#' @param title [character] (*optional*):
#' Plot title. Set `title = "auto"` for an automatically generated title.
#'
#' @param filename [character] (*optional*):
#' File name or path to save the plot as image. If just a file name is given, the image is
#' saved in the working directory. The image type is chosen by the file ending. Both, vector images
#' as well as pixel images are possible. Allowed are `.pdf`, `.eps`, `.svg` (vector graphics),
#' `.jpg`, `.png`, `.bmp` (pixel graphics) and more, see [ggplot2::ggsave].
#'
#'
#' @return
#' A [data.frame] of the average signal curve is returned, containing two columns: `$time` and `$signal`.
#'
#'
#' @section Last updates:
#'
#' 2026-02-25, DM: Revised and refactored whole code
#' * Function won't accept data sets with varying x-axes and measurement lengths anymore. Now, ONLY records with x-axes IDENTICAL to those of the first record ARE INCLUDED.
#' * Accepts list of data.frames as input objects
#' * Changed naming and default values of some secondary arguments
#' * Faster calculation
#' * Plotting is less prone of errors and warnings
#'
#' @author
#' Dirk Mittelstraß, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite this package, including its version number. Enter the command `citation("OSLdecomposition")` to generate the correct reference.
#'
#' @seealso [fit_OSLcurve], [RLum.OSL_correction], [RLum.OSL_global_fitting]
#'
#'
#' @examples
#'
#' # 'FB_10Gy' is a dose recovery test with the Fontainebleau quartz
#' # measured in a lexsyg research with green LED stimulation
#' data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
#' data_set <- Luminescence::read_BIN2R(data_path, fastForward = TRUE)
#'
#' # Give average CW-OSL curve back
#' average_curve <- sum_OSLcurves(data_set)
#'
#' @md
#' @export
sum_OSLcurves <- function(
  object,
  record_type = "OSL",
  selection = NULL,
  Y_offset = 0,
  verbose = TRUE,
  output.plot = FALSE,
  theme.set = ggplot2::theme_classic(),
  plot.first = FALSE,
  title = NULL, # "auto"
  filename = NULL
){

  # Changelog:
  # * 2018-05-23, DM: first version
  # * 2019-03-15, DM: rewritten for ggplots and new data format, renamed into sum_OSLcurves
  # * 2020-08-30, DM: Overworked plotting; Expanded roxygen documentation
  # * 2026-02-23, DM: Accepts list of data.frames as input objects
  # * 2026-02-25, DM: Refactored whole code: Faster and less prone of errors and warnings
  # * 2026-02-25, DM: Changed naming and default values of some secondary arguments
  # * 2026-02-25, DM: Function won't accept data sets with varying x-axes and measurement lengths anymore.
  #                   Now, ONLY records with x-axes IDENTICAL to those of the first record ARE INCLUDED.

  #  ToDo:
  # * Like in the old version, allow curves of different length. Shorten to the shortest curve
  # * Add options for plotting, like for plot_OSLcurve
  # * Faster plotting
  # * add more options for record selection (e.g. dose)
  # * add legend to plot
  # * add info box with number of OSL curves to plot
  # * add 'hide.plot' parameter

  # prove if object is a list of aliquots or just a single aliquot
  if (is.list(object)) {

    # prove if aliquot selection is given. If not, take all aliquots of the data set
    if (is.null(selection)
        || is.na(selection)
        || length(selection) > length(object)) {
      selection <- c(1:length(object))
    }
  } else {

    object <- list(object)
    selection <- c(1)
  }
  check_verbose <- FALSE

  # ---------------------------------------------------------------- #
  # ------ Step 1: Build big data.frame with all curves in it ------ #
  # ---------------------------------------------------------------- #

  skipped <- 0
  big_df <- data.frame()
  master_curve <- NULL
  is_first_record <- TRUE

  for(obj in object){

    # Case 1: List entry is a RLum.Analysis object
    if (inherits(obj, "RLum.Analysis")) {

      for (record in obj@records) {
        if (check_RLum.Data(record, record_type = record_type,
                            curve_template = master_curve, verbose = check_verbose)) {

          if (is_first_record) {

            big_df <- data.frame(time = record@data[, 1], signal = record@data[, 2])
            master_curve <- record
            is_first_record <- FALSE

          }else{
            big_df <- cbind(big_df, record@data[, 2])}

        }else{
          skipped <- skipped + 1}
      }

      # Case 2: List entry is a RLum.Data.Curve object
    } else if (check_RLum.Data(obj, record_type = record_type, curve_template = master_curve, verbose = check_verbose)){

      if (is_first_record) {

        big_df <- data.frame(time = obj@data[, 1], signal = obj@data[, 2])
        master_curve <- obj
        is_first_record <- FALSE
      }else{

        big_df <- cbind(big_df, obj@data[, 2])}


      # Case 3: List entry is a data.frame
    } else if (inherits(obj, "data.frame") && ncol(obj) >= 2
               && is.numeric(obj[, 1]) && is.numeric(obj[, 2])) {

      if(is_first_record){

        big_df <- data.frame(time = obj[, 1], signal = obj[, 2])
        is_first_record <- FALSE

      } else if (all(obj[, 1] == big_df[, 1])){ # Are the x-axes identical?

        big_df <- cbind(big_df, obj[, 2])

      } else {
        skipped <- skipped + 1}


    } else {
      skipped <- skipped + 1}
  }

  if (verbose && skipped > 0) cat(skipped, "records were not suitable.\n")

  # ---------------------------------------------------------------- #
  # ------ Step 2: Build big data.frame with all curves in it ------ #
  # ---------------------------------------------------------------- #

  no_records <- ncol(big_df) - 1
  if (no_records < 1) {
    if (verbose) cat("No suitable record found in the input object.\n")
    return(NULL)

  } else if(no_records == 1) {

    if (verbose) cat("Only one suitable record found in the input object.\n")
    mean_curve <- first_curve <- big_df - Y_offset

  } else {
    if (verbose) cat("Built arithmetic mean curve from", no_records, record_type, "records\n")
     mean_values <- rowMeans(big_df[, -1]) - Y_offset
     mean_curve <- data.frame(time = big_df$time, signal = mean_values)
     first_curve <- big_df[, c(1,2)] #data.frame(time = big_df$time, big_df$signal)
  }

  # Leave function if no plotting is wanted
  if (!output.plot) return(mean_curve)

  # --------------------------------- #
  # ------ Step 3: Do the plot ------ #
  # --------------------------------- #

  # Formatting text elements and overall appearance
  ggplot2::theme_set(theme.set)
  text_format <-
    ggplot2::theme(axis.title = ggplot2::element_text(size = 8),
                   plot.title = ggplot2::element_text(size = 9, face = "bold"))

  # Make points lighter with increasing number of points
  alpha.value <- round(10 / no_records, digits = 3)
  if (alpha.value > 0.4) alpha.value <- 0.4
  if (alpha.value < 0.01) alpha.value <- 0.01

  ## The following code is optimized for performance and was suggested by ChatGPT ##

  # Extract all signal columns as matrix
  signal_matrix <- as.matrix(big_df[, -1])

  # Create long vectors WITHOUT pivoting (memory efficient)
  x_vals <- rep(big_df$time, times = ncol(signal_matrix))
  y_vals <- as.vector(signal_matrix)

  if (length(y_vals) > 42000)
    if (verbose) cat("Plotting may take a while because of many data points...\n")

  # ---------- Linear plot --------------
  p.lin <- ggplot2::ggplot() +

    # Grey/black data cloud
    ggplot2::geom_point(
      ggplot2::aes(x = x_vals, y = y_vals),
      color = "black", alpha = alpha.value, size = 0.5) +

    # Global average curve
    ggplot2::geom_line(
      ggplot2::aes(x = mean_curve$time, y = mean_curve$signal),
      color = "red", linewidth = 0.5) +

    # Scaling and formatting
    ggplot2::coord_cartesian(ylim = c(0, round(max(mean_curve$signal) * 1.5))) +
    ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
    ggplot2::labs(title = "CW-OSL", x = "Time (s)", y = "Signal (cts)") +
    text_format

  # ---------- LM-OSL plot --------------

  # Double the stimulation time
  P = 2*max(mean_curve$time)


  # Transform data points
  y_LM <- y_vals * sqrt(2 * x_vals / P)
  x_LM <- sqrt(2 * P * x_vals)

  # Transform mean curve
  LMcurve <- mean_curve
  LMcurve$signal <- LMcurve$signal * sqrt(2 * LMcurve$time / P)
  LMcurve$time <- sqrt(2 * P * LMcurve$time)

  p.LM <- ggplot2::ggplot() +

    # Grey/black data cloud
    ggplot2::geom_point(
      ggplot2::aes(x = x_LM, y = y_LM),
      color = "black", alpha = alpha.value, size = 0.5) +

    # Global average curve
    ggplot2::geom_line(
      ggplot2::aes(x = LMcurve$time, y = LMcurve$signal),
      color = "red", linewidth = 0.5) +

    # Scaling and formatting
    ggplot2::coord_cartesian(ylim = c(0, round(max(LMcurve$signal) * 1.5))) +
    ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
    ggplot2::scale_y_continuous(labels = scales::label_number_auto()) +
    ggplot2::labs(title = "pseudoLM-OSL", x = "Ramping time (s)", y = "Signal (cts)") +
    text_format

  # ---------- Plot first curve --------------
  if (plot.first) {
    p.lin <- p.lin +
      ggplot2::geom_line(
        ggplot2::aes(x = first_curve$time, y = first_curve$signal),
        size = 0.5, color = "blue")

    LMfirst <- first_curve
    LMfirst$signal <- LMfirst$signal * sqrt(2 * LMfirst$time / P)
    LMfirst$time <- sqrt(2 * P * LMfirst$time)

    p.LM <- p.LM +
      ggplot2::geom_line(
        ggplot2::aes(x = LMfirst$time, y = LMfirst$signal),
        size = 0.5, color = "blue")
    }

  # ---------- Arrange plots --------------

  if (!is.null(title) && title == "auto") {

    title <- paste0("Data points and average curve of ", no_records, " ", record_type, " records")}

  # use grid function from gridExtra package to display linear and log plot side by side
  plot_object <- gridExtra::arrangeGrob(p.lin, p.LM, nrow = 1, top = title)

  # ---------- Output --------------

  # save plot as file
  if (!is.null(filename)) {
    try(suppressMessages(ggplot2::ggsave(filename, plot = plot_object, units = "cm")),
        silent = FALSE)}

  # show plot
  gridExtra::grid.arrange(plot_object)

  invisible(mean_curve)
}
