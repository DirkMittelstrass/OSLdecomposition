#' Combine RLum OSL records to one global average curve
#'
#' The function adds up CW-OSL records saved in [RLum.Analysis-class] objects and
#' calculates the arithmetic mean signal from all records for each channel.
#' This is useful to create on global average curve with sufficient signal-to-noise ratio
#' for OSL components identification with [fit_OSLcurve] or to create one signal background
#' reference curve.
#'
#' Function supports currently just input objects of type `RLum.Analysis` from
#' the [Luminescence] package. Major overhaul necessary, see ToDo list
#' in the source code
#'
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' data set of one or multiple aliquots containing CW-OSL records
#'
#' @param record_type [character] (*with default*):
#' type of records which are selected from the input `object`,
#' for example: `"OSL"`,`"SGOSL"` or `"IRSL"`
#'
#' @param aliquot_selection [numeric] vector (*optional*):
#' vector specifying the indicies of elements (aliquots) of a list of [RLum.Analysis-class] objects
#' which shall be included
#'
#' @param offset_value [numeric] (*with default*):
#' signal offset (background) which will be substracted from each record
#'
#' @param verbose [logical] (*with default*):
#' enables console text output
#'
#' @param output.plot [logical] (*with default*):
#' returns a plot with all data points of all records and the average curve
#'
#' @param theme.set [ggplot2] object (*with default*):
#' sets the graphical theme of the output plot.
#' See [ggplot2::theme_bw] for available themes
#'
#' @param plot.first [logical] (*with default*):
#' plot includes additional drawing of first `record_type` record of first `object` list element
#'
#' @param title [character] (*with default*):
#' plot title. Set `title = "default"` for an automatically generated title.
#' Set `title = NULL` for no title
#'
#' @param filename [character] (*optional*):
#' file name or path to save the plot as image. If just a name is given, the image is
#' saved in the working directory. The image type is chosen by the file ending.
#' Allowed are `.pdf`, `.eps`, `.svg` (vector graphics) and `.jpg`, `.png`, `.bmp` (pixel graphics)
#' and more, see [ggplot2::ggsave] for details
#'
#'
#' @return
#' A [data.frame] of the average CW-OSL curve is returned, containing two columns: `$time` and `$signal`
#'
#'
#' @section Last updates:
#'
#' 2020-10-30, DM: Overworked plotting; Expanded roxygen documentation
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite the package the following way:
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' Automated identification and separation of quartz CW-OSL signal components with R, *in preparation*.
#'
#' @seealso [fit_OSLcurve], [RLum.OSL_correction], [RLum.OSL_global_fitting]
#'
#'
#' @examples
#'
#' # 'FB_10Gy' is a dose recovery test with the La Fontainebleau quartz
#' # measured in a lexsyg research with green LED stimulation
#' data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
#' data_set <- Luminescence::read_BIN2R(data_path, fastForward = TRUE)
#'
#  # Plot all data points and give average CW-OSL curve back
#' average_curve <- sum_OSLcurves(data_set)
#'
#' @md
#' @export
sum_OSLcurves <- function(
  object,
  record_type = "OSL",
  aliquot_selection = NULL,
  offset_value = 0,
  verbose = TRUE,
  output.plot = TRUE,
  theme.set = ggplot2::theme_classic(),
  plot.first = FALSE,
  title = "default",
  filename = NULL
){

  # Changelog:
  # * 2018-05-23, DM: first version
  # * 2019-03-15, DM: rewritten for ggplots and new data format, renamend into sum_OSLcurves
  # * 2020-08-30, DM: Overworked plotting; Expanded roxygen documentation
  #
  #  ToDo:
  # * add more options for record selection (e.g. dose)
  # * interpolate x-axis if it doesn't match
  # * accept list of data.frames as input objects
  # * test if first value is zero
  # * add legend to plot
  # * add info box with number of OSL curves to plot
  # * add 'hide.plot' parameter

  # prove if object is a list of aliquots or just a single aliquot
  if (is.list(object)) {

    # prove if aliquot selection is given. If not, take all aliquots of the data set
    if (is.null(aliquot_selection)
        || is.na(aliquot_selection)
        || length(aliquot_selection) > length(object)) {
      aliquot_selection <- c(1:length(object))
    }
  } else {

    object <- list(object)
    aliquot_selection <- c(1)
  }

  # Solve a devtools::check problem were it complains about ggplot2 syntax ("no visible binding for global variable")
  time <- signal <- record <- NULL


  ##============================================================================##
  # CALC MEAN CURVE
  ##============================================================================##

  mean.values <- c(0)
  first.values <- c(0)
  all.values <- NULL
  shortest_record_length <- Inf
  n <- 0
  old_x_axis <- c(0)
  new_x_axis <- c(0)


  for (j in aliquot_selection) {
    if (j < 1 || j > length(object)) {
      if (verbose) warning("Item ", j," is not a part of the data set. Item skipped")

    } else {
      if (class(object[[j]]) != "RLum.Analysis") {
        if (verbose) warning("Item ", j," is not of class RLum.Analysis. Item skipped")

      } else {######### LOOP BEGINS ########

        records <- object[[j]]@records
        for (i in c(1:length(records))) {

          if (records[[i]]@recordType == record_type) {

            mean.values <- mean.values + records[[i]]@data[,2] - offset_value

            # add record to one big data.frame for the multiple line plot
            all.values <- rbind(all.values, cbind(data.frame(records[[i]]@data), factor(paste0(j,"-",i))))

            # remember the length of the shortest record to cut the superposition curve later
            if (length(records[[i]]@data[,2]) < shortest_record_length) {
              shortest_record_length <- length(records[[i]]@data[,2])
            }

            # prove the x-axis for consistency
            new_x_axis <- records[[i]]@data[,1][1:shortest_record_length]
            if (length(old_x_axis)>1) {
              if (sum(new_x_axis) != sum(old_x_axis[1:shortest_record_length])) {
                if (verbose) warning("Record ", i, " in list item ", j, ": x-axis values differ from previous records")
              }
            }
            old_x_axis <- new_x_axis

            # count the records
            n <- n + 1

            # save the first curve separatly for displaying purposes
            if (n == 1) first.values <- mean.values
          }
        }
      }
    } ######### LOOP ENDS ########
  }

  # cut curve length to shortest record and calc mean
  mean.curve <- data.frame(time = new_x_axis, signal = mean.values[1:shortest_record_length] / n)

  if (verbose) cat(paste0("Built global average curve from arithmetic means from first ", shortest_record_length, " data points of all ", n, " ", record_type, " records\n"))

  ##============================================================================##
  # PLOT
  ##============================================================================##

  if (output.plot == TRUE) {

    # supress warnings
    options(warn = -1)

    #ggplot2::theme_set(ggplot2::theme_bw())
    ggplot2::theme_set(theme.set)

    text_format <- ggplot2::theme(axis.title = ggplot2::element_text(size = 8),
                                  plot.subtitle = ggplot2::element_text(size = 9, face = "bold"))

    #theme_set(theme_classic())

    colnames(all.values) <- c("time", "signal", "record")

    curve <- data.frame(mean.curve,
                        record = factor("mean"))

    first.curve <- data.frame(time = new_x_axis,
                              signal = first.values[1:shortest_record_length],
                              record = factor("first"))

    alpha.value <- round(10 / n, digits = 3)

    if (alpha.value > 0.4) alpha.value <- 0.4
    if (alpha.value < 0.01) alpha.value <- 0.01

    # plot a multiple line summary using the ggplot2 library
    p.lin <- ggplot2::ggplot(all.values, ggplot2::aes(x=time, y=signal, group=record)) +
      ggplot2::geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE, color = "black") +
      ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
      ggplot2::scale_y_continuous(limits = c(0, round(max(curve$signal) * 1.5)),
                                  labels = scales::label_number_auto()) +
      ggplot2::labs(subtitle = "CW-OSL", x = "Time (s)", y ="Signal (cts)") +
      text_format

    ######################## pseudoLM PLOT #########################################################

    # transform data
    P = 2*max(curve$time)

    LMdata <- all.values
    LMdata$signal <- LMdata$signal * (2 * LMdata$time / P)^0.5
    LMdata$time <- (2 * P * LMdata$time)^0.5

    LMcurve <- curve
    LMcurve$signal <- LMcurve$signal * (2 * LMcurve$time / P)^0.5
    LMcurve$time <- (2 * P * LMcurve$time)^0.5

    LMfirst.curve <- first.curve
    LMfirst.curve$signal <- LMfirst.curve$signal * (2 * LMfirst.curve$time / P)^0.5
    LMfirst.curve$time <- (2 * P * LMfirst.curve$time)^0.5


    p.LM <- ggplot2::ggplot(LMdata, ggplot2::aes(x=time, y=signal, group=record)) +
      ggplot2::geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE, color = "black") +
      ggplot2::scale_y_continuous(limits = c(0, round(max(LMcurve$signal) * 1.5)),
                                  labels = scales::label_number_auto()) +
      ggplot2::scale_x_continuous(labels = scales::label_number_auto()) +
      ggplot2::labs(subtitle = "pseudoLM-OSL", x = "Ramping time (s)", y ="Signal (cts)") +
      text_format

    if (plot.first) {

      p.lin <- p.lin + ggplot2::geom_line(data = first.curve, size = 0.5, color = "blue")
      p.LM <- p.LM + ggplot2::geom_line(data = LMfirst.curve, size = 0.5, color = "blue")}

    # plot global average curve
    p.lin <- p.lin + ggplot2::geom_line(data = curve, size = 0.5, color = "red")
    p.LM <- p.LM + ggplot2::geom_line(data = LMcurve, size = 0.5, color = "red")


    #plot title?
    if (!is.null(title) && title == "default") {

      title <- paste0("Data points and average curve of ", n, " ", record_type, " records")}

    # use grid function from gridExtra package to display linear and log plot side by side
    plot_object <- gridExtra::arrangeGrob(p.lin, p.LM, nrow = 1, top = title)

    # save plot as file
    if (!is.null(filename)) {
      try(ggplot2::ggsave(filename, plot = plot_object, units = "cm"), silent = FALSE)}

    # show plot
   gridExtra::grid.arrange(plot_object)

  }

  invisible(mean.curve)
}
