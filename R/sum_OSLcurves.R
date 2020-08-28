#' Combine RLum OSL records to a global mean curve
#'
#' The function adds up RLum records of one specific type
#' it is useful to create background records and enable component fitting
#' at sets of records with low signal-to-noise ratio
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' data set of one or multiple aliquots containing containing luminescence measurements
#'
#' @param record_type [character] (*with default*):
#' record type ("OSL","TL","IRSL") which shall be summed up
#'
#' @param aliquot_selection [c] (*optional*):
#' vector specifying which items (aliquots) of a list of [RLum.Analysis-class] objects shall be combined
#'
#' @param offset_value [numeric] (*with default*):
#' signal offset (background) which will be substracted from each record before combining
#'
#' @param output.plot [logical] (*with default*):
#' returns a plot of the created superposition curve
#' @param plot.first
#' @param plot.global
#' @param verbose
#' @param title
#' @param return.plot
#'
#' @return
#' A [data.frame] is returned, containing the created mean curve
#'
#' @section Changelog:
#' * 2018-05-23, DM: first version
#' * 2019-03-15, DM: rewritten for ggplots and new data format, renamend into sum_OSLcurves
#'
#' @section ToDo:
#' * add more options for record selection (e.g. dose)
#' * interpolate x-axis if it doesn't match
#' * accept data.frames as input objects
#' * test if first value is zero
#' * add legend to plot
#' * add info box with number of OSL curves to plot
#'
#' @section Last changed: 2020-04-20
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [fit_OSLcurve], [RLum.OSL_correction], [RLum.OSL_global_fitting]
#'
#' @references
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' Automated identification and separation of quartz CW-OSL signal components with R, *in preparation*.
#'
#' @export
#'
#' @examples
#'

sum_OSLcurves <- function(
  object,
  record_type = "OSL",
  aliquot_selection = NULL,
  offset_value = 0,
  output.plot = TRUE,
  plot.first = TRUE,
  plot.global = TRUE,
  verbose = TRUE,
  title = "default",
  return.plot = FALSE
){

  ##============================================================================##
  # Function Setup and Input check
  ##============================================================================##

  library(Luminescence)


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

  if (verbose) writeLines(paste0("Built global average curve from arithmetic means from first ", shortest_record_length, " data points of all ", n, " ", record_type, " records"))

  ##============================================================================##
  # PLOT
  ##============================================================================##

  if (output.plot == TRUE) {

    library(gridExtra)
    library(ggplot2)
    theme_set(theme_bw())
    #theme_set(theme_classic())

    colnames(all.values) <- c("time", "signal", "record")

    curve <- data.frame(mean.curve,
                        record = factor("mean"))

    first.curve <- data.frame(time = new_x_axis,
                              signal = first.values[1:shortest_record_length],
                              record = factor("first"))

    alpha.value <- round(8 / n, digits = 3)

    if (alpha.value > 0.5) alpha.value <- 0.5
    if (alpha.value < 0.01) alpha.value <- 0.01

    # plot a multiple line summary using the ggplot2 library
    p.lin <- ggplot(all.values, aes(x=time, y=signal, group=record))+
      geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE) +
      ylim(0, round(max(curve$signal) * 1.5)) +
      ylab("signal (cts)") + xlab("time (s)") +
      theme(axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            axis.title = element_text(size = 12),
            legend.position = "none")



    ### plot log-log scale ###

    # set y-axis minimum
    log_limits <- c(1, max(all.values$signal))
    if (min(all.values$signal) > 0) {
      log_limits[1] <- 10^floor(log10(min(all.values$signal)))
    }

    breaks <- 10^(-10:10)
    minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

    p.log <- ggplot(all.values, aes(x=time, y=signal, group=record))+
      geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE) +
      scale_y_log10(limits = log_limits, position = "right",
                    breaks = breaks,
                    minor_breaks = minor_breaks) +
      scale_x_log10(limits = c(new_x_axis[2] - new_x_axis[1], max(new_x_axis)), position = "right",
                    breaks = breaks,
                    minor_breaks = minor_breaks) +
      ylab("signal (cts)") + xlab("time (s)") +
     # ylab("log(signal (cts))") + xlab("log(time (s))") +
      theme(axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            axis.title = element_text(size = 12),
            legend.position = "none")

    if (plot.first) {

      p.lin <- p.lin + geom_line(data = first.curve, size = 0.5, color = "red")
      p.log <- p.log + geom_line(data = first.curve, size = 0.5, color = "red")
    }
    if (plot.global) {

      p.lin <- p.lin + geom_line(data = curve, size = 1, color = "blue")
      p.log <- p.log + geom_line(data = curve, size = 1, color = "blue")
    }

    #plot title?
    if (!is.null(title)) {
      if (title == "default") {
        title <- paste0("Global average curve and signal scattering of ", n, " ", record_type, " records")
      }
    }

    # use grid function from gridExtra package to display linear and log plot side by side

    if (return.plot) {

      return(p.lin)
    } else {
      grid.arrange(p.lin, p.log, nrow = 1, top = title)
    }
  }

  invisible(mean.curve)
}
