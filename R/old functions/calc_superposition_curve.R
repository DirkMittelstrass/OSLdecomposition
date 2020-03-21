#' Combine RLum records
#'
#' The function adds up RLum records of one specific type
#' it is useful to create background records and enable component fitting
#' at sets of records with low signal-to-noise ratio#'
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' data set of one or multiple aliquots containing containing luminescence measurements
#'
#' @param record_type [character] (*with default*):
#' record type ("OSL","TL","IRSL") which shall be summed up
#'
#' @param mean.values [logical] (*with default*):
#' if TRUE, an arithmetic mean curve is given back; if FALSE, a sum curve is given back
#'
#' @param aliquot_selection [c] (*optional*):
#' vector specifying which items (aliquots) of a list of [RLum.Analysis-class] objects shall be combined
#'
#' @param offset_value [numeric] (*with default*):
#' signal offset (background) which will be substracted from each record before combining
#'
#' @param output.plot [logical] (*with default*):
#' returns a plot of the created superposition curve
#'
#' @return
#' A [RLum.Analysis-class] is returned, containing the created superposition curve
#'
#' @section Changelog:
#' * 2018-05-23, DM: first version
#'
#' @section ToDo:
#' * calculate standard deviation for each data point and include error in plot
#' * add more options for record selection (e.g. dose type)
#' * interpolate x-axis if it doesn't match
#' * accept data.frames as input objects
#'
#' @section Function version: 0.1.1
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [calc_CWOSLcomponents], [analyse_SAR.CWOSL.deconvolved], [calc_deconvolution.intervals]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for dose evaluation within the R luminescence package
#'
#' @md
#' @export

calc_superposition.curve <- function(
  object,
  record_type = "OSL",
  aliquot_selection = NULL,
  offset_value = 0,
  output.plot = TRUE
){


  ##============================================================================##
  # Function Setup and Input check
  ##============================================================================##



  #library(Luminescence)


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
  all.values <- NULL
  shortest_record_length <- Inf
  n <- 0
  old_x_axis <- c(0)


  for (j in aliquot_selection) {
    if (j < 1 || j > length(object)) {
      message("Warning: Item ", j," is not a part of the data set. Item skipped")

    } else {
      if (class(object[[j]]) != "RLum.Analysis") {
        message("Warning: Item ", j," is not of class RLum.Analysis. Item skipped")

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

            # prove if the x-axis for consistency
            new_x_axis <- records[[i]]@data[,1][1:shortest_record_length]
            if (length(old_x_axis)>1) {
              if (sum(new_x_axis) != sum(old_x_axis[1:shortest_record_length])) {
                message("Warning (Record ", i, " in list item ", j, ": x-axis values differ from previous records")
              }
            }
            old_x_axis <- new_x_axis

            # count the records
            n <- n + 1

          }
        }
      }
    } ######### LOOP ENDS ########
  }

  # cut curve length to shortest record and calc mean
  #mean.curve <- cbind(new_x_axis, mean.values[1:shortest_record_length] / n)

  mean.curve <- data.frame(time = new_x_axis, signal = mean.values[1:shortest_record_length] / n)

  writeLines(paste0("#1: Arithmetic mean curve with length ", shortest_record_length, " from ", n, " ", record_type, " records"))

  ##============================================================================##
  # PLOT MEAN CURVE
  ##============================================================================##

  if (output.plot == TRUE) {

    library(gridExtra)
    library(ggplot2)
    theme_set(theme_bw())

    colnames(all.values) <- c("time", "signal", "record")
    #curve <- data.frame(time = mean.curve[,1], signal = mean.curve[,2], record = factor("mean"))
    curve <- data.frame(mean.curve, record = factor("mean"))

    alpha.value <- 10 / n
    if (alpha.value > 0.5) alpha.value <- 0.5

    # plot a multiple line summary using the ggplot2 library
    p.lin <- ggplot(all.values, aes(x=time, y=signal, group=record))+
      geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE)+
      geom_line(data = curve, size = 1.5)+
      ylim(0, round(max(curve$signal) * 2))


    #plot log scale
    p.log <- ggplot(all.values, aes(x=time, y=signal, group=record))+
      geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE)+
      scale_y_log10(limits = c(1, max(all.values$signal)))+
      geom_line(data = curve, size = 1.5)+
      ylab("log(signal)")

    # use grid function from gridExtra package to display linear and log plot side by side
    grid.arrange(p.lin, p.log, nrow = 1, top = paste0("Global mean curve and signal scattering of ", n, " ", record_type, " records"))

  }




  return(mean.curve)




'  ##============================================================================##
  # INPUT OBJECTS & INTEGRITY CHECKS
  ##============================================================================##

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

  values <- c(0)
  shortest_record_length <- Inf
  n <- 0
  old_x_axis <- c(0)


  ##============================================================================##
  # MAIN CODE
  ##============================================================================##

  for (j in aliquot_selection) {
    if (j < 1 || j > length(object)) {
      message("Warning: Item ", j," is not a part of the data set. Item skipped")
    } else {
      if (class(object[[j]]) != "RLum.Analysis") {
        message("Warning: Item ", j," is not of class RLum.Analysis. Item skipped")
      } else {######### LOOP BEGINS ########

        records <- object[[j]]@records
        for (i in c(1:length(records))) {

          if (records[[i]]@recordType == record_type) {

            values <- values + records[[i]]@data[,2] - offset_value

            # remember the length of the shortest record to cut the superposition curve later
            if (length(records[[i]]@data[,2]) < shortest_record_length) {
              shortest_record_length <- length(records[[i]]@data[,2])
            }

            # prove if the x-axis for consistency
            new_x_axis <- records[[i]]@data[,1][1:shortest_record_length]
            if (length(old_x_axis)>1) {
              if (sum(new_x_axis) != sum(old_x_axis[1:shortest_record_length])) {
               message("Warning (Record ", i, " in list item ", j, ": x-axis values differ from previous records")
              }
            }
            old_x_axis <- new_x_axis

            # count the records
            n <- n + 1

          }
        }
      }
    } ######### LOOP ENDS ########
  }

  info_text <- "Sum"
  if (mean.values) {
    values <- values/n
    info_text <- "Arithmetic mean"}
  message(info_text, " curve with length ", shortest_record_length, " from ", n, " ", record_type, " records calculated")


  ##============================================================================##
  # RETURN VALUES
  ##============================================================================##

  results <- set_RLum(class  = "RLum.Data.Curve",
           recordType = record_type,
           curveType = info_text,
           data = cbind(new_x_axis, values[1:shortest_record_length]),
           info = list())

  if(output.plot){
    title_text <- paste(info_text, " of ", n, record_type, " records")
    plot_RLum.Data.Curve(results, main = title_text)
  }

  return(results)'
}
