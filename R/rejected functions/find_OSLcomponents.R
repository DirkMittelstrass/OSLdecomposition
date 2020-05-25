find_OSLcomponents <- function(
  object,
  aliquot_selection = NA,
  n.max = 3,
  record_type = "OSL",
  offset_value = 0,
  add_to_data = FALSE,
  output_plot = TRUE
){

  #object <- Risoe.BINfileData2RLum.Analysis(read_BIN2R(choose.files(multi = FALSE)))

  ##============================================================================##
  # INTEGRITY CHECKS
  ##============================================================================##

  library(ggplot2)
  theme_set(theme_bw())

  library(Luminescence)
  library(gridExtra)

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
  mean.curve <- cbind(new_x_axis, mean.values[1:shortest_record_length] / n)
  writeLines(paste0("#1: Arithmetic mean curve with length ", shortest_record_length, " from ", n, " ", record_type, " records"))

  ##============================================================================##
  # PLOT MEAN CURVE
  ##============================================================================##


  colnames(all.values) <- c("time", "signal", "record")
  mean.curve <- data.frame(time = mean.curve[,1], signal = mean.curve[,2], record = factor("mean"))
  alpha.value <- 10 / n
  if (alpha.value > 0.5) alpha.value <- 0.5

  # plot a multiple line summary using the ggplot2 library
  p.lin <- ggplot(all.values, aes(x=time, y=signal, group=record))+
    geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE)+
    geom_line(data = mean.curve, size = 1.5)+
    ylim(0, round(max(mean.curve$signal) * 2))


  #plot log scale
  p.log <- ggplot(all.values, aes(x=time, y=signal, group=record))+
    geom_point(size = 0.5, alpha = alpha.value, na.rm = TRUE)+
    scale_y_log10()+
    geom_line(data = mean.curve, size = 1.5)+
    ylab("log(signal)")

  # use grid function from gridExtra package to display linear and log plot side by side
  grid.arrange(p.lin, p.log, nrow = 1, top = paste0("Global mean curve and signal scattering of ", n, " ", record_type, " records"))


    # transform mean curve into RLum object
'  mean.curve <- set_RLum(class  = "RLum.Data.Curve",
                      recordType = record_type,
                      curveType = paste0("mean.", record_type),
                      data = mean.curve,
                      info = list())'
  #return(mean.curve)

  ##============================================================================##
  # FIT COMPONENTS
  ##============================================================================##

  # find components via fitting
  fit.results <- Luminescence::fit_CWCurve(mean.curve,
                             n.components.max = n.max,
                             fit.calcError = TRUE,
                             plot =  TRUE)

  return(fit.results)


  f_fast <- fit_results@data[["data"]][["lambda1"]]
  f_medium <- fit_results@data[["data"]][["lambda2"]]
  f_slow <- fit_results@data[["data"]][["lambda3"]]

  # if just one component could be found
  if (is.na(f_medium) && is.na(f_slow)) message("Just one component found, defined as 'fast' component by default")

  # if no third component could be found, set second component as 'slow'
  if (!is.na(f_medium) && is.na(f_slow)) {
    f_slow <- f_medium
    f_medium <- NA
    message("Just two signal components found, defined as 'fast' and 'slow' component by default")
  }


  ##============================================================================##
  # CALC INTERVALS
  ##============================================================================##

  # get the deconvolution intervals
  intervals <- calc_deconvolution.intervals(f.fast = f_fast,
                                            f.medium = f_medium,
                                            f.slow = f_slow,
                                            channel.width = superposition_curve[2,1] - superposition_curve[1,1],
                                            channel.number = length(superposition_curve[,1]))

  results$decay.parameter <- rbind(results$decay.parameter,
                                   data.frame(f.fast = f_fast,
                                              f.medium = f_medium,
                                              f.slow = f_slow,
                                              intervals))

  ##============================================================================##
  # ADD TO LUM DATA
  ##============================================================================##


}
