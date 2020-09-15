#' Advanced plot function for CW-OSL curves and CW-OSL signal components
#'
#' @param curve
#' @param components
#' @param display
#' Choices: "detailed" (default), "components", "raw" ,"compare_lin", "compare_log", "lin-log"
#'
#' @param show.intervals
#' Draws vertical lines into the residual plot showing the signal bin intervals for Step 2 CW-OSL decomposition
#'
#' @param title
#' Plot title
#'
#' @param hide.plot
#' If true, plot is not drawn but can be catched by 'A <- plot_OSLcurve(...)' or saved as file.
#' If catched, the plot ca be displayed later for example by [gridExtra::grid.arrange]
#'
#' @param filename
#' Filename or path to save the diagram as file. If just a name is chosen, the file is
#' saved in the working directory. The graphic type is chosen according to the file ending.
#' Allowed are .pdf, .eps, .svg (vector graphics) and .jpg, .png, .bmp (pixel graphics)
#' and more, see [ggplot2::ggsave]
#'
#' @section Notes:
#'
#' pseudoLM-OSL curves are created using the transformation described by Bulur (2000).
#' The stimulation ramp duration is chosen double the CW-OSL duration.
#' See also Bos and Wallinga (2012) for a detailed explanation and discussion
#'
#' @references
#' Bos, A. J. J. and Wallinga, J.: How to visualize quartz OSL signal components, Radiation Measurements, 47(9), 752–758, doi:10.1016/j.radmeas.2012.01.013, 2012.
#'
#' Bulur, E.: A simple transformation for converting CW-OSL curves to LM-OSL curves, Radiation Measurements, 32(2), 141–145, doi:10.1016/S1350-4487(99)00247-4, 2000.
#'
#'
#' @section Last changed. 2020-09-13
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' @references
#' Bos, A. J. J. and Wallinga, J.: How to visualize quartz OSL signal components,
#' Radiation Measurements, 47(9), 752–758, doi:10.1016/j.radmeas.2012.01.013, 2012.
#'
#' Bulur, E.: A simple transformation for converting CW-OSL curves to LM-OSL curves,
#' Radiation Measurements, 32(2), 141–145, doi:10.1016/S1350-4487(99)00247-4, 2000.
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' Automated identification and separation of quartz CW-OSL signal components with R, *in preparation*.
#'
#' @export
#'
#' @examples
#'
#' # Set some reasonable parameter for a weak quartz CW-OSL decay
#' components <- data.frame(name = c("fast", "medium", "slow"), lambda = c(2, 0.5, 0.02), n = c(1000, 1000, 10000))
#'
#' # Simulate the CW-OSL curve and add some signal noise
#' curve <- simulate_OSLcurve(components, simulate.curve = TRUE, add.poisson.noise = TRUE)
#'
#' # Display the simulated curve
#' plot_OSLcurve(curve, components)
#'
plot_OSLcurve <- function(curve = NULL,
                          components,
                          display = "detailed",
                          show.intervals = FALSE,
                          title = NULL,
                          hide.plot = FALSE,
                          filename = NULL){

#' Changelog:
#' * 2019-03-06, DM: First reasonable version
#' * 2019-04-03, DM: Rebuild whole function
#' * 2019-10-02, DM: Added background component support and some little tweaks
#' * 2020-04-22, DM: Enabled hidden output to draw later
#' * 2020-06-19, DM: Added pseudoLM-OSL diagrams
#' * 2020-08-04, DM: Added subtitles and RSS info
#' * 2020-09-02, DM: Added graphic saving with [ggplot2::ggsave]
#'
#' ToDo:
#' * REFACTORIZE CODE
#' * ! Put the ggplot building (or at least its style options) in its own sub-function !
#' * Add input checks and data conversions from [decompose_OSLcurve]
#' * When drawing components without curve, skip residual curve
#' * Display fitting formula
#' * Get rid of library 'ggpubr' to decrease dependencies
#' * Change from library("XXX") to XXX::
#' * Add argument ggsave.control() to give direct control about image saving
#' * Add table column for the lambda error


  # Hidden parameters
  zoom <- 1

  #library(gridExtra)
  library(ggplot2)
  library(ggpubr)
  library(scales)

  ########## Checks and data preperations ###########

  # Set color and line themes
  theme_set(theme_minimal())
  graph.colors <- c("darkgrey", "black","red3","green3","royalblue3","darkorchid","gold","brown","pink")
  graph.sizes <- c(1, 1, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  graph.shapes <-  c(16, 32, 32, 32, 32, 32, 32, 32, 32)
  graph.lines <- c("blank", "solid","solid","solid","solid","solid","solid","solid","solid")

  text_format <- theme(axis.title = element_text(size = 8),
                       plot.subtitle = element_text(size = 9, face = "bold"))

  # Create curve from component table if not given
  if (is.null(curve) && !is.null(components)) {

    channel.width <- 1 / (2*max(components$lambda))
    channel.number <- ceiling(1 / (min(components$lambda) * channel.width))
    curve <- simulate_OSLcurve(components,
                               channel.width = channel.width,
                               channel.number = channel.number,
                               simulate.curve = TRUE,
                               add.poisson.noise = FALSE)
  } else {

    # Check input curve
    if(is(curve, "RLum.Data.Curve") == FALSE & is(curve, "data.frame") == FALSE & is(curve, "matrix") == FALSE){
      stop("[plot_OSLcurve()] Error: Input object is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")
    }

    if (!("time" %in% colnames(curve)) ||
        !("signal" %in% colnames(curve))) {
      curve <- data.frame(time = curve[,1],
                          signal = curve[,2])
    }
  }

  time <- curve$time
  channel.width <- time[2] - time[1]

  # Set x-axis zoom
  X_limits <- c(0, max(time)*zoom)

  # Transform data into factor-based long data set
  data <- data.frame(time = time, signal = curve$signal, graph = factor("meas"))



  if ((display == "components") || (display == "detailed") ||
      (display == "compare_lin") || (display == "compare_log") ||
      (display == "return.lin") || (display =="lin-log")) {

    # if a component table is given, add or overwrite the component signal columns in the curve data.frame
    if (!is.null(components)) {
      #|| !("residual" %in% colnames(curve))
      curve <- simulate_OSLcurve(components, curve, simulate.curve = FALSE)
    }

    if (length(curve[1,]) < 3) {
      warning("No components given neither in 'components' nor in 'curve'")
      return(NULL)
    }

    # set the columns which shall be displayed
    X <- c(1:length(curve[1,]))
    X <- X[colnames(curve) != "residual"]
    X <- X[3:length(X)]

    # set legend text
    graph.labels <- c("", colnames(curve)[X])
    graph.labels[1] <- "measured"
    graph.labels[2] <- "fitted"



    # rearrange data to work in ggplot-function
    for (i in X) {

      data.append <- data.frame(time = time, signal = curve[,i], graph = colnames(curve)[i])
      data <- rbind(data, data.append)
    }

    ######################## LINEAR PLOT #########################################################

    p.lin <- ggplot(data, aes(x=time, y=signal, colour=graph, size=graph, shape=graph, linetype=graph)) +
      geom_point(na.rm = TRUE) + geom_line(na.rm = TRUE) +
      scale_y_continuous(limits = c(0, max(curve$signal) + 1)) +
      scale_x_continuous(limits = X_limits) +
      scale_colour_manual(values = graph.colors, labels = graph.labels, guide = FALSE) +
      scale_size_manual(values = graph.sizes, labels = graph.labels, guide = FALSE) +
      scale_shape_manual(values = graph.shapes, labels = graph.labels, guide = FALSE) +
      scale_linetype_manual(values = graph.lines, labels = graph.labels, guide = FALSE) +
      labs(subtitle = "CW-OSL", x = "Time (s)", y ="Signal (cts)") +
      text_format



    ######################## pseudoLM PLOT #########################################################

    # transform data
    P = 2*max(time)
    LMdata <- data
    LMdata$signal <- LMdata$signal * (2 * LMdata$time / P)^0.5
    LMdata$time <- (2 * P * LMdata$time)^0.5


    p.LM <- ggplot(LMdata, aes(x=time, y=signal, colour=graph, size=graph, shape=graph, linetype=graph)) +
      geom_point(na.rm = TRUE) + geom_line(na.rm = TRUE) +
      scale_y_continuous(limits = c(0, max(LMdata$signal) + 1)) +
      scale_x_continuous(limits = c(0, X_limits[2]*2)) +
      scale_colour_manual(values = graph.colors, labels = graph.labels, guide = FALSE) +
      scale_size_manual(values = graph.sizes, labels = graph.labels, guide = FALSE) +
      scale_shape_manual(values = graph.shapes, labels = graph.labels, guide = FALSE) +
      scale_linetype_manual(values = graph.lines, labels = graph.labels, guide = FALSE) +
      labs(subtitle = "pseudoLM-OSL", x = "Ramping time (s)", y ="Signal (cts)") +
      text_format



    ######################## RESIDUAL PLOT #########################################################

    res <- curve$residual
    res_text <- paste0("Residual (RSS = ", formatC(sum(res^2, na.rm = TRUE)), ")")

    # set y-axis
    res.max <- ceiling(abs(max(res))) + 1
    if (abs(min(res)) > res.max) res.max <- ceiling(abs(min(res))) + 1

    # if the intervals are given, draw them into the residual curve
    res.intervals <- list(NULL)
    if (("t.start" %in% colnames(components))
        & (zoom == 1)
        & show.intervals) {

      for (i in c(1:nrow(components))) {

        res.intervals[[i*2 - 1]] <- annotate("segment",
                                             x = components$t.start[i],
                                             xend = components$t.start[i],
                                             y = - res.max, yend = res.max, colour = "black")
        res.intervals[[i*2]] <- annotate("segment",
                                         x = components$t.end[i],
                                         xend = components$t.end[i],
                                         y = - res.max, yend = res.max, colour = "black")
      }
      t.times <- c(components$t.start[1], components$t.end)
      t.labels <- formatC(t.times)
      if (t.times[1] == 0) t.labels <- c("", formatC(components$t.end))
      if (t.labels[length(t.labels)] == time[length(time)]) {
        t.labels[length(t.labels)] <- ""
      }

      scale.intervals <- scale_x_continuous(breaks = t.times, limits = X_limits,
                                            labels = t.labels)#, position = "top")

    } else {

      scale.intervals <- scale_x_continuous(limits = X_limits)#, position = "top")
    }


    ########## Plot residual curve ###########

    p.res <- ggplot(curve, aes(x=time, y=residual)) +
      #error.ribbon +
      geom_point(size = 1, shape =  16, color = "darkgrey", na.rm = TRUE) +
      scale_y_continuous(limits = c(- res.max, res.max)) +
      labs(subtitle = res_text, x = "Time (s)", y ="Signal (cts)") +
      annotate("segment", x = 0, xend = max(curve$time), y = 0, yend = 0, colour = "black", size = 1) +
      text_format +
      scale.intervals +
      res.intervals




    ######################## LOG LOG PLOT #########################################################
    if (display == "detailed" || (display == "compare_lin") ||
        (display == "compare_log") || (display == "lin-log")) {

      # set y-axis minimum
      log_limits <- c(1, max(data$signal))
      if (min(curve$signal) > 0) {
        log_limits[1] <- 10^floor(log10(min(curve$signal)))
      }

      breaks <- 10^(-10:10)
      minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))


      # draw plot
      p.log <- ggplot(data, aes(x=time, y=signal, colour=graph, size=graph, shape=graph, linetype=graph)) +
        geom_point(na.rm = TRUE) + geom_line(na.rm = TRUE) +
        scale_y_log10(limits = log_limits, position = "right",
                      breaks = breaks,
                      minor_breaks = minor_breaks) +
        scale_x_log10(limits = c(channel.width, max(time)), position = "right",
                      breaks = breaks,
                      minor_breaks = minor_breaks) +
        scale_colour_manual(values = graph.colors, labels = graph.labels, guide = FALSE) +
        scale_size_manual(values = graph.sizes, labels = graph.labels, guide = FALSE) +
        scale_shape_manual(values = graph.shapes, labels = graph.labels, guide = FALSE) +
        scale_linetype_manual(values = graph.lines, labels = graph.labels, guide = FALSE) +
        ylab("Signal (cts)") + xlab("Time (s)")  +
        text_format


      ######################## LIN LOG PLOT #########################################################
      if (display == "lin-log") {

        # set y-axis minimum
        log_limits <- c(1, max(data$signal))
        if (min(curve$signal) > 0) {
          log_limits[1] <- 10^floor(log10(min(curve$signal)))
        }

        breaks <- 10^(-10:10)
        minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))


        # draw plot
        p.linlog <- ggplot(data, aes(x=time, y=signal, colour=graph, size=graph, shape=graph, linetype=graph)) +
          geom_point(na.rm = TRUE) + geom_line(na.rm = TRUE) +
          scale_y_log10(limits = log_limits, position = "right",
                        breaks = breaks, minor_breaks = minor_breaks,
                        labels = scientific) +
          scale_x_continuous(limits = X_limits, labels = scientific) +
          scale_colour_manual(values = graph.colors, labels = graph.labels, guide = FALSE) +
          scale_size_manual(values = graph.sizes, labels = graph.labels, guide = FALSE) +
          scale_shape_manual(values = graph.shapes, labels = graph.labels, guide = FALSE) +
          scale_linetype_manual(values = graph.lines, labels = graph.labels, guide = FALSE) +
          ylab("Signal (cts)") + xlab("Time (s)")  +
          text_format

      }
      ######################## TABLE #########################################################

      if (!is.null(components) ) {

        # Table styles can be found at:
      # https://rpkgs.datanovia.com/ggpubr/reference/ggtexttable.html

      p.colnames <- c("      ",
                      expression(italic(lambda) ~~ (s^-1)),
                      expression(italic(n)))

      #print.lambda <- components$lambda
      print.lambda <- formatC(components$lambda, digits = 3)
      print.lambda[print.lambda == "NA"] <- ""

      p.table <- data.frame(name = components$name,
                            lambda = print.lambda,
                            n = prettyNum(round(components$n)))

      # add "prettyNum()" to shorten big n numbers
      if ("n.error" %in% colnames(components)) {
        p.table <- data.frame(p.table,
                              n.error = prettyNum(round(components$n.error)))
        p.colnames <- c(p.colnames,
                        expression(sigma[italic(n)]))
      }

      #if ("n.residual" %in% colnames(components)) {
      #  p.table <- data.frame(p.table,
      #                        n.residual = prettyNum(round(components$n.residual)))
      #  p.colnames <- c(p.colnames,
      #                  expression(tail[italic(n)]))
      #  p.table$n.residual[is.na(p.table$n.residual)] <- ""
      #}


      colnames(p.table) <- p.colnames

      print.size = 8

      p.tab <-  ggtexttable(p.table,
                            rows = NULL,
                            theme = ttheme(base_style = "default",
                                           colnames.style = colnames_style(fill = "white",
                                                                           face = "bold",
                                                                           size = print.size,
                                                                           parse = TRUE),
                                           tbody.style = tbody_style(fill = c("white","white"),
                                                                     size = print.size)))

      for (z in 1:nrow(p.table)) {
        p.tab <- table_cell_bg(p.tab, row = 1 + z, column = 1, fill = graph.colors[z+2],
                               color = "white", linewidth = 5)
      }


      } else {

        p.tab <- ggplot() + geom_blank()
      }
    }
    ##################################################################################################
  } else {

    ########## Create plot without components ###########

    # plot using the ggplot2 library
    p.lin <- ggplot(curve, aes(x=time, y=signal)) +
      geom_line(na.rm = TRUE, color = "black") +
      scale_y_continuous(limits = c(0, max(curve$signal))) +
      scale_x_continuous(limits = X_limits) +
      # geom_point(size = 1, shape =  16, color = "darkgrey", na.rm = TRUE) +
      #  ylim(0, round(max(curve$signal) * 1.1)) +
      ylab("signal (cts)") + xlab("time (s)")

  }


  ########## Arrange and display plots ###########


  if (display == "return.lin") {

    plot_object <- p.lin

  } else if (display == "components") {

    # display layout matrix:
    lay <- cbind(c(1,1,2))

    # use grid function from gridExtra package to display linear and log plot side by side
    plot_object <- gridExtra::arrangeGrob(p.lin, p.res, layout_matrix = lay, top = title)

  } else if (display == "detailed") {

    'lay <- rbind(c(1,1,1,3,3),
                 c(1,1,1,3,3),
                 c(1,1,1,3,3),
                 c(1,1,1,4,4),
                 c(2,2,2,4,4),
                 c(2,2,2,4,4))'

    lay <- rbind(c(1,3),
                 c(1,4),
                 c(2,4))

    #nrow(p.table)
    tab.shift <- 0
    if (nrow(p.table) > 4) tab.shift <- nrow(p.table) - 4

    lay <- matrix(c(rep(1,10),
                    rep(2,5),
                    rep(3,10 - tab.shift),
                    rep(4,5 + tab.shift))
                  ,ncol = 2)

    #grid.arrange(p.lin, p.res, p.log, p.tab, layout_matrix = lay, top = title)
    plot_object <- gridExtra::arrangeGrob(p.lin, p.res, p.LM, p.tab, layout_matrix = lay, top = title)

  } else if (display == "compare_lin") {

    lay <- rbind(c(1),
                 c(1),
                 c(2),
                 c(3))

    plot_object <- gridExtra::arrangeGrob(p.lin, p.res, p.tab, layout_matrix = lay, top = title)

  } else if (display == "compare_log") {

    lay <- rbind(c(1),
                 c(1),
                 c(2),
                 c(3))

    plot_object <- gridExtra::arrangeGrob(p.log, p.res, p.tab, layout_matrix = lay, top = title)

  } else if (display == "lin-log") {


    lay <- rbind(c(1,1,1,3),
                 c(1,1,1,3),
                 c(2,2,2,3))

    plot_object <- gridExtra::arrangeGrob(p.linlog, p.res, p.tab, layout_matrix = lay, top = title)


  } else {
    plot_object <- gridExtra::arrangeGrob(p.lin, nrow = 1, top = title)}

  # save plot as file
  if (!is.null(filename)) {
    try(ggplot2::ggsave(filename, plot = plot_object, units = "cm"), silent = FALSE)}

  # show plot
  if (!hide.plot) gridExtra::grid.arrange(plot_object)

  # return plot object
  invisible(plot_object)
}
