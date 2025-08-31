#' @title Advanced plot function for component resolved CW-OSL curves
#'
#' @description This function is used for plotting CW-OSL curves and its signal components. It can handle data returned
#' by [fit_OSLcurve] or [decompose_OSLcurve]. Besides CW-OSL curves, pseudoLM-OSL curves and
#' residual plots can also be plotted.
#'
#' **Change graph types with parameter:** `display`
#'
#' \tabular{ll}{
#'  `"detailed"` \tab (default) Output plot consists of: Linear CW-OSL plot, pseudoLM-OSL plot, residual curve and component table \cr
#'  `"lin"` \tab Linear CW-OSL plot only \cr
#'  `"compare_lin"` \tab Linear CW-OSL plot with residual curve below and component table on bottom. Useful if two CW-OSL measurements shall be compared side by side \cr
#'  `"log"` \tab CW-OSL plot with logarithmic y-Axis and linear x-Axis \cr
#'  `"compare_log"` \tab CW-OSL plot with logarithmic y-Axis with residual curve below and component table on bottom. Useful if two CW-OSL measurements shall be compared side by side \cr
#'  `"loglog"` \tab Double-logarithmic CW-OSL plot \cr
#'  `"LM"` \tab PseudoLM-OSL plot \cr
#'  `"res"` \tab Plot of residual curve: Measurement minus fitting model \cr
#'  `"tab"` \tab Table of component parameters as image \cr
#'  `"raw"` \tab Raw x-y plot without further data
#' }
#'
#' PseudoLM-OSL curves are created using the transformation described by Bulur (2000).
#' The stimulation ramp duration is twice the CW-OSL duration.
#' See Bos and Wallinga (2012) for a detailed explanation and discussion.
#'
#' This function was black-box tested prior release.
#' These tests as well as code examples are available at:
#' https://luminescence.de/OSLdecomposition/module_tests/Test_plot_OSLcurve.html
#'
#'
#' @param curve [data.frame] or [matrix] or [Luminescence::RLum.Data.Curve-class] (*optional*):
#' CW-OSL curve x-Axis: `$time` or first column as measurement time (must have constant time intervals);
#' y-Axis: `$signal` or second column as luminescence signal.
#' Other columns will be plotted as component curves, in case no input object `components` is defined.
#' If no input is given, a CW-OSL curve will be simulated with the parameters of `components`
#'
#' @param components [data.frame] (*optional*):
#' Table with OSL component parameters. The parameters are used to approximate separate signal decay curves
#' for each component. Need to have at least the columns: `$names`, `$lambda` and `$n`.
#' If an insufficient or no input object is provided, the `curve`` object will be searched for
#' component-related signal values.
#'
#' @param display [character] (*with default*):
#' Sets the arrangement of graphs, see section **Details**.
#'
#' @param show.legend [logical] (*with default*):
#' Draws a legend in the top right corner of the first graph.
#'
#' @param show.intervals [logical] (*with default*):
#' Draws vertical lines into the residual plot showing the signal bin intervals (if available) for the CW-OSL
#' decomposition with [decompose_OSLcurve].
#'
#' @param show.crosssec [logical] (*with default*):
#' Displays photoionisation cross section values in the component table (if available).
#'
#' @param show.initial [logical] (*with default*):
#' Displays signal share at the first channel in the component table (if available).
#'
#' @param title [character] (*with default*):
#' Plot title. Overwrites automatic titles but affects just the first (upper left)
#' diagram in case of multi-plot display setting.
#' Set `title = NULL` for auto-title and `title = ""` for no title.
#'
#' @param graph.colors [character] vector (*optional*):
#' Color for the graphs/stacked areas are defined in the following order: 1. Measurement,
#' 2. Model, 3. Component 1, 4. Component 2, etc. The color vector is allowed to
#' be shorter than the needed colors. For missing colors, the default colors will be used
#'
#' @param theme.set [ggplot2::ggplot2-package] object (*with default*):
#' Graphical theme of the output plot. This argument is forwarded to [ggplot2::theme_set].
#' Recommended themes are `ggplot2::theme_minimal()`, `ggplot2::theme_classic()` and `ggplot2::theme_bw()`,
#' see [ggplot2::theme_bw] or [here](https://ggplot2.tidyverse.org/reference/ggtheme.html) for
#' a full list.
#'
#' @param hide.plot [logical] (*with default*):
#' If true, plot is not drawn but can still be saved as file or caught by `A <- plot_OSLcurve(...)`.
#' If caught, the plot can be drawn manually for example by using [gridExtra::grid.arrange].
#'
#' @param filename [character] (*optional*):
#' File name or path to save the plot as image. If just a file name is given, the image is
#' saved in the working directory. The image type is chosen by the file ending. Both, vector images
#' as well as pixel images are possible. Allowed are `.pdf`, `.eps`, `.svg` (vector graphics), `.jpg`, `.png`, `.bmp` (pixel graphics)
#' and more, see [ggplot2::ggsave].
#'
#' @return
#' An invisible [ggplot2::ggplot] object containing the diagram will be returned. "Invisible" means, the no value
#' will be returned (e.g. no console printout) if the function is not assigned to a variable via `<-`.
#' If the function is assigned, the returned object can be further manipulated by [ggplot2::ggplot2-package] methods
#' or manually drawn by various functions like for example [gridExtra::grid.arrange].
#'
#' @section Last updates:
#'
#' 2025-08-29, DM: Function is now a wrapper of [plot_MultiExponential] which makes it more stable and universal.
#' 2025-08-29, DM: Changed from component graphs to colored areas. Changed also default colors..
#'
#' @author
#' Dirk Mittelstraß, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite the package the following way:
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' R package OSLdecomposition: Automated identification and separation of quartz CW-OSL signal components, *in preparation*.
#'
#' @seealso [fit_OSLcurve], [RLum.OSL_decomposition], [RLum.OSL_global_fitting], [simulate_OSLcomponents]
#'
#' @references
#' Bos, A. J. J. and Wallinga, J., 2012. How to visualize quartz OSL signal components,
#' Radiation Measurements, 47(9)
#'
#' Bulur, E., 2000. A simple transformation for converting CW-OSL curves to LM-OSL curves,
#' Radiation Measurements, 32(2)
#'
#'
#' @examples
#'
#' # Set some arbitrary decay parameter for a dim CW-OSL measurement of quartz
#' components <- data.frame(name = c("fast", "medium", "slow"),
#'                          lambda = c(2, 0.5, 0.02),
#'                          n = c(1000, 1000, 10000))
#'
#' # Simulate a CW-OSL curve including some signal noise
#' curve <- simulate_OSLcomponents(components, simulate.curve = TRUE, add.poisson.noise = TRUE)
#'
#' # Display the simulated curve
#' plot_OSLcurve(curve, components)
#'
#' @md
#' @export
plot_OSLcurve <- function(curve = NULL,
                          components = NULL,
                          display = "detailed",
                          show.legend = FALSE,
                          show.intervals = FALSE,
                          show.crosssec = FALSE,
                          show.initial = FALSE,
                          title = NULL,
                          graph.colors = NULL,
                          theme.set = ggplot2::theme_classic(),
                          hide.plot = FALSE,
                          filename = NULL){

  # Changelog:
  # * 2019-03-06, DM: First reasonable version
  # * 2019-04-03, DM: Rebuild whole function
  # * 2019-10-02, DM: Added background component support and some little tweaks
  # * 2020-04-22, DM: Enabled hidden output to draw later
  # * 2020-06-19, DM: Added pseudoLM-OSL diagrams
  # * 2020-08-04, DM: Added subtitles and RSS info
  # * 2020-09-02, DM: Added graphic saving with [ggplot2::ggsave]
  # * 2020-11-03, DM: Refactored code; Changed parameter `display` choices and added parameters `theme.set` and `show.legend`
  # * 2020-11-11, DM: Improved parameter table
  # * 2020-11-16, DM: Added parameter `show.crosssec` and `show.initial`
  # * 2020-11-19, DM: RLum.Data.Curves can now be displayed without 'components' if RLum.OSL_decomposition was already performed at data set
  # * 2021-03-29, DM: Hidden output objects are now [ggplot2] objects if the plot is not a composite diagram
  # * 2025-08-29, DM: Function is now a wrapper of [plot_MultiExponential] which makes it more stable and universal.
  # * 2025-08-29, DM: Changed from component graphs to colored areas. Changed also default colors..
  #
  # ToDo:
  # * When drawing components without curve, skip residual curve
  # * Display fitting formula
  # * Add argument ggsave.control() to give direct control about image saving
  # * Add lambda error in the component table
  # * Draw bin intervals directly into the diagram (optionally)

  ###################### INPUT DATA  ######################

  if (is.null(curve) && is.null(components)) stop("[plot_OSLcurve()]: Either 'curve' or 'components' must be defined")

  # Create curve from component table if not given
  if (is.null(curve) && !is.null(components)) {

    channel.width <- 1 / (2*max(components$lambda))
    channel.number <- ceiling(1 / (min(components$lambda) * channel.width))
    curve <- simulate_OSLcomponents(components,
                               channel.width = channel.width,
                               channel.number = channel.number,
                               simulate.curve = TRUE,
                               add.poisson.noise = FALSE)
  } else {

    # Check input curve
    if(!inherits(curve, c("RLum.Data.Curve", "data.frame", "matrix"))){
      stop("[plot_OSLcurve()] Error: Input object 'curve' is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")}

    if(inherits(curve, "RLum.Data.Curve")) {

      # if no component table is given, check if the data set was already decomposed by RLum.OSL_decomposition
      if (is.null(components)) {

        if ("COMPONENTS" %in% names(curve@info)) {

          components <- curve@info$COMPONENTS

        }else{
          # switch to simple plot
          display <- "raw"}}

      # now transform the RLum object into a simple table. There is no need to keep all the extra info
      curve <- as.data.frame(Luminescence::get_RLum(curve))
    }

    if (!("time" %in% colnames(curve)) |
        !("signal" %in% colnames(curve))) {

      curve <- data.frame(time = curve[,1], signal = curve[,2])}
  }

  time <- curve$time

  # if a component table is given, add or overwrite the component signal columns in the curve data.frame
  if (!is.null(components))
    curve <- simulate_OSLcomponents(components, curve, simulate.curve = FALSE)

  # Does 'curve' now contain component information?
  if (length(curve[1,]) < 3) {

    warning("No component curve given neither in 'components' nor in 'curve'. Set 'display' = 'raw'")
    display <- "raw"
  }

  ###################### PLOT DESIGN ######################

  # Set color and line themes
  ggplot2::theme_set(theme.set)

  guide <- "none"
  if (show.legend) guide <- "right"

  # Old colors of OSLdecomposition 1.0.0
  # graph.colors <- c("darkgrey", "black","red3","green3","royalblue3","darkorchid","gold","brown","pink")

  # New colors
  graph_colors <- c("grey30", "darkblue", "skyblue2","orchid","cyan2","orange","red2","pink3","brown2")
  if (length(graph.colors) > 0)
    graph_colors[1:length(graph.colors)] <- graph.colors

  # The colors are NOT assigned to the components iterativly
  # see plot_MultiExponential for the rules. However, we have to resort
  # the colors to be in line with plot_MultiExponential
  if(display != "raw" && !is.null(components) && nrow(components) > 0){

    component_colors <- rep("black", nrow(components))
    col_i <- 3
    for(i in nrow(components):1){

      if (components$n[i] > 0) {
        component_colors[i] <- graph_colors[col_i]
        col_i <- col_i + 1
      }
    }
    for(i in nrow(components):1){

      if (components$n[i] < 0) {
        component_colors[i] <- graph_colors[col_i]
        col_i <- col_i + 1
      }
    }
  }



  ######################## LINEAR PLOT #########################################################

  if ((display == "detailed") | (display == "lin") | (display == "compare_lin")) {

    if (is.null(title)) title <- "CW-OSL"

    p.lin <- plot_MultiExponential(curve, components,
                                   main = title,
                                   xlab = "Time (s)", ylab = "Signal (cts)",
                                   theme.set = theme.set,
                                   font.size = 8,
                                   legend.position = guide,
                                   hide.plot = TRUE)

    # Deactivate legend for all following plots
    guide <- "none"
    title <- NULL

  }

  ######################## LOG PLOT #########################################################

  if ((display == "log") | (display == "compare_log")) {

    if (is.null(title)) title <- "CW-OSL (logarithmic scale)"

    p.log  <- plot_MultiExponential(curve, components,
                                    main = title,
                                    theme.set = theme.set,
                                    ylog = TRUE,
                                    font.size = 8,
                                    legend.position = guide,
                                    xlab = "Time (s)", ylab = "Signal (cts)",
                                    hide.plot = TRUE)
    guide <- "none"
    title <- NULL
  }

  ######################## LOG-LOG PLOT #########################################################

  if (display == "loglog") {

    if (is.null(title)) title <- "CW-OSL (double logarithmic scale)"

    p.loglog  <- plot_MultiExponential(curve, components,
                                       main = title,
                                       theme.set = theme.set,
                                       ylog = TRUE, xlog = TRUE,
                                       font.size = 8,
                                       legend.position = guide,
                                       xlab = "Time (s)", ylab = "Signal (cts)",
                                       hide.plot = TRUE)
    guide <- "none"
    title <- NULL
  }


  ######################## pseudoLM PLOT #########################################################

  if ((display == "detailed") | (display == "LM")) {

    if (is.null(title)) title <- "pseudoLM-OSL"

    p.LM  <- plot_MultiExponential(curve, components,
                                   main = title,
                                   theme.set = theme.set,
                                   linear.modulation = TRUE,
                                   xlab = "Time (s)", ylab = "Signal (cts)",
                                   font.size = 8,
                                   legend.position = guide,
                                   hide.plot = TRUE)
    guide <- "none"
    title <- NULL
  }


  ######################## RESIDUAL #########################################################

  if ((display == "detailed") | (display == "compare_lin") |
      (display == "compare_log") | (display == "res")) {

    res <- curve$residual
    res_curve <- data.frame(time = curve$time, signal = res)

    if (is.null(title)) title <- "Residual"
    # IMPORTANT: Deactivated RSS display for the moment as there is a mismatch
    #            between RSS here and those returned by fit_OSLcurve()
    # if (is.null(title))
    #   title <- paste0("Residual (RSS = ", prettyNum(sum(res^2, na.rm = TRUE), digits = 1), ")")

    # Plot residual curve
    p.res  <- plot_MultiExponential(res_curve,
                                    main = title,
                                    theme.set = theme.set,
                                    xlab = "Time (s)", ylab = "Signal (cts)",
                                    font.size = 8,
                                    legend.position = guide,
                                    hide.plot = TRUE)

    # if the integration intervals are given, draw them into the residual curve
    if (("t.start" %in% colnames(components))
        & show.intervals) {

      res.intervals <- list(NULL)

      # set y-axis
      res.max <- ceiling(abs(max(res))) + 1
      if (abs(min(res)) > res.max) res.max <- ceiling(abs(min(res))) + 1

      for (i in c(1:nrow(components))) {

        res.intervals[[i*2 - 1]] <- ggplot2::annotate("segment",
                                                      x = components$t.start[i],
                                                      xend = components$t.start[i],
                                                      y = - res.max, yend = res.max, colour = "black")

        res.intervals[[i*2]] <- ggplot2::annotate("segment",
                                                  x = components$t.end[i],
                                                  xend = components$t.end[i],
                                                  y = - res.max, yend = res.max, colour = "black")}

      t.times <- c(components$t.start[1], components$t.end)
      t.labels <- formatC(t.times)
      if (t.times[1] == 0) t.labels <- c("", formatC(components$t.end))
      if (t.labels[length(t.labels)] == time[length(time)]) t.labels[length(t.labels)] <- ""

      p.res <- p.res + res.intervals +
        ggplot2::scale_x_continuous(breaks = t.times, limits = c(0, max(time)), labels = t.labels)
    }
  }

  ######################## TABLE #########################################################

  if ((display == "detailed") | (display == "compare_lin") |
      (display == "compare_log") | (display == "tab")) {

    if (!is.null(components) ) {

      # Table styles can be found at:
      # Link: rpkgs.datanovia.com/ggpubr/reference/ggtexttable.html

      ### Build Names and Lambda column
      p.colnames <- c("      ",
                      expression(italic(lambda) ~~ (s^-1)))

      print.lambda <- prettyNum(components$lambda, digits = 3)
      print.lambda[print.lambda == "NA"] <- ""

      p.table <- data.frame(name = components$name,
                            lambda = print.lambda)

      ### Are photo cross-sections given?
      if (show.crosssec & ("cross.section" %in% colnames(components))){
        p.table <- data.frame(p.table,
                              cross = prettyNum(components$cross.section, digits = 3))
        p.colnames <- c(p.colnames, expression(italic(sigma) ~~ (cm^2)))}

      ### Add Intensity column
      p.table <- data.frame(p.table,
                            n = prettyNum(round(components$n)))
      p.colnames <- c(p.colnames, expression(italic(n)))

      if ("n.error" %in% colnames(components)) {
        p.table$n <- paste0(p.table$n, " \u00b1 ", prettyNum(round(components$n.error)))}

      ### Are shares of initial signal given?
      if (show.initial & ("initial.signal" %in% colnames(components))){
        p.table <- data.frame(p.table,
                              initial = scales::percent(components$initial.signal,
                                                        suffix = "%",
                                                        accuracy = 1))
        p.colnames <- c(p.colnames, expression(italic(I)[0]))}


      colnames(p.table) <- p.colnames
      print.size = 8

      # draw component table
      p.tab <-  ggpubr::ggtexttable(p.table,
                                    rows = NULL,
                                    theme = ggpubr::ttheme(base_style = "default",
                                                           colnames.style = ggpubr::colnames_style(fill = "white",
                                                                                                   face = "bold",
                                                                                                   size = print.size,
                                                                                                   parse = TRUE),
                                                           tbody.style = ggpubr::tbody_style(fill = c("white","white"),
                                                                                             size = print.size)))

      # assign graph colors to component names
      for (z in 1:nrow(p.table)) {
        p.tab <- ggpubr::table_cell_bg(p.tab, row = 1 + z, column = 1, fill = component_colors[z],
                                       color = "white", linewidth = 5)}

    } else {

      # if no component table is given, draw blank image
      p.table <- data.frame(NULL)
      p.tab <- ggplot2::ggplot() + ggplot2::geom_blank()}}


  ######################## RAW PLOT ########################################################

  if (display == "raw") {

    # Create plot without components
    p.raw <- plot_MultiExponential(curve,
                                   main = title,
                                   theme.set = theme.set,
                                   xlab = "Time (s)", ylab = "Signal (cts)",
                                   font.size = 8,
                                   legend.position = guide,
                                   hide.plot = TRUE)
  }


  #########################################################################################
  ######                              Arrange plots                                  ######
  #########################################################################################


  if (display == "detailed") {

    # Create layout matrix, see ?gridExtra::arrangeGrob for details
    # As the component table needs more space with increasing number of components,
    # pane 4 is increased in height
    tab.shift <- 0
    if (nrow(p.table) > 4) tab.shift <- nrow(p.table) - 4

    lay <- matrix(c(rep(1, 10),
                    rep(2, 5),
                    rep(3, 10 - tab.shift),
                    rep(4, 5 + tab.shift)),
                  ncol = 2)

    plot_object <- gridExtra::arrangeGrob(p.lin, p.res, p.LM, p.tab, layout_matrix = lay)


  } else if (display == "compare_lin") {

    lay <- rbind(c(1), c(1), c(2), c(3))

    plot_object <- gridExtra::arrangeGrob(p.lin, p.res, p.tab, layout_matrix = lay)

  } else if (display == "compare_log") {

    lay <- rbind(c(1), c(1), c(2), c(3))

    plot_object <- gridExtra::arrangeGrob(p.log, p.res, p.tab, layout_matrix = lay)

  } else if (display == "lin") {

    plot_object <- p.lin

  } else if (display == "log") {

    plot_object <- p.log

  } else if (display == "loglog") {

    plot_object <- p.loglog

  } else if (display == "LM") {

    plot_object <- p.LM

  } else if (display == "res") {

    plot_object <- p.res

  } else if (display == "tab") {

    plot_object <- p.tab

  } else if (display == "raw") {

    plot_object <- p.raw

  } else {

    stop("No valid argument 'display'")}



  # save plot as file
  if (!is.null(filename)) {

    try(suppressMessages(ggplot2::ggsave(filename, plot = plot_object, units = "cm")),
        silent = FALSE)}

  # show plot
  if (!hide.plot) gridExtra::grid.arrange(plot_object)

  # return plot object
  invisible(plot_object)
}
