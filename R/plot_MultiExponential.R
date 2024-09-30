#' Advanced plot function for displaying component resolved signal curves
#'
#' @description This function plots multi-exponentially decaying measurements and its signal components.
#' [plot_OSLcurve] is a wrapper for this function.
#'
#' This function was black-box tested prior release.
#' These tests as well as many code examples are available at:
#' https://luminescence.de/OSLdecomposition/module_tests/Test_decompose_OSLcurve.html
#'
#'
#' @param curve [data.frame] or [matrix] or [Luminescence::RLum.Data.Curve-class] (*optional*):
#' Measured signal curve. First column is used as x-axis, second column is used as y-axis.
#' Further columns are ignored. If this argument is `NULL` but a component table is given,
#' signal components and a modeled multi-exponential signal curve will be displayed nonetheless.
#'
#' @param components [data.frame] or [numeric] vector (**optional**)
#' Either table with the signal component parameters or numeric vector with decay rates.
#' The component parameter table is usually given by [fit_OSLcurve] or [decompose_OSLcurve].
#' If handmade, it needs the columns `$names`, `$lambda` and `$n`.T
#' This argument also accepts numeric vectors, the component intensity will be calculated automatically
#' using [decompose_OSLcurve]. If the vector elements are named, these names will be used as component names.
#'
#' @param fill.components [logical] (*with default*):
#' If `TRUE`, signal components are displayed ad stacked areas. Otherwise they are displayed as line graphs.
#'
#' @param linear.modulation [logical] (*with default*):
#' If `TRUE`, all graphs are transformed to linear modulated curves, a peak-like representation for decay curves.
#' See Bulur (2000) and Bos and Wallinga (2012) for details.
#'
#' @param xlog [logical] (*with default*):
#' If `TRUE`, x-axis is transformed to logarithmic scale.
#'
#' @param ylog [logical] (*with default*):
#' If `TRUE`, y-axis is transformed to logarithmic scale.
#'
#' @param main [character] (*optional*):
#' Plot title, drawn at the top left of the diagram.
#'
#' @param xlab [character] (*optional*):
#' Axis title of the x-axis.
#'
#' @param ylab [character] (*optional*):
#' Axis title of the y-axis.
#'
#' @param xlim [numeric] vector (*optional*):
#' Minimum and maximum `c(min, max)` of the x-scale.
#'
#' @param ylim [numeric] vector (*optional*):
#' Minimum and maximum `c(min, max)` of the y-scale.
#'
#' @param font.size [numeric] (*with default*):
#' Scale factor for all text elements. Legend title and main title are one bigger
#'
#' @param graph.colors [color] vector (*optional*):
#' Color for the graphs/stacked areas are defined in the following order: 1. Measurement,
#' 2. Model, 3. Component 1, 4. Component 2, etc. The color vector is allowed to
#' be shorter than the needed colors. For missing colors, the default colors will be used
#'
#' @param graph.names [character] vector (*optional*):
#' Alternative graph names which shall be displayed in the legend. The names are defined
#' in the following order: 1. Measurement, 2. Model, 3. Residual, 4. Component 1, 5. Component 2, etc..
#' For missing names, the default names will be used.
#'
#'
#'
#'
#'
#'
#' @return
#' Returns an invisible [ggplot2::ggplot] object containing the plot "Invisible" means, the no value
#' will be returned (e.g. no console printout) if the function is not assigned to a variable via `<-`.
#' If the function is assigned, the returned object can be further manipulated by [ggplot2-package] methods
#' or manually drawn by various functions like for example [gridExtra::grid.arrange].
#'
#' @section Last update:
#'
#' 2024-08-29, DM: Forked this function from plot_OSLcurve
#'
#' @author
#' Dirk Mittelstraß, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please add the following references to your publication: Your currently used package version, obtained by
#' `citation("OSLdecomposition")`
#'
#' Mittelstraß, D., Schmidt, C., Kreutzer, S., Beyer, J., Straessner, A. and Heitmann, J.:
#' R package OSLdecomposition: Automated signal component analysis of multi-exponential decays for optically stimulated luminescence applications., *in preparation*.
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
#' @md
#' @export

plot_MultiExponential <- function(curve = NULL,
                             components = NULL,
                             fill.components = TRUE,
                             linear.modulation = FALSE,
                             xlog = FALSE,
                             ylog = FALSE,
                             main = NULL,
                             xlab = "Time",
                             ylab = "Signal",
                             xlim = NULL,
                             ylim = NULL,
                             font.size = 10,
                             graph.colors = NULL,
                             graph.names = NULL,
                             legend.position = "right",
                             legend.justification = NULL,
                             theme.set = ggplot2::theme_classic(),
                             hide.plot = FALSE){

  # Changelog:
  # * 2024-08-29, DM: First version of the new function

  #### INPUT CHECKS #### -------------------------------------------------------

  if (is.null(curve) && is.null(components)) stop("Either 'curve' or 'components' must be defined")

  # Check input components
  if (!is.null(components)) {

    if (inherits(components, "numeric")) {

      if (is.null(curve))
        stop("An input 'curve' is needed if 'components' input is just a numeric vector")

      components <- decompose_OSLcurve(curve,
                                       components,
                                       algorithm = "det+nls",
                                       error.estimation = "none",
                                       verbose = FALSE)

    } else if (inherits(components, "data.frame")) {

      if(!("name" %in% colnames(components)) ||
         !("lambda" %in% colnames(components)) ||
         !("n" %in% colnames(components)))
        stop("Input data.frame for 'components' needs at least the columns 'name', 'lambda' and 'n'")

      # Create curve from component table if not given
      if (is.null(curve)) {
        channel_width <- 1 / (2*max(components$lambda))
        channel_number <- ceiling(1 / (min(components$lambda) * channel_width))
        curve <- simulate_OSLcomponents(components,
                                        channel.width = channel_width,
                                        channel.number = channel_number,
                                        simulate.curve = TRUE,
                                        add.poisson.noise = FALSE)
      }

    } else {
      stop("Argument 'components' is not of type numeric or data.frame")
    }
  }

  # Check input curve
  if(!inherits(curve, c("RLum.Data.Curve", "data.frame", "matrix"))){
    stop("Input object 'curve' is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")}

  if(inherits(curve, "RLum.Data.Curve")) {

    # if no component table is given, check if the data set was already decomposed by RLum.OSL_decomposition
    if (is.null(components)) {
      if ("COMPONENTS" %in% names(curve@info)) components <- curve@info$COMPONENTS}

    # now transform the RLum object into a simple table. There is no need to keep all the extra info
    curve <- as.data.frame(Luminescence::get_RLum(curve))
  }

  if(inherits(curve, "matrix"))
    curve <- data.frame(time = curve[,1], signal = curve[,2])

  #### CREATE COMPONENT GRAPHS #### --------------------------------------------

  # Calculate components data points (signal values will not be overwritten)
  if (!is.null(components)){
    curve <- simulate_OSLcomponents(components, curve, simulate.curve = FALSE)
  } else {
    if (ncol(curve) > 2) {
      message("Input object 'curve' has more than two columns. Just the first two are used to define time and signal")
      curve <- curve[,1:2]
    }
  }

  K <- 0
  if (ncol(curve) > 4) K <- ncol(curve) - 4

  # Transform to linearly modulated curves here in accordance to Bulur (2000)
  if (linear.modulation) {

    P = 2*max(curve[,1])
    for (i in 2:ncol(curve))
      curve[,i] <- curve[,i] * sqrt(2 * curve[,1] / P)

    curve[,1] <- sqrt(2 * P * curve[,1])
  }


  # Build a new well defined object which will contain all graphs
  graphs <- data.frame(time = curve[,1],
                       signal = curve[,2])

  # Is there a "sum" curve?
  if (ncol(curve) > 2) {
    graphs <- cbind(graphs, curve[,3])
  }  else {
   # graphs <- cbind(graphs, NA)
  }

  # Is there a "residual" curve?
  if (ncol(curve) > 3) {
    graphs <- cbind(graphs, curve[,4])
  }  else {
  #  graphs <- cbind(graphs, NA)
  }

  # We will need these for the legend later
  graph_names <- c("Measurement", "Model", "Residual")

  # Stacked plots need special treatment to be displayed properly
  if (K > 0) {

    positive_comps <- negative_comps <- data.frame(NULL)
    if (K == 1) {

      if (components$n >= 0) {
        positive_comps <- as.data.frame(curve[,5])
        colnames(positive_comps) <- colnames(curve)[5]
      } else {
        negative_comps <- as.data.frame(curve[,5])
        colnames(negative_comps) <- colnames(curve)[5]
      }

    } else {

      # Separate components into two types: Increasing and decreasing
      #positive_comps <- curve[,-1:-4][, components$n >= 0]
      #negative_comps <- curve[,-1:-4][, components$n < 0]
      positve_ns <- components$n >= 0
      comps <- curve[,-1:-4]

      positive_comps <- as.data.frame(comps[, positve_ns])
      colnames(positive_comps) <- colnames(comps)[positve_ns]

      negative_comps <- as.data.frame(comps[, !positve_ns])
      colnames(negative_comps) <- colnames(comps)[!positve_ns]
    }

    # Re-iterate K to account for components which cannot be displayed
    K <- ncol(positive_comps)

    # Now stack up the components with positive n (decreasing graphs)
    # but put the SLOWEST component first. This way, it will always have
    # the same color, no matter if later fittings add some new
    # component or not

    if (ncol(positive_comps) > 0) {

      if (fill.components) {
        added_signal <- rep(0, nrow(curve))
        for (i in ncol(positive_comps):1) {

          # ymin boundary
          graphs <- cbind(graphs, added_signal)

          # ymax boundary
          added_signal <- added_signal + positive_comps[,i]
          graphs <- cbind(graphs, added_signal)
        }
      } else { # for "line" plots
        # rev() will revert the column order
        graphs <- cbind(graphs, rev(positive_comps))
      }
      graph_names <- c(graph_names, rev(colnames(positive_comps)))
    }

    # Do the same for components with n below zero (increasing graphs)
    if (ncol(negative_comps) > 0 && !ylog) {

    #  are_there_negative_values <- TRUE
      K <- K + ncol(negative_comps)

      if (fill.components) {
        previous_signal <- rep(0, nrow(curve))
        for (i in ncol(negative_comps):1) {

          # ymin boundary
          added_signal <- previous_signal + negative_comps[, i]
          graphs <- cbind(graphs, added_signal)

          # ymax boundary
          graphs <- cbind(graphs, previous_signal)
          previous_signal <- added_signal
        }
      } else {
        graphs <- cbind(graphs, rev(negative_comps))
      }
      graph_names <- c(graph_names, rev(colnames(negative_comps)))
    } else if(ncol(negative_comps) > 0 && ylog){
      message(paste("Input data contains signal components with negative intensity (increasing signals).",
                    "These can not be displayed with logarithmic y-axis and are therefore omitted."))
    }
  }


  #### CHECK AND SET AXES LIMITS #### ------------------------------------------

  if (!is.null(xlim)) {

    if (xlim[2] <= xlim[1])
      stop("X-axis minimum is larger then X-axis maximum.")

    if (xlog && (xlim[1] <= 0 || xlim[2] <= 0))
      stop("X-axis limits must be larger than zero when using logarithmic axis.")

    # Remove too small values
    graphs <- graphs[graphs$time >= xlim[1] ,]

    # Remove too large values
    graphs <- graphs[graphs$time <= xlim[2] ,]

    if (nrow(curve) < 1) stop("X-axis limits are too small for this data.")
  } else if(is.null(xlim) && xlog){

    # Check if there are values which won't work on logarithmic axis and remove them
    invalid_x <- curve[,1] <= 0
    if (sum(invalid_x) > 0) {
      message("X-axis contains values below/equal zero. Those are removed to enable logaritmic X-axis.")
      graphs <- graphs[!invalid_x,]
    }
  }

  if (!is.null(ylim)) {
    if (ylim[2] <= ylim[1])
      stop("Y-axis minimum is larger then Y-axis maximum.")

    if (ylog && (ylim[1] <= 0 || ylim[2] <= 0))
      stop("Y-axis limits must be larger than zero when using logarithmic axis.")

  } else if (is.null(ylim) && ylog){

    # The auto-limits extend too far to lower values to fully show all components.
    # Thus, we set some reasonable limits manually

    # Take also the sum curve into account, if one is defined
    signal_min <- min(c(curve[,2], curve$sum), na.rm = TRUE)
    signal_max <- max(c(curve[,2], curve$sum), na.rm = TRUE)

    if (signal_max <= 0)
      stop("Signal curve with just negative values can not be displayed with logarithmic Y-axis.")

    if (signal_min <= 0) {
      message("Signal curve contains value equal/lower zero. Values below 0.1 will not be displayed.")
      signal_min <- 0.1
    }

    ylim <- c(signal_min * 0.2, signal_max * 1.5)
  }

  #### DESIGN CHOICES #### -----------------------------------------------------

  library(ggplot2)

  # Set color and line themes
  ggplot2::theme_set(theme.set)

  # Set colors and legend names
  if (length(graph.names) > 0)
    graph_names[1:length(graph.names)] <- graph.names

  graph_colors <- c("grey30", "darkblue", "skyblue2","orchid","cyan2","orange","red2","pink3","brown2")

  if (length(graph.colors) > 0)
    graph_colors[1:length(graph.colors)] <- graph.colors

  # Use signal color also for residual graph
  graph_colors <- append(graph_colors, graph_colors[1], after = 2)

  graph_colors <- graph_colors[1:length(graph_names)]
  names(graph_colors) <- graph_names

  # Do this to prevent generic col names which might throw an exception in ggplot2
  if(ncol(graphs) > 3)
    colnames(graphs)[3:ncol(graphs)] <- letters[1:(ncol(graphs) - 2)]


  # Ensure that all titles begin with a upper case character
  #for (i in 1:length(curve_params)) {
  #  substr(curve_params[i], 1, 1) <- toupper(substr(curve_params[i] , 1, 1))}


  #### ADD COMPONENT PLOTS #### ------------------------------------------------

  # We start with an empty plot to keep full control over everything
  p <- ggplot(graphs, aes(x = time))

  # Are there any columns besides Signal, Time, Sum and Residual?
  # Then these are the components. Draw them first!
  if (K > 0) {

    if (K >= 8)
      warning("Graphs with more than 7 components are not supported. Only the first 7 are displayed.")

    if (fill.components) {

      # Yes, the code looks weird. However, it seems like ggplot does not copy object, instead it
      # of the input data. Instead it works address pointer. Thus, dynamical approaches like increasing
      # indices won't work.

      if (K >= 1) p <- p +
          geom_ribbon(aes(ymin = graphs[,5],  ymax = graphs[,6], fill = graph_names[4]))

      if (K >= 2) p <- p +
          geom_ribbon(aes(ymin = graphs[,7],  ymax = graphs[,8], fill = graph_names[5]))

      if (K >= 3) p <- p +
          geom_ribbon(aes(ymin = graphs[,9],  ymax = graphs[,10], fill = graph_names[6]))

      if (K >= 4) p <- p +
          geom_ribbon(aes(ymin = graphs[,11],  ymax = graphs[,12], fill = graph_names[7]))

      if (K >= 5) p <- p +
          geom_ribbon(aes(ymin = graphs[,13],  ymax = graphs[,14], fill = graph_names[8]))

      if (K >= 6) p <- p +
          geom_ribbon(aes(ymin = graphs[,15],  ymax = graphs[,16], fill = graph_names[9]))

      if (K >= 7) p <- p +
          geom_ribbon(aes(ymin = graphs[,17],  ymax = graphs[,18], fill = graph_names[10]))

    } else { # line graphs

      if (K >= 1) p <- p +
          geom_line(aes(y = graphs[,5], color = graph_names[4]), linewidth = 0.75)

      if (K >= 2) p <- p +
          geom_line(aes(y = graphs[,6], color = graph_names[5]), linewidth = 0.75)

      if (K >= 3) p <- p +
          geom_line(aes(y = graphs[,7], color = graph_names[6]), linewidth = 0.75)

      if (K >= 4) p <- p +
          geom_line(aes(y = graphs[,8], color = graph_names[7]), linewidth = 0.75)

      if (K >= 5) p <- p +
          geom_line(aes(y = graphs[,9], color = graph_names[8]), linewidth = 0.75)

      if (K >= 6) p <- p +
          geom_line(aes(y = graphs[,10], color = graph_names[9]), linewidth = 0.75)

      if (K >= 7) p <- p +
          geom_line(aes(y = graphs[,11], color = graph_names[10]), linewidth = 0.75)
    }
  }


  #### BASIC PLOT #### ---------------------------------------------------------

  # If there is any value with negative sign, draw a zero line. Time and residual does not count
  if (any(graphs[,c(-1, -4)] < 0, na.rm = TRUE)) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
  }

  # Signal curve
 # p <- p + geom_point(aes(y = graphs[, 2], color = "Measurement"), size = 1)

  p <- p + geom_line(aes(y = graphs[, 2], color = graph_names[1]), linewidth = 0.5)
  # Summed up components
  if (ncol(graphs) > 2) {
    p <- p + geom_line(aes(y = graphs[, 3], color = graph_names[2]), linewidth = 0.75)
  }


  #### SET SCALES #### ---------------------------------------------------------

  if (ylog) {
    p <- p + scale_y_log10(limits = ylim)
  } else {
    p <- p + scale_y_continuous(limits = ylim)
  }

  if (xlog) {
    p <- p + scale_x_log10(limits = xlim)
  } else {
    p <- p + scale_x_continuous(limits = xlim)
  }

  #### APPLY DESIGN SETTINGS #### ----------------------------------------------

  # Apply axis labels and design choices
  if (is.null(xlab)) xlab <- "Time"
  if (is.null(ylab)) ylab <- "Signal"

  # Set legend, axis labels and design choices
  p <- p +
    scale_color_manual(values = graph_colors) +
    scale_fill_manual(values = graph_colors) +
    labs(color = NULL, #"Multi-exponential curve",
         fill = "Signal components",
         subtitle = main, x = xlab, y = ylab) +
    theme(axis.title = element_text(size = font.size),
          element_text(size = font.size + 1, face = "bold"),
          legend.position = legend.position,
          legend.justification = legend.justification,
          legend.text = element_text(size = font.size))

  #### RETURN OBJECTS #### -----------------------------------------------------

  # show plot
  if (!hide.plot){

    # ggplot will throw a lot of annoying warnings when data points are outside
    # of the limits, thus we suppress them.
    # However, better would be to set them NA
    suppressWarnings(print(p))
  }

  # return plot object
  invisible(p)
}
