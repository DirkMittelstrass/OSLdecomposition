#' Advanced plot function for displaying component resolved signal curves
#'
#'
#'
#'
#' @return
#' An invisible [ggplot2::ggplot] object containing the diagram will returned. "Invisible" means, the no value
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
#'
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
                             xlab = NULL,
                             ylab = NULL,
                             xlim = NULL,
                             ylim = NULL,
                             font.size = 10,
                             component.colors = NULL,
                             component.names = NULL,
                             legend.position = "right",
                             legend.justification = NULL,
                             hide.plot = FALSE,
                             theme.set = ggplot2::theme_classic()){

  # Changelog:
  # * 2024-08-29, DM: First version of the new function

  #### INPUT CHECKS #### -------------------------------------------------------

  if (is.null(curve) && is.null(components)) stop("Either 'curve' or 'components' must be defined")

  if (!is.null(components)) {

    # Add decay rate input here
  }

  # Create curve from component table if not given
  if (is.null(curve) && !is.null(components)) {

    channel_width <- 1 / (2*max(components$lambda))
    channel_number <- ceiling(1 / (min(components$lambda) * channel_width))
    curve <- simulate_OSLcomponents(components,
                                    channel.width = channel_width,
                                    channel.number = channel_number,
                                    simulate.curve = TRUE,
                                    add.poisson.noise = FALSE)
  } else {

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
  }

  #### CREATE COMPONENT GRAPHS #### --------------------------------------------

  if (!is.null(components))
    curve <- simulate_OSLcomponents(components, curve, simulate.curve = FALSE)

  K <- 0
  if (ncol(curve) > 4) K <- ncol(curve) - 4

  # Transform to linearly modulated curves here in accordance to Bulur (2000)
  if (linear.modulation) {
    # transform data
    P = 2*max(curve[,1])
    for (i in 2:ncol(curve)) {
      curve[,i] <- curve[,i] * sqrt(2 * curve[,1] / P)
    }
    curve[,1] <- sqrt(2 * P * curve[,1])
  }




  # Build a new well defined object which will contain all graphs
  graphs <- data.frame(time = curve[,1],
                       signal = curve[,2])

  # Is there a "sum" curve?
  if (ncol(curve) > 2) {
    graphs <- cbind(graphs, curve[,3])
  }  else {
    graphs <- cbind(graphs, NA)
  }

  # Is there a "residual" curve?
  if (ncol(curve) > 3) {
    graphs <- cbind(graphs, curve[,4])
  }  else {
    graphs <- cbind(graphs, NA)
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
  #  curve <- curve[curve[,1] >= xlim[1] ,]
    graphs <- graphs[graphs$time >= xlim[1] ,]

    # Remove too large values
  #  curve <- curve[curve[,1] <= xlim[2] ,]
    graphs <- graphs[graphs$time <= xlim[2] ,]

    if (nrow(curve) < 1) stop("X-axis limits are too small for this data.")
  } else if(is.null(xlim) && xlog){

    # Check if there are values which won't work on logarithmic axis and remove them
    invalid_x <- curve[,1] <= 0
    if (sum(invalid_x) > 0) {
      message("X-axis contains values below/equal zero. Those are removed to enable logaritmic X-axis.")
  #    curve <- curve[!invalid_x,]
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
  if (length(component.names) > 0)
    graph_names[4:(3 + length(component.names))] <- component.names

  #graph_names <- c("Measurement", "Model", "Residual", "C1", "C2", "C3", "C4", "C5", "C6", "C7")
  graph_colors <- c("grey40", "black", "grey40", "skyblue2","orchid","cyan2","orange","red2","pink3","brown2")

  if (length(component.colors) > 0)
    graph_colors[4:length(component.colors)] <- component.colors

  graph_colors <- graph_colors[1:length(graph_names)]
  names(graph_colors) <- graph_names

  # Do this to prevent generic col names which might throw an exception in ggplot2
  colnames(graphs)[3:length(graphs)] <- letters[1:(length(graphs)-2)]


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

      #return(graphs)

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
          geom_line(aes(y = graphs[,5], color = graph_names[4]), size = 0.75)

      if (K >= 2) p <- p +
          geom_line(aes(y = graphs[,6], color = graph_names[5]), size = 0.75)

      if (K >= 3) p <- p +
          geom_line(aes(y = graphs[,7], color = graph_names[6]), size = 0.75)

      if (K >= 4) p <- p +
          geom_line(aes(y = graphs[,8], color = graph_names[7]), size = 0.75)

      if (K >= 5) p <- p +
          geom_line(aes(y = graphs[,9], color = graph_names[8]), size = 0.75)

      if (K >= 6) p <- p +
          geom_line(aes(y = graphs[,10], color = graph_names[9]), size = 0.75)

      if (K >= 7) p <- p +
          geom_line(aes(y = graphs[,11], color = graph_names[10]), size = 0.75)
    }
  }


  #### BASIC PLOT #### ---------------------------------------------------------

  # If there is any value with negative sign, draw a zero line. Time and residual does not count
  if (any(graphs[,c(-1, -4)] < 0, na.rm = TRUE)) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
  }

  # Signal curve
  p <- p + geom_point(aes(y = graphs[, 2], color = "Measurement"), size = 1)

  # Summed up components
  if (ncol(graphs) > 2) {
    p <- p + geom_line(aes(y = graphs[, 3], color = "Model"), size = 0.75)
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
  if (is.null(xlab)) xlab <- graph_names[1]
  if (is.null(ylab)) ylab <- graph_names[2]

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
