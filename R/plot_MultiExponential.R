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
                             scaling.type = "linear", #" log", "loglog", "LM"
                             component.type = "fill", # "line", "pattern"
                             fill.components = TRUE,
                             linear.modulation = FALSE,
                             show.legend = TRUE,
                             show.residuals = TRUE,
                             colors = NULL,
                             font.size = 10,
                             main = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             xlim = NULL,
                             ylim = NULL,
                             xlog = FALSE,
                             ylog = FALSE,
                             hide.plot = FALSE,
                             filename = NULL,
                             theme.set = ggplot2::theme_classic()){

  # Changelog:
  # * 2024-08-29, DM: First version of the new function

  # Hidden parameters, might be removed
  show.intervals <- FALSE
  show.crosssec <- FALSE

  #### INPUT CHECKS #### -------------------------------------------------------

  if (is.null(curve) && is.null(components)) stop("[plot_MultiExponential()]: Either 'curve' or 'components' must be defined")

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
      stop("[plot_MultiExponential()] Error: Input object 'curve' is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")}

    if(inherits(curve, "RLum.Data.Curve")) {

      # if no component table is given, check if the data set was already decomposed by RLum.OSL_decomposition
      if (is.null(components)) {
        if ("COMPONENTS" %in% names(curve@info)) components <- curve@info$COMPONENTS}

      # now transform the RLum object into a simple table. There is no need to keep all the extra info
      curve <- as.data.frame(Luminescence::get_RLum(curve))
    }

    if (!("time" %in% colnames(curve)) | !("signal" %in% colnames(curve))) {
      curve <- data.frame(time = curve[,1], signal = curve[,2])}
  }

  #### CREATE COMPONENT GRAPHS #### --------------------------------------------

  if (!is.null(components)) {
    curve <- simulate_OSLcomponents(components, curve, simulate.curve = FALSE)
  } else {
    # warning("[plot_MultiExponential()] When setting argument 'component.type' to 'fill', the argument 'components' must also be set. Thus, changed  'component.type' to 'line' ")
    # component.type <- "line"
  }

  # channel.width <- curve$time[2] - curve$time[1]
  curve_params <- colnames(curve)
  K <- 0
  if (ncol(curve) > 4) K <- ncol(curve) - 4

  # Transform to linearly modulated curves here in accordance to Bulur (2000)
  if (linear.modulation) {
    # transform data
    P = 2*max(curve$time)
    for (i in 2:ncol(curve)) {
      curve[,i] <- curve[,i] * sqrt(2 * curve$time / P)
    }
    curve$time <- sqrt(2 * P * curve$time)
  }


  # Stacked plots need special treatment to be displayed properly
  c_plots <- data.frame(time = curve$time)
  are_there_negative_ns <- FALSE
  if (K > 0) {

    # Separate components into two types: Increasing and decreasing
    comps_only <- curve[,5:ncol(curve)]
    positive_comps <- as.data.frame(comps_only[, components$n >= 0])
    negative_comps <- as.data.frame(comps_only[, components$n < 0])

    # Now stack up the components with positive n (decreasing graphs)
    # but put the SLOWEST component first. This way, it will always have
    # the same color, no matter if later fittings add some new
    # component or not
    if (ncol(positive_comps) > 0) {

      if (fill.components) {
        added_signal <- rep(0, nrow(curve))
        for (i in ncol(positive_comps):1) {

          # ymin boundary
          c_plots <- cbind(c_plots, added_signal)

          # ymax boundary
          added_signal <- added_signal + positive_comps[,i]
          c_plots <- cbind(c_plots, added_signal)
        }
      } else { # for "line" plots
        # rev() will revert the column order
        c_plots <- cbind(c_plots, rev(positive_comps))
      }
    }

    # Do the same for components with n below zero (increasing graphs)
    if (ncol(negative_comps) > 0) {
      are_there_negative_ns <- TRUE

      if (fill.components) {
        previous_signal <- rep(0, nrow(curve))
        for (i in ncol(negative_comps):1) {

          # ymin boundary
          added_signal <- previous_signal + negative_comps[, i]
          c_plots <- cbind(c_plots, added_signal)

          # ymax boundary
          c_plots <- cbind(c_plots, previous_signal)
          previous_signal <- added_signal
        }
      } else {
        c_plots <- cbind(c_plots, rev(negative_comps))
      }
    }
  }


  #### CHECK AND SET AXES LIMITS #### ------------------------------------------

  if (!is.null(xlim)) {
    if (xlim[2] <= xlim[1]) stop("X-axis minimum is larger then X-axis maximum ")
    if (scaling.type == "loglog" &&
        (xlim[1] <= 0 || xlim[2] <= 0)) stop("Negative X-axis limits are not allowed when using logarithmic axis")

    # Remove too small values
    curve <- curve[curve$time >= xlim[1] ,]
    c_plots <- c_plots[c_plots$time >= xlim[1] ,]

    # Remove too large values
    curve <- curve[curve$time <= xlim[2] ,]
    c_plots <- c_plots[c_plots$time <= xlim[2] ,]

    if (nrow(curve) < 1) stop("X-axis limits are too small for this data")
  }



  #### DESIGN CHOICES #### -----------------------------------------------------

  library(ggplot2)

  # Set color and line themes
  ggplot2::theme_set(theme.set)
  comp_col <- c("skyblue2","orchid","cyan2","orange","red2","pink3","brown2")
  if (length(colors) > 0) comp_col[1:length(colors)] <- colors

  # Reduce font size of titles; also set the legend design
  text_format <- ggplot2::theme(axis.title = ggplot2::element_text(size = font.size),
                                plot.subtitle = ggplot2::element_text(size = font.size + 1, face = "bold"),
                                legend.position.inside = c(1, 1), legend.justification = c("right", "top"),
                                legend.title = ggplot2::element_blank(),
                                legend.text = ggplot2::element_text(size = font.size))





  #### ADD COMPONENT PLOTS #### ------------------------------------------------

  # We start with an empty plot to keep full control over everything
  p <- ggplot()

  # Are there any columns besides Signal, Time, Sum and Residual?
  # Then these are the components. Draw them first!
  if (K > 0) {

    if (K >= 8) warning("[plot_MultiExponential()] Graphs with more than 7 components are not supported. Only the first 7 are displayed.")

    col_index <- length(curve_params) - 4

    if (fill.components) {

      # Yes, the code looks weird. However, it seems like ggplot does not copy object, instead it
      # of the input data. Instead it works address pointer. Thus, dynamical approaches like increasing
      # indices won't work.

      if (K >= 1) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,2],  ymax = c_plots[,3]), colour = comp_col[1], fill = comp_col[1])

      if (K >= 2) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,4],  ymax = c_plots[,5]), colour = comp_col[2], fill = comp_col[2])

      if (K >= 3) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,6],  ymax = c_plots[,7]), colour = comp_col[3], fill = comp_col[3])

      if (K >= 4) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,8],  ymax = c_plots[,9]), colour = comp_col[4], fill = comp_col[4])

      if (K >= 5) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,10],  ymax = c_plots[,11]), colour = comp_col[5], fill = comp_col[5])

      if (K >= 6) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,12],  ymax = c_plots[,13]), colour = comp_col[6], fill = comp_col[6])

      if (K >= 7) p <- p +
          geom_ribbon(aes(x = c_plots[,1], ymin = c_plots[,14],  ymax = c_plots[,15]), colour = comp_col[7], fill = comp_col[7])

    } else { # line graphs

      if (K >= 1) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,2]), color = comp_col[1], linewidth = 0.75)

      if (K >= 2) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,3]), color = comp_col[2], linewidth = 0.75)

      if (K >= 3) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,4]), color = comp_col[3], linewidth = 0.75)

      if (K >= 4) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,5]), color = comp_col[4], linewidth = 0.75)

      if (K >= 5) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,6]), color = comp_col[5], linewidth = 0.75)

      if (K >= 6) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,7]), color = comp_col[6], linewidth = 0.75)

      if (K >= 7) p <- p +
          geom_line(aes(x = c_plots[,1], y = c_plots[,8]), color = comp_col[7], linewidth = 0.75)
    }
  }


  #### BASIC PLOT #### ---------------------------------------------------------

  if (are_there_negative_ns) {
   # p <- p + geom_hline(yintercept = 0, color = "white", linewidth = 0.5)
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5)
  }

  # Signal curve
  p <- p + geom_point(aes(x = curve[, 1], y = curve[, 2]), color = "grey40", size = 1)

  # Summed up components
  if (ncol(curve) > 2) {
    p <- p + geom_line(aes(x = curve[, 1], y = curve[, 3]), color = "black", linewidth = 0.75)
  }

  #### RESIDUAL PLOT #### ------------------------------------------------------

  if (show.residuals && ncol(curve) >= 4) {




  }

  #### SET SCALES #### ---------------------------------------------------------

  if (scaling.type == "log") {
    p <- p + scale_y_log10()
  }


  #### APPLY DESIGN SETTINGS #### ----------------------------------------------


  # Ensure that all titles begin with a upper case character
  for (i in 1:length(curve_params)) {
    substr(curve_params[i], 1, 1) <- toupper(substr(curve_params[i] , 1, 1))}

  if (is.null(xlab)) xlab <- curve_params[1]
  if (is.null(ylab)) ylab <- curve_params[2]

  p <- p +
    labs(subtitle = main, x = xlab, y = ylab) +
    text_format

  #### RETURN OBJECTS #### -----------------------------------------------------

  # save plot as file
  if (!is.null(filename)) {

    try(suppressMessages(ggplot2::ggsave(filename, plot = p, units = "cm")),
        silent = FALSE)}

  # show plot
  if (!hide.plot) print(p)

  # return plot object
  invisible(p)
}
