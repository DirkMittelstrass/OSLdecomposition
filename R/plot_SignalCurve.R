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

plot_SignalCurve <- function(curve = NULL,
                             components = NULL,
                             scaling.type = "linear", #" log", "loglog", "LM"
                             component.type = "fill", # "line", "pattern"
                             show.legend = TRUE,
                             show.residuals = TRUE,
                             colors = NULL,
                             font.size = 10,
                             main = NULL,
                             xlab = NULL,
                             ylab = NULL,
                             xlim = NULL,
                             ylim = NULL,
                             hide.plot = FALSE,
                             filename = NULL,
                             theme.set = ggplot2::theme_classic()){

  # Changelog:
  # * 2024-08-29, DM: Created function

  # Hidden parameters, might be removed
  show.intervals <- FALSE
  show.crosssec <- FALSE

  #### INPUT CHECKS #### -------------------------------------------------------

  if (is.null(curve) && is.null(components)) stop("[plot_SignalCurve()]: Either 'curve' or 'components' must be defined")

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
      stop("[plot_SignalCurve()] Error: Input object 'curve' is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")}

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

  #### ADD COMPONENT SIGNALS #### ----------------------------------------------

  if (!is.null(components)) {
    curve <- simulate_OSLcomponents(components, curve, simulate.curve = FALSE)
  } else {
    warning("[plot_SignalCurve()] When setting argument 'component.type' to 'fill', the argument 'components' must also be set. Thus, changed  'component.type' to 'line' ")
    component.type <- "line"
  }


  # channel.width <- curve$time[2] - curve$time[1]

  curve_params <- colnames(curve)

  # Stacked plots need special treatment to be displayed properly
  positive_comps <- negative_comps <- data.frame(NULL)
  if (component.type == "fill" && length(curve_params) > 4) {

    # Separate components into two types: Increasing and decreasing
    comps_only <- curve[,5:ncol(curve)]
    positive_comps <- comps_only[, components$n >= 0]
    negative_comps <- comps_only[, components$n < 0]

    # Now stack up the components
    if (ncol(positive_comps) > 1) {
      for (i in ncol(positive_comps):2) {
        positive_comps[,i - 1] <- positive_comps[,i - 1] + positive_comps[,i]
      }

      # Add zero line which is needed by ggplot_ribbon()
      positive_comps$zero <- 0
    }

    #return(positive_comps)

    # Doe the same for the components below zero
    if (ncol(negative_comps) > 1) {
      for (i in ncol(negative_comps):2) {
        negative_comps[,i - 1] <- negative_comps[,i - 1] + negative_comps[,i]
      }
      negative_comps$zero <- 0
    }
  }



  #### DESIGN CHOICES #### -----------------------------------------------------

  library(ggplot2)

  # Set color and line themes
  ggplot2::theme_set(theme.set)
  comp_col <- c("orchid","royalblue2","red2","orange","green3","pink2","gold","brown2")
  if (length(colors) > 0) comp_col[1:length(colors)] <- colors

  # Reduce font size of titles; also set the legend design
  text_format <- ggplot2::theme(axis.title = ggplot2::element_text(size = font.size),
                                plot.subtitle = ggplot2::element_text(size = font.size + 1, face = "bold"),
                                legend.position.inside = c(1, 1), legend.justification = c("right", "top"),
                                legend.title = ggplot2::element_blank(),
                                legend.text = ggplot2::element_text(size = font.size))





  #### ADD COMPONENT PLOTS #### ------------------------------------------------

  p <- ggplot(curve)

  # Are there any columns besides Signal, Time, Sum and Residual?
  # Then these are the components. Draw them first!
  if (length(curve_params) > 4) {

    col_index <- length(curve_params) - 4

    if (component.type == "fill") {

      # QUESTION: Must  ymin/ymax be part of the input data?
      if (ncol(positive_comps) > 1) {
        for (i in ncol(positive_comps):2) {
          p <- p + geom_ribbon(aes(x = time,
                                   ymin = positive_comps[, i],
                                   ymax = positive_comps[, i - 1]),
                               color = comp_col[col_index], fill = comp_col[col_index])
          col_index <- col_index - 1
        }
      }

      if (ncol(negative_comps) > 1) {
        for (i in ncol(negative_comps):2) {
          p <- p + geom_ribbon(aes(x = time,
                                   ymax = negative_comps[, i],
                                   ymin = negative_comps[, i - 1]),
                               color = comp_col[col_index], fill = comp_col[col_index])
          col_index <- col_index - 1
        }
      }





    } else {

      for (i in 5:length(curve_params)) {
        p <- p + geom_line(aes(x = time, y = curve[, i]), color = comp_col[col_index], linewidth = 0.5)
        col_index <- col_index - 1
      }
    }
  }


  #### BASIC PLOT #### ---------------------------------------------------------

  p <- p + geom_line(aes(x = time, y = signal), color = "grey40", linewidth = 0.5)


  #### SET SCALES #### ---------------------------------------------------------

  if (scaling.type == "log") {
    p <- p + scale_y_log10()
  }


  #### APPLY DESIGN SETTINGS #### ----------------------------------------------

  p <- p + text_format

  # Ensure that all titles begin with a upper case character
  for (i in 1:length(curve_params)) {
    substr(curve_params[i], 1, 1) <- toupper(substr(curve_params[i] , 1, 1))}

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
