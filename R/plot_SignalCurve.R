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

  # channel.width <- curve$time[2] - curve$time[1]

  #### PLOT DESIGN #### --------------------------------------------------------

  library(ggplot2)

  # Set color and line themes
  ggplot2::theme_set(theme.set)
  comp_col <- c("royalblue3","green3","red3","gold","orchid","orange","pink2","brown2")
  if (length(colors) > 0) comp_col[1:length(colors)] <- colors

  # Reduce font size of titles; also set the legend design
  text_format <- ggplot2::theme(axis.title = ggplot2::element_text(size = 8),
                                plot.subtitle = ggplot2::element_text(size = 9, face = "bold"),
                                legend.position.inside = c(1, 1), legend.justification = c("right", "top"),
                                legend.title = ggplot2::element_blank(),
                                legend.text = ggplot2::element_text(size = 8))

  #### BASIC PLOT #### ---------------------------------------------------------

  p <- ggplot(curve, aes(x = time, y = signal)) +
    geom_line(color = "darkgrey", linewidth = 0.5)

  #### SET SCALES #### ---------------------------------------------------------

  if (scaling.type == "log") {
    p <- p + scale_y_log10()
  }


  #### APPLY DESIGN SETTINGS #### ----------------------------------------------

  p <- p + text_format

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
