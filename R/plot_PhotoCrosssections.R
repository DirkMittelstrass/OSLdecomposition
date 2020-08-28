#' Plot comparison of the global OSL curve photo-crosssections with literature values
#'
#' @param fit.list
#'
#' @section Changelog:
#' * 2019-10-08 Seperated from fit_OSLcurve()
#' * 2020-04-09 Changed from decay rates to photoionisation cross-sections; cleaned code
#'
#' @section ToDo:
#' * Documentation
#' * Sometimes the red rectangle is not drawn (fit_OSLcurve(Batagai, K.max = 6)). Why?
#' * add literature values for other than just 470nm values
#' * remove library(XXXX)
#'
#' @section Last changed. 2020-04-09
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' @return
#' @export
#'
#' @examples
#'
plot_PhotoCrosssections <- function(
  fit.list,
  stimulation.intensity,
  stimulation.wavelength,
  title = NULL
){

  ##### Load graphic libraries and set output theme #####
  library(gridExtra)
 # library(scales)
  library(ggplot2)
  theme_set(theme_bw())

  # supress warnings in the whole script
  options( warn = -1 )

  K <- fit.list$K.selected
  plot_data <- fit.list$plot.data
  x <- nrow(fit.list$F.test)

  # Calc photon energy: E = h*v  [W*s^2 * s^-1 = W*s = J]
  E <-6.62606957e-34 * 299792458 * 10^9 / stimulation.wavelength

  # Calc photon flux of stimulation light: Flux = I / E  [W/cm^2 / W*s = 1/s*cm^2]
  Flux <- stimulation.intensity / (E * 1000)

  # Transform decay rates into photoionisation cross-sections
  plot_data$lambda <- plot_data$lambda / Flux
  plot_data$lambda.low <- plot_data$lambda.low / Flux
  plot_data$lambda.up <- plot_data$lambda.up / Flux

  # Plot literature values just in case of blue light stimulation
  plot_literature <- FALSE
  if ((stimulation.wavelength >= 465) && (stimulation.wavelength <= 480)) plot_literature <- TRUE

  if(plot_literature) {

    # set drawing borders for fast component (1:2) and medium component (3:4)
    y.components <- c(1.9e-17, 3.1e-17,
                      3e-18, 9e-18,
                      1e-18, 1.85e-18,
                      1.1e-19, 4e-19,
                      1e-20, 4.67e-20)
    y.means <- c((y.components[1] - y.components[2])/ log(y.components[1] / y.components[2]),
                 (y.components[3] - y.components[4])/ log(y.components[3] / y.components[4]),
                 (y.components[5] - y.components[6])/ log(y.components[5] / y.components[6]),
                 (y.components[7] - y.components[8])/ log(y.components[7] / y.components[8]),
                 (y.components[9] - y.components[10])/ log(y.components[9] / y.components[10]))

     # Increase number of plotting rows
    x <- x + 1

    # Durcan and Duller 2011
    plot_data <- rbind(plot_data,
                       data.frame(lambda = c(2.6e-17, 4.28e-18, 1.09e-18, 3.04e-19, 3.39e-20, 9.06e-21),
                                  lambda.low = c(2.54e-17, 3.93e-18, 1.04e-18, 2.58e-19, 2.73e-20, 8.31e-21),
                                  lambda.up = c(2.66e-17, 4.62e-18, 1.14e-18, 3.5e-19, 4.03e-20, 9.81e-21),
                                  name = factor("Durcan & Duller (2011)"),
                                  x = x))
    x <- x + 1

    # Jain et al. 2003
    plot_data <- rbind(plot_data,
                       data.frame(lambda = c(2.9e-16, 2.32e-17, 5.59e-18, 1.33e-18, 2.08e-19, 2.06e-20, 2.76e-21),
                                  lambda.low = c(2.9e-16,2.16e-17, 5.15e-18, 1.07e-18, 1.62e-19, 1.9e-20, 2.59e-21),
                                  lambda.up = c(2.9e-16,2.48e-17, 6.03e-18, 1.59e-18, 2.54e-19, 2.22e-20, 2.93e-21),
                                  name = factor("Jain et al. (2003)"),
                                  x = x))
    x <- x + 1

    # Singarayer and Bailey 2003
    plot_data <- rbind(plot_data,
                       data.frame(lambda = c(2.5e-17, 5.9e-18, 2.1e-19, 1.2e-20, 1.9e-21),
                                  lambda.low = c(2.2e-17, 3.9e-18, 1.6e-19, 1.0e-20, 0.7e-21),
                                  lambda.up = c(2.8e-17, 7.9e-18, 2.6e-19, 1.4e-20, 4.7e-21),
                                  name = factor("Singarayer & Bailey (2003)"),
                                  x = x))
  }

  # set x-Axis
  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
    }

  breaks <- 10^(-25:-15)
  minor_breaks <- rep(1:9, 21)*(10^rep(-25:-15, each=9))
  xmax <- x

  # set plotting band width
  ymin <- min(plot_data$lambda.low) / 1.5
  ymax <- max(plot_data$lambda.up) * 1.5

  ##### build F-test plot #####

  # Draw decay rates from F-table as diagram
  p <- ggplot(plot_data, aes(x = x,
                             y = lambda)) +
    geom_point(size = 2.5) + geom_errorbar(aes(ymin = lambda.low, ymax = lambda.up), width = .3) +
    scale_y_continuous(trans = reverselog_trans(10),
                       limits = c(ymax, ymin),
                       breaks = breaks,
                       minor_breaks = minor_breaks) +
    scale_x_reverse(labels = levels(plot_data$name),
                    breaks = 1:xmax,
                    minor_breaks = 1:(xmax + 1) - 0.5)

  # Draw component definition ranges from literature
  if(plot_literature) {p <- p +  annotate("rect", xmin = 0.5, xmax = xmax + 1.5,  #################
             ymin = y.components[1], ymax = y.components[2], alpha = 0.1, fill = "red3") +
    annotate("text", x = xmax + 1, y = y.means[1], label = "Fast", colour = "red3") +
    annotate("rect", xmin = 0.5, xmax = xmax + 1.5, #################
             ymin = y.components[3], ymax = y.components[4], alpha = 0.1, fill = "green3") +
    annotate("text", x = xmax + 1, y = y.means[2], label = "Medium", colour = "green3") +
    annotate("rect", xmin = 0.5, xmax = xmax + 1.5, #################
             ymin = y.components[5], ymax = y.components[6], alpha = 0.1, fill = "blue3") +
    annotate("text", x = xmax + 1, y = y.means[3], label = "Slow1", colour = "blue3") +
    annotate("rect", xmin = 0.5, xmax = xmax + 1.5, #################
             ymin = y.components[7], ymax = y.components[8], alpha = 0.1, fill = "darkorchid") +
    annotate("text", x = xmax + 1, y = y.means[4], label = "Slow2", colour = "darkorchid") +
    annotate("rect", xmin = 0.5, xmax = xmax + 1.5, #################
             ymin = y.components[9], ymax = y.components[10], alpha = 0.1, fill = "gold") +
    annotate("text", x = xmax + 1, y = y.means[5], label = "Slow3", colour = "gold") +
    geom_vline(mapping = aes(xintercept = xmax - 2.5), color="black")
  }

  # Draw red rectangle which indicates the selected case
  p <- p + geom_rect(mapping = aes(xmin = K - 0.5, xmax = K + 0.5,
                                   ymin = ymin, ymax = ymax),
                     colour = "red", size = 1, fill = NA)

  # rotate diagramm and delete unnecessary visual elements
  p <- p + coord_flip() +
    ylab(expression("photoionisation cross-section (cm^2)")) +
    theme(panel.border = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(colour = "black", size = 10),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.y = element_blank())

  # show plot
  grid.arrange(p, nrow = 1, top = title)
}
