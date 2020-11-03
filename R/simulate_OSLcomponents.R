#' Simulates signal component decay curves and whole CW-OSL curves
#'
#' The function builds CW-OSL component decay curves or a whole CW-OSL curve from OSL component parameters.
#' Therewith it supports [fit_OSLcurve], [decompose_OSLcurve] and [plot_OSLcurve] by providing model and residual curves.
#'
#'
#' @param components [data.frame] (**required**):
#' Table with component parameters. Need to have columns `$names`, `$lambda` and `$n`, see section **Examples**
#'
#' @param curve [data.frame] (*optional*):
#' CW-OSL curve serving as template for the time axis. Need to have a column `$time`.
#' If no input object is given or the object contains no column `$signal`,
#' `simulate.curve` will be set `TRUE`
#'
#' @param channel.width [numeric] (*with default*):
#' channel width in seconds if no template `curve` is given
#'
#' @param channel.number [numeric] (*with default*):
#' number of channels resp. data points if no template `curve` is given
#'
#' @param simulate.curve [logical] (*with default*):
#' if `FALSE`, the output curve will take over the column `$signal` from the input `curve`.
#' If `TRUE`, a new column `$signal` will be created which is the sum of all component curves
#'
#' @param add.poisson.noise [logical] (*with default*):
#' adds poisson distributed shot noise to `$signal` if `simulate.curve = TRUE`
#'
#' @param add.gaussian.noise [numeric] (*with default*):
#' standard deviation of the detector noise in **cts/s**, added to `$signal` if `simulate.curve = TRUE`
#'
#' @param add.background [numeric] (*with default*):
#' signal background level in **cts/s**, added to `$signal` if `simulate.curve = TRUE`
#'
#' @param round.values [logical] (*with default*):
#' rounds `$signal` values to integers if `simulate.curve = TRUE`
#'
#' @return
#'
#' A [data.frame]  of a CW-OSL curve with the columns: `$time`, `$signal`, `$residual`, `$sum`
#' and a signal decay curve for each single component named after the input object `components$names`
#'
#' @section Last updates:
#'
#' 2020-10-30, DM: Renamed from *simulate_OSLcurve* to *simulate_OSLcomponents*;
#' Renamed argument from *template.curve* to *curve*; Rewrote roxygen documentation
#'
#' @author
#'
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' Please cite the package the following way:
#'
#' Mittelstraß, D., Schmidt, C., Beyer, J., Heitmann, J. and Straessner, A.:
#' Automated identification and separation of quartz CW-OSL signal components with R, *in preparation*.
#'
#' @references
#'
#' Mittelstraß, D., 2019. Decomposition of weak optically stimulated luminescence signals and its application in retrospective dosimetry at quartz (Master thesis). TU Dresden, Dresden.
#'
#'
#' @seealso [fit_OSLcurve], [decompose_OSLcurve], [plot_OSLcurve]
#'
#' @examples
#'
#' # Set some reasonable parameter for a weak quartz CW-OSL decay
#' components <- data.frame(name = c("fast", "medium", "slow"), lambda = c(1.5, 0.5, 0.1), n = c(1000, 1000, 10000))
#'
#' # Simulate the CW-OSL curve and add some signal noise
#' curve <- simulate_OSLcomponents(components, simulate.curve = TRUE, add.poisson.noise = TRUE)
#'
#' # Display the simulated curve
#' plot_OSLcurve(curve, components)
#'
#' @md
#' @export
simulate_OSLcomponents <- function(components,
                              curve = NULL,
                              channel.width = 0.1,
                              channel.number = 400,
                              simulate.curve = FALSE,
                              add.poisson.noise = TRUE,
                              add.gaussian.noise = 0,
                              add.background = 0,
                              round.values = TRUE){

  # Changelog:
  # * 2019-03-04, DM: combined old code fragments from 2014 and created this function
  # * 2019-03-06, DM: works now well with [plot_OSLcurve]
  # * 2019-10-02, DM: added background component
  # * 2020-07-24, DM: value rounding when curve is simulated is now optional; little tweaks to increase performance
  # * 2020-10-30, DM: Renamed from simulate_OSLcurve to simulate_OSLcomponents; Renamed argument from template.curve to curve; Rewrote roxygen documentation
  #
  # ToDo:
  # * Check literature for exact noise model and put references into roxygen docu
  # * Add not-first-order case
  # * Optimize computing time
  # * Use lambda.error to actually vary lambda when simulating curves
  # * Add input checks

  signal <- 0
  K <- nrow(components)

  # Which component curves shall be calculated? Can be a single number or a selection. If NA, all are added together.
  component.selection <- NA
  X <- component.selection
  if (is.na(X)) X <- 1:K

  # is there a template for the x-axis provided?
  if (!is.null(curve)) {

    time <- curve$time
    channel.width <- time[2] - time[1]

    if (!("signal" %in% colnames(curve))) {

      simulate.curve <- TRUE

    }else{

      signal <- curve$signal}

  } else  {

    simulate.curve <- TRUE

    # build x-axis
    time <- c(1:channel.number) * channel.width}

  # Create another time vector, which includes zero. You will see why
  #time0 <- c(0, time[1:(length(time)-1)])
  time0 <- c(0, time)

  # create a data table
  data <- data.frame(time = time, signal = signal, sum = 0, residual = 0)




  ##============================================================================##
  ## Building the components
  ##============================================================================##

  # add components
  for (i in X) {

    lambda <- components[i,]$lambda
    n <- components[i,]$n
    component <- data.frame(A = time)

    if(is.na(lambda)) {
      if (is.na(n)) {

        component$A <- rep(0, length(time))
      } else {

        component$A <- rep(n * channel.width, length(time))}

    } else {

      # Speed up things here by calculating just "exp(-lambda*time)" and subtracting
      # it from
      channel_probs <- exp(-lambda * time0)
      component$A <- - n * diff(channel_probs)
      #component$A <- n * (exp(-lambda*(time - channel.width)) - exp(-lambda*time))
    }


    data$sum <- data$sum + component$A

    colnames(component) <- components[i,]$name
    data <- cbind(data, component)}



  ##============================================================================##
  ## Curve simulation
  ##============================================================================##

  if (simulate.curve) {

    signal <- data$sum

    if (add.poisson.noise == TRUE){
      for (i in c(1:length(signal))) signal[i] <- rpois(1, signal[i])}

    if ((add.gaussian.noise > 0) || (add.background != 0)) {

      stddev <- sqrt(channel.width) * add.gaussian.noise
      offset <- add.background * channel.width

      signal <- signal + rnorm(length(signal), mean = offset, sd = stddev)}

    if (round.values) signal <- round(signal)

    data$signal <- signal

  }
  #+++ non-first order kinetics +++

  # decay formula:  f(t) = A * (1 + B * t)^C
  # therefore:      F(t) = (A  / ((C + 1) * B) * (1 + B * t)^(C + 1)
  #
  # parameters from Yukihara & McKeever (2011):
  #   A <- (n0^b * f) / (N^(b - 1))
  #   B <- (b - 1) * (n0 / N)^(b - 1) * f
  #   C <- -b / (b - 1)

  # for(i in c(1:channels)){
  #   I1 <- -A / ((C + 1) * B) * (1 + B * (i - 1)*channel.time)^(C + 1)
  #   I2 <- -A  / ((C + 1) * B) * (1 + B * i*channel.time)^(C + 1)
  #   signal[i] <- I1 - I2}



  ### Build residual curve
  data$residual <- data$signal - data$sum

  return(data)
}
