#' Builds a CW-OSL curve from the parameters of its components or adds component columns to an existing curve
#'
#' @param components [data.frame] (**required**):
#' Table with the component parameters. It needs the rows "lambda" and "n" containing numeric values
#'
#' @param template.curve [data.frame] (*with default*):
#' OSL curve providing its x-axis (time axis) as template. The signal values of this curve are NOT used.
#'
#' @param channel.width [numeric] (*with default*):
#' channel width in seconds
#'
#' @param channel.number [integer] (*with default*):
#' number of channels resp. data points
#'
#' @param component.selection [vector] (*with default*):
#' Which component curves shall be calculated? Can be a single number or a selection. If NA, all are added together.
#'
#' @param add.poisson.noise [logical] (*with default*):
#' adds poisson distributed noise to input signal
#'
#' @param add.gaussian.noise [numeric] (*with default*):
#' standard deviation of the detection background in **cts/s**
#'
#' @param add.background [numeric] (*with default*):
#' background caused signal offset in **cts/s**
#'
#' @return
#'
#' **CW-OSL curve**
#' A [data.frame] with the X = time and Y = signal
#'
#' @section Changelog:
#' * 2019-03-04, DM: combining old functions (simulate_...) from post-bachelor-time (2014) to this one
#' * 2019-03-06, DM: works now well with plot_OSLcurve
#' * 2019-10-02, DM: added background component
#' * 2020-07-24, DM: value rounding when curve simulated now optional; little tweaks to increase performance
#'
#' @section ToDo:
#' * Check literature for exact noise model
#' * Add not-first-order case
#' * correct and expand Roxygen help
#' * Check handling of background component
#'
#' @section Last changed: 2020-07-24
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @export
#'
#' @examples
#' components <- data.frame(name = c("fast","medium","slow"), lambda = c(1.5,0.5,0.1), n = c(1000,1000,10000))
#' curve <- simulate_OSLcurve(components, simulate.curve = TRUE, add.poisson.noise = TRUE, add.background = 20)
#' plot(curve)
#'
simulate_OSLcurve <- function(components,
                              template.curve = NULL,
                              channel.width = 0.1,
                              channel.number = 400,
                              component.selection = NA,
                              simulate.curve = FALSE,
                              round.values = TRUE,
                              add.poisson.noise = TRUE,
                              add.gaussian.noise = 0,
                              add.background = 0){

  ##============================================================================##
  ## Preparations
  ##============================================================================##

  signal <- 0
  K <- nrow(components)

  # is there a template for the x-axis provided?
  if (!is.null(template.curve)) {

    time <- template.curve$time
    channel.width <- time[2] - time[1]

    signal <- template.curve$signal

  } else  {

    # build x-axis
    time <-  c(1:channel.number)*channel.width
  }

  # Create another time vector, which includes zero. You will see why
  #time0 <- c(0, time[1:(length(time)-1)])
  time0 <- c(0, time)


  # create a data table
  data <- data.frame(time = time, signal = signal, sum = 0, residual = 0)

  X <- component.selection
  if (is.na(X)) X <- 1:K

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

        component$A <- rep(n * channel.width, length(time))
      }


    } else {

      # Speed up things here by calculating just "exp(-lambda*time)" and subtracting
      # it from
      channel_probs <- exp(-lambda * time0)
      component$A <- - n * diff(channel_probs)
      #component$A <- n * (exp(-lambda*(time - channel.width)) - exp(-lambda*time))
    }


    data$sum <- data$sum + component$A

    colnames(component) <- components[i,]$name
    data <- cbind(data, component)

  }

  #+++ non-first order kinetics +++

    # decay formula:  f(t) = A * (1 + B * t)^C
    # therefore:      F(t) = (A  / ((C + 1) * B) * (1 + B * t)^(C + 1)
    #
    # parameters from Yukihara & McKeever (2011):
 #   A <- (n0^b * f) / (N^(b - 1))
 #   B <- (b - 1) * (n0 / N)^(b - 1) * f
 #   C <- -b / (b - 1)

  #  for(i in c(1:channels))
  #  {
  #    I1 <- -A / ((C + 1) * B) * (1 + B * (i - 1)*channel.time)^(C + 1)
   #   I2 <- -A  / ((C + 1) * B) * (1 + B * i*channel.time)^(C + 1)
   #   signal[i] <- I1 - I2}

  ##============================================================================##
  ## Curve simulation
  ##============================================================================##

  if (simulate.curve) {

    signal <- data$sum

    if (add.poisson.noise == TRUE){
      for (i in c(1:length(signal))){
        signal[i] <- rpois(1, signal[i])
      }
    }

    if ((add.gaussian.noise > 0) || (add.background != 0)) {

      stddev <- sqrt(channel.width) * add.gaussian.noise
      offset <- add.background * channel.width

      signal <- signal + rnorm(length(signal), mean = offset, sd = stddev)
    }


    if (round.values) {
      signal <- round(signal)
    }

    data$signal <- signal
  }

  ##============================================================================##

  ### Build residual curve
  data$residual <- data$signal - data$sum

  return(data)

}
