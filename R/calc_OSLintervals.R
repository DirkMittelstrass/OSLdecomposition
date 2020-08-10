#' Calculate integration intervals for CW-OSL deconvolution
#'
#' The function provides the integration intervals for CW-OSL component separation with [decompose_OSLcomponents]
#'
#' @param components [data.frame] (**required**):
#' Table containing the decay constants of the signal components
#'
#' @param curve [data.frame] (*optional*):
#' OSL signal curve. The x-axis (time axis) will be used ot define channel.width and channel.number
#'
#' @param channel.width [numeric] (*optional*):
#' channel width in seconds. Necessary if *curve* is not given
#'
#' @param channel.number [integer] (*optional*):
#' number of channels resp. data points. Necessary if *curve* is not given
#'
#' @param t.start [numeric] (*with default*):
#' starting point of the first interval, per default the start of the measurement
#'
#' @param t.end [numeric] (*optional*):
#' end point of the last interval, per default the end of the measurement
#'
#' @return
#' The input table *components* [data.frame] will be returned with four additional columns:
#' *$t.start*, *$t.end* defining the interval borders in time; *$ch.start*, *$ch.end* defining the intervals as channels
#'
#'
#' @section Changelog:
#' * 2018-04-05, DM: first running version
#' * 2018-06-16, DM: added 2-component and 1-component case
#' * 2018-06-19, DM: changed t0 and t3 parameter to t.start and t.end and added full increment proof
#' * 2018-06-24, DM: changed data structure to get static tables for [analyse_SAR.CWOSL.deconvolved]
#' * 2019-03-21, DM: Rewritten for arbitrary component numbers and changed data structure
#' * 2019-10-02, DM: added Background-fitting interval determination
#'
#' @section ToDo:
#' * rename function into "optimise_OSLintervals" or similar
#' * replace own algorithm with DEoptim()
#' * optional RLum.Data.Curve as input template
#' * add more comments
#'
#' @section Last changed: 2019-10-02
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [decompose_OSLcomponents], [...], [..]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for dose evaluation within the R luminescence package
#'
#' @keywords CWOSL OSL decomposition deconvolution components
#'
#' @examples
#' ...
#'
#' @md
#' @export


calc_OSLintervals <- function(
  components,
  curve = NULL,
  background.fitting = TRUE,
  channel.width = NA,
  channel.number = NA,
  t.start = 0,
  t.end = NA,
  verbose = TRUE
){

  ########## is a template curve given? ###########

  if (!is.null(curve)) {

    dt <- curve$time[2] - curve$time[1]
    n <- length(curve$time)

  } else if ((!is.na(channel.width)) && (!is.na(channel.number))) {

    dt <- channel.width
    n <- channel.number

  } else {
    warning("No template curve nor channel parameters given")
    return(components)
  }

  # round start and end point to full increments
  t0 <- floor(t.start / dt) * dt

  ## set t.end if not preset
  if (is.na(t.end) || (t.end > n*dt) || (t.end < 3*dt)) {
    t.end <- n*dt
  } else {
    t.end <- floor(t.end / dt) * dt
  }

  # If the background also shall be fitted, add another component wit lambda = 0
  if (background.fitting) {


    if (!is.na(components$lambda[nrow(components)])) {
      new.row <- components[nrow(components),]
      rownames(new.row) <- "Background"
      new.row[1:length(new.row)] <- NA
      new.row[1]  <- "Background"
      components <- rbind(components, new.row)
    }
  } else {

    if (is.na(components$lambda[nrow(components)])) {
      components <- components[1:(nrow(components) - 1),]
    }
  }

  component.number <- nrow(components)
  lambdas <- components$lambda


  ########## is there just 1 components? ###########

  if (component.number == 1) {

    components$t.start <- t0
    components$t.end <- t.end
    components$ch.start <- 1 + floor(t0 /  dt)
    components$ch.end <- ceiling(t.end / dt)

    return(components)
  }

  ########## Create matrix ###########

  ## channels: just the end of an interval. For example: 1:3 / 4:5 / 6:7 becomes c(0,3,5,7)
  calc_determinant <- function(lambdas, channels,component.number, dt) {

    M <- matrix(0,component.number,component.number)
    for (i in c(1:component.number)) {

      for (j in c(1:component.number)) {

        #P <- exp(-t0 * f.fast) - exp(-t1 * f.fast)
        if (is.na(lambdas[j])) {

          M[i, j] <- channels[i + 1] * dt - channels[i] * dt

        } else {

          M[i, j] <- exp(- channels[i] * dt * lambdas[j]) - exp(- channels[i + 1] * dt * lambdas[j])
        }

      }
    }

    # calc determinante
    return(abs(det(M)))
  }


  ########## Search for determinant maximum ###########

  max_found <- FALSE
  max_iteration <- 50
  iterations <- 0
  iterations.stop <- 20 * 2^component.number
  iterations.bests <- 2 * 2^component.number
  min.ch <- NULL
  max.ch <- NULL
  A <- data.frame(NULL)

  # start parameters for the choosable time intervals
  for (i in c(1:component.number - 1)) {
    min.ch <- c(min.ch, i)
    max.ch <- c(max.ch, n - component.number + i)
  }

  while ((max_found == FALSE) && (iterations < max_iteration)) {

    # create random set of numbers
    interval_set <- c(1,1)
    while (any(duplicated(interval_set))) {
      interval_set <- NULL
      for (i in c(1:(component.number - 1))) {
        interval_set <- c(interval_set, round(runif(1, min = min.ch[i], max = max.ch[i])))
      }
    # interval_set <- round(runif(component.number - 1, min = t0 / dt + 1, max = t.end / dt - 1))
    }

    # ... and sort it
    interval_set <- interval_set[order(interval_set)]

    # add start and end value
    #interval_set <- c(t0 / dt, interval_set, t.end / dt)

    # now calc the determinante
    det.value <- calc_determinant(lambdas,
                                  c(t0 / dt, interval_set, t.end / dt),
                                  component.number,
                                  dt)

    # append the parameters to a data.frame
    A <- rbind(A, c(det.value, interval_set))

    if (nrow(A) == iterations.stop) {

      # order list elements by determinant value
      A <- A[order(A[,1]),]

      # reverse order to get large values first and reduce to the best 10 rows
      A <- A[nrow(A):(iterations.stop - iterations.bests),]

      # check if all elements are the same, then
      if (all(A[,1] == A[1,1])) {

        interval_set <- A[1,2:component.number]
        max_found <- TRUE

        if (verbose) {
          writeLines(paste0("Maximum determinant = ", round(A[1,1], digits = component.number + 2),
                            " with interval breaking channels [", paste0(interval_set, collapse = ", " ),
                            "] found after ", iterations * iterations.stop, " iterations"))
        }

      } else {

        # redefine start parameters
        for (i in c(1:(component.number - 1))) {
          min.ch[i] <- min(A[,i + 1])
          max.ch[i] <- max(A[,i + 1])
        }
        iterations <- iterations + 1
        A <- data.frame(NULL)
      }

    }

  }

'  components <- cbind(components,
                      data.frame(t.start = c(t0, interval_set * dt),
                                 t.end = c(interval_set * dt, t.end),
                                 ch.start = c(1 + floor(t0 /  dt), interval_set + 1),
                                 ch.end = c(interval_set, ceiling(t.end / dt))))'

  components$t.start <- unlist(c(t0, interval_set * dt))
  components$t.end <- unlist(c(interval_set * dt, t.end))
  components$ch.start <- unlist(c(1 + floor(t0 /  dt), interval_set + 1))
  components$ch.end <- unlist(c(interval_set, ceiling(t.end / dt)))

  return(components)

}
