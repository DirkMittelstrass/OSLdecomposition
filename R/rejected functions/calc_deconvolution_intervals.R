#' Calculate integration intervals for CW-OSL deconvolution
#'
#' The function provides the integration intervals for CW-OSL component separation with [calc_CWOSLcomponents]
#'
#' @param f.fast [numeric] (**required**):
#' fast component decay parameter
#'
#' @param f.medium [numeric] (*optional*):
#' medium component decay parameter
#'
#' @param f.slow [numeric] (*optional*):
#' slow component decay parameter
#'
#' @param channel.width [numeric]
#' channel width in seconds
#'
#' @param channel.number [integer]
#' number of channels resp. data points
#'
#' @param t.start [numeric] (*with default*):
#' starting point of the first interval, per default the start of the measurement
#'
#' @param t.end [numeric] (*optional*):
#' end point of the last interval, per default the end of the measurement
#'
#' @return
#' A [data.frame] object is returned containing the following elements:
#'
#' \item{t0 ... t3}{[numeric] The interval boundaries. Be aware: t1 = end of first interval and start of second, and so on}
#' \item{determinant.value}{[numeric] determinant value of the release probability matrix of the equation system and denominator in cramers rule; serves as optimization target; transfered to [calc_CWOSLcomponents] a bit of computing time can be saved}
#'
#' @section Changelog:
#' * 2018-04-05, DM: first running version
#' * 2018-06-16, DM: v0.2
#' * 2018-06-16, DM: added 2-component and 1-component case
#' * 2018-06-19, DM: changed t0 and t3 parameter to t.start and t.end and added full increment proof
#' * 2018-06-24, DM: changed data structure to get static tables for [analyse_SAR.CWOSL.deconvolved]
#'
#' @section ToDo:
#' * optional data.frame with result-matrix
#' * optional RLum.Data.Curve as input and error calculation as result-matrix
#' * more comments
#'
#' @section Function version: 0.2.2
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [calc_CWOSLcomponents], [analyse_SAR.CWOSL.deconvolved], [calc_superposition.curve]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for dose evaluation within the R luminescence package
#'
#' @keywords CWOSL deconvolution components
#'
#' @examples
#' calc_deconvolution.intervals(1.5, 0.5, 0.1, channel.width = 0.1, channel.number = 100)
#'
#' @md
#' @export


calc_deconvolution.intervals <- function(
  f.fast,
  f.medium = NA,
  f.slow = NA,
  channel.width,
  channel.number,
  t.start = 0,
  t.end = NA
){

  ##============================================================================##
  # INPUT OBJECTS & INTEGRITY CHECKS
  ##============================================================================##

  # check if fast component is given
  if (is.na(f.fast)) {
    warning("Decay parameter for fast component missing")
    return(NULL)
  } else {
    if (!(f.fast > 0)) {
      warning("Decay parameter for fast component must be greater than zero")
      return(NULL)
    }
  }

  dt <- channel.width
  n <- channel.number

  ## decrease number of channels if necessary to increase performance
  while (n > 500) {
    n <- floor(n/2)
    dt <- dt*2
  }


  # round start and end point to full increments
  t0 <- floor(t.start / dt) * dt

  ## set t.end if not preset
  if (is.na(t.end) || (t.end > n*dt) || (t.end < 3*dt)) {
    t.end <- n*dt
  } else {
    t.end <- floor(t.end / dt) * dt
  }

  results <- data.frame(t0 = t0,
                        t1 = NA,
                        t2 = NA,
                        t3 = NA,
                        determinant.value = 0)

  ##============================================================================##
  # Check if there are less components than 3
  ##============================================================================##

  ### no f.medium and f.slow are set (1-component case):
  if (is.na(f.medium) && is.na(f.slow)) {
    message("Just fast component given, therefore just one interval calculated (tO,t1)")

    t1 <- t.end
    results$t1 <- t1
    P1.fast <- (exp(-t0 * f.fast) - exp(-t1 * f.fast))
    results$determinant.value <- P1.fast
  }

  ### f.medium OR f.slow is missing, but not both:
  if ((is.na(f.medium) || is.na(f.slow)) && !(is.na(f.medium) && is.na(f.slow))) {

    ##============================================================================##
    #   2 - COMPONENT CASE
    ##============================================================================##

    message("Just two components given, therefore just two intervals calculated (tO,t1,t2)")

    if (is.na(f.slow)) {
      f.slow <- f.medium
      f.medium <- NA
    }

    t2 <- t.end
    results$t2 <- t2

    ## precalculate values related to t0 and t2
    t0.fast <- exp(-t0 * f.fast)
    t0.slow <- exp(-t0 * f.slow)

    t2.fast <- exp(-t2 * f.fast)
    t2.slow <- exp(-t2 * f.slow)

    # vary t1
    for (i in 2:(n-1)) {

      t1 <- i*dt
      t1.fast <- exp(-t1 * f.fast)
      t1.slow <- exp(-t1 * f.slow)

      # calculate the release propabilities during an interval
      P1.fast <- (t0.fast - t1.fast)
      P1.slow <- (t0.slow - t1.slow)

      P2.fast <- (t1.fast - t2.fast)
      P2.slow <- (t1.slow - t2.slow)

      # create the release propability matrix ..
      A <- matrix(c(P1.fast, P1.slow, P2.fast, P2.slow), nrow = 2)

      # ... and calculate the determinant
      value <- abs(det(A))

      if (value > results$determinant.value) {
        results$determinant.value <- value
        results$t1 <- t1
      }
    }
  }

  if (!is.na(f.medium) && !is.na(f.slow)) {
    ##============================================================================##
    #   3 - COMPONENT CASE
    ##============================================================================##

    t3 <- t.end
    results$t3 <- t3

    ## precalculate values related to t0 and t3
    t0.fast <- exp(-t0 * f.fast)
    t0.medium <- exp(-t0 * f.medium)
    t0.slow <- exp(-t0 * f.slow)

    t3.fast <- exp(-t3 * f.fast)
    t3.medium <- exp(-t3 * f.medium)
    t3.slow <- exp(-t3 * f.slow)

    # vary t2
    for (i in 2:(n-1)) {

      t2 <- i*dt
      t2.fast <- exp(-t2 * f.fast)
      t2.medium <- exp(-t2 * f.medium)
      t2.slow <- exp(-t2 * f.slow)

      # vary t1
      for (j in 1:(i-1)) {

        t1 <- j*dt
        t1.fast <- exp(-t1 * f.fast)
        t1.medium <- exp(-t1 * f.medium)
        t1.slow <- exp(-t1 * f.slow)

        #formula (7-3) bachelor thesis Mittelstrass 2013
        P1.fast <- (t0.fast - t1.fast)
        P1.medium <- (t0.medium - t1.medium)
        P1.slow <- (t0.slow - t1.slow)

        P2.fast <- (t1.fast - t2.fast)
        P2.medium <- (t1.medium - t2.medium)
        P2.slow <- (t1.slow - t2.slow)

        P3.fast <- (t2.fast  - t3.fast)
        P3.medium <- (t2.medium - t3.medium)
        P3.slow <- (t2.slow - t3.slow)

        # create the matrix ..
        A <- matrix(c(P1.fast, P1.medium, P1.slow, P2.fast, P2.medium, P2.slow, P3.fast, P3.medium, P3.slow), nrow = 3)

        # ... and calculate the determinant
        value <- abs(det(A))

        if (value > results$determinant.value) {
          results$determinant.value <- value
          results$t1 <- t1
          results$t2 <- t2
        }
      }
    }
  }

  return(results)


  # DOCUMENTATION - INLINEDOC LINES -----------------------------------------

  ##details<<
  ## use this section for details on your function

  ##value<<
  ## \item{plot}{(optional)

  ##references<<
  ##


  ##note<<
  ## use this section for further notes

  ##seealso<<
  ## \code{\link{fit_LMCurve}}, \code{\link{plot}},
  ## \code{\linkS4class{RLum.Data.Curve}}, \code{\linkS4class{RLum.Results}},
  ## \code{\link{get_RLum.Results}}

  ##keyword<<
  ## dplot
  ## models

}
