#' Classic CW-OSL signal calculation by late light or early light background substraction
#'
#' The function calculates the CW-OSL signal by late light or early light background substraction.
#' The intervals are set by the rules formulated in Galbraith and Roberts (2012) and Cunningham and Wallinga (2010).
#' Error calculation is performed by using formula (3) in Galbraith and Roberts (2012), also used in
#' Analyst v3.24 see Duller (2007) and as "poisson" approach in [calc_OSLLxTxRatio]
#'
#' Late light background approach:
#' * Signal: First 1 sec
#' * Background: Last 5 sec
#' * Round up if intervals are not divisible through channel time
#'
#' Early light background approach:
#' * Signal: First 0.4 sec
#' * Background: Following 1 sec (resp. 2.5x signal time)
#' * Round up if intervals are not divisible through channel time
#'
#' @param values [RLum.Data.Curve-class] or [data.frame] (**required**):
#' CW-OSL record. First row (x-axis) must contain time marks (seconds), Second row (y-axis) must contain signal values
#'
#' @param type [character] (*with default*): "late" for late light background substraction, "early" for early light background substraction
#'
#'
#' @return
#' A [data.frame] object is returned containing the following elements:
#' \item{signal}{[numeric] background corrected signal value}
#' \item{noise}{[numeric] standard deviation of the background corrected signal}
#' \item{signal.raw}{[numeric] signal value before background substraction}
#' \item{signal.min}{[numeric] start channel signal integral}
#' \item{signal.max}{[numeric] end channel signal integral}
#' \item{signal.duration}{[numeric] time interval signal (sec)}
#' \item{background.raw}{[numeric] background integral without normalisation}
#' \item{background.min}{[numeric] start channel background integral}
#' \item{background.max}{[numeric] end channel background integral}
#' \item{background.duration}{[numeric] time interval background (sec)}
#'
#' @section Changelog:
#' * 2018-08-08, DM: first version
#'
#' @section ToDo:
#' * Add non-poisson error calculation
#'
#' @section Function version: 0.0.1
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [calc_deconvolution.intervals], [analyse_SAR.CWOSL.deconvolved], [calc_superposition.curve]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for dose evaluation within the R luminescence package
#'
#' @keywords CWOSL deconvolution components
#'
#' @md
#' @export
#'
#' @examples
#'
calc_classic.OSL.signal <- function(
  values,
  type = "late"
){

  # INTEGRITY CHECKS -------------------------------------------------------

  ##INPUT OBJECTS
  if(is(values, "RLum.Data.Curve") == FALSE & is(values, "data.frame") == FALSE & is(values, "matrix") == FALSE){
    stop("[calc_CWOSLcomponents] Error: Input object is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")
  }


  if(is(values, "RLum.Data.Curve") == TRUE) values <- as.data.frame(get_RLum(values))

  x <- values[,1]
  y <- values[,2]

  # Channel width
  dt <- x[2] - x[1]
  # we will use 'dt' to build a x-axis!
  # Why? Is the first x-Value in the provided data = 0.0 sec (start of counting)
  # or is it 0 + dt sec (end of counting)?




  # preset result table, so the data format is static, even if not all cells are filled
  ooo <- data.frame(signal = 0,
                    noise = 0,
                    signal.raw = NA,
                    signal.min = NA,
                    signal.max = NA,
                    signal.duration = NA,
                    background.raw = NA,
                    background.min = NA,
                    background.max = NA,
                    background.duration = NA)

  if (type == "late") {
    ooo$signal.min <- 1
    ooo$signal.max <- ceiling(1 / dt)
    ooo$signal.duration <- ooo$signal.max * dt
    ooo$signal.raw = sum(y[ooo$signal.min:ooo$signal.max])

    ooo$background.max <- length(y)
    ooo$background.min <- length(y) - ceiling(5/dt)
    if (ooo$background.min <= ooo$signal.max) ooo$background.min <- ooo$signal.max + 1
    if (ooo$background.min > ooo$background.max) stop("data curve too short")
    ooo$background.duration <- (ooo$background.max - ooo$background.min + 1) * dt
    ooo$background.raw <- sum(y[ooo$background.min:ooo$background.max])


  } else {
    # "early"
    ooo$signal.min <- 1
    ooo$signal.max <- ceiling(0.4 / dt)
    ooo$signal.duration <- ooo$signal.max * dt
    ooo$signal.raw = sum(y[ooo$signal.min:ooo$signal.max])

    ooo$background.min <- ooo$signal.max + 1
    ooo$background.max <- ooo$signal.max + ceiling(1 / dt)

    if (ooo$background.max > length(y)) ooo$background.max <- length(y)
    if (ooo$background.min  > ooo$background.max) stop("data curve too short")
    ooo$background.duration <- (ooo$background.max - ooo$background.min + 1) * dt
    ooo$background.raw <- sum(y[ooo$background.min:ooo$background.max])
  }

  # poisson noise (absolute)
  k <-  ooo$background.duration / ooo$signal.duration
  ooo$signal <- ooo$signal.raw - ooo$background.raw / k
  ooo$noise <- (ooo$signal.raw + ooo$background.raw / k^2)^0.5


    # extra-poisson noise
    # b.mean <- background / length(ooo$background.min:ooo$background.max)
    #background.noise <- var(y(ooo$background.min:ooo$background.max))


  return(ooo)
}
