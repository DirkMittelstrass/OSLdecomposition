#' classic CW-OSL signal calculation
#'
#' Early background calculation after [1], error calculation after (3) in [2].
#' If the stimulation intensity differs from 35 mW /cm^-2 [1], the intervals are increased/decreased accordingly
#'
#' [1]A. C. Cunningham and J. Wallinga, ‘Selection of integration time intervals for quartz OSL decay curves’, Quaternary Geochronology, vol. 5, no. 6, pp. 657–666, Dec. 2010.
#' [2]R. F. Galbraith and R. G. Roberts, ‘Statistical aspects of equivalent dose and error calculation and display in OSL dating: An overview and some recommendations’, Quaternary Geochronology, vol. 11, pp. 1–27, Aug. 2012.
#'
#'
#' @return
#' The input table **components** [data.frame] will be returned with added/overwritten columns
#'
#' @section Changelog:
#' *
#
#' @section ToDo:
#' * Update documentation
#' *
#'
#' @section Last changed: 2019-10-11
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [decompose_OSLcurve], [simulate_OSLcurve]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for natural dose evaluation
#'
#' @keywords
#'
#' @md
#' @export
#'
#' @examples
#'
calc_classicOSLsignal <- function(
  curve,
  components = NULL,
  algorithm = "late", # "early"
  signal.end.late = 0.8,
  signal.end.early = 0.4,
  stimulation.intensity = 35,
  verbose = TRUE
){

  ########## Input checks ###########

  if(is(curve, "RLum.Data.Curve") == FALSE & is(curve, "data.frame") == FALSE & is(curve, "matrix") == FALSE){
    stop("[decompose_OSLcurve()] Error: Input object is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")
  }

  if(is(curve, "RLum.Data.Curve") == TRUE) curve <- as.data.frame(get_RLum(curve))

  if (!("time" %in% colnames(curve)) ||
      !("signal" %in% colnames(curve))) {
    curve <- data.frame(time = curve[,1],
                        signal = curve[,2])
  }

  channel.width <- curve$time[2] - curve$time[1]

  stretcher <- 35 / stimulation.intensity
  signal.end.late <- stretcher * signal.end.late
  signal.end.early <- stretcher * signal.end.early

  # check if time beginns with zero and add channel.width if the case
  if (curve$time[1] == 0)  curve$time <- curve$time + channel.width


  ########## Set parameters ###########

  signal <- curve$signal
  time <- curve$time

  dt <- time[2] - time[1]
  channel.number <- length(time)

  components <- components

  if (is.null(components)) {
    components <- data.frame(name = c("Signal","Background"),
                             lambda = c(NA,NA))

    if (algorithm == "late") {
      ######## LATE BACKGROUND SUBSTRACTION ###############

      # where end the c(signal, background) intervals?
      ch.end <- c(ceiling(signal.end.late / dt),
                  channel.number)

      # and the matching time stamps?
      t.end <- c(ch.end[1] * dt,
                 channel.number * dt)

      # where start the c(signal, background) intervals?
      ch.start <- c(1,
                    channel.number - ch.end[1] * 10)

      # shorten background interval, if necessary
      if (ch.start[2] <= ch.end[1]) ch.start[2] <- ch.end[1] + 1

      # and the matching time stamps?
      t.start <- c(0,
                   ch.start[2] * dt)

    } else {
      ######## EARLY BACKGROUND SUBSTRACTION ###############

      # where start the c(signal, background) intervals?
      ch.start <- c(1,
                    ceiling(signal.end.early / dt) + 1)

      # and the matching time stamps?
      t.start <- c(0,
                   (ch.start[2] - 1) * dt)

      # where end the c(signal, background) intervals?
      ch.end <- c(ch.start[2] - 1,
                  ch.start[2] - 1 + ceiling(2.5*(ch.start[2] - 1)))

      # shorten background interval, if necessary
      if (ch.end[2] > channel.number) ch.end[2] <- channel.number

      # and the matching time stamps?
      t.end <- c(ch.end[1] * dt,
                 ch.end[2] * dt)

    }

    components$t.start <- t.start
    components$t.end <- t.end
    components$ch.start <- ch.start
    components$ch.end <- ch.end

  } else if ((nrow(components) != 2) || !((components$name[1] != "Signal") || (components$name[1] != "signal"))) {

    stop("Wrong input format of component table")
  } else {

    ch.start <- components$ch.start
    ch.end <- components$ch.end

  }

  # calc bins
  bin <- c(sum(signal[ch.start[1]:ch.end[1]]),
           sum(signal[ch.start[2]:ch.end[2]]))

  # calc signal
  offset <- bin[2] * (length(ch.start[1]:ch.end[1]) /
                      length(ch.start[2]:ch.end[2]))

  components$n <- c(bin[1] - offset, offset)

  # equation (3) from Galbraith (2012), also used in Galbraith (2002) and Duller (2007)
  components$n.error <- c(sqrt(bin[1] + offset), sqrt(offset))

  components$n.residual <- rep(0, 2)

  components$bin <- bin
  components$bin.error <- sqrt(bin)

  if (verbose) print(components)

  return(components)
}
