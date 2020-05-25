#' Alternativy methods for OSL decomposition
#'
#' @param curve
#' @param components
#' @param algorithm
#' @param offset.value
#' @param negative.values.to.zero
#'
#' @return
#' @export
#'
#' @examples
decompose_OSLalternatively <- function(
  curve,
  components,
  algorithm = "nls",
  offset.value = 0,
  negative.values.to.zero = FALSE
){
  # Last changed. 2019-09-11

  ## Work-Log Dirk
  #   2019-02-14 Started coding
  #   2019-03-15 Seperated from fit_OSLcurve
  #   2019-05-07 Guess starting n values from 'calc_OSLintervals()' results
  #   2019-09-11 Fixed bug: 1-component case gave 'n.error = NA'
  #
  ## TODO
  #   Add background value handling
  #   Add Shen & lang approach
  #   Documentation

  ### INTEGRITY CHECKS ###

  if(is(curve, "RLum.Data.Curve") == FALSE & is(curve, "data.frame") == FALSE & is(curve, "matrix") == FALSE){
    stop("[decompose_OSLalternatively()] Error: Input object is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")
  }

  if(is(curve, "RLum.Data.Curve") == TRUE) curve <- as.data.frame(get_RLum(curve))

  if (!("time" %in% colnames(curve)) ||
      !("signal" %in% colnames(curve))) {
    curve <- data.frame(time = curve[,1],
                        signal = curve[,2])
  }

  curve$signal <- curve$signal - offset.value
  channel.width <- curve$time[2] - curve$time[1]

  # check if time beginns with zero and add channel.width if the case
  if (curve$time[1] == 0)  curve$time <- curve$time + channel.width


  component.number <- length(components[,1])
  X <- 1:component.number
  n <- NA[X]

  ## START PARAMETER ESTIMATION
  if (is.null(components$ch.start) || any(is.na(components$ch.start)) ||
      is.null(components$ch.end) || any(is.na(components$ch.end))) {

    # if no preinformation is given, divide the curve in n areas of the same size and guess the area values as start n
    n <- rep(sum(curve) / component.number, component.number)

  } else {

    # take the signal integrals, where a component is dominant as first guess
    for (i in X) {

      n[i] <-  sum(curve$signal[components$ch.start[i]:components$ch.end[i]])
    }
  }


  ##============================================================================##
  ## FITTING
  ##============================================================================##


 if (algorithm == "nls") {

    ################### NLS standard with fixed lambdas ###########################


    # Create formula
    n.names <- paste0("n.",1:component.number)
    #decays <- paste(n.names," * ", components$lambda, "* exp(-",components$lambda," * time)", collapse=" + ")
    decays <- paste(n.names," * (exp(-",components$lambda," * (time - ", channel.width,")) - exp(-",components$lambda," * time))"
                    , collapse=" + ")
    fit.formula <- as.formula(paste0("signal ~ ", decays))

    names(n) <- n.names

    fit <- try(nls(fit.formula,
                   data = curve,
                   start = c(n)),
               silent = FALSE)

    if (attr(fit,"class") == "try-error") {

      components$n <- NA
      components$n.error <- 0
    } else {

      components$n <- coef(fit)
      #components$n.error <- summary(fit)$parameters[, "Std. Error"][paste0("n.",1:component.number)]
      components$n.error <- summary(fit)$parameters[, "Std. Error"][X]
    }

#return(fit)



  } else if (algorithm == "lm") {

    ################### LM with fixed lambdas ###########################
  }


  stim.end <- curve$time[length(curve$time)]
  for (i in 1:component.number) {
    if (!is.na(components$n[i])) {

      if (negative.values.to.zero && (components$n[i] < 0)) components$n[i] <- 0
      components$n.residual[i] <- round(components$n[i] * exp(- stim.end * components$lambda[i]))
    }

  }

  return(components)
}


