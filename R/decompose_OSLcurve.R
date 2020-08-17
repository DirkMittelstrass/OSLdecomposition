#' multi-exponential CW-OSL decomposition
#'
#' The function calculates the CW-OSL component amplitudes by a determinant-based algorithm.
#' It also estimates the standard deviation of the amplitudes by using the error propagation method.
#'
#' @param components [data.frame] or [numeric] vector (**required**):
#' Either a vector containing the decay parameters of the CW-OSL components or a table (data.frame), usually the table returned by [fit_OSLcurve].
#' In case of a vector: It is recommended to use less than 7 parameters. The parameters will be sorted in decreasing order.
#' In case of a data.frame. One column must be named *$lambda*.
#' It is recommended to provide also integration interval parameters (columns *t.start, t.end, ch.start, ch.end*),
#' which can be found by applying [calc_OSLintervals] on the global mean curve, calculated by [sum_OSLcurves].
#' If one or more column is missing, a simple interval definition algorithm is run automatically (see **Notes**)..
#'
#' @param curve [RLum.Data.Curve-class] or [data.frame] (**required**):
#' CW-OSL record. `$time` or first column (x-axis) for the measurement time (must have constand time intervals),
#' `$signal` or second column (y-axis) for the signal values. Further columns will be ignored
#'
#' @param error.calculation (*with default*):
#' integral error estimation approach, either **"empiric"** or **"poisson"** or [numeric] value;
#' Per default the data of *curve$residual* provided by [simulate_OSLcurve] is used to calculate an **empiric** standard error for each integral which will be processed in the error propagation formula.
#' Alternatively the integral standard error can be calculated by assuming a **poisson** distributed signal error, known as *Shot noise*.
#' This is suitable if the lack of data points on the x-axis circumvent an empiric error estimation, like with spatial or spectral resolved CCD measurements.
#' Also the parameter can be set to a [numeric] value; which will be handled as standard deviation per channel and added to the poisson distributed *Shot noise*.
#' If no error calculation is wished, set the argument to **"none"**. This reduces the necessary computing time heavily.
#'
#' @return
#' The input table **components** [data.frame] will be returned with added/overwritten columns: *$n, $n.error, $n.residual, $I, $I.error*
#'
#' @section Changelog:
#' * winter 2012/13: first version, used for Bachelor thesis DM
#' * autumn 2013   : added empiric error estimation, shown in germanLED Freiberg 2013
#' * 2014-11-??, SK: formated into Rluminecence package standard
#' * 2014-11-07, DM: Binomial error propagation added
#' * 2018-05-04, DM: added residuals for n values (necessary for slow component dosimetry) and many little tweaks
#' * 2018-06-22, DM: added decomposition of data sets with just 1 or 2 components
#' * 2018-06-22, DM: added negative.values.to.zero and set it on TRUE as default (on request of C. Schmidt) (later removed)
#' * 2018-07-05, DM: overworked error estimation; replaced binomial with poisson error approach; added auto-switch to poisson if integral length = 1; integrated background.noise into error
#' * 2019-03-29, DM: Rewritten function for several purposes: 1. working now with any number of components  2. shorter and more elegant code 3. data format more suitable for markdown/shiny applications.
#' * 2019-09-09, DM: Added anticorrelating covariance terms into error estimation (later removed)
#' * 2019-09-25, DM: Merged function with decompose_OSLalternatively() and added algorithm argument
#' * 2019-09-25, DM: Deleted unnecessary functions (negative.values.to.zero, offset, anticorrelation)
#' * 2019-10-02, DM: Added optional background fitting
#' * 2020-04-06, DM: Added 'initial.signal' column in output data.frame; cleaned print output
#' * 2020-07-20, DM: Added algorithm for fast interval definition based on logarithmic means; More input data checks
#
#' @section ToDo:
#' * Update documentation (example, notes)
#' * Replace Cramers rule equations with 'solve()' to increase performance
#' * In some very rare cases, negative values for n.error are returned. Why?
#' * Test and expand interval determination algorithm in case of very few (N ~ K) data points
#'
#' @section Last changed: 2020-08-17
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [calc_OSLintervals], [fit_OSLcurve], [sum_OSLcurves]
#'
#' @references
#'
#' @return
#' @export
#'
#' @examples
#'
decompose_OSLcurve <- function(
  curve,
  components,
  background.fitting = FALSE,
  algorithm = "det", # "det", "nls", "det+nls"
  error.calculation = "empiric",   # "poisson", "empiric", "nls", numeric value
  verbose = TRUE
){

  ########## Input checks ###########

  if(!(is(curve, "RLum.Data.Curve") | is(curve, "data.frame") | is(curve, "matrix"))){
   stop("[decompose_OSLcurve()] Error: Input object 'curve' is not of type 'RLum.Data.Curve' or 'data.frame' or 'matrix'!")
  }

  if(is(curve, "RLum.Data.Curve") == TRUE) curve <- as.data.frame(Luminescence::get_RLum(curve))

  if (!("time" %in% colnames(curve)) ||
      !("signal" %in% colnames(curve))) {
    curve <- data.frame(time = curve[,1],
                        signal = curve[,2])}

  if ((algorithm == "nls") &! (error.calculation == "nls")) {
    if (verbose) warning("When algorithm 'nls' is chosen, error.calculation must be also 'nls'. Argument changed to error.calculation='nls'")
    error.calculation <- "nls"}

  # What is the channel width?
  dt <- curve$time[2] - curve$time[1]

  # check if time beginns with zero and add dt if the case
  if (curve$time[1] == 0)  curve$time <- curve$time + dt

  # Check if 'components' is of valid type
  if (class(components)=="numeric") {
    components <- data.frame(names = paste0("Component ", 1:length(components)),
                             lambda = sort(components, decreasing = TRUE))

  }else if(class(components)=="data.frame"){
    if (!("lambda" %in% colnames(components))) {
      stop("[decompose_OSLcurve()] Error: Input object 'components' contains no column '$lambda'!")}

  }else{
    stop("[decompose_OSLcurve()] Error: Input object 'components' is not of type 'numeric vector' or 'data.frame' !")}

  # if background.fitting = FALSE (recommended), remove last row
  # this removes also the last integration interval (which is good)
  if (is.na(components$lambda[nrow(components)]) && (background.fitting==FALSE)) {
    components <- components[1:(nrow(components)-1),]}


  lambda <- components$lambda
  K <- nrow(components)
  X <- c(1:K)

  if (K > nrow(curve)) {
    stop("[decompose_OSLcurve()] Error: Number of decay rates in 'components' exceeds number of data points in 'curve'!")}


  # are the integration intervals given?
  if (!("t.start" %in% colnames(components)) ||
      !("t.end" %in% colnames(components)) ||
      !("ch.start" %in% colnames(components)) ||
      !("ch.end" %in% colnames(components))) {
    #if (verbose) warning("Integration intervals not provided. calc_OSLintervals() executed")

    # Define the K = 1 case first:
    ch.start <- 1
    ch.end <- nrow(curve)

    if (K > 1) {

      # Calc the logarithmic means between following lambdas
      intervals <- diff(log(lambda)) / diff(lambda)

      # Test if each interval starts before k/K
      intervals <- pmin(intervals, curve$time[nrow(curve)] * c(1:(K-1)) / K)

      # Round values up to full channels
      ch.end <- ceiling(intervals / dt)
      ch.start <- c(1, ch.end + 1)
      ch.end <- c(ch.end, nrow(curve))

      # Test if each interval is at least one channel wide
      for (i in 1:(K-1)) {
        if ((ch.end[i] - ch.start[i]) < 1) {
          ch.end[i] <- ch.start[i] + 1
          ch.start[i + 1] <- ch.end[i] + 1}}

      # In the very unlikely event that the last interval is shifted out of the measurement
      if (ch.start[K] >= ch.end[K]) {
        stop("[decompose_OSLcurve()] Error: Last interval is shifted out of the measurement.")}
      }

    t.start <- (ch.start - 1) * dt
    t.end <- ch.end * dt
    #if (verbose) cat("Intervals set to: ")

    components$t.start <- t.start
    components$t.end <- t.end
    components$ch.start <- ch.start
    components$ch.end <- ch.end

  } else {

    t.start <- components$t.start
    t.end <- components$t.end
    ch.start <- components$ch.start
    ch.end <- components$ch.end

  }




  ########## Set parameters ###########

  signal <- curve$signal[1:components$ch.end[K]]
  time <- curve$time[1:components$ch.end[K]]

  components$n <- rep(NA, K)
  #components$n.error <- rep(NA, K)
  components$n.residual <- rep(NA, K)


  ### calculate integrals  ###

  I <- NULL
  for (i in X) {
    I <- c(I, sum(signal[c(ch.start[i]:ch.end[i])]))
  }
  components$bin <- I
  #components$bin.error <- rep(NA, K)

  n <- NULL

  ######################### DET ###########################

  if ((algorithm == "det")||(algorithm == "det+nls")) {


    # ToDo: Substitute Cramers rule with solve()
    ### define matrices ###

    # Build denominator matrix
      D <- matrix(0, K, K)
      for (i in X) {
        for (j in X) {

          if (is.na(lambda[j])) {

            D[i, j] <- t.end[i] - t.start[i]
          } else {

            D[i, j] <- exp(-t.start[i] * lambda[j]) - exp(- t.end[i] * lambda[j])
          }


        }
      }

    # Build enumerator matrices
    A <- list(NULL)
    for (j in X) {

      A.temp <- D
      A.temp[,j] <- I
      A[[j]] <- A.temp
    }

    ### Calculate component amplitudes ###
    for (i in X) {

      n.temp <- det(A[[i]])/det(D)
      n <- c(n, n.temp)
    }

    components$n <- n

  }  ########### end DET ############

  ######################### NLS ###########################

  if ((algorithm == "nls")||(algorithm == "det+nls")) {

    # use outcome from DET as start parameters. If not given, use integral values
    if(is.null(n)) n <- I

    ### Create fit formula ###
    n.names <- paste0("n.",1:K)

    if (is.na(components$lambda[K])) {

      lambda <- components$lambda[1:(K - 1)]
      decays <- paste(n.names[1:(K - 1)],
                      " * (exp(-",lambda," * (time - ", dt,")) - exp(-",lambda," * time))"
                    , collapse=" + ")
      decays <- paste0(decays, " + ", n.names[K], " * ",dt)


    } else {

      decays <- paste(n.names," * (exp(-",components$lambda," * (time - ", dt,")) - exp(-",components$lambda," * time))"
                      , collapse=" + ")
    }



    fit.formula <- as.formula(paste0("signal ~ ", decays))

    names(n) <- n.names

    ### try Gauss-Newton fit ###
    fit <- try(nls(fit.formula,
                   data = curve,
                   start = c(n)),
               silent = TRUE)

    if (attr(fit,"class") == "try-error") {

      if (algorithm == "nls") {

        warning("nls-fit failed. Input component table returned")
        #if (is.na(components$lambda[K])) {warning("A") } else { warning("B")}
        return(components)
      } else {

        if (verbose) warning("nls-fit failed. Falling back to det-results")
        #if (is.na(components$lambda[K])) {  warning("X") } else { warning("Y")}
        algorithm <- "det-fallback"
      }

    } else {

      n <- coef(fit)
      components$n <- n

      # add error estimations of fit as default and 'error.calculation=nls'-result
      components$n.error <- summary(fit)$parameters[, "Std. Error"][X]
    }
  } ########### end NLS ############



  ################### ERROR CALC ##################


  if ((error.calculation == "empiric")
      || (error.calculation == "poisson")
      || is.numeric(error.calculation)) {

    ### Calculate signal bin variances  ###
    I.err <- NULL
    if (error.calculation == "empiric") {

      # Calc reconstructed noise-free curve
      curve <- simulate_OSLcurve(components, curve, simulate.curve = FALSE)

      # Calc corrected sample variance
      for (i in X) {

        if (ch.start[i] == ch.end[i]) {

          # if signal bin consists just of one channel, assume Poisson statistics:
          I.err <- I[i]^0.5
        } else {

          # in all other cases: Use the corrected sample variance formula
          korrektor <- length(ch.start[i]:ch.end[i]) / (length(ch.start[i]:ch.end[i]) - 1)
          I.err <- c(I.err,
                     (korrektor * sum(curve$residual[ch.start[i]:ch.end[i]]^2))^0.5)
        }
      }
    } else {

      # Use poisson approach, add instrumental noise if defined
      if (!is.numeric(error.calculation)) error.calculation <- 0

      for (i in X) {

        I.err[i] <- (I[i] + length(ch.start[i]:ch.end[i]) * error.calculation^2 )^0.5
      }
    }

    components$bin.error <- I.err

    ### Propagation of uncertainty ###
    for (k in X) {
      sum.err <- 0

      for (i in X) {

        A.k <- A[[k]]

        # Differate the determinant term after I[j]
        A.k[i,] <- 0
        A.k[,k] <- 0
        A.k[i,k] <- 1
        sum.err <- sum.err + (det(A.k)*I.err[i])^2
      }

      components$n.error[k] <- sum.err^0.5 / det(D)
    }
  } ############ end ERROR CALC ############


  ########## component residuals  ###########

  # set the end of the record as the end of stimulation. Need not to be the same value as t.end
  stim.end <- curve$time[length(curve$time)]
  for (i in X) {

    components$n.residual[i] <- round(n[i] * exp(- stim.end * lambda[i]))
  }


  # Calculate average share at initial signal
  first_ch_signal <- n * (1 - exp(- lambda * dt))
  components$initial.signal <- round(first_ch_signal / sum(first_ch_signal), digits = 4)


  if (verbose) {
    col_set <- c("name", "n", "n.error", "n.residual", "initial.signal")
    col_set <- col_set[col_set %in% colnames(components)]
    print(subset(components, select = col_set))}

invisible(components)
}
