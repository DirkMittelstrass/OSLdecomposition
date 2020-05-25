#' Algebraic CW-OSL deconvolution
#'
#' The function provides a CW-OSL deconvolution for up to 3 first order OSL signal components.
#'
#' @param values [RLum.Data.Curve-class] or [data.frame] (**required**):
#' CW-OSL record. First row (x-axis) must contain time marks (seconds), Second row (y-axis) must contain signal values
#'
#' @param f.fast [numeric] (**required**):
#' fast component decay parameter (s^-1)
#'
#' @param f.medium [numeric] (*optional*):
#' medium component decay parameter (s^-1)
#'
#' @param f.slow [numeric] (*optional*):
#' slow component decay parameter (s^-1)
#'
#' @param t0 [numeric] (*with default*): point of time where the first integration interval starts
#'
#' @param t1 [numeric] (*optional*):
#' point of time dividing first and second integration interval; if not defined, value is determined by [calc_deconvolution.intervals]
#'
#' @param t2[numeric] (*optional*):
#' point of time dividing second and third integration interval; if not defined, value is determined by [calc_deconvolution.intervals]
#'
#' @param t3 [numeric] (*optional*):
#' point of time ending the third integration interval; if not defined, value is determined by [calc_deconvolution.intervals]
#'
#' @param error [character] (*with default*): integral error estimation approach, either **"empiric"** or **"poisson"** or a **numeric value**;
#' Per default the sample variance (see [var]) is used to calculate an **empiric** standard error for each integral which will be processed in the error propagation.
#' Alternatively the integral standard error can be calculated by assuming a **poisson** distributed signal error, known as *Shot noise*.
#' This is suitable if the lack of data points on the x-axis circumvent an empiric error estimation, like with spatial or spectral resolved CCD measurements.
#' Also the parameter can be set to a **numerical value** which will be handled as standard deviation per channel of the *Background noise*
#' and added to the poisson distributed *Shot noise*
#'
#' @param offset_value [numeric] (*with default*): signal offset per channel, usually detection or stimulation background or electronic offset to prevent negative counts.
#' Necessary if an accurate slow component output is wished
#'
#' @param negative.values.to.zero [logical] (*with default*): if TRUE negative n values are overwritten and set to n = 0;
#' negative n values are mathematically correct but physically debatable and can cause problems with [plot_GrowthCurve]
#'
#' @param output.plot [logical] (*with default*): returns a plot of the fitted curve

#' @param cex.global [numeric] (*with default*): global scaling factor
#'
#' @param sample_code [character] (*with default*): sample code used as subtitle for the plot
#'
#' @param ... further arguments and graphical parameters passed to [plot]
#'
#' @return
#' A [data.frame] object is returned containing the following elements:
#' \item{n}{[numeric] value equal to the number of filled traps at the start of the OSL stimulation}
#' \item{std.dev}{[numeric] standard deviation of n}
#' \item{n.residual}{[numeric] value equal to the number of remaining filled traps at the end of the OSL stimulation}
#' \item{sum}{[numeric] integral values}
#' \item{std.deviation}{[numeric] standard deviation of the integral values}
#' \item{signal.contribution}{[numeric] portion of signal contribution of one component to one integral}
#'
#' @section Changelog:
#' * winter 2012/13: first version, used for Bachelor thesis DM
#' * autumn 2013   : added empiric error estimation, shown in germanLED Freiberg 2013
#' * 2014-11-??, SK: formated into Rluminecence package standard
#' * 2014-11-07, DM: Binomial error propagation added
#' * 2018-05-01, DM: added interval time arguments and many little tweaks
#' * 2018-05-04, DM: added residuals for n values (necessary for slow component dosimetry)
#' * 2018-06-22, DM: added deconvolution of data sets with 1 or 2 components
#' * 2018-06-22, DM: added negative.values.to.zero and set it on TRUE as default (on request of C. Schmidt)
#' * 2018-07-05, DM: overworked error estimation; replaced binomial with poisson error approach; added auto-switch to poisson if integral length = 1; integrated background.noise into error
#'
#' @section ToDo:
#' * (!) generalised algorithmn which can deconvolve more than 3 components (!)
#' * allow for mixed empiric and poisson error estimations
#' * optional pseudoLM-OSL plotting
#'
#' @section Function version: 0.4.2
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
calc_CWOSLcomponents <- function(
  values,
  f.fast,
  f.medium = NA,
  f.slow = NA,
  t0 = 0,
  t1 = NA,
  t2 = NA,
  t3 = NA,
  error = "empiric",
  offset_value = 0,
  negative.values.to.zero = TRUE,
  output.plot = TRUE,
  cex.global = 0.6,
  sample_code = " ",
  ...
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

  # Set plot format parameters -----------------------------------------------
  # grep further paramters from the '...' argument
  extraArgs <- list(...) # read out additional arguments list

  main <- if("main" %in% names(extraArgs)) {extraArgs$main} else
  {"CW-OSL component seperation of quartz"}

  xlab <- if("xlab" %in% names(extraArgs)) {extraArgs$xlab} else
  {"time (s)"}

  ylab <- if("ylab" %in% names(extraArgs)) {extraArgs$ylab} else
  {"signal (cts)"}

  log <- if("log" %in% names(extraArgs)) {extraArgs$log} else
  {""}


  # preset result table, so the data format is static, even if not all cells are filled
  results <- data.frame(component = c("fast", "medium", "slow"),
                        n = c(NA, NA, NA),
                        std.dev = c(NA, NA, NA),
                        n.residual = c(NA, NA, NA),
                        Integration.parameters = c("1st integral", "2nd integral", "3rd integral"),
                        channels = c(NA, NA, NA),
                        t.start  = c(NA, NA, NA),
                        t.end = c(NA, NA, NA),
                        sum  = c(NA, NA, NA),
                        std.deviation = c(NA, NA, NA),
                        signal.contribution.fast = c(NA, NA, NA),
                        signal.contribution.medium  = c(NA, NA, NA),
                        signal.contribution.slow = c(NA, NA, NA))

  #### CHECK NUMBER OF COMPONENTS ####

  # no f.medium and f.slow are set (1-component case):
  if (is.na(f.medium) && is.na(f.slow)) {

    ##============================================================================##
    ## 1 COMPONENT
    ##============================================================================##

    if (is.na(t1)) t1 <- dt * length(x)

    ## calculating which channels shall be integrated:

    # first interval
    m1 <- 1 + floor(t0 / dt)
    m2 <- ceiling(t1 / dt)

    # recalculating time marks
    t0 <- (m1 - 1) * dt
    t1 <- m2 * dt

    # write interval properties in result table
    results$channels <- c(paste(m1, ":", m2), NA , NA)
    results$t.start  <- c(t0, NA, NA)
    results$t.end <- c(t1, NA, NA)

    # calculate release probability during interval
    P1.fast <- (exp(-t0 * f.fast) - exp(-t1 * f.fast))

    #calculate integral
    I1 <- sum(y[m1:m2]) - length(m1:m2) * offset_value

    # ... and solve:
    n.fast <- I1 / P1.fast
    if ((n.fast < 0) && negative.values.to.zero) n.fast <- 0

    ### component decay curve reconstruction (needed for numeric error calculation and the plots)
    #based on formula (7-2) with integration width = channel width:

    y.fast <- c(1: length(x))
    for(i in 1: length(x)){
      y.fast[i] <- n.fast * (exp(-(i-1)*dt*f.fast) - exp(-i*dt*f.fast))
    }
    y.virtual <- y.fast + offset_value

    ### ERROR CALCULATION ###

    # First step: numeric or analytic approach to get uncertainty of the integrals
    if(error=="empiric"){
      if (length(m1:m2) == 1){
        warning("interval contains just one channel, error estimation switched to Poisson type")
        error <- "poisson"
      } else {
        sigma1.square <- var(y[m1:m2] - y.virtual[m1:m2]) * length(m1:m2)}}

    if(error!="empiric"){
      if (!is.numeric(error)) error <- 0
      sigma1.square <- I1 + error^2 }

    # solution from gaussian error probagation method
    sigma.fast <- sqrt(sigma1.square) * n.fast / I1

    # calculate the remaining signal after the OSL stimulation ended
    residue.fast <- round(n.fast * exp(-dt * length(x) * f.fast))

    # write results
    results$n <- c(n.fast, NA, NA)
    results$std.dev <- c(sigma.fast, NA, NA)
    results$n.residual <- c(residue.fast, NA, NA)
    results$sum <- c(I1, NA, NA)
    results$std.deviation <- c(sigma1.square, NA, NA)
    results$signal.contribution.fast <- c(sum(y.fast[m1:m2])/I1, NA, NA)

  }

  ### f.medium or f.slow is missing (2-component case):
  if ((is.na(f.medium) || is.na(f.slow)) && !(is.na(f.medium) && is.na(f.slow)))  {

    ##============================================================================##
    ## 2 COMPONENTS
    ##============================================================================##

    # convert medium into slow and convert back later
    is_medium_given <- FALSE
    if (is.na(f.slow)) {
      f.slow <- f.medium
      is_medium_given <- TRUE
    }

    # calculate intervals
    if (is.na(t1)) {
      message("No integration intervals given. Call calc_deconvolution.intervals()")
      intervals <- calc_deconvolution.intervals(f.fast,
                                                NA,
                                                f.slow,
                                                channel.width = dt,
                                                channel.number = length(x),
                                                t.start = t0,
                                                t.end = t2)
      t1 <- intervals$t1
    }
    if (is.na(t2)) t2 <- dt * length(x)

    ## calculating which channels shall be integrated

    # first interval
    m1 <- 1 + floor(t0 / dt)
    m2 <- ceiling(t1 / dt)

    # second interval
    m3 <- m2 + 1
    m4 <- ceiling(t2 / dt)

    # recalculating time marks
    t0 <- (m1 - 1) * dt
    t1 <- m2 * dt
    t2 <- m4 * dt

    # write interval properties in result table
    results$channels = c(paste(m1, ":", m2),
                         paste(m3, ":", m4), NA)
    results$t.start  = c(t0, t1, NA)
    results$t.end = c(t1, t2, NA)


    #probability of luminescence photon appearing of each channel
    #formula (7-3) bachelor thesis Mittelstrass 2013
    P1.fast <- (exp(-t0 * f.fast) - exp(-t1 * f.fast))
    P1.slow <- (exp(-t0 * f.slow) - exp(-t1 * f.slow))

    P2.fast <- (exp(-t1 * f.fast) - exp(-t2 * f.fast))
    P2.slow <- (exp(-t1 * f.slow) - exp(-t2 * f.slow))

    #calculate the integrals
    I1 <- sum(y[m1:m2]) - length(m1:m2) * offset_value
    I2 <- sum(y[m3:m4]) - length(m3:m4) * offset_value
    I.vector <- c(I1,I2)

    #creating the solving matrices after Cramers rule
    #formulas (7-4) to (7-9) bachelor thesis Mittelstrass 2013
    A <- matrix(c(P1.fast, P1.slow, P2.fast, P2.slow), nrow = 2)
    A.fast <- A
    A.fast[1,] <- I.vector
    A.slow <- A
    A.slow[2,] <- I.vector

    # ... and solve:
    n.fast <- det(A.fast)/det(A)
    n.slow <- det(A.slow)/det(A)

    # check for negative values:
    if ((n.fast < 0) && negative.values.to.zero) n.fast <- 0
    if ((n.slow < 0) && negative.values.to.zero) n.slow <- 0

    ### component decay curve reconstruction (needed for numeric error calculation and the plots)
    y.fast <- c(1: length(x))
    y.slow <- c(1: length(x))

    #based on formula (7-2) with integration width = channel width:
    for(i in 1: length(x)){
      y.fast[i] <- n.fast * (exp(-(i-1)*dt*f.fast) - exp(-i*dt*f.fast))
      y.slow[i] <- n.slow * (exp(-(i-1)*dt*f.slow) - exp(-i*dt*f.slow))
    }
    y.virtual <- y.fast + y.slow + offset_value


    ### ERROR CALCULATION ###

    # First step: numeric or analytic approach to get uncertainty of the integrals
    if(error=="empiric"){
      if ((length(m1:m2) == 1) || (length(m3:m4) == 1) ){
        warning("intervals contain just one channel, error estimation switched to Poisson type")
        error <- "poisson"
      } else {
        # how is the noise-inflicted difference between the virtual and the real signal curve?
        # we calculate the variance for every integration area
        sigma1.square <- var(y[m1:m2] - y.virtual[m1:m2]) * length(m1:m2)
        sigma2.square <- var(y[m3:m4] - y.virtual[m3:m4]) * length(m3:m4)}}

    if(error!="empiric"){
      if (!is.numeric(error)) error <- 0
      sigma1.square <- I1 + error^2
      sigma2.square <- I2 + error^2}

    # Second step: gaussian error propagation to get uncertainty of the signal components

    #the partial differentials of the matrices A1, A2, A3:
    dA.matrix <- matrix(c(1:4), nrow = 2)
    for(i in 1:2){
      for(j in 1:2){
        A.temp <- A
        A.temp[i,] <- 0
        A.temp[i,j] <- 1
        dA.matrix[i,j] <- det(A.temp)^2
      }
    }

    # the big Gauss term:
    sigma.fast <- sqrt(dA.matrix[1,1]*sigma1.square + dA.matrix[1,2]*sigma2.square)/det(A)
    sigma.slow <- sqrt(dA.matrix[2,1]*sigma1.square + dA.matrix[2,2]*sigma2.square)/det(A)

    # calculate the remaining signal after the OSL stimulation ended
    residue.fast <- round(n.fast * exp(-dt * length(x) * f.fast))
    residue.slow <- round(n.slow * exp(-dt * length(x) * f.slow))

    # write results
    if (is_medium_given) {

      results$n <- c(n.fast, n.slow, NA)
      results$std.dev <- c(sigma.fast, sigma.slow, NA)
      results$n.residual <- c(residue.fast, residue.slow, NA)
      results$signal.contribution.medium <- c(sum(y.slow[m1:m2])/I1,
                                              sum(y.slow[m3:m4])/I2, NA)
    } else {

      results$n <- c(n.fast, NA, n.slow)
      results$std.dev <- c(sigma.fast, NA, sigma.slow)
      results$n.residual <- c(residue.fast, NA, residue.slow)
      results$signal.contribution.slow <- c(sum(y.slow[m1:m2])/I1,
                                            sum(y.slow[m3:m4])/I2, NA)
    }

    results$sum <- c(I1, I2, NA)
    results$std.deviation <- c(sigma1.square, sigma2.square, NA)
    results$signal.contribution.fast <- c(sum(y.fast[m1:m2])/I1,
                                          sum(y.fast[m3:m4])/I2, NA)

  }

if (!(is.na(f.fast) || is.na(f.medium) || is.na(f.slow))) {
  ##============================================================================##
  ## 3 COMPONENTS
  ##============================================================================##

  # calculate intervals
  if (is.na(t1) | is.na(t2)) {
    message("No integration intervals given. Call calc_deconvolution.intervals()")
    intervals <- calc_deconvolution.intervals(f.fast,
                                              f.medium,
                                              f.slow,
                                              channel.width = dt,
                                              channel.number = length(x),
                                              t0,
                                              t3)
    t1 <- intervals$t1
    t2 <- intervals$t2
  }
  if (is.na(t3)) t3 <- dt * length(x)

  ## calculating which channels shall be integrated

  # first interval
  m1 <- 1 + floor(t0 / dt)
  m2 <- ceiling(t1 / dt)

  # second interval
  m3 <- m2 + 1
  m4 <- ceiling(t2 / dt)

  # third interval
  m5 <- m4 + 1
  m6 <- floor(t3 / dt)

  # recalculating time marks
  t0 <- (m1 - 1) * dt
  t1 <- m2 * dt
  t2 <- m4 * dt
  t3 <- m6 * dt

  # write interval properties in result table
  results$channels = c(paste(m1, ":", m2),
                       paste(m3, ":", m4),
                       paste(m5, ":", m6))
  results$t.start  = c(t0, t1, t2)
  results$t.end = c(t1, t2, t3)

  #probability of luminescence photon appearing of each channel
  #formula (7-3) bachelor thesis Mittelstrass 2013
  P1.fast <- (exp(-t0 * f.fast) - exp(-t1 * f.fast))
  P1.medium <- (exp(-t0 * f.medium) - exp(-t1 * f.medium))
  P1.slow <- (exp(-t0 * f.slow) - exp(-t1 * f.slow))

  P2.fast <- (exp(-t1 * f.fast) - exp(-t2 * f.fast))
  P2.medium <- (exp(-t1 * f.medium) - exp(-t2 * f.medium))
  P2.slow <- (exp(-t1 * f.slow) - exp(-t2 * f.slow))

  P3.fast <- (exp(-t2 * f.fast) - exp(-t3 * f.fast))
  P3.medium <- (exp(-t2 * f.medium) - exp(-t3 * f.medium))
  P3.slow <- (exp(-t2 * f.slow) - exp(-t3 * f.slow))

  #calculate the integrals
  I1 <- sum(y[m1:m2]) - length(m1:m2) * offset_value
  I2 <- sum(y[m3:m4]) - length(m3:m4) * offset_value
  I3 <- sum(y[m5:m6]) - length(m3:m4) * offset_value
  I.vector <- c(I1,I2,I3)

  #creating the solving matrices after Cramers rule
  #formulas (7-4) to (7-9) bachelor thesis Mittelstrass 2013
  A <- matrix(c(P1.fast, P1.medium, P1.slow, P2.fast, P2.medium, P2.slow, P3.fast, P3.medium, P3.slow), nrow = 3)
  A.fast <- A
  A.fast[1,] <- I.vector
  A.medium <- A
  A.medium[2,] <- I.vector
  A.slow <- A
  A.slow[3,] <- I.vector

  # ... and solve:
  n.fast <- det(A.fast)/det(A)
  n.medium <- det(A.medium)/det(A)
  n.slow <- det(A.slow)/det(A)

  # check for negative values:
  if ((n.fast < 0) && negative.values.to.zero) n.fast <- 0
  if ((n.medium < 0) && negative.values.to.zero) n.medium <- 0
  if ((n.slow < 0) && negative.values.to.zero) n.slow <- 0

  ### component decay curve reconstruction (needed for numeric error calculation and the plots)
  y.fast <- c(1: length(x))
  y.medium <- c(1: length(x))
  y.slow <- c(1: length(x))

  #based on formula (7-2) with integration width = channel width:
  for(i in 1: length(x)){
    y.fast[i] <- n.fast * (exp(-(i-1)*dt*f.fast) - exp(-i*dt*f.fast))
    y.medium[i] <- n.medium * (exp(-(i-1)*dt*f.medium) - exp(-i*dt*f.medium))
    y.slow[i] <- n.slow * (exp(-(i-1)*dt*f.slow) - exp(-i*dt*f.slow))
  }
  y.virtual <- y.fast + y.medium + y.slow + offset_value


  ### ERROR CALCULATION ###


  # First step: numeric or analytic approach to get uncertainty of the integrals
  if(error=="empiric"){
    if ((length(m1:m2) == 1) || (length(m3:m4) == 1) || (length(m5:m6) == 1)){
      warning("intervals contain just one channel, error estimation switched to Poisson type")
      error <- "poisson"
    } else {
      # how is the noise-inflicted difference between the virtual and the real signal curve?
      # we calculate the variance for every integration area
      sigma1.square <- var(y[m1:m2] - y.virtual[m1:m2]) * length(m1:m2)
      sigma2.square <- var(y[m3:m4] - y.virtual[m3:m4]) * length(m3:m4)
      sigma3.square <- var(y[m5:m6] - y.virtual[m5:m6]) * length(m5:m6)}}

  if(error!="empiric"){
    if (!is.numeric(error)) error <- 0
    sigma1.square <- I1 + error^2
    sigma2.square <- I2 + error^2
    sigma3.square <- I3 + error^2}

  # Binomial approach: Another alternative of error estimation
  # sigma1.square <- n.fast*P1.fast*(1 - P1.fast) + n.medium*P1.medium*(1 - P1.medium) + n.slow*P1.slow*(1 - P1.slow) + (t1 - t0)*background.noise^2
  # sigma2.square <- n.fast*P2.fast*(1 - P2.fast) + n.medium*P2.medium*(1 - P2.medium) + n.slow*P2.slow*(1 - P2.slow) + (t2 - t1)*background.noise^2
  # sigma3.square <- n.fast*P3.fast*(1 - P3.fast) + n.medium*P3.medium*(1 - P3.medium) + n.slow*P3.slow*(1 - P3.slow) + (t3 - t2)*background.noise^2


  ## Second step: gaussian error propagation to get uncertainty of the signal components

  # the partial differentials of the matrices A1, A2, A3:
  dA.matrix <- matrix(c(1:9), nrow = 3)
  for(i in 1:3){
    for(j in 1:3){
      A.temp <- A
      A.temp[i,] <- 0
      A.temp[i,j] <- 1
      dA.matrix[i,j] <- det(A.temp)^2
    }
  }

  # the big Gauss term:
  sigma.fast <- sqrt(dA.matrix[1,1]*sigma1.square + dA.matrix[1,2]*sigma2.square + dA.matrix[1,3]*sigma3.square)/det(A)
  sigma.medium <- sqrt(dA.matrix[2,1]*sigma1.square + dA.matrix[2,2]*sigma2.square + dA.matrix[2,3]*sigma3.square)/det(A)
  sigma.slow <- sqrt(dA.matrix[3,1]*sigma1.square + dA.matrix[3,2]*sigma2.square + dA.matrix[3,3]*sigma3.square)/det(A)

  # calculate the remaining signal after the OSL stimulation ended
  residue.fast <- round(n.fast * exp(-dt * length(x) * f.fast))
  residue.medium <- round(n.medium * exp(-dt * length(x) * f.medium))
  residue.slow <- round(n.slow * exp(-dt * length(x) * f.slow))

  # write results
  results$n <- c(n.fast, n.medium, n.slow)
  results$std.dev <- c(sigma.fast, sigma.medium, sigma.slow)
  results$n.residual <- c(residue.fast, residue.medium, residue.slow)
  results$sum <- c(I1, I2, I3)
  results$std.deviation <- c(sigma1.square, sigma2.square, sigma3.square)
  results$signal.contribution.fast <- c(sum(y.fast[m1:m2])/I1,
                                        sum(y.fast[m3:m4])/I2,
                                        sum(y.fast[m5:m6])/I3)
  results$signal.contribution.medium <- c(sum(y.medium[m1:m2])/I1,
                                          sum(y.medium[m3:m4])/I2,
                                          sum(y.medium[m5:m6])/I3)
  results$signal.contribution.slow <- c(sum(y.slow[m1:m2])/I1,
                                        sum(y.slow[m3:m4])/I2,
                                        sum(y.slow[m5:m6])/I3)

}


  ##==========================================================================##
  ## PLOTTING
  ##==========================================================================##
  if(output.plot==TRUE){

    ##set colors gallery to provide more colors
    col<-unlist(colors())
    col<-col[c(261,552,51,62,76,151,451,474,654)]

    ##set plot frame
    layout(matrix(c(1,2,3),3,1,byrow=TRUE),c(1.6,1,1), c(1,0.3,0.4),TRUE)
    par(oma=c(1,1,1,1),mar=c(0,4,3,0),cex=cex.global)

    ##==uppper plot==##
    ##open plot area

    plot(NA,NA,
         xlim=c(min(x),max(x)),
         ylim=if(log=="xy"){c(1,max(y))}else{c(0,max(y))},
         xlab="",
         xaxt="n",
         ylab=if(missing(ylab)==TRUE){
           paste("OSL [cts/",length(x)/max(x)," s]",sep="")}else{
             ylab
           },
         main=main,
         sub="",
         log=log)

    ##plotting measured signal
    points(x,y,pch=20, col="grey")

    ##add additional labeling (fitted function)
    mtext(side=3,sample_code,cex=0.7*cex.global)

    ##plot decay functions
    lines(x,y.virtual, col=col[1])
    lines(x,y.fast, col=col[2])
    if(!is.na(f.medium)) lines(x,y.medium, col=col[3])
    if(!is.na(f.slow))lines(x,y.slow, col=col[4])

    ##plot legend
    #legend.caption <- c("sum of components","fast component", "medium component", "slow components")
    legend.caption <- c("reconstructed OSL decay curve",
                        paste("fast component:   ",
                              "n = ", format(n.fast, digits = 3, scientific = TRUE),
                              " +- ", format(sigma.fast, digits = 3, scientific = TRUE)))
    if(!is.na(f.medium)) legend.caption <- c(legend.caption, paste("medium component: ",
                                                   "n = ", format(n.medium, digits = 3, scientific = TRUE),
                                                   " +- ", format(sigma.medium, digits = 3, scientific = TRUE)))
    if(!is.na(f.slow)) legend.caption <- c(legend.caption, paste("slow component:   ",
                                                                 "n = ", format(n.slow, digits = 3, scientific = TRUE),
                                                                 " +- ", format(sigma.slow, digits = 3, scientific = TRUE)))


    legend("topright",legend.caption, lty=rep(1,length(legend.caption),NA),lwd=1.5, col=col[1:length(legend.caption)], bty="n")

      ##==lower plot==##
      ##plot residuals
      par(mar=c(4.2,4,0,0))
      plot(x,y - y.virtual,
           xlim=c(min(x),max(x)),
           xlab=if(missing(xlab)==TRUE){if(log=="x" | log== "xy"){"log time [s]"}else{"time [s]"}}else{xlab},
           type="l",
           col="darkgrey",
           ylab="residual [a.u.]",
           lwd=2,
           log=if(log=="x" | log=="xy"){log="x"}else{""}
      )

      ##add 0 line
      abline(h=0)

#
#       ##------------------------------------------------------------------------##
#       ##++component to sum contribution plot ++##
#       ##------------------------------------------------------------------------##
#
#       ##plot component contribution to the whole signal
#       #open plot area
#       par(mar=c(4,4,3.2,0))
#       plot(NA,NA,
#            xlim=c(min(x),max(x)),
#            ylim=c(0,100),
#            ylab="contribution [%]",
#            xlab=if(missing(xlab)==TRUE){if(log=="x" | log=="xy"){"log time [s]"}else{"time [s]"}}else{xlab},
#            main="Component Contribution To Sum Curve",
#            log=if(log=="x" | log=="xy"){log="x"}else{""})
#
#       stepping <- seq(3,length(component.contribution.matrix),2)
#
#       for(i in 1:length(I0)){
#
#         polygon(c(component.contribution.matrix[,1],
#                   component.contribution.matrix[,2]),
#                 c(component.contribution.matrix[,stepping[i]],
#                   component.contribution.matrix[,stepping[i]+1]),
#                 col = col[i+1])
#       }
#       rm(stepping)


#    }#end if try-error for fit
  }

##============================================================================##
## Return Values
##============================================================================##

  return(results)
}
