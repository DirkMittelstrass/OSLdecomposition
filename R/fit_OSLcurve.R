#' Fit global mean curve to get decay constants
#'
#' @param curve
#' @param K.max
#' @param F.threshold
#' @param stimulation.intensity
#' @param stimulation.wavelength
#' @param verbose
#'
#' @section Changelog:
#' * 2019-02-14, DM: First version
#' * 2019-03-15, DM: Seperated 'decompose_OSLalternatively()'
#' * 2019-04-29, DM: Added Blucszs & Adamiec-like approach using numOSL::decomp (Peng et al. 2014)
#' * 2019-05-14, DM: Added "fit_OSLLifeTimes" approach from Luminescence package 0.9; Corrected and improved numOSL approach; Deleted nls.default approach
#' * 2019-06-28, DM: Deleted "fit_OSLLifeTimes" approach. Added stretched exponentials for testing. Added overview plot
#' * 2019-10-07, DM: Streamlined function; added optional background fitting
#' * 2019-10-08, DM: Seperated plotting to plot_PhotoCrosssections()
#' * 2020-04-04, DM: Extended output list (curve & arguments)
#' * 2020-04-06, DM: Extended print output and made some  tweaks. Replaced 'SAR.compatible' with 'fully.bleached'
#' * 2020-05-05, DM: Replaced bolean 'fully.bleached' with numeric 'bleaching.grade'
#' * 2020-08-05, DM: Added DEoptim + nlsLM algorithm
#' * 2020-08-10, DM: Optional parallel computing enabled
#'
#' @section ToDo:
#' * Write documentation
#' * Reactivate optional background level fitting
#' * Enable list of DE and LM control parameters as argument
#' * Enable optional weighted fitting and give out reduced Chi²
#'
#' @section Last changed. 2020-08-10
#'
#' @author
#' Dirk Mittelstrass, \email{dirk.mittelstrass@@luminescence.de}
#'
#' @return
#' @export
#'
#' @examples
fit_OSLcurve <- function(
  curve,
  K.max = 3,
  F.threshold = 150,
  stimulation.intensity = 30,
  stimulation.wavelength = 470,
  algorithm = "DE+LM", # "DE", "DE+LM", "numOSL"
  parallel.computing = TRUE, # set to FALSE before release!
  verbose = TRUE,
  output.complex = TRUE
){
  # measure computing time
  time.start <- Sys.time() # Delete before release

  ################### Prepare input data ###########################

  if(is(curve, "RLum.Data.Curve") == FALSE & is(curve, "data.frame") == FALSE){
    stop("[fit_CWCurve()] Input object is not of type 'RLum.Data.Curve' or 'data.frame'!", call. = FALSE)
  }

  if(is(curve, "RLum.Data.Curve") == TRUE){

    time <- curve@data[,1]
    signal <- curve@data[,2]

  }else{

    ##set x and y values
    time <- curve$time
    signal <- curve$signal
  }

  data <- data.frame(time = time, signal = signal)

  channel_width <- time[2] - time[1]

  # check if time beginns with zero and add channel_width if the case
  if (time[1] == 0)  time <- time + channel_width


  # check if any signal = 0 or smaller and set them to 0.1
  if (any(signal <= 0)) {

    if (verbose) warning("[fit_OSLcurve] Signal values equal or smaller than zero detected. Replaced with 0.1")
    signal[signal <= 0] <- 0.1
  }

  # Some presets ...
  K_selected <- 0
  #x <- 1
  X <- c(1:K.max)
  RSS <- NA
  RSS_old <- Inf
  fittings <- list(NULL)
  component_tables <- list(NULL)
  F_table <- data.frame(NULL)
  F_table_print <- data.frame(NULL)
  plot_data <- data.frame(NULL)
  lambda <- c(NULL)
  info_text <- ""

  # supress warnings in the whole script
  options(warn = -1)

  # prepare printed table
  if (verbose) cat("Cycle \t", paste(paste0("       f_", X), collapse = "  "),
                   "        RSS     F-value\n")

  ###################### Subfunction for DE minimisation ###########################################
  calc_RSS <- function(lambda_vector, curve = data, decompose_algorithm = "det"){

    # The linear part of HELA (see Bluszcz & Adamiec 2006) is performed by decompose_OSLcurve()
    RSScomponents <- decompose_OSLcurve(curve,
                                        lambda_vector,
                                        algorithm = decompose_algorithm,
                                        error.calculation = "none",
                                        verbose = FALSE)

    # Now add the residual curve to the input curve ...
    curve <- simulate_OSLcurve(RSScomponents,
                               template.curve = curve,
                               simulate.curve = FALSE)

    # ... and calculate the residual sum of squares (RSS)
    RSS <- sum(curve$residual^2)
    if (is.na(RSS) || (RSS <= 0)) RSS <- Inf
    return(RSS)}

  ###################### Reduced Chi² ###############################################################
  calc_Chi2 <- function(components, curve = data, N_curves = M, detec_noise = 0, K = K){

    curve <- simulate_OSLcurve(components,
                                    template.curve = curve)

    RS <- curve$residual^2 * N_curves / (abs(curve$sum) + detec_noise)
    Chi2 <- sum(RS) / (length(RS) - K * 2)
    return(Chi2)}

  ###################### Photo-Ionisation Crosssections #############################################
  build_component_table <- function(lambda_vector, lambda_err, curve = data){

    Y <- 1:length(lambda_vector)

    # Calc photon energy: E = h*v  [W*s^2 * s^-1 = W*s = J]
    E <-6.62606957e-34 * 299792458 * 10^9 / stimulation.wavelength

    # Calc photon flux of stimulation light: Flux = I / E  [W/cm^2 / W*s = 1/s*cm^2]
    Flux <- stimulation.intensity / (E * 1000)

    # Calc cross-sections: Sigma = lambda / Flux  [s^-1 / 1/s*cm^2 = cm^2]
    cross.section <- lambda_vector / Flux
    cross.relative <- round(cross.section / cross.section[1], digits=4)

    ### NAME COMPONENTS ###
    # default names:
    name <- paste0("Component ", Y)

    if ((stimulation.wavelength >= 460) && (stimulation.wavelength <= 485)  ) {

      # Rename components according to literature
      for (i in Y) {

        c <- cross.section[i]

        # Autonaming uses Table 1 in Durcan & Duller (2011)
        # Minimum = lowest value - 2-sigma
        # Maximum = highest value + 2-sigma
        # Exception are Ultrafast and Slow4 which are not well defined in literature and guessed freely
        if ((c > 1e-16) && (c < 1e-15)) name[i] <- "Ultrafast"
        if ((c > 1.9e-17) && (c < 3.1e-17)) name[i] <- "Fast"
        if ((c > 3e-18) && (c < 9e-18)) name[i] <- "Medium"
        if ((c > 1e-18) && (c < 1.85e-18)) name[i] <- "Slow1"
        if ((c > 1.1e-19) && (c < 4e-19)) name[i] <- "Slow2"
        if ((c > 1e-20) && (c < 4.67e-20)) name[i] <- "Slow3"
        if ((c > 1e-21) && (c < 1e-20)) name[i] <- "Slow4"}}


    # Check for double-naming
    name[duplicated(name)] <- paste0(name[duplicated(name)], ".a")

    # Check again
    name[duplicated(name)] <- paste0(substr(name[duplicated(name)], 1, nchar(name[duplicated(name)]) - 2), ".b")

    # And again
    name[duplicated(name)] <- paste0(substr(name[duplicated(name)], 1, nchar(name[duplicated(name)]) - 2), ".c")

    # How much is the component bleached during stimulation?
    bleaching.grade <- round(1 - exp(- lambda_vector * max(time)), digits = 4)

    # Decay with zero or negative error had no correct error estimation
    lambda_err[lambda_err <= 0] <- NA

    ##### Build result table #####
    components <- data.frame(name = name,
                             lambda = lambda_vector,
                             lambda.error = lambda_err,
                             cross.section = cross.section,
                             cross.relative = cross.relative,
                             bleaching.grade = bleaching.grade)

    row.names(components) <- Y

    # Go the easy way to extract additional information from the fitting
    components <- decompose_OSLcurve(curve = curve,
                                     components = components,
                                     verbose = FALSE)

  } ### END building tables for the various cases ###



  #----------------------------------------------------------------------------------------------------#
  #------------------------------------- K = K + 1 cycle ----------------------------------------------#
  #----------------------------------------------------------------------------------------------------#
  for (K in X) {

    lambda_error <- rep(0, K)

    if (algorithm == "numOSL") {
      ########################################### numOSL ##############################################
      fit <- numOSL::decomp(cbind(time, signal),
                            ncomp = K,
                            plot = FALSE,
                            weight = FALSE)

      if (fit$message == 0){

        lambda <- fit$LMpars[,3][1:K]
        RSS <- fit$value
        fit_sucessful <- TRUE
      } else {
        if (verbose) cat("Left loop. numOSL::decomp() fitting failed at K =", K, "\n")

        # leave loop
        break
      }

    } else {

      ########################################### DEoptim ##############################################

      # Divide the DE parameter space a the decay values of the previous cycle
      # Additional constraints:
      # - no negative values (decay >= 0)
      # - no superfast decays, that the channel frequency couldn't resolve it (decay <= 3 / channel_width)
      lower_lambda <- c(lambda, 0)
      upper_lambda <- c(3 / channel_width, lambda)

      # Perform differential evolution (DE). As minimisation function, use calc_RSS()
      DE_min <- try(DEoptim::DEoptim(
        fn = calc_RSS,
        lower = lower_lambda,
        upper = upper_lambda,
        control = DEoptim::DEoptim.control(
          NP = K * 15,
          strategy = 2,
          itermax = 100,
          c = 0.2,
          reltol = 1e-4,
          steptol = 10,
          parallelType = parallel.computing,
          packages = c("OSLdecomposition"),
          parVar = c("data"),
          trace = FALSE)))

      # Did the DE algorithm break?
      if (attr(DE_min, "class") == "try-error") {
        if (verbose) cat("Differential evolution failed at K =", K, ". Algorithm stopped.\n")

        # leave loop and set flag to delete unused columns
        fit_failed <- TRUE
        break}

      # Otherwise, extract results
      fit <- list(NULL)
      fit[["DE"]] <- DE_min
      lambda <- DE_min$optim$bestmem
      RSS <- DE_min$optim$bestval

      ########################################### nlsLM ##############################################

      if (algorithm == "DE+LM") {

        # We need the signal intensities of the components as start values for the LM fitting
        DE_components <- decompose_OSLcurve(data,
                                            lambda,
                                            algorithm = "det",
                                            error.calculation = "none",
                                            verbose = FALSE)
        n <- DE_components$n

        ### Create fit formula ###
        n.names <- paste0("n.",1:K)
        lambda.names <- paste0("lambda.",1:K)

        # now creat the optimization formula
        fit.formula <- as.formula(paste0("signal ~ ",
                                         paste(n.names," * (exp(-", lambda.names," * (time - ", channel_width,")) - exp(-", lambda.names," * time))",
                                               collapse=" + ")))

        # Name the vectors to allow the correct value allocation
        names(n) <- n.names
        names(lambda) <- lambda.names

        ### Apply LM algorithm  ###
        LM_fit <- try(minpack.lm::nlsLM(fit.formula,
                                        data = data,
                                        start = c(n, lambda)),
                      silent = TRUE)

        if (attr(LM_fit, "class") == "try-error") {

          if (verbose) cat("Levenberg-Marquardt fitting failed at K =", K, ". Differential evolution result are used:\n")

        } else {

          fit[["LM"]] <- LM_fit
          lambda <- summary(LM_fit)$parameters[paste0("lambda.", 1:K),"Estimate"]
          RSS <- LM_fit$m$deviance()
          lambda_error <- summary(LM_fit)$parameters[paste0("lambda.", 1:K),"Std. Error"]}
      }
    }

    # save the raw fitting objects
    fittings[[K]] <- fit

    # identify the components and build the Component table
    component_tables[[K]] <- build_component_table(lambda, lambda_error)

    # Add values to [plot_Photocrosssections()] ploting table
    plot_data <- rbind(plot_data,
                       data.frame(lambda = lambda,
                                  lambda.low = lambda - lambda_error,
                                  lambda.up = lambda + lambda_error,
                                  name = factor(paste0("Fit with K = ", K)),
                                  x = K))
    #x <- x + 1

    ### F-test ###
    F_value <- 0.5*(RSS_old - RSS) / (RSS / (length(signal) - 2 * K))
    RSS_old <- RSS

    # Create live console output
    table_row <- c(rep("          ", K.max),
                   formatC(RSS, digits = 4, width = 10),
                   formatC(F_value, digits = 4, width = 10))
    table_row[1:length(lambda)] <- formatC(lambda, digits = 4, width = 10)
    if (verbose) cat("K =", K, "\t", paste(table_row, collapse = "  "), "\n")

    # Build output F-test table for [report_Step2.rmd] script
    F_table_print <- rbind(F_table_print, table_row, stringsAsFactors = FALSE)
    F_table <- rbind(F_table, (c(lambda[X], RSS, F_value)))

    # Stop fitting if K - 1 was the correct model
    if (F_value <= F.threshold) {
      info_line <- paste0("Left loop because F-test value (F = ", formatC(F_value),
                          ") fell below threshold value (F = ", F.threshold,")\n")
      info_text <- paste0(info_text, info_line)
      if (verbose) cat(info_line)
      break}

    # If the current iteration cycle succeded until this point, it must be the best fit so far. Therefore:
    K_selected <- K

    if (K == K.max) {
      info_line <- paste0("Left loop because maximum number of allowed components K is reached\n")
      info_text <- paste0(info_text, info_line)
      if (verbose) cat(info_line)}

  } #---------------------------------------- End cycle -----------------------------------------------

  if ((K_selected == 0)||(nrow(F_table_print) == 0)) stop("[fit_OSLcurve] No sucessful fit")

  # Give F tables approbiate headers
  colnames(F_table) <- c(paste0("f_", X),"RSS","F-value")
  colnames(F_table_print) <- c(paste0("f_", X),"RSS","F-value")

  # Delete unused columns
  if (nrow(F_table) < K.max) {
    F_table <- F_table[,-c((nrow(F_table) + 1):K.max)]
    F_table_print <- F_table_print[,-c((nrow(F_table) + 1):K.max)]}

  # Standard set of components is the on chosen by the F-test
  components <- component_tables[[K_selected]]

  ######### Further console output ######
  if (verbose) {

    cat(paste0("->  ", K_selected,"-component model choosen"), "\n")

    cat("\nSignal component parameter:\n")
    # Reduce amount of information in the table to avoid user irritation, also round some values
    print_components <- subset(components,
                               select = c(name, lambda, lambda.error,
                                          n, n.error, cross.section, initial.signal, bleaching.grade))

     for (col in 2:ncol(print_components)) {
      print_components[,col] <- formatC(print_components[,col], digits = 4)}

    print.data.frame(print_components, row.names = FALSE)

    # Give advice which components are suited for further analysis
    bleach <- components$bleaching.grade
    if (any(bleach < 0.99)) {

      if (sum(bleach < 0.99) == nrow(components)) {
        cat("WARNING: No component was fully bleached during stimulation. Check your experimental settings!\n")

      } else if (sum(bleach < 0.99)  == 1) {

        cat(paste0(components$name[bleach < 0.99],
                   " was not fully bleached during stimulation and is not recommended for further dose evaluation\n"))
      } else {

        cat(paste0(paste(components$name[bleach < 0.99], collapse = ", "),
                   " were not fully bleached during stimulation and are not recommended for further dose evaluation\n"))}}

     # print computing time
    cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")

  }

  ######################### Return values #######################

  if (output.complex) {

    output <- list(decay.rates = components$lambda,
                   K.selected = K_selected,
                   F.test = F_table,
                   F.test.print = F_table_print,
                   info.text = info_text,
                   component.tables = component_tables,
                   curve = curve,
                   fit.results = fittings,
                   plot.data = plot_data,
                   parameters = list(K.max = K.max,
                                     F.threshold = F.threshold,
                                     stimulation.intensity = stimulation.intensity,
                                     stimulation.wavelength = stimulation.wavelength,
                                     algorithm = algorithm))
    invisible(output)

  } else {
    invisible(components)}

}
