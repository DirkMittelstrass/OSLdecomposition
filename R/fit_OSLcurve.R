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
#' * 2019-02-14 First version
#' * 2019-03-15 Seperated 'decompose_OSLalternatively()'
#' * 2019-04-29 Added Blucszs & Adamiec-like approach using numOSL::decomp (Peng et al. 2014)
#' * 2019-05-14 Added "fit_OSLLifeTimes" approach from Luminescence package 0.9; Corrected and improved numOSL approach; Deleted nls.default approach
#' * 2019-06-28 Deleted "fit_OSLLifeTimes" approach. Added stretched exponentials for testing. Added overview plot
#' * 2019-10-07 Streamlined function; added optional background fitting
#' * 2019-10-08 Seperated plotting to plot_PhotoCrosssections()
#' * 2020-04-04 Extended output list (curve & arguments)
#' * 2020-04-06 Extended print output and made some  tweaks. Replaced 'SAR.compatible' with 'fully.bleached'
#'
#' @section ToDo:
#' * ! Replace numOSL::decomp by own code. Then enable sigma.lambda calculation and replace "easy way" of intensity calculation !
#' * Documentation
#' * Test background value fitting
#'
#' @section Last changed. 2020-04-06
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @return
#' @export
#'
#' @examples
fit_OSLcurve <- function(
  curve,
  K.max = 5,
  F.threshold = 50,
  stimulation.intensity = 30,
  stimulation.wavelength = 470,
  applied.time.cut = FALSE,
  weight.Chi = FALSE,
  background.fitting = FALSE,
  verbose = TRUE,
  output.complex = TRUE,
  output.plot = TRUE
){

  ################### Prepare input data ###########################

  library(numOSL)

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

  channel.width <- time[2] - time[1]

  # check if time beginns with zero and add channel.width if the case
  if (time[1] == 0)  time <- time + channel.width


  # check if any signal = 0 or smaller and set them to 0.1
  if (any(signal <= 0)) {

    if (verbose) warning("[fit_OSLcurve] Signal values equal or smaller than zero detected. Replaced with 0.1")
    signal[signal <= 0] <- 0.1
  }

  plot.data <- data.frame(NULL)
  K.selected <- 0
  x <- 1

  # supress warnings in the whole script
  options( warn = -1 )

  ################### Genetic algorithm fitting ###########################


  X <- c(1:K.max)

  Chi2_old <- Inf
  best_fit <- 0
  fit.list <- list(NULL)
  results <- data.frame(NULL)

  for (i in X) {

    fit <- decomp(cbind(time, signal),
                  ncomp = i,
                  constant = background.fitting,
                  plot = FALSE,
                  weight = weight.Chi)

    fit.list[[i]] <- fit

    # was fit sucessfull?
    if (fit$message == 0) {

      Chi2 <- fit$value
      F_value <- 0.5*(Chi2_old - Chi2) / (Chi2 / (length(signal) - 2*i))

      Chi2_old <- Chi2

      if (F_value > F.threshold) K.selected <- i

      # Add values to ploting table
      plot.data <- rbind(plot.data,
                         data.frame(lambda = fit$LMpars[,3],
                                    lambda.low = fit$LMpars[,3],
                                    lambda.up = fit$LMpars[,3],
                                    name = factor(paste0("Fit with K = ", i)),
                                    x = x))
      x <- x + 1

      # Build overview table
      result_vector <- fit$LMpars[,3][X]
      if (background.fitting == TRUE) result_vector <- c(result_vector, fit$constant[1])
      result_vector <- c(result_vector,
                         fit$value,
                         fit$FOM,
                         F_value)

      results <- rbind(results, result_vector)

      }
  }

  if (K.selected == 0) stop("no sucessful fit")

  ###### Build component tables names and estimate photo-ionisation crosssections #####

  C.list <- list(NULL)
  for (k in 1:nrow(results)) {

    lambda <- fit.list[[k]][["LMpars"]][,3]

    #n <- NA[1:k]
    n.error <- NA[1:k]
    ##### Calc photoionisation crosssections ######

    # Calc photon energy: E = h*v  [W*s^2 * s^-1 = W*s = J]
    E <-6.62606957e-34 * 299792458 * 10^9 / stimulation.wavelength

    # Calc photon flux of stimulation light: Flux = I / E  [W/cm^2 / W*s = 1/s*cm^2]
    Flux <- stimulation.intensity / (E * 1000)

    # Calc crosssections: Sigma = lambda / Flux  [s^-1 / 1/s*cm^2 = cm^2]
    cross.section <- lambda / Flux
    cross.relative <- round(cross.section / cross.section[1], digits=4)



    ### NAME COMPONENTS ###

    # default names:
    name <- paste0("Component ",1:k)

    if ((stimulation.wavelength >= 460) && (stimulation.wavelength <= 485)  ) {

      # Rename components according to literature
      for (i in c(1:k)) {

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
        if ((c > 1e-21) && (c < 1e-20)) name[i] <- "Slow4"
      }
    }


    # Check for double-naming
    name[duplicated(name)] <- paste0(name[duplicated(name)], ".a")

    # Check again
    name[duplicated(name)] <- paste0(substr(name[duplicated(name)], 1, nchar(name[duplicated(name)]) - 2), ".b")

    # And again
    name[duplicated(name)] <- paste0(substr(name[duplicated(name)], 1, nchar(name[duplicated(name)]) - 2), ".c")


    # Is the component fully bleached (less than 1% residual) during stimulation?
    fully.bleached <- rep(1,k)
    for (y in 1:k) if(exp(- lambda[y] * max(time)) > 0.01) fully.bleached[y] <- 0

    ##### Build result table #####
    components <- data.frame(name = name,
                             lambda = lambda,
                             cross.section = cross.section,
                             cross.relative = cross.relative,
                             fully.bleached = fully.bleached)

    row.names(components) <- 1:k

    # Go the easy way to extract additional information from the fitting
    components <- decompose_OSLcurve(curve = curve,
                                     components = components,
                                     background.fitting = background.fitting,
                                     verbose = FALSE)

    C.list[[k]] <- components
  } ### END building tables for the various cases ###

  # Standard set of components is the on chosen by the F-test
  components <- C.list[[K.selected]]

  # Reduce amount of information in the table to avoid user irritation, also round some values
  components <- subset(components, select = c(name, lambda, cross.section, initial.signal, n, fully.bleached))
  N.not_fully_bleached <- (nrow(components) - sum(components$fully.bleached))
  components$fully.bleached[components$fully.bleached == 1] <- "true"
  components$fully.bleached[components$fully.bleached == 0] <- "false"
  components$n <- round(components$n)
  components$cross.section <- formatC(components$cross.section)

  ##### Format F-tables #####
  # create a better looking table for publishing purposes

  # add "K=" column
  # results <- cbind(data.frame(K = 1:nrow(results)), results)


  #results.print <- results
  # Chi² in Rmarkdown: $\\chi^2$

  if (background.fitting == TRUE) {

    colnames(results) <- c(paste0("k_", X),"background","Chi2","FOM","F.value")
    #colnames(results.print) <- c("  K  ", paste0("$\\lambda_", X,"$ $(s^{-1})$"),"background","RSS","FOM","$F_K$")
  } else {

    colnames(results) <- c(paste0("k_", X),"Chi2","FOM","F.value")
    #colnames(results.print) <- c("  K  ", paste0("$\\lambda_", X,"$ $(s^{-1})$"), "RSS","FOM","$F_K$")
  }

  results$Chi2 <- formatC(results$Chi2)
  results$F.value <- formatC(results$F.value)

  results[is.na(results)] <- ""
  results[results == "Inf"] <- ""


 # for (y in 1:ncol(results)) {
#
 #   results[grepl("NA", results[, y]), y] <- ""
 #   results[grepl("Inf", results[, y]), y] <- ""
 #}

 # results.print <- subset(results.print, select = -FOM)

  # Print F-table and Component table
  if (verbose) {

    writeLines("F-test table:")
    print(subset(results, select = -FOM))
    writeLines(paste0("-->  ", K.selected,"-component model choosen"))

    writeLines("")
    writeLines("Signal component table:")
    print(components)

    # Give advice which components are suited for further analysis
    if (N.not_fully_bleached > 0) {

      if (N.not_fully_bleached == nrow(components)) {
        writeLines("WARNING: No component was fully bleached during stimulation. Check your experimental settings!")

      } else if (N.not_fully_bleached == 1) {

        writeLines(paste0(components$name[components$fully.bleached == "false"],
                          " was not fully bleached during stimulation and is not recommended for further dose evaluation"))
      } else {

        writeLines(paste0(paste(components$name[components$fully.bleached == "false"], collapse = ", "),
                          " were not fully bleached during stimulation and are not recommended for further dose evaluation"))
      }
    }
  }

  ######################### Output #######################

  if (output.complex) {

    output <- list(curve = curve,
                   K.selected = K.selected,
                   components = components,
                   F.test = results,
                #   F.test.print = results.print,
                   fit.results = fit.list,
                   case.tables = C.list,
                   plot.data = plot.data,
                   fit.arguments = list(K.max = K.max,
                                        F.threshold = F.threshold,
                                        stimulation.intensity = stimulation.intensity,
                                        stimulation.wavelength = stimulation.wavelength,
                                        applied.time.cut = applied.time.cut,
                                        weight.Chi = weight.Chi,
                                        background.fitting = background.fitting))

    invisible(output)

  } else {

    invisible(components)
  }

}
