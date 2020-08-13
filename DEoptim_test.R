library(OSLdecomposition)
#set.seed(1234)

verbose <- TRUE
M <- 500


components <- data.frame(lambda = c(2, 0.5, 0.05),
                              n = c(2000, 2000, 100000))
K_max <- length(components$lambda) + 1

cat("input n:      \t", paste(formatC(components$n, digits = 6), collapse = "\t"), "\n")

cat("input lambdas:\t", paste(formatC(components$lambda, digits = 6), collapse = "\t"), "\n")


C <- simulate_OSLcurve(components, simulate.curve = TRUE)
for (i in 2:M) {
  temp <- simulate_OSLcurve(components, simulate.curve = TRUE)
  C$signal <- C$signal + temp$signal
}
C$signal <- C$signal / M
channel_width <- C$time[2] - C$time[1]




########################################################################


#plot(curve$time, curve$signal)
#cat(Chi2(c(3,0.1,0.01), C), "\n")

# measure computing time
time.start <- Sys.time()

best_lambdas <- c(NULL)

for (K in 1:K_max) {
  lower_lambda <- c(best_lambdas, 0)
  upper_lambda <- c( 3 / channel_width, best_lambdas)

  DE_min <- DEoptim::DEoptim(
    fn = calc_RSS,
    lower = lower_lambda,
    upper = upper_lambda,
    control = DEoptim::DEoptim.control(NP = K * 15,
                                       itermax = 100,
                                       c = 0.2,
                                       reltol = 1e-4,
                                       steptol = 10,
                                       trace = FALSE))

  #Get results
  best_lambdas <- DE_min$optim$bestmem
  best_RSS <- DE_min$optim$bestval

  cat("\nK =", K,"\t(nfeval", DE_min$optim$nfeval, ", iter", DE_min$optim$iter,")",
      "\nDEoptim:    \t", paste(round(best_lambdas, digits = 6), collapse = "\t"),
      "\tRSS =", formatC(best_RSS), "\n")

  components <- decompose_OSLcurve(C, best_lambdas, verbose = FALSE)
  plot_OSLcurve(C, components)
  weights <- simulate_OSLcurve(components,
                               template.curve = C,
                               simulate.curve = TRUE,
                               add.poisson.noise = FALSE)

  weights <- M / weights$signal
  #weights <- rep(1, length(weights$signal))


  # Calculate lambda error
 # for (i in 1:length(best_lambdas)) {
#
#    lambda_deviation <- 0.1
#    Chi2 <- calc_Chi2(components)

#    deviated_components <- components
#    deviated_components$lambda[i] <-  deviated_components$lambda[i] * (1 - lambda_deviation)
#    Chi1 <- calc_Chi2(deviated_components)

#    deviated_components <- components
#    deviated_components$lambda[i] <-  deviated_components$lambda[i] * (1 + lambda_deviation)
#    Chi3 <- calc_Chi2(deviated_components)

#    lambda_error <- sqrt((components$lambda[i] * lambda_deviation)^2 /
#                           ((Chi1 + Chi3)/2 - Chi2))

#    if(verbose) cat("Component", i,
#                    " Chi²'s: ", formatC(Chi1), " > ", formatC(Chi2), " < ", formatC(Chi3),
#                    " --> rel. error = ", formatC(lambda_error / components$lambda[i]), " \n")}

  # what are the numOSL results? ###############################################

  fit <- numOSL::decomp(cbind(C$time, C$signal),
                ncomp = K,
                plot = FALSE,
                weight = FALSE)

  if (fit$message == 0){

    lambda <- fit$LMpars[,3][1:K]
    cat("numOSL:     \t", paste(round(lambda, digits = 6), collapse = "\t"), "\tRSS =", formatC(fit$value), "\n")

  } else {
    cat("!LM-fit failed!\n")
  }


  # and now try nls and compare ###############################################

  ### Create fit formula ###
  n.names <- paste0("n.",1:K)
  lambda.names <- paste0("lambda.",1:K)

  decays <- paste(n.names," * (exp(-", lambda.names," * (time - ", channel_width,")) - exp(-", lambda.names," * time))"
                  , collapse=" + ")

  fit.formula <- as.formula(paste0("signal ~ ", decays))

  n <- components$n
  names(n) <- n.names

  lambda <- components$lambda
  names(lambda) <- lambda.names

  ### try Gauss-Newton fit ###
  fit2 <- try(minpack.lm::nlsLM(fit.formula, #nls
                 data = C,
                 start = c(n, lambda)), # weights = weights
             silent = FALSE)

  if (attr(fit2, "class") == "try-error") {

    cat("!nls-fit failed!\n")
  } else {

    #summary(fit)$parameters[paste0("lambda.", 1:4),2]
    lambda_results <- summary(fit2)$parameters[paste0("lambda.", 1:K),"Estimate"]
    lambda_error <- summary(fit2)$parameters[paste0("lambda.", 1:K),"Std. Error"]
    cat("nls:        \t", paste(round(lambda_results, digits = 6), collapse = "\t"),
        "\tChi² =", formatC(fit2$m$deviance()), "\n")
    cat("rel. errors:\t", paste(round(lambda_error / lambda, digits = 6) , collapse = "\t"), "\n")
  }

  # minpack.lm::nlsLM  ########################################################

  ### try Gauss-Newton fit ###
  fit3 <- try(minpack.lm::nlsLM(fit.formula,
                  data = C,
                  start = c(n, lambda),
                  weights = weights),
              silent = FALSE)

  if (attr(fit3, "class") == "try-error") {

    cat("!nlsLM-fit failed!\n")
  } else {

    #summary(fit)$parameters[paste0("lambda.", 1:4),2]
    lambda_results <- summary(fit3)$parameters[paste0("lambda.", 1:K),"Estimate"]
    lambda_error <- summary(fit3)$parameters[paste0("lambda.", 1:K),"Std. Error"]
    cat("nlsLM:      \t", paste(round(lambda_results, digits = 6), collapse = "\t"),
        "\tChi² =", formatC(fit3$m$deviance()), "\n")
    cat("rel. errors:\t", paste(round(lambda_error / lambda, digits = 6) , collapse = "\t"), "\n")
  }


}






# print computing time
cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")

