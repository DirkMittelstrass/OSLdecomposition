#' Analyse SAR CW-OSL measurements using OSL decomposition
#'
#' Calculates separate growth curves for each OSL signal component.
#' Alternative for [Luminescence::analyse_SAR.CWOSL]
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' data set of one or multiple aliquots containing containing SAR CW-OSL measurements
#'
#' @param aliquot_selection [c] (*optional*):
#' vector specifying which items (aliquots) of a list of [RLum.Analysis-class] objects shall be combined
#'
#' @param components [data.frame] (**required**):
#' Table containing the decay parameters of the OSL curve. One column must be named *$lambda*.
#' It is recommended to provide also the integration interval parameters (columns *t.start, t.end, ch.start, ch.end*),
#' which can be found by applying [calc_OSLintervals] on the global mean curve, calculated by [sum_OSLcurves].
#' If one or more column is missing, [calc_OSLintervals] is run automatically.
#'
#' @param record_type [character] (*with default*):
#' record type ("OSL","SGOSL","IRSL") which shall be evaluated
#'
#' @param recuperation_rate [numeric] (*with default*):
#' check aliquot data for the recuperation rate as rejection criteria; set to 'NA' for ignoring
#'
#' @param recycling_ratio [numeric] (*with default*):
#' check aliquot data for the recycling ratio as rejection criteria; set to 'NA' for ignoring
#'
#' @param dose.points [c] (*optional*):
#' a numeric vector containg the dose points values Using this argument
#' overwrites dose point values in the signal curves.
#'
#' @param negative_values_to_zero [logical] (*with default*): if TRUE negative n values are overwritten and set to n = 0;
#' negative n values are mathematically correct but physically debatable and can cause problems with [Luminescence::plot_GrowthCurve]
#'
#' @param data_has_Tx [logical] (*with default*):
#' if TRUE it is assumed that the sequence have a test dose steps to normalize the regeneration dose points.
#' If FALSE, every dose point is used as regeneration dose point
#'
#' @return
#' A [list] is returned containing the following elements:
#'
#' \item{[[i]]}{**list entry** RLum object related to the fast OSL decay component}
#' \item{...$data}{[data.frame] contains De results from [plot_GrowthCurve]}
#' \item{...$LnLxTnTx.table}{[list]([data.frame]) LxTx tables for each [RLum.Analysis-class] object (aliquot)}
#' \item{...$rejection.criteria}{[list]([data.frame]) rejection criteria tables for each [RLum.Analysis-class] object (aliquot)}
#' \item{...$Formula}{[list]([character]) growth curve fit formula for each [RLum.Analysis-class] object (aliquot)}
#'
#' @section Changelog:
#' * 2018-05-18, DM: first running version
#' * 2018-05-27, DM: first version with good chance of suceeding (for Bayreuth meeting 2018-05-28)
#' * 2018-07-02, DM: try(plot_GrowthCurve(...)); more integrity checks for rejection critera calculation
#' * 2019-04-30, DM: deleted fit-parameter-evaluation stuff; works now for any number of components; renamed function
#' * 2019-10-11, DM: Tweaked for calc_classicOSLsignals() and Rmarkdown-script
#'
#' @section ToDo:
#' * streamline code
#' * add capability to include background fitting
#'
#' @section Last changed: 2019-10-27
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [...], [...], [...]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for dose evaluation within the R luminescence package
#'
#' @keywords OSL CW-OSL CWOSL decompostion components SAR
#'
#' @md
#' @export
#'
#' @examples
#'

decompose_SARdata <- function(
  object,
  aliquot_selection = NULL,
  components,
  record_type = "OSL",
  algorithm = "det+nls", #"nls"
  recuperation_rate = 0.05,
  recycling_ratio = 0.1,
  dose.rate = 1,
  dose.points = NULL,
  applied.pre.dose = 0,
  fit.method = "EXP",
  negative_values_to_zero = TRUE,
  data_has_Tx = TRUE,
  verbose = TRUE
){

  ########## Input checks ###########

  # prove if object is a list of aliquots or just a single aliquot
  if (is.list(object)) {

    # prove if aliquot selection is given. If not, take all aliquots of the data set
    if (is.null(aliquot_selection)
        || is.na(aliquot_selection)
        || length(aliquot_selection) > length(object)) {
      aliquot_selection <- c(1:length(object))
    }
  } else {

    object <- list(object)
    aliquot_selection <- c(1)
  }


  K <- nrow(components)
  if ((components$name[K] == "background") || (components$name[K] == "Background")) {
    K <- K - 1
  }

  if ("SAR.compatible" %in% colnames(components)) {

    K <- sum(components$SAR.compatible, na.rm = TRUE)
  }
  X <- c(1:K)

  ########## Create data structure ###########
  results <- list(NULL)

  for (k in X) {
    results[[k]] <- set_RLum(class = "RLum.Results",
                             data = list(data = data.frame(),
                                         LnLxTnTx.table = list(),
                                         rejection.criteria = list(),
                                         Formula = list()))
  }

  # TODO: creating parameters for background substraction

  # predefine some intern parameters
  n <- 0
  n_record <- 0



  ##============================================================================##
  # MAIN CODE
  ##============================================================================##

  for (j in aliquot_selection) {
    if (j < 1 || j > length(object)) {
      message("Warning: Item ", j," is not a part of the data set. Item skipped")
    } else {
      if (class(object[[j]]) != "RLum.Analysis") {
        message("Warning: Item ", j," is not of class RLum.Analysis. Item skipped")
      } else {

        ##============================================================================##
        # DECOMPOSE DATA
        ##============================================================================##
        if (verbose) writeLines(" ", sep = "\n")
        if (verbose) writeLines(c("Decompose aliquot: ", j), sep = "\t")
        if (verbose) writeLines(" ", sep = "\n")

        records <- object[[j]]@records
        is_it_Lx <- TRUE # marker if the current record is an Lx or Tx record

        LxTx_table <- list(NULL)
        for (k in X) {
          LxTx_table[[k]] <- data.frame(NULL)
        }

        LxTx_row <- data.frame(NULL)
        Xn_residuals <- rep(0, K)
        n_Lx <- 1

        record_index <- NULL

        ########## Perform decomposition ###########

        for (i in c(1:length(records))) {
          if (records[[i]]@recordType == record_type) {

            record_index <- c(record_index, i)

            if ((algorithm == "late") || (algorithm == "early")) {
              decomposition <- calc_classicOSLsignal(records[[i]]@data,
                                                      components,
                                                      algorithm = algorithm,
                                                      verbose = FALSE)
            } else {
              decomposition <- decompose_OSLcurve(records[[i]]@data,
                                                  components,
                                                  algorithm = algorithm,
                                                  background.fitting = FALSE,
                                                  verbose = FALSE)
            }

            Xn <-  decomposition$n[X] - Xn_residuals
            Xn_residuals <- decomposition$n.residual[X]
            Xn_error <-  decomposition$n.error[X]

            n_record <- n_record + 1

            ########## Build LxTx table ###########

            if (is_it_Lx) {

              ## get dose point
              if (is.null(dose.points)){


                dose <- records[[i]]@info[["IRR_TIME"]]

                dose_rate <- 1
                # is the dose rate user given?
                if (dose.rate != 1) {

                  dose_rate <- dose.rate
                } else {

                  dose_rate <- records[[i]]@info[["IRR_DOSERATE"]]
                }

                # does the data contain dose rate information?
                if (is.finite(dose_rate) && dose_rate > 0)  dose <- dose * dose_rate

              }else{
                dose <- dose.points[n_Lx]
              }
              n_Lx <- n_Lx + 1


              for (k in X) {
                LxTx_row <- data.frame(Dose = dose,
                                       LxTx = Xn[k],
                                       LxTx.Error = Xn_error[k],
                                       LnLx = Xn[k],
                                       LnLx.Error = Xn_error[k],
                                       LnLx.Residual = Xn_residuals[k],
                                       TnTx = as.numeric(!data_has_Tx), # if there is no Tx although there should be one, set Tx = 0
                                       TnTx.Error = 0,
                                       TnTx.Residual = 0,
                                       Test.dose = 0)
                LxTx_table[[k]] <- rbind(LxTx_table[[k]], LxTx_row)
              }


              # in case testdoses are provided ...
              if (data_has_Tx) is_it_Lx <- FALSE

            } else {

              ### Testdose step ###
              # overwrite default LxTx and TnTx values
              for (k in X) {
                r <- length(LxTx_table[[k]]$Dose)
                LxTx_table[[k]]$LxTx[r] <- LxTx_table[[k]]$LxTx[r] / Xn[k]

                LxTx.Error <- abs(LxTx_table[[k]]$LxTx[r]
                               * ((LxTx_table[[k]]$LnLx.Error[r] / LxTx_table[[k]]$LnLx[r])^2
                                  + (Xn_error[k] / Xn[k])^2)^0.5)
                if (is.nan(LxTx.Error)) LxTx.Error <- 0

                ###########

                # ToDo: Refactorize this:

                test_dose <- records[[i]]@info[["IRR_TIME"]]

                dose_rate <- 1
                # is the dose rate user given?
                if (dose.rate != 1) {

                  dose_rate <- dose.rate
                } else {

                  dose_rate <- records[[i]]@info[["IRR_DOSERATE"]]
                }

                # does the data contain dose rate information?
                if (is.finite(dose_rate) && dose_rate > 0)  test_dose <- test_dose * dose_rate


                #################

                LxTx_table[[k]]$LxTx.Error[r] <- LxTx.Error
                LxTx_table[[k]]$Test.dose[r] <- test_dose
                LxTx_table[[k]]$TnTx[r] <- Xn[k]
                LxTx_table[[k]]$TnTx.Error[r] <- Xn_error[k]
                LxTx_table[[k]]$TnTx.Residual[r] <- Xn_residuals[k]
              }
              is_it_Lx <- TRUE
            }


          } # record if #
        } # record loop #


        ########## Analyze LxTx tables ###########
        for (k in X) {

          if (data_has_Tx) {

            # remove row with negative or zero Tx values
            LxTx_table[[k]] <- LxTx_table[[k]][LxTx_table[[k]]$TnTx >= 1,]
          }
          dose_values <- LxTx_table[[k]]$Dose
          LxTx_values <- LxTx_table[[k]]$LxTx

          # Check for negative LxTx values
          if (negative_values_to_zero) {

            LxTx_values[LxTx_values < 0] <- 0
            LxTx_table[[k]]$LxTx <- LxTx_values
          }

          ########## Rejection criteria ###########
          RC_status <- "FAILED"
          RC_table <- data.frame(recuperation.rate = NA,
                                 recuperation.rate.threshold = recuperation_rate,
                                 recycling.ratio = NA,
                                 recycling.ratio.threshold = recycling_ratio)

          try(
            if (!is.na(LxTx_values[1]) && (LxTx_values[1] > 0)) {
              RC_status <- "OK"

              # cycle through the LxTx table rows to find the recuperation and the recycling step
              for (l in c(2:length(LxTx_values))) {

                # recuperation rate
                if (dose_values[l] == 0) {
                  RC_table$recuperation.rate <- LxTx_values[l] / LxTx_values[1]
                }

                # recycling ratio
                for (m in c(3:length(LxTx_values))) {
                  if ((dose_values[l] == dose_values[m]) && (l != m) && (LxTx_values[m] > 0)) {
                    RC_table$recycling.ratio <- LxTx_values[l] / LxTx_values[m]
                    break
                  }
                }
              }

              if (!(is.na(RC_table$recuperation.rate) || is.na(recuperation_rate))){
                if(RC_table$recuperation.rate > recuperation_rate) RC_status <- "FAILED"}

              if (!(is.na(RC_table$recycling.ratio) || is.na(recycling_ratio))){
                if(abs(RC_table$recycling.ratio - 1) > recycling_ratio) RC_status <- "FAILED"}

              if (is.nan(RC_table$recuperation.rate) || is.nan(RC_table$recycling.ratio)) RC_status <- "FAILED"

            })


          ########## De calculation ###########

          # try to calculate De
          De_results <- try(plot_GrowthCurve(LxTx_table[[k]],
                                             fit.method = fit.method,
                                             output.plot = FALSE,
                                             verbose = verbose,
                                             fit.weights = TRUE),
                            silent = TRUE)

          # has it worked?
          if (class(De_results) == "RLum.Results") {

            # we are just interested in the results data.frame and the fitting formula
            De <- De_results@data$De
            Fit_formula <- De_results@data$Formula[[1]]

            if (is.finite(De$De)) {

              De$De <- De$De - applied.pre.dose
              if (De$De < 0) De$De <- NA
            }

            if (is.finite(De$De.MC)) {

              De$De.MC <- De$De.MC - applied.pre.dose
              if (De$De.MC < 0) De$De.MC <- NA
            }
          } else {

            # even if fitting failed, an empty data item is necessary to maintain data structure
            Fit_formula <- "fitting failed"
            De <- data.frame(De = NA,
                             De.Error = NA,
                             D01 = NA,
                             D01.ERROR = NA,
                             D02 = NA,
                             D02.ERROR = NA,
                             De.MC = NA,
                             Fit = "")

            RC_status <- "FAILED"

          }



          De <- tryCatch({

            De <- data.frame(De,
                             RC.Status = RC_status,
                             recuperation.rate = RC_table$recuperation.rate,
                             recycling.ratio = RC_table$recycling.ratio,
                             Ln = LxTx_table[[k]]$LnLx[1],
                             LnTn = LxTx_values[1])
          }, error = function(error_condition) {

            De <- data.frame(De,
                             RC.Status = "FAILED",
                             recuperation.rate = NA,
                             recycling.ratio = NA,
                             Ln = NA,
                             LnTn = NA)
          })

          # add Rejection criteria status

          if (class(De) == "data.frame") {

            # fill the RLum.object with data
            results[[k]]@data$data <- rbind(results[[k]]@data$data, De)
            results[[k]]@data$Formula <- c(results[[k]]@data$Formula, Fit_formula)
            results[[k]]@data$LnLxTnTx.table <- c(results[[k]]@data$LnLxTnTx.table, list(LxTx_table[[k]]))
            results[[k]]@data$rejection.criteria <- c(results[[k]]@data$rejection.criteria, list(RC_table))
            #results[[k]]$name <- components$name[k]
          } else {

            warning("Result data.frame couldn'be created")
          }
        }

        # count the analyzed aliquots
        n <- n + 1

      } # aliquot if 2 #
    } # aliquot if 1#
  } # aliquot loop #

  #################################################

  if (verbose) writeLines(paste0("Calculated equivalent doses of ", n, " aliquots with ", n_record ," OSL records"))

  results$index <- record_index
  results$name <- components$name
  results$lambda <- components$lambda
  if ((algorithm == "late") || (algorithm == "early")) {
    if (algorithm == "late") {

      results$table.header <- "late b."
    } else {

      results$table.header <- "early b."
    }

    #results$table.header <- paste0(algorithm, " background")
    results$comment <- paste0("Calculated with ", algorithm, " light background subtraction")
  } else {

    results$table.header <- paste0("K = ", length(components$lambda))
    results$comment <- paste0("Calculated with K = ", length(components$lambda), " decompositon")
  }
  if (verbose) writeLines(result$comment)
  return(results)
}
