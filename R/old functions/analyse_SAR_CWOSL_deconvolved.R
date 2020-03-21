#' Analyse SAR CW-OSL measurements using OSL deconvolution
#'
#' Calculates separate growth curves for each OSL signal component.
#' Alternative for [analyse_SAR.CWOSL]
#'
#' @param object [RLum.Analysis-class] or [list](RLum.Analysis) (**required**):
#' data set of one or multiple aliquots containing containing SAR CW-OSL measurements
#'
#' @param aliquot_selection [c] (*optional*):
#' vector specifying which items (aliquots) of a list of [RLum.Analysis-class] objects shall be combined
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
#' @param single_aliquot_fitting [logical] (*with default*):
#' ifTRUE the superposition curve and the decay parameter fitting is done for each [RLum.Analysis-class] object (e.g. aliquot) seperatly.
#' If FALSE (default), the whole data set is analyzed with the same superposition curve and the same OSL decay parameters
#'
#' @param data_has_Tx [logical] (*with default*):
#' if TRUE it is assumed that the sequence have a test dose steps to normalize the regeneration dose points.
#' If FALSE, every dose point is used as regeneration dose point
#'
#' @param offset_value [numeric] (*with default*):
#' signal offset (background) which will be substracted from each record; value will be transfered to
#' [calc_CWOSLcomponents] and [calc_superposition.curve]
#'
#' @param output.plot [logical] (*with default*):
#' returns a plot of the growth curves and the fitted superposition curve
#'
#' @param output.plot.extended [logical] (*with default*):
#' plots also every curve with drawn-in OSL components
#'
#' @return
#' A [list] is returned containing the following elements:
#'
#' \item{fast.component}{[RLum.Results-class] RLum object related to the fast OSL decay component}
#' \item{medium.component}{[RLum.Results-class] similar to fast.component; empty if just one or two components was found}
#' \item{slow.component}{[RLum.Results-class] similar to fast.component; empty if just one component was found}
#' \item{...$data}{[data.frame] contains De results from [plot_GrowthCurve]}
#' \item{...$LnLxTnTx.table}{[list]([data.frame]) LxTx tables for each [RLum.Analysis-class] object (aliquot)}
#' \item{...$rejection.criteria}{[list]([data.frame]) rejection criteria tables for each [RLum.Analysis-class] object (aliquot)}
#' \item{...$Formula}{[list]([character]) growth curve fit formula for each [RLum.Analysis-class] object (aliquot)}
#' \item{indices}{[c] vector of selected aliquot; similar to input parameter aliquot_selection but after integrity check; useful to associate result index with data set item}
#' \item{superposition.curve}{[list] of [RLum.Data.Curve-class], return from [calc_superposition.curve]}
#' \item{decay.parameter}{[data.frame] results from [fit_CWCurve] in a nutshell}
#'
#' @section Changelog:
#' * 2018-05-18, DM: first running version
#' * 2018-05-27, DM: first version with good chance of suceeding (for Bayreuth meeting 2018-05-28)
#' * 2018-07-02, DM: try(plot_GrowthCurve(...)); more integrity checks for rejection critera calculation
#'
#' @section ToDo:
#' * (!) More than 3 components possible
#' * (!) Add classic OSL signal calculation (Galbraith & Roberts 2012)
#' * overview table records
#' * RecordType other than OSL possible
#' * possibility to define one (or more) aliquots as background aliquot. Their aritmetic mean curve shall be used as background
#' * rejection criteria and equal to classic analyze_SAR
#' * improve performance by transmitting the release propability determinant from calc_intervals to calc_components
#' * inherit rejection criteria and other stuff from classic analyze_SAR
#' * component seperation plot of superposition curve, including pseudoLM-OSL plot
#'
#' @section Function version: 0.1.2
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#' @seealso [calc_CWOSLcomponents], [calc_deconvolution.intervals], [calc_superposition.curve]
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for dose evaluation within the R luminescence package
#'
#' @keywords CWOSL deconvolution components SAR
#'
#' @md
#' @export
#'
#' @examples
#'

analyse_SAR.CWOSL.deconvolved <- function(
  object,
  aliquot_selection = NULL,
  recuperation_rate = 0.1,
  recycling_ratio = 0.1,
  dose.points = NULL,
  single_aliquot_fitting = FALSE,
  data_has_Tx = TRUE,
  offset_value = 0,
  output.plot = TRUE,
  output.plot.extended = FALSE
){

  ##============================================================================##
  # INPUT OBJECTS & INTEGRITY CHECKS
  ##============================================================================##

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

  ###  creating data structure
  results <- list(late.background = set_RLum(class = "RLum.Results",
                                            data = list(data = data.frame(),
                                                        LnLxTnTx.table = list(),
                                                        rejection.criteria = list(),
                                                        Formula = list())),

                  fast.component = set_RLum(class = "RLum.Results",
                                            data = list(data = data.frame(),
                                                        LnLxTnTx.table = list(),
                                                        rejection.criteria = list(),
                                                        Formula = list())),
                  medium.component = set_RLum(class = "RLum.Results",
                                              data = list(data = data.frame(),
                                                          LnLxTnTx.table = list(),
                                                          rejection.criteria = list(),
                                                          Formula = list())),
                  slow.component = set_RLum(class = "RLum.Results",
                                            data = list(data = data.frame(),
                                                        LnLxTnTx.table = list(),
                                                        rejection.criteria = list(),
                                                        Formula = list())),
                  early.background = set_RLum(class = "RLum.Results",
                                              data = list(data = data.frame(),
                                                          LnLxTnTx.table = list(),
                                                          rejection.criteria = list(),
                                                          Formula = list())),
                  indices = aliquot_selection,
                  superposition.curve = list(),
                 # fit.result = list(),
                  decay.parameter = data.frame(NULL))

  ### creating parameters for background substraction

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
        # GET OSL DECAY PARAMETERS
        ##============================================================================##

        if (single_aliquot_fitting || (n == 0)) {

          # shall every aliquot be analyzed with its own OSL components?
          if (single_aliquot_fitting) {
            superpos_selection <- c(j)
          }else{
            superpos_selection <- aliquot_selection
          }

          # calc arithmetic mean curve
          superposition_curve <- calc_superposition.curve(object,
                                                          record_type = "OSL",
                                                          aliquot_selection = superpos_selection,
                                                          offset_value = offset_value,
                                                          output.plot = output.plot.extended)
          results$superposition.curve <- c(results$superposition.curve,
                                              superposition_curve)

          # find components via fitting
          fit_results <- fit_CWCurve(superposition_curve,
                                              n.components.max = 3,
                                              plot =  output.plot)
          f_fast <- fit_results@data[["data"]][["lambda1"]]
          f_medium <- fit_results@data[["data"]][["lambda2"]]
          f_slow <- fit_results@data[["data"]][["lambda3"]]

          # if just one component could be found
          if (is.na(f_medium) && is.na(f_slow)) message("Just one component found, defined as 'fast' component by default")

          # if no third component could be found, set second component as 'slow'
          if (!is.na(f_medium) && is.na(f_slow)) {
            f_slow <- f_medium
            f_medium <- NA
            message("Just two signal components found, defined as 'fast' and 'slow' component by default")
          }

          # get the deconvolution intervals
          intervals <- calc_deconvolution.intervals(f.fast = f_fast,
                                                    f.medium = f_medium,
                                                    f.slow = f_slow,
                                                    channel.width = superposition_curve[2,1] - superposition_curve[1,1],
                                                    channel.number = length(superposition_curve[,1]))

          results$decay.parameter <- rbind(results$decay.parameter,
                                          data.frame(f.fast = f_fast,
                                                     f.medium = f_medium,
                                                     f.slow = f_slow,
                                                     intervals))
        }

        ##============================================================================##
        # DECONVOLVE DATA
        ##============================================================================##
        writeLines(" ", sep = "\n")
        writeLines(c("Deconvolve aliquot: ", aliquot_selection[j]), sep = "\t")
        writeLines(" ", sep = "\n")

        records <- object[[j]]@records
        is_it_Lx <- TRUE # marker if the current record is an Lx or Tx record
        LxTx_table <- list(fast = data.frame(NULL),
                           medium = data.frame(NULL),
                           slow = data.frame(NULL))
        LxTx_row <- data.frame(NULL)
        Xn_residuals <- c(0,0,0)
        n_Lx <- 1

        for (i in c(1:length(records))) {
          if (records[[i]]@recordType == "OSL") {

            if (is_it_Lx) {
              plot_text <- "Lx_"
            } else{
              plot_text <- "Tx_"
            }
            plot_text <- paste(plot_text, n_Lx - 1, " of aliquot ", j)

            deconvolution <- calc_CWOSLcomponents(records[[i]]@data,
                                                  f.fast = f_fast,
                                                  f.medium = f_medium,
                                                  f.slow = f_slow,
                                                  t0 = intervals$t0,
                                                  t1 = intervals$t1,
                                                  t2 = intervals$t2,
                                                  t3 = intervals$t3,
                                                  offset_value = offset_value,
                                                  output.plot = output.plot.extended,
                                                  info = TRUE,
                                                  main = plot_text)

            Xn <-  deconvolution$n - Xn_residuals
            Xn_error <-  deconvolution$std.dev
            Xn_residuals <- deconvolution$n.residual
            n_record <- n_record + 1


            if (is_it_Lx) {

              ## get dose point
              if (is.null(dose.points)){

                dose <- records[[i]]@info[["IRR_TIME"]]
                dose_rate <- records[[i]]@info[["IRR_DOSERATE"]]
                if (is.finite(dose_rate) && dose_rate > 0) {
                  # data contains dose rate informations
                  dose <- dose * dose_rate
                }
              }else{
                dose <- dose.points[n_Lx]
              }
              n_Lx <- n_Lx + 1


              for (k in c(1:3)) {
                LxTx_row <- data.frame(Dose = dose,
                                       LxTx = Xn[k],
                                       LxTx.Error = Xn_error[k],
                                       LnLx = Xn[k],
                                       LnLx.Error = Xn_error[k],
                                       LnLx.Residual = Xn_residuals[k],
                                       TnTx = 1,
                                       TnTx.Error = 0,
                                       TnTx.Residual = 0)
                LxTx_table[[k]] <- rbind(LxTx_table[[k]], LxTx_row)
              }


              # in case testdoses are provided ...
              if (data_has_Tx) is_it_Lx <- FALSE

            } else {

              ### Testdose step ###
              # overwrite default LxTx and TnTx values
              for (k in c(1:3)) {
                r <- length(LxTx_table[[k]]$Dose)
                LxTx_table[[k]]$LxTx[r] <- LxTx_table[[k]]$LxTx[r] / Xn[k]
                LxTx_table[[k]]$LxTx.Error[r] <- (LxTx_table[[k]]$LxTx[r]
                                                  * ((LxTx_table[[k]]$LnLx.Error[r] / LxTx_table[[k]]$LnLx[r])^2
                                                     + (Xn_error[k] / Xn[k])^2)^0.5)
                LxTx_table[[k]]$TnTx[r] <- Xn[k]
                LxTx_table[[k]]$TnTx.Error[r] <- Xn_error[k]
                LxTx_table[[k]]$TnTx.Residual[r] <- Xn_residuals[k]
              }
              is_it_Lx <- TRUE
            }


          } # record if #
        } # record loop #


        # make the next step seperatly for each component
        for (k in 1:3) {

          ##============================================================================##
          # REJECTION CRITERIA
          ##============================================================================##

          dose_values <- LxTx_table[[k]]$Dose
          LxTx_values <- LxTx_table[[k]]$LxTx
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
          })





          ##============================================================================##
          # RESULTS
          ##============================================================================##


         # if(is.na(LxTx_table[[k]]$LxTx[1])){
         #   results[[k]]@data$data <- rbind(results[[k]]@data$data, NA)
         #   results[[k]]@data$Formula <- c(results[[k]]@data$Formula, "")
         # } else {

            # try to calculate De
            De_results <- try(plot_GrowthCurve(LxTx_table[[k]],
                                               fit.method = "EXP OR LIN",
                                               output.plot = output.plot),
                              silent = TRUE)

            # has it worked?
            if (class(De_results) == "RLum.Results") {

              # we are just interested in the results data.frame and the fitting formula
              De <- De_results@data$De
              Fit_formula <- De_results@data$Formula[[1]]

              # Check if De is infinite
              if (is.infinite(De$De)) De$De <- NaN

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
            }

            # add Rejection criteria status
            De <- data.frame(De,
                             RC.Status = factor(RC_status),
                             recuperation.rate = RC_table$recuperation.rate,
                             recycling.ratio = RC_table$recycling.ratio)

            # fill the RLum.object with data
            results[[k]]@data$data <- rbind(results[[k]]@data$data, De)
            results[[k]]@data$Formula <- c(results[[k]]@data$Formula, Fit_formula)
            results[[k]]@data$LnLxTnTx.table <- c(results[[k]]@data$LnLxTnTx.table, list(LxTx_table[[k]]))
            results[[k]]@data$rejection.criteria <- c(results[[k]]@data$rejection.criteria, list(RC_table))

                      }

        # count the analyzed aliquots
        n <- n + 1

      } # aliquot if 2 #
    } # aliquot if 1#
  } # aliquot loop #



  ##============================================================================##
  # RETURN VALUES
  ##============================================================================##

  message("Calculated component seperated equivalent doses of ", n, " aliquots with ", n_record ," OSL records")
  return(results)


}
