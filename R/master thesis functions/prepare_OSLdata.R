#' CW-OSL preperation
#'
#' The function prepares raw CW-OSL data sets for OSL decomposition process.
#'
#'
#' @section Changelog:
#'
#' @section ToDo:
#' * Give back data set structure
#' * Correct x-axis, if it begins with Zero
#' * Delete preset aliquots
#' * Get vector with OSL record indicies
#' * test sequence on SAR compability for each aliquot
#'
#' @section Last changed: 2019-10-25
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'
#'
#' @references
#' Dirk Mittelstrass, Christoph Schmidt, Sebastian Kreutzer, ... .*in preperation*. Algebraic CW-OSL
#' signal component decomposition of quartz and its use for natural dose evaluation
#'
#' @keywords OSL CWOSL CW-OSL deconvolution decomposition component components
#'
#' @md
#' @export
#'
#' @examples
#'
prepare_OSLdata <- function(
  object,
  record.type = "OSL",
  cut.time = NA,
  background.aliquot = NA,
  background.plot = TRUE,
  tailor.single.grain = FALSE
){

  ########## Input checks ###########

  if(is.na(cut.time)) cut.time <- Inf

  # prove if object is a list of aliquots or just a single aliquot
  if (!is.list(object)) stop("Error: Object is not a list")

  #### Background sum ####
  if (!is.na(background.aliquot)) {

    # create background curve
    B.curve <- sum_OSLcurves(object, record.type,
                             aliquot_selection = background.aliquot,
                             output.plot = background.plot,
                             plot.first = FALSE,
                             plot.global = TRUE,
                             verbose = FALSE,
                             title = NULL)

    # delete background aliquots
    object[background.aliquot] <- NULL
  }

#### Single grain sum ####
  if (tailor.single.grain) {
    S.curve <- sum_OSLcurves(object, record.type,
                             output.plot = FALSE,
                             verbose = FALSE)
    time <- S.curve$time
    signal <- S.curve$signal
    d1 <- diff(S.curve$signal)
    d2 <- diff(d1)
    half <- floor(length(d2)/2)
    cut1 <- which.min(d2[1:half]) +  1
    cut2 <- which.min(d2[half:length(d2)])
  }


  for (j in 1:length(object)) {

    record_index <- NULL

    for (i in c(1:length(object[[j]]@records))) {
      if (object[[j]]@records[[i]]@recordType == record.type) {

        record_index <- c(record_index, i)

        time <- object[[j]]@records[[i]]@data[,1]
        signal <- object[[j]]@records[[i]]@data[,2]
        channel.width <- time[2] - time[1]

        ##### remove background #####
        if (!is.na(background.aliquot)) {

          signal <- signal - B.curve$signal
        }

        ##### tailor curve ####
        if (tailor.single.grain) {

          signal <- signal[cut1:cut2]
          time <- c(1:length(signal)) * channel.width
        }

        ##### cut curve, if too long #####
        if (max(time) > cut.time) {

          time <- time[1:ceiling(cut.time/channel.width)]
          signal <- signal[1:ceiling(cut.time/channel.width)]
        }

        object[[j]]@records[[i]]@data <- matrix(c(time, signal), ncol = 2)
      }

    }

  }

  return(object)
}
