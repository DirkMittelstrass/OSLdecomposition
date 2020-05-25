#' Various tools for checking and manipulation CW-OSL curves in RLum.Analysis data sets
#'
#' @last_changed 2020-05-24
#'
#' @param object [RLum.Analysis-class] (**required**):
#' input object containing data for analysis, alternatively a [list] of
#' [RLum.Analysis-class] objects can be provided.
#'
#' @param record_type
#' @param background_sequence
#' @param subtract_offset
#' @param check_consistency
#' @param check_signal_level **Inactive**
#' @param tailor_records **Inactive**
#' @param record_cut
#' @param correct_PMTsaturation **Inactive**
#' @param verbose
#' @param report
#'
#'
#' @Notes
#' The order of processing is:
#' check_consistency -> tailor_records -> record_cut -> correct_PMTsaturation -> background_sequence -> check_signal_level -> subtract_offset
#'
#' @return
#' @examples

RLum.OSL_correction <- function(
  object,
  record_type = "OSL",
  background_sequence = NULL,
  subtract_offset = 0,
  check_consistency = TRUE,
  check_signal_level = FALSE,
  tailor_records = FALSE,
  record_cut = 20,
  correct_PMTsaturation = NA,
  verbose = TRUE,
  report = FALSE
){
  ### ToDo's
  # - Check for Zero as first value at the time axis
  # - enhance 'record_sorting' to accept vectors of @info-arguments, include LPOWER and LIGHTSOURCE per default and print arguments
  # - enhance 'background' to accept whole RLum objects
  # - deploy Luminescence::verify_SingleGrainData() for 'check_single_grain_signal'
  # - handle previous CORRECTION steps
  # - new 'reduce data' argument to delete all unnecessary data (like TL curves etc.)

  library(Luminescence)

  object_name <- deparse(substitute(object))

  # define new list object to safely ignore incompatible list elements
  data_set <- list()
  data_set_overhang <- list()

  # Build overview list
  cor_data <- list(parameters = list(record_type = record_type,
                                     check_consistency = check_consistency,
                                     tailor_records = tailor_records,
                                     background_sequence = background_sequence,
                                     check_signal_level = check_signal_level,
                                     subtract_offset = subtract_offset))

  # Test if object is a list. If not, create a list
  if (is.list(object)) {

    for (i in 1:length(object)) {

      if (class(object[[i]]) == "RLum.Analysis") {

        data_set[[length(data_set) + 1]] <- object[[i]]
      } else {

        element_name <- names(object)[i]
        if (element_name=="CORRECTION") {
          warning("Data set was already manipulated by [RLum.OSL_correction()]. Old information in $CORRECTION were overwritten")

          } else {

            data_set_overhang[[element_name]] <- object[[i]]
            warning("List element ", i, " is not of type 'RLum.Analysis' and was not included in fitting procedure")}}}


  } else {

    data_set <- list(object)
    warning("Input was not of type list, but output is of type list")}

  if (length(data_set) == 0) stop("Input object contains no RLum.Analysis data")

  if (!(check_consistency) && !(is.null(background_sequence))) {
    stop("Background correction requires consistent data! Please set 'check_consistency=TRUE' and try again.")}

  correction_step <- 0

  if (check_consistency) {
    correction_step <- correction_step + 1 #########################################################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Check records for consistency in the detection settings -----\n")

    # measure computing time
    time.start <- Sys.time()

    # Characteristics vector. Will be extentable later
    Cvector <- c("Channels", "Channel width")

    ### Build table ###
    Ctable <- data.frame(NULL)
    for (j in 1:length(data_set)) {
      for (i in 1:length(data_set[[j]]@records)) {
        if (data_set[[j]]@records[[i]]@recordType == record_type) {

          channels <- length(data_set[[j]]@records[[i]]@data[,1])
          channel_width <- data_set[[j]]@records[[i]]@data[,1][2] - data_set[[j]]@records[[i]]@data[,1][1]

          new_line <- data.frame(sequence = j,
                        record = i,
                        record_type = as.character(record_type),
                        channels = channels,
                        channel_width = channel_width)


          ####### Add code for @info parameters here

          # build parameter string
          new_line <- cbind(new_line,
                            data.frame(string = paste(new_line[4:length(new_line)], collapse = ", ")))

          # add measurement to overview table
          Ctable <- rbind(Ctable, new_line)
        }
      }
    } #### end loop

    # now get parameter statistics
    Cstats <- as.data.frame(table(Ctable$string))

    N <- nrow(Cstats)
    if (N > 1) {

      # sort rows for frequency
      Cstats <- Cstats[order(-Cstats$Freq),]

      # Add record_type column and numerate them
      Cstats$record_type <- record_type
      Cstats$record_type[2:N] <- paste0(record_type, 2:N)

      # Insert new levels in OSL record type colection
      levels(Ctable$record_type) <- Cstats$record_type


      sequences_with_replacements <- c(NULL)
      for (i in 2:N) {
        # Replace record_types in the condition table
        Ctable$record_type[Ctable$string == Cstats$Var1[i]] <- Cstats$record_type[i]

        # Remember the measurement sequences
        sequences_with_replacements <- c(sequences_with_replacements,
                                         Ctable$sequence[Ctable$string == Cstats$Var1[i]])
      }
      # delete doublettes
      sequences_with_replacements <- unique(sequences_with_replacements)

      # Rename the recordType properties in the RLum.curve objects
      for (i in 1:nrow(Ctable)) {

        data_set[[Ctable$sequence[i]]]@records[[Ctable$record[i]]]@recordType <- as.character(Ctable$record_type[i])
      }

      # Display statistics
      colnames(Cstats) <- c("settings", "frequency", "record_type")
      Cvector
      if(verbose) cat(paste0("Frequency table of different sets of detection settings (",
                             paste(Cvector, collapse = ", "),"):\n"))
      if(verbose) print(Cstats)
      if(verbose) cat("RLum.Data.Curve@RecordType changed to",
                      paste(paste0(record_type, 2:N), collapse = " or "),
                      "in sequence:", paste(sequences_with_replacements, collapse = ", "), "\n")
      if(verbose) cat("Further data manipulations are performed just on", record_type,"records\n")

    } else {
      if(verbose) cat("All", record_type, "records have the same detection settings\n")
    }

    cor_data <- c(cor_data, list(measurement_characteristics = Ctable,
                            character_statistics = Cstats))


    # print needed computing time
    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
  }

  if (tailor_records) {
    correction_step <- correction_step + 1 ########################### TAILOR ##############################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Cut not-stimulated measurement parts -----\n")
    time.start <- Sys.time()


    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
  }


  if (is.numeric(record_cut) && (record_cut > 0)) {
    correction_step <- correction_step + 1 ########################## CUT ###############################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Cut records at specific time -----\n")
    time.start <- Sys.time()
    records_changed <- 0

    for (j in 1:length(data_set)) {
      for (i in c(1:length(data_set[[j]]@records))) {
        if (data_set[[j]]@records[[i]]@recordType == record_type) {

          time <- data_set[[j]]@records[[i]]@data[,1]
          signal <- data_set[[j]]@records[[i]]@data[,2]
          channel_width <- time[2] - time[1]

          ##### cut curve, if too long #####
          if (max(time) > record_cut) {

            # reduce channel number and write record back into the RLum object
            time <- time[1:ceiling(record_cut/channel_width)]
            signal <- signal[1:ceiling(record_cut/channel_width)]
            data_set[[j]]@records[[i]]@data <- matrix(c(time, signal), ncol = 2)
            records_changed <- records_changed + 1
          }
        }
      }
    }

    if (records_changed==0) {
      if(verbose) cat("No record was shorter than", record_cut, "s\n")
    } else{
      if(verbose) cat("Measurement duration of", records_changed, "records reduced to", record_cut, "s\n")}

    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
  }

  if (check_signal_level) {
    correction_step <- correction_step + 1 ######################### PMT SATURATION ################################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Correct for PMT saturation effects -----\n")
    time.start <- Sys.time()

    if(verbose) cat("THIS FUNCTION IS CURRENTLY NOT AVAILABLE")

    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
  }

  if (!(is.null(background_sequence) || is.na(background_sequence))
      && is.vector(background_sequence)
      && (all(background_sequence> 0) ) && (all(background_sequence<= length(data_set)))) {
    correction_step <- correction_step + 1 ######################### BACKGROUND ################################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Subtract background measurement -----\n")
    time.start <- Sys.time()


    N <- 0
    # create background curve
    background_curve <- sum_OSLcurves(data_set, record_type,
                             aliquot_selection = background_sequence,
                             output.plot = FALSE,
                             verbose = TRUE)

    # rename background OSL curves
    for (j in background_sequence) {
      for (i in c(1:length(data_set[[j]]@records))) {
        if (data_set[[j]]@records[[i]]@recordType == record_type) {

          N <- N + 1
          data_set[[j]]@records[[i]]@recordType <- paste0(record_type, "background")}}}

    #if(verbose) cat("Göobal background curve created from", N, record_type,"records of sequence \n")
    if(verbose) cat("RLum.Data.Curve@RecordType is changed to",
                    paste0(record_type, "background"),"for those records\n")

    # subtract background curve
    for (j in 1:length(data_set)) {
      for (i in c(1:length(data_set[[j]]@records))) {
        if (data_set[[j]]@records[[i]]@recordType == record_type) {

          time <- data_set[[j]]@records[[i]]@data[,1]
          signal <- data_set[[j]]@records[[i]]@data[,2]

          signal <- signal - background_curve$signal

          N <- N + 1
          data_set[[j]]@records[[i]]@data <- matrix(c(time, signal), ncol = 2)}}}

    cor_data <- c(cor_data, list(background_curve = background_curve))

    if(verbose) cat("Background saved at @CORRECTION@background_curve\n")
    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")

  } else {
    if (!(is.null(background_sequence))) {

      if(verbose) cat("Invalid 'background' argument. Step skipped\n")
      warning("Invalid 'background' argument")}}


  if (check_signal_level) {
    correction_step <- correction_step + 1 ######################### NOISE CHECK ################################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Check measurements for sufficient signal levels -----\n")
    time.start <- Sys.time()

    if(verbose) cat("THIS FUNCTION IS CURRENTLY NOT AVAILABLE")

    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
  }


  if (subtract_offset != 0) {
    correction_step <- correction_step + 1 ######################### OFFSET ################################
    if(verbose) cat("CORRECTION STEP", correction_step,"----- Subtract offset value -----\n")
    time.start <- Sys.time()

    for (j in 1:length(data_set)) {
      for (i in c(1:length(data_set[[j]]@records))) {
        if (data_set[[j]]@records[[i]]@recordType == record_type) {

          time <- data_set[[j]]@records[[i]]@data[,1]
          signal <- data_set[[j]]@records[[i]]@data[,2]
          channel_width <- time[2] - time[1]

          signal <- signal - subtract_offset * channel_width
          data_set[[j]]@records[[i]]@data <- matrix(c(time, signal), ncol = 2)
        }
      }
    }
    if(verbose) cat("Offset of", subtract_offset, "counts per second subtracted from every", record_type, "record\n")
    if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
  }



  ################################ REPORT  ################################

  if (FALSE) {
    if(("rmarkdown" %in% rownames(installed.packages())) && ("kableExtra" %in% rownames(installed.packages()))) {

      if(verbose) cat("CORRECTION SUMMARY ----- Create report -----\n")
      if(verbose) cat("This process can take a few minutes...\n")
      library(rmarkdown)

      time.start <- Sys.time()

      # the RMD script has to be located in the "/inst" folder of the project
      # then, it will be installed with the package
      try({

        report_format <- "html"
        # for test purposes:
        rmd_path <- "C:\\Users\\mitte\\Desktop\\R\\OSLdecomposition\\inst\\rmd\\report_Correction.Rmd"
        #rmd_path <- system.file("rmd", "report_Step1.Rmd", package = "OSLdecomposition")
        output_file <- paste0(getwd(), "/", "report_Correction.", report_format)

        rmarkdown::render(rmd_path,
                          params = list(dec_data = dec_data,
                                        data_set = data_set,
                                        object_name = object_name),
                          output_file = output_file,
                          output_format = paste0(report_format,"_document"),
                          quiet = TRUE)

        cat("Save", toupper(report_format), "report to:", output_file, "\n")

        # ToDo: Replace the following try() outside the big try
        try({
          browseURL(output_file)
          cat("Open", toupper(report_format), "report in the systems standard browser\n")
        })


        if(verbose) cat("(time needed:", round(as.numeric(difftime(Sys.time(), time.start, units = "s")), digits = 2),"s)\n\n")
      })

    } else {

      warning("Packages 'rmarkdown' and 'kableExtra' are needed to create reports. One or both are missing.")
    }
  }

  # Return decomposed data
  object <- c(data_set, data_set_overhang, CORRECTION = list(cor_data))
  invisible(object)


}
