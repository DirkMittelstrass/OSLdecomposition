#' Creates result table of the SAR analysis
#'
#'
#' @section Changelog:
#' * 2019-10-14, DM: first version
#
#' @section ToDo:
#' * nice output with ggtexttable or formattable or kable
#' * refactorize, streamline and comment whole funciton
#' * Roxygen documentation
#'
#' @section Last changed: 2019-10-24
#'
#' @author
#' Dirk Mittelstrass, TU Dresden (Germany), \email{dirk.mittelstrass@@luminescence.de}
#'

table_SARdata <- function(
  SAR.list,
  data.type = "De",
  criterion = "median", # "mean", "CAM", "MAM"
  apply.RC = FALSE,
  age_rate = NA,
  unit = "",
  digits = 3,
  add.number = FALSE,
  add.number.prefix = "(", # "of"
  expected_value = NA,
  output.plot = FALSE,
  title = NULL
){


  if(is.na(age_rate)) {
    age_rate <- 1
    if (unit == "ka") unit <- "Gy"
  }
  x.length <- length(SAR.list)
  y.length <- max(lengths(SAR.list)) - 5
  Table <- data.frame(matrix(NA, nrow = y.length, ncol = x.length))
  Plot.data <- data.frame(NULL)

  for (x in 1:x.length) {
    for (y in 1:y.length) {

      Table[y,x] <- "-"

      if (class(SAR.list[[x]][[y]]) == "RLum.Results") {


        Data <- SAR.list[[x]][[y]]@data$data[[data.type]]


        if(data.type != "RC.Status"){

          RC.ok <- TRUE
          if (apply.RC) {

            RC <- SAR.list[[x]][[y]]@data$data[["RC.Status"]]
            RC.ok <- (RC == "OK")
            Data <- Data[RC.ok]
          }
        }
        Data <- Data[is.finite(Data)]

        added_value <- length(Data)

        if (length(Data) > 0) { ############################

        ### CRITERIONS ###

        if (criterion == "number") {

          Table[y,x] <- length(Data)

          if (data.type == "RC.Status") {


            Table[y,x] <- length(Data[(Data == "OK") &
                                        is.finite(SAR.list[[x]][[y]]@data$data[["De"]])])
            Data <- Data[is.finite(SAR.list[[x]][[y]]@data$data[["De"]])]
          }

        } else if (criterion == "median") {



          Table[y,x] <- round(median(Data), digits = digits)

        } else if (criterion == "mean") {

          Mean <- round(mean(Data), digits = digits)

          if (length(Data) > 1) {

            Sigma <- round(sd(Data), digits = digits)
            Table[y,x] <- paste0(Mean, " $\\pm$ ", Sigma)
          } else {

            Table[y,x] <- Mean
          }


          ### AGE MODELS ###

        } else if ((criterion == "CAM")||(criterion == "MAM")) {

          De <- SAR.list[[x]][[y]]@data$data$De
          De.Error <- SAR.list[[x]][[y]]@data$data$De.Error

          Selection <- is.finite(De) & is.finite(De.Error) & RC.ok

          Input <- data.frame(De = De[Selection],
                              De.Error = De.Error[Selection])

          Data <- De[Selection]


          if (criterion == "CAM") {

            if (length(Selection) > 1) {
              dose <- try(calc_CentralDose(Input,
                                           plot = FALSE,
                                           verbose = FALSE),
                          silent = TRUE)
            } else {

              Table[y,x] <- "-"
            }

          } else {

            dose <- try(calc_MinDose(Input,
                                 plot = FALSE,
                                 verbose = FALSE,
                                 sigmab = 0.2),
                        silent = TRUE)
          }

          if (attr(dose, "class") == "try-error") {

            Table[y,x] <- "-"
          } else if (is.na(dose@data[["summary"]][["de"]]) ||
                     is.na(dose@data[["summary"]][["de_err"]])) {

            Table[y,x] <- "-"
          } else if (length(Selection) > 1){

            return
            age <- round(dose@data[["summary"]][["de"]] * age_rate, digits = digits)
            age.error <- round(dose@data[["summary"]][["de_err"]] * age_rate, digits = digits)

            if (criterion == "CAM") added_value <- round(dose@data[["summary"]][["OD"]]) / 100

            Table[y,x] <- paste0(age, " $\\pm$ ", age.error)
          } else {

            Table[y,x] <- "-"
          }
        }


        ### UNIT ###

        if(Table[y,x] != "-") Table[y,x] <- paste0(Table[y,x], " ", unit)

        ### ADD NUMBER ###

        if (add.number && (Table[y,x] != "-")) {

          if ((add.number.prefix == "(") || (add.number.prefix == "()")) {

            Table[y,x] <- paste0(Table[y,x], " (", added_value, ")")
          } else {

            Table[y,x] <- paste0(Table[y,x],  add.number.prefix, " ", added_value)
          }
        }

        ### PLOT DATA ###

        if ((Table[y,x] != "-") && (length(Data) > 0)) {
          Box.style <- factor("classic")
          if (x > 2) Box.style <- factor(paste0("component ", y))
          Plot.data <- rbind(Plot.data,
                             data.frame(X = SAR.list[[x]]$table.header,
                                        Y = Data,
                                        G = Box.style))
        }

        #################################
        } else { Table[y,x] <- "-" }

      } else if (y <= length(SAR.list[[x]]$name)) {


        Table[y,x] <- "-"
        if (SAR.list[[x]]$name[y] == "Background") Table[y,x] <- ""

      } else {

        Table[y,x] <- ""
      }
    } # end for y

   colnames(Table)[x] <- SAR.list[[x]]$table.header
  } # end for x

  ### Add component indicies ###

  First.col <- data.frame(comp. = c(1:nrow(Table)))
  Table <- cbind(First.col, Table)

############################################

  if (output.plot) {

    library(gridExtra)
    library(scales)
    library(ggplot2)
    #theme_set(theme_bw())
    theme_set(theme_classic())

   # return(Plot.data)
    #Upper.limit <- median(Plot.data$Y, na.rm = TRUE) * 2
    Upper.limit <- quantile(Plot.data$Y, 0.9, na.rm = TRUE) * 1.25

    if (!is.na(expected_value) && !is.null(expected_value)) {

      if (1.2 * expected_value > Upper.limit) Upper.limit <- 1.2 * expected_value
    }


    p <- ggplot(Plot.data, aes(x=X, y=Y, fill=G)) +
      geom_boxplot() +
      scale_fill_manual(values=c("lightgrey", "red3","green3","blue3","darkorchid","gold","brown","pink")) +
      scale_y_continuous(limits = c(0, Upper.limit)) +
      ylab(paste0("Equivalent dose (", unit, ")")) + xlab("Signal calculation method") +
      labs(fill='source signals')

    if (!is.na(expected_value) && !is.null(expected_value)) {
      p <- p + geom_hline(yintercept = expected_value,
                          colour="grey", linetype = "longdash")
    }

     # geom_dotplot(binaxis='y', stackdir='center') + position=position_dodge(0)

    grid.arrange(p, nrow = 1, top = title)

  }


  return(Table)
}
