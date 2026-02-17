#' @title Check if an object is a valid RLum.Data.Curve record for use in RLum.OSL functions
#'
#' @description
#' The input object is checked for the following properties:
#' * Is the object of class [Luminescence::RLum.Data.Curve-class] ?
#' * Has the object a slot `recordType`?
#' * Does the objects record type match with this functions argument `record_type`?
#' * Is the record not just a XSYG metadata object (marked by '_' before the record type name)?
#' * Contains the record a curve, thus has the object a slot `data`?
#' * Is the curve of type XY, thus has it a 2 x N dimension?
#' * If a `curve_template` is provided, the input object is also checked if it matches number of data points, x-axis and the `info` parameters "LPOWER", "LTYPE" and "TEMPERATURE".
#'
#' If all checks are positive, the input object is regarded as suitable for the functions
#' [RLum.OSL_correction], [RLum.OSL_global_fitting], [RLum.OSL_decomposition] and other functions if their
#' input curve is of type [Luminescence::RLum.Data.Curve-class].
#'
#' @param object [Luminescence::RLum.Data.Curve-class] (**required**):
#' Input object which shall be tested.
#'
#' @param record_type [character] (*with default*):
#' Expected type of record of the input `object`,
#' for example: `"OSL"`,`"SGOSL"` or `"IRSL"`.
#'
#' @param curve_template [Luminescence::RLum.Data.Curve-class] (*optional*):
#' Curve to check x-axis and some measurement parameter against.
#'
#' @param verbose [logical] (*with default*):
#' Enables console text output.
#'
#' @returns
#' A bolean value: `TRUE` or `FALSE`.
#'
#' @section Last updates:
#'
#' 2026-02-17, DM: Created function.
#'
#' @author
#' Dirk Mittelstraß, \email{dirk.mittelstrass@@luminescence.de}
#'
#' @examples
#'
#' # Load example data
#' data_path <- system.file("examples", "FB_10Gy_SAR.bin", package = "OSLdecomposition")
#' data_set <- Luminescence::read_BIN2R(data_path, fastForward = TRUE)
#'
#' # Test if record is of type OSL
#' check_RLum.Data(data_set[[5]]@records[[1]])
#'
#' @md
#' @export
check_RLum.Data <- function(
    object,
    record_type = "OSL",
    curve_template = NULL,
    verbose = TRUE
){

  obj_class <- class(object)[1]

  if (obj_class != "RLum.Data.Curve") {
    if (verbose) cat("Object is not of class 'RLum.Data.Curve' but of class",
                     paste0("'", obj_class, "'.\n"))
    return(FALSE)}

  obj_slots <- methods::slotNames(object)

  if (!("recordType" %in% obj_slots)) {
    if (verbose) cat("Object does not contain a slot 'recordType'.\n")
    return(FALSE)}

  obj_type <- object@recordType

  if (!grepl(paste0(record_type, "\\s\\([a-zA-Z]+\\)"), obj_type, fixed = FALSE)) {
    if (verbose) cat("Record is not of type", paste0("'", record_type, "'"),
                     "but of type", paste0("'", obj_type, "'."), "\n")
    return(FALSE)}

  if (startsWith(obj_type, "_")) {
    if (verbose) cat("Record consists only of XSYG metadata.\n")
    return(FALSE)}

  if (!("data" %in% obj_slots)) {
    if (verbose) cat("Record does not contain a slot 'data'.\n")
    return(FALSE)}

  if (length(dim(object@data)) != 2 || ncol(object@data) != 2) {
    if (verbose) cat("Curve data is no XY data.\n")
    return(FALSE)}

  # If everything is right so far, compare the record with the template
  if (!is.null(curve_template)) {

    # Test also the curve template
    if (!check_RLum.Data(curve_template, verbose = FALSE)) {
      cat("Template curve is invalid: ")
      check_RLum.Data(curve_template, verbose = TRUE)
      stop("Invalid value of argument 'curve_template'.")}

    if (nrow(curve_template@data) != nrow(object@data)) {
      if (verbose) cat("Number of data points differ between record and template.\n")
      return(FALSE)}

    # Check if the x-axes match. Allow small deviations (0.1 %)
    if (all(abs(curve_template@data[,1] - object@data[,1]) / abs(object@data[,1]) <= 0.001)) {
      if (verbose) cat("X-axes do not match between record and template.\n")
      return(FALSE)}

    # Check if measurement settings match
    info_params <- c("LPOWER", "LTYPE", "TEMPERATURE")
    for (param in info_params) {
      if ("info" %in% obj_slots &&
          "info" %in% methods::slotNames(curve_template) &&
          param %in% names(object@info) &&
          param %in% names(curve_template@info) &&
          object@info[[param]] != curve_template@info[[param]]) {

        if (verbose) cat("Value of parameter", paste0("'", param, "'"),
                         "does not match between record and template.\n")
        return(FALSE)}
    }
  }

  if (verbose) cat("Object is a valid record.\n")
  return(TRUE)
}
