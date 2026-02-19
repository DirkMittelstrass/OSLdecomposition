test_that("check functionality", {
  testthat::skip_on_cran()

  #library(Luminescence)

  expect_output(expect_false(check_RLum.Data(matrix())),
                "Object is not of class 'RLum.Data.Curve' but of class 'matrix'")
  empty <- set_RLum("RLum.Data.Curve")
  expect_output(expect_false(check_RLum.Data(empty)),
                "Record type of object is not set")

  empty@recordType <- "_OSL"
  expect_output(expect_false(check_RLum.Data(empty)),
                "Record consists only of XSYG metadata")

  empty@recordType <- "OSL"
  empty@data <- matrix(1:10, ncol = 1)
  expect_output(expect_false(check_RLum.Data(empty)),
                "Curve data is no XY data")

  empty@recordType <- "OSL (PMT)"
  empty@data <- matrix(1:10, ncol = 2)
  expect_output(expect_true(check_RLum.Data(empty)))

  expect_output(expect_false(check_RLum.Data(empty, record_type = "IRSL")))
  expect_output(expect_false(check_RLum.Data(empty, record_type = "OSL2")))

  template <- empty
  expect_output(expect_true(check_RLum.Data(empty, curve_template = template)))

  template@data <- matrix(1:20, ncol = 1)
  expect_output(expect_error(check_RLum.Data(empty, curve_template = template),
                             "Invalid value of argument 'curve_template'"),
                "Template curve is invalid: Curve data is no XY data")

  template@data <- matrix(1:20, ncol = 2)
  expect_output(expect_false(check_RLum.Data(empty, curve_template = template)),
                "Number of data points differ between record and template")

  template@data <- matrix(1:10, ncol = 2)
  template@data[,1] <- 0:4
  expect_output(expect_false(check_RLum.Data(empty, curve_template = template)),
                "X-axes do not match between record and template")

  template@data[,1] <- 1:5
  template@info$TEMPERATURE <- 125
  empty@info$TEMPERATURE <- 125
  expect_output(expect_true(check_RLum.Data(empty, curve_template = template)))

  empty@info$TEMPERATURE <- 70
  expect_output(expect_false(check_RLum.Data(empty, curve_template = template)),
                "Value of parameter 'TEMPERATURE' does not match")
})
