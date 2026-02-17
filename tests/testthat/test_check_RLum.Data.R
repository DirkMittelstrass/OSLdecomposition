test_that("check functionality", {
  testthat::skip_on_cran()

  expect_output(expect_false(check_RLum.Data(matrix())),
                "Object is not of class 'RLum.Data.Curve' but of class 'matrix'")
  empty <- set_RLum("RLum.Data.Curve")
  expect_output(expect_false(check_RLum.Data(empty)),
                "Record is not of type 'OSL' but of type 'NA'")

  empty@recordType <- "_OSL"
  expect_output(expect_false(check_RLum.Data(empty)),
                "Record is not of type 'OSL' but of type '_OSL'")

  empty@recordType <- "OSL"
  empty@data <- matrix(1:10, ncol = 1)
  expect_output(expect_false(check_RLum.Data(empty)),
  ## should say "Curve data is no XY data"
                "Record is not of type 'OSL' but of type 'OSL'")
})
