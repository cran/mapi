library("parallel")
library("testthat")
library("mapi")
test_that("ignore.weights+MAPI_Plot2", {
  data("samples")
  data("metric")
  grid <- MAPI_GridHexagonal(samples, crs=3857, hw=250)
  set.seed(1234)  # for reproductibility
  results <- MAPI_RunOnGrid(samples, metric, grid=grid, nbPermuts=50, nbCores=1, ignore.weights=TRUE) # serial computation for reproductibility!
  expect_equal(sum(results$proba), 969.44)
  tails <- MAPI_Tails(results)
  expect_equal(nrow(tails[tails$tail=='upper',]), 2)
  expect_equal(nrow(tails[tails$tail=='lower',]), 1)
  expect_no_error(MAPI_Plot2(results, tails, samples=samples))
})
