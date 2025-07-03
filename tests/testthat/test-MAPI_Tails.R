library("parallel")
library("testthat")
library("mapi")
test_that("MAPI_RunOnGrid", {
  data("samples")
  data("metric")
  grid <- MAPI_GridHexagonal(samples, crs=3857, hw=250)
  set.seed(1234)  # for reproductibility
  results <- MAPI_RunOnGrid(samples, metric, grid=grid, nbPermuts=50, nbCores=1) # serial computation for reproductibility!
  expect_equal(sum(results$proba), 286.34)
  tails <- MAPI_Tails(results)
  expect_equal(nrow(tails[tails$tail=='upper',]), 1)
  expect_equal(nrow(tails[tails$tail=='lower',]), 3)
})
