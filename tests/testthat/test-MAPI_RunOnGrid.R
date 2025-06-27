library("testthat")
library("mapi")
test_that("MAPI_RunOnGrid", {
  data("samples")
  data("metric")
  grid <- MAPI_GridHexagonal(samples, crs=3857, hw=250)
  results <- MAPI_RunOnGrid(samples, metric, grid=grid, nbPermuts=0, nbCores=1)
  expect_equal(sum(results$avg_value), 103.419112384369)
  # ignore.weights=TRUE
  results2 <- MAPI_RunOnGrid(samples, metric, grid=grid, nbPermuts=0, nbCores=1, ignore.weights=TRUE)
  expect_true(all(results2$nb_ell==results2$sum_wgts))
  expect_equal(sum(results2$avg_value), 204.278229490566)
})
