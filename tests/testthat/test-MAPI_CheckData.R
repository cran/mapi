library("testthat")
library("mapi")
test_that("MAPI_CheckData", {
  data("samples")
  data("metric")
  clean.data <- mapi::MAPI_CheckData(samples, metric)
  expect_equal(nrow(clean.data[['samples']]), 200)
  expect_equal(nrow(clean.data[['metric']]), 39800)
})
