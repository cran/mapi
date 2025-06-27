library("testthat")
library("mapi")
test_that("MAPI_EstimateHalfwidth", {
    data("samples")
    hw50 <- MAPI_EstimateHalfwidth(samples, crs=3857, beta=0.5)
    expect_equal(hw50, 286.835888252914)
    hw25 <- MAPI_EstimateHalfwidth(samples, crs=3857, beta=0.25)
    expect_equal(hw25, 143.417944126457)
})
