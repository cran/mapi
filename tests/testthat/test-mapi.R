library("testthat")
library("mapi")
test_that("mapi", {
    expect_equal(nrow(mapi::samples), 200)
    expect_equal(nrow(mapi::metric), 19900)
})
