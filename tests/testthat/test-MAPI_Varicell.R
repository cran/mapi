library("testthat")
library("mapi")
test_that("MAPI_Varicell", {
    keep <- c(2,3,6,9,10,11,16,18,19,20,21,23,26,27,29,31,33,34,38,41,46,50,54,58,61,63,65,71,
    72,73,76,78,79,81,85,92,93,94,97,98,99,101,103,113,115,119,120,121,124,127,130,134,142,143,
    151,152,159,160,176,185,189,191,195,196,197,198,199) # for reproductibility
    samples2 <- mapi::samples[keep,]
    grid <- MAPI_GridHexagonal(samples2, crs=3857, hw=250)
    expect_equal(nrow(grid), 1045)
    expect_equal(as.integer(sum(sf::st_area(grid))), 169686852)
    set.seed(1234) # for reproductibility
    grid.var <- MAPI_Varicell(grid, samples2, buf=250, var.coef=20)
    expect_equal(nrow(grid.var), 173)
    expect_equal(as.integer(sum(sf::st_area(grid.var))), 189126625)
    expect_equal(sd(sf::st_area(grid.var)), 963339.563502512)
})
