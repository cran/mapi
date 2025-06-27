#' @title Function MAPI_EstimateHalfwidth
#' @export
#' @description This function computes the side length (= halfwidth) of the hexagonal cells. Halfwidth value can be further used to build a MAPI grid.
#' 
#' @param samples a data.frame with names and geographical coordinates of samples. Column names must be: 'ind', 'x', 'y'.  
#'   Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude).
#' @param crs coordinate reference system: integer with the EPSG code, or character with proj4string. 
#'   The coordinates system must be a projection, not latitude/longitude.
#' @param beta A value depending on spatial regularity of sampling: 0.5 for regular sampling, 0.25 for random sampling (Hengl, 2006).
#' 
#' @details 
#' \eqn{h_w = \frac{\beta \sqrt{A/N}}{\sqrt{2.5980}}}{hw = p / sqrt(2.5980) with p = beta * sqrt(A/N)}, where A is the study area (convex hull of sampling points) and N the number of samples. 
#' Parameter beta allows to respect the Nyquist-Shannon sampling theorem depending on sampling regularity.
#'
#' @return halfwidth cell value (side length of hexagonal cells). 
#' 
#' @examples
#' data(samples)
#' # Computes hexagonal cell halfwidth for the 'samples' dataset using beta=0.5
#' hw <- MAPI_EstimateHalfwidth(samples, crs=3857, beta=0.5)
#'
#' @references
#' Hengl, T. (2006) Finding the right pixel size. Computers & Geosciences, 32, 1283--1298.
#' 

MAPI_EstimateHalfwidth <- function(samples, crs, beta=0.25) {
    my.crs <- sf::st_crs(crs)
    samples2 <- sf::st_as_sf(samples, coords=c("x", "y"), crs=my.crs)
    if (sf::st_is_longlat(samples2)) {
        message("NOTE: coordinates are in angular units, note that x stands for longitude and y for latitude.")
        my.ch <- s2::s2_convex_hull_agg(sf::st_as_s2(samples2$geometry))
        # whether convex hull exists or not...
        if (s2::s2_is_empty(my.ch)) {
            # return earth area
            # Classical calculation : A <- 4 * pi * s2::s2_earth_radius_meters() * s2::s2_earth_radius_meters()
            # S2 calculation should be (and actually is) the same! :-D
            A <- s2::s2_area(s2::as_s2_geography(TRUE))
        } else {
            # return area of the convex hull
            A <- s2::s2_area(my.ch)
        }
    } else {
        A <- sf::st_area(sf::st_convex_hull(sf::st_union(samples2$geometry)))
    }
    A <- as.double(A)
    N <- nrow(samples)
    p <- beta*sqrt(A/N)
    hw <- p/sqrt(2.5980) # NOTE: Should be 2.598076211 but rounded for compatibility with old SQL version
    message(sprintf("Estimated halfwidth: %f", hw))
    return(hw)
}
