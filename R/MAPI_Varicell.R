#' @title Function MAPI_Varicell
#' @export
#' 
#' @description Alter grid by changing polygons size according to samples density.
#' 
#' @param grid a spatial object of class 'sf' with the geometry of each cell. 
#'   When using your own grid, please check that the object structure is the same as returned by 
#'   \code{\link{MAPI_GridAuto}} or \code{\link{MAPI_GridHexagonal}}.
#'   Note that the computation of the cartogram may shrink the grid close to samples; therefore the grid 
#'   should be built with a very large buffer around the sampling points.
#' @param samples a data.frame with names and geographical coordinates of samples. Column names must be: 'ind', 'x', 'y'.  
#'   Optional column 'errRad' with an error radius for sample locations (eg. GPS uncertainty). 
#'   Coordinates must be projected (not latitude/longitude) and in the same projection as the grid.
#' @param var.coef optional, default value=50. The expected ratio of size between largest and smallest cells in final grid.
#' @param buf This parameter allows to clip the grid around samples by a number of units in 
#'   the same reference system as the grid geographical coordinates (0 by default). A value 3-5 times smaller
#'   than the size used to build the grid may be OK.
#' 
#' @return a spatial object of class 'sf' including cell area, the x and y coordinates of cell centroids, 
#'   cell geometry (polygons) and cell id (gid) plus the computed densities and weights.
#' 
#' 
#' @examples
#' \dontrun{
#' # Dummy example!
#' require(data.table, sf, mapi, spatstat.geom, spatstat.core)
#' data("samples")
#' keep <- c(2,3,6,9,10,11,16,18,19,20,21,23,26,27,29,31,33,34,38,41,46,50,54,58,61,63,65,71,72,73,
#' 76,78,79,81,85,92,93,94,97,98,99,101,103,113,115,119,120,121,124,127,130,134,142,143,151,152,159,
#' 160,176,185,189,191,195,196,197,198,199)
#' samples2 <- samples[keep, ] # keep only one third in order to create discontinuities
#' set.seed(1234)
#' # Builds a grid of hexagonal cells according to samples coordinates (columns x and y) 
#' # using the EPSG:3857 projection and an halfwidth cell value of hw=250m.
#' grid <- MAPI_GridHexagonal(samples2, crs=3857, hw=250, buf=1500)
#' grid.var <- MAPI_Varicell(grid, samples2, buf=250, var.coef=20)
#' ggplot() + 
#'     geom_sf(data=grid.var, aes(fill=dens), size=0.1) + 
#'     geom_sf(data=st_as_sf(samples2, coords=c("x","y"), crs=3857), aes())
#' }
#'

MAPI_Varicell <- function(grid, samples, buf=0, var.coef=50) {
    
    # Convert samples into spatial sf object
    samples.g <- sf::st_as_sf(samples, coords=c("x", "y"), crs=sf::st_crs(grid))
    
    # ROIs
    study.area <- sf::st_buffer(sf::st_convex_hull(sf::st_union(grid$geometry)), buf)
 
    # mean area then halfwidth
    A <- mean(as.numeric(sf::st_area(grid))) # average cell area
    hw.mean <- sqrt(2*A/(3*sqrt(3))) # regular hexagon formula
    zz <- (3*sqrt(3))/2 # NOTE: corrction factor such as the grid with var.coef=1 has more or less the same resolution as the original one
    
    # build the mesh
    mesh <- fmesher::fm_mesh_2d_inla(sf::st_coordinates(samples.g),
                       boundary = list(study.area),
                       max.edge = c(hw.mean*zz*sqrt(var.coef)), cutoff=hw.mean*zz/sqrt(var.coef))
    # get mesh points
    my.ppp <- data.table::as.data.table(mesh$loc)
    colnames(my.ppp) <- c("X","Y","Z")
    my.ppp.g <- sf::st_as_sf(my.ppp, coords=c("X","Y"), crs=sf::st_crs(grid))
    
    # build tessellation from mesh points
    vor0 <- sf::st_voronoi(sf::st_union(my.ppp.g), envelope=study.area)
    # convert to sf polygons
    vor1 <- sf::st_as_sf(sf::st_collection_extract(vor0))
    vor1 <- sf::st_set_geometry(vor1, "geometry")
    # cut within area
    vor <- sf::st_intersection(vor1, study.area)
    # add centroids
    coor <- data.table::as.data.table(sf::st_coordinates(sf::st_centroid(vor)))
    colnames(coor) <- tolower(colnames(coor))
    grid.out <- cbind(vor, coor)
    # cell area
    grid.out$area <- as.numeric(sf::st_area(vor))
    # add gid column
    grid.out$gid <- seq(1,nrow(grid.out))
    # Returns the result
    return(grid.out)
}
