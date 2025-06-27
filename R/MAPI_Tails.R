#' @title Function MAPI_Tails
#' @export
#' 
#' @description Determine significant continuous and discontinuous areas from the result 
#'   of a MAPI analysis when run with permutations.
#' 
#' @param resu A spatial object of class 'sf' resulting from a MAPI analysis done using either
#'    \code{\link{MAPI_RunAuto}} or \code{\link{MAPI_RunOnGrid}}.
#' @param minQ Threshold under which cells with the smallest sum-of-weights percentile (range 1 .. 100) are discarded (default value = 0). 
#'   This parameter allows to discard cells for which the average value of the pairwise metric is computed 
#'   using either a small number and/or only long-distance ellipses.
#' @param alpha Significance level (default=0.05)
#' 
#' @return a spatial object of class 'sf' with the area and geometry of the polygons delineating the significant areas. 
#'    A column provides the tail for each polygon (upper or lower).
#'
#' @details When permutations are performed, in \code{\link{MAPI_RunOnGrid}} for each cell, the proportion of permuted values that are smaller or greater 
#'    than the observed value provides a lower-tailed (ltP) and upper-tailed (utP) test p-value.
#'    A false discovery rate (FDR) procedure (Benjamini and Yekutieli, 2001) is applied to account for multiple 
#'    testing (number of cells) under positive dependency conditions (spatial autocorrelation). An adjusted
#'    p-value is computed for each cell using the function \code{p.adjust} from the 'stats' package with the method 'BY'.
#'    The significance level at which FDR is controlled is set through the parameter alpha. For example, when alpha is
#'    set to 0.05, this means that 5\% of the cells detected as significant can be false positives.
#'    
#'    Significant cells belonging to the lower (or upper) tail that are spatially connected are aggregated 
#'    together to form the significant areas with the lowest (or greater) average values of the pairwise metric analyzed.
#' 
#' @examples
#' \dontrun{
#' data("metric")
#' data("samples")
#' # Run MAPI computation
#' resu <- MAPI_RunAuto(samples, metric, crs=3857, nbPermuts=1000)
#' # Discards the 10% cells with the smallest sum-of-weights 
#' #    and aggregates adjacent cells belonging to the same tail 
#' #    at a 5% significance level
#' tails <- MAPI_Tails(resu, minQ=10, alpha=0.05)
#' }
#'
#' @references 
#' Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165–1188.


MAPI_Tails <- function(resu, minQ=0, alpha=0.05) {
    # discard empty cells, if any
    resu <- resu[!is.na(resu$avg_value) , ]
    
    # check data availability
    if (any(colnames(resu)=="ltP") || any(colnames(resu)=="utP")) {
        message(sprintf("Significant areas aggregation, with percentile filter minQ > %f and alpha = %f ...", minQ, alpha))
    } else {
        stop("Error: no adjusted probabiblites for lower- and upper-tail columns (\"ltP\", \"utP\") found in results table provided.")
    }
    
    # Check confidence level value
    if (alpha < 0.0 || alpha > 0.5) {
        stop(sprintf("Incorrect value for parameter alpha: %f (Range: [0, 0.5])", minQ))
    }
    
    # Allow sum-of-weights percentile filtering
    if (minQ > 0) {
        my.resu <- resu[ resu$swQ > minQ , ]
    } else if (minQ <= 100) {
        my.resu <- resu
    } else {
        stop(sprintf("Incorrect value for parameter minQ: %f (Range: [0, 100])", minQ))
    }
    
    # TODO: check for consistence between alpha and nbPermuts!
    
    if (sf::st_is_longlat(my.resu)) {
        # We need to reconsider default snap distance, otherwise holes appear in merged polygons.
        # This is an alternative to the small buffer used in planar grids.
        # Let's estimate snap as small fraction of average cell radius
        snap.dist <- mean(sqrt(s2::s2_area(sf::st_geometry(my.resu)) /pi)) * 1e-5
        # IMPORTANT! Convert to radians
        snap.dist <- snap.dist / s2::s2_earth_radius_meters()
    }
    
    # get upper tail
    anyUpper <- any(c(my.resu$utP <= alpha, FALSE), na.rm=TRUE)
    if (anyUpper) {
        if (sf::st_is_longlat(my.resu)) {
            tails.up.g <- s2::s2_union_agg(sf::st_geometry(my.resu[my.resu$utP <= alpha, ]), options=s2::s2_options(model='closed', snap=s2::s2_snap_distance(snap.dist)))
            tails.up <- data.frame(tail=rep("upper", length(tails.up.g)))
            tails.up$area <- s2::s2_area(tails.up.g)
            tails.up$geometry <- tails.up.g
            tails.up <- sf::st_as_sf(tails.up, sf_column_name="geometry")
        } else {
            tails.up.g <- sf::st_cast(sf::st_buffer(sf::st_union(my.resu$geometry[my.resu$utP <= alpha]), 0.0001), 'POLYGON')
            tails.up <- data.frame(tail=rep("upper", length(tails.up.g)))
            sf::st_geometry(tails.up) <- tails.up.g
            tails.up$area <- sf::st_area(tails.up)
        }
    }
    # get lower tail
    anyLower <- any(c(my.resu$ltP <= alpha, FALSE), na.rm=TRUE)
    if (anyLower) {
        if (sf::st_is_longlat(my.resu)) {
            tails.low.g <- s2::s2_union_agg(sf::st_geometry(my.resu[my.resu$ltP <= alpha, ]), options=s2::s2_options(model='closed', snap=s2::s2_snap_distance(snap.dist)))
            tails.low <- data.frame(tail=rep("lower", length(tails.low.g)))
            tails.low$area <- s2::s2_area(tails.low.g)
            tails.low$geometry <- tails.low.g
            tails.low <- sf::st_as_sf(tails.low, sf_column_name="geometry")
        } else {
            tails.low.g <- sf::st_cast(sf::st_buffer(sf::st_union(my.resu$geometry[my.resu$ltP <= alpha]), 0.0001), 'POLYGON')
            tails.low <- data.frame(tail=rep("lower", length(tails.low.g)))
            sf::st_geometry(tails.low) <- tails.low.g
            tails.low$area <- sf::st_area(tails.low)
        }
    }
    # merge and returns tails (if any)
    if (!anyUpper && !anyLower) { 
        # geometries empty (ie. no tails)
        tails <- sf::st_cast(sf::st_buffer(sf::st_union(my.resu$geometry[my.resu$utP <= alpha]), 0.0001), 'POLYGON') # returns an empty geometry
        message("... no significant area")
    } else if (anyUpper && anyLower) { 
        # both geometries exists
        tails <- rbind(tails.up, tails.low)
        message(sprintf("... %d upper-tail and %d lower-tail significant areas returned", nrow(tails.up), nrow(tails.low)))
    } else if (anyUpper) { 
        # only upper tail exists
        tails <- tails.up
        message(sprintf("... %d upper-tail significant areas returned", nrow(tails.up)))
    } else if (anyLower) { 
        # only lower tail exists
        tails <- tails.low
        message(sprintf("... %d lower-tail significant areas returned", nrow(tails.low)))
    }
    # add an id if needed
    if (anyUpper || anyLower) {
        tails$id <- 1:nrow(tails)
    }
    return(tails)
}

