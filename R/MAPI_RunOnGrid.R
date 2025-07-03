#' @title Function MAPI_RunOnGrid
#' @export
#' @description Launch a MAPI analysis for a given grid computed with \code{\link{MAPI_GridAuto}} 
#'   or \code{\link{MAPI_GridHexagonal}} or provided by users.
#' 
#' @param samples a data.frame with names and geographical coordinates of samples. 
#'   Column names must be: 'ind', 'x', 'y'.  Optional column 'errRad' with an error radius 
#'   for sample locations (eg. GPS uncertainty).
#'   Coordinates should be projected (not latitude/longitude) for a local dataset, angular coordinates 
#'   (ie. latitude/longitude) are accepted for worldwide datasets.
#'   In this case, the longitude will be in the "x" column and the latitude in the "y" column.
#' @param metric a data.frame or a square matrix with the pairwise metric computed for all pairs of samples. 
#'   If data.frame, column names must be: 'ind1', 'ind2', 'value'.  
#'   If matrix, sample names must be the row- and column names.
#' @param grid a spatial object of class 'sf' with the geometry of each cell. 
#'   When using your own grid, please check that the object structure is the same as returned by 
#'   \code{\link{MAPI_GridAuto}} or \code{\link{MAPI_GridHexagonal}}.
#' @param isMatrix Boolean. Depends on the 'metric' data:\cr
#'   TRUE if 'metric' is a square matrix with column names = row names and standing for sample names.\cr
#'   FALSE if 'metric is a three columns data.frame ('ind1', 'ind2', 'value'). \cr
#'   The default value is determined using a "matrix" class detection for metric as well as identity between row and column number.
#' @param ecc ellipse eccentricity value (0.975 by default).
#' @param errRad global error radius for sample locations (same radius for all samples, 10 by default). 
#'   Units are in the same reference system as the sample geographical coordinates.
#'   To use different error radius values for sample locations, add a column 'errRad' in the 'sample' data (see \code{\link{mapi}}).
#' @param nbPermuts number of permutations of sample locations (0 by default).
#' @param dMin minimum distance between individuals. 0 by default.
#' @param dMax maximal distance between individuals. +Inf by default.
#' @param use_s2 transform grid and coordinates in latitude,longitude in order to use s2 library functions (default FALSE).
#'   When set to TRUE mapi is able to process worldwide grids and ellipses and generate results on the sphere.
#' @param nbCores number of CPU cores you want to use during parallel computation. 
#'   The default value is estimated as the number of available cores minus 1, suitable for a personal computer. 
#'   On a cluster you might have to set it to a reasonable value (eg. 8) in order to keep resources for other tasks. 
#' @param N number of points used per quarter of ellipse, 8 by default. 
#'   Don't change it unless you really know what you are doing.
#' @param ignore.weights allow to compute a simple mean instead of the default weighted mean.
#'   This implies that all ellipses bear the same weight.
#' 
#' @return a spatial object of class 'sf' providing for each cell: \cr
#'  - gid: Cell ID \cr
#'  - x and y coordinates of cell center \cr
#'  - nb_ell: number of ellipses used to compute the weighted mean \cr
#'  - avg_value: weighted mean of the pairwise metric \cr
#'  - sum_wgts: sum of weights of ellipses used to compute the weighted mean \cr
#'  - w_stdev: weighted standard deviation of the pairwise metric \cr
#'  - swQ: percentile of the sum of weights \cr
#'  - geometry \cr
#'  When permutations are performed: \cr
#'  - proba: proportion of the permuted weighted means below the observed weighted mean \cr
#'  - ltP: lower-tail p-value adjusted using the FDR procedure of Benjamini and Yekutieli \cr
#'  - utP: upper-tail p-value adjusted using the FDR procedure of Benjamini and Yekutieli \cr
#' 
#' @details
#' To test whether the pairwise metric values associated with the ellipses are independent of the sample locations, those are permuted 'nbPermuts' times. 
#'   At each permutation, new cell values are computed and stored to build a cumulative null distribution for each cell of the grid. 
#'   Each cell value from the observed data set is then ranked against its null distribution.
#'   For each cell, the proportion of permuted values that are smaller or greater 
#'   than the observed value provides a lower-tailed (ltP) and upper-tailed (utP) test p-value.
#'   
#' A false discovery rate (FDR) procedure (Benjamini and Yekutieli, 2001) is applied to account for multiple 
#'    testing (number of cells) under positive dependency conditions (spatial autocorrelation).  An adjusted
#'    p-value is computed for each cell using the function \code{p.adjust} from the 'stats' package with the method 'BY'.
#'    
#' @references 
#' Benjamini, Y. and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165–1188.
#' 
#' @examples
#' \dontrun{
#' data(metric)
#' data(samples)
#' my.grid <- MAPI_GridHexagonal(samples, crs=3857, 500) # 500m halfwidth
#'
#' # Note: 10 permutations is only for test purpose, increase to >=1000 in real life!
#' my.results <- MAPI_RunOnGrid(samples, metric, grid=my.grid, nbPermuts=10, nbCores=1)
#' 
#' # eg. Export results to shapefile "myFirstMapiResult" in current directory
#' sf::st_write(my.results, dsn=".", layer="myFirstMapiResult", driver="ESRI Shapefile")
#' }
#' 


MAPI_RunOnGrid <- function(samples, metric, grid, isMatrix=FALSE, ecc=0.975, errRad=10, nbPermuts=0, dMin=0, dMax=Inf, use_s2=FALSE, nbCores=parallel::detectCores()-1, N=8, ignore.weights=FALSE) {
    message("MAPI COMPUTATION STARTED")
    if (!("sf" %in% class(grid))) {stop("grid parameter is not a 'sf' object!")}
    if (base::requireNamespace("doParallel", quietly=TRUE) && base::requireNamespace("foreach", quietly=TRUE)) {
        `%dopar%` <- foreach::`%dopar%`
    } else {
        nbCores <- 1
        `%dopar%` <- NA
    }
    tot <- system.time({
        my.data.in <- MAPI_CheckData(samples, metric, isMatrix=isMatrix)
        my.samples <- my.data.in[[1]]
        my.metric <- my.data.in[[2]]
        # prepare index
        data.table::setkey(my.metric, "ind1", "ind2")
        
        # fill missing errRad with parameter value
        if (base::any(base::colnames(my.samples)=="errRad")) {
            my.samples[base::is.na(my.samples$errRad), "errRad"] <- errRad
        } else {
            my.samples$errRad <- errRad
        }
        # get crs value from grid geometry
        crs0 <- sf::st_crs(grid) ; crs <- crs0

        # add geometry (point) to samples
        my.samples.geom <- sf::st_as_sf(my.samples, coords=c("x", "y"), crs=crs0, remove=FALSE)
        
        # force using s2 functions
        if (use_s2) {
            crs <- "EPSG:4326"
            my.samples.geom <- sf::st_transform(my.samples.geom, crs)
            my.samples.geom$x <- as.numeric((sf::st_coordinates(my.samples.geom)[,1]))
            my.samples.geom$y <- as.numeric((sf::st_coordinates(my.samples.geom)[,2]))
            grid <- sf::st_transform(grid, crs)
        }
        
        # Flag to true if angular coordinates and spherical geometry applies
        mapi.spherical <- sf::st_is_longlat(grid)

        if (mapi.spherical) {
            ## Companion function for worldwide ellipses
            mkEll_s2 <- function(i0, i1, loc.pts, err.circ, N=8, ecc=0.975, fracs=NULL) {
                if (i0 == i1) {
                    # Same location & same error circle
                    return(err.circ[i0])
                } else {
                    # Compute distance fractions if not provided
                    if (is.null(fracs)) {
                        N <- 2 * N
                        epsilon <- 0.01 # avoids "The three spheres do not intersect!" errors and still close enough of the apex
                        angles <- c(epsilon, base::seq(pi/N, pi-pi/N, by=pi/N), pi-epsilon)
                        fracs <- (cos(angles) + 1) / 2
                    }
                    # Get points
                    p0 <- loc.pts[i0]
                    p1 <- loc.pts[i1]
                    
                    if (s2::s2_equals(p0, p1)) {
                        # Returns largest error circle if two points are superposed with different error radius
                        c0 <- err.circ[i0]
                        c1 <- err.circ[i1]
                        return(s2::s2_union_agg(c(c0, c1)))
                    } else {
                        # Get points coordinates
                        p0x <- s2::s2_x(p0)
                        p0y <- s2::s2_y(p0)
                        p1x <- s2::s2_x(p1)
                        p1y <- s2::s2_y(p1)
                        # Great circle distance between points
                        d <- s2::s2_distance(p0, p1)
                        # Ellipse distance 
                        dd <- d / ecc
                        # compute point pairs
                        pts1 <- NULL
                        pts2 <- NULL 
                        for(fr in fracs) {
                            # great circle distance to ellipse border (first focus)
                            gc0 <- d * fr + (dd-d)/2
                            # complementary great circle distance to ellipse border (second focus)
                            gc1 <- dd - gc0

                            # trilateration of:
                            #   earth, 
                            #   1st sphere based on first focus and radius from first great circle distance,
                            #   2nd sphere based on second focus and radius from second great circle distance
                            tp.cpp <- .trilaterate_cpp(p0x, p0y, p1x, p1y, gc0, gc1)

                            # get coordinates for the two points
                            tp1.cpp <- tp.cpp[1:2]
                            tp2.cpp <- tp.cpp[3:4]

                            # build the two halves, one in reverse order in order to loop in the same direction
                            if (! base::any(base::is.na(tp1.cpp)))
                                pts1 <- base::rbind(pts1, tp1.cpp)
                            if (! base::any(base::is.na(tp2.cpp)))
                                pts2 <- base::rbind(tp2.cpp, pts2)
                        }
                        # asssemble the two halves
                        pts.s2 <- base::rbind(pts1, pts2)
                        # builds ellipse polygon from points coordinates
                        ell.s2 <- s2::s2_make_polygon(pts.s2[,1], pts.s2[,2])
                        # tryCatch({ ell.s2 <- s2::s2_make_polygon(pts.s2[,1], pts.s2[,2]) }, error=function(e) { print(pts.s2) ; stop(e) })
                        # if first error circle is wider than ellipse then builds the convex hull of both
                        c0 <- err.circ[i0]
                        if (!s2::s2_contains(ell.s2, c0)) {
                            ell.s2 <- s2::s2_convex_hull_agg(c(ell.s2, c0))
                        }
                        # if second error circle is wider than ellipse (or previous convex hull) then builds the convex hull of both
                        c1 <- err.circ[i1]
                        if (!s2::s2_contains(ell.s2, c1)) {
                            ell.s2 <- s2::s2_convex_hull_agg(c(ell.s2, c1))
                        }
                        # done; let's return the result
                        return(ell.s2)
                    }
                }
            }
        }
        
        ## Let's compute ellipses
        t <- system.time({
            # create locality code from hex representation of geometry and error radius
            my.samples.geom$locCode <- base::as.character(sf::st_as_binary(my.samples.geom$geometry, hex=TRUE))
            my.samples.geom$locCode <- base::apply(my.samples.geom, 1, function(r) { base::paste(r["locCode"], r["errRad"]) })
            
            # get distinct localities with geometry
            distinct.locations <- base::unique(base::as.data.frame(my.samples.geom[, c("locCode", "geometry", "errRad", "x", "y")]))
            distinct.locations$dlid <- 1:base::nrow(distinct.locations)
            
            # precompute error circles once for all
            if (mapi.spherical) {
                location.pts <- s2::s2_lnglat(distinct.locations$x, distinct.locations$y)
                error.circles <- s2::s2_geography()
                for (i in 1:length(location.pts)) {
                    # simplify buffer in order to speedup computations using the same tolearance as ellipses (1/4N --> 32 segments by circle for N=8)
                    error.circles[i] <- s2::s2_simplify( s2::s2_buffer_cells(location.pts[i], distance=distinct.locations$errRad[i]) , tolerance=2*pi*(distinct.locations$errRad[i])/(4*N))
                }
            }
            
            # prepare pairs of points as cartesian product of distinct localities ...
            ellipses <- data.table::as.data.table(expand.grid(loc1=distinct.locations$locCode, loc2=distinct.locations$locCode))
            # ... and keep only the half-matrix, including diagonal as two distinct individuals may share a same locality
            ellipses <- ellipses[base::as.character(ellipses$loc1) <= base::as.character(ellipses$loc2) , ]
            
            # merge half-matrix with locations
            ellipses <- data.table::merge.data.table(ellipses, distinct.locations, by.x="loc2", by.y="locCode")
            ellipses <- data.table::merge.data.table(ellipses, distinct.locations, by.x="loc1", by.y="locCode", suffixes=c("2","1"))
            
            # computes distance between the two localities (Orthodromic or Pythagore)
            if (mapi.spherical) {
                if (nbCores > 1) {
                    doParallel::registerDoParallel(nbCores)
                    ellipses$dist <- foreach::foreach(i=1:nrow(ellipses), .combine=base::c) %dopar% {
                        i1 <- ellipses$dlid1[i]
                        i2 <- ellipses$dlid2[i]
                        p1 <- location.pts[i1]
                        p2 <- location.pts[i2]
                        d <- s2::s2_distance(p1, p2)
                        return(d)
                    }
                    doParallel::stopImplicitCluster()
                } else {
                    ellipses$dist <- base::apply(ellipses[ , c("dlid1","dlid2")], 1, function(r) {
                        return(s2::s2_distance(location.pts[r[1]], location.pts[r[2]]))
                    })
                }
            } else {
                ellipses$dist <- base::sqrt( (ellipses$x1 - ellipses$x2)^2 + (ellipses$y1 - ellipses$y2)^2 )
            }
            
            # distance filtering out (if any filter set)
            if (!base::is.na(dMin) && dMin > 0)                             { ellipses <- ellipses[ellipses$dist >= dMin, ] }
            if (!base::is.infinite(dMax) && !base::is.na(dMax) && dMax > 0) { ellipses <- ellipses[ellipses$dist <= dMax, ] }
            
            base::message(base::sprintf("Building %s elliptical polygons...", base::nrow(ellipses)))
            
            # sort the completed half-matrix [NOTE: why?] and add rowid
            ellipses <- ellipses[base::order(ellipses$dlid1, ellipses$dlid2) , ]
            ellipses$rowid <- 1:base::nrow(ellipses)

            if (mapi.spherical) {
                # get distinct locations: ids of each pair
                ellipses.matrix <- base::as.matrix(ellipses[,c("dlid1", "dlid2")])
                
                # precomute length fractions for ellipses tracing
                epsilon <- 0.01
                angles <- base::c(epsilon, base::seq(pi/(2*N), pi-pi/(2*N), by=pi/(2*N)), pi-epsilon)
                fracs <- (base::cos(angles) + 1) / 2
                
                # build ellipses on sphere using companion functions
                if (nbCores > 1) { 
                    # parallel computation
                    doParallel::registerDoParallel(nbCores)
                    doIt <- function(todo) {
                        foreach::foreach(i=todo, .export=c("mkEll_s2"), .combine=c) %dopar% {
                            # return as text as pointers vanishes from cluster
                            s2::s2_as_text(mkEll_s2(ellipses.matrix[i,1], ellipses.matrix[i,2], location.pts, error.circles, ecc=ecc, fracs=fracs))
                        }
                    }
                    ePolysT <- doIt(1:nrow(ellipses.matrix))
                    doParallel::stopImplicitCluster()
                    ePolys <- s2::as_s2_geography(ePolysT) # Therefore, we converts back to s2 geography
                } else { 
                    # serial computation
                    ePolys <- s2::s2_geography()
                    for (i in (1:nrow(ellipses.matrix))) {
                        r <- ellipses.matrix[i,]
                        ePolys[i] <- mkEll_s2(r[1], r[2], location.pts, error.circles, ecc=ecc, fracs=fracs)
                    }
                }
            } else {
                # get distinct locations: x,y coordinates and error radius of each point
                ellipses.matrix <- as.matrix(ellipses[,c("x1", "y1", "x2", "y2", "errRad1", "errRad2")])
                # call to MAPI c++ function .mkP4st_cpp for building the polygon (convex hull of ellipse and error circles)
                ePolys <- apply(ellipses.matrix, 1, function(r) {
                    sf::st_polygon(list(.mkP4st_cpp(r, N, ecc)))
                })
            }
            # simplify ellipses table
            ellipses <- ellipses[,c("rowid", "loc1", "loc2", "dist")]
            # compute weight
            # set geometry to ellipses table (same order !)
            if (mapi.spherical) {
                # compute area
                ellipses$area <- s2::s2_area(ePolys)
                ellipses$geometry <- ePolys
            } else {
                sf::st_geometry(ellipses) <- sf::st_sfc(ePolys, crs=crs)
                # compute area
                ellipses$area <- sf::st_area(ellipses$geometry) # area in "units"
            }
            if (ignore.weights) {
                ellipses$weight <- 1.0
            } else {
                ellipses$weight <- 1.0 / as.numeric(ellipses$area)
            }
            
            # free memory
            rm(ellipses.matrix)
            
            # For each ellipse we keep a simplified table with its rowid, weight and the two locality codes
            ells <- data.table::data.table(loc1=ellipses$loc1, loc2=ellipses$loc2, rowid=ellipses$rowid, weight=ellipses$weight)
            data.table::setkey(ells, "rowid")
        })
        tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
        # NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
        message(sprintf("... %s ellipses.    [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", as.character(nrow(ellipses)), tv[1]+tv[4], tv[2]+tv[5], tv[3]))

        ## Spatial intersect between ellipses and grid
        # NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
        message(sprintf("Computing spatial intersection between %s grid cells and %s ellipses...", as.character(nrow(grid)), as.character(nrow(ellipses))))
        t <- system.time({
            if (mapi.spherical) {
                g2 <- sf::st_as_s2(grid)
                if (nbCores > 1) { 
                    # parallel computation
                    doParallel::registerDoParallel(nbCores)
                    inter0 <- foreach::foreach(i=1:length(g2)) %dopar% {
                        which(s2::s2_intersects(g2[i], ellipses$geometry))
                    }
                    doParallel::stopImplicitCluster()
                } else { 
                    # serial computation
                    inter0 <- base::lapply(g2, FUN=function(cs2){
                        which(s2::s2_intersects(cs2, ellipses$geometry))
                    })
                }
            } else {
                inter0 <- sf::st_intersects(grid, ellipses, sparse=TRUE, prepared=TRUE)
            }
        })
        tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
        message(sprintf("... done.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
        # We are finished with the ellipses geometry, let's free some memory
        rm(ellipses) # free memory
        
        
        ## Build sample pairs
        message("Building sample pairs...")
        t <- system.time({
            # build a simple table with sample code and locality code only
            my.sampleCode.locCode <- data.table::data.table(ind=my.samples.geom$ind, locCode=my.samples.geom$locCode)
            # cartesian product of samples
            sampX <- data.table::as.data.table(expand.grid(ind1=my.sampleCode.locCode$ind,ind2=my.sampleCode.locCode$ind))
            # complete matrix *without* diagonal
            sampX <- sampX[as.character(sampX$ind1) != as.character(sampX$ind2) , ]
            sampX <- sampX[order(sampX$ind1, sampX$ind2) , ]
            # merge for getting locality codes for both samples
            sampX <- merge(sampX, my.sampleCode.locCode, by.x="ind2", by.y="ind")
            sampX <- merge(sampX, my.sampleCode.locCode, by.x="ind1", by.y="ind", suffixes=c("2","1"))
            # merge with ellipses for weight
            sampX_12 <- merge(sampX, ells, by.x=c("locCode1", "locCode2"), by.y=c("loc1", "loc2"))
            sampX_21 <- merge(sampX, ells, by.x=c("locCode1", "locCode2"), by.y=c("loc2", "loc1"))
            sampX <- unique(rbind(sampX_12, sampX_21))
            # merge with metric (already symmetrized) for getting value
            sampX <- merge(sampX, my.metric, by=c("ind1", "ind2"), all.x=TRUE)
            # sort and compute a rowid
            sampX <- sampX[order(sampX$ind1, sampX$ind2) , ]
            sampX$id <- 1:nrow(sampX)
            rm(my.sampleCode.locCode, sampX_12, sampX_21) # free memory
        })
        tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
        # NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
        message(sprintf("... %s sample pairs.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", as.character(nrow(sampX)), tv[1]+tv[4], tv[2]+tv[5], tv[3]))
        
        
        ## Expand results from intersection between grid and ellipses into a clean list of vectors of sample pairs ids
        message("Matching grid cells with sample pairs...")
        t <- system.time({
            # sample pairs with symmetrized localities for faster merges
            sampX2 <- rbind(sampX[ , c("id", "locCode1", "locCode2")], sampX[ , c("id", "locCode2", "locCode1")])
            data.table::setkey(sampX2, "locCode1", "locCode2")
            # iterate on intersection
            inter <- lapply(inter0, function(r) {
                # get matched ellipses row numbers
                ids <- as.vector(r)
                # get matched ellipses
                myElls <- ells[ids, ]
                # join with symmetrized sample pairs
                tmp1 <- merge(myElls, sampX2, by.x=c("loc1", "loc2"), by.y=c("locCode1", "locCode2"))
                # extract unique sample pair ids
                v <- as.vector(unique(tmp1$id))
                return(v)
            })
            # Fast count
            nbMatches <- .countMatches_cpp(inter)
        })
        tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
        if (is.na(nbMatches)) { nbMatches <- 0 }
        # NOTE: %s used for stringified integers due to overflow for very large datasets (thanks to Simon Dellicour)
        message(sprintf("... %s matches found.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", as.character(nbMatches), tv[1]+tv[4], tv[2]+tv[5], tv[3]))
        
        
        ## Get results for unpermuted MAPI analysis
        # For each cell, we compute the number of intersecting ellipses, the sum of their weightd, the weighted mean and weighted standard deviation.
        message("Computing values for grid cells...")
        t <- system.time({
            # Call to C++ function parseInter_cpp which iterates on grid, intersetion and sample pairs for computing values
            resu <- data.table::as.data.table(.parseInter_cpp(grid$gid, inter, sampX$weight, sampX$value))
            colnames(resu) <- c("gid", "nb_ell", "avg_value", "sum_wgts", "w_stdev")
            data.table::setkey(resu, "gid")
            # computes sum-of-weights percentile
            resu$swQ <- unclass(.bincode(resu$sum_wgts, stats::quantile(resu$sum_wgts, 0:100/100.0, na.rm=TRUE), include.lowest=TRUE))
            # merge to grid for getting geometry
            resu <- merge(grid, resu, by="gid", all.x=TRUE)
        })
        tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
        message(sprintf("... done.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
        

        
        ## Let's start permutations block (if any)
        ##########################################################################################################################################################
        if(nbPermuts > 0) {
            t <- system.time({
                message(sprintf("Starting %d permutations...", nbPermuts))
                
                # prepare permuted samples table
                my.samples_perm <- data.table::data.table(ind=my.samples$ind)
                data.table::setkey(my.samples_perm, "ind")
                
                # run permutation function (parallel or not)
                if (nbCores > 1) {
                    # parallelized
                    message(sprintf("... parallelized over %d cores ...", nbCores))
                    doParallel::registerDoParallel(nbCores)
                    doIt <- function(todo) {
                        resu2 <- foreach::foreach( i=todo, .export=c("my.samples_perm", "my.samples", "sampX", "my.metric", "inter", "grid"), .verbose=FALSE, .combine=cbind ) %dopar% {
                            cat(sprintf("%d ",i))
                            .doPerm(i, my.samples_perm, my.samples, sampX, my.metric, inter, grid)
                        }
                    }
                    resu2 <- doIt(1:nbPermuts)
                    cat("\n")
                    doParallel::stopImplicitCluster()
                } else {
                    # serially
                    message(sprintf("... serially ..."))
                    pbs <- utils::txtProgressBar(max=nbPermuts, style=3, width=80)
                    resu2 <- data.table::as.data.table(lapply(1:nbPermuts, function(i){
                        pp <- .doPerm(i, my.samples_perm, my.samples, sampX, my.metric, inter, grid)
                        utils::setTxtProgressBar(pbs, i)
                        return(pp)
                    }))
                    close(pbs) # end of progress bar
                }
                
                message("... computing probabilities from permuted values ...")
                resu$proba <- NA
                resu2 <- as.matrix(resu2) # faster for getting a row ;-)
                for(i in 1:nrow(resu)) {
                    val <- resu$avg_value[i]
                    perms <- as.numeric(resu2[i,])
                    resu$proba[i] <- ( sum(perms < val) / nbPermuts )
                }
                
                # Benjamini-Yekutieli correction of probabilities for upper and lower tails
                resu$ltP <- stats::p.adjust(resu$proba, method="BY")
                resu$utP <- stats::p.adjust((1.0 - resu$proba), method="BY")
                
            })
            tv <- as.vector(t) ; tv[is.na(tv)] <- 0.0
            message(sprintf("... done.  [user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
        }
        ##########################################################################################################################################################
    })
    
    # revert results to initial crs
    if (use_s2) {
        resu <- sf::st_transform(resu, crs0)
    }
    
    tv <- as.vector(tot) ; tv[is.na(tv)] <- 0.0
    message(sprintf("MAPI COMPUTATION ENDED.  [TOTAL: user: %0.3f, system: %0.3f, elapsed: %0.3f seconds]", tv[1]+tv[4], tv[2]+tv[5], tv[3]))
    return(resu)
}



#### Internal function

#' @noRd
# prepare permutation function
.doPerm <- function(permut, my.samples_perm, my.samples, sampX, my.metric, inter, grid) {
    # shuffle individuals
    my.samples_perm$indP <- base::sample(my.samples$ind)

    # replace unpermuted individuals by permuted ones
    sampX.P <- data.table::merge.data.table(sampX,   my.samples_perm, by.x="ind1", by.y="ind")
    sampX.P <- data.table::merge.data.table(sampX.P, my.samples_perm, by.x="ind2", by.y="ind", suffixes=c("1","2"))

    # merge permuted samplepairs with metric
    sampX.P <- data.table::merge.data.table(sampX.P, my.metric, all.x=TRUE, by.x=c("indP1","indP2"), by.y=c("ind1","ind2"))
    data.table::setkey(sampX.P, "id")

    # compute permuted result
    return( .parseInterPerm_cpp(grid$gid, inter, sampX.P$weight, sampX.P$value.y) )
}
