#' Map ne120pg2 grid points to closest ne30pg2 points
#'
#' @param dist distance metric between grid points; defaults to
#'     plotutils::dist.gc(), expected signature is (lon1, lon2, lat1, lat2)
#' @param ... arguments passed to ddply
#'
#' @examples
#' df.120.map <- map_fine_to_coarse(df.120, df.30, .parallel = TRUE, .progress = "text")
#' df.120.map <- map_fine_to_coarse(df.120, df.30,
#'     dist = function(lon1, lon2, lat1, lat2) sqrt((lon1 - lon2)^2 + (lat1 - lat2)^2), 
#'     .parallel = TRUE, .progress = "text")
#' @export
map_fine_to_coarse <- function(fine, coarse, dist = plotutils::dist.gc,
                               ...) {
    plyr::ddply(fine, ~ ncol, function(df ) {
        coarse %>%
            dplyr::mutate(dist = dist(df$lon, lon,
                                      df$lat, lat)) %>%
            filter(dist == min(dist)) %>%
            transmute(dist,
                      ncol.coarse = ncol)
    }, ...)
}

#' Map ne120pg2 grid points to closest ne30pg points 
#'
#' @export
dist_fine_to_coarse <- function(fine, coarse, ...) {
    refine <- nrow(fine) / nrow(coarse)
    fine %>%
        dplyr::mutate(ncol.coarse = floor(ncol / refine) + 1,
                      dist = plotutils::dist.gc(coarse$lon[ncol.coarse], lon,
                                                coarse$lat[ncol.coarse], lat))
}
