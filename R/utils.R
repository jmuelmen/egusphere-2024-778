#' @export
#' 
between_quantiles <- function(x, min, max) {
    qs <- quantile(x, c(min, max))
    dplyr::between(x, qs[1], qs[2])
}

#' @export
#' 
esat <- function(T) {
    ## Tetens formula
    T <- T - 273.15
    0.61078e3 * exp((17.27 * T) / (T + 237.3))
}

#' @export
#' 
qsat <- function(T, p) {
    es <- esat(T)
    0.622 * es / (p - 0.378 * es)
}

#' @export
#' 
equivalent_factor <- function(T, p) {
    qs <- qsat(T, p)
    exp((Lv * qs) / (cp * T))
}

#' Gamma is the adiabatic condensation rate in kg m^-4
#' 
#' @export
#' 
adiabatic_nd_from_tau_re <- function(tau, re, Gamma = 2e-6) {
    sqrt(10) / (4 * pi * sqrt(1e3)) * sqrt(Gamma) * sqrt(tau) * re ^ (-5/2)
}

#' @export
#' 
adiabatic_lwp_from_tau_re <- function(tau, re) {
    5 * 1e3 * tau * re / 9
}

#' @export
#'
regionalize_zhang <- function(lon, lat) {
    ifelse(
        between(lon, -140, -120) & between(lat, 15, 35), "NEP",
    ifelse(
        between(lon, -90, -70) & between(lat, -30, -10), "SEP",
    ifelse(
        between(lon, -10, 10) & between(lat, -25, -5), "SEA",
    ifelse(
        between(lon, -40, -20) & between(lat, 10, 30), "NEA",
    ifelse(
        between(lon, 90, 110) & between(lat, -35, -15), "AUS",
        NA
           )))))
}

#' @export
#'
regionalize_seasonal <- function(month, lon, lat, fscu.annual, landfrac) {
    ifelse(landfrac > 0 | fscu.annual <= 0.3 | abs(lat) > 45,
           NA,
    ifelse(lat > 0 & month %in% c("Jun", "Jul", "Aug"),
    ifelse(lon < 30 | lon > 330,
           "NEA",
           "NEP"),
    ifelse(lat < 0 & month %in% c("Oct", "Nov", "Dec", "Jan", "Feb"),
    ifelse(lon < 30 | lon > 330,
           "SEA",
    ifelse(lon > 60 & lon < 120,
           "AUS",
           "SEP")),
    NA)))
}
