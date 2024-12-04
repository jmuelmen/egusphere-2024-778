#' ECMWF L137 hybrid sigma coordinate coefficients
#'
#' A dataset containing the hybrid sigma coordinate coefficients used
#' in ECMWF 137-level model configurations.  The pressure at level
#' \eqn{i} is: \eqn{p_i = a_i + b_i p_\text{surf}}.
#'
#' @format A data frame with 3 variables:
#' \describe{
#'   \item{lev}{level; 0 is model top, 137 is surface}
#'   \item{a}{coefficient \eqn{a_i} (units: Pa)}
#'   \item{b}{coefficient \eqn{b_i} (unitless)}
#' }
#' 
#' @source Text file: https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
#'     Rda file: \code{sigma_l137 <- generate_ecmwf_sigma("inst/l137.txt");
#'     devtools::use_data(sigma_l137, pkg = ".", internal = FALSE)}
"sigma_l137"

#' E3SM ne30pg2 grid
#'
#' A dataset containing the latitude, longitude, cell area, and land fraction of the E3SM ne30pg2 grid
#'
#' @format A data frame with 3 variables:
#' \describe{
#'   \item{ncol}{grid cell ("column") index}
#'   \item{lon}{grid cell midpoint longitude}
#'   \item{lat}{grid cell midpoint latitude}
#'   \item{area}{grid cell area}
#'   \item{LANDFRAC}{grid cell fractional land area}
#' }
#' 
#' @source Any E3SMv2 (other versions not tested) ne30pg2 output file
#'     containing a single timestep of the LANDFRAC variable (i.e.,
#'     nfilt=1 namelist option).  The specific file used is given
#'     below.
#'
#' Creation:
#' \code{nc <- nc_open("scratch//e3sm_scratch/v2_aer_nudged/v2_1850aer_nudged/run/v2_1850aer_nudged.eam.h0.2012-12.nc");
#'    ne30pg2 <- dplyr::full_join(plotutils::nc.to.df(nc, c("lon", "lat", "area")),
#'                                plotutils::nc.to.df(nc, c("LANDFRAC")) %>% dplyr::select(-time))}
#'
#' Use:
#'     devtools::use_data(ne30pg2, pkg = ".", internal = FALSE)}
"ne30pg2"

