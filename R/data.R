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

