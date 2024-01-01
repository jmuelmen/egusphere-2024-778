#' @export
generate_ecmwf_sigma <- function(fname = "inst/l137.txt") {
    read.table(file = fname, header = TRUE, na.strings = "-") %>%
        dplyr::select(lev = n, a = a..Pa., b) 
}
