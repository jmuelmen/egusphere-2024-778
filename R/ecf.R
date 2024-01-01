#' @export
#' 
ecf.standard <- function() {
    c("qflx",
      "sh.cp",
      "delta.F.cp",
      "log10.delta.R",
      "h",
      "cdnc_ic",
      "dthetacore.0",
      "dqcore.0",
      "delta.theta",
      "delta.q",
      "omega.pbl")
}

#' @export
#' 
ecf.cdnc <- function() {
    
}

#' @export
#' 
ecf.precip <- function() {
    
}

#' @export
#' 
ecf.top <- function() {
    
}

#' @export
#' 
ecf.plot <- function(df,
                     E.group = c("E_theta", "E_q"),
                     ecf.group = ecf.standard(),
                     quantile.filter = 5e-3,
                     median.line = TRUE,
                     raster = TRUE,
                     ...) {
    df.ecf <- df %>%
        select(any_of(c(E.group, ecf.group, "f.minus"))) %>%
        tidyr::gather(E, E_val, dplyr::any_of(E.group)) %>%
        group_by(E, .add = TRUE) %>%
        filter(between_quantiles(E_val, quantile.filter, 1 - quantile.filter)) %>%
        ungroup(E) %>%
        gather(ecf, ecf_val, dplyr::any_of(ecf.group)) %>%
        mutate(ecf = factor(ecf, unique(ecf))) %>%
        filter(ifelse(grepl("_ic$", ecf), f.minus > 0.1, TRUE)) %>%
        ## filter(ifelse(ecf == "delta.R", ecf_val > 1e-5 / 86400 , TRUE)) %>%
        group_by(ecf, .add = TRUE) %>%
        filter(between_quantiles(ecf_val, 0.02, 0.98)) %>%
        ungroup(ecf) %>%
        mutate(ecf = revalue(ecf, tikz_replacements_unitful()),
               E = revalue(E, tikz_replacements_unitful()))

    g <- df.ecf %>%
        ggplot(aes(x = ecf_val, y = E_val, ...)) +
        ## geom_point(size = 0) +
        ## stat_bin2d(bins = 100, geom = "raster") +
        ## geom_smooth() +
        scale_fill_distiller(palette = "Greys", trans = "log",
                             direction = 1, breaks = c(1,10,100)) +
        facet_grid(E ~ ecf, scales = "free_x") +
        scale_x_continuous(labels = tikz_sanitize_sparse) +
        labs(x = "", y = "")

    if (raster) {
        g <- g + stat_bin2d(bins = c(100, 100), geom = "raster") 
    }
    
    if (median.line) {
        ext.groups <- groups(df.ecf)
        print(c("E", "ecf", as.character(ext.groups)))
        df.ecf %>%
            ddply(c("E", "ecf", as.character(groups(df.ecf))), discretize, ecf_val, 25, equal_contents = TRUE) %>%
            str
        median.ecf <- df.ecf %>%
            ddply(c("E", "ecf", as.character(ext.groups)), discretize, ecf_val, 25, equal_contents = TRUE) %>%
            group_by_(as.character(ext.groups)) %>%
            group_by(E, ecf, ecf_val, .add = TRUE) %>%
            summarize(E_val = median(E_val)) %>%
            ungroup(E, ecf, ecf_val)
        str(median.ecf)
    
        g <- g + geom_line(data = median.ecf)
    } else {
        g <- g + geom_smooth()
    }

    return(g)
}
