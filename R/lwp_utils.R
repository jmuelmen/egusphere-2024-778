#' Load the fields required for the LWP-CDNC analysis
#'
#' @export
load_lwp <- function(fname, radiation = FALSE) {
    gc()
    nc <- nc_open(fname)
    ## nbdate <- ncdf4::ncvar_get(nc, "nbdate")
    nbdate <- ncatt_get(nc, "time")$"units"
    df.lonlat <- nc.to.df(nc, c("lon", "lat"))
    df.ttop <- nc.to.df(nc, c("ttop", "icc"))
    mask.ttop <- with(df.ttop, icc < 1e-3 & ttop > 273.15)
    rm(df.ttop)
    df.tend <- nc.to.df(nc, c("lwp", "cdnc")) %>%
        group_by(ncol) %>%
        transmute(time = time,
                  dlwp = lead(lwp) - lwp,
                  dcdnc = lead(cdnc) - cdnc) %>%
        ungroup()
    df.global.ttop <- nc.to.df(nc, grep(sprintf("area|lcc|lwp|cdnc|ttop|icc|PRECC|PRECL|OMEGA[57]00|TH7001000|LANDFRAC|rwp%s",
                                                if (radiation) "|FSNTOA|FSUTOAC|SOLIN" else ""),
                                        names(nc$var), value = TRUE), mask = mask.ttop)
    nc_close(nc)
    ## names(df) <- gsub("_210e_to_245e_15n_to_35n", "", names(df))
    df.global.ttop %>%
        left_join(df.lonlat, by = "ncol") %>%
        left_join(df.tend, by = c("ncol", "time")) %>%
        rename_with(tolower) %>%
        dplyr::mutate(time = as.POSIXct(nbdate, format = "days since %Y-%m-%d %H:%M:%S") + time * 86400) %>%
        mutate(is.scu = omega500 * 86400 > 10e2 & omega500 * 86400 > 10e2 & th7001000 > 18.55)
}

#' Load from GISS
#'
#' @export
load_lwp_giss <- function(fname) {
    gc()
    Sys.setenv(TZ = "UTC")
    nc <- nc_open(fname)
    ## nbdate <- ncdf4::ncvar_get(nc, "nbdate")
    nbdate <- ncatt_get(nc, "time")$"units"
    df.iwp <- nc.to.df(nc, c("iwp"))
    mask.iwp <- df.iwp$iwp == 0
    rm(df.iwp)
    ## df.tend <- nc.to.df(nc, c("lwp", "cdnc")) %>%
    ##     group_by(ncol) %>%
    ##     transmute(time = time,
    ##               dlwp = lead(lwp) - lwp,
    ##               dcdnc = lead(cdnc) - cdnc) %>%
    ##     ungroup()
    df.global.iwp <- nc.to.df(nc, grep("lwp|cLWPss|cldss_2d|ssct_ncl", names(nc$var), value = TRUE), mask = mask.iwp)
    nc_close(nc)
    df.global.iwp %>%
        filter(cldss_2d > 0) %>%
        mutate(cdnc_ic = ssct_ncl,
               lwp_ic = cLWPss / cldss_2d) %>%
        dplyr::mutate(time = as.POSIXct(nbdate, format = "hours since %Y-%m-%d %H:%M UTC") + time * 3600) %>%
        dplyr::filter(is.finite(cdnc_ic))
    ## ## names(df) <- gsub("_210e_to_245e_15n_to_35n", "", names(df))
    ## df.global.ttop %>%
    ##     left_join(df.lonlat, by = "ncol") %>%
    ##     left_join(df.tend, by = c("ncol", "time")) %>%
    ##     rename_with(tolower) %>%
    ##     dplyr::mutate(time = as.POSIXct(nbdate, format = "days since %Y-%m-%d %H:%M:%S") + time * 86400) %>%
    ##     mutate(is.scu = omega500 * 86400 > 10e2 & omega500 * 86400 > 10e2 & th7001000 > 18.55)
}

#' Load the fields required for the LWP-CDNC analysis (Andrew Gettelman's CAM run)
#'
#' @export
load_lwp_cam <- function(fname) {
    gc()
    cl <- snow::makeCluster(rep("localhost", 8), type = "SOCK", outfile = "/dev/null")
    on.exit({
        snow::stopCluster(cl)
    })
    doSNOW::registerDoSNOW(cl)
    snow::clusterExport(cl, "fname", environment())
    print(snow::clusterEvalQ(cl, ls()))
    snow::clusterEvalQ(cl, {
        library(ncdf4)
        library(plotutils)
        library(magrittr)
        library(dplyr)

        nc <- nc_open(fname)
    })
    print(snow::clusterEvalQ(cl, ls()))
    
    nc <- nc_open(fname)
    ## nbdate <- ncdf4::ncvar_get(nc, "nbdate")
    time <- ncvar_get(nc, "time")
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")
    nbdate <- ncatt_get(nc, "time")$"units"
    ## df.tend <- nc.to.df(nc, c("TGCLDLWP", "ACTNL")) %>%
    ##     group_by(lon, lat) %>%
    ##     rename(lwp = TGCLDLWP, cdnc = ACTNL) %>%
    ##     transmute(time = time,
    ##               dlwp = lead(lwp) - lwp,
    ##               dcdnc = lead(cdnc) - cdnc) %>%
    ##     ungroup()
    df.global.iwp <- plyr::ldply(1 : length(lat), function(step) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))
        df.global.iwp <- nc.to.df(nc, grep("^[A-Z]+", names(nc$var), value = TRUE),
                                  start = c(1, step, 1), count = c(-1, 1, -1)) %>%
            filter(TGCLDIWP == 0) %>%
            filter(ACTNL != 0) %>%
            select(-TGCLDIWP)
    }, .parallel = TRUE)
    str(df.global.iwp)
    gc()
    nc_close(nc)
    ## names(df) <- gsub("_210e_to_245e_15n_to_35n", "", names(df))
    df.global.iwp %>%
        ## left_join(df.tend, by = c("lon", "lat")) %>%
        ## rename_with(tolower) %>%
        dplyr::mutate(time = as.POSIXct(nbdate, format = "days since %Y-%m-%d %H:%M:%S") + time * 86400) %>%
        mutate(is.scu = NA)
}

#' @export
#' 
load_gfdl <- function() {
    readRDS("../am4.rds") %>%
        filter(is.finite(cldtop)) %>%
        filter(LWP > 0) %>%
        filter(IWP <= 0) %>%
        filter(ttop > 273) %>%
        mutate(lwp_ic = 1e2 * LWP / lclmodis) %>% ## until we can get the 2D cloud cover
        mutate(cdnc_ic = cdnc / cldtop) %>%
        mutate(lts = t700 * 0.7 ^ -(2/7) - ts,
               is.scu = lts > 18.55 & omega500 * 864 > 10 & omega700 * 864 > 10)
}

#' @export
#' 
load_aerocom <- function(path = "/global/homes/j/jmuelmen", pattern = ".*_IND3_lwp_cdnc.rds") {
    library(magrittr)
    library(dplyr)
    lf <- list.files(path, pattern, full.names = TRUE)
    plyr::ldply(lf, function(fname) {
        ## parallel::mccollect(parallel::mcparallel( {
        df <- readRDS(fname) %>%
            dplyr::filter(is.finite(cdnc),
                          is.finite(lwp),
                          is.finite(tcc),
                          cdnc > 0,
                          lwp > 0,
                          tcc > 0,
                          iwp < 1e-3) %>%
            dplyr::mutate(cdnc = ifelse(grepl("HadGEM3", model), cdnc * 1e6, cdnc),
                          cdnc_ic = ifelse(grepl("SPRINTARS", model), cdnc, cdnc / tcc),
                          lwp_ic = lwp / tcc) 
        df.cldclass <- dplyr::bind_rows(df %>%
                                        dplyr::mutate(class = "all"),
                                        df %>%
                                        dplyr::filter(tcc > 0.25) %>%
                                        dplyr::mutate(class = "cloudy"),
                                        df %>%
                                        dplyr::filter(tcc > 0.9) %>%
                                        dplyr::mutate(class = "ovc"))
        dplyr::bind_rows(df.cldclass %>%
                         dplyr::mutate(iwp.limit = 1e-3),
                         df.cldclass %>%
                         dplyr::filter(iwp == 0) %>%
                         dplyr::mutate(iwp.limit = 0)) %>%
            dplyr::group_by(model, class, iwp.limit) %>%
            lwp.cdnc.conditional(lwp.bins = exp(seq(log(1e-4), log(3), length.out = 25 * 2)),
                                 cdnc.bins = exp(seq(log(1e6), log(1000e6), length.out = 345)))
        ## }))[[1]]
    })
}

#' @export
#' 
filter_warm <- function(df) {
    df %>%
        dplyr::filter(icc < 1e-3, ttop > 273.15, lcc > 0.1) 
}

#' Like filter_warm but without the ttop requirement
#'
#' @export
filter_liquid <- function(df) {
    df %>%
        dplyr::filter(icc < 1e-3, lcc > 0.1) 
}

#' @export
#' 
filter_ocean <- function(df) {
    df %>%
        dplyr::filter(landfrac < 1e-3) 
}

#' Apply all filters used in the standard analysis (warm cloud, ocean, Sc dynamical regime)
#' @export
filter_standard <- function(df) {
    df %>%
        filter_warm() %>%
        filter_ocean() %>%
        filter_sc()
}

#' Select standard variables for L-N analysis
#' @export
select_standard <- function(df) {
    dplyr::select(df,
                  dplyr::any_of(c("lon", "lat", "time",
                                  "is.scu", "period", "model",
                                  "lcc", "lwp_ic", "cdnc_ic"))) 
}

#' @export
#' 
calculate_incloud <- function(df) {
    df %>%
        dplyr::mutate(lwp_ic = lwp / lcc,
                      cdnc_ic = cdnc / lcc)
}

#' @export
#' 
calculate_albedo_from_lcc <- function(df) {
    ## ## probably the cleanest way to do temp variables for efficiency:
    ## solin_recip = with(df, ifelse(solin == 0, NA, 1 / solin))
    df %>%
        dplyr::mutate(fsutoa = solin - fsntoa,
                      fsutoa_cld = (fsutoa - (1 - lcc) * fsutoac) / lcc,
                      alpha_cld = fsutoa_cld / solin)
}

#' @export
#' 
discretize_lwp <- function(df, bins = exp(seq(log(0.01), log(0.3), length.out = 25 * 2))) {
    df %>%
        plotutils::discretize(lwp_ic, bins) 
}

#' @export
#' 
discretize_cdnc <- function(df, bins = exp(seq(log(1e6), log(300e6), length.out = 37 * 2))) {
    df %>%
        plotutils::discretize(cdnc_ic, bins)
}

#' @export
#' 
summarize_lwp <- function(df, ...) {
    df %>%
        dplyr::group_by(cdnc_ic, .add = TRUE) %>%
        dplyr::summarize(lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE)),
                         ...) %>%
        dplyr::ungroup(cdnc_ic)
}

#' @export
#' 
summarize_pcond <- function(df) {
    df %>%
        dplyr::group_by(cdnc_ic, lwp_ic, .add = TRUE) %>%
        dplyr::summarize(p = n()) %>%
        ungroup(lwp_ic) %>%
        mutate(p = p / sum(p)) %>%
        dplyr::ungroup(cdnc_ic)
}

#' @export
#' 
summarize_weighted_mean <- function(df) {
    df %>%
        group_by(cdnc_ic, .add = TRUE) %>%
        mutate(lwp_ic_bar = exp(sum(log(lwp_ic) * count.lwp, na.rm = TRUE) /
                                sum(is.finite(lwp_ic) * count.lwp, na.rm = TRUE))) %>%
        ungroup(cdnc_ic) %>%
        group_by(lwp_ic, .add = TRUE) %>%
        mutate(cdnc_ic_bar = exp(sum(log(cdnc_ic) * count.cdnc, na.rm = TRUE) /
                                sum(is.finite(cdnc_ic) * count.cdnc, na.rm = TRUE))) %>%
        ungroup(lwp_ic) %>%
        ## ## note: when 2D histo contains empty CDNC-LWP bins, those rows/columns will then be missing elements of the CDNC or LWP weighted mean
        ## ## note 2: since there are usually not many missing bins, taking the median of the weighted means will usually be fine; better fix might be to discard the estimate if the vector has missing elements; even better would be to calculate the global means in lwp.cdnc.conditional(), but that would require rerunning
        mutate(cdnc_ic_bar = median(cdnc_ic_bar),
               lwp_ic_bar = median(lwp_ic_bar)) 
}

#' Parameterization of the "inverted v" Nd-LWP relationship
#'
#' Note: the parameterization is piecewise linear in log-log space; in
#' the notation below, the log is implicit
#'
#' @param N really log(N)
#' @param Np really the log of the location of the peak of the inverted v
#' @param Lp really the log of the location of the peak of the inverted v
#' @param ml the slope below the peak
#' @param mh the slope above the peak
#'
#' @return log(L)
#'
#' @export
inverted.v <- function(N, Np, Lp, ml, mh) {
    ifelse(N < Np, Lp + ml * (N - Np), Lp + mh * (N - Np))
}

#' @export
#' 
logify <- function(df) {
    df %>%
        dplyr::transmute(N = log(cdnc_ic),
                         L = log(lwp_ic)) 
}

#' @export
#' 
unlogify <- function(df) {
    df %>%
        dplyr::mutate(Np = exp(Np),
                      Lp = exp(Lp))
}

#' @export
#' 
fit_inverted.v <- function(df, summarize = TRUE) {
    fit <- df %>%
        logify() %>%
        stats::nls(L ~ inverted.v(N, Np, Lp, ml, mh), .,
                   start = list(Np = log(1e7), Lp = log(0.1), ml = 0.5, mh = -0.5))
    if (summarize) {
        summ <- fit %>%
            summary() %$%
            coefficients[,1] %>%
            rbind() %>%
            transform() %>%
            unlogify()
    } else {
        fit
    }
}

#' @export
#' 
summarize_inverted.v <- function(df) {
    df %>%
        dplyr::summarize(Nmax = cdnc_ic[which.max(lwp_ic)],
                         Lmax = max(lwp_ic),
                         L100 = lwp_ic[which.min(abs(cdnc_ic - 1e8))],
                         Lratio = L100 / Lmax)
}

#' Use inverted v fit to predict LWP from Nd
#' 
#' @export
#' 
predict_inverted.v <- function(df.train, df.apply) {
    fit.v <- fit_inverted.v(df.train, summarize = FALSE)
    df.log <- df.apply %>% logify()
    df.apply %>% dplyr::mutate(log.err = predict(fit.v, newdata = df.log) - df.log$L,
                               lwp_ic.pred = exp(predict(fit.v, newdata = df.log)))
}

#' Use full P(LWP | Nd) to predict LWP from Nd
#'
#' @param conditional.train A data.frame produced by
#'     lwp.cdnc.conditional() that constitutes the P(L|N) to be used
#'     for prediction
#' @param conditional.apply A data.frame produced by
#'     lwp.cdnc.conditional that constitutes the perturbed P(N) to
#'     which the P(L|N) is to be applied
#'
#' Note: the conditionals shall have idential cdnc_ic bins
#' 
#' Note: the conditionals shall not be grouped
#' 
#' @export
predict_lwp_cdnc_conditional <- function(conditional.train, conditional.apply) {
    ## summarizing CDNC counts will produce the wrong answer if there are any groups
    stopifnot(length(groups(conditional.train)) == 0)
    stopifnot(length(groups(conditional.apply)) == 0)

    ## count CDNC in lieu of calculating perturbed P(N)
    df.cdnc.apply <- conditional.apply %>%
        dplyr::group_by(cdnc_ic) %>%
        dplyr::summarize(count.cdnc = count.cdnc[1]) %>%
        dplyr::ungroup()

    df.lwp.pred <- plyr::ddply(df.cdnc.apply, ~cdnc_ic, function(df) {
        df.lwp <- conditional.train %>%
            filter(cdnc_ic == df$cdnc_ic) %>%
            mutate(p = p * df$count.cdnc)
    }, .progress = "text") %>%
        group_by(lwp_ic) %>%
        summarize(p = sum(p)) %>%
        ungroup() %>%
        mutate(p = p / sum(p))
}

#' Parameterization of the "inverted v" Nd-LWP relationship
#'
#' Note: the parameterization is piecewise linear in log-log space; in
#' the notation below, the log is implicit
#'
#' @param N really log(N)
#' @param Np really the log of the location of the peak of the inverted v
#' @param Lp really the log of the location of the peak of the inverted v
#' @param ml the slope below the peak
#' @param mh the slope above the peak
#'
#' @return log(L)
#'
#' @export
inverted_v_path_data <- function(Np, Lp, ml, mh, Nmin = 1e6, Nmax = 3e8, n.points = 2) {
    ret <- data.frame(cdnc_ic = c(seq(Nmin, Np, length.out = n.points),
                                  seq(Np, Nmax, length.out = n.points)[-1])) ## don't visit Np twice
    ret %>%
        mutate(lwp_ic = exp(inverted.v(log(cdnc_ic), log(Np), log(Lp), ml, mh)))
}

#' @export
#' 
lwp.cdnc.conditional <- function(df, cdnc.bins, lwp.bins) {
    if (missing(cdnc.bins)) {
        df %<>% discretize_cdnc()
    } else {
        df %<>% discretize_cdnc(cdnc.bins)
    }
    df.cdnc <- df %>% summarize_lwp(count.cdnc = n())
    if (missing(lwp.bins)) {
        df.lwp.cdnc <- df %>% discretize_lwp()
    } else {
        df.lwp.cdnc <- df %>% discretize_lwp(lwp.bins)
    }
    df.lwp <- df.lwp.cdnc %>% 
        dplyr::group_by(lwp_ic, .add = TRUE) %>%
        dplyr::summarize(count.lwp = n()) %>%
        dplyr::ungroup(lwp_ic)
    df.lwp.cdnc %<>% summarize_pcond()
    full_join(df.lwp.cdnc, df.lwp) %>%
        full_join(df.cdnc %>% rename(lwp_ic_mean = lwp_ic))
}

#' @export
#'
lwp.cdnc.plot <- function(df, heatmap = FALSE, marginal = FALSE,
                          labels = labs(x = "$N_d$ (in cloud, m$^{-3}$)",
                                        y = "$\\mathcal{L}$ (in cloud, kg~m$^{-2}$)"),
                          legend.position = c(0.4, 0.25),
                          cdnc.density.limits = c(0, 0.01),
                          expand = TRUE,
                          ...) {
    g <- df %>%
        group_by(cdnc_ic, .add = TRUE) %>%
        summarize(lwp_ic = lwp_ic_mean[1]) %>%
        ungroup(cdnc_ic) %>%
        ggplot(aes(x = cdnc_ic, y = lwp_ic, ...))

    if (heatmap) {
        g <- g + geom_raster(aes(fill = p), df)
    }

    g <- g +
        geom_line() +
        scale_x_log10(limits = range(df$cdnc_ic, na.rm = TRUE)) +
        scale_y_log10(limits = if (marginal) range(df$lwp_ic, na.rm = TRUE) else NULL) +
        labels +
        coord_cartesian(expand = expand) +
        theme(legend.position = legend.position,
              legend.background = element_rect(fill = NA),
              legend.spacing.y = unit(0, "lines"),
              legend.direction = "vertical", legend.box = "horizontal") 
    
    if (marginal) {
        g.lwp <- df %>%
            filter(is.finite(lwp_ic)) %>%
            group_by(lwp_ic, .add = TRUE) %>%
            summarize(count.lwp = count.lwp[1]) %>%
            ungroup(lwp_ic) %>%
            mutate(p = count.lwp / sum(count.lwp)) %>%
            ggplot(aes(x = lwp_ic, y = p, ...)) +
            geom_line() +
            ## coord_flip(xlim = range(df$lwp_ic_mean)) +
            coord_flip(expand = expand) +
            scale_x_log10() +
            theme_void() +
            theme(legend.position="none")
        g.cdnc <- df %>%
            filter(is.finite(cdnc_ic)) %>%
            group_by(cdnc_ic, .add = TRUE) %>%
            summarize(count.cdnc = count.cdnc[1]) %>%
            ungroup(cdnc_ic) %>%
            mutate(p = count.cdnc / sum(count.cdnc)) %>%
            ggplot(aes(x = cdnc_ic, y = p, ...)) +
            geom_line() +
            scale_x_log10() +
            theme_void() +
            coord_cartesian(ylim = cdnc.density.limits, expand = expand) +
            theme(legend.position="none")
        wrap_plots(list(marginal.cdnc = g.cdnc, plot_spacer(), correlation = g, marginal.lwp = g.lwp),
                   widths = c(1, 0.1), heights = c(0.1, 1),
                   ncol = 2, nrow = 2, byrow = TRUE)
    } else {
        g
    }
    
}

