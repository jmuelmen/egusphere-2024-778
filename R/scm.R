#' @export
#' 
load_scm <- function(fname) {
    nc <- ncdf4::nc_open(fname)
    ## calculate 3d pressure from hybrid vertical coordinate
    ps <- ncdf4::ncvar_get(nc, "PS")
    p0 <- ncdf4::ncvar_get(nc, "P0")
    hyam <- ncdf4::ncvar_get(nc, "hyam")
    hybm <- ncdf4::ncvar_get(nc, "hybm")
    dhyai <- diff(ncdf4::ncvar_get(nc, "hyai"))
    dhybi <- diff(ncdf4::ncvar_get(nc, "hybi"))
    p.3d <- aperm(outer(ps, hybm) + outer(array(p0, dim(ps)), hyam), c(2,1))
    dp.3d <- aperm(outer(ps, dhybi) + outer(array(p0, dim(ps)), dhyai), c(2,1))

    df.3d <- plotutils::nc.to.df(nc, grep("^(T|Q|CLDLIQ|QRL|QRS|OMEGA|OMEGAT|DCQ|DTCOND|DTCORE|DQCORE|VT|Z3|CCN3|CCN4|NUMLIQ|CLOUD)$",
                                          names(nc$var), value = TRUE)) %>%
        dplyr::mutate(QRS = ifelse(is.finite(QRS), QRS, 0)) ## QRS is NA if shortwave radiative transfer was not run
    df.2d <- plotutils::nc.to.df(nc, grep("^(SHFLX|LHFLX|QFLX|PRECL|PRECC|CLDLOW|cdnc|lcc|lwp)$",
                                          names(nc$var), value = TRUE))
    df.3d %<>% dplyr::left_join(df.2d, by = c("time")) %>%
        dplyr::rename_with(tolower) %>%
        dplyr::rename(qv = q,
                      ql = cldliq) %>%
        dplyr::group_by(time) %>%
        dplyr::mutate(jk = 1:length(lev)) %>%
        dplyr::ungroup()

    df.3d %<>% 
        dplyr::mutate(p = as.vector(p.3d),
                      dp = as.vector(dp.3d),
                      p0 = as.vector(p0)) %>%
        dplyr::mutate(prec = precc + precl) %>%
        calculate_thermodynamics()

    ncdf4::nc_close(nc)

    return(df.3d)
}

#' @export
#' 
load_scm_giss <- function(fname) {
    nc <- ncdf4::nc_open(fname)
    df.3d <- plotutils::nc.to.df(nc, grep("^(p_3d|t|th|q|qcl|dth_rad|omega|z|nclic|cf|ccn0p2)$", names(nc$var), value = TRUE))
    df.2d <- plotutils::nc.to.df(nc, c("shflx", "lhflx", "prec", "prsurf", "lwp", "cLWPss", "cldss_2d", "ssct_ncl")) %>%
        dplyr::mutate(qflx = lhflx / Lv)
    df.3d %<>% dplyr::left_join(df.2d, by = c("time", "lon", "lat")) %>%
        dplyr::mutate(time = time / 24) %>%
        dplyr::arrange(time, p) %>%
        dplyr::group_by(time) %>%
        dplyr::mutate(jk = 1 : length(p)) %>%
        dplyr::mutate(dp = 100 * prsurf / 1000 * (p - lag(p))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(dtcore = 0,
                      dqcore = 0)

    ncdf4::nc_close(nc)

    return(df.3d)
}

#' @export
#' 
scm_budget <- function(df.3d, inv.type = "jk.tinv") {
    df.budget <- df.3d %>%
        calculate_thermodynamics() %>%
        group_by(time) %>%
        find_inversion(inv.type) %>%
        calculate_h() %>%
        group_by(lev) %>%
        calculate_tendencies() %>%
        mutate(dh.dt_lagged = (h - lag(h)) / dt) %>%
        transform_dT_dtheta() %>%
        filter_invalid_qrad() %>%
        group_by(time) %>%
        calculate_budgets() %>%
        ungroup()
}
