nc_getstep <- function(nc, varname, stride, i, fname = varname) {
    nc <- nc[[fname]]
    ndm1 <- nc$var[[varname]]$ndims - 1 ## number of non-time dimensions
    if (length(ndm1) == 0)
        ndm1 <- 0
    ncdf4::ncvar_get(nc, varname,
                     c(rep(1, ndm1), (i - 1) * stride + 1),
                     c(rep(-1, ndm1), stride),
                     collapse_degen = FALSE)
}

#' Calculate horizontal scalar advection from u and v
#'
#' @param s Scalar quantity (2- or higher-dimensional array)
#' @param u Zonal wind component (2- or higher-dimensional array,
#'     dimensionality matching scalar)
#' @param v Meridional wind component
#' @param dx Zonal grid spacing in units of m; two-dimensional (lon-lat)
#' @param dy Meridional grid spacing in units of m; one-dimensional (lat)
#'
#' @export
advection <- function(s, u, v, dx, dy) {
    dsx <- plyr::aaply(s, 2:length(dim(s)),
                       function(x) as.vector(stats::filter(x, filter = 1:-1, circular = TRUE))) %>%
        aperm(c(length(dim(s)), 1:(length(dim(s)) - 1)))
    dsy <- plyr::aaply(s, (1:length(dim(s)))[-2],
                       function(x) as.vector(stats::filter(x, filter = 1:-1))) %>%
        aperm(c(1, length(dim(s)), 2:(length(dim(s)) - 1)))
    ## extend dx and dy to match dimensions of s
    dx <- outer(dx, array(1, dim(s)[-(1:2)]))
    dy <- outer(outer(array(1, dim(s)[1]), dy),
                array(1, dim(s)[-(1:2)]))
    return(u * dsx/dx + v * dsy/dy)
}

process.omega <- function() {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    doParallel::registerDoParallel(12)
    
    
    lf = list.files("/global/cscratch1/sd/plma/e3sm_scratch/EAMv1_test/latlon/", "EAMv1_test.cam.h1.2013-.*.nc", full.names = TRUE)
    
    arr = laply(lf, function(fname) { nc = nc_open(fname); on.exit(nc_close(nc)); ncvar_get(nc, "OMEGA", start = c(26, 26, 53, 1), count = c(1,1,1,-1)); }, .parallel = TRUE)
    
    df.omega <- data.frame(time = as.POSIXct("2013-01-01") + 3600 * 1 : (365 * 24), omega.700 = as.vector(t(arr)) * 864) %>%
        (function(df) {
            ldply(c(1, 12, 24, 72, 168), function(fil) {
                df %>%
                    mutate(omega.700 = stats::filter(omega.700, rep(1 / fil, fil)) %>%
                               as.vector,
                           filter = fil)
            })
        })(.)

    saveRDS(df.omega, "~/omega.rds")
    
    df.omega %>%
        ## filter(time >= "2013-06-01", time < "2013-08-01") %>%
        ggplot(aes(time, omega.700, color = factor(filter))) +
        geom_line() +
        scale_color_brewer("Filter (h)"## , palette = "Dark2"
                           ) +
        geom_hline(yintercept = 10, lty = "dashed", col = "grey") +
        labs(x = "", y = "hPa d$^{-1}$") +
        scale_x_datetime(date_labels = "%Y--%m--%d") +
        theme_bw()
}

process.pblh <- function() {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    doParallel::registerDoParallel(12)
    
    lf <- list.files("/global/cscratch1/sd/plma/e3sm_scratch/EAMv1_test/latlon/",
                    "EAMv1_test.cam.h1.2013-.*.nc",
                    full.names = TRUE)
    
    ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        ## fields required to calculate 3D pressure field
        ps <- ncvar_get(nc, "PS", start = c(26, 26, 1), count = c(1,1,-1))
        p0 <- ncdf4::ncvar_get(nc, "P0")
        hyam <- ncdf4::ncvar_get(nc, "hyam")
        hybm <- ncdf4::ncvar_get(nc, "hybm")

        p.3d <- aperm(outer(ps, hybm) + outer(array(p0, dim(ps)), hyam),
                      c(2,1))

        ## other 2D fields
        pblh <- ncvar_get(nc, "PBLH", start = c(26, 26, 1), count = c(1,1,-1))
        precc <- ncvar_get(nc, "PRECC", start = c(26, 26, 1), count = c(1,1,-1))
        precl <- ncvar_get(nc, "PRECL", start = c(26, 26, 1), count = c(1,1,-1))
        qflx <- ncvar_get(nc, "QFLX", start = c(26, 26, 1), count = c(1,1,-1))
        shflx <- ncvar_get(nc, "SHFLX", start = c(26, 26, 1), count = c(1,1,-1))

        ## other 3D fields
        z3 <- ncvar_get(nc, "Z3", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        t <- ncvar_get(nc, "T", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        qv <- ncvar_get(nc, "Q", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        ql <- ncvar_get(nc, "CLDLIQ", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        omega <- ncvar_get(nc, "OMEGA", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        prao <- ncvar_get(nc, "PRAO", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        prco <- ncvar_get(nc, "PRCO", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        qrl <- ncvar_get(nc, "QRL", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        qrs <- ncvar_get(nc, "QRS", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        cloud <- ncvar_get(nc, "CLOUD", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        cloudfrac_clubb <- ncvar_get(nc, "CLOUDFRAC_CLUBB", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        numliq <- ncvar_get(nc, "NUMLIQ", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))

        ## 3D fields on full levels
        wp2_clubb <- ncvar_get(nc, "WP2_CLUBB", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        wp3_clubb <- ncvar_get(nc, "WP3_CLUBB", start = c(26, 26, 1, 1), count = c(1,1,-1,-1))
        
        df.3d <- expand.grid(lev = 100 * ncvar_get(nc, "lev") %>% as.vector,
                             time = as.POSIXct("2013-01-01") + 86400 * ncvar_get(nc, "time") %>% as.vector) %>%
            mutate(z3   = as.vector(z3),
                   p    = as.vector(p.3d),
                   qv   = as.vector(qv),
                   ql   = as.vector(ql),
                   omega = as.vector(omega),
                   prao =  as.vector(prao),
                   prco =  as.vector(prco),
                   qrl =  as.vector(qrl),
                   qrs =  as.vector(qrs),
                   cloud = as.vector(cloud),
                   cloudfrac_clubb = as.vector(cloudfrac_clubb),
                   numliq = as.vector(numliq),
                   t    = as.vector(t)) %>%
            mutate(theta = t * (p0 / p) ^ (2/7),
                   theta_l = theta - 2.5e3 * ql,
                   qt = qv + ql) %>%
            ## total kluge, but hopefully level numbers will be static for a while
            group_by(time) %>%
            mutate(jk = 1:72) %>%
            ungroup()

        df.3d.ilev <- expand.grid(ilev = 100 * ncvar_get(nc, "ilev") %>% as.vector,
                             time = as.POSIXct("2013-01-01") + 86400 * ncvar_get(nc, "time") %>% as.vector) %>%
            mutate(wp2_clubb   = as.vector(wp2_clubb),
                   wp3_clubb   = as.vector(wp3_clubb)) %>%
            ## total kluge, but hopefully level numbers will be static for a while
            group_by(time) %>%
            mutate(jk = 0:72 + 0.5) %>%
            ungroup()
        
        df.2d <- expand.grid(time = as.POSIXct("2013-01-01") + 86400 * ncvar_get(nc, "time") %>% as.vector) %>%
            mutate(pblh         = as.vector(pblh),
                   qflx =  as.vector(qflx),
                   shflx =  as.vector(shflx),
                   precc        = as.vector(precc),
                   precl        = as.vector(precl))

        df <- bind_rows(df.3d, df.3d.ilev) %>%
            left_join(df.2d, by = "time")
    }, .parallel = TRUE) %>%
    saveRDS("profiles.rds")
}

#' Take a look at different PBL height definitions
#'
#' @export
process.pblh.inversion <- function(path = "/qfs/people/muel306/scratch/e3sm_scratch/EAMv1_entrain_NEP_compy_nudged/run",
                                   pattern = "EAMv1_entrain_.*.cam.h1.20.0-.*.nc") {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    doParallel::registerDoParallel(12)

    Sys.setenv(TZ = "UTC")

    lf <- list.files(path, pattern,
                     full.names = TRUE)

    ## print(lf)

    plyr::ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        df.1d <- plotutils::nc.to.df(nc, names(nc$var)[grepl("lon|lat", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))
        
        ## 3D fields
        df.3d <- plotutils::nc.to.df(nc,
                                     names(nc$var)[grepl("T|Q|CLDLIQ", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        df.2d <- df.3d %>%
            dplyr::group_by(ncol, time) %>%
            dplyr::mutate(jk = 1 : length(lev)) %>%
            dplyr::filter(lev > 600) %>%
            dplyr::mutate(theta = T * (1000 / lev) ^ (2/7),
                          qt = Q + CLDLIQ) %>%
            dplyr::summarize(jk.tinv = max(which(diff(T) < 0)),
                             jk.theta = which.min(diff(theta) / diff(lev)), ## theta jump is negative in increasing lev direction
                             jk.qt = which.max(diff(qt) / diff(lev)), ## qt jump is positive in increasing lev direction
                             theta.jump = diff(theta)[jk.theta],
                             qt.jump = diff(qt)[jk.qt]) %>% 
            dplyr::ungroup()

        ## print(summary(df.2d))
        
        df.2d %<>%
            dplyr::left_join(df.1d, by = "ncol")
        
    }, .parallel = TRUE) %>%
    saveRDS("pblh.rds")
 
}

#' Pluck the entrainment debugging variables from the grid point with
#' max Sc occurrence from the EAMv1_entrain_NEP run
#'
#' By default, only July, August, September are used
#'
#' In some runs, the correct pattern is ...cam.h2...
#'
#' @export
process.NEP.entrain <- function(path = "/qfs/people/muel306/scratch/e3sm_scratch/EAMv1_entrain_NEP_compy_nudged/run",
                                pattern = "EAMv1_entrain_.*.cam.h1.20.0-0[7-9]-.*.nc",
                                icol = 291) {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    doParallel::registerDoParallel(12)
    
    Sys.setenv(TZ = "UTC")

    lf <- list.files(path, pattern,
                     full.names = TRUE)

    ## print(lf)

    fname.out <- sprintf("profiles_NEP_entrain_run-%s_icol-%s.rds",
                         ## run name is the second-to-last component of the path (the last is "/run/")
                         rev(strsplit(path, "/")[[1]])[2],
                         ## icol could be a vector, so use the unparsed expression
                         deparse(substitute(icol)))

    print(sprintf("fname.out: %s", fname.out))

    plyr::ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        df.entrain <- plotutils::nc.to.df(nc,
                                          grep("^ENTRAIN_.*",
                                               names(nc$var), value = TRUE)) %>%
            dplyr::rename_with(. %>% gsub("_[0-9]+e_to.*$", "", .)) %>%
            dplyr::filter(ncol %in% icol)

        df.entrain.2d <- df.entrain %>%
            dplyr::filter(is.na(lev)) %>%
            dplyr::select(-lev) %>%
            tidyr::gather(var, val, -c(ncol, time)) %>%
            dplyr::group_by(var) %>%
            dplyr::filter(!all(is.na(val))) %>%
            dplyr::ungroup() %>%
            tidyr::spread(var, val)

        ## str(df.entrain.2d)
        
        df.entrain.3d <- df.entrain %>%
            dplyr::filter(!is.na(lev)) %>%
            tidyr::gather(var, val, -c(ncol, time, lev)) %>%
            dplyr::group_by(var) %>%
            dplyr::filter(!all(is.na(val))) %>%
            dplyr::ungroup() %>%
            tidyr::spread(var, val)

        ## str(df.entrain.3d)

        df.1d <- plotutils::nc.to.df(nc, names(nc$var)[grepl("lon|lat", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        nbdate <- ncdf4::ncvar_get(nc, "nbdate")
        p0 <- ncdf4::ncvar_get(nc, "P0")
        df.hym <- plotutils::nc.to.df(nc, c("hyam", "hybm"))
        df.dhyi <- with(plotutils::nc.to.df(nc, c("hyai", "hybi")),
                        ## layer thicknesses
                        data.frame(lev = df.hym$lev,
                                   dhyai = diff(hyai),
                                   dhybi = diff(hybi)))

        ## 2D fields
        df.2d <- plotutils::nc.to.df(nc,
                                     grep("^PRECC_|^PRECL_|^QFLX_|^SHFLX_|^OMEGA500_|^OMEGA700_|^TH7001000_|^PS_",
                                          names(nc$var), value = TRUE)) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .)) %>%
            dplyr::filter(ncol %in% icol) 

        ## 3D fields
        df.3d <- plotutils::nc.to.df(nc,
                                     grep("^Z3_|^T_|^Q_|^CLDLIQ_|^OMEGA_|^PRAO_|^PRCO_|^QRL_|^QRS_|^CLOUD_|^NUMLIQ_|^CCN4_|^DTCORE_|^DQCORE_|^DTCOND_",
                                          names(nc$var), value = TRUE)) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .)) %>%
            dplyr::filter(ncol %in% icol) 

        ret <- df.3d %>%
            dplyr::left_join(df.2d, by = c("time", "ncol")) %>%
            dplyr::left_join(df.entrain.3d, by = c("time", "ncol", "lev")) %>%
            dplyr::left_join(df.entrain.2d, by = c("time", "ncol")) %>%
            dplyr::left_join(df.1d, by = "ncol") %>%
            dplyr::left_join(df.hym, by = "lev") %>%
            dplyr::left_join(df.dhyi, by = "lev") %>%
            dplyr::mutate(p0 = p0) %>%
            dplyr::rename_with(tolower) %>%
            dplyr::mutate(time = as.POSIXct(sprintf("%d", nbdate), format = "%Y%m%d") + time * 86400)

        str(ret)

        ret
    }, .parallel = TRUE) %>%
    saveRDS(fname.out)
}

#' Extract the fields needed for the entrainment calculation from all
#' columns satisfying the uniform PBL height among nearest neighbors
#' condition
#'
#' @export
process.NEP.uniform.neighbors <- function(path = "/qfs/people/muel306/scratch/e3sm_scratch/EAMv1_entrain_NEP_compy_nudged/run",
                                          pattern = "EAMv1_entrain_.*.cam.h1.20.0-.*.nc",
                                          icol.select = 237:346) {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    doParallel::registerDoParallel(12)
    
    Sys.setenv(TZ = "UTC")

    lf <- list.files(path, pattern,
                     full.names = TRUE)

    ## print(lf)

    plyr::ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        df.1d <- plotutils::nc.to.df(nc, names(nc$var)[grepl("lon|lat", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        ## find nearest neighbors
        df.neighbors <- df.1d %>%
            select(ncol, lon, lat) %>%
            neighbors()

        ## find PBL heights
        df.pblh <- plotutils::nc.to.df(nc,
                                       names(nc$var)[grepl("^T_|^Q_|^CLDLIQ_", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .)) %>%
            dplyr::group_by(ncol, time) %>%
            dplyr::mutate(jk = 1 : length(lev)) %>%
            dplyr::filter(lev > 600) %>%
            dplyr::mutate(theta = T * (1000 / lev) ^ (2/7),
                          qt = Q + CLDLIQ) %>%
            dplyr::summarize(jk.tinv = max(which(diff(T) < 0)),
                             jk.theta = which.min(diff(theta) / diff(lev)), ## theta jump is negative in increasing lev direction
                             jk.qt = which.max(diff(qt) / diff(lev)), ## qt jump is positive in increasing lev direction
                             theta.jump = diff(theta)[jk.theta],
                             qt.jump = diff(qt)[jk.qt]) %>% 
            dplyr::ungroup() %>%
            dplyr::left_join(df.1d, by = "ncol")

        ## select cases of uniform PBL height in all nearest neighbors
        df.uniform <- df.pblh %>%
            dplyr::select(ncol, lon, lat, time, jk.tinv : jk.qt) %>%
            tidyr::gather(type, jk, jk.tinv : jk.qt) %>%
            dplyr::filter(type == "jk.tinv") %>%  ## remove this filter to get other inversion types
            plyr::ddply(~ type, function(x) {
                dplyr::inner_join(x, df.neighbors, by = c("ncol" = "ncol1")) %>%
                    dplyr::left_join(x %>% select(ncol, time, jk),
                                     by = c("ncol2" = "ncol", "time")) %>%
                    dplyr::group_by(ncol, lon, lat, time) %>%
                    dplyr::filter(is.finite(jk.x) & all(jk.x == jk.y)) %>%
                    dplyr::group_by(ncol, lon, lat, time) %>%
                    dplyr::summarize(jk = jk.x[1]) %>%
                    dplyr::ungroup()
            })

        nbdate <- ncdf4::ncvar_get(nc, "nbdate")
        p0 <- ncdf4::ncvar_get(nc, "P0")
        df.hym <- plotutils::nc.to.df(nc, c("hyam", "hybm"))
        df.dhyi <- with(plotutils::nc.to.df(nc, c("hyai", "hybi")),
                        ## layer thicknesses
                        data.frame(lev = df.hym$lev,
                                   dhyai = diff(hyai),
                                   dhybi = diff(hybi)))

        ## 2D fields
        df.2d <- plotutils::nc.to.df(nc,
                                     grep("^PRECC_|^PRECL_|^QFLX_|^SHFLX_|^OMEGA500_|^OMEGA700_|^TH7001000_|^PS_",
                                          names(nc$var), value = TRUE)) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        ## 3D fields
        df.3d <- plotutils::nc.to.df(nc,
                                     grep("^Z3_|^T_|^Q_|^CLDLIQ_|^OMEGA_|^PRAO_|^PRCO_|^QRL_|^QRS_|^CLOUD_|^NUMLIQ_|^CCN4_|^DTCORE_|^DQCORE_|^DTCOND_",
                                          names(nc$var), value = TRUE)) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        ret <- df.3d %>%
            dplyr::left_join(df.2d, by = c("time", "ncol")) %>%
            dplyr::right_join(df.uniform %>%
                              dplyr::select(ncol, time, type),
                              by = c("time", "ncol")) %>%
            dplyr::left_join(df.1d, by = "ncol") %>%
            dplyr::left_join(df.hym, by = "lev") %>%
            dplyr::left_join(df.dhyi, by = "lev") %>%
            dplyr::mutate(p0 = p0) %>%
            dplyr::rename_with(tolower) %>%
            dplyr::mutate(time = as.POSIXct(sprintf("%d", nbdate), format = "%Y%m%d") + time * 86400)

        ret %<>% dplyr::filter(ncol %in% icol.select)

        str(ret)

        ## if (!missing(icol.select)) {
        ##     ret %>% dplyr::filter(ncol %in% icol.select)
        ## } else {
        ##     ret
        ## }

        ret
    }, .parallel = TRUE) %>%
    saveRDS("profiles_NEP_uniform_pblh.rds")
}

#' Calculate Sc occurrence from the EAMv1_entrain_NEP run
#'
#' @export
process.NEP.fscu <- function(path = "/global/cscratch1/sd/jmuelmen/e3sm_scratch/EAMv1_entrain_NEP/run/",
                             pattern = "EAMv1_entrain_NEP.cam.h1.2000-.*.nc") {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    doParallel::registerDoParallel(12)
    
    lf <- list.files(path, pattern,
                     full.names = TRUE)

    ncvar_get_NEP <- function(nc, varname, ...) {
        ncdf4::ncvar_get(nc, sprintf("%s_210e_to_245e_15n_to_35n", varname), ...)
    }
    
    plyr::ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        df.1d <- plotutils::nc.to.df(nc, names(nc$var)[grepl("lon|lat", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))
        
        ## 2D fields
        df.2d <- plotutils::nc.to.df(nc,
                                     names(nc$var)[grepl("OMEGA500|OMEGA700|TH7001000", names(nc$var))]) %>%
            dplyr::rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        df.2d %<>%
            dplyr::left_join(df.1d, by = "ncol")
        
    }, .parallel = TRUE) %>%
    dplyr::group_by(ncol, lon, lat) %>%
    dplyr::summarize(fscu = mean(OMEGA500 * 86400 > 10e2 & OMEGA700 * 86400 > 10e2 & TH7001000 > 18.55)) %>%
    dplyr::ungroup() %>%
    saveRDS("profiles_NEP_fscu.rds")
}

#' Pluck the grid point with max Sc occurrence from the EAMv1_entrain_NEP run
#'
#' @export
process.NEP <- function(path = "/global/cscratch1/sd/jmuelmen/e3sm_scratch/EAMv1_entrain_NEP/run/",
                        pattern = "EAMv1_entrain_NEP.cam.h1.2000-.*.nc", icol = 316) {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    doParallel::registerDoParallel(12)
    
    lf <- list.files(path, pattern,
                     full.names = TRUE)

    ncvar_get_NEP <- function(nc, varname, ...) {
        ncdf4::ncvar_get(nc, sprintf("%s_210e_to_245e_15n_to_35n", varname), ...)
    }
    
    ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        ## fields required to calculate 3D pressure field
        ps <- ncvar_get_NEP(nc, "PS", start = c(icol, 1), count = c(1,-1))
        p0 <- ncvar_get(nc, "P0")
        hyam <- ncvar_get(nc, "hyam")
        hybm <- ncvar_get(nc, "hybm")

        p.3d <- aperm(outer(ps, hybm) + outer(array(p0, dim(ps)), hyam),
                      c(2,1))

        ## other 2D fields
        pblh <- ncvar_get_NEP(nc, "PBLH", start = c(icol, 1), count = c(1,-1))
        precc <- ncvar_get_NEP(nc, "PRECC", start = c(icol, 1), count = c(1,-1))
        precl <- ncvar_get_NEP(nc, "PRECL", start = c(icol, 1), count = c(1,-1))
        qflx <- ncvar_get_NEP(nc, "QFLX", start = c(icol, 1), count = c(1,-1))
        shflx <- ncvar_get_NEP(nc, "SHFLX", start = c(icol, 1), count = c(1,-1))

        ## other 3D fields
        z3 <- ncvar_get_NEP(nc, "Z3", start = c(icol, 1, 1), count = c(1,-1,-1))
        t <- ncvar_get_NEP(nc, "T", start = c(icol, 1, 1), count = c(1,-1,-1))
        qv <- ncvar_get_NEP(nc, "Q", start = c(icol, 1, 1), count = c(1,-1,-1))
        ql <- ncvar_get_NEP(nc, "CLDLIQ", start = c(icol, 1, 1), count = c(1,-1,-1))
        omega <- ncvar_get_NEP(nc, "OMEGA", start = c(icol, 1, 1), count = c(1,-1,-1))
        prao <- ncvar_get_NEP(nc, "PRAO", start = c(icol, 1, 1), count = c(1,-1,-1))
        prco <- ncvar_get_NEP(nc, "PRCO", start = c(icol, 1, 1), count = c(1,-1,-1))
        qrl <- ncvar_get_NEP(nc, "QRL", start = c(icol, 1, 1), count = c(1,-1,-1))
        qrs <- ncvar_get_NEP(nc, "QRS", start = c(icol, 1, 1), count = c(1,-1,-1))
        cloud <- ncvar_get_NEP(nc, "CLOUD", start = c(icol, 1, 1), count = c(1,-1,-1))
        cloudfrac_clubb <- ncvar_get_NEP(nc, "CLOUDFRAC_CLUBB", start = c(icol, 1, 1), count = c(1,-1,-1))
        numliq <- ncvar_get_NEP(nc, "NUMLIQ", start = c(icol, 1, 1), count = c(1,-1,-1))
        ccn <- ncvar_get_NEP(nc, "CCN4", start = c(icol, 1, 1), count = c(1,-1,-1))

        ## 3D fields on full levels
        wp2_clubb <- ncvar_get_NEP(nc, "WP2_CLUBB", start = c(icol, 1, 1), count = c(1,-1,-1))
        wp3_clubb <- ncvar_get_NEP(nc, "WP3_CLUBB", start = c(icol, 1, 1), count = c(1,-1,-1))
        
        df.3d <- expand.grid(lev = 100 * ncvar_get(nc, "lev") %>% as.vector,
                             time = as.POSIXct("1999-10-01", tz = "UTC") + 86400 * ncvar_get(nc, "time") %>% as.vector) %>%
            mutate(z3   = as.vector(z3),
                   p    = as.vector(p.3d),
                   qv   = as.vector(qv),
                   ql   = as.vector(ql),
                   omega = as.vector(omega),
                   prao =  as.vector(prao),
                   prco =  as.vector(prco),
                   qrl =  as.vector(qrl),
                   qrs =  as.vector(qrs),
                   cloud = as.vector(cloud),
                   cloudfrac_clubb = as.vector(cloudfrac_clubb),
                   numliq = as.vector(numliq),
                   ccn = as.vector(ccn),
                   t    = as.vector(t)) %>%
            mutate(theta = t * (p0 / p) ^ (2/7),
                   theta_l = theta - 2.5e3 * ql,
                   qt = qv + ql) %>%
            ## total kluge, but hopefully level numbers will be static for a while
            group_by(time) %>%
            mutate(jk = 1:72) %>%
            ungroup()

        df.3d.ilev <- expand.grid(ilev = 100 * ncvar_get(nc, "ilev") %>% as.vector,
                             time = as.POSIXct("1999-10-01", tz = "UTC") + 86400 * ncvar_get(nc, "time") %>% as.vector) %>%
            mutate(wp2_clubb   = as.vector(wp2_clubb),
                   wp3_clubb   = as.vector(wp3_clubb)) %>%
            ## total kluge, but hopefully level numbers will be static for a while
            group_by(time) %>%
            mutate(jk = 0:72 + 0.5) %>%
            ungroup()
        
        df.2d <- expand.grid(time = as.POSIXct("1999-10-01", tz = "UTC") + 86400 * ncvar_get(nc, "time") %>% as.vector) %>%
            mutate(pblh         = as.vector(pblh),
                   qflx =  as.vector(qflx),
                   shflx =  as.vector(shflx),
                   precc        = as.vector(precc),
                   precl        = as.vector(precl))

        df <- bind_rows(df.3d, df.3d.ilev) %>%
            left_join(df.2d, by = "time")
    }, .parallel = TRUE) %>%
    saveRDS("profiles_NEP.rds")
}

#' Pluck the profiles that meet Sc critera from the EAMv1_entrain_CDNC-LWP run
#'
#' @export
process.Sc <- function(path = "/global/cscratch1/sd/jmuelmen/e3sm_scratch/EAMv1_entrain_CDNC-LWP/run",
                       pattern = "EAMv1_entrain_CDNC-LWP.cam.h1.2010-.*.nc") {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(plotutils)
    library(magrittr)
    doParallel::registerDoParallel(12)
    
    lf <- list.files(path, pattern, full.names = TRUE)
    
    ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        df.1d <- plotutils::nc.to.df(nc, c("lon", "lat"))
        
        ## 2D fields
        df.2d <- plotutils::nc.to.df(nc,
                                     c("PRECC",
                                       "PRECL",
                                       "cdnc",
                                       "lcc",
                                       "icc",
                                       "lwp",
                                       "iwp",
                                       "OMEGA500",
                                       "TH7001000",
                                       "LANDFRAC",
                                       "ttop"))
        
        ## 3D fields
        df.3d <- plotutils::nc.to.df(nc,
                                     c("T",
                                       "CLOUD",
                                       "CLDLIQ",
                                       "NUMLIQ"))

        df.3d %<>%
            group_by(ncol,time) %>%
            mutate(ilev = 1 : length(lev)) %>%
            ungroup()

        df.2d %<>%
            filter(icc < 1e-3 & ttop > 273.15 & OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55 & LANDFRAC < 1e-3) %>%
            left_join(df.1d, by = "ncol")

        df <- df.3d %>%
            filter(lev >= 690) %>%
            right_join(df.2d, by = c("time", "ncol"))

        df
    }, .parallel = TRUE) %>%
    saveRDS("profiles_Sc.rds")
}

#' Save PBL depth, Sc criteria, CCN
#'
#' @export
process.ccn <- function(path = "/global/cscratch1/sd/jmuelmen/e3sm_scratch/EAMv1_entrain_NEP/run",
                        pattern = "EAMv1_entrain.*cam.h1.2000-.*.nc") {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(plotutils)
    library(magrittr)
    doParallel::registerDoParallel(12)
    
    lf <- list.files(path, pattern, full.names = TRUE)
    
    ldply(lf, function(fname) {
        nc <- nc_open(fname)
        on.exit(nc_close(nc))

        gc()

        df.1d <- plotutils::nc.to.df(nc,
                                     sprintf("%s_210e_to_245e_15n_to_35n",
                                             c("lon", "lat"))) %>%
            rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))
        
        ## 2D fields
        df.2d <- plotutils::nc.to.df(nc,
                                     sprintf("%s_210e_to_245e_15n_to_35n",
                                             c("OMEGA500",
                                               "OMEGA700",
                                               "TH7001000",
                                               ## "LANDFRAC",
                                               "PS"))) %>%
            rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))
        
        ## 3D fields
        df.3d <- plotutils::nc.to.df(nc,
                                     sprintf("%s_210e_to_245e_15n_to_35n",
                                             c("T",
                                               "U",
                                               "V",
                                               "CCN4",
                                               "CLOUD",
                                               "CLDLIQ",
                                               "NUMLIQ"))) %>%
            rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .))

        df.3d %<>%
            group_by(ncol, time) %>%
            mutate(ilev = 1 : length(lev)) %>%
            summarize(jk.pbl = ilev[max(which(diff(T) < 0))],
                      jk.700 = ilev[which.min(abs(lev - 700))],
                      T.sfc = T[72], U.sfc = U[72], V.sfc = V[72], CCN.sfc = CCN4[72],
                      T.ft = T[jk.pbl - 1], U.ft = U[jk.pbl - 1], V.ft = V[jk.pbl - 1], 
                      T.700 = T[jk.700], U.700 = U[jk.700], V.700 = V[jk.700], 
                      CLOUD.top = CLOUD[jk.pbl + 1], NUMLIQ.top = NUMLIQ[jk.pbl + 1]) %>%
            ungroup()
            

        df.2d %<>%
            ## filter(icc < 1e-3 & ttop > 273.15 & OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55 & LANDFRAC < 1e-3) %>%
            left_join(df.1d, by = "ncol")

        df <- df.3d %>%
            right_join(df.2d, by = c("time", "ncol"))

        df
    }, .parallel = TRUE) %>%
    saveRDS("profiles_ccn.rds")
}

#' Worth a try, but doesn't fit in memory for the 3D variables for a year
postprocess.cdnc.lwp.cloud3d <- function(path = "~/scratch/e3sm_scratch/EAMv1_entrain_CDNC-LWP/postprocessed") {
    library(ncdf4)
    library(plotutils)
    library(plyr)
    library(dplyr)
    library(magrittr)

    nc <- nc_open(sprintf("%s/EAMv1_entrain_CDNC-LWP_cloud.nc", path))
    df.ttop <- nc.to.df(nc, c("ttop", "icc"))
    mask.ttop <- with(df.ttop, icc < 1e-3 & ttop > 273.15)
    rm(df.ttop)
    df.mask <- nc.to.df(nc, c("OMEGA500", "TH7001000", "LANDFRAC"))
    mask <- with(df.mask, mask.ttop & OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55 & LANDFRAC < 1e-3)
    rm(df.mask)
    df.global.tmp <- nc.to.df(nc, names(nc$var)[grepl("lcc|lwp|cdnc|ttop|icc|OMEGA[57]00|TH7001000|LANDFRAC", names(nc$var))], mask = mask)
    ## names(df) <- gsub("_210e_to_245e_15n_to_35n", "", names(df))
    nc_close(nc)
    nc.3d <- nc_open(sprintf("%s/EAMv1_entrain_CDNC-LWP_cloud_3d.nc", path))
    df.global.3d <- nc.to.df(nc.3d, names(nc.3d$var)[grepl("NUMLIQ", names(nc.3d$var))])
}

#' Add Medeiros and Stevens (2011)-like Sc flag to file
add.scu.flag <- function (fname = "../vignettes/EAMv1_entrain_NEP.cam.h1.2000_scu.nc") {
    library(ncdf4)
    library(plotutils)
    library(plyr)
    library(dplyr)
    library(magrittr)

    nc <- nc_open(fname)
    df <- nc.to.df(nc, names(nc$var)[-(1:3)])
    names(df) <- gsub("_210e_to_245e_15n_to_35n", "", names(df))
    df.lon.lat <- nc.to.df(nc, names(nc$var)[2:3])
    names(df.lon.lat) <- gsub("_210e_to_245e_15n_to_35n", "", names(df.lon.lat))
    df %<>% mutate(is.scu = OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55)
    nc_close(nc)

    df.mean <- ddply(df, ~ncol, colMeans) %>%
        full_join(df.lon.lat)

    df.monthly.mean <- ddply(df %>%
                             mutate(time = as.POSIXct("1999-10-01", tz = "UTC") + 86400 * time,
                                    month = factor(months(time),
                                                   levels = unique(months(time)),
                                                   ordered = TRUE)) %>%
                             select(-time),
                             ~ ncol + month, function(x) { colMeans(x %>% select(-month)) }) %>%
        full_join(df.lon.lat)

    df.jjas.mean <- ddply(df %>%
                          mutate(time = as.POSIXct("1999-10-01", tz = "UTC") + 86400 * time) %>%
                          filter(time >= as.POSIXct("2000-06-01", tz = "UTC"),
                                 time < as.POSIXct("2000-10-01", tz = "UTC")) %>%
                          select(-time),
                          ~ncol, colMeans) %>%
        full_join(df.lon.lat)

    ggplot(df.monthly.mean %>%
           mutate(lon = ifelse(lon > 180, lon - 360, lon)),
           aes(lon, lat, col = is.scu)) + geom_point(size = 5) +
        scale_color_distiller(palette = "Spectral") +
        geom_world_polygon(highres = TRUE) +
        scale_x_geo() +
        scale_y_geo() +
        facet_wrap(~ month) +
        coord_cartesian(xlim = c(-160, -110), ylim = c(10,40)) +
        theme_bw()

    df.max.scu <- df %>%
        mutate(time = as.POSIXct("1999-10-01", tz = "UTC") + 86400 * time) %>%
        filter(ncol == with(df.jjas.mean, ncol[which.max(is.scu)]))

    df.max.scu %>%
        mutate(OMEGA500 = OMEGA500 * 864,
               OMEGA700 = OMEGA700 * 864) %>%
        gather(var, val, OMEGA500, OMEGA700, TH7001000) %>%
        ggplot(aes(x = time, y = val)) +
        geom_line() +
        geom_hline(aes(yintercept = val), data = data.frame(var = c("OMEGA500", "OMEGA700", "TH7001000"), val = c(10,10,18.55)),
                   lty = "dashed", color = "darkgrey") +
        facet_grid(var ~ ., scales = "free_y") +
        theme_bw()

    append_scu_flag <- FALSE
    if (append_scu_flag) {
        nc <- nc_open(fname, write = TRUE)
        var.is_scu <- ncvar_def("is_scu_210e_to_245e_15n_to_35n", "", list(nc$dim$ncol_210e_to_245e_15n_to_35n, nc$dim$time), prec = "float")
        nc <- ncvar_add(nc, var.is_scu)
        ncvar_put(nc, var.is_scu, df$is.scu)
        nc_close(nc)
    }
}

plot_lwp_cdnc <- function() {
    library(ncdf4)
    library(plotutils)
    library(dplyr)
    library(plyr)
    library(ggplot2)
    nc <- nc_open("EAMv1_entrain_NEP_CDNC-LWP_cloud.nc")
    df <- nc.to.df(nc, names(nc$var)[!grepl("bnds|lat|lon", names(nc$var))])
    names(df) <- gsub("_210e_to_245e_15n_to_35n", "", names(df))
    df %<>% mutate(is.scu = OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55)

    df %>%
        filter(is.scu) %>%
        ggplot(aes(lcc)) + geom_histogram(binwidth = 0.001)

    df %>%
        filter(is.scu) %>%
        ggplot(aes(lwp)) + geom_histogram(binwidth = 0.001)

    df %>%
        ## filter(is.scu) %>%
        ggplot(aes(log10(iwp / icc))) + geom_histogram(binwidth = 0.001)

    df %>%
        filter(lcc > 0.1) %>%
        mutate(lwp_ic = lwp / lcc,
               cdnc_ic = cdnc / lcc) %>%
        filter(lwp_ic < 1) %>%
        discretize(lwp_ic, 100) %>%
        discretize(cdnc_ic, 100) %>%
        group_by(cdnc_ic, lwp_ic) %>%
        dplyr::summarize(p = n()) %>%
        ungroup() %>%
        group_by(cdnc_ic) %>%
        mutate(p = p / sum(p)) %>%
        ungroup() %>%
        ggplot(aes(cdnc_ic, lwp_ic)) +
        ## stat_bin2d(aes(fill = after_stat(count))) # , binwidth = c(3,1))
        geom_raster(aes(fill = p)) + geom_smooth()
    
    df %>%
        ## filter(is.scu) %>%
        filter(iwp < 1e-10) %>%
        filter(lcc > 0.1) %>%
        mutate(lwp_ic = lwp / lcc,
               cdnc_ic = cdnc / lcc) %>%
        filter(lwp_ic < 1) %>%
        discretize(lwp_ic, 100) %>%
        discretize(cdnc_ic, 100) %>%
        group_by(cdnc_ic, lwp_ic) %>%
        dplyr::summarize(p = n()) %>%
        ungroup() %>%
        group_by(cdnc_ic) %>%
        dplyr::summarize(lwp_ic = sum(p * lwp_ic) / sum(p)) %>%
        ungroup() %>%
        ggplot(aes(cdnc_ic, lwp_ic)) +
        geom_line()

    df %>%
        ## filter(is.scu) %>%
        filter(iwp < 1e-10) %>%
        filter(lcc > 0.1) %>%
        mutate(lwp_ic = lwp / lcc,
               cdnc_ic = cdnc / lcc) %>%
        filter(lwp_ic < 1) %>%
        discretize(lwp_ic, 100) %>%
        discretize(cdnc_ic, 100) %>%
        group_by(cdnc_ic) %>%
        dplyr::summarize(lwp_ic = mean(lwp_ic)) %>%
        ungroup() %>%
        ggplot(aes(cdnc_ic, lwp_ic)) +
        geom_line()
    
    df %>%
        ## filter(is.scu) %>%
        filter(iwp < 1e-10) %>%
        filter(lcc > 0.1) %>%
        mutate(lwp_ic = lwp,
               cdnc_ic = cdnc) %>%
        filter(lwp_ic < 1) %>%
        discretize(lwp_ic, 100) %>%
        discretize(cdnc_ic, 100) %>%
        group_by(cdnc_ic, lwp_ic) %>%
        dplyr::summarize(p = n()) %>%
        ungroup() %>%
        group_by(cdnc_ic) %>%
        dplyr::summarize(lwp_ic = sum(p * lwp_ic) / sum(p)) %>%
        ungroup() %>%
        ggplot(aes(cdnc_ic, lwp_ic)) +
        geom_line()

    df %>%
        ## filter(is.scu) %>%
        filter(iwp < 1e-10) %>%
        filter(lcc > 0.1) %>%
        discretize(cdnc, 100) %>%
        discretize(lcc, 100) %>%
        group_by(lcc) %>%
        dplyr::summarize(cdnc = mean(cdnc)) %>%
        ungroup() %>%
        ggplot(aes(lcc, cdnc)) +
        geom_line()
    
    df %>%
        filter(is.scu) %>%
        ## filter(lcc != 0) %>%
        ## ggplot(aes(lcc)) + geom_histogram(binwidth = 0.001)
        filter(lcc > 0.1) %>%
        ggplot(aes(cdnc / lcc, lwp / lcc)) +
        stat_bin2d(aes(fill = after_stat(count))) # , binwidth = c(3,1))
    
    str(df)

}

process.acme.ppe <- function(path) {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(plotutils)
    library(magrittr)
    library(mltools)
    doParallel::registerDoParallel(32)

    df.runs <- expand.grid(period = c("pi", "pd"),
                           model = sprintf("e%02d", 0:27))

    ddply(df.runs, ~ period + model, function(df.run) {
        gc()
        out.summary <- with(df.run, sprintf("scratch/acmev1_%s_%s_summary.rds", model, period))
        ret.summary <- try(readRDS(out.summary))
        if ("data.frame" %in% class(ret.summary))
            return(ret.summary)
        ## otherwise regenerate
        mccollect(mcparallel(with(df.run, {
            lf <- list.files(sprintf("%s/acmev1_%s_%s/ne30", path, model, period),
                             "e.._2007_.._.._cld_p..nc", full.names = TRUE)
            out.name <- sprintf("scratch/acmev1_%s_%s.rds", model, period)
            ## ret <- try(readRDS(out.name))
            ## if (!("data.frame" %in% class(ret))) {
            ##   print(sprintf("(re)generating %s", out.name))
                 ret <- ldply(lf, function(fname) {
                    nc <- nc_open(fname)
                    on.exit(nc_close(nc))
                    
                    gc()
                    
                    df.1d <- plotutils::nc.to.df(nc, c("lon", "lat"))
                    
                    ## 2D fields
                    df.2d <- plotutils::nc.to.df(nc,
                                                 c("cdnc_cldtop", "lwp", "cldfra", "PRECL", "FCTI")) %>%
                        dplyr::filter(FCTI == 0, is.finite(lwp))
                    
                    df.2d %>%
                        left_join(df.1d, by = "ncol")
                    
                }, .parallel = TRUE)
                mcparallel(saveRDS(ret, out.name), detached = TRUE)
            ## }
            ret %<>%
                dplyr::filter(cldfra > 0.9) %>%
                dplyr::mutate(lwp_ic = 1e-3 * lwp / cldfra, cdnc_ic = 1e6 * cdnc_cldtop / cldfra, model = model, period = toupper(period)) %>%
                dplyr::group_by(model, period) %>%
                lwp.cdnc.conditional()
            mcparallel(saveRDS(ret, out.summary), detached = TRUE)
            return(ret)
        })))
    }, .parallel = FALSE)
}

#' Process Po-Lun's EAGLES v2 PPElet (explicit list of runs is in the source)
#'
#' This function can (but doesn't have to) be run as a two-stage
#' workflow, first producing summary files for each E3SM run, then
#' combining them
#'
#' @return data.frame with conditional probability summaries for all
#'     runs; I saved these to eaglesv2_ppe.rds
#' 
#' @export
#' 
process.eagles.ppe <- function(path = "/global/cscratch1/sd/plma/e3sm_scratch") {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(plotutils)
    library(magrittr)
    library(mltools)
    doParallel::registerDoParallel(32)

    df.runs <- data.frame(model = c("v2r.rain.wup.auto.gccn.act.f2010",
                                    "v2r.rain.wup.auto.gccn.act.f2010.a1850",
                                    "v2r.rain.wup.auto.gccn.f2010.a1850",
                                    "v2r.rain.wup.auto.gccn.f2010",
                                    "v2r.rain.wup.auto.f2010.a1850",
                                    "v2r.rain.wup.auto.f2010",
                                    "v2r.rain.wup.f2010.a1850",
                                    "v2r.rain.wup.f2010",
                                    "v2r.rain.f2010.a1850",
                                    "v2r.rain2.f2010",
                                    "v2r.f2010.a1850",
                                    "v2r.f2010",
                                    "v2.f2010.a1850",
                                    "v2.f2010"))
    
    ddply(df.runs, ~ model, function(df.run) {
        gc()
        out.summary <- with(df.run, sprintf("scratch/%s_summary.rds", model))
        ret.summary <- try(readRDS(out.summary))
        if ("data.frame" %in% class(ret.summary))
            return(ret.summary)
        ## otherwise regenerate
        mccollect(mcparallel(with(df.run, {
            lf <- list.files(sprintf("%s/%s/run", path, model),
                             ".*eam.h2.2010-..-..-.*.nc", full.names = TRUE)
            out.name <- sprintf("scratch/%s.rds", model)
            ## ret <- try(readRDS(out.name))
            ## if (!("data.frame" %in% class(ret))) {
            ##   print(sprintf("(re)generating %s", out.name))
                 ret <- ldply(lf, function(fname) {
                    nc <- nc_open(fname)
                    on.exit(nc_close(nc))
                    
                    gc()
                    
                    df.1d <- plotutils::nc.to.df(nc, c("lon", "lat"))
                    
                    ## 2D fields
                    df.2d <- plotutils::nc.to.df(nc,
                                                 c("cdnc", "lwp", "lcc", "PRECL", "FCTI")) %>%
                        dplyr::filter(FCTI == 0, is.finite(lwp))
                    
                    df.2d %>%
                        left_join(df.1d, by = "ncol")
                    
                }, .parallel = TRUE)
                mcparallel(saveRDS(ret, out.name), detached = TRUE)
            ## }
            ret %<>%
                dplyr::filter(lcc > 0.9) %>%
                dplyr::mutate(lwp_ic = lwp / lcc,
                              cdnc_ic = cdnc / lcc,
                              period = ifelse(grepl("a1850", model), "PI", "PD"),
                              model = gsub(".a1850", "", model)) %>%
                dplyr::group_by(model, period) %>%
                lwp.cdnc.conditional()
            mcparallel(saveRDS(ret, out.summary), detached = TRUE)
            return(ret)
        })))
    }, .parallel = FALSE)
}

#' Process 3D profiles in h9; combine with time-step output of T, Q,
#' CLDLIQ in h8 for tendencies; combine with ACI diagnostics in h7
#'
#' @export
process.F2010A.entrain <- function(path = "~/scratch/e3sm_scratch/v2_F2010A0_00/run",
                                   pattern = ".*.eam.h9..*.nc",
                                   time.min = "1990-01-01",
                                   time.max = "2010-01-01",
                                   h.entrain = "eam.h9",
                                   h.tendenc = "eam.h8",
                                   h.aerocom = "eam.h7",
                                   parallel = TRUE) {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    library(mltools)
    doParallel::registerDoParallel(24)
    
    Sys.setenv(TZ = "UTC")

    ## h9, h8, h7 happened to be the tapes I used for these purposes
    ## in in v2, hence the numbers in the variable names
    lf.9 <- list.files(path, pattern, full.names = TRUE)
    lf.8 <- list.files(path, sub(h.entrain, h.tendenc, pattern), full.names = TRUE)
    lf.7 <- list.files(path, sub(h.entrain, h.aerocom, pattern), full.names = TRUE)

    stopifnot(length(lf.8) == length(lf.9))
    stopifnot(length(lf.7) == length(lf.9))

    ## print(lf)

    fname.out <- sprintf("entrain_F2010A_run-%s.rds",
                         ## run name is the second-to-last component of the path (the last is "/run/")
                         rev(strsplit(path, "/")[[1]])[2])
    fname.out.filtered <- sprintf("entrain_F2010A_run-%s_filtered.rds",
                         ## run name is the second-to-last component of the path (the last is "/run/")
                         rev(strsplit(path, "/")[[1]])[2])
    fname.out.fscu <- sprintf("entrain_F2010A_run-%s_fscu.rds",
                         ## run name is the second-to-last component of the path (the last is "/run/")
                         rev(strsplit(path, "/")[[1]])[2])

    print(sprintf("fname.out: %s", fname.out))

    ## print(lf)

    df <- plyr::ldply(1 : length(lf.9), function(i) {
        fname.out <- sub(h.entrain, "eam.entrain", lf.9[i]) %>%
            sub("nc$", "rds", .)
        print(sprintf("fname.out: %s", fname.out))
        ret <- try(readRDS(fname.out))
        if (!("try-error" %in% class(ret))) {
            if (!parallel)
                return(ret)
            else
                return(NULL)
        }
        
        nc.9 <- nc_open(lf.9[i])
        nc.8 <- nc_open(lf.8[i])
        nc.7 <- nc_open(lf.7[i])
        
        on.exit(nc_close(nc.9))
        on.exit(nc_close(nc.8), add = TRUE)
        on.exit(nc_close(nc.7), add = TRUE)

        gc()

        ## figure out how time steps align
        time.h9 <- ncvar_get(nc.9, "time")
        time.h8 <- ncvar_get(nc.8, "time")
        time.h7 <- ncvar_get(nc.7, "time")

        ## orient ourselves in time and space
        df.1d <- plotutils::nc.to.df(nc.9, c("lon", "lat"))
        nbdate <- ncdf4::ncvar_get(nc.9, "nbdate")
        p0 <- ncdf4::ncvar_get(nc.9, "P0")
        df.hym <- plotutils::nc.to.df(nc.9, c("hyam", "hybm"))
        df.dhyi <- with(plotutils::nc.to.df(nc.9, c("hyai", "hybi")),
                        ## layer thicknesses
                        data.frame(lev = df.hym$lev,
                                   dhyai = diff(hyai),
                                   dhybi = diff(hybi)))
        
        ret <- plyr::ldply(1 : length(time.h9), function(j) {
            gc()
            
            ## 2D fields
            df.2d <- plotutils::nc.to.df(nc.9,
                                         grep("^PRECC$|^PRECL$|^QFLX$|^SHFLX$|^OMEGA500$|^OMEGA700$|^TH7001000$",
                                              names(nc.9$var), value = TRUE),
                                         start = c(1, j), count = c(-1, 1))
            ## 3D fields
            df.3d <- plotutils::nc.to.df(nc.9,
                                         grep("^Z3$|^OMEGA$|^PRAO$|^PRCO$|^QRL$|^QRS$|^CLOUD$|^NUMLIQ$|^CCN4$|^DTCORE$|^DQCORE$|^DTCOND$",
                                              names(nc.9$var), value = TRUE),
                                         start = c(1, 1, j), count = c(-1, -1, 1))

            ## 3D fields from h8 that need forward time differencing; select time.h8 and following time step
            j.h8 <- which(time.h8 == time.h9[j])
            df.2d.h8 <- plotutils::nc.to.df(nc.8,
                                            grep("^PS$",
                                                 names(nc.8$var), value = TRUE),
                                            start = c(1, j.h8), count = c(-1, 2))
            df.3d.h8 <- plotutils::nc.to.df(nc.8,
                                            grep("^T$|^Q$|^CLDLIQ$",
                                                 names(nc.8$var), value = TRUE),
                                            start = c(1, 1, j.h8), count = c(-1, -1, 2))

            ## 2D ACI diagnostics from h7
            j.h7 <- which(time.h7 == time.h9[j])
            df.2d.h7 <- plotutils::nc.to.df(nc.7,
                                            grep("^cdnc$|^ccn$|^lwp$|^lcc$|^iwc$|^icc$|^ttop$",
                                                 names(nc.7$var), value = TRUE),
                                            start = c(1, j.h7), count = c(-1, 1))
            
            ## combine the various streams
            df <- df.3d.h8 %>%
                dplyr::left_join(df.2d.h8, by = c("time", "ncol")) %>%
                dplyr::left_join(df.3d, by = c("time", "ncol", "lev")) %>%
                dplyr::left_join(df.2d, by = c("time", "ncol")) %>%
                ## dplyr::left_join(df.entrain.3d, by = c("time", "ncol", "lev")) %>%
                ## dplyr::left_join(df.entrain.2d, by = c("time", "ncol")) %>%
                dplyr::left_join(df.1d, by = "ncol") %>%
                dplyr::left_join(df.hym, by = "lev") %>%
                dplyr::left_join(df.dhyi, by = "lev") %>%
                dplyr::mutate(p0 = p0) %>%
                dplyr::rename_with(tolower) 

            calculate_h_cheaply <- function(df) {
                df %>%
                    dplyr::mutate(
                               ## h on the vertical grid
                               h = z3[jk == jk.pbl],
                               ## subgrid h calculated from the theta_l profile
                               h_theta = NA, # reconstruct_inversion(p, z3, theta_l, debug.theta),
                               ## subgrid h calculated from the qt profile
                               h_q = NA, # reconstruct_inversion(p, z3, qt, debug.q),
                               ## same, but applied to pressure
                               p.pbl_theta = NA, # reconstruct_inversion_pressure(p, z3, theta_l, debug.theta),
                               p.pbl_q = NA) # reconstruct_inversion_pressure(p, z3, qt, debug.q))
            }

            calculate_aligned_forward_tendencies <- function(df) {
                ## *assuming* time is the slowest-varying dimension,
                ## calculate forward tendencies by subtracting the first
                ## half of the data.frame from the second
                block.diff <- function(var) {
                    len <- length(var)
                    mid <- len / 2
                    stopifnot(len %% 1 == 0)
                    c(var[(mid + 1) : len] - var[1 : mid], rep(NA, mid))
                }
                
                df %<>% dplyr::mutate(dt = block.diff(time))
                
                if ("difftime" %in% class(df$dt)) {
                    ## let difftime do the conversion to seconds
                    df %<>% dplyr::mutate(dt = as.double(dt, units = "secs"))
                } else {
                    ## assume [time] == days
                    df %<>% dplyr::mutate(dt = 86400 * dt)
                }
                
                df %>%
                    dplyr::mutate(dtheta_l.dt = block.diff(theta_l) / dt,
                                  dq_t.dt = block.diff(qt) / dt,
                                  dq_l.dt = block.diff(ql) / dt,
                                  dh.dt = block.diff(h) / dt,
                                  dh_theta.dt = NA,
                                  dh_q.dt = NA,
                                  dp.pbl_theta.dt = NA,
                                  dp.pbl_q.dt = NA)
            }

            ## calculate mixed-layer budgets
            ret <- df %>%
                rename_fields() %>%
                group_by(time, ncol) %>%
                mutate(jk = 1 : length(lev)) %>%
                ungroup() %>%
                calculate_p3d() %>%
                calculate_thermodynamics() %>%
                group_by(time, ncol) %>%
                find_inversion() %>%
                calculate_h_cheaply() %>%
                ungroup() %>%
                calculate_aligned_forward_tendencies() %>%
                slice(1 : (n() / 2)) %>%
                transform_dT_dtheta() %>%
                group_by(time, ncol) %>%
                calculate_budgets() %>%
                ungroup() %>%
                dplyr::select(-(cdnc:lcc)) %>%
                dplyr::left_join(df.2d.h7, by = c("time", "ncol")) %>%
                dplyr::mutate(time = as.POSIXct(sprintf("%d", nbdate), format = "%Y%m%d") + time * 86400)

            str(ret)

            ret
        })

        saveRDS(ret, fname.out)

        if (!parallel)
            ret
        else
            NULL
    }, .parallel = parallel) 

    ## saveRDS(df, fname.out)
    gc()

    df.fscu <- df %>%
        filter(time >= time.min, time < time.max) %>%
        mutate(is.scu = lts > 18.55 & omega500 * 864 > 10 & omega700 * 864 > 10) %>%
        mutate(month = months(time, abbreviate = TRUE),
               month = factor(month, unique(month))) %>%
        group_by(month, ncol) %>%
        summarize(fscu = mean(is.scu)) %>%
        ungroup() %>%
        group_by(ncol) %>%
        mutate(fscu.annual = mean(fscu)) %>%
        ungroup()

    df.fscu %>%
        saveRDS(fname.out.fscu)

    df.filtered <- df %>%
        filter_sc() %>%
        filter(time >= time.min, time < time.max) %>%
        mutate(month = months(time, abbreviate = TRUE),
               month = factor(month, unique(month))) %>%
        left_join(df.fscu)

    df.filtered %>%
        saveRDS(fname.out.filtered)

    
    ##     df.filtered
}


#' Select Sc regions from the output of process.F2010A.entrain()
#'
#' @export
regionalize.F2010A.entrain <- function(file = "~/scratch/", grid = "~/scratch/ne30pg2.rds") {
    library(dplyr)
    library(mltools)

    fname.regional <- gsub("filtered.rds$", "regional.rds", file)
    fname.golden <- gsub("filtered.rds$", "golden.rds", file)

    df.entrain <- readRDS(file) %>%
        left_join(readRDS(grid)) %>%
        filter(LANDFRAC == 0) %>%
        filter(fscu.annual > 0.3) %>%
        mutate(region = regionalize_seasonal(month, lon, lat, fscu.annual, LANDFRAC)) %>%
        filter(!is.na(region))
    saveRDS(df.entrain, fname.regional)

    df.entrain.golden <- df.entrain %>%
        filter(between(73 - jk.pbl, 6, 15)) %>%
        filter(abs(lat) < 40) %>%
        filter(icc == 0)
    saveRDS(df.entrain.golden, fname.golden)

}
                                       
#' Cobble on PRECC, which is only included in the process.F2010A.entrain() output as of mltools >= 0.0.14
#'
#' @export
cobble.on.precc <- function(path = "~/scratch/e3sm_scratch/v2_F2010A0_00/run",
                            pattern = ".*.eam.h9..*.nc",
                            time.min = "1990-01-01",
                            time.max = "2010-01-01",
                            h.entrain = "eam.h9",
                            parallel = TRUE) {

    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    library(mltools)
    doParallel::registerDoParallel(12)
    
    Sys.setenv(TZ = "UTC")

    ## h9, h8, h7 happened to be the tapes I used for these purposes
    ## in in v2, hence the numbers in the variable names
    lf.9 <- list.files(path, pattern, full.names = TRUE)

    fname.out <- sprintf("entrain_F2010A_run-%s_preccmask.rds",
                         ## run name is the second-to-last component of the path (the last is "/run/")
                         rev(strsplit(path, "/")[[1]])[2])

    print(sprintf("fname.out: %s", fname.out))

    ## print(lf)

    df <- plyr::ldply(1 : length(lf.9), function(i) {
        
        nc.9 <- nc_open(lf.9[i])
        ## nc.8 <- nc_open(lf.8[i])
        ## nc.7 <- nc_open(lf.7[i])
        
        on.exit(nc_close(nc.9))
        ## on.exit(nc_close(nc.8), add = TRUE)
        ## on.exit(nc_close(nc.7), add = TRUE)

        gc()

        ## figure out how time steps align
        nbdate <- ncdf4::ncvar_get(nc.9, "nbdate")
        df.2d <- plotutils::nc.to.df(nc.9, "PRECC") %>%
            dplyr::mutate(time = as.POSIXct(sprintf("%d", nbdate), format = "%Y%m%d") + time * 86400) %>%
            dplyr::filter(time >= time.min, time < time.max) %>%
            dplyr::filter(PRECC == 0) %>%
            dplyr::select(time, ncol)
    }, .parallel = parallel)

    saveRDS(df, fname.out)
}

#' Extract lon and lat for unstructured grids
#'
#' @export
e3sm.grid <- function(fname) {
    nc <- ncdf4::nc_open(fname)
    on.exit(ncdf4::nc_close(nc))
    dplyr::select(dplyr::left_join(plotutils::nc.to.df(nc, c("lon", "lat")),
                                   plotutils::nc.to.df(nc, "LANDFRAC")),
                  -time)
}

#'
#'
#' @export
merge.entrainment <- function(df, path, pattern = ".*.eam.hent.%04d-%02d.*") {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    library(mltools)
    doParallel::registerDoParallel(8)
    
    Sys.setenv(TZ = "UTC")

    timesteps_per_month <- function(month, year, dt) {
        month.p1 <- (month %% 12) + 1  ## because remainder is 0-based
        year.p1 <- year + month %/% 12
        difftime(as.POSIXct(sprintf("%d-%02d-01", year.p1, month.p1)),
                 as.POSIXct(sprintf("%d-%02d-01", year, month)),
                 units = "secs") / dt
    }
    
    summarize_entrainment <- function(df, trim = 0.05) {
        .time <- df$time %>% as.POSIXlt()
        df %>%
            mutate(year = .time$year + 1900,
                   month = .time$mon + 1) %>%
            filter(fscu.annual > 0.3) %>%
            filter(between(73 - jk.pbl, 6, 15)) %>%
            group_by(year, month, ncol) %>%
            summarize(n = sum(!is.na(E_theta)) * (1 - 2 * trim),
                      ntot = as.vector(timesteps_per_month(month[1], year[1], 27 * 3600)),
                      E_theta = mean(E_theta, trim, na.rm = TRUE),
                      E_q = mean(E_q, trim, na.rm = TRUE)) %>%
            ungroup() %>%
            mutate(sc_frac = n / ntot)
    }

    append_entrainment <- function(df, fname, template.var = "FSNT") {
        nc <- nc_open(fname, write = TRUE)
        nbdate <- ncdf4::ncvar_get(nc, "nbdate")
        df.to.match <- plotutils::nc.to.df(nc, template.var) %>%
            dplyr::mutate(time = as.POSIXct(sprintf("%d", nbdate), format = "%Y%m%d") + time * 86400)
        df %<>% left_join(df.to.match, ., by = c("ncol")) %>%
            mutate(sc_frac = replace(sc_frac, is.na(sc_frac), 0))
        fillval <- ncatt_get(nc, template.var)$missing_value
        var.E_theta <- ncvar_def("E_theta", "", list(nc$dim$ncol, nc$dim$time), prec = "float", missval = fillval)
        var.E_q <- ncvar_def("E_q", "", list(nc$dim$ncol, nc$dim$time), prec = "float", missval = fillval)
        var.sc_frac <- ncvar_def("sc_frac", "", list(nc$dim$ncol, nc$dim$time), prec = "float")
        nc %<>% ncvar_add(var.E_theta)
        nc %<>% ncvar_add(var.E_q)    
        nc %<>% ncvar_add(var.sc_frac)
        with(df, {
            ncvar_put(nc, var.E_theta, E_theta)
            ncvar_put(nc, var.E_q,     E_q)
            ncvar_put(nc, var.sc_frac, sc_frac)
        })
        nc_close(nc)
    }

    df %<>% summarize_entrainment()
    plyr::d_ply(df, ~ year + month, function(df) {
        fname <- list.files(path,
                            df %>% slice(1) %$% sprintf(pattern, year, month),
                            full.names = TRUE)
        append_entrainment(df, fname)
    })
}

#' Process a single run into Nd-LWP histograms
#'
#' @return none; data.frame with conditional probability summaries is
#'     saved to file
#' 
#' @export
#' 
process_e3sm_aer <- function(path = "/compyfs/dtn/muel306/EAMv1_nudged_2010/run",
                             pattern = "EAMv1_nudged_2010.cam.h2..*.nc",
                             time.min = "2010-01-01",
                             time.max = "2015-01-01",
                             version = 1) {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    library(mltools)
    library(plotutils)
    doParallel::registerDoParallel(8)
    
    Sys.setenv(TZ = "UTC")

    fname.out <- sprintf("e3smv%d_aer-%s.rds",
                         version,
                         ## run name is the second-to-last component of the path (the last is "/run/")
                         rev(strsplit(path, "/")[[1]])[2])

    lf <- list.files(path, pattern,
                          full.names = TRUE)
    df <- plyr::ldply(lf, mltools::load_lwp, .parallel = TRUE) %>%
        filter(time >= time.min, time < time.max)
    saveRDS(df, fname.out)
}

process_e3smv2_aer <- function(radiation = FALSE) {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    library(mltools)
    library(plotutils)
    doParallel::registerDoParallel(8)
    
    Sys.setenv(TZ = "UTC")

    lf.2010 <- list.files("~/scratch/e3sm_scratch/v2_aer_nudged/v2_2010aer_nudged/run", "v2_.*aer_nudged.eam.h7..*.nc", full.names = TRUE)
    df.2010 <- plyr::ldply(lf.2010, mltools::load_lwp, radiation = radiation, .parallel = TRUE) %>%
        filter(time >= "2010-01-01", time < "2015-01-01")
    saveRDS(df.2010, sprintf("v2_2010aer_nudged%s.rds", if (radiation) "_rad" else ""))
    
    lf.1850 <- list.files("~/scratch/e3sm_scratch/v2_aer_nudged/v2_1850aer_nudged/run", "v2_.*aer_nudged.eam.h7..*.nc", full.names = TRUE)
    df.1850 <- plyr::ldply(lf.1850, mltools::load_lwp, radiation = radiation, .parallel = TRUE) %>%
        filter(time >= "2010-01-01", time < "2015-01-01")
    saveRDS(df.1850, sprintf("v2_1850aer_nudged%s.rds", if (radiation) "_rad" else ""))

}
