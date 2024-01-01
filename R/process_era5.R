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

## path <- "/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag"
## pattern <- "amip-tdiag_.*.01_tdiag_t.nc" ## arbitrary *unique* basename
## lf <- ## "/tmp/CER-NEWS_CCCM_Aqua-FM3-MODIS-CAL-CS_RelB1_905906.20071226.hdf"
##     list.files(path = path, pattern = pattern, recursive = TRUE, full.names = TRUE)

#' @export
summary_2d_era5 <- function(df, df.2d) {
    df %>%
        ## make 2D summaries of the 3D data
        dplyr::group_by(lon, lat, time) %>%
        dplyr::summarize(
            ## add regime-identifying variables (omega.500, LTS)
            omega.500 = omega[which.min(abs(p - 500e2))],
            omega.700 = omega[which.min(abs(p - 700e2))],
            LTS = theta[which.min(abs(p - 700e2))] - theta[jk.surf],
            t.surf = t[jk.surf],
            ## identify boundary layer top
            jk.pbl = max(which(diff(t) < 0)),
            p.pbl = p[jk.pbl],
            ## boundary-layer averages
            xl.0 = mean(xl[lev > jk.pbl]),
            theta.0 = mean(theta[lev > jk.pbl]),
            theta_l.0 = mean(theta_l[lev > jk.pbl]),
            q.0 = mean(q[lev > jk.pbl]),
            q_t.0 = mean(q_t[lev > jk.pbl]),
            ## boundary-layer-average advection
            adv.theta_l.0 = mean(adv.theta_l[lev > jk.pbl]),
            adv.q_t.0 = mean(adv.q_t[lev > jk.pbl]),
            ## free-tropospheric values
            xl.plus = xl[lev == jk.pbl],
            theta.plus = theta[lev == jk.pbl],
            theta_l.plus = theta_l[lev == jk.pbl],
            q.plus = q[lev == jk.pbl],
            q_t.plus = q_t[lev == jk.pbl],
            ## boundary-layer cloudiness (assuming max overlap)
            f.cloud = max(aclcac[lev > jk.pbl]),
            ## radiative cooling
            rad.cool =
                (flx_uplw[jk.pbl + 0] - flx_uplw[jk.surf]) -
                (flx_dnlw[jk.pbl + 0] - flx_dnlw[jk.surf]) +
                (flx_upsw[jk.pbl + 0] - flx_upsw[jk.surf]) -
                (flx_dnsw[jk.pbl + 0] - flx_dnsw[jk.surf])) %>%
        ## add other 2D fields
        dplyr::ungroup() %>%
        dplyr::left_join(df.2d, by = c("lon", "lat", "time")) %>%
        dplyr::ungroup() %>%
        ## pre-filter (using variables that do not need full time series)
        dplyr::filter(## LTS > 18.55, 
                      ## omega.500 * 86400 > 10e2, ## 10 hPa / day
                      abs(lat) < 45) ##,
                      ## p.pbl > 850e2,
                      ## jk.pbl < 46,
                      ## t.surf > 283.15)
}

#' @export
summary_3d.inversion_era5 <- function(df, df.2d, jk.surf) {
    df.scu <- df %>%
        dplyr::left_join(df.2d, by = c("lon", "lat", "time")) %>%
        dplyr::group_by(lon, lat, time) %>%
        dplyr::mutate(
                   LTS = theta[which.min(abs(p - 700e2))] - theta[jk.surf]) %>%
        ## select SCu; conditions from Medeiros and Stevens (2009)
        dplyr::filter(LTS > 18.55, ## 18.55 K
                      omega.500 * 86400 > 10e2) %>% ## 10 hPa / day
        dplyr::ungroup()
    
    ## select and plot theta and q in the vicinity of the PBL top
    df.inversion <- df.scu %>%
        dplyr::group_by(lon, lat, time) %>%
        dplyr::mutate(zinv = max(which(diff(t) < 0))) %>%
        dplyr::mutate(lev = lev - 99) %>% ## lev starts at 100; shift so it can be used as index
        dplyr::mutate(delta.zpbl = zpbl - zinv) %>%
        tidyr::gather(ytype, jk.pbl, ##zpbl,
                      zinv) %>%
        dplyr::mutate(delta.k = lev - jk.pbl,
                      delta.p = p - p[jk.pbl],
                      p.pbl = p[jk.pbl],
                      t.surf = t[jk.surf]) %>%
        dplyr::filter(abs(lat) < 25,
                      p.pbl > 850e2,
                      t.surf > 283.15) %>%
        dplyr::filter(abs(delta.k) <= 10) %>%
        dplyr::ungroup()
    
    df.inversion.relative <- df.inversion %>%
        tidyr::gather(xtype, x, q, theta, xl, q_t, theta_l, aclcac, tt_lw, tt_sw) %>%
        dplyr::group_by(lon, lat, time, xtype, ytype) %>%
        dplyr::mutate(x.pbl = mean(x[delta.k > 0]),
                      delta.x = x - ifelse(xtype %in% c("q", "theta", "q_t", "theta_l"),
                                           x.pbl, 0),
                      max.delta.k = max(delta.k),
                      k.norm = delta.k / max(delta.k),
                      p.norm = delta.p / (aps - p.pbl)) %>%
        dplyr::ungroup()
}

#' @export
process_era5 <- function(in.name,
                         out.name,
                         ncores,
                         sum.fun = summary_2d,
                         max.step = -1,
                         stride = 8,
                         filter = FALSE) {

    ## filename stubs for all tdiag files 
    ## prefix <- gsub("([0-9]*)_[a-z_]*\\.nc$", "\\1_tdiag_.*.nc", in.name)

    ## do not descend into subdirectories (which could be, e.g., daily means)
    lf2 <- list.files(path = in.name, pattern = ".*_ml.nc", 
                      recursive = FALSE, full.names = TRUE)


    ## cl <- snow::makeCluster(rep("localhost", ncores), type = "SOCK", outfile = "/dev/null")
    ## on.exit({
    ##     snow::stopCluster(cl)
    ## })
    ## doSNOW::registerDoSNOW(cl)
    
    ## ## snow::clusterEvalQ(cl, ls())
    
    ## snow::clusterExport(cl, "lf2", environment())
    ## snow::clusterEvalQ(cl, {
    ##     library(plyr)
    ##     library(dplyr)
        
    ##     nc_getstep <- function(nc, varname, stride, i, fname = varname) {
    ##         ncdf4::ncvar_get(nc, varname)
    ##     }
        
    ##     nc <- plyr::llply(lf2, ncdf4::nc_open)
    ##     names(nc) <- gsub("[^_]*_[0-9.]*_tdiag_([a-z0-9_]*)\\.nc$", "\\1", basename(lf2))
        
    ##     lon  <- ncdf4::ncvar_get(nc[["xl"]], "lon")
    ##     lat  <- ncdf4::ncvar_get(nc[["xl"]], "lat")
    ##     lev  <- ncdf4::ncvar_get(nc[["xl"]], "lev")
    ##     hyam <- ncdf4::ncvar_get(nc[["xl"]], "hyam")
    ##     hybm <- ncdf4::ncvar_get(nc[["xl"]], "hybm")
    ##     jk.surf <- length(hyam)
    ##     dx <- 6.371e6 * pi * outer(rep(median(stats::filter(lon, 1:-1), na.rm = TRUE),
    ##                                    length(lon)),
    ##                                cos(pi * lat / 180)) / 180
    ##     dy <- pi * array(6.371e6 * stats::filter(lat, 1:-1)) / 180
    ## })
    
    fname <- lf2[1]
    nc.ml <- nc_open(fname)
    nc.sf <- nc_open(gsub("_ml\\.nc", "_sf.nc", fname))
    nc.pl <- nc_open(gsub("_ml\\.nc", "_pl.nc", fname))

    lon  <- ncdf4::ncvar_get(nc.ml, "longitude")
    lat  <- ncdf4::ncvar_get(nc.ml, "latitude")
    lev  <- ncdf4::ncvar_get(nc.ml, "level")
    pl   <- ncvar_get(nc.pl, "level")
    time <- nc_getstep(nc.ml, "time", stride, i, "t")
    jk.surf <- length(lev)
    dx <- 6.371e6 * pi * outer(rep(median(stats::filter(lon, 1:-1), na.rm = TRUE),
                                   length(lon)),
                               cos(pi * lat / 180)) / 180
    dy <- pi * array(6.371e6 * stats::filter(lat, 1:-1)) / 180


    ## tm1 <- ncvar_get(nc, "tm1", c(1,1,1,(0 * 192) + 1), c(-1,-1,-1,192))
    ## qm1 <- ncvar_get(nc, "qm1", c(1,1,1,(0 * 192) + 1), c(-1,-1,-1,192))
    ## omega <- ncvar_get(nc, "omega", c(1,1,1,(0 * 192) + 1), c(-1,-1,-1,192))
    ## lh <- ncvar_get(nc, "ahflwac", NA, c(-1,-1,192))
    ## sh <- ncvar_get(nc, "ahfswac", NA, c(-1,-1,192))
    ## zpbl <- ncvar_get(nc, "zpbl", NA, c(-1,-1,192))
    ## aps <- ncvar_get(nc, "aps", NA, c(-1,-1,192))
    ahflwac <- nc_getstep(nc.sf, "slhf", stride, i) / 3600 ## accumulated over the last hour (J m^-2)
    ahfswac <- nc_getstep(nc.sf, "sshf", stride, i) / 3600 ## accumulated over the last hour (J m^-2)
    zpbl <- nc_getstep(nc.sf, "blh", stride, i) ## boundary layer height (m)
    aps <- nc_getstep(nc.sf, "sp", stride, i) ## surface pressure (Pa)
    aprl <- nc_getstep(nc.sf, "lsp", stride, i) / 3600 ## large-scale precip (m)
    aprc <- nc_getstep(nc.sf, "cp", stride, i) / 3600 ## convective precip (m)

    ## 2D variables
    df.2d <- expand.grid(lon = as.vector(lon),
                         lat = as.vector(lat),
                         time = as.vector(time)) %>%
        dplyr::mutate(ahflwac = as.vector(ahflwac),
                      ahfswac = as.vector(ahfswac),
                      zpbl = as.vector(zpbl),
                      aps = as.vector(aps),
                      aprl = as.vector(aprl),
                      aprc = as.vector(aprc))

    ## pressure-level variables (only omega currently); convert to 2d
    df.pl <- expand.grid(lon = as.vector(lon),
                         lat = as.vector(lat),
                         pl = as.vector(pl),
                         time = as.vector(time)) %>%
        dplyr::mutate(omega = as.vector(nc_getstep(nc.pl, "w", stride, i))) %>%
        tidyr::spread(pl, omega, sep = ".") %>%
        rename(omega.500 = pl.500, omega.700 = pl.700)

    df.2d %<>%
        dplyr::left_join(df.pl)
    
    data(sigma_l137)
    sigma_l137 %<>%
        dplyr::right_join(data.frame(lev = as.vector(lev)), by = "lev")
    hyam <- sigma_l137$a
    hybm <- sigma_l137$b

    p.3d <- aperm(outer(aps, hybm) + outer(array(1, dim(aps)), hyam),
                  c(1,2,4,3))
    
    ## u and v for advection
    u <- nc_getstep(nc.ml, "u", stride, i, "uv")
    v <- nc_getstep(nc.ml, "v", stride, i, "uv")

    ## scalars for which advection is to be calculated
    q <- nc_getstep(nc.ml, "q", stride, i)
    xl <- nc_getstep(nc.ml, "clwc", stride, i)
    t <- nc_getstep(nc.ml, "t", stride, i, "t")
        
    res <- expand.grid(lon = as.vector(lon),
                       lat = as.vector(lat),
                       lev = as.vector(lev),
                       time = as.vector(time)) %>%
        ## read in 3D fields
        dplyr::mutate(## adv.q  = as.vector(mltools::advection(q , u, v, dx, dy)),
                   ## adv.xl = as.vector(mltools::advection(xl, u, v, dx, dy)),
                   ## adv.t  = as.vector(mltools::advection(t , u, v, dx, dy)),
                   q  = as.vector(q ),
                   xl = as.vector(xl),
                   t  = as.vector(t ),
                   aclcac = as.vector(nc_getstep(nc, "cc", stride, i)),
                   tt_sw = as.vector(nc_getstep(nc.ml, "mttswr", stride, i)),
                   tt_lw = as.vector(nc_getstep(nc.ml, "mttlwr", stride, i)),
                   p = as.vector(p.3d)) %>%
        ## calculate thermodynamic variables
        dplyr::mutate(theta = t * (1e5 / p) ^ (2/7),
                          theta_l = theta - 2.5e3 * xl,
                          q_t = q + xl,
                          adv.theta = adv.t * (1e5 / p) ^ (2/7),
                          adv.theta_l = adv.theta - 2.5e3 * adv.xl,
                          adv.q_t = adv.q + adv.xl) %>%
            sum.fun(df.2d)

    gc()

    res
    ## }, .progress = "text", .parallel = TRUE, .inform = TRUE, stride = stride)
    
    
    res %<>%
        ## calculate tendencies
        dplyr::group_by(lon, lat) %>%
        dplyr::mutate(dq_t.0.dt = (q_t.0 - dplyr::lag(q_t.0)) / 3600,
                      dtheta_l.0.dt = (theta_l.0 - dplyr::lag(theta_l.0)) / 3600,
                      lh = ahflwac, ## in ECHAM, ahf{l,s}wac are accumulated since the beginning of the run
                      sh = ahfswac) %>%
        dplyr::ungroup() ##  %>%
        ## dplyr::filter(LTS > 18.55, 
        ##               omega.500 * 86400 > 10e2, ## 10 hPa / day
        ##               abs(lat) < 45,
        ##               p.pbl > 850e2,
        ##               jk.pbl < 46,
        ##               t.surf > 283.15)

    if (filter) {
        ## low-pass filter tendencies, LH, SH, precip (1 h; can afford to lose 1 h per month)
        filt <- rep(1/8, 8)
        res %<>%
            dplyr::group_by(lon, lat) %>%
            dplyr::mutate(dq_t.0.dt     = as.vector(stats::filter(dq_t.0.dt, filt)),
                          dtheta_l.0.dt = as.vector(stats::filter(dtheta_l.0.dt, filt)),
                          lh   = as.vector(stats::filter(lh, filt)),
                          sh   = as.vector(stats::filter(sh, filt)),
                          aprl = as.vector(stats::filter(aprl, filt)),
                          aprc = as.vector(stats::filter(aprc, filt))) %>%
            dplyr::ungroup()
        ## replace omega with filtered omega (168 h, so needs to be processed a year at a time)
        lf.omega <- list.files(path = sprintf("%s/filt168h", dirname(prefix)),
                               pattern = basename(prefix), 
                               recursive = FALSE, full.names = TRUE)
        nc.omega <- ncdf4::nc_open(lf.omega)
        df.omega <- expand.grid(lon =  as.vector(ncdf4::ncvar_get(nc.omega, "lon")),
                                lat =  as.vector(ncdf4::ncvar_get(nc.omega, "lat")),
                                lev =  as.vector(ncdf4::ncvar_get(nc.omega, "lev")),
                                time = as.vector(ncdf4::ncvar_get(nc.omega, "time"))) %>%
            dplyr::mutate(omega = as.vector(ncdf4::ncvar_get(nc.omega, "omega"))) %>%
            dplyr::mutate(lev = lev * 1e-2) %>%
            tidyr::spread(lev, omega, sep = ".")

        df.omega %<>%
            dplyr::rename(omega.700.filt = lev.700,
                          omega.500.filt = lev.500) %>%
            dplyr::mutate(time = round(time * 24 * 60 / 7.5 - 0.5) / (24 * 60 / 7.5))
        
        res %<>%
            dplyr::mutate(time = round(time * 24 * 60 / 7.5) / (24 * 60 / 7.5)) %>%
            dplyr::left_join(df.omega, by = c("lon", "lat", "time"))

        res
    }
    
    saveRDS(res, out.name)

    res
}

## postprocess used to perform work (calculate tendencies, cut on Scu
## conditions, rename variables), but now it just combines months into
## an annual file
postprocess <- function() {
    library(dplyr)
    lf <- list.files("/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag",
                     ".*tdiag.rds", full.names = TRUE)
    doParallel::registerDoParallel(12)
    plyr::ldply(lf, function(fname) {
        readRDS(fname) %>%
            dplyr::filter(LTS > 18.55, 
                          omega.500 * 86400 > 10e2, ## 10 hPa / day
                          abs(lat) < 45,
                          p.pbl > 850e2,
                          jk.pbl < 46,
                          t.surf > 283.15) 
    }, .parallel = TRUE) %>%
        saveRDS("ml-1979.rds")
}

## filtered-field postprocessing
postprocess.filt <- function() {
    library(dplyr)
    lf <- list.files("/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag",
                     ".*tdiag_filtered.rds", full.names = TRUE)
    doParallel::registerDoParallel(12)
    plyr::ldply(lf, function(fname) {
        readRDS(fname) %>%
            dplyr::filter(LTS > 18.55, 
                          omega.500 * 86400 > 10e2, ## 10 hPa / day
                          abs(lat) < 45,
                          p.pbl > 850e2,
                          jk.pbl < 46,
                          t.surf > 283.15) 
    }, .parallel = TRUE) %>%
        saveRDS("ml-filt-1979.rds")
}

## daily mean version of postprocess; no other data reduction is required
postprocess.daymean <- function() {
    library(magrittr)
    lf <- list.files("/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag/daymean",
                     ".*tdiag_daymean.rds", full.names = TRUE)
    doParallel::registerDoParallel(12)
    plyr::ldply(lf, readRDS, .parallel = TRUE) %>%
        dplyr::filter(LTS > 18.55, 
                      omega.500 * 86400 > 10e2, ## 10 hPa / day
                      abs(lat) < 45,
                      p.pbl > 850e2,
                      jk.pbl < 46,
                      t.surf > 283.15) %>%
            saveRDS("ml-daymean-1979.rds")
}

## omega is special, since we may want to band-pass filter it and
## therefore need the entire time series at every grid point
postprocess.omega <- function() {
    library(dplyr)
    lf <- list.files("/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag",
                     ".*tdiag.rds", full.names = TRUE)
    doParallel::registerDoParallel(12)
    plyr::ldply(lf, function(fname) {
        readRDS(fname) %>%
            select(lon, lat, time, omega.500, omega.700)
    }, .parallel = FALSE) %>%
        saveRDS("ml-omega-1979.rds")
}

filter.omega <- function() {
    library(magrittr)
    doParallel::registerDoParallel(12)
    df.omega <- readRDS("ml-omega-1979.rds")
    df.omega %<>%
        dplyr::arrange(lon, lat, time)
    gc()
    df.omega %<>%
        dplyr::mutate(omega.500 = as.vector(stats::filter(omega.500, rep(1 / (8 * 168), 8 * 168))))
    gc()
    df.omega %<>%
        dplyr::mutate(omega.700 = as.vector(stats::filter(omega.700, rep(1 / (8 * 168), 8 * 168))))
    gc()
    saveRDS(df.omega, "ml-omega-filt-1979.rds")
        plyr::ddply(~ lon + lat, dplyr::mutate,
                    omega.500.filt7day = stats::filter(omega.500, rep(1 / (8 * 168), 8 * 168)),
                    omega.700.filt7day = stats::filter(omega.700, rep(1 / (8 * 168), 8 * 168)), .progress = "text")

    dplyr::left_join(readRDS("ml-1979.rds") %>%
                     dplyr::select(-c(omega.500, omega.700)),
                     readRDS("ml-omega-filt-1979.rds")) %>%
        saveRDS("ml-filt-1979.rds")
}

## select the full time series of omega for analysis of noisiness and
## filtering
select.omega <- function() {
    df <- readRDS("~/ml-omega-1979.rds")
    lat.max <- unique(df$lat)[15]
    df.omega <- filter(df, lon == 360 - 82.5, lat == lat.max)
    saveRDS(df.omega, "~/ml-omega-fscumax-1979.rds")
}
