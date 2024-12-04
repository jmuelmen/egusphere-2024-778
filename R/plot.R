#' @export
plot.sample.era5 <- function() {
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(plotutils)
    library(ncdf4)

    nc.ml <- nc_open("/work/bb1036/b380126/era5/era5_20171114_0600_ml.nc")
    nc.sf <- nc_open("/work/bb1036/b380126/era5/era5_20171114_0600_sf.nc")
    nc.pl <- nc_open("/work/bb1036/b380126/era5/era5_20171114_0600_pl.nc")
    lon <- ncvar_get(nc.ml, "longitude")
    lat <- ncvar_get(nc.ml, "latitude")
    lev <- ncvar_get(nc.ml, "level")
    pl <- ncvar_get(nc.pl, "level")

    nc_getstep <- function(nc, varname, stride, i, fname = varname) {
        ## don't need to bother with all the stuff we need in the
        ## enormous ECHAM files
        ncdf4::ncvar_get(nc, varname)
    }

    ## 2D fields
    ahflwac <- nc_getstep(nc.sf, "slhf", stride, i) / 3600 ## accumulated over the last hour (J m^-2)
    ahfswac <- nc_getstep(nc.sf, "sshf", stride, i) / 3600 ## accumulated over the last hour (J m^-2)
    zpbl <- nc_getstep(nc.sf, "blh", stride, i) ## boundary layer height (m)
    aps <- nc_getstep(nc.sf, "sp", stride, i) ## surface pressure (Pa)
    aprl <- nc_getstep(nc.sf, "lsp", stride, i) / 3600 ## large-scale precip (m)
    aprc <- nc_getstep(nc.sf, "cp", stride, i) / 3600 ## convective precip (m)

    ## ## save TOA and surf fluxes separately
    ## flx_uplw_toa = as.vector(nc_getstep(nc, "flx_uplw", stride, i)[,,1,]),
    ## flx_dnlw_toa = as.vector(nc_getstep(nc, "flx_dnlw", stride, i)[,,1,]),
    ## flx_upsw_toa = as.vector(nc_getstep(nc, "flx_upsw", stride, i)[,,1,]),
    ## flx_dnsw_toa = as.vector(nc_getstep(nc, "flx_dnsw", stride, i)[,,1,]),
    ## ## save TOA and surf fluxes separately
    ## flx_uplw_surf = as.vector(nc_getstep(nc, "flx_uplw", stride, i)[,,48,]),
    ## flx_dnlw_surf = as.vector(nc_getstep(nc, "flx_dnlw", stride, i)[,,48,]),
    ## flx_upsw_surf = as.vector(nc_getstep(nc, "flx_upsw", stride, i)[,,48,]),
    ## flx_dnsw_surf = as.vector(nc_getstep(nc, "flx_dnsw", stride, i)[,,48,]))

    ## u and v for advection
    u <- nc_getstep(nc.ml, "u", stride, i, "uv")
    v <- nc_getstep(nc.ml, "v", stride, i, "uv")

    ## scalars for which advection is to be calculated
    q <- nc_getstep(nc.ml, "q", stride, i)
    xl <- nc_getstep(nc.ml, "clwc", stride, i)
    t <- nc_getstep(nc.ml, "t", stride, i, "t")
    
    omega = as.vector(nc_getstep(nc.pl, "w", stride, i)) ## omega (Pa s^-1)
    aclcac = as.vector(nc_getstep(nc, "cc", stride, i))
    ## ## for fluxes (on half levels), discard TOA,
    ## ## treat the rest as if on full level above the half level
    ## flx_uplw = as.vector(nc_getstep(nc, "flx_uplw", stride, i)[,,-1,]),
    ## flx_dnlw = as.vector(nc_getstep(nc, "flx_dnlw", stride, i)[,,-1,]),
    ## flx_upsw = as.vector(nc_getstep(nc, "flx_upsw", stride, i)[,,-1,]),
    ## flx_dnsw = as.vector(nc_getstep(nc, "flx_dnsw", stride, i)[,,-1,]),
    ## ## get radiative heating directly
    tt_sw = as.vector(nc_getstep(nc.ml, "mttswr", stride, i))
    tt_lw = as.vector(nc_getstep(nc.ml, "mttlwr", stride, i))

}

#' @export
plot.sample <- function() {
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(plotutils)
    library(ncdf4)


    nc3d_getstep <- function(nc, varname, i) {
        ncvar_get(nc, varname,
                  c(1,1,1,(i - 1) * 192 + 1),
                  c(-1,-1,-1,192))
    }

    nc <- nc_open("/scratch/b/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag/amip-tdiag_197901.01_tdiag.nc")
    lon <- ncvar_get(nc, "lon")
    lat <- ncvar_get(nc, "lat")
    lev <- ncvar_get(nc, "lev")
    hyam <- ncvar_get(nc, "hyam")
    hybm <- ncvar_get(nc, "hybm")
    time <- ncvar_get(nc, "time", NA, 192)

    ## tm1 <- ncvar_get(nc, "tm1", c(1,1,1,(0 * 192) + 1), c(-1,-1,-1,192))
    ## qm1 <- ncvar_get(nc, "qm1", c(1,1,1,(0 * 192) + 1), c(-1,-1,-1,192))
    ## omega <- ncvar_get(nc, "omega", c(1,1,1,(0 * 192) + 1), c(-1,-1,-1,192))
    lh <- ncvar_get(nc, "ahflwac", NA, c(-1,-1,192))
    sh <- ncvar_get(nc, "ahfswac", NA, c(-1,-1,192))
    zpbl <- ncvar_get(nc, "zpbl", NA, c(-1,-1,192))
    aps <- ncvar_get(nc, "aps", NA, c(-1,-1,192))

    nc.T <- nc_open("/work/bb0839/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag/amip-tdiag_197901.01_tdiag_t.nc")

    jk.surf <- length(hyam)
    
    df.2d <- expand.grid(lon = as.vector(lon),
                         lat = as.vector(lat),
                         time = as.vector(time)) %>%
        mutate(lh = as.vector(lh),
               sh = as.vector(sh),
               zpbl = as.vector(zpbl),
               aps = as.vector(aps))

    p.3d <- aperm(outer(aps, hybm) + outer(array(1, dim(aps)), hyam),
                  c(1,2,4,3))

    df.3d <- expand.grid(lon = as.vector(lon),
                         lat = as.vector(lat),
                         lev = as.vector(lev),
                         time = as.vector(time)) %>%
        mutate(q = as.vector(nc3d_getstep(nc, "q", 1)),
               xl = as.vector(nc3d_getstep(nc, "xl", 1)),
               t = as.vector(nc3d_getstep(nc.T, "st", 1)),
               omega = as.vector(nc3d_getstep(nc, "omega", 1)),
               p = as.vector(p.3d))  %>%
        mutate(theta = t * (1e5 / p) ^ (2/7),
               theta.l = theta - 2.5e3 * xl,
               qt = q + xl) %>%
        left_join(df.2d, by = c("lon", "lat", "time")) %>%
        ## add regime-identifying variables (omega.500, LTS)
        group_by(lon, lat, time) %>%
        mutate(omega.500 = omega[which.min(abs(p - 500e2))],
               omega.700 = omega[which.min(abs(p - 700e2))],
               LTS = theta[which.min(abs(p - 700e2))] - theta[jk.surf]) %>%
        ungroup()

    gc()

    df.3d %>%
        filter(lev > 35) %>%
        group_by(lat, lev) %>%
        summarize(theta = mean(theta)) %>%
        ggplot(aes(x = lat, y = lev, fill = theta)) +
        geom_raster() +
        scale_fill_distiller() +
        scale_y_reverse() +
        theme_bw(24)

    o500 <- df.3d %>% group_by(lon, lat) %>% summarize(omega.500 = mean(omega.500))

    mutate(df.3d,
           scu = LTS > 18.55 & ## 18.55 K
               omega.500 * 86400 > 10e2) %>% ## 10 hPa / day
        group_by(lon, lat) %>%
        summarize(fscu = mean(scu)) %>%
        ungroup() %>%
        mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
        ggplot(aes(x = lon, y = lat, fill = fscu)) +
        geom_raster() +
        scale_fill_distiller() +
        scale_x_geo() + scale_y_geo() + geom_world_polygon() +
        theme_bw(24)

    df.3d %>%
        filter(time == time[1]) %>%
        group_by(lon, lat, time) %>%
        summarize(theta.zpbl = theta[zpbl[1]],
                  theta.surf = theta[47],
                  t.surf = t[47]) %>%
        ungroup() %>%
        mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
        ggplot(aes(x = lon, y = lat, fill = t.surf > 273.15)) +
        geom_raster() +
        scale_fill_brewer(palette = "Spectral") +
        scale_x_geo() + scale_y_geo() + geom_world_polygon() +
        theme_bw(24)


    ## select SCu; conditions from Medeiros and Stevens (2009)
    df.scu <- filter(df.3d,
                     LTS > 18.55, ## 18.55 K
                     omega.500 * 86400 > 10e2) ## 10 hPa / day

    ## select and plot theta and q in the vicinity of the PBL top
    df.inversion <- df.scu %>%
        group_by(lon, lat, time) %>%
        mutate(zinv = max(which(diff(t) < 0))) %>%
        mutate(delta.zpbl = zpbl - zinv) %>%
        tidyr::gather(ytype, jk.pbl, zpbl, zinv) %>%
        mutate(delta.k = lev - jk.pbl,
               delta.p = p - p[jk.pbl],
               p.pbl = p[jk.pbl],
               t.surf = t[jk.surf]) %>%
        filter(abs(lat) < 25,
               p.pbl > 850e2,
               t.surf > 283.15) %>%
        filter(abs(delta.k) < 15) %>%
        ungroup()

    df.inversion %>%
        ggplot(aes(x = zpbl)) +
        geom_histogram(aes(fill = ytype), position = "dodge") +
        theme_bw(24)

    df.inversion %>%
        ggplot(aes(x = delta.zpbl)) +
        geom_histogram() +
        theme_bw(24)

    ## plot theta and q inversions 
    df.inversion %>%
        tidyr::gather(xtype, x, q, theta, xl, qt, theta.l) %>%
        group_by(delta.k, xtype, ytype) %>%
        summarize(sd = sd(x, na.rm = TRUE),
                  med.x = median(x),
                  x.25 = quantile(x, 0.25),
                  x.75 = quantile(x, 0.75),
                  x = mean(x),
                  delta.p = median(delta.p)) %>%
        ##     ggplot(aes(y = x, x = delta.p * 1e-2, ymin = x - sd, ymax = x + sd)) +
        ggplot(aes(y = med.x, x = delta.p * 1e-2, ymin = x.25, ymax = x.75)) +
        geom_ribbon(fill = "gray", alpha = 0.5) +
        geom_line() +
        geom_point() +
        facet_grid(ytype ~ xtype, scales = "free_x") +
        coord_flip() +
        scale_x_reverse() +
        theme_bw(24)

    df.inversion %>%
        ungroup() %>%
        filter(time == time[1]) %>%
        tidyr::gather(xtype, x, q, theta) %>%
        group_by(lon, lat, time) %>%
        ggplot(aes(y = x, x = delta.k, ##delta.p * 1e-2,
                   alpha = paste(time, lon, lat, sep = "_"))) +
        scale_alpha_discrete(range = c(0.1,0.1), guide = "none") +
        geom_line(col = "gray", lwd = 2) +
        facet_grid(. ~ xtype, scales = "free_x") +
        coord_flip() +
        scale_x_reverse() +
        theme_bw(24)

    ## plot theta and q profiles relative to boundary-layer mean
    df.inversion.relative <- df.inversion %>%
        tidyr::gather(xtype, x, q, theta, xl, qt, theta.l) %>%
        group_by(lon, lat, time, xtype, ytype) %>%
        mutate(x.pbl = mean(x[delta.k > 0]),
               delta.x = x - ifelse(xtype == "xl",
                                    0, x.pbl),
               max.delta.k = max(delta.k),
               k.norm = delta.k / max(delta.k),
               p.norm = delta.p / (aps - p.pbl)) %>%
        ungroup()

    set.seed(0xbeef)
    times <- sample(unique(df.inversion.relative$time), 10)

    tikz_sanitize <- function(x) cbasetools::sanitize.numbers(x, "latex", TRUE, TRUE)
    tikz_sanitize_sparse <- function(x) {
        ret <- tikz_sanitize(x)
        str(ret)
        ret[seq(2, length(labels), by = 2)] <- ""
        ret
    }

                                        # change the default scales
    scale_x_continuous <- function(..., labels=tikz_sanitize)
        ggplot2::scale_x_continuous(..., labels=labels)

    scale_x_reverse <- function(..., labels=tikz_sanitize)
        ggplot2::scale_x_reverse(..., labels=labels)

    scale_y_continuous <- function(..., labels=tikz_sanitize_sparse)
        ggplot2::scale_y_continuous(..., labels=labels)

    df.inversion.relative %>%
        ungroup() %>%
        filter(jk.pbl < 46) %>%
        filter(time %in% times) %>%
        mutate(xtype = revalue(xtype, c(q = "$q$",
                                        qt = "$q_t$",
                                        theta = "$\\theta$",
                                        theta.l = "$\\theta_l$",
                                        xl = "$r_l$")),
               ytype = revalue(ytype, c(zinv = "$z_\\text{inv}$",
                                        zpbl = "$z_\\text{pbl}$"))) %>%
        ggplot(aes(y = delta.x,  x = delta.k,
                   ## k.norm, ##
                   ## x = p.norm,
                   ## x = delta.p * 1e-2,
                   ## col = p.pbl,
                   col = factor(jk.pbl),
                   ## col = max.delta.k,
                   ## col = t.surf,
                   ## col = LTS,
                   ## col = omega.500,
                   ## col = omega.700 * 86400 * 1e-2,
                   alpha = paste(time, lon, lat, sep = "_"))) +
        scale_alpha_discrete(range = c(0.1,0.1), guide = "none") +
        ## scale_color_distiller(palette = "Spectral") +
        scale_color_brewer(palette = "Spectral") +
        guides(color = guide_legend(override.aes = list(lwd = 1, alpha = 1))) +
        geom_line(##col = "gray",
            lwd = 0.25) +
        facet_grid(ytype ~ xtype, scales = "free_x") +
        ## coord_flip(xlim = c(-2.5, 1)) +
        coord_flip(xlim = c(-10, 5)) +
        scale_x_reverse() +
        theme_bw(24)

    df.inversion.relative %>%
        ungroup() %>%
        ## filter(time == time[1]) %>%
        tidyr::gather(xtype, x, delta.q, delta.theta) %>%
        group_by(lon, lat, time) %>%
        ## sample_frac(0.1) %>%
        ggplot(aes(y = x, x = ## delta.k,
                              p.norm,
                   ## k.norm, ##delta.p * 1e-2,
                   ## col = p.pbl,
                   ## col = max.delta.k,
                   ## col = t.surf,
                   ## col = zpbl,
                   ## col = LTS,
                   ## col = omega.500,
                   ## col = omega.700 * 86400 * 1e-2,
                   alpha = paste(time, lon, lat, sep = "_"))) +
        scale_alpha_discrete(range = c(0.01,0.01), guide = "none") +
        scale_color_distiller(palette = "Spectral") +
        geom_line(##col = "gray",
            lwd = 0.1) +
        facet_grid(. ~ xtype, scales = "free_x") +
        coord_flip(xlim = c(-2.5, 1)) +
        ## coord_flip(xlim = c(-10, 5)) +
        scale_x_reverse() +
        theme_bw(24)

    df.inversion.relative %>%
        ungroup() %>%
        ## filter(time == time[1]) %>%
        mutate(omega = omega * 86400 * 1e-3) %>%
        tidyr::gather(xtype, x, delta.q, delta.theta, omega) %>%
        group_by(lon, lat, time) %>%
        ## sample_frac(0.1) %>%
        ggplot(aes(y = x, x = delta.k,
                   ## p.norm,
                   ## k.norm, ##delta.p * 1e-2,
                   ## col = p.pbl,
                   ## col = max.delta.k,
                   col = lat,
                   ## col = factor(zpbl),
                   ## col = LTS,
                   ## col = omega.500,
                   ## col = omega.700 * 86400 * 1e-2,
                   alpha = paste(time, lon, lat, sep = "_"))) +
        scale_alpha_discrete(range = c(0.01,0.01), guide = "none") +
        scale_color_distiller(palette = "Spectral") +
        geom_line(##col = "gray",
            lwd = 0.1) +
        facet_grid(. ~ xtype, scales = "free_x") +
        ## coord_flip(xlim = c(-2.5, 1)) +
        coord_flip(xlim = c(-10, 5)) +
        scale_x_reverse() +
        theme_bw(24)



    nc.w500 <- nc_open("/scratch/b/b380126/mpiesm-1.2.00p1/src/echam/experiments/amip-tdiag/amip-tdiag_197901.01_omega.nc")
    w500 <- ncvar_get(nc.w500, "omega", NA, c(-1,-1,-1,192))

    rm(tm1)
    rm(qm1)
    rm(omega)
}

#' @export
plot.sample2 <- function() {
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(plotutils)

    df <- readRDS("mltools/ml-1979.rds")

    df %>%
        group_by(lon, lat) %>%
        summarize(fscu = n() / (365 * 24 * 8)) %>%
        ungroup() %>%
        discretize(fscu, seq(0, 1, 0.1), as_factor = TRUE) %>%
        mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
        ggplot(aes(x = lon, y = lat, fill = fscu)) +
        geom_raster() +
        geom_world_polygon() +
        scale_x_geo() +
        scale_y_geo() +
        coord_fixed(1, ylim = c(-45, 45)) +
        scale_fill_brewer(palette = "Blues", direction = -1) +
        ## scale_color_brewer(palette = "Blues", direction = -1) +
        ## geom_contour(aes(z = fscu, lty = factor(..level..)), binwidth = 0.2) +
        ## scale_linetype("SCu fraction") +
        theme_bw(24)

    df %>%
        group_by(lon, lat) %>%
        summarize(fscu = n() / (365 * 24 * 8)) %>%
        ungroup() %>%
        mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
        ## discretize(fscu, seq(0, 1, 0.1), as_factor = TRUE) %>%
        ggplot(aes(x = lon, y = lat, fill = fscu)) +
        ## geom_raster() +
        geom_world_polygon() +
        scale_x_geo() +
        scale_y_geo() +
        coord_fixed(1, ylim = c(-45, 45)) +
        ## scale_fill_brewer(palette = "Blues", direction = -1) +
        ## scale_color_brewer(palette = "Blues", direction = -1) +
        geom_contour(aes(z = fscu, lty = factor(..level..)), binwidth = 0.2) +
        scale_linetype("SCu fraction") +
        theme_bw(24)

    entrain <- df %>%
        ungroup() %>%
        filter(lh != 0, sh != 0) %>%
        ## approximate boundary layer thickness (m)
        mutate(H = (aps - p.pbl) * 0.1) %>%
        ## calculate entrainment for theta_l and for q_t
        transmute(lon, lat, time,
                  E_theta = (1.3 * H * dtheta_l.0.dt + sh / 1e3) / (theta_l.plus - theta_l.0),
                  E_q = (1.3 * H * dq_t.0.dt + lh / 2.5e6) / (q_t.plus - q_t.0))

    entrain %>%
        tidyr::gather(var, E, E_theta, E_q) %>%
        group_by(lon, lat, var) %>%
        summarize(E = median(E, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(E = pmin(pmax(E, -0.01), 0.01)) %>%
        mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
        ggplot(aes(x = lon, y = lat, fill = E)) +
        geom_raster() +
        geom_world_polygon() +
        scale_x_geo() +
        scale_y_geo() +
        coord_cartesian(ylim = c(-60, 60)) +
        scale_fill_distiller(palette = "RdYlBu") +
        facet_grid(. ~ var) +
        theme_bw(24)


    doParallel::registerDoParallel(12)
    entrain <- df %>%
        plyr::ddply(~ lon + lat + time, function(x) {
            x %>%
                ## approximate boundary layer thickness (m)
                mutate(H = (aps - p.pbl) * 0.1) %>%
                ## calculate entrainment for theta_l and for q_t
                transmute(lon, lat, time,
                          E_theta = (H * dtheta_l.0.dt + sh / 1e3) / (theta_l.plus - theta_l.0),
                          E_q = (H * dq_t.0.dt + lh / 2.5e6) / (q_t.plus - q_t.0))
        }, .parallel = TRUE)
    
}
