#' 
#'
#' @export
process.gfdl.lwp.cdnc <- function(fname = "~/scratch/cms/am4/nc4/atmos.2010010100-2010123123.%s.nc") {
    library(ncdf4)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(doParallel)
    library(magrittr)
    doParallel::registerDoParallel(12)
    
    Sys.setenv(TZ = "UTC")

    ## preliminaries
    nc <- nc_open(sprintf(fname, "LWP"))
    ## lat <- ncvar_get(nc, "lat")
    time <- ncvar_get(nc, "time")
    start.date <- ncatt_get(nc, "time")$"units"
    nc_close(nc)

    getstep <- function(name, step, mask, is.3d) {
        gc()
        nc <- nc_open(sprintf(fname, name))
        on.exit(nc_close(nc))
        plotutils::nc.to.df(nc, name,
                            spread = FALSE,
                            mask = mask,
                            start = if (is.3d) c(1, 1, 1, step) else c(1, 1, step),
                            count = if (is.3d) c(-1, -1, -1, 1) else c(-1, -1, 1))
    }
    
    df <- plyr::ldply(1 : length(time), function(step) {
        df.2d <- plyr::ldply(c("IWP",    
                               "lclmodis",
                               "lremodis",
                               "ltaumodis",
                               "LWP",    
                               "precip", 
                               "ps",     
                               "u_ref",  
                               "v_ref",  
                               "z_pbl",  
                               "u700",   
                               "v700")[1:6],
                             getstep, step, NULL, FALSE) %>%
            tidyr::spread(name, .var_val)

        mask <- df.2d %$% (IWP <= 0 & is.finite(lclmodis) & lclmodis > 0)

        str(df.2d)

        df.3d <- plyr::ldply(c("cld_amt",
                               "droplets",
                               "liq_wat",
                               "omega",
                               "temp"),
                             getstep, step, NULL, TRUE) %>%
            tidyr::spread(name, .var_val) %>%
            group_by(lon, lat, time) %>%
            ## mutate(jk = 1 : length(pfull)) %>%
            summarize(jk.500 = which.min(abs(pfull - 500)),
                      jk.700 = which.min(abs(pfull - 700)),
                      jk.surface = which.max(pfull),
                      jk.liqtop = min(which(cld_amt > 0.1 & droplets > 0 & liq_wat > 0)),
                      omega500 = omega[jk.500],
                      omega700 = omega[jk.700],
                      t700 = temp[jk.700],
                      ts = temp[jk.surface],
                      ttop = temp[jk.liqtop],
                      cdnc = droplets[jk.liqtop],
                      cldtop = cld_amt[jk.liqtop],
                      liqtop = liq_wat[jk.liqtop]) %>%
            ungroup()

        str(df.3d)

        left_join(df.2d %>% filter(IWP <= 0 & is.finite(lclmodis) & lclmodis > 0),
                  df.3d,
                  by = c("lon", "lat", "time"))
        
    }, .parallel = TRUE)

}        
