#' Calculate start and count values for blockwise traversal of EAMxx (i.e., very
#' large) NetCDF files
#'
#' Only intended for 3D fields with time, ncol, lev dimension for now
#'
#' @return data.frame with start and count for time, ncol, lev
eamxx_step <- function(nc, varname, blk.time = 1, blk.ncol = 1024^2, blk.lev = 128) {
    sizes <- nc$var[[varname]]$size
    start.time <- seq(1, sizes[3], blk.time)
    start.ncol <- seq(1, sizes[2], blk.ncol)
    start.lev <- seq(1, sizes[1], blk.lev)

    expand.grid(start.time = start.time, start.ncol = start.ncol, start.lev = start.lev,
                count.time = blk.time, count.ncol = blk.ncol, count.lev = blk.lev)
}

#' Calculate start and count values for blockwise traversal of EAMxx (i.e., very
#' large) NetCDF files 
#'
#' Only intended for 3D fields with time, ncol, lev dimension for now
#'
#' @param ncols Vector of ncols to select; each entry turns into a size-1 block
#'
#' @return data.frame with start and count for time, ncol, lev
eamxx_step_ncol_select <- function(nc, varname, blk.time = 1, ncols, blk.lev = 128) {
    sizes <- nc$var[[varname]]$size
    start.time <- seq(1, sizes[3], blk.time)
    start.lev <- seq(1, sizes[1], blk.lev)

    expand.grid(start.time = start.time, start.ncol = ncols, start.lev = start.lev,
                count.time = blk.time, count.ncol = 1, count.lev = blk.lev,
                KEEP.OUT.ATTRS = FALSE)
}

#' Reduce fields stepwise
eamxx_reduce <- function(nc, varnames.3d, varnames.2d = NULL, eamxx.steps, reduce) {
    plyr::ddply(eamxx.steps, ~ start.time + start.ncol + start.lev, function(step) {
        ## print(step)
        df.3d <- step %$% plotutils::nc.to.df(nc, varnames.3d,
                                              start = c(start.lev, start.ncol, start.time),
                                              count = c(count.lev, count.ncol, count.time))
        df.2d <- if (!is.null(varnames.2d)) {
                     step %$% plotutils::nc.to.df(nc, varnames.2d,
                                                  start = c(start.ncol, start.time),
                                                  count = c(count.ncol, count.time))
                 } else {
                     NULL
                 }
        reduce(df.3d, df.2d)
    }, .progress = "text")
}

#' Reduce fields stepwise
#'
#' Supply a file name instead of an ncdf4 object, since each thread will need
#' its own file handle
eamxx_reduce_parallel <- function(nc, varnames.3d, varnames.2d = NULL, eamxx.steps, reduce) {
    plyr::ddply(eamxx.steps, ~ start.time + start.ncol + start.lev, function(step) {
        print(step)
        nc <- ncdf4::nc_open(nc)
        df.3d <- step %$% plotutils::nc.to.df(nc, varnames.3d,
                                              start = c(start.lev, start.ncol, start.time),
                                              count = c(count.lev, count.ncol, count.time))
        df.2d <- if (!is.null(varnames.2d)) {
                     step %$% plotutils::nc.to.df(nc, varnames.2d,
                                                  start = c(start.ncol, start.time),
                                                  count = c(count.ncol, count.time))
                 } else {
                     NULL
                 }
        nc_close(nc)
        reduce(df.3d, df.2d)
    }, .parallel = TRUE)
}

process_eamxx_profiles <- function() {
    library(ncdf4)
    library(plotutils)
    library(magrittr)
    library(plyr)
    library(dplyr)
    library(ggplot2)
    library(mltools)
    
    nc <- nc_open("~/scratch.pm/ForACIdiags/output.scream.daily.AVERAGE.ndays_x1.2019-08-04-00000.nc")
    hybm <- nc.to.df(nc, c("hyam", "hybm"))
    hybi <- nc.to.df(nc, c("hyai", "hybi"))
    hybm %<>% mutate(dhyai = hybi %$% diff(hyai),
                     dhybi = hybi %$% diff(hybi))
    lonlat <- nc.to.df(nc, c("lon", "lat"))
    
    lev.500 <- with(hybm, lev[which.min(abs(hybm - 0.5))])
    lev.700 <- with(hybm, lev[which.min(abs(hybm - 0.7))])

    ## find columns meeting Medeiros & Stevens Sc criteria
    step.daily <- eamxx_step(nc, "omega")
    df.sc <- try(readRDS("~/scratch.pm/ForACIdiags/scream.daily.sc.rds"))
    if (any(class(df.sc) == "try-error")) {
        df.sc <- eamxx_reduce(nc, c("omega", "T_mid"), NULL, step.daily,
                              function(df, df.2d) { ## only 3D fields in this pass
                                  str(df)
                                  df %<>%
                                      group_by(time, ncol) %>%
                                      summarize(omega.500 = omega[lev.500],
                                                omega.700 = omega[lev.700],
                                                lts = T_mid[lev.700] * 0.7^(-2/7) - T_mid[128]) %>%
                                      ungroup() %>%
                                      filter(omega.500 * 864 > 10,
                                             omega.700 * 864 > 10,
                                             lts > 18.55)
                                  str(df)
                                  df
                              })
        saveRDS(df.sc, "~/scratch.pm/ForACIdiags/scream.daily.sc.rds")
    }
    df.sc %<>% select(ncol)

    ## find inversion based on daily T and q_v profiles
    step.daily.sc <- eamxx_step_ncol_select(nc, "T_mid", ncols = df.sc$ncol)
    step.daily.short <- eamxx_step(nc, "omega", blk.ncol = 1024)
    doParallel::registerDoParallel(cores = 12)
    df.profiles.daily <- try(readRDS("~/scratch.pm/ForACIdiags/scream.daily.profiles.rds"))
    if (any(class(df.profiles.daily) == "try-error")) {
        df.profiles.daily <-
            eamxx_reduce_parallel("~/scratch.pm/ForACIdiags/output.scream.daily.AVERAGE.ndays_x1.2019-08-04-00000.nc",
                                  c("qc", "qv", "T_mid"), c("ps"), step.daily,
                                  function(df.3d, df.2d) {
                                      str(df.2d)
                                      df <- df.sc %>%
                                          inner_join(df.3d) %>%
                                          left_join(df.2d, by = c("time", "ncol")) %>%
                                          left_join(hybm, by = c("lev"))
                                      str(df)
                                      df %<>%
                                          filter(ps > 975e2) %>%
                                          group_by(time, ncol) %>%
                                          mutate(p0 = 1e5) %>%
                                          calculate_p3d() %>%
                                          rename_fields_eamxx() %>%
                                          calculate_thermodynamics() %>%
                                          find_inversion_eamxx() %>%
                                          ungroup() %>%
                                          filter(lev >= pmin(jk.theta, jk.qt) - 20) %>%
                                          select(-c(ps : p0))
                                      str(df)
                                      df
                                  })
        saveRDS(df.profiles.daily, "~/scratch.pm/ForACIdiags/scream.daily.profiles.rds")
    }
    df.profiles.daily %<>%
        select(-c(start.time : start.lev, time)) %>% 
        mutate(jk.tinv = as.integer(jk.tinv),
               jk.theta = as.integer(jk.theta),
               jk.qt = as.integer(jk.qt)) %>%
        select(-jk.pbl)

    df.profiles.jk.pbl <- df.profiles.daily %>%
        filter(abs(jk.qt - jk.theta) < 5) %>%
        mutate(delta.k.qt = lev - jk.qt) %>%
        mutate(delta.k.theta = lev - jk.theta) %>%
        select(ql, qv, theta, theta_l, delta.k.qt, delta.k.theta, jk.theta, jk.qt) %>%
        tidyr::gather(vert, delta.k, delta.k.qt, delta.k.theta) %>%
        tidyr::gather(var, val, ql, qv, theta, theta_l) 

    gc()

    df.daily.statistics <- df.profiles.jk.pbl %>%
        mutate(jk.pbl = ifelse(vert == "delta.k.qt", jk.qt, jk.theta)) %>%
        group_by(vert, var, delta.k, jk.pbl) %>%
        mutate(val.nonzero = ifelse(val == 0, NA, val)) %>%
        summarize(across(c(val, val.nonzero),
                         list(mean = ~ mean(.x, na.rm = TRUE),
                              sd = ~ sd(.x, na.rm = TRUE),
                              q = ~ {
                                  ret <- transform(quantile(.x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                                                            na.rm = TRUE)) 
                                  mutate(ret, names = sprintf("p_%s", row.names(ret))) %>%
                                      tidyr::spread(names, X_data)
                              },
                              n = ~ sum(!is.na(.x))),
                         .unpack = TRUE)) %>%
        ungroup() %>%
        tidyr::unnest_wider(q) %>%
        rename_with(. %>% gsub("^([0-9]*)%$", "p\\1", .))
    
    parallel::mcparallel(df.profiles.daily %>%
                         filter(lev == lev[1]) %>%
                         select(ncol, jk.tinv, jk.theta, jk.qt, qt.jump, theta.jump) %>%
                         saveRDS("~/scratch.pm/ForACIdiags/scream.daily.inv.rds"))

    df.profiles.daily %>%
        filter(lev == lev[1]) %$%
        summary(jk.qt - jk.theta)
    
    df.profiles.daily %>%
        filter(lev == lev[1]) %>%
        ## filter(abs(jk.theta - jk.qt) < 5) %>%
        summary()

    df.profiles.daily %>%
        filter(lev == lev[1]) %>%
        filter(abs(jk.theta - jk.qt) < 5) %$%
        summary(jk.theta)

    df.profiles.daily %>%
        filter(lev == lev[1]) %>%
        summary(jk.theta)
    
    df.profiles.daily %>%
        filter(lev == lev[1]) %$%
        summary(theta.jump)

    df.profiles.daily %>%
        filter(lev == lev[1]) %$%
        summary(qt.jump)

    fname.nd <- "~/scratch.pm/ForACIdiags/output.scream.NcCldfrac_liq.INSTANT.nhours_x3.2019-08-04-00000.nc"
    fname.ql <- "~/scratch.pm/ForACIdiags/output.scream.QcQi.INSTANT.nhours_x3.2019-08-04-00000.nc"

    nc.3h <- nc_open(fname.nd)
    step.3h <- eamxx_step(nc.3h, "cldfrac_liq")

    doParallel::registerDoParallel(cores = 8)
    df.profiles.daily %<>% filter(between(jk.theta, 95, 118),
                                  between(jk.qt, 95, 118),
                                  abs(jk.qt - jk.theta) < 5)
    gc()
    df.nd.f.3h <- eamxx_reduce_parallel(fname.nd,
                              c("cldfrac_liq", "nc"), NULL, step.3h,
                              function(df.3d, df.2d) {
                                  inner_join(df.3d, df.profiles.daily, by = c("ncol", "lev"))
                              })
    rm(df.profiles.daily)
    gc()
    mcparallel(saveRDS(df.nd.f.3h, "~/scratch.pm/ForACIdiags/scream.3h.nd.f.rds"))
    df.ncol.lev <- readRDS("~/scratch.pm/ForACIdiags/scream.daily.profiles.rds") %>%
        filter(between(jk.theta, 95, 118),
               between(jk.qt, 95, 118),
               abs(jk.qt - jk.theta) < 5) %>%
        select(ncol, lev)
    gc()
    df.ql.qi.3h <- eamxx_reduce_parallel(fname.ql,
                                         c("qc", "qi"), NULL, step.3h,
                                         function(df.3d, df.2d) {
                                             inner_join(df.3d, df.ncol.lev, by = c("ncol", "lev"))
                                            }) %>%
        select(-c(start.time : start.lev))
    saveRDS(df.ql.qi.3h, "~/scratch.pm/ForACIdiags/scream.3h.ql.qi.rds")
    df.nd.f.3h <- readRDS("~/scratch.pm/ForACIdiags/scream.3h.nd.f.rds") %>%
        select(-c(start.time : start.lev))
    gc()
    df.profiles.3h <- inner_join(readRDS("~/scratch.pm/ForACIdiags/scream.3h.ql.qi.rds"),
                                 readRDS("~/scratch.pm/ForACIdiags/scream.3h.nd.f.rds") %>%
                                 select(-c(start.time : start.lev)))
    saveRDS(df.profiles.3h, "~/scratch.pm/ForACIdiags/scream.3h.profiles.rds")


    df.nd.f.3h <- readRDS("~/scratch.pm/ForACIdiags/scream.3h.nd.f.rds")
    df.nd.f.3h %<>% select(ncol, lev, time, cldfrac_liq, nc)
    df.ql.qi.3h <- readRDS("~/scratch.pm/ForACIdiags/scream.3h.ql.qi.rds")
    df.ql.qi.3h %<>% select(ncol, lev, time, qc, qi)
    gc()
    df.3h <- full_join(df.ql.qi.3h, df.nd.f.3h)
    rm(df.ql.qi.3h)
    rm(df.nd.f.3h)
    gc()
    df.profiles.daily <- readRDS("~/scratch.pm/ForACIdiags/scream.daily.profiles.rds")
    df.3h %<>% full_join(df.profiles.daily, by = c("ncol", "lev"))

    df.3h %<>% select(-c(start.time : time.y)) %>% rename(time = time.x)
    df.3h %<>% mutate(jk.tinv = as.integer(jk.tinv),
                      jk.theta = as.integer(jk.theta),
                      jk.qt = as.integer(jk.qt)) %>%
        select(-jk.pbl)
    
    parallel::mcparallel(saveRDS(df.3h, "~/scratch.pm/ForACIdiags/scream.3h.profiles.rds"))

    df.3h <- readRDS("~/scratch.pm/ForACIdiags/scream.3h.profiles.rds")

    df.3h %<>%
        filter(abs(jk.qt - jk.theta) < 5) %>%
        mutate(delta.k.qt = lev - jk.qt) %>%
        mutate(delta.k.theta = lev - jk.theta) %>%
        select(qc, qi, cldfrac_liq, nc, delta.k.qt, delta.k.theta, jk.theta, jk.qt)
    gc()
    saveRDS(df.3h, "~/scratch.pm/ForACIdiags/scream.3h.profiles.short.rds")

    df.3h <- readRDS("~/scratch.pm/ForACIdiags/scream.3h.profiles.short.rds")

    df.3h %<>% tidyr::gather(vert, delta.k, delta.k.qt, delta.k.theta)
    ## df.3h %<>% tidyr::gather(var, val, qc, qi, cldfrac_liq, nc) 

    df.3h %>%
        mutate(jk.pbl = ifelse(vert == "delta.k.qt", jk.qt, jk.theta)) %>%
        group_by(vert, delta.k, jk.pbl) %>%
        ## mutate(val.nonzero = ifelse(val == 0, NA, val)) %>%
        summarize(across(c(qc, qi, cldfrac_liq, nc ## val ## , val.nonzero
                           ),
                         list(mean = ~ mean(.x, na.rm = TRUE),
                              sd = ~ sd(.x, na.rm = TRUE),
                              q = ~ {
                                  ret <- transform(quantile(.x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                                                            na.rm = TRUE)) 
                                  mutate(ret, names = sprintf("p%s", row.names(ret))) %>%
                                      tidyr::spread(names, X_data)
                              },
                              n = ~ sum(!is.na(.x))),
                         .unpack = "{outer}.{inner}")) %>%
        ungroup() %>%
    saveRDS("~/scratch.pm/ForACIdiags/3h.statistics.rds")

    df.3h %>%
        mutate(jk.pbl = ifelse(vert == "delta.k.qt", jk.qt, jk.theta)) %>%
        group_by(vert, delta.k, jk.pbl) %>%
        ## mutate(val.nonzero = ifelse(val == 0, NA, val)) %>%
        mutate(across(c(qc, qi, cldfrac_liq, nc),
                      ~ ifelse(.x != 0, .x, NA))) %>%
        summarize(across(c(qc, qi, cldfrac_liq, nc ## val ## , val.nonzero
                           ),
                         list(mean = ~ mean(.x, na.rm = TRUE),
                              sd = ~ sd(.x, na.rm = TRUE),
                              q = ~ {
                                  ret <- transform(quantile(.x, c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                                                            na.rm = TRUE)) 
                                  mutate(ret, names = sprintf("p%s", row.names(ret))) %>%
                                      tidyr::spread(names, X_data)
                              },
                              n = ~ sum(!is.na(.x))),
                         .unpack = "{outer}.{inner}")) %>%
        ungroup() %>%
    saveRDS("~/scratch.pm/ForACIdiags/3h.nonzero.statistics.rds")
    
    
    df.3h %>%
        filter(lev == lev[1]) %$%
        table(jk.theta) %>%
        (function(x) cumsum(x) / sum(x))

    df.3h %>%
        filter(lev == lev[1]) %$%
        table(jk.qt) %>%
        (function(x) cumsum(x) / sum(x))

    df.3h %>%
        filter(lev == lev[1]) %$%
        table(jk.theta) %>%
        (function(x) x / sum(x))
    
    df.3h %>%
        filter(lev == lev[1]) %$%
        summary(jk.theta)    

    df.3h %>%
        filter(lev == lev[1]) %$%
        summary(jk.qt)    

    df.3h %>%
        filter(lev == 100) %$%
        summary(p)    

}

    ## eamxx_reduce(nc, c("qc", "qv", "qi", "T_mid"), c("ps"), step.test,
    ##                                   function(df.3d, df.2d) {
    ##                                       str(df.2d)
    ##                                       df <- df.sc %>%
    ##                                           inner_join(df.3d) %>%
    ##                                           left_join(df.2d, by = c("time", "ncol")) %>%
    ##                                           left_join(hybm, by = c("lev"))
    ##                                       str(df)
    ##                                       df %<>%
    ##                                           group_by(time, ncol) %>%
    ##                                           mutate(p0 = 1e5) %>%
    ##                                           calculate_p3d() %>%
    ##                                           rename_fields_eamxx() %>%
    ##                                           calculate_thermodynamics() %>%
    ##                                           find_inversion_eamxx() %>%
    ##                                           ungroup() %>%
    ##                                           filter(lev >= 100) %>%
    ##                                           select(-c(ps : p0))
    ##                                       str(df)
    ##                                       df
    ##                                   })
    
