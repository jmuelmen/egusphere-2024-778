#' @export
#' 
process_aerocom <- function(path = "/global/homes/j/jmuelmen/scratch/cms/aerocom") {
    library(ncdf4)
    library(plotutils)
    library(magrittr)
    library(dplyr)
    library(parallel)

    gc()
    Sys.setenv(TZ = "UTC")

    doParallel::registerDoParallel(4)
    
    models <- list.files(path, recursive = FALSE)
    models <- rev(models)
    plyr::l_ply(models, function(model) {
        fname <- function(path, model, var) {
            sprintf("%s/%s/%s_%s.nc",
                    path,
                    model,
                    var,
                    gsub("_IND3", "", model))
        }
        nc_get_df <- function(var, model, spread = TRUE, mask = NULL, step) {
            gc()
            nc <- try(nc_open(fname(path, model, var)))
            if ("try-error" %in% class(nc))
                return(NULL)
            on.exit(nc_close(nc))
            nbdate <- ncatt_get(nc, "time")$"units"
            ret <- nc.to.df(nc, var, spread = spread, mask = mask, start = c(1, step, 1), count = c(-1, 1, -1)) %>%
                dplyr::mutate(time = as.POSIXct(nbdate, format = "days since %Y-%m-%d %H:%M:%S") + time * 86400)
            ## str(ret)
            ret
        }
        print(model)
        ## mccollect(mcparallel( {
            nc <- try(nc_open(fname(path, model, "iwp")))
            lat <- ncvar_get(nc, "lat")
            ret <- plyr::ldply(1 : length(lat), function(step) {
                if (model != "HadGEM3-renamed_IND3") {
                    mccollect(mcparallel( {
                        print(lat[step])
                        df.iwp <- nc_get_df("iwp", model, step = step)
                        mask.iwp <- df.iwp$iwp < 1e-3
                        ## print(sum(mask.iwp))
                        ret <- plyr::ldply(c("cdnc", "lwp", "iwp", "tcc"), nc_get_df, model = model, spread = FALSE, mask = mask.iwp, step = step) %>%
                            tidyr::spread(name, .var_val) %>%
                            mutate(model = model)
                        ret
                    }))[[1]]
                } else {
                    print(lat[step])
                    ret <- plyr::ldply(c("cdnc", "lwp", "iwp", "tcc"), nc_get_df, model = model, spread = FALSE, step = step)
                    ## str(ret)
                    ## print(ret %$% table(time, name))
                    ret %<>%
                        group_by(time, name) %>%
                        mutate(n = n()) %>%
                        group_by(time) %>%
                        mutate(valid = all(n == n[1])) %>%
                        ungroup() %>%
                        filter(valid) %>%
                        tidyr::spread(name, .var_val) %>%
                        filter(is.finite(iwp)) %>%
                        filter(iwp < 1e-3) %>%
                        mutate(model = model)
                    str(ret)
                    ret
                }
            }, .parallel = TRUE)
            str(ret)
            parallel::mcparallel(saveRDS(ret, sprintf("%s_lwp_cdnc.rds", model)), detached = TRUE)
        ## }))
    }, .parallel = FALSE)
}
