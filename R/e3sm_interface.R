#' @useDynLib mltools entrainment_diags_eam

#' @export
entrainment.diags <- function(pver,
                              p, t, qv, ql,
                              dpdt, dtdt, dqvdt, dqldt, 
                              dtdt_adv, dqvdt_adv, dqldt_adv, 
                              dtdt_radheat, 
                              qauto, qacc, 
                              hflux_srf, qvflx_srf, 
                              kpblh) {
    ## res <- .Fortran("entrainment_diags",
    ##                 as.integer(pver),
    ##                 as.double(p), as.double(t), as.double(qv), as.double(ql),
    ##                 as.double(dpdt), as.double(dtdt), as.double(dqvdt), as.double(dqldt), 
    ##                 as.double(dtdt_adv), as.double(dqvdt_adv), as.double(dqldt_adv), 
    ##                 as.double(dtdt_radheat), 
    ##                 as.double(qauto), as.double(qacc), 
    ##                 as.double(hflux_srf), as.double(qvflx_srf), 
    ##                 as.integer(kpblh),
    ##                 Etheta = double(1),
    ##                 Eq = double(1))
    return(NULL)
}

