## @knitr aerocom-load
df.aerocom <- readRDS("../aerocom.rds") %>%
    mutate(model = gsub("_", "\\\\_", model))

## ---- cam-load -------------------------------------------
df.cam <- load_lwp_cam("../cam_cmsaci.002.cam.h1.2015_2d.nc") %>%
    rename(cdnc = ACTNL,
           lcc = CLDTOT,
           prec = PRECT,
           lwp = TGCLDLWP) %>%
    filter(cdnc > 0) %>%
    calculate_incloud()

## ---- cam-filter -----------------------------------------

## ---- giss-load ------------------------------------------
df.giss <- expand.grid(period = c("PD", "PI"),
                       experiment = c("jmu"## , "e25s0"
                                      )) %>%
    ddply(~ period + experiment, function(df) {
        fname <- df %$%
            sprintf("../ModelE3_2x2p5deg_ACI_2010.allmerget705ml_amip_2010nudge_%s_%s.rds",
                    experiment, tolower(period))
        readRDS(fname) %>%
            mutate(lcc = 1e-2 * cldss_2d) %>%
            mutate(cdnc_ic = 1e6 * cdnc_ic,
                   lwp_ic = 1e-3 / 1e-2 * lwp_ic)
    })

## ---- giss-filter ----------------------------------------

## ---- e3sm-load ------------------------------------------
df.e3sm <- bind_rows(
    load_lwp("../EAMv1_entrain_NEP_nudged_pd_cdnc-lwp.cam.h1.2010.nc") %>%
    mutate(period = "PD"),
    load_lwp("../EAMv1_entrain_NEP_nudged_pi_cdnc-lwp.cam.h1.2010.nc") %>%
    mutate(period = "PI")) %>%
    rename(lts = th7001000) %>%
    filter_warm() %>%
    filter_ocean() %>%
    calculate_incloud()

## ---- e3sm-load-precip ------------------------------------------
df.e3sm.precip <- readRDS("../e3smv1_aer-EAMv1_nudged_2010.rds") %>%
    rename(lts = th7001000) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower)) %>%
    filter_warm() %>%
    filter_ocean() %>%
    calculate_incloud()

## ---- e3smv2-load ------------------------------------------
df.e3smv2 <- bind_rows(
    readRDS("../v2_2010aer_nudged.rds") %>%
    mutate(period = "PD"),
    readRDS("../v2_1850aer_nudged.rds") %>%
    mutate(period = "PI")) %>%
    rename(lts = th7001000) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower)) %>%
    filter_warm() %>%
    filter_ocean() %>%
    calculate_incloud()

## ---- e3smv2-load-zhang ------------------------------------------
df.zhang <- readRDS("../v2_2010aer_nudged_rad.rds") %>%
    rename(lts = th7001000) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower)) %>%
    filter_warm() %>%
    filter_ocean() %>%
    calculate_incloud() %>%
    ## rough approximation to SZA > 65 degrees
    filter(solin > cos(65 * pi / 180) * 1360) %>% 
    calculate_albedo_from_lcc()

## @knitr e3sm-filter ----------------------------------------
warning("Test")

## ---- gfdl-load ------------------------------------------
df.gfdl <- load_gfdl() %>%
    mutate(cdnc_ic = cdnc_ic * 1e6) %>%
    ## also use better lcc eventually
    mutate(lcc = 1e-2 * lclmodis)

## ---- gfdl-filter ----------------------------------------

## ---- multimodel-setup -----------------------------------
df.multimodel <- bind_rows(
    df.cam %>%
    select_standard() %>%
    mutate(model = "CAM6"),
    df.giss %>%
    filter(experiment == "jmu") %>%
    select_standard() %>%
    ## discretize GISS twice (coarser first --> idempotent) until more data is available
    ## discretize_cdnc(bins = exp(seq(log(1e6), log(300e6), length.out = 25))) %>%
    mutate(model = "GISS"),
    df.gfdl %>%
    select_standard() %>%
    select(-time) %>%
    mutate(model = "GFDL"),
    df.e3sm %>%
    select_standard() %>%
    mutate(model = "E3SM")## ,
    ## df.e3smv2 %>%
    ## select_standard() %>%
    ## mutate(model = "E3SMv2")
) %>%
    filter(lcc > 0.9) %>%
    mutate(period = ifelse(is.na(period), "PD", period)) %>%
    mutate(model = factor(model, levels = c("GISS", "CAM6", "E3SM", "GFDL")))

## @knitr multimodel-plot ------------------------------------
df.multimodel %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional() %>%
    lwp.cdnc.plot(color = model, lty = period, marginal = TRUE)

## @knitr multimodel-plot-lcc-0.1 ------------------------------------
ggplot.multimodel.lcc.0.1 <- bind_rows(
    ## multimodel with lcc > 0.1 threshold
    bind_rows(
        df.cam %>%
        select_standard() %>%
        mutate(model = "CAM6"),
        df.giss %>%
        filter(experiment == "jmu") %>%
        select_standard() %>%
        ## discretize GISS twice (coarser first --> idempotent) until more data is available
        ## discretize_cdnc(bins = exp(seq(log(1e6), log(300e6), length.out = 25))) %>%
        mutate(model = "GISS"),
        df.gfdl %>%
        select_standard() %>%
        select(-time) %>%
        mutate(model = "GFDL"),
        df.e3sm %>%
        select_standard() %>%
        mutate(model = "E3SM")## ,
        ## df.e3smv2 %>%
        ## select_standard() %>%
        ## mutate(model = "E3SMv2")
    ) %>%
    mutate(period = ifelse(is.na(period), "PD", period)) %>%
    mutate(model = factor(model, levels = c("GISS", "CAM6", "E3SM", "GFDL"))) %>%
    filter(period == "PD") %>%
    filter(lcc > 0.1) %>%
    mutate(lcc.lim = 0.1),
    ## multimodel with lcc > 0.9 threshold
    df.multimodel %>%
    filter(period == "PD") %>%
    mutate(lcc.lim = 0.9)) %>%
    group_by(model, lcc.lim) %>%
    lwp.cdnc.conditional() %>%
    mutate(lcc.lim = sprintf("$f_\\text{liq} > %g$", lcc.lim)) %>%
    replace_with_conventional_units() %>%
    lwp.cdnc.plot(color = model, alpha = lcc.lim, cdnc.density.limits = c(0, 0.25),
                  marginal = FALSE, expand = FALSE) +
    scale_color_uscms("") +
    scale_alpha_manual("", values = c(1, 0.33)) +
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"))

## @knitr multimodel-selected-pd-pi-plot ------------------------------------
ggplot.multimodel.selected.pd.pi <- df.multimodel %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional(lwp.bins = exp(seq(log(0.95e-2), log(0.31), length.out = 25 * 2))) %>%
    lwp.cdnc.plot(cdnc.density.limits = c(0, 0.1), color = model, lty = period, marginal = TRUE, expand = FALSE)

## @knitr multimodel-selected-pd-pi-plot-conventional
ggplot.multimodel.selected.pd.pi <- df.multimodel %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional(lwp.bins = exp(seq(log(0.95e-2), log(0.31), length.out = 25 * 2))) %>%
    replace_with_conventional_units() %>%
    lwp.cdnc.plot(cdnc.density.limits = c(0, 0.1), color = model, lty = period, marginal = TRUE, expand = FALSE)

## @knitr multimodel-pd-plot ------------------------------------
ggplot.multimodel.pd <- df.multimodel %>%
    filter(period == "PD") %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional(lwp.bins = exp(seq(log(1e-4), log(3), length.out = 25 * 2))) %>%
    lwp.cdnc.plot(cdnc.density.limits = c(0, 0.1), color = model, marginal = TRUE)

## @knitr multimodel-pd-plot-conventional ------------------------------------
ggplot.multimodel.pd <- df.multimodel %>%
    filter(period == "PD") %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional(lwp.bins = exp(seq(log(1e-4), log(3), length.out = 25 * 2))) %>%
    replace_with_conventional_units() %>%
    lwp.cdnc.plot(cdnc.density.limits = c(0, 0.1), color = model, marginal = TRUE)

## @knitr multimodel-plot-heatmap ------------------------------------
df.multimodel %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional() %>%
    filter(period == "PD") %>%
    lwp.cdnc.plot(color = model, lty = period, heatmap = TRUE) +
    facet_wrap(~model) +
    scale_fill_distiller("$p(\\lwp|\\nd)$", palette = "Spectral", limits = c(0, 0.2), trans = "sqrt") +
    coord_cartesian(ylim = c(0.01, 0.3))

## @knitr multimodel-plot-lwp
df.multimodel %>%
    ## filter(is.scu) %>%
    filter(lwp_ic > 1e-4) %>%
    ggplot(aes(lwp_ic, col = model, lty = period)) +
    geom_freqpoly(aes(y = ..density..), ## aes(y = abs(x) * ..density..), 
        bins = 100) +
    ## stat_ecdf() +
    scale_x_log10() + ## + facet_grid(factor(sed) ~ .)
    labs(x = "$\\mathcal{L}$ (in cloud, kg~m$^{-2}$)") +
    theme(legend.position = c(0.4, 0.25), legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(0, "lines"),
          legend.direction = "vertical", legend.box = "horizontal") 

## @knitr multimodel-selected-plot-lwp
df.multimodel %>%
    ## filter(is.scu) %>%
    filter(lwp_ic > 1e-4) %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    ggplot(aes(lwp_ic, col = model, lty = period)) +
    geom_freqpoly(aes(y = ..density..), ## aes(y = abs(x) * ..density..), 
        bins = 100) +
    ## stat_ecdf() +
    scale_x_log10() + ## + facet_grid(factor(sed) ~ .)
    labs(x = "$\\mathcal{L}$ (in cloud, kg~m$^{-2}$)") +
    theme(legend.position = c(0.4, 0.25), legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(0, "lines"),
          legend.direction = "vertical", legend.box = "horizontal") 

## @knitr multimodel-selected-ecdf-lwp
df.multimodel %>%
    ## filter(is.scu) %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    ggplot(aes(lwp_ic, col = model, lty = period)) +
    ## geom_freqpoly(aes(y = ..density..), ## aes(y = abs(x) * ..density..), 
    ##     bins = 100) +
    stat_ecdf(n = 200, geom = "line") +
    scale_x_log10() + ## + facet_grid(factor(sed) ~ .)
    labs(x = "$\\mathcal{L}$ (in cloud, kg~m$^{-2}$)") +
    coord_cartesian(xlim = c(1e-3, 1)) +
    theme(legend.position = c(0.4, 0.25), legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(0, "lines"),
          legend.direction = "vertical", legend.box = "horizontal") 

## @knitr multimodel-selected-stats
df.multimodel %>%
    ## filter(is.scu) %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    group_by(model, period) %>%
    summarize(lwp_ic = exp(mean(log(lwp_ic)))) %>%
    group_by(model) %>%
    mutate(ratio = lwp_ic[period == "PD"] / lwp_ic[period == "PI"]) %>%
    ungroup() 

## @knitr multimodel-selected-stats
df.multimodel %>%
    ## filter(is.scu) %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    group_by(model, period) %>%
    summarize(lwp_ic = mean(lwp_ic)) %>%
    group_by(model) %>%
    mutate(ratio = lwp_ic[period == "PD"] / lwp_ic[period == "PI"]) %>%
    ungroup() 

## @knitr e3sm-plot ------------------------------------
df.multimodel %>%
    filter(grepl("E3SM", model)) %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional() %>%
    lwp.cdnc.plot(color = model, lty = period, marginal = TRUE)

## @knitr aerocom-plot-iwp-limits ------------------------------------
df.aerocom %>%
    group_by(model, class, iwp.limit) %>%
    filter(class == "ovc") %>%
    ## adjust axis ranges
    filter(between(lwp_ic, 10e-3, 1/3), between(lwp_ic_mean, 10e-3, 1/3), between(cdnc_ic, 3e6, 300e6)) %>%
    lwp.cdnc.plot(color = model, lty = factor(iwp.limit), marginal = TRUE)

## @knitr aerocom-plot-cloudiness ------------------------------------
df.aerocom %>%
    group_by(model, class, iwp.limit) %>%
    filter(iwp.limit == 1e-3) %>%
    lwp.cdnc.plot(color = model, lty = class, marginal = TRUE)

## @knitr aerocom-plot-heatmap ------------------------------------
df.aerocom %>%
    group_by(model, class, iwp.limit) %>%
    filter(iwp.limit == 1e-3) %>%
    lwp.cdnc.plot(color = model, lty = class, heatmap = TRUE) +
    facet_grid(class ~ model) +
    scale_fill_distiller(palette = "BuPu", direction = 1, limits = c(0, 0.2))

## @knitr aerocom-line-plot ------------------------------------
ggplot.aerocom.line.plot <- df.aerocom %>%
    mutate(model = gsub("\\\\_IND3$", "", model)) %>%
    mutate(model = gsub("-ETHZ$|^NCAR-|-renamed", "", model)) %>%
    mutate(model = gsub("-CLUBB", "-CL", model)) %>%
    group_by(model, class, iwp.limit) %>%
    filter(iwp.limit == 1e-3) %>%
    filter(class == "ovc") %>%
    filter(cdnc_ic < 3e8) %>%
    lwp.cdnc.plot(color = model, marginal = TRUE) +
    scale_fill_distiller(palette = "BuPu", direction = 1, limits = c(0, 0.2))

## @knitr aerocom-line-plot-conventional ------------------------------------
ggplot.aerocom.line.plot <- df.aerocom %>%
    mutate(model = gsub("\\\\_IND3$", "", model)) %>%
    mutate(model = gsub("-ETHZ$|^NCAR-|-renamed", "", model)) %>%
    mutate(model = gsub("-CLUBB", "-CL", model)) %>%
    group_by(model, class, iwp.limit) %>%
    filter(iwp.limit == 1e-3) %>%
    filter(class == "ovc") %>%
    filter(cdnc_ic < 3e8) %>%
    replace_with_conventional_units() %>%
    lwp.cdnc.plot(color = model, marginal = TRUE) +
    scale_fill_distiller(palette = "BuPu", direction = 1, limits = c(0, 0.2))

## @knitr giss-ppe-plot
df.giss %>%
    filter(lcc > 0.9) %>%
    group_by(experiment, period) %>%
    lwp.cdnc.conditional() %>%
    filter(between(cdnc_ic, 3e6, 300e6)) %>%
    lwp.cdnc.plot(color = experiment, lty = period, marginal = TRUE) &
    scale_color_brewer(palette = "Dark2")

## @knitr giss-hist
df.giss %>%
    filter(lcc > 0.9) %>%
    ggplot(aes(lwp_ic, color = experiment, lty = period)) +
    geom_density() +
    scale_x_log10()

## @knitr giss-stats
df.giss %>%
    ## filter(lcc > 0.9) %>%
    group_by(experiment, period) %>%
    summarize(lwp.mean = mean(lwp_ic),
              lwp.geomean = exp(mean(log(lwp_ic))),
              lwp.media = median(lwp_ic))

## @knitr multimodel-plot-pd-pi
df.multimodel %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional() %>%
    lwp.cdnc.plot(color = model, lty = period, marginal = TRUE, legend.position = c(0.5, 0.25))

## @knitr e3sm-load-profiles
df.profiles.fscu <- readRDS("../profiles_Sc.rds") %>%
    group_by(ncol) %>%
    mutate(fscu = n() / ((diff(range(ilev)) + 1) * 365 * 8)) %>%
    ungroup() %>%
    filter(fscu > 0.2)

df.profiles.jk.pbl <- df.profiles.fscu %>%
    rename_with(tolower) %>%
    mutate(numliq = numliq * 100 * lev / (287 * t)) %>% ## convert mixing ratio to concentration
    group_by(time, ncol) %>%
    mutate(jk.pbl = ilev[max(which(diff(t) < 0))]) %>%
    mutate(delta.k = ilev - jk.pbl) %>%
    ## filter(all(cloud[delta.k < -1 ] < 1e-3)) %>%
    ungroup() %>%
    mutate(jk.pbl = 72 - jk.pbl) %>% ## cut(72 - jk.pbl, c(0,5,10,20,30,40,80), right = FALSE)) %>%
    filter(jk.pbl < 20)

df.profiles <- df.profiles.jk.pbl %>%
    filter(delta.k >= -20) %>%
    mutate(f_cdnc = numliq / cdnc,
           numliq_ic = ifelse(cloud > 0.1, numliq / cloud, NA),
           f_cdnc_ic = numliq_ic / ifelse(lcc > 0.1, cdnc / lcc, NA),
           ql_ic = ifelse(cloud > 0.1, cldliq / cloud, NA)) %>%
    mutate(t.c = t - 273.15) %>%
    select(jk.pbl : t.c, cloud) %>%
    gather(xtype, x, f_cdnc, numliq_ic, f_cdnc_ic, ql_ic, cloud, t.c) %>%
    mutate(xtype = revalue(xtype, c(t.c = "$T$ (\\textdegree C)",
                                    f_cdnc = "$N_d / N_d^\\text{top}$ (grid mean)",
                                    f_cdnc_ic = "$N_d / N_d^\\text{top}$ (in cloud)",
                                    ql_ic = "$q_l\\ (\\text{kg kg}^{-1})$ (in cloud)",
                                    cloud = "$f$",
                                    numliq_ic = "$N_d$ (in cloud, m$^{-3}$)"))) 

df.profiles.nep <- readRDS("../profiles_NEP.rds")

df.2d <- df.profiles.nep %>%
    dplyr::group_by(time) %>%
    dplyr::summarize(
               ## add regime-identifying variables (omega.500, LTS)
               omega.500 = omega[which.min(abs(p - 500e2))],
               omega.700 = omega[which.min(abs(p - 700e2))],
               LTS = theta[which.min(abs(p - 700e2))] - theta[71],
               ## jk.pbl =  which.min(abs(z3 - pblh)), ## turns out not to work well
               precc = precc[1],
               precl = precl[1]) %>%
    ungroup()

df.scu.e3sm <- df.2d %>%
    filter(LTS > 18.55 & omega.500 * 864 > 10) %>%
    left_join(df.profiles.nep, by = "time", multiple = "all")

## @knitr profiles-ccn
df.scu.e3sm %>%
    group_by(time) %>%
    mutate(jk.pbl = max(which(diff(t) < 0))) %>%
    ungroup() %>%
    mutate(delta.k = jk - jk.pbl) %>%
    filter(delta.k >= -20) %>%
    gather(xtype, x, ccn) %>%
    group_by(time, xtype) %>%
    mutate(delta.x = x - ifelse(xtype %in% c("ccn", "numliq", "numliq_ic", "ql_ic", "t", "ql", "omega", "wp2_clubb", "wp3_clubb", "cloud", "cloudfrac_clubb", "qrl", "qrs", "skw", "prao", "prco"), 0, mean(x[delta.k > 2], na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(xtype = revalue(xtype, tikz_replacements_unitful())) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>% ## cut(72 - jk.pbl, c(0,5,10,20,30,40,80), right = FALSE)) %>%
    filter(jk.pbl %in% 10:15) %>%
    left_join(df.pblh.stats.v1 %>%
              mutate(across(where(is.double), round)),
              by = c("jk.pbl")) %>%
    mutate(pblh = sprintf("%d ($%d\\pm%d$ m)", jk.pbl, mean, sd)) %>%
    group_by(xtype, delta.k, pblh) %>%
    ## ## summarize(q25 = quantile(x, 0.25), q50 = quantile(x, 0.5), q75 = quantile(x, 0.75)) %>%
    summarize(q50 = mean(delta.x, na.rm = TRUE), q25 = mean(delta.x) - 2 * sd(delta.x), q75 = mean(delta.x) + 2 * sd(delta.x), n = n()) %>%
    ungroup() %>%
    filter(!is.na(q50)) %>%
    ggplot(aes(# y = n,
        y = ifelse(n > 100, q50, NA),
        x = delta.k,
        ## x =  delta.k, y = delta.x,
        ## group = time,
        ## k.norm, ##
        ## x = p.norm,
        ## x = delta.p * 1e-2,
        ## col = p.pbl,
        col = pblh,
        fill = pblh)) +
    ## scale_color_distiller(palette = "Spectral") +
    scale_y_continuous(labels = tikz_sanitize_sparse ) +
    ## scale_color_discrete("$k_\\text{sfc} - k_\\text{pbl}$"## , palette = "Spectral"
    ##                      ) +
    ## scale_fill_discrete("$k_\\text{sfc} - k_\\text{pbl}$") +
    scale_color_brewer("PBL depth (model levels)", palette = "Spectral", drop = FALSE,
                       guide = guide_legend(ncol = 2)) +
    ## scale_fill_brewer("$k_\\text{sfc} - k_\\text{pbl}$", palette = "Spectral") +
    ## guides(color = guide_legend("$k_\\text{pbl}$", override.aes = list(lwd = 1, alpha = 1))) +
    geom_line(lwd = 1) +
    ## geom_ribbon(aes(ymin = q25, ymax = q75, col = NULL), alpha = 0.1) +
    facet_wrap(~ xtype, scales = "free_x", nrow = 2, strip.position = "bottom") +
    ## coord_flip(xlim = c(-2.5, 1)) +
    coord_flip(xlim = c(16, -10)## , ylim = 5e-2 * c(-1, 1)
               ) +
    ## scale_x_reverse() +
    labs(x = "$k - k_\\text{pbl}$", y = "") + ## , y = "kg kg$^{-1}$ d$^{-1}$") +
    geom_hline(yintercept = 0, lty = "dashed", col = "grey") +
    geom_vline(xintercept = 0, lty = "dashed", col = "grey") +
    theme(axis.title.x = element_blank(),
          strip.placement = "outside", strip.background = element_blank()) +
    theme(legend.position = c(0.475, 0.125), legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(0, "lines")) +
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"))

## @knitr e3sm-load-cosp
nc <- nc_open("../EAMv1_entrain_NEP_CDNC-LWP.nc")
df.cosp <- nc.to.df(nc, grep("CLWMODIS|CLIMODIS|IWPMODIS|LWPMODIS|TAUWMODIS|OMEGA|TH7001000", names(nc$var), value = TRUE))
df.cosp %<>% rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .)) 
df.cosp %<>% mutate(is.scu = OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55)
df.cosp %<>% mutate(re = REFFCLWMODIS / CLWMODIS,
                    tau = TAUWMODIS / CLWMODIS,
                    lwp_ic = LWPMODIS / CLWMODIS,
                    cdnc_ic = sqrt(10) / (4 * pi * sqrt(1e3)) * sqrt(2e-6) * sqrt(tau) * re ^ (-5/2))
df.cosp %<>% filter(CLIMODIS == 0)
nc_close(nc)

nc <- nc_open("../EAMv1_entrain_NEP_CDNC-LWP_cloud.nc")
df.nep <- nc.to.df(nc, names(nc$var)[!grepl("bnds|lat|lon", names(nc$var))]) %>%
    rename_with(. %>% gsub("_210e_to_245e_15n_to_35n", "", .)) %>%
    mutate(is.scu = OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55)
nc_close(nc)

df.cosp.nep <- left_join(df.nep %>%
                         filter_warm() %>%
                         calculate_incloud(),
                         df.cosp,
                         by = c("time", "ncol"),
                         suffix =  c(".native", ".cosp"))


## @knitr cosp-tau-re-plot
## calculate r_e and tau isolines in log N-log L space
df.tau.isolines <- data.frame(tau = c(1,10,100)) %>%
    ## 
    mutate(slope = -2/5,
           intercept = (6/5) * log10(5 * 1e3 * tau / 9) + (1/5) * log10(9 * 2e-6 / (8 * pi^2 * 1e6)),
           label = sprintf("$\\tau=%g$", tau)) %>%
    mutate(y = log10(adiabatic_lwp_from_tau_re(tau, re = 30 * sqrt(3) * 1e-6)),
           x = log10(adiabatic_nd_from_tau_re(tau, re = 30 * sqrt(3) * 1e-6)))
df.re.isolines <- data.frame(re = c(3,10,30) * 1e-6) %>%
    ## L = (re^6 * (4 pi rho_w)^2) / (18 Gamma_eff) * N_eff^2
    mutate(slope = 2,
           intercept = 6 * log10(re) + 2 * log10(4 * pi * 1e3) - log10(18 * 2e-6),
           label = sprintf("$r_e=%g~\\upmu\\text{m}$", re * 1e6)) %>%
    mutate(y = log10(adiabatic_lwp_from_tau_re(tau = sqrt(0.1), re)),
           x = log10(adiabatic_nd_from_tau_re(tau = sqrt(0.1), re)))
ggplot(df.cosp.nep,
       aes(cdnc_ic.cosp %>% log10, lwp_ic.cosp %>% log10)) +
    stat_bin2d() +
    geom_smooth() +
    scale_fill_distiller(palette = "Spectral", trans = "log10", guide = "none") +
    ## geom_abline(slope = 2, intercept = c(-14, -20)) + ## geom_abline(slope = -1/2, intercept = c(1, 4), lty = "dashed") +
    ## geom_abline(slope = -2/5, intercept = c(1 / 5, 16 / 5)) +
    coord_fixed(xlim = c(5.5, 9.5)) +
    labs(x = "MODIS in-cloud $\\log_{10} (\\nd~\\text{m}^{3})$",
         y = "MODIS in-cloud $\\log_{10} (\\lwp~\\text{kg}^{-1}~\\text{m}^{2})$") +
    geom_abline(aes(slope = slope, intercept = intercept), data = df.tau.isolines, lty = "dashed") +
    geom_label(aes(x = x, y = y, label = label), data = df.tau.isolines) +
    geom_abline(aes(slope = slope, intercept = intercept), data = df.re.isolines, lty = "dashed") +
    geom_label(aes(x = x, y = y, label = label), data = df.re.isolines)


## @knitr e3sm-precip-setup
df.e3sm.precip.bins <- df.e3sm.precip %>%
    left_join(readRDS("../data/perlmutter/entrain_F2010A_run-EAMv1_nudged_2010_fscu.rds") %>%
              left_join(readRDS("../ne30pg2.rds") %>%
                        rename_with(tolower)) %>%
              filter(landfrac == 0) %>%
              group_by(lat, lon) %>%
              summarize(fscu.annual = fscu.annual[1]) %>%
              ungroup()) %>%
    filter_sc() %>%
    filter(precc == 0) %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    mutate(cloudsat = factor(ifelse(precip < 1e-2,
                                    "$R < 10^{-2}$~mm~d$^{-1}$",
                                    "$R\\geq 10^{-2}$~mm~d$^{-1}$"),
                             levels = c("$R < 10^{-2}$~mm~d$^{-1}$",
                                        "$R\\geq 10^{-2}$~mm~d$^{-1}$"))) %>%
    discretize(precip, 6, equal_contents = TRUE, as_factor = TRUE) 

## @knitr precip-lwp-cdnc
  ggplot.precip.bins <- df.e3sm.precip.bins %>%
      mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip))))) %>%
      group_by(precip) %>%
      lwp.cdnc.conditional(cdnc.bins =  exp(seq(log(3e6), log(3e8), length.out = 20))) %>%
      replace_with_conventional_units() %>%
      lwp.cdnc.plot(marginal = TRUE,
                    col = precip, lty = "Mediated by $R$",
                    expand = FALSE, cdnc.density.limits = c(0, 0.25))

  ggplot.precip.bins[[3]] <- ggplot.precip.bins[[3]] +
      geom_point(aes(col = precip),
                 data = df.e3sm.precip.bins %>%
                     group_by(precip) %>%
                     summarize(cdnc_ic = exp(mean(log(cdnc_ic), na.rm = TRUE)),
                               lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
                     replace_with_conventional_units() %>%
                     ungroup() %>%
                     mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip)))))) +
      geom_line(aes(lty = "Mediated by $R$"), col = "black",
                data = df.e3sm.precip.bins %>%
                    group_by(precip) %>%
                    summarize(cdnc_ic = exp(mean(log(cdnc_ic), na.rm = TRUE)),
                              lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
                    replace_with_conventional_units() %>%
                    ungroup() %>%
                    mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip)))))) +
      geom_line(aes(lty = "Model native relationship"),
                col = "black",
                data = df.e3sm.precip.bins %>%
                    discretize_cdnc(bins =  exp(seq(log(3e6), log(3e8), length.out = 20))) %>%
                    mutate(period = "PD") %>%
                    group_by(period) %>%
                    summarize_lwp() %>%
                    replace_with_conventional_units() %>%
                    ungroup()) +
      labs_nd_lwp_conventional() +
      theme(plot.margin = margin(0,0,6,0)) +
      theme(legend.position = c(0.45, 0.8), legend.background = element_rect(fill = NA),
             legend.spacing.y = unit(0, "lines"),
             legend.margin = margin(7, unit = "lines"),
             legend.direction = "vertical", legend.box = "vertical")

  ggplot.precip.bins[[1]] <- ggplot.precip.bins[[1]] + theme(plot.margin = margin(6,0,0,0))
  ggplot.precip.bins[[2]] <- ggplot.precip.bins[[2]] + theme(plot.margin = margin(0,0,0,0))
  ggplot.precip.bins[[4]] <- ggplot.precip.bins[[4]] + theme(plot.margin = margin(0,6,0,0))
  
  ggplot.precip.bins &
      scale_color_brewer("$R$~(mm~d$^{-1}$)", palette = "BrBG", guide = guide_legend(order = 2, ncol = 2)) &
      scale_linetype_discrete("", guide = guide_legend(order = 1)) &
      theme(legend.text = element_text(size = 7),
            legend.title = element_text(size = 8),
            legend.key.size = unit(0.8, "lines"))
  
  ## cowplot::plot_grid(
  ##              df.e3sm.precip.bins %>%
  ##              ## filter(exp == "c1 = 2.4") %>%
  ##              ## discretize_cdnc(bins =  exp(seq(log(1e6), log(3e8), length.out = 20))) %>%
  ##              ## ungroup() %>%
  ##              mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip))))) %>%
  ##              arrange(precip, cdnc_ic) %>%
  ##              ggplot(aes(cdnc_ic, color = precip)) +
  ##              geom_density() +
  ##              ## geom_bar() +
  ##              scale_x_log10() +
  ##              ## facet_grid(. ~ precip) +
  ##              labs(x = "$\\nd$ (in cloud, m$^{-3}$)", ) +
  ##              theme_void() +
  ##              theme(legend.position = "none") +
  ##              theme(
  ##                  strip.background = element_blank(),
  ##                  strip.text.x = element_blank()
  ##              ) +
  ##              coord_cartesian(xlim = c(5e6, 2.5e8), expand = TRUE),
  ##              df.e3sm.precip.bins %>%
  ##                    ## filter(exp == "c1 = 2.4") %>%
  ##     ## discretize(cdnc_ic, 30, equal_contents = TRUE) %>%
  ##     discretize_cdnc(bins =  exp(seq(log(3e6), log(3e8), length.out = 20))) %>%
  ##     group_by(precip) %>%
  ##     summarize_lwp() %>%
  ##     ungroup() %>%
  ##     mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip))))) %>%
  ##     arrange(precip, cdnc_ic) %>%
  ##     ggplot(aes(cdnc_ic, lwp_ic)) +
  ##     geom_line(aes(lty = "Mediated by $R$", col = precip)) +
  ##     geom_point(aes(col = precip),
  ##                data = df.e3sm.precip.bins %>%
  ##                    group_by(precip) %>%
  ##                    summarize(cdnc_ic = exp(mean(log(cdnc_ic), na.rm = TRUE)),
  ##                              lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
  ##                    ungroup() %>%
  ##                    mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip)))))) +
  ##     geom_line(aes(lty = "Mediated by $R$"),
  ##               data = df.e3sm.precip.bins %>%
  ##                    group_by(precip) %>%
  ##                    summarize(cdnc_ic = exp(mean(log(cdnc_ic), na.rm = TRUE)),
  ##                              lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
  ##                    ungroup() %>%
  ##                    mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip)))))) +
  ##     geom_line(aes(lty = "Model native relationship"),
  ##               data = df.e3sm.precip.bins %>%
  ##                   discretize_cdnc(bins =  exp(seq(log(3e6), log(3e8), length.out = 20))) %>%
  ##                   mutate(period = "PD") %>%
  ##                   group_by(period) %>%
  ##                   summarize_lwp() %>%
  ##                   ungroup()) +
  ##     scale_x_log10() +
  ##     scale_y_log10() +
  ##     scale_color_discrete("$R$~(mm~d$^{-1}$)", guide = guide_legend(order = 2, ncol = 1)) +
  ##     scale_linetype_discrete("", guide = guide_legend(order = 1)) +
  ##     ## facet_grid(. ~ precip) +
  ##     labs(x = "$\\nd$ (in cloud, m$^{-3}$)",
  ##          y = "$\\mathcal{L}$ (in cloud, kg~m$^{-2}$)") +
  ##     theme(plot.margin = margin(0,0,12,12)) +
  ##     theme(legend.position = c(0.225, 0.65), legend.background = element_rect(fill = NA),
  ##            legend.spacing.y = unit(0, "lines"),
  ##            legend.margin = margin(6, unit = "lines"),
  ##            legend.direction = "vertical", legend.box = "vertical"
  ##      ) +
  ##     coord_cartesian(xlim = c(5e6, 2.5e8), expand = FALSE)
  ##  , ## ) / (
  ##  nrow = 2, axis = "lr", align = "v", rel_heights = c(0.1, 1)
  ## )  

## @knitr precip-lwp-cdnc-fscu
ggplot.precip.fscu.bins <- df.e3sm.precip.bins %>%
    mutate(precip = factor(precip, labels = (sprintf("%s", levels(precip))))) %>%
    discretize(fscu.annual, seq(0, 1, by = 0.1), as_factor = TRUE) %>%
    group_by(precip, fscu.annual) %>%
    lwp.cdnc.conditional(cdnc.bins =  exp(seq(log(3e6), log(3e8), length.out = 20))) %>%
    replace_with_conventional_units() %>%
    lwp.cdnc.plot(marginal = FALSE,
                  col = precip, lty = "Mediated by $R$",
                  expand = FALSE, cdnc.density.limits = c(0, 0.25)) +
    labs_nd_lwp_conventional() +
    facet_wrap(~ fscu.annual, ncol = 3) +
    theme(plot.margin = margin(0,0,6,0)) +
    theme(legend.position = c(0.45, 0.8), legend.background = element_rect(fill = NA),
          legend.spacing.y = unit(0, "lines"),
          legend.margin = margin(7, unit = "lines"),
          legend.direction = "vertical", legend.box = "vertical")

ggplot.precip.fscu.bins &
    scale_color_brewer("$R$~(mm~d$^{-1}$)", palette = "BrBG", guide = guide_legend(order = 2, ncol = 2)) &
    scale_linetype_discrete("", guide = guide_legend(order = 1)) &
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"))


## @knitr e3sm-precip-multilin
lm.e3sm.precip <- df.e3sm.precip %>%
    filter_sc() %>%
    filter(precc == 0) %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    filter(precip > 0) %>% ## for logarithmic precip
    filter(cdnc_ic > 20e6) %>%
    lm(log10(lwp_ic) ~ log10(cdnc_ic) + log10(precip), ., na.action = na.omit)

lm.e3sm.precip.nd <- df.e3sm.precip %>%
    filter_sc() %>%
    filter(precc == 0) %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    filter(precip > 0) %>%  
    filter(cdnc_ic > 20e6) %>%
    lm(log(lwp_ic) ~ log(cdnc_ic), ., na.action = na.omit)

lm.e3sm.precip.r <- df.e3sm.precip %>%
    filter_sc() %>%
    filter(precc == 0) %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    filter(precip > 0) %>%  
    filter(cdnc_ic > 20e6) %>%
    lm(log(lwp_ic) ~ log(precip), ., na.action = na.omit)

lm.e3sm.precip.binned <- df.e3sm.precip %>%
    filter_sc() %>%
    filter(precc == 0) %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    discretize(precip, 6, equal_contents = TRUE, as_factor = TRUE) %>%
    lm(log10(lwp_ic) ~ log10(cdnc_ic):precip, ., na.action = na.omit)

## @knitr e3sm-precip-exps-setup
## same as e3sm-precip-setup but for multiple experiments
df.e3sm.precip.exps.bins <- df.e3sm.precip.exps %>%
    filter_sc() %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    mutate(cloudsat = factor(ifelse(precip < 1e-2,
                                    "$R < 10^{-2}$~mm~d$^{-1}$",
                                    "$R\\geq 10^{-2}$~mm~d$^{-1}$"),
                             levels = c("$R < 10^{-2}$~mm~d$^{-1}$",
                                        "$R\\geq 10^{-2}$~mm~d$^{-1}$"))) %>%
    discretize(precip, 6, equal_contents = TRUE, as_factor = TRUE) 

## @knitr e3sm-ccn-met-setup
df.ccn <- readRDS("../profiles_ccn.rds") %>%
    mutate(is.scu = OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    group_by(ncol, lon, lat) %>%
    mutate(fscu = mean(is.scu)) %>%
    ungroup() 

df.ccn.met <- df.ccn %>%
    group_by(time) %>%
    mutate(jk.max = jk.pbl[fscu == max(fscu)]) %>%
    ungroup() %>%
    filter(jk.max %in% 10:15) %>%
    pivot_longer(c(T.sfc : V.700, NUMLIQ.top, PS)) %>%
    group_by(lon, lat, jk.max, name) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    pivot_wider() %>%
    mutate(spd.sfc = sqrt(U.sfc^2 + V.sfc^2),
           spd.700 = sqrt(U.700^2 + V.700^2),
           dir.sfc = atan2(V.sfc, U.sfc),
           dir.700 = atan2(V.700, U.700))

df.ccn.max.fscu <- df.ccn %>%
    filter(fscu == max(fscu)) %>%
    slice(1)

## @knitr e3sm-ccn-met-plot
df.ccn.met %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower),
              by = c("lon", "lat")) %>%
    filter(landfrac == 0) %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
    left_join(df.pblh.stats.v1 %>%
              mutate(across(where(is.double), round)),
              by = c("jk.max"="jk.pbl")) %>%
    mutate(pblh = sprintf("{\\small $h=%d$ levels ($%d\\pm%d$ m)}", jk.max, mean, sd)) %>%
    ggplot(aes(lon, lat)) +
    geom_point(aes(col = CCN.sfc), size = 1.25) +
    geom_world_polygon() +
    geom_text(data = df.ccn.max.fscu %>% mutate(lon = ifelse(lon > 180, lon - 360, lon)), col = "black", label = "$\\star$") +
    metR::geom_arrow(aes(mag = spd.700, angle = 180 * dir.700 / pi), 
                     data = . %>% filter(ncol %% 4 == 0),
                     size = 0.4, col = "black", alpha = 0.25) +
    scale_x_geo() +
    scale_y_geo() +
    metR::scale_mag("$|v|_{700\\text{~hPa}}$~(m~s$^{-1}$)", max_size = 0.75, max = 15) +
    coord_fixed(1, c(-150, -115), c(14, 36)) +
    scale_color_distiller("CCN (cm$^{-3}$)", palette = "YlOrBr", direction = 1) +
    facet_wrap(~ pblh, ncol = 3) + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank()) +
    theme(legend.position = "bottom", legend.box = "horizontal") 

## @knitr e3sm-pblh-plot
df.e3sm.pblh <- df.profiles.jk.pbl %>%
    filter(ilev == 72, jk.pbl %in% 10:15) %>%
    filter(precc == 0) %>%
    filter_warm() %>%
    filter_ocean() %>%
    filter(lcc > 0.9) %>%
    calculate_incloud()

ggplot.pblh.bins <- df.e3sm.pblh %>%
    mutate(type = "") %>%
    left_join(df.pblh.stats.v1 %>%
              mutate(across(where(is.double), round)),
              by = c("jk.pbl")) %>%
    mutate(pblh = sprintf("%d ($%d\\pm%d$ m)", jk.pbl, mean, sd)) %>%
    group_by(pblh, type) %>%
    lwp.cdnc.conditional(cdnc.bins = exp(seq(log(3e6), log(3e8), length.out = 50)),
                         lwp.bins = exp(seq(log(0.95e-2), log(0.31), length.out = 25))) %>%
    replace_with_conventional_units() %>%
    lwp.cdnc.plot(marginal = TRUE,
                  col = pblh, ##factor(jk.pbl, levels = 10:15), ## group = jk.pbl,
                  lty = type,
                  expand = FALSE, cdnc.density.limits = c(0, 0.55),
                  legend.position = c(0.475, 0.725)) 

df.e3sm.pblh.summary <- df.e3sm.pblh %>%
    mutate(lwp_ic = ifelse(lcc > 0.9, lwp / lcc, NA),
           cdnc_ic = ifelse(lcc > 0.9, cdnc / lcc, NA)) %>%
    left_join(df.pblh.stats.v1 %>%
              mutate(across(where(is.double), round)),
              by = c("jk.pbl")) %>%
    mutate(pblh = sprintf("%d ($%d\\pm%d$ m)", jk.pbl, mean, sd)) %>%
    group_by(jk.pbl, pblh) %>%
    summarize(lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE)),
              cdnc_ic = exp(mean(log(cdnc_ic), na.rm = TRUE))) %>%
    replace_with_conventional_units() %>%
    ungroup() %>%
    mutate(type = "Mediated by PBL thickness")

ggplot.pblh.bins[[3]] <- ggplot.pblh.bins[[3]] +
    geom_path(data = df.e3sm.pblh.summary, col = "black") +
    geom_path(data = df.e3sm.pblh.summary, aes(group = jk.pbl %/% 2), col = "black", arrow = grid::arrow(type = "closed", angle = 15, length = unit(0.05, "inches")), show.legend = FALSE) +
    geom_path(data = df.e3sm.pblh.summary, aes(group = (jk.pbl + 1) %/% 2), col = "black", arrow = grid::arrow(type = "closed", angle = 15, length = unit(0.05, "inches")), show.legend = FALSE) +
    geom_point(data = df.e3sm.pblh.summary, aes(col = pblh)) +  ## factor(jk.pbl, levels = 10:15))) +
    geom_line(data = df.profiles.fscu %>%
                  filter(ilev == 72) %>%
                  filter(lcc > 0.9, ttop > 273.15) %>%
                  ## filter(is.scu) %>%
                  filter(icc < 1e-3) %>%
                  filter(LANDFRAC < 1e-3) %>%
                  mutate(lwp_ic = lwp / lcc,
                         cdnc_ic = cdnc / lcc) %>%
                  discretize(cdnc_ic, exp(seq(log(7e6), log(120e6), length.out = 25))) %>%
                  group_by(cdnc_ic) %>%
                  summarize(lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
                  replace_with_conventional_units() %>%
                  ungroup() %>%
                  mutate(type = "Model native relationship"),
              col = "black") +
    scale_linetype_manual("", values = c("blank", "solid", "dotted")) +
    labs_nd_lwp_conventional() +
    coord_cartesian(c(9, 1.1e2),
                    c(0.04e3, 0.18e3),
                    expand = FALSE) +
    theme(plot.margin = margin(0,0,6,0)) 

ggplot.pblh.bins[[1]] <- ggplot.pblh.bins[[1]] +
    coord_cartesian(c(9, 1.1e2), c(0, 0.15), expand = FALSE) +
    theme(plot.margin = margin(6,0,0,0))
ggplot.pblh.bins[[2]] <- ggplot.pblh.bins[[2]] + theme(plot.margin = margin(0,0,0,0))
ggplot.pblh.bins[[4]] <- ggplot.pblh.bins[[4]] +
    coord_flip(c(0.04e3, 0.18e3), expand = FALSE) +
    theme(plot.margin = margin(0,6,0,0))

ggplot.pblh.bins &
    scale_color_brewer("PBL depth (model levels)", palette = "Spectral", ## drop = FALSE,
                       guide = guide_legend(ncol = 2)) &
    theme(legend.direction = "vertical", legend.box = "vertical",
          legend.margin = margin(4, unit = "lines"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.key.size = unit(0.8, "lines"))

    ## geom_line(aes(cdnc_ic, lwp_ic, col = factor(jk.pbl, levels = 10:15)),
    ##                          df.profiles.jk.pbl %>%
    ##                          filter(ilev == 72, jk.pbl %in% 10:15) %>%
    ##                          filter(precc == 0) %>%
    ##                          mutate(lwp_ic = ifelse(lcc > 0.9, lwp / lcc, NA),
    ##                                 cdnc_ic = ifelse(lcc > 0.9, cdnc / lcc, NA)) %>%
    ##                          discretize_cdnc(bins =  exp(seq(log(3e6), log(3e8), length.out = 20))) %>%
    ##                          group_by(jk.pbl) %>%
    ##                          summarize_lwp() %>%
    ##                          ungroup() %>%
    ##                          mutate(type = "Mediated by PBL thickness")) +
                   ## geom_line(data = df.global.ttop %>% 
                   ##               filter(OMEGA500 * 86400 > 10e2 & TH7001000 > 18.55) %>%
                   ##               filter(lcc > 0.9, ttop > 273.15) %>%
                   ##               filter(icc < 1e-3) %>%
                   ##               filter(LANDFRAC < 1e-3) %>%
                   ##               mutate(lwp_ic = lwp / lcc,
                   ##                      cdnc_ic = cdnc / lcc) %>%
                   ##               discretize(cdnc_ic, exp(seq(log(1e6), log(300e6), length.out = 37 * 2))) %>%
                   ##               group_by(cdnc_ic) %>%
                   ##               dplyr::summarize(lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
                   ##               ungroup() %>%
                   ##               mutate(type = "Full model")) +
                   ## geom_line(data = df.profiles.fscu %>%
                   ##               filter(ilev == 72) %>%
                   ##               filter(lcc > 0.9, ttop > 273.15) %>%
                   ##               ## filter(is.scu) %>%
                   ##               filter(icc < 1e-3) %>%
                   ##               filter(LANDFRAC < 1e-3) %>%
                   ##               mutate(lwp_ic = lwp / lcc,
                   ##                      cdnc_ic = cdnc / lcc) %>%
                   ##               discretize(cdnc_ic, exp(seq(log(7e6), log(120e6), length.out = 25))) %>%
                   ##               group_by(cdnc_ic) %>%
                   ##               summarize(lwp_ic = exp(mean(log(lwp_ic), na.rm = TRUE))) %>%
                   ##               ungroup() %>%
                   ##               mutate(type = "Model native relationship")) +
                   ## scale_x_log10(labels = NULL) +
                   ## scale_y_log10() +
                   ## labs(x = "$N_d$ (in cloud, m$^{-3}$)",
                   ##      y = "$\\mathcal{L}$ (in cloud, kg~m$^{-2}$)") +
                   ## coord_cartesian(xlim = c(1e7, 1e8), ylim = c(6.5e-2, 11e-2)) +
                   ## theme(axis.title.x = element_blank()) +
                   ## theme(legend.position = c(0.3, 0.275), legend.background = element_rect(fill = NA),
                   ##       legend.spacing.y = unit(0, "lines")),
                   ## df.profiles.jk.pbl %>%
                   ## filter(ilev == 72, jk.pbl %in% 10:15) %>%
                   ## filter(lcc > 0.9, ttop > 273.15) %>%
                   ## ## filter(is.scu) %>%
                   ## filter(icc < 1e-3) %>%
                   ## ## filter(LANDFRAC < 1e-3) %>%
                   ## mutate(lwp_ic = lwp / lcc,
                   ##        cdnc_ic = cdnc / lcc) %>%
                   ## ggplot(aes(cdnc_ic, col = factor(jk.pbl, levels = 10:15))) +
                   ## geom_density() +
                   ## scale_x_log10() +
                   ## scale_y_continuous(labels = tikz_sanitize_sparse) +
                   ## scale_color_brewer("PBL depth (model levels)", palette = "Spectral", drop = FALSE,
                   ##                    guide = "none") +
                   ## labs(x = "$N_d$ (in cloud, m$^{-3}$)") +
                   ## coord_cartesian(xlim = c(1e7, 1e8)),
                   ## ncol = 1, align = "v", rel_heights = c(3, 1))

## @knitr e3sm-pblh-multilin
lm.e3sm.pblh <- df.e3sm.pblh %>%
    filter(precc == 0) %>%
    left_join(df.pblh.stats.v1 %>%
              mutate(across(where(is.double), round)),
              by = c("jk.pbl")) %>%
    mutate(pblh = mean) %>%
    filter(cdnc_ic > 20e6) %>%
    lm(log(lwp_ic) ~ log(cdnc_ic) + log(pblh), ., na.action = na.omit)

lm.e3sm.pblh.nd <- df.e3sm.pblh %>%
    filter(precc == 0) %>%
    left_join(df.pblh.stats.v1 %>%
              mutate(across(where(is.double), round)),
              by = c("jk.pbl")) %>%
    mutate(pblh = mean) %>%
    filter(cdnc_ic > 20e6) %>%
    lm(log(lwp_ic) ~ log(cdnc_ic), ., na.action = na.omit)

lm.e3sm.pblh.binned <- df.e3sm.pblh %>%
    filter(precc == 0) %>%
    filter(cdnc_ic > 20e6) %>%
    lm(log(lwp_ic) ~ log(cdnc_ic) : factor(jk.pbl), ., na.action = na.omit)

## @knitr e3sm-pblh-stats
df.pblh.summary.v1 <- df.entrain.aer.v1 %>%
    mutate(jk.pbl = 73 - jk.pbl) %>%
    plyr::ddply(~ jk.pbl, . %$% summary(h)) 

df.pblh.stats.v1 <- df.entrain.aer.v1 %>%
    mutate(jk.pbl = 73 - jk.pbl) %>%
    group_by(jk.pbl) %>%
    summarize(mean = mean(h),
              sd = sd(h)) %>%
    ungroup()

## @knitr multimodel-delta-pred-actual-setup
df.multimodel.conditional <- df.multimodel %>%
    filter(model %in% c("E3SM", "GISS")) %>%
    group_by(model, period) %>%
    lwp.cdnc.conditional(lwp.bins = exp(seq(log(0.95e-2), log(0.31), length.out = 25 * 2))) %>%
    ungroup()

weighted.mean.fun <- weighted.mean
df.pred <- ddply(df.multimodel.conditional, ~ model, function(df) {
    predict_lwp_cdnc_conditional(df %>% filter(period == "PD"),
                                 df %>% filter(period == "PI"))
}) %>%
    group_by(model) %>%
    summarize(lwp_ic = weighted.mean.fun(lwp_ic, p, na.rm = TRUE),
              period = "PI_pred")

df.multimodel.pred.actual <- df.multimodel.conditional %>%
    group_by(model, period) %>%
    summarize(lwp_ic = weighted.mean.fun(lwp_ic, count.lwp, na.rm = TRUE),
              cdnc_ic = weighted.mean.fun(cdnc_ic, count.cdnc, na.rm = TRUE)) %>%
    bind_rows(df.pred) %>%
    group_by(model) %>%
    mutate(cdnc_ic = replace(cdnc_ic, is.na(cdnc_ic), cdnc_ic[period == "PI"]))
  
## @knitr multimodel-delta-pred-actual-plot
df.multimodel.pred.actual %>%
    replace_with_conventional_units() %>%
    pivot_wider(names_from = period, values_from = c(lwp_ic, cdnc_ic)) %>%
    ggplot(aes(col = model)) +
    geom_segment(aes(x = cdnc_ic_PD, y = lwp_ic_PD, xend = cdnc_ic_PI, yend = lwp_ic_PI, lty = "PI actual"), arrow = grid::arrow(type = "closed", angle = 15, length = unit(0.15, "inches")), show.legend = FALSE) +
    geom_segment(aes(x = cdnc_ic_PD, y = lwp_ic_PD, xend = cdnc_ic_PI_pred, yend = lwp_ic_PI_pred, lty = "PI predicted"), arrow = grid::arrow(type = "closed", angle = 15, length = unit(0.15, "inches")), show.legend = FALSE) +
    geom_segment(aes(x = cdnc_ic_PD, y = lwp_ic_PD, xend = cdnc_ic_PI, yend = lwp_ic_PI, lty = "PI actual")) + 
    geom_segment(aes(x = cdnc_ic_PD, y = lwp_ic_PD, xend = cdnc_ic_PI_pred, yend = lwp_ic_PI_pred, lty = "PI predicted")) + 
    geom_segment(aes(x = cdnc_ic_PD, y = lwp_ic_PD, xend = cdnc_ic_PD, yend = lwp_ic_PD, lty = "PD")) + 
    geom_point(aes(cdnc_ic_PI, lwp_ic_PI, shape = "PI actual")) + 
    geom_point(aes(cdnc_ic_PI_pred, lwp_ic_PI_pred, shape = "PI predicted")) + 
    geom_point(aes(cdnc_ic_PD, lwp_ic_PD, shape = "PD")) +
    scale_x_log10() +
    scale_y_log10() +
    scale_shape("") +
    scale_color_uscms_selected("") +
    scale_linetype_manual("", values = c("blank", "solid", "dotted")) +
    labs_nd_lwp_conventional() +
    theme(legend.position = c(0.75, 0.75), legend.background = element_rect(fill = NA))

## @knitr multimodel-lon-lat-susc-setup
doParallel::registerDoParallel(8)
df.lon.lat.slopes <- df.multimodel %>%
    mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
    mutate(lon = 5 * (lon %/% 5),
           lat = 5 * (lat %/% 5)) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              mutate(lon = ifelse(lon > 180, lon - 360, lon)) %>%
              mutate(lon = 5 * (lon %/% 5),
                     lat = 5 * (lat %/% 5)) %>%
              group_by(lon, lat) %>%
              summarize(landfrac = mean(LANDFRAC))) %>%
    filter(landfrac == 0) %>%
    plyr::ddply(~ model + period + lon + lat, function(df) {
        ret <- data.frame(n = nrow(df))
        lm.fit <- try(lm(log(lwp_ic) ~ log(cdnc_ic), df))
        if (!any(class(lm.fit) == "try-error")) {
            ret %<>% mutate(lwp_susc = lm.fit$coefficients[2])
            sum.lm.fit <- try(summary(lm.fit)$coefficients[2,2])
            if (!any(class(sum.lm.fit) == "try-error")) 
                ret %<>% mutate(lwp_susc.se = sum.lm.fit)
        }
        ret
    }, .parallel = TRUE) 

## @knitr multimodel-lon-lat-susc-plot
df.lon.lat.slopes %>%
    filter(period == "PD") %>%
    filter(between(lwp_susc, -1, 1)) %>%
    ## filter(between(lwp_susc.se, -0.1, 0.1)) %>%
    ggplot(aes(lon, lat, fill = lwp_susc)) +
    geom_raster() + 
    geom_world_polygon(highres = FALSE) +
    scale_x_geo() +
    scale_y_geo() +
    coord_fixed() +
    facet_wrap( ~ model) +
    scale_fill_distiller("$d\\log\\lwp/d\\log\\nd$", palette = "RdBu", labels = tikz_sanitize) +
    guides(fill = guide_colorbar(direction = "horizontal", title.vjust = 0.75, barwidth = unit(5, "cm"))) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
