## @knitr e3sm-scm-aer-ind3-setup
doParallel::registerDoParallel(8)
df.e3sm.scm.aer.ind3 <- expand.grid(aer = 2^(-3:3), sed = c(0, 1, 2, 4), c1 = c(2.4)) %>%
    ddply(~ aer + sed + c1, function(df) {
        res <- try(df %$% load_scm(sprintf("../scm_aer/e3sm_scm_DYCOMSrf02_aer%g_sed%g_c1%g_ind3.eam.h0.1999-07-11-00000.nc", aer, sed, c1)) %>%
                   scm_budget())
        if (any(class(res) == "try-error")) {
            NULL 
        } else {
            res
        }
    }, .parallel = TRUE) %>%
    mutate(cdnc_ic_ind3 = ifelse(lcc < 0.75, NA, cdnc / lcc),
           lwp_ic_ind3 = ifelse(lcc < 0.75, NA, lwp / lcc)) %>%
    mutate(numliq_ic = ifelse(f.minus < 0.75, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.75, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R))

## @knitr e3sm-scm-timeseries
df.e3sm.scm.aer.ind3 %>%
    filter(time < 1) %>%
    filter(sed < 4) %>%
    ## filter(f.minus > 0.1) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    gather(var, val, lcc, f.minus, lwp_ic_ind3, cdnc_ic_ind3, lwp_ic, cdnc_ic) %>%
    ## mutate(var = revalue(var, tikz_replacements_unitful())) %>%
    mutate(var = gsub("_", "\\\\_", var)) %>%
    ggplot(aes(time, val, ## col = factor(aer),
               lty = factor(sed))) +
    scale_color_brewer(palette = "Spectral", direction = -1) +
    geom_line() +
    facet_grid(var ~ aer, scales = "free_y")

## @knitr e3sm-scm-aer-scans-setup
doParallel::registerDoParallel(8)

df.aer.all <- expand.grid(aer = 2^(-3:3), sed = c(0, 1, 2, 4), c1 = c(1.335, 2.4, 3.6)) %>%
    ddply(~ aer + sed + c1, function(df) {
        res <- try(df %$% load_scm(sprintf("../scm_aer/e3sm_scm_DYCOMSrf02_aer%g_sed%g_c1%g.eam.h0.1999-07-11-00000.nc", aer, sed, c1)) %>%
                   rename(ccn4 = ccn3) %>%
                   scm_budget())
        if (class(res) != "try-error") {
            res
        } else {
            NULL
        }
    }, .parallel = TRUE) %>%
    mutate(numliq_ic = ifelse(f.minus < 0.75, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.75, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R))

df.aer.sed <- df.aer.all %>% filter(c1 == 2.4) ## v2 default c1

## @knitr e3sm-scm-aer-swoff-setup
df.aer.swoff <- expand.grid(aer = 2^(-3:3), sed = c(0, 1, 2, 4), c1 = c(2.4)) %>%
    ddply(~ aer + sed + c1, function(df) {
        res <- try(df %$% load_scm(sprintf("../scm_aer/e3sm_scm_DYCOMSrf02_aer%g_sed%g_c1%g_swoff.eam.h0.1999-07-11-00000.nc", aer, sed, c1)) %>%
                   rename(ccn4 = ccn3) %>%
                   scm_budget())
        if (class(res) != "try-error") {
            res
        } else {
            NULL
        }
    }, .parallel = TRUE) %>%
    mutate(numliq_ic = ifelse(f.minus < 0.75, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.75, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R))

## @knitr e3sm-scm-aer-sed-timeseries-plot
df.aer.swoff %>%
    mutate(numliq_ic = ifelse(f.minus < 0.75, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.75, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R),
           precip = delta.R) %>%
    filter(time < 1) %>%
    filter(sed == 1) %>%
    ## filter(f.minus > 0.1) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    mutate(h = h - 6254.913 / 9.8) %>% ## kluge due to nonzero surface geopotential in SCM runs
    replace_with_conventional_units() %>%
    gather(var, val,
           h, f.minus, cdnc_ic,
           E_q, # E_theta, E_q.delta.q, E_theta.delta.theta,
           ## q_t.plus, q_t.0, theta_l.plus, theta_l.0
           ## delta.theta, delta.q
           lwp_ic, precip
        ) %>%
    mutate(var = revalue(var, c(lwp_ic = "lwp_ic.conventional",
                                cdnc_ic = "cdnc_ic.conventional"))) %>%
    mutate(var = revalue(var, tikz_replacements_unitful()),
           var = factor(var, levels = unique(var))) %>%
    ggplot(aes(time, val, col = factor(aer))) +
    scale_color_brewer("$\\na \\times$", palette = "Spectral", direction = -1) +
    geom_line() +
    facet_grid(var ~ ., scales = "free_y", switch = "y") +
    ggh4x::facetted_pos_scales(y = list(scale_y_continuous(), scale_y_continuous(), scale_y_log10(),
                                        scale_y_continuous(), scale_y_continuous(), scale_y_log10())) +
    labs(x = "Time (d)") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside", strip.background.y = element_blank())

## @knitr giss-scm-na-setup
doParallel::registerDoParallel(8)

df.giss.na <- expand.grid(aer = 10 * 2^(1:4)) %>%
    ddply(~ aer, function(df) {
        res <- try(df %$% load_scm_giss(sprintf("../ModelE3_2x2p5deg_ACI/SCM_old/baseline/allsteps.allmergeSCM_RF02_jmu_na%d.nc",
                                                aer)) %>%
                   mutate(ccn0p2 = NA,
                          precl = prec,
                          precc = 0) %>%
                   rename_fields_giss() %>%
                   scm_budget(inv.type = "jk.qt"))
        if ("try-error" %in% class(res)) {
            NULL
        } else {
            res
        }
    }, .parallel = TRUE) %>%
    mutate(log10.delta.R = log10(delta.R))

## @knitr giss-scm-ppe-setup
doParallel::registerDoParallel(8)

df.giss.ppe <- expand.grid(aer = c(20, 30, 40, 60, 80, 120, 160, 320), experiment = c("t705ml", "t705ml_e25s0")) %>%
    ddply(~ aer + experiment, function(df) {
        res <- try(df %$% load_scm_giss(sprintf("../ModelE3_2x2p5deg_ACI/SCM_RF02/%s/allsteps.allmergeSCM_RF02_%s.nc",
                                                gsub("^t705ml", "t705ml_dp10", experiment),
                                                gsub("^t705ml", sprintf("t705ml_na%d_dp10", aer), experiment))) %>%
                   mutate(ccn0p2 = NA,
                          precl = prec,
                          precc = 0) %>%
                   rename_fields_giss() %>%
                   scm_budget(inv.type = "jk.qt"))
        if ("try-error" %in% class(res)) {
            NULL
        } else {
            res
        }
    }, .parallel = TRUE) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    mutate(experiment = gsub("_", "\\\\_", experiment)) 

## @knitr giss-scm-def-setup
doParallel::registerDoParallel(8)

df.giss.ppe <- ## df.giss.scm.def <-
    expand.grid(aer = c(20, 30, 40, 60, 80, 120, 160, 320, 640), experiment = c("def", "t705ml")) %>%
    ddply(~ aer + experiment, function(df) {
        res <- try(df %$% load_scm_giss(sprintf("../data/giss/scm_rf02_nudge0_L62scm/allsteps.allmergerf02_%s_%s.nc",
                                                experiment,
                                                sprintf("na%d", aer))) %>%
                   mutate(ccn0p2 = NA,
                          precl = prec,
                          precc = 0) %>%
                   rename_fields_giss() %>%
                   scm_budget(inv.type = "jk.qt", gradient.type = "p", lev.type = "half"))
        if ("try-error" %in% class(res)) {
            NULL
        } else {
            res
        }
    }, .parallel = TRUE) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    mutate(experiment = gsub("_", "\\\\_", experiment)) %>%
    mutate(time = time - 0.25/24)  ## Andy thinks we should use the middle of
                                   ## the time step rather than the end

## @knitr dharma-setup
nc.dharma <- nc_open("../data/giss/Muelmenstadt/dharma_rf02_na190_newgrid_scalars.nc")
df.dharma <- nc.to.df(nc.dharma, c("zi", "LWP", "cc", "Nc", "prec_sfc")) %>%
        rename_fields_dharma() 
nc_close(nc.dharma)

## @knitr dharma-aer-setup
df.dharma.aer <- 
    expand.grid(aer = c(95, 190, 380, 760, 1520, 3040)) %>%
    ddply(~ aer, function(df) {
        on.exit(nc_close(nc.dharma))
        nc.dharma <- nc_open(sprintf("../data/giss/Muelmenstadt/dharma_runs/dharma_rf02_na%d_scalars.nc",
                                     df$aer))
        df.dharma <- nc.to.df(nc.dharma, c("zi", "LWP", "cc", "Nc", "prec_sfc")) %>%
            rename_fields_dharma() 
    }, .parallel = TRUE) 

## @knitr giss-dharma-setup
df.giss.dharma <- ## df.giss.scm.def <-
    expand.grid(experiment = c("na60", "na60_gissLW")) %>%
    ddply(~ experiment, function(df) {
        res <- try(df %$% load_scm_giss(sprintf("../data/giss/Muelmenstadt/allsteps.allmergerf02_def_%s.nc",
                                                experiment)) %>%
                   mutate(ccn0p2 = NA,
                          omega = NA,
                          precl = prec,
                          precc = 0) %>%
                   mutate(ssct_ncl = weighted.mean(nclic, qcl)) %>%
                   mutate(cLWPss = lwp + cLWPmc + pLWPmc) %>%
                   rename_fields_giss() %>%
                   scm_budget(inv.type = "jk.qt", gradient.type = "p", lev.type = "half"))
        if ("try-error" %in% class(res)) {
            NULL
        } else {
            res
        }
    }, .parallel = FALSE) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    mutate(experiment = gsub("_", "\\\\_", experiment)) %>%
    mutate(time = ifelse(time %% (0.5 / 24) == 0,
                         time - 0.25 / 24,
                         time))  ## Andy thinks we should use the middle of
                                   ## the time step rather than the end

## @knitr e3sm-scm-aer-sed-lwpsusc-plot
df.aer.swoff %>%
    mutate(numliq_ic = ifelse(f.minus < 0.9, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.9, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    filter(between(time, 2/24, 1/2)) %>%
    filter(sed == 1) %>%
    mutate(sed = sprintf("$Q_\\text{sed} \\times %d$", sed)) %>%
    group_by(time, sed, c1) %>%
    mutate(delta.L = lwp_ic - lwp_ic[aer == 1],
           delta.log.L = log(lwp_ic / lwp_ic[aer == 1]),
           delta.log.N = log(cdnc_ic / cdnc_ic[aer == 1]),
           S.L = delta.log.L / delta.log.N) %>%
    ungroup() %>%
    ## filter(f.minus > 0.1) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    ggplot(aes(x = factor(""), y = S.L, fill = factor(aer))) +
    scale_fill_brewer("$\\na \\times$", palette = "Spectral", direction = -1, drop = TRUE) +
    ## geom_boxplot() +
    geom_violin(draw_quantiles = 0.5) +
    labs(x = "", y = "$d\\log\\lwp/d\\log\\nd$") +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(legend.position = c(0.12, 0.3))

## @knitr e3sm-scm-aer-sed-esusc-plot
df.aer.swoff %>%
    mutate(numliq_ic = ifelse(f.minus < 0.9, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.9, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    filter(between(time, 2/24, 1/2)) %>%
    filter(sed < 4) %>%
    filter(sed == 1) %>%
    mutate(sed = sprintf("$Q_\\text{sed} \\times %d$", sed)) %>%
    group_by(time, sed, c1) %>%
    mutate(delta.L = lwp_ic - lwp_ic[aer == 1],
           delta.log.L = log(lwp_ic / lwp_ic[aer == 1]),
           delta.log.N = log(cdnc_ic / cdnc_ic[aer == 1]),
           delta.E_q = E_q - E_q[aer == 1],
           delta.E_theta = E_theta - E_theta[aer == 1],
           delta.log.E_q = log(E_q / E_q[aer == 1]),
           delta.log.E_theta = log(E_theta / E_theta[aer == 1]),
           S.L = delta.log.L / delta.log.N,
           S.Eq = delta.log.E_q / delta.log.N,
           S.Etheta = delta.log.E_theta / delta.log.N) %>%
    ungroup() %>%
    mutate(S.E = S.Eq) %>%
    ## filter(f.minus > 0.1) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    ggplot(aes(x = "", y = S.E, fill = factor(aer))) +
    scale_fill_brewer("$\\na \\times$", palette = "Spectral", direction = -1, drop = FALSE) +
    scale_linetype("") +
    ## geom_boxplot() +
    geom_violin(draw_quantiles = 0.5) +
    coord_cartesian(ylim = c(-0.25, 1)) +
    labs(x = "", y = "$d\\log E/d\\log\\nd$") +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    theme(legend.position = "none")

## @knitr aer-sed-1e-present-2
df.aer.swoff %>%
    mutate(numliq_ic = ifelse(f.minus < 0.9, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.9, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    filter(between(time, 2/24, 1/2)) %>%
    filter(sed < 4) %>%
    mutate(sed = sprintf("$Q_\\text{sed} \\times %d$", sed)) %>%
    group_by(time, sed, c1) %>%
    mutate(delta.L = lwp_ic - lwp_ic[aer == 1],
           delta.log.L = log(lwp_ic / lwp_ic[aer == 1]),
           delta.log.N = log(cdnc_ic / cdnc_ic[aer == 1]),
           delta.E_q = E_q - E_q[aer == 1],
           delta.E_theta = E_theta - E_theta[aer == 1],
           delta.log.E_q = log(E_q / E_q[aer == 1]),
           delta.log.E_theta = log(E_theta / E_theta[aer == 1]),
           S.L = delta.log.L / delta.log.N,
           S.Eq = delta.log.E_q / delta.log.N,
           S.Etheta = delta.log.E_theta / delta.log.N) %>%
    ungroup() %>%
    ## filter(f.minus > 0.1) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    ggplot(aes(x = factor(sed), y = S.L)) + ## , fill = factor(aer))) +
    ## scale_fill_brewer("$N_a \\times$", palette = "Spectral", direction = -1, drop = TRUE) +
    ## geom_boxplot() +
    geom_violin(draw_quantiles = 0.5) +
    labs(x = "", y = "$d\\log\\lwp/d\\log\\nd$") +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    theme(legend.position = c(0.75, 0.8))

## @knitr aer-sed-swoff-1e
df.aer.swoff %>%
    mutate(numliq_ic = ifelse(f.minus < 0.75, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = ifelse(f.minus < 0.75, NA, lwp / f.minus)) %>%
    mutate(log10.delta.R = log10(delta.R)) %>%
    filter(between(time, 2/24, 1/2)) %>%
    filter(sed == 1) %>%
    group_by(time, sed, c1) %>%
    mutate(delta.L = lwp_ic - lwp_ic[aer == 1],
           delta.log.L = log(lwp_ic / lwp_ic[aer == 1]),
           delta.log.N = log(cdnc_ic / cdnc_ic[aer == 1]),
           S.L = delta.log.L / delta.log.N) %>%
    ungroup() %>%
    ## filter(f.minus > 0.1) %>%
    mutate(jk.pbl = 72 - jk.pbl) %>%
    ggplot(aes(x = "", y = S.L, fill = factor(aer))) +
    scale_fill_brewer("$\\na \\times$", palette = "Spectral", direction = -1, drop = TRUE) +
    geom_violin(draw_quantiles = 0.5) +
    labs(x = "", y = "$d\\log\\lwp/d\\log\\nd$") +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    theme(legend.position = c(0.8, 0.35))

## @knitr e3sm-scm-rename
df.aer.sed.renamed <- df.aer.swoff %>%
    filter(sed < 4) %>%
    mutate(lcc = f.minus) %>%
    mutate(numliq_ic = ifelse(f.minus < 0.75, NA, numliq.minus / f.minus)) %>%
    mutate(cdnc_ic = numliq_ic * rho.minus) %>%
    mutate(lwp = q_l.0 * rhoH) %>%
    mutate(lwp_ic = lwp / f.minus) %>%
    mutate(experiment = factor(sed))

## @knitr e3sm-scm-means-setup
df.aer.sed.means <- df.aer.sed.renamed %>%
    filter(between(time, 2/24, 12/24)) %>%
    filter(lcc > 0.9) %>%
    group_by(aer, experiment) %>%
    summarize(across(c(lwp_ic, cdnc_ic, E_theta, E_q, E_theta.delta.theta, E_q.delta.q, delta.R),
                     .fns = list(mean = ~ mean(.x, na.rm = TRUE),
                                 sd = . %>% sd(na.rm = TRUE))),
              n = n()) %>%
    ungroup()  %>%
    ## conventional units
    mutate(lwp_ic_mean_conv = lwp_ic_mean * 1000,
           lwp_ic_sd_conv = lwp_ic_sd * 1000) 

## @knitr e3sm-scm-endpoints-setup
df.aer.sed.endpoints <- df.aer.sed.renamed %>%
    filter(between(time, 10/24, 1/2)) %>%
    filter(lcc > 0.9) %>%
    group_by(aer, experiment) %>%
    summarize(across(c(lwp_ic, cdnc_ic, E_theta, E_q, delta.R),
                     .fns = list(mean = ~ mean(.x, na.rm = TRUE),
                                 sd = . %>% sd(na.rm = TRUE))),
              n = n()) %>%
    ungroup()

## @knitr e3sm-scm-tendency-setup
df.aer.sed.tendency <- df.aer.sed.renamed %>%
    filter(between(time, 2/24, 12/24)) %>%
    filter(lcc > 0.9) %>%
    plyr::ddply(~ experiment + aer, function(df) {
        ret <- data.frame(n = nrow(df))
        lm.fit <- try(lm(lwp_ic ~ time, df %>% mutate(time = time * 86400)))
        sum.lm.fit <- try(summary(lm.fit)$coefficients[2,2])
        if (!any(class(lm.fit) == "try-error")) 
            ret %<>% mutate(lwp_ic_tend = lm.fit$coefficients[2])
        if (!any(class(sum.lm.fit) == "try-error")) 
            ret %<>% mutate(lwp_ic_tend.se = sum.lm.fit)
    }) 

## @knitr e3sm-scm-ndeps-setup
df.e3sm.scm.ndeps <-
    bind_rows(
        ## dL/dt vs Nd
        df.aer.sed.tendency %>%
        full_join(df.aer.sed.means, by = c("experiment", "aer")) %>%
        ## filter(experiment == 1) %>%
        transmute(experiment,
                  x = cdnc_ic_mean,
                  y = lwp_ic_tend,
                  ymin = lwp_ic_tend - lwp_ic_tend.se,
                  ymax = lwp_ic_tend + lwp_ic_tend.se,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = "$\\partial\\lwp/\\partial t$~(kg~m$^{-2}$~s$^{-1}$)"),
        ## E vs Nd
        df.aer.sed.tendency %>%
        full_join(df.aer.sed.means) %>%
        ## filter(experiment == 1) %>%
        transmute(experiment,
                  x = cdnc_ic_mean,
                  y = E_theta.delta.theta_mean,
                  ymin = E_theta.delta.theta_mean - E_theta.delta.theta_sd,
                  ymax = E_theta.delta.theta_mean + E_theta.delta.theta_sd,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["E_theta.delta.theta.hatted"]),
        ## E vs Nd
        df.aer.sed.tendency %>%
        full_join(df.aer.sed.means) %>%
        ## filter(experiment == 1) %>%
        transmute(experiment,
                  x = cdnc_ic_mean,
                  y = E_q.delta.q_mean,
                  ymin = E_q.delta.q_mean - E_q.delta.q_sd,
                  ymax = E_q.delta.q_mean + E_q.delta.q_sd,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["E_q.delta.q.hatted"]),
        ## E vs Nd
        df.aer.sed.tendency %>%
        full_join(df.aer.sed.means) %>%
        ## filter(experiment == 1) %>%
        transmute(experiment,
                  x = cdnc_ic_mean,
                  y = E_q_mean,
                  ymin = E_q_mean - E_q_sd,
                  ymax = E_q_mean + E_q_sd,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["E"]),
        ## L vs Nd
        df.aer.sed.tendency %>%
        full_join(df.aer.sed.means) %>%
        ## filter(experiment == 1) %>%
        transmute(experiment,
                  x = cdnc_ic_mean,
                  y = lwp_ic_mean,
                  y = lwp_ic_mean_conv,
                  ymin = lwp_ic_mean_conv - lwp_ic_sd_conv,
                  ymax = lwp_ic_mean_conv + lwp_ic_sd_conv,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["lwp_ic.conventional"])) %>%
    mutate(model = "E3SM")

## @knitr e3sm-scm-ndeps-plot
df.e3sm.scm.ndeps %>%
    mutate(x = x * 1e-6) %>% ## conventional Nd units
    filter(experiment == 1) %>%
    ggplot(aes(x, y)) +
    geom_line(aes(lty = model)) +
    geom_pointrange(aes(shape = model,
                        ymin = ymin,
                        ymax = ymax)) +
    scale_x_log10() +
    ## scale_y_log10() +
    facet_grid(ytype ~ ., scales = "free_y", switch = "y") +
    labs(x = tikz_replacements_unitful()["cdnc_ic.conventional"], y = "") +
    guides(shape = "none", linetype = "none") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside", strip.background.y = element_blank())

## @knitr giss-ppe-timeseries
df.giss.ppe %>% ## mutate(aer = as.character(aer)) %>%
    ## bind_rows(df.giss %>% mutate(aer = "default")) %>%
    mutate(model = "GISS") %>%
    filter(experiment == "def") %>%
    ## filter(1 - lcc < 1e-7) %>%
    rename(precip = delta.R,
           cloud = lcc) %>%
    replace_with_conventional_units() %>%
    gather(var, val, h, cloud, ## delta.theta, delta.q, theta_l.plus, q_t.plus, theta_l.0, q_t.0, delta.F.cp,
           cdnc_ic, ## E_theta,
           E_q, ## E_theta.delta.theta, E_q.delta.q,
           lwp_ic, precip) %>%
    mutate(var = revalue(var, c(lwp_ic = "lwp_ic.conventional",
                                cdnc_ic = "cdnc_ic.conventional"))) %>%
    mutate(var = revalue(var, tikz_replacements_unitful()),
           var = factor(var, levels = unique(var))) %>%
    group_by(var) %>%
    ## mutate(val = as.vector(stats::filter(val, c(0.5, 0.5)))) %>%
    ungroup() %>%
    filter(between(time, 0/24, 23/24)) %>%
    ggplot(aes(time, val, 
               col = factor(aer))) +
    scale_color_brewer("$\\na$~(cm$^{-3}$)", palette = "Spectral", direction = -1) +
    geom_line() +
    facet_grid(var ~ ., scales = "free_y", switch = "y") +
    ggh4x::facetted_pos_scales(y = list(scale_y_continuous(), scale_y_continuous(), scale_y_log10(),
                                        scale_y_continuous(), scale_y_continuous(), scale_y_log10())) +
    labs(x = "Time (d)") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside", strip.background.y = element_blank())

## @knitr dharma-aer-timeseries
df.dharma.aer %>% ## mutate(aer = as.character(aer)) %>%
    mutate(model = "GISS") %>%
    mutate(time = time / 86400) %>%
    rename(precip = prec) %>%
    gather(var, val, h, cloud, ## delta.theta, delta.q, theta_l.plus, q_t.plus, theta_l.0, q_t.0, delta.F.cp,
           cdnc_ic, ## E_theta,
           lwp_ic, precip) %>%
    mutate(var = revalue(var, tikz_replacements_unitful()),
           var = factor(var, levels = unique(var))) %>%
    group_by(var) %>%
    ungroup() %>%
    filter(between(time, 0/24, 23/24)) %>%
    ggplot(aes(time, val, 
               col = factor(aer))) +
    scale_color_brewer("$\\na$~(cm$^{-3}$)", palette = "Spectral", direction = -1) +
    geom_line() +
    facet_grid(var ~ ., scales = "free_y", switch = "y") +
    ggh4x::facetted_pos_scales(y = list(scale_y_continuous(), scale_y_continuous(), scale_y_log10(),
                                        scale_y_continuous(), scale_y_continuous(), scale_y_log10())) +
    labs(x = "Time (d)") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside", strip.background.y = element_blank())

## @knitr giss-dharma-timeseries
bind_rows(df.dharma %>%
          mutate(experiment = "LES") %>%
          mutate(time = time / 86400) %>%
          rename(precip = prec),
          df.giss.dharma %>%
          mutate(experiment = revalue(experiment, c("na60" = "SCM (Beer's Law LW)",
                                                    "na60\\_gissLW" = "SCM (Native LW)"))) %>%
          rename(precip = delta.R,
                 cloud = lcc)) %>%
    mutate(experiment = factor(experiment, unique(experiment))) %>%
    replace_with_conventional_units() %>%
    gather(var, val, h, cloud, ## delta.theta, delta.q, theta_l.plus, q_t.plus, theta_l.0, q_t.0, delta.F.cp,
           cdnc_ic, ## E_theta,
           lwp_ic, precip) %>%
    mutate(var = revalue(var, c(lwp_ic = "lwp_ic.conventional",
                                cdnc_ic = "cdnc_ic.conventional"))) %>%
    mutate(var = revalue(var, tikz_replacements_unitful()),
           var = factor(var, levels = unique(var))) %>%
    group_by(var) %>%
    ungroup() %>%
    filter(between(time, 0/24, 23.75/24)) %>%
    ggplot(aes(time, val, 
               col = experiment, lty = experiment)) +
    scale_color_brewer("", palette = "Dark2", direction = -1) +
    scale_linetype_discrete("") +
    geom_line() +
    facet_grid(var ~ ., scales = "free_y", switch = "y") +
    ggh4x::facetted_pos_scales(y = list(scale_y_continuous(), scale_y_continuous(limits = c(0.9, 1)),
                                        scale_y_log10(limits = c(50, 70)),
                                        scale_y_continuous(), scale_y_log10(limits = c(1e-7,1e-4)))) +
    labs(x = "Time (d)") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside", strip.background.y = element_blank()) +
    theme(legend.position = c(0.75, 0.7)) 



## @knitr giss-scm-susc-setup
df.giss.ppe.susc <- df.giss.ppe %>%
    filter(between(time, 3/24, 1/2)) %>%
    filter(1 - lcc < 1e-7) %>%
    group_by(time, experiment) %>%
    mutate(delta.L = lwp_ic - lwp_ic[aer == 80],
           delta.log.L = log(lwp_ic / lwp_ic[aer == 80]),
           delta.log.N = log(cdnc_ic / cdnc_ic[aer == 80]),
           delta.E_q = E_q - E_q[aer == 80],
           delta.E_theta = E_theta - E_theta[aer == 80],
           delta.log.E_q = log(E_q / E_q[aer == 80]),
           delta.log.E_theta = log(E_theta / E_theta[aer == 80]),
           delta.log.R = log(delta.R / delta.R[aer == 80]),
           S.L = delta.log.L / delta.log.N,
           S.Eq = delta.log.E_q / delta.log.N,
           S.Etheta = delta.log.E_theta / delta.log.N,
           S.R = delta.log.R / delta.log.N) %>%
    ungroup()

## @knitr giss-scm-means-setup
df.giss.ppe.means <- df.giss.ppe %>%
    ## only average E until first deepening
    mutate(across(c(E_q, E_theta, E_q.delta.q, E_theta.delta.theta),
                  ~ ifelse(jk.pbl == max(jk.pbl), .x, NA))) %>%
    filter(between(time, 2/24, 12/24)) %>%
    filter(lcc > 0.9) %>%
    group_by(aer, experiment) %>%
    summarize(across(c(lwp_ic, cdnc_ic, E_theta, E_q, E_theta.delta.theta, E_q.delta.q, delta.R),
                     .fns = list(mean = ~ mean(.x, na.rm = TRUE),
                                 sd = ~ sd(.x, na.rm = TRUE)))) %>%
    ungroup() %>%
    ## conventional units
    mutate(lwp_ic_mean_conv = lwp_ic_mean * 1000,
           lwp_ic_sd_conv = lwp_ic_sd * 1000,
           cdnc_ic_mean_conv = cdnc_ic_mean * 1e-6,
           cdnc_ic_sd_conv = cdnc_ic_sd * 1e-6) 

## @knitr giss-scm-endpoints-setup
df.giss.ppe.endpoints <- df.giss.ppe %>%
    filter(between(time, 10/24, 1/2)) %>%
    filter(lcc > 0.9) %>%
    group_by(aer, experiment) %>%
    summarize(across(c(lwp_ic, cdnc_ic, E_theta, E_q, delta.R),
                     .fns = list(mean = ~ mean(.x, na.rm = TRUE),
                                 sd = ~ sd(.x, na.rm = TRUE)))) %>%
    ungroup()

## @knitr giss-scm-tendency-setup
df.giss.ppe.tendency <- df.giss.ppe %>%
    filter(between(time, 2/24, 12/24)) %>%
    filter(lcc > 0.9) %>%
    plyr::ddply(~ experiment + aer, function(df) {
        ret <- data.frame(n = nrow(df))
        lm.fit <- try(lm(lwp_ic ~ time, df %>% mutate(time = time * 86400)))
        sum.lm.fit <- try(summary(lm.fit)$coefficients[2,2])
        if (!any(class(lm.fit) == "try-error")) 
            ret %<>% mutate(lwp_ic_tend = lm.fit$coefficients[2])
        if (!any(class(sum.lm.fit) == "try-error")) 
            ret %<>% mutate(lwp_ic_tend.se = sum.lm.fit)
    }) 

## @knitr giss-scm-ndeps-setup
df.giss.scm.ndeps <-
    bind_rows(
        ## dL/dt vs Nd
        df.giss.ppe.tendency %>%
        full_join(df.giss.ppe.means) %>%
        filter(experiment == "def") %>%
        transmute(x = cdnc_ic_mean,
                  y = lwp_ic_tend,
                  ymin = lwp_ic_tend - lwp_ic_tend.se,
                  ymax = lwp_ic_tend + lwp_ic_tend.se,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = "$\\partial\\lwp/\\partial t$~(kg~m$^{-2}$~s$^{-1}$)"),
        ## E vs Nd
        df.giss.ppe.tendency %>%
        full_join(df.giss.ppe.means) %>%
        filter(experiment == "def") %>%
        transmute(x = cdnc_ic_mean,
                  y = E_theta.delta.theta_mean,
                  ymin = E_theta.delta.theta_mean - E_theta.delta.theta_sd,
                  ymax = E_theta.delta.theta_mean + E_theta.delta.theta_sd,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["E_theta.delta.theta.hatted"]),
        ## E vs Nd
        df.giss.ppe.tendency %>%
        full_join(df.giss.ppe.means) %>%
        filter(experiment == "def") %>%
        transmute(x = cdnc_ic_mean,
                  y = E_q.delta.q_mean,
                  ymin = E_q.delta.q_mean - E_q.delta.q_sd,
                  ymax = E_q.delta.q_mean + E_q.delta.q_sd,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["E_q.delta.q.hatted"]),
        ## E vs Nd
        df.giss.ppe.tendency %>%
        full_join(df.giss.ppe.means) %>%
        filter(experiment == "def") %>%
        transmute(x = cdnc_ic_mean,
                  y = E_q_mean,
                  ymin = E_q_mean - E_q_sd,
                  ymax = E_q_mean + E_q_sd,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["E"]),
        ## L vs Nd
        df.giss.ppe.tendency %>%
        full_join(df.giss.ppe.means) %>%
        filter(experiment == "def") %>%
        transmute(x = cdnc_ic_mean,
                  y = lwp_ic_mean_conv,
                  ymin = lwp_ic_mean_conv - lwp_ic_sd_conv,
                  ymax = lwp_ic_mean_conv + lwp_ic_sd_conv,
                  xtype = "\\nd~(m$^{-3}$)",
                  ytype = tikz_replacements_unitful()["lwp_ic.conventional"])) %>%
    mutate(model = "GISS")

## @knitr giss-scm-ndeps-plot
bind_rows(df.giss.scm.ndeps) %>%
    mutate(x = x * 1e-6) %>% ## conventional Nd units
    ggplot(aes(x, y)) +
    geom_line(aes(lty = model)) +
    geom_pointrange(aes(shape = model,
                        ymin = ymin,
                        ymax = ymax)) +
    scale_x_log10() +
    facet_grid(ytype ~ ., scales = "free_y", switch = "y") +
    labs(x = tikz_replacements_unitful()["cdnc_ic.conventional"], y = "") +
    guides(shape = "none", linetype = "none") +
    theme(axis.title.y = element_blank(),
          strip.placement = "outside", strip.background.y = element_blank())

## @knitr dharma-means-setup
df.dharma.means <- df.dharma.aer %>%
    filter(between(time / 86400, 2/24, 12/24)) %>%
    filter(cloud > 0.9) %>%
    group_by(aer) %>%
    summarize(across(c(lwp_ic, cdnc_ic),
                     .fns = list(mean = ~ mean(.x, na.rm = TRUE),
                                 sd = ~ sd(.x, na.rm = TRUE)))) %>%
    ungroup() %>%
    ## conventional units
    mutate(lwp_ic_mean_conv = lwp_ic_mean * 1000,
           lwp_ic_sd_conv = lwp_ic_sd * 1000,
           cdnc_ic_mean_conv = cdnc_ic_mean * 1e-6,
           cdnc_ic_sd_conv = cdnc_ic_sd * 1e-6)

## @knitr giss-scm-les-nd-lwp
bind_rows(df.giss.ppe.means %>%
          mutate(experiment = sprintf("SCM (%s)", experiment)),
          df.dharma.means %>%
          mutate(experiment = "LES")) %>%
    ggplot(aes(cdnc_ic_mean_conv, lwp_ic_mean_conv,
               col = experiment, pch = experiment, lty = experiment)) +
    geom_line() +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_brewer("", palette = "Dark2", direction = -1) +
    scale_shape("") +
    scale_linetype("") +
    labs_nd_lwp_conventional() +
    theme(legend.position = c(0.7, 0.5))
