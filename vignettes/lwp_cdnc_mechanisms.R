## @knitr aer-sc-setup-v2
df.entrain.aer <- expand.grid(aer = c(1850, 2010),
                              exp = c("default", "sed0", "sed2")) %>%
    ddply(~ aer + exp, function(df) {
        gc()
        ## ret <- try(
        ##     with(df, 
        ##          readRDS(sprintf("../entrain_F2010A_run-v2_%daer%s_nudged_filtered2010-2015_JJA.rds",
        ##                          aer,
        ##                          ifelse(exp == "default", "", sprintf("_%s", exp)))) %>% 
        ##          left_join(readRDS("../ne30pg2.rds")) %>%
        ##          filter(LANDFRAC == 0) %>%
        ##          filter(fscu.annual > 0.3) %>%
        ##          mutate(region = regionalize(month, lon, lat, fscu.annual, LANDFRAC)) %>%
        ##          filter(!is.na(region))))
        ## if (!("try-error" %in% class(ret)))
        ##     return(ret)
        ret <- try(
            with(df, 
                 readRDS(sprintf("../entrain_F2010A_run-v2_%daer%s_nudged_filtered.rds",
                                 aer,
                                 ifelse(exp == "default", "", sprintf("_%s", exp)))) %>% 
                 left_join(readRDS("../ne30pg2.rds")) %>%
                 filter(LANDFRAC == 0) %>%
                 filter(fscu.annual > 0.3) %>%
                 mutate(region = regionalize_seasonal(month, lon, lat, fscu.annual, LANDFRAC)) %>%
                 filter(!is.na(region))))
        if (!("try-error" %in% class(ret)))
            return(ret)
        else
            return(NULL)
    }) 
        
df.entrain.aer.golden <- df.entrain.aer %>%
    filter(between(73 - jk.pbl, 10, 15)) %>%
    filter(abs(lat) < 40)

## @knitr aer-sc-noprecsup-setup-v2
df.entrain.aer.noprecsup.golden <-
    expand.grid(aer = c("D", "I"),
                sed = c("D", "N", "G")) %>%
    ddply(~ aer + sed, function(df) {
        gc()
        with(df, 
             readRDS(sprintf("../data/perlmutter/entrain_F2010A_run-rEn0aNc%sp%so2_golden.rds",
                             sed, aer)))
    }) %>%
    filter(between(73 - jk.pbl, 10, 15)) %>%
    filter(abs(lat) < 40) %>%
    mutate(aer = sprintf("P%s", aer),
           sed = revalue(sed, c("D" = 1,
                                "N" = 0,
                                "G" = 2)))

## @knitr aer-noprecsup-setup-v2
df.aer.noprecsup <-
    expand.grid(aer = c("D", "I"),
                sed = c("D")) %>% ## sed = c("D", "N", "G")) %>%
    ddply(~ aer + sed, function(df) {
        gc()
        with(df, 
             readRDS(sprintf("../data/perlmutter/e3smv2_aer-rEn0aNc%sp%so2.rds",
                             sed, aer)))
    }) %>%
    mutate(aer = sprintf("P%s", aer),
           sed = revalue(sed, c("D" = 1,
                                "N" = 0,
                                "G" = 2))) %>%
    rename(lts = th7001000) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower)) %>%
    filter_warm() %>%
    filter_ocean() %>%
    calculate_incloud()

## @knitr aer-noprecsup-v2-lwpratio
aer.noprecsup.lwp.pd.pi <- with(df.aer.noprecsup %>%
                                filter(lcc == 1, precc == 0, icc == 0, ttop > 273) %>%
                                group_by(aer) %>%
                                summarize(across(c(lwp, lwp_ic, cdnc, cdnc_ic, lcc), mean)) %>%
                                ungroup(),
                                log(lwp_ic[aer == "PD"]) - log(lwp_ic[aer == "PI"])) %>%
    sprintf("%.1e", .) %>%
    tikz_sanitize() %>%
    gsub("\\$", "", .)

## @knitr aer-sc-setup-v1
df.entrain.aer.v1 <- expand.grid(aer = 2010,
                                 exp = c("default", "eul", "noprecsup")) %>%
    ddply(~ aer + exp, function(df) {
        gc()
        with(df, 
             readRDS(sprintf("../entrain_F2010A_run-EAMv1_nudged_2010%s_golden.rds",
                             ifelse(exp == "default", "", sprintf("_%s", exp)))) %>% 
             left_join(readRDS("../ne30pg2.rds")) %>%
             filter(LANDFRAC == 0) %>%
             filter(fscu.annual > 0.3) %>%
             mutate(region = regionalize_seasonal(month, lon, lat, fscu.annual, LANDFRAC)) %>%
             filter(!is.na(region)))
    }) 
        
df.entrain.aer.golden.v1 <- df.entrain.aer.v1 %>%
    filter(between(73 - jk.pbl, 10, 15)) %>%
    filter(abs(lat) < 40)

## @knitr precsup-conditional
df.e3sm.precsup.cond <- bind_rows(readRDS("../e3smv1_aer-EAMv1_nudged_2010.rds") %>%
                                 mutate(exp = "E3SM"),
                                 readRDS("../e3smv1_aer-EAMv1_nudged_2010_noprecsup.rds") %>%
                                 mutate(exp = "E3SM-noprecsup")) %>%
    rename(lts = th7001000) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower)) %>%
    filter_warm() %>%
    filter_ocean() %>%
    calculate_incloud() %>%
    group_by(exp) %>%
    lwp.cdnc.conditional() %>%
    ungroup()

## @knitr precsup-precip-conditional
df.e3sm.precsup.precip.cond <- bind_rows(readRDS("../e3smv1_aer-EAMv1_nudged_2010.rds") %>%
                                 mutate(exp = "E3SM"),
                                 readRDS("../e3smv1_aer-EAMv1_nudged_2010_noprecsup.rds") %>%
                                 mutate(exp = "E3SM-noprecsup")) %>%
    rename(lts = th7001000) %>%
    left_join(readRDS("../ne30pg2.rds") %>%
              rename_with(tolower)) %>%
    filter_warm() %>%
    filter_ocean() %>%
    filter_sc() %>%
    calculate_incloud() %>%
    mutate(precip = 1e3 * 86400 * (precc + precl)) %>% ## m s^-1 --> mm d^-1
    mutate(cloudsat = factor(ifelse(precip < 1e-2,
                                    "$R < 10^{-2}$~mm~d$^{-1}$",
                                    "$R\\geq 10^{-2}$~mm~d$^{-1}$"),
                             levels = c("$R < 10^{-2}$~mm~d$^{-1}$",
                                        "$R\\geq 10^{-2}$~mm~d$^{-1}$"))) %>%
    discretize(precip, 6, equal_contents = TRUE, as_factor = TRUE) 
    group_by(exp, precip) %>%
    lwp.cdnc.conditional() %>%
    ungroup()

## @knitr entrain-susc-v2-noprecsup-linear-setup
df.summary <- df.entrain.aer.noprecsup.golden %>%
    filter(precc == 0) %>%
    mutate(period = aer) %>%
    ## filter(period == 2010) %>%
    ## filter(aer == "PD") %>%
    mutate(exp = sprintf("sed%s", sed)) %>%
    calculate_incloud() %>%
    calculate_fluxes() %>%
    convert_omega() %>%
    ## mutate(log10.delta.R = log10(delta.R)) %>%
    ## mutate(discretized_R = 86400 * delta.R) %>% ## m s^-1 --> mm d^-1
    ## discretize(discretized_R, 4, as_factor = TRUE, equal_contents = TRUE) %>%
    filter(is.finite(lwp_ic)) %>%
    discretize_lwp(seq(0, 1, 0.025)) %>%
    ##       discretize_lwp(exp(seq(log(0.001), log(1), length.out = 25))) %>%
    ## discretize_cdnc(exp(seq(log(1e7), log(3e8), length.out = 10))) %>%
    mutate(E = E_q) %>%
    mutate(period = factor(period)) %>%
    filter(lcc == 1, icc == 0, ttop > 273)

df.susc <- df.summary %>% ## delta.theta > 5) %>%
    ## filter(E > 0) %>%
    plyr::ddply(~ exp + period + lwp_ic + region, function(df) {
        ret <- data.frame(n = nrow(df))
        lm.fit <- try(lm(log(E) ~ log(cdnc_ic), df))
        sum.lm.fit <- try(summary(lm.fit)$coefficients[2,2])
        if (!any(class(lm.fit) == "try-error")) 
            ret %<>% mutate(S.E = lm.fit$coefficients[2])
        if (!any(class(sum.lm.fit) == "try-error")) 
            ret %<>% mutate(S.E.se = sum.lm.fit)
        ret %>% mutate(E = mean(df$E, trim = 0.05),
                       cdnc_ic = mean(df$cdnc_ic, trim = 0.05))
    })

df.lwp.susc <- df.summary %>% ## delta.theta > 5) %>%
    plyr::ddply(~ exp + period + lwp_ic + region, function(df) {
        ret <- data.frame(n = nrow(df))
        lm.fit <- try(lm(rhoH.dq_t.0.dt ~ log(cdnc_ic), df))
        sum.lm.fit <- try(summary(lm.fit)$coefficients[2,2])
        if (!any(class(lm.fit) == "try-error")) 
            ret %<>% mutate(S.L = lm.fit$coefficients[2])
        if (!any(class(sum.lm.fit) == "try-error")) 
            ret %<>% mutate(S.L.se = sum.lm.fit)
        ret %>% mutate(E = mean(df$E, trim = 0.05),
                       cdnc_ic = mean(df$cdnc_ic, trim = 0.05))
    })

df.lwp.susc.rh <- df.summary %>% ## delta.theta > 5) %>%
    plyr::ddply(~ exp + period + region, function(df) {
        ## calculate RH just above the inversion
        df %>%
            mutate(rh.plus = q_t.plus / qsat(theta_l.plus * (p.pbl / 1e5) ^ (2/7), p.pbl)) %>%
            filter(is.finite(rh.plus), between(rh.plus, 0, 1)) %>%
            discretize(rh.plus, 3, equal_contents = TRUE)
    }) %>%
    plyr::ddply(~ exp + period + lwp_ic + region + rh.plus, function(df) {
        ret <- data.frame(n = nrow(df))
        lm.fit <- try(lm(rhoH.dq_t.0.dt ~ log(cdnc_ic), df))
        sum.lm.fit <- try(summary(lm.fit)$coefficients[2,2])
        if (!any(class(lm.fit) == "try-error")) 
            ret %<>% mutate(S.L = lm.fit$coefficients[2])
        if (!any(class(sum.lm.fit) == "try-error")) 
            ret %<>% mutate(S.L.se = sum.lm.fit)
        ret %>% mutate(E = mean(df$E, trim = 0.05),
                       cdnc_ic = mean(df$cdnc_ic, trim = 0.05))
    })

df.diff <- df.summary %>%
    group_by(exp, region, period, lwp_ic) %>%
    summarize(E = mean(E, trim = 0.05),
              cdnc = mean(cdnc_ic, trim = 0.05),
              n = n()) %>%
    ungroup() %>%
    gather(var, val, E : n) %>%
    spread(period, val) %>%
    ##    mutate(dlog = PD / PI - 1) %>%
    mutate(dlog = log(PD) - log(PI)) %>%
    select(-c(PD, PI)) %>%
    spread(var, dlog) %>%
    mutate(S.E = E / cdnc)
