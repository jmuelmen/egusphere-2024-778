## @knitr E-nd-L-lineplot
df.entrain.aer.noprecsup.golden %>%
    filter(precc == 0) %>%
    mutate(period = aer) %>%
    filter(period == "PD") %>%
    filter(sed == 1) %>%
    mutate(exp = sed) %>%
    ## filter(time > "2010-01-01", time < "2011-01-01") %>%
    calculate_incloud() %>%
    calculate_fluxes() %>%
    convert_omega() %>%
    filter(is.finite(cdnc_ic),
           is.finite(lwp_ic)) %>%
    ## mutate(log10.delta.R = log10(delta.R)) %>%
    ## mutate(discretized_R = 86400 * delta.R) %>% ## m s^-1 --> mm d^-1
    ## discretize(discretized_R, 4, as_factor = TRUE, equal_contents = TRUE) %>%
    discretize_lwp(seq(0, 0.25, 25e-3)) %>%
    discretize(cdnc_ic, 6, equal_contents = TRUE) %>% ## (seq(0, 3e8, 25e6)) %>%
    mutate(E = E_q) %>%
    filter(is.finite(E),
           is.finite(lwp_ic)) %>%
    filter(between(E, quantile(E, 0.02), quantile(E, 0.98))) %>%
    mutate(period = factor(period)) %>%
    filter(region == "NEP") %>%
    filter(lcc > 0.9, icc == 0, ttop > 273) %>% ## delta.theta > 5) %>%
    group_by(exp, period, cdnc_ic, lwp_ic, region) %>%
    filter(n() > 100) %>%
    summarize(n = n(),
              se.E = sd(E, na.rm = TRUE) / sqrt(n),
              E = mean(E, na.rm = TRUE)) %>%
    ungroup(cdnc_ic) %>%
    mutate(ntot = sum(n)) %>%
    ggplot(aes(cdnc_ic, E, col = lwp_ic, group = lwp_ic)) +
    ## lwp.cdnc.conditional(lwp.bins = seq(-0.1, 0.1, length.out = 50 * 2),
    ##                      cdnc.bins = seq(0, 3e8, 25e6)) %>% ## exp(seq(log(1e6), log(300e6), length.out = 25))) %>%
    ## filter(count.cdnc > 25) %>%
    ## lwp.cdnc.plot(marginal = FALSE) +
    geom_line() +
    geom_point(aes(size = n)) + 
    geom_linerange(aes(ymin = E - 1.96 * se.E,
                       ymax = E + 1.96 * se.E)) + 
    ## scale_x_log10() +
    ## scale_y_log10() +
    ## scale_y_continuous(breaks = (3:5) * 1e-3) +
    coord_cartesian(xlim = c(10e6, 175e6),
                    expand = FALSE) +
    labs(x = "\\nd~(m$^{-3}$)",
         y = "$E$ (kg m$^{-2}$ s$^{-1}$)") +
    scale_size("$n$") +
    scale_color_distiller("\\lwp~(kg~m$^{-2}$)", palette = "PuRd", direction = -1) +
    guides(color = guide_colorbar(direction = "horizontal", title.vjust = 0.75, barwidth = 8)) +
    theme(legend.direction = "horizontal",
          legend.position = c(0.65, 0.35)) 

## @knitr nd-lwp-v1
df.entrain.aer.golden.v1 %>%
    ## filter(precc == 0) %>%
    filter(exp %in% c("default", "noprecsup")) %>%
    mutate(period = aer) %>%
    filter(period == 2010) %>%
    mutate(aer = sprintf("%d", aer),
           period = sprintf("%d", period)) %>%
    mutate(version = "v1") %>%
    calculate_incloud() %>%
    calculate_fluxes() %>%
    convert_omega() %>%
    ## filter(time > "2010-01-01", time < "2011-01-01") %>%
    filter(region == "NEP") %>%
    filter(lcc > 0.9, icc == 0) %>% ## , ttop > 273) %>% ## delta.theta > 5) %>%
    ## mutate(log10.delta.R = log10(delta.R)) %>%
    ## mutate(discretized_R = 86400 * delta.R) %>% ## m s^-1 --> mm d^-1
    ## discretize(discretized_R, 4, as_factor = TRUE, equal_contents = TRUE) %>%
    filter(is.finite(lwp_ic),
           is.finite(cdnc_ic)) %>%
    group_by(exp, period) %>%
    lwp.cdnc.conditional(quantile(cdnc_ic, seq(0, 1, 0.04))) %>%
    lwp.cdnc.plot(col = exp)

## @knitr nd-lwp-v2
bind_rows(df.entrain.aer.golden.v1 %>%
          ## filter(precc == 0) %>%
          mutate(period = aer) %>%
          filter(period == 2010) %>%
          filter(exp == "noprecsup") %>%
          mutate(aer = sprintf("%d", aer),
                 period = sprintf("%d", period)) %>%
          mutate(version = "v1"),
          df.entrain.aer.noprecsup.golden %>%
          ## filter(precc == 0) %>%
          mutate(period = aer) %>%
          filter(period == "PD") %>%
          filter(sed == 1) %>%
          mutate(exp = sed) %>%
          mutate(version = "v2")) %>%
    calculate_incloud() %>%
    calculate_fluxes() %>%
    convert_omega() %>%
    ## filter(time > "2010-01-01", time < "2011-01-01") %>%
    filter(region == "NEP") %>%
    filter(lcc > 0.9, icc == 0) %>% ## , ttop > 273) %>% ## delta.theta > 5) %>%
    ## mutate(log10.delta.R = log10(delta.R)) %>%
    ## mutate(discretized_R = 86400 * delta.R) %>% ## m s^-1 --> mm d^-1
    ## discretize(discretized_R, 4, as_factor = TRUE, equal_contents = TRUE) %>%
    filter(is.finite(lwp_ic),
           is.finite(cdnc_ic)) %>%
    group_by(version, period) %>%
    lwp.cdnc.conditional(quantile(cdnc_ic, seq(0, 1, 0.04))) %>%
    lwp.cdnc.plot(col = version)

## @knitr e3sm-noprecsup-v2-entrainment-ecdfs
df.entrain.aer.noprecsup.golden %>%
    filter(sed == 1) %>%
    mutate(period = aer) %>%
    filter(precc == 0) %>%
    filter(lcc == 1, icc == 0, ttop > 273) %>% ## delta.theta > 5) %>%
    calculate_incloud() %>%
    calculate_fluxes() %>%
    filter_sc() %>%
    mutate(rh.plus = q_t.plus / qsat(theta_l.plus * (p.pbl / 1e5) ^ (2/7), p.pbl)) %>%
    filter(rh.plus <= 1, rh.plus >= 0) %>%
    group_by(sed, region) %>%
    ## discretize(rh.plus, 3, equal_contents = TRUE) %>%
    mutate(rh.plus = cut(rh.plus, 4, include.lowest = TRUE)) %>% ## , labels = FALSE, dig.lab = 3)) %>%
    ungroup() %>%
    mutate(lwp.bins = cut(lwp_ic, c(0,0.02,0.05,0.1,0.2,1), include.lowest = TRUE)) %>% ## , labels = FALSE, dig.lab = 3)) %>%
    mutate(cdnc.bins = cut(cdnc_ic, c(0,0.02,0.05,0.1,0.2) * 1e9, include.lowest = TRUE)) %>% ## , labels = FALSE, dig.lab = 3)) %>%
    filter(between(E_q, -1e-1, 1e-1)) %>%
    filter(!is.na(cdnc.bins)) %>%
    ggplot(aes(E_q)) +
    stat_ecdf(aes(col = lwp.bins, lty = period)) +
    ## or
    ## stat_ecdf(aes(col = cdnc.bins, lty = period)) +
    facet_grid(. ~ region) +
    scale_x_continuous(limits = c(-0.01, 0.01))

## @knitr e3sm-noprecsup-v2-susceptibility.maps
df.entrain.aer.noprecsup.golden %>%
    filter(sed == 1) %>%
    mutate(period = aer) %>%
    filter(precc == 0) %>%
    filter(lcc == 1, icc == 0, ttop > 273) %>% ## delta.theta > 5) %>%
    calculate_incloud() %>%
    calculate_fluxes() %>%
    filter_sc() %>%
    mutate(rh.plus = q_t.plus / qsat(theta_l.plus * (p.pbl / 1e5) ^ (2/7), p.pbl)) %>%
    filter(rh.plus <= 1, rh.plus >= 0) %>%
    ## discretize(rh.plus, 3, equal_contents = TRUE) %>%
    mutate(rh.plus = cut(rh.plus, 4, include.lowest = TRUE)) %>% ## , labels = FALSE, dig.lab = 3)) %>%
    mutate(lon = 10 * lon %/% 10,
           lat = 10 * lat %/% 10) %>%
    group_by(sed, region, period, lon, lat, rh.plus) %>%
    summarize(lwp = mean(log(lwp_ic), trim = 0.05),
              cdnc = mean(log(cdnc_ic), trim = 0.05),
              E_theta = mean(E_theta[lwp_ic > 0.05], trim = 0.05),
              E_q = mean(E_q[lwp_ic > 0.05], trim = 0.05),
              n = n()) %>%
    ungroup() %>%
    gather(var, val, lwp : n) %>%
    spread(period, val) %>%
    ##    mutate(dlog = PD / PI - 1) %>%
    mutate(dlog = PD - PI) %>%
    select(-c(PD, PI)) %>%
    spread(var, dlog) %>%
    mutate(S.L = lwp / cdnc) %>%
    ggplot(aes(lon, lat, fill = sign(lwp))) +
    geom_raster() +
    facet_grid(rh.plus ~ .) +
    scale_fill_distiller(palette = "RdBu") # , limits = c(-1, 1) * 1e-2)
