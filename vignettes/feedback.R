## @knitr regionalize
regionalize <- function(month, lon, lat, fscu.annual, landfrac) {
    ifelse(landfrac > 0 | fscu.annual <= 0.3 | abs(lat) > 45,
           NA,
    ifelse(lat > 0 & month %in% c("Jun", "Jul", "Aug"),
    ifelse(lon < 30 | lon > 330,
           "NEA",
           "NEP"),
    ifelse(lat < 0 & month %in% c("Oct", "Nov", "Dec", "Jan", "Feb"),
    ifelse(lon < 30 | lon > 330,
           "SEA",
    ifelse(lon > 60 & lon < 120,
           "AUS",
           "SEP")),
    NA)))
}

## @knitr amip-sc-setup
df.amip.fscu <- left_join(readRDS("../entrain_F2010A_run-v2_F2010A4_00_fscu.rds"),
                          readRDS("../ne30pg2.rds"))

df.entrain.amip <- expand.grid(amip = c(0, 4),
                               exp = c(0, 8)) %>%
    ddply(~ amip + exp, function(df) {
        gc()
        with(df, 
             readRDS(sprintf("../entrain_F2010A_run-v2_F2010A%d_%02d_filtered.rds", amip, exp)) %>%
             left_join(readRDS("../ne30pg2.rds")) %>%
             filter(LANDFRAC == 0) %>%
             filter(fscu.annual > 0.3) %>%
             mutate(region = regionalize(month, lon, lat, fscu.annual, LANDFRAC)) %>%
             filter(!is.na(region)) %>%
             mutate(period = if (amip == 0) "AMIP" else sprintf("AMIP+%d K", amip)))
    }) 

df.entrain.amip.golden <- df.entrain.amip %>%
    filter(between(73 - jk.pbl, 6, 15)) %>%
    filter(abs(lat) < 40) %>%
    filter(icc == 0)

## @knitr amip-sc-monthly-stats
df.entrain.amip.golden.monthly.stats <- df.entrain.amip.golden %>%
    ## only do default parameter model run
    ## filter(exp == 0) %>%
    ## calculate year
    mutate(year = as.POSIXlt(time)$year + 1900) %>%
    ## do a bit of renaming
    mutate(theta_adv = rhoH.dthetacore.0 / rhoH,
           q_adv = rhoH.dqcore.0 / rhoH) %>%
    calculate_fluxes() %>%
    convert_omega() %>%
    ## 
    group_by(region, period, exp, month, year) %>%
    summarize(across(c(E_theta, E_q, E_theta.delta.theta, E_q.delta.q,
                       lcc, lwp, omega500.hPa.d, omega700.hPa.d, lts,
                       theta_adv, q_adv,
                       theta_l.0, q_t.0,
                       delta.q, delta.theta,
                       sh, lh, delta.F, delta.R, SST_cpl),
                     list(median = median,
                          trimmean = ~ mean(.x, trim = 0.05))),
              n = n()) %>%
    ungroup() 

## @knitr amip-sc-monthly-stats-json-output
library(jsonlite)
write_json(df.entrain.amip.golden.monthly.stats, "sc_entrain_monthly.json", digits = NA, pretty = TRUE)
# write_json(as.data.frame(df.entrain.amip.golden.monthly.stats %>% filter(exp == 0)), "sc_entrain_monthly.json")
write.table(df.entrain.amip.golden.monthly.stats %>% filter(exp == 0), "sc_entrain_monthly.txt", quote = TRUE, row.names = FALSE, sep = ",")
saveRDS(df.entrain.amip.golden.monthly.stats, "sc_entrain_monthly.rds" )
