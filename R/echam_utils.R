#' Rename variables
#'
#' @export
rename_fields_echam <- function(df) {
    gc()
    df %>%
        dplyr::rename(lts = LTS,
                      omega500 = omega.500,
                      omega700 = omega.700)
}

#' Scale variables with unexpected units
#'
#' @export
rescale_fields_echam <- function(df) {
    gc()
    df %>%
        dplyr::mutate(dthetadt.rad.0 = dthetadt.rad.0 / 86400,
                      dthetadt.hdiff.0 = dthetadt.hdiff.0 / 86400) 
}

#' Replace certain fields with the canonical combinations used in our
#' chosen form of the budget equations
#'
#' @export
canonize_fields_echam <- function(df) {
    gc()
    df %>%
        dplyr::mutate(rhoH.dtheta_l.0.dt = rhoHg * dtheta_l.0.dt / g,
                      delta.F.cp = -dthetadt.rad.0 * rhoHg / g,
                      sh.cp = -sh / cp,
                      delta.R = (aprl + aprc),
                      rhoH.dq_t.0.dt = rhoHg * dq_t.0.dt / g,
                      qflx = -lh / Lv,
                      rhoH.dthetadt.hdiff.0 = rhoHg * dthetadt.hdiff.0 / g,
                      rhoH.adv.theta_l.0 = rhoHg * adv.theta_l.0 /g,
                      rhoH.adv.q_t.0 = rhoHg * adv.q_t.0 /g) %>%
        dplyr::select(-c(dtheta_l.0.dt, dthetadt.rad.0, dthetadt.hdiff.0, adv.theta_l.0,
                         sh, lh, aprl, aprc,
                         dq_t.0.dt, adv.q_t.0,
                         rhoHg))
}

#' @export
calculate_budgets_echam <- function(df) {
    df %>%
        dplyr::mutate(delta.theta = theta_l.plus - theta_l.0,
                      delta.q = q_t.plus - q_t.0,
                      ## without advection, otherwise following Kalmus et al. eq. (10)--(11)
                      E_theta.delta.theta.noadv = rhoH.dtheta_l.0.dt + delta.F.cp - sh.cp - delta.R * Lv / cp,
                      E_q.delta.q.noadv = rhoH.dq_t.0.dt - qflx + delta.R,
                      E_theta.noadv = E_theta.delta.theta.noadv / delta.theta,
                      E_q.noadv = E_q.delta.q.noadv / delta.q,
                      ## with horizontal advection and diffusion (note dthetadt.hdiff is of
                      ## opposite sign to the advective terms)
                      E_theta.delta.theta = E_theta.delta.theta.noadv + rhoH.adv.theta_l.0, ## - rhoH.dthetadt.hdiff.0,
                      E_q.delta.q = E_q.delta.q.noadv + rhoH.adv.q_t.0,
                      E_theta = E_theta.delta.theta / delta.theta,
                      E_q = E_q.delta.q / delta.q)
}
