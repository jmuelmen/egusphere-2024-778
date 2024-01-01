#' @export
#' 
tikz_minus <- function(x) sub("^-", "$-$", format(x))

#' @export
#' 
tikz_sanitize <- function(x) plotutils::sanitize.numbers(x, "latex", TRUE, TRUE)

#' @export
#' 
tikz_sanitize_sparse <- function(x) {
    ## print(x)
    ret <- tikz_sanitize(x)
    print(ret)
    if (length(na.omit(ret)) > 2) {
        if (length(na.omit(ret)) %% 2) {
            ret[seq(2 + is.na(ret[1]), length(ret), by = 2)] <- ""
        } else {
            ret[seq(1 + is.na(ret[1]), length(ret), by = 2)] <- ""
        }
    }
    print(ret)
    ret
}

# change the default scales
#' @export
#' 
scale_x_continuous <- function(..., labels=tikz_sanitize)
    ggplot2::scale_x_continuous(..., labels=labels)

#' @export
#' 
scale_x_reverse <- function(..., labels=tikz_sanitize)
    ggplot2::scale_x_reverse(..., labels=labels)

#' @export
#' 
scale_x_log10 <- function(..., labels=tikz_sanitize)
    ggplot2::scale_x_log10(..., labels=labels)

#' @export
#' 
scale_y_continuous <- function(..., labels=tikz_sanitize)
    ggplot2::scale_y_continuous(..., labels=labels)

#' @export
#' 
scale_y_log10 <- function(..., labels=tikz_sanitize)
    ggplot2::scale_y_log10(..., labels=labels)

#' @export
#' 
tikz_replacements_unitful <- function() {
    c(delta.F = "$\\Delta F$ (W~m$^{-2}$)",
      lh = "$\\text{LH}$ (W~m$^{-2}$)",
      sh = "$\\text{SH}$ (W~m$^{-2}$)",
      ctcool = "$\\Delta F$ (in-cloud, W~m$^{-2}$)",
      t = "$T$ (\\textdegree C)",
      qv = "$\\Delta q\\ (\\text{kg kg}^{-1})$",
      ql = "$q_l\\ (\\text{g kg}^{-1})$",
      q_l.0 = "$\\hat{q}_l\\ (\\text{g kg}^{-1})$",
      q_t.0 = "$\\hat{q}_t\\ (\\text{g kg}^{-1})$",
      q_t.plus = "$q_t^+$",
      ql_ic = "$q_l/f\\ (\\text{kg kg}^{-1})$",
      qn_ic = "$q_N/f\\ (\\text{kg kg}^{-1})$",
      qt = "$\\Delta q_t\\ (\\text{g kg}^{-1})$",
      lts = "LTS (K)",
      omega = "$\\omega$ (Pa~s$^{-1}$)",
      omega500.hPa.d = "$\\omega_{500}$ (hPa d$^{-1}$)",
      omega700.hPa.d = "$\\omega_{700}$ (hPa d$^{-1}$)",
      omega.pbl.hPa.d = "$\\omega_\\text{PBL}$ (hPa d$^{-1}$)",
      theta = "$\\Delta \\theta\\ (\\text{K})$",
      theta_l = "$\\Delta \\theta_l\\ (\\text{K})$",
      theta_l.0 = "$\\hat{\\theta}_l\\ (\\text{K})$",
      theta_l.plus = "$\\theta_l^+\\ (\\text{K})$",
      theta_e = "$\\theta_e\\ (\\text{K})$",
      aclcac = "$f$",
      cloud = "$f$",
      f.minus = "$f$",
      h = "$h$ (m)",
      cloudfrac_clubb = "$f_\\text{CLUBB}$",
      numliq = "$q_N$",
      numliq_ic = "$q_N / f$",
      cdnc = "$\\nd$ (m$^{-3}$)",
      cdnc_ic = "$\\nd$ (m$^{-3}$)",
      cdnc_ic.conventional = "$\\nd$ (cm$^{-3}$)",
      lwp = "$\\lwp$ (kg m$^{-2}$)",
      lwp_ic = "$\\lwp$ (kg m$^{-2}$)",
      lwp_ic.conventional = "$\\lwp$ (g m$^{-2}$)",
      ccn = "CCN ($S=0.2\\%$, cm$^{-3}$)",
      adv.xl = "${\\bf v}\\cdot\\nabla q_l$",
      adv.q = "${\\bf v}\\cdot\\nabla q$",
      wp2_clubb = "$\\overline{w'^2}\\ (\\text{m}^2\\ \\text{s}^{-2})$",
      wp3_clubb = "$\\overline{w'^3}\\ (\\text{m}^3\\ \\text{s}^{-3})$",
      skw = "Sk$_w$",
      delta.theta = "$\\Delta\\theta_l$ (K)",
      delta.q = "$\\Delta q_t$ (kg~kg$^{-1}$)",
      E = "$E$ (kg~m$^{-2}$~s$^{-1}$)",
      E_theta = "$E_\\theta$ (kg~m$^{-2}$~s$^{-1}$)",
      E_q = "$E_q$ (kg~m$^{-2}$~s$^{-1}$)",
      E_theta.delta.theta = "$E_\\theta(\\theta_l^+ - \\hat\\theta_l)$ (kg~K~m$^{-2}$~s$^{-1}$)",
      E_q.delta.q = "$E_q(q_t^+ - \\hat q_t)$ (kg~m$^{-2}$~s$^{-1}$)",
      E_theta.noadv = "$E_\\theta$ (no adv)",
      E_q.noadv = "$E_q$ (no adv)",
      E_h_theta = "$E_{h_\\theta}$",
      E_h_q = "$E_{h_q}$",
      E_p_theta = "$E_{p_\\theta}$",
      E_p_q = "$E_{p_q}$",
      dh.dt_lagged = "lagged $\\ddp{h}{t}$",
      ## EAMxx additions
      cldfrac_liq = "$f_\\text{liq}$",
      nc = "$N_c$",
      ql = "$q_l$~(kg~kg$^{-1}$)",
      qc = "$q_c$~(kg~kg$^{-1}$)",
      qi = "$q_i$~(kg~kg$^{-1}$)",
      delta.R = "$\\Delta R$~(kg m$^{-2}$ s$^{-1}$)",
      precip = "$R$~(kg m$^{-2}$ s$^{-1}$)"
      )
}

#' @export
#' 
tikz_replacements_unitless <- function() {
    c(delta.F.cp = "$\\Delta F/c_p$",
      ## Dtheta_l.0.Dt = "$\\DDP{\\hat\\theta_l}{t}$",
      rhoH.Dtheta_l.0.Dt = "$\\rho H(\\DDP{\\hat\\theta_l}{t})$",
      rhoH.dtheta_l.0.dt = "$\\rho H(\\ddp{\\hat\\theta_l}{t})$",
      rhoH.dthetacore.0 = "$-\\rho H({\\bf v}\\cdot\\nabla\\hat\\theta_l)$",
      rhoH.adv.theta_l.0 = "$\\rho H({\\bf v}\\cdot\\nabla\\hat\\theta_l)$",
      rhoH.dthetadt.hdiff.0 = "$\\rho H(\\ddp{\\hat\\theta_l}{t})|_{\\text{diff}}$",
      ## rhoH.dtheta_l.0.dt = "$\\rho H(\\ddp{\\hat\\theta_l}{t})$",
      sh.cp = "$\\text{SH}/c_p$",
      ## delta.theta = "$\\theta_l^+ - \\hat\\theta_l$",
      E_theta.delta.theta = "$E_\\theta(\\theta_l^+ - \\hat\\theta_l)$",
      ## E_theta = "$E_\\theta$"))
      delta.R = "$\\Delta R$",
      ## Dq_t.0.Dt = "$\\DDP{\\hat q_t}{t}$",
      rhoH.Dq_t.0.Dt = "$\\rho H(\\DDP{\\hat q_t}{t})$",
      rhoH.dq_t.0.dt = "$\\rho H(\\ddp{\\hat q_t}{t})$",
      rhoH.dqcore.0 = "$-\\rho H({\\bf v}\\cdot\\nabla\\hat q_t)$",
      rhoH.adv.q_t.0 = "$\\rho H({\\bf v}\\cdot\\nabla\\hat q_t)$",
      ## rhoH.dq_t.0.dt = "$\\rho H(\\ddp{\\hat q_t}{t})$",
      qflx = "$\\text{LH}/L_v$",
      ## delta.q = "$q_t^+ - \\hat q_t$",
      E_q.delta.q = "$E_q(q_t^+ - \\hat q_t)$",
      ## E_q = "$E_q$"))
      E_theta = "$E_\\theta$",
      E_q = "$E_q$",
      E_theta.noadv = "$E_\\theta$ (no adv)",
      E_q.noadv = "$E_q$ (no adv)",
      E_h_theta = "$E_{h_\\theta}$",
      E_h_q = "$E_{h_q}$",
      E_p_theta = "$E_{p_\\theta}$",
      E_p_q = "$E_{p_q}$")
}
