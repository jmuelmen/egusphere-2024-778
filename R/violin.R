#' @export
#' 
violin.theta <- function(df.budget, geom = geom_violin(draw_quantiles = 0.5, col = PNNLonyx, fill = PNNLcopper)) {
    df.budget %>%
        tidyr::gather(term, x, any_of(c("delta.F.cp",
                                        "rhoH.dtheta_l.0.dt",
                                        "sh.cp",
                                        "rhoH.dthetacore.0",
                                        "rhoH.adv.theta_l.0",
                                        "rhoH.vadv.theta_l.0",
                                        ## rhoH.dthetadt.hdiff.0,
                                        "E_theta.delta.theta"))) %>%
        mutate(
            variable = "$\\theta_l$ budget terms",
            ## ## preserve order
            ## term = factor(term, levels = 
            ##                         names(df.theta)[which(names(df.theta) == "delta.F.cp") : 
            ##                                         which(names(df.theta) == "E_theta.delta.theta")]), 
            ## make labels fancy
            term = revalue(term, 
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
                             E_theta.delta.theta = "$E_\\theta(\\theta_l^+ - \\hat\\theta_l)$" ##,
                             ## E_theta = "$E_\\theta$"))
                             ))) %>%
        ## group_by(term) %>%
        filter(between(x, quantile(x, 0.01, na.rm = TRUE), quantile(x, 0.99, na.rm = TRUE))) %>%
        ## ungroup() %>%
        ggplot(aes(x = term, y = x)) +
        ## geom_violin(aes(fill = factor(2 * ceiling((72 - jk.pbl) / 2))), draw_quantiles = 0.5) + ## , col = PNNLonyx, fill = PNNLcopper) +
        ## geom_violin(draw_quantiles = 0.5, col = PNNLonyx, fill = PNNLcopper) +
        geom + 
        ## geom_boxplot() +
        ##scale_y_continuous(labels = tikz_sanitize) +
        ## scale_y_continuous(labels = tikz_sanitize_sparse, limits = c(-0.25, 0.25)) +
        labs(x = "", y = "kg K m$^{-2}$ s$^{-1}$") +
        ##ylim(-0.2, 0.2) +
        geom_hline(yintercept = 0, lty = "dashed", col = "lightgrey") +
        facet_grid(variable ~ .) +
        theme_bw() 
}

#' @export
#' 
violin.q <- function(df.budget, geom = geom_violin(draw_quantiles = 0.5, col = PNNLonyx, fill = PNNLcopper),
                     clip.delta.R = 1e-6, clip.quantiles = 1e-2) {
    if (!is.null(clip.delta.R)) {
        df.budget %<>%
            mutate(delta.R = replace(delta.R, delta.R < clip.delta.R, NA))
    }

    df.budget %<>%
        tidyr::gather(term, x, any_of(c("delta.R",
                                        "rhoH.dq_t.0.dt",
                                        "qflx",
                                        "rhoH.dqcore.0",
                                        "rhoH.adv.q_t.0",
                                        "rhoH.vadv.q_t.0",
                                        "E_q.delta.q"))) %>%
        ## group_by(term) %>%
        ## mutate(p.min = quantile(x, 0.02, na.rm = TRUE), 
        ##        p.max = quantile(x, 0.98, na.rm = TRUE)) %>%
        ## filter(x >= p.min, x <= p.max) %>%
        ## mutate(z = (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)) %>%
        ## ungroup() %>%
        mutate(
            variable = "$q_t$ budget terms",
            ## ## preserve order
            ## term = factor(term, levels = 
            ##                         names(df.q)[which(names(df.q) == "delta.R") : 
            ##                                     which(names(df.q) == "E_q.delta.q")]), 
            ## make labels fancy
            term = revalue(term, 
                           c(delta.R = "$\\Delta R$",
                             ## Dq_t.0.Dt = "$\\DDP{\\hat q_t}{t}$",
                             rhoH.Dq_t.0.Dt = "$\\rho H(\\DDP{\\hat q_t}{t})$",
                             rhoH.dq_t.0.dt = "$\\rho H(\\ddp{\\hat q_t}{t})$",
                             rhoH.dqcore.0 = "$-\\rho H({\\bf v}\\cdot\\nabla\\hat q_t)$",
                             rhoH.adv.q_t.0 = "$\\rho H({\\bf v}\\cdot\\nabla\\hat q_t)$",
                             ## rhoH.dq_t.0.dt = "$\\rho H(\\ddp{\\hat q_t}{t})$",
                             qflx = "$\\text{LH}/L_v$",
                             ## delta.q = "$q_t^+ - \\hat q_t$",
                             E_q.delta.q = "$E_q(q_t^+ - \\hat q_t)$"#,
                             ## E_q = "$E_q$"))
                             )))

    if (!is.null(clip.quantiles)) {
        df.budget %<>%
            ## group_by(term) %>%
            filter(between(x, quantile(x, clip.quantiles, na.rm = TRUE), quantile(x, 1 - clip.quantiles, na.rm = TRUE)))
        ## ungroup() %>%
    }

    df.budget %>%
        ggplot(aes(x = term, y = x)) +
        ## geom_violin(aes(fill = factor(2 * ceiling((72 - jk.pbl) / 2))), draw_quantiles = 0.5) + ## , col = PNNLonyx, fill = PNNLcopper) +
        ## geom_violin(draw_quantiles = 0.5, col = PNNLonyx, fill = PNNLcopper) +
        geom + 
        ## geom_boxplot() +
        ##scale_y_continuous(labels = tikz_sanitize) +
        ## scale_y_continuous(labels = tikz_sanitize, limits = c(-0.15, 0.15) * 1e-3) +
        labs(x = "", y = "kg m$^{-2}$ s$^{-1}$") +
        ##ylim(-0.2, 0.2) +
        geom_hline(yintercept = 0, lty = "dashed", col = "lightgrey") +
        facet_grid(variable ~ .) +
        theme_bw() 
}
