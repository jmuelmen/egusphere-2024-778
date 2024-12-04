## @knitr giss-ppe-1e-lwp
df.giss.ppe.susc %>%
    ggplot(aes(x = factor(""), y = S.L, fill = factor(aer), lty = experiment)) +
    scale_fill_brewer("$N_a$", palette = "Spectral", direction = -1, drop = FALSE) +
    geom_violin(draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    labs(x = "", y = "$d\\log\\lwp/d\\log\\nd$")

## @knitr giss-ppe-1e-E
df.giss.ppe.susc %>%
    gather(Etype, S.E, S.Etheta, S.Eq) %>%
    mutate(Etype = revalue(Etype, c(S.Etheta = "$E_\\theta$", S.Eq = "$E_q$"))) %>%
    ggplot(aes(x = "$E$", y = S.E, fill = factor(aer), lty = experiment)) +
    scale_fill_brewer("$N_a$", palette = "Spectral", direction = -1, drop = FALSE) +
    scale_linetype("") +
    geom_violin(draw_quantiles = 0.5) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    facet_grid(Etype ~ .) +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    labs(x = "", y = "$d\\log E/d\\log\\nd$")

## @knitr giss-ppe-1e-R
df.giss.ppe.susc %>%
    ggplot(aes(x = factor(""), y = S.R, fill = factor(aer), lty = experiment)) +
    scale_fill_brewer("$N_a$", palette = "Spectral", direction = -1, drop = FALSE) +
    geom_violin(draw_quantiles = 0.5) +
    geom_hline(yintercept = 0, lty = "dashed", col = "darkgrey") +
    labs(x = "", y = "$d\\log R/d\\log\\nd$")
