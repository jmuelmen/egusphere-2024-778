#'
#' @export
entrainment.scatterplot <- function(df, heatmap = TRUE, marginal = FALSE,
                                    labels = labs(x = "$E_\\theta$ (kg~m$^{-2}$~s$^{-1}$)",
                                                  y = "$E_q$ (kg~m$^{-2}$~s$^{-1}$)"),
                                    ...) {
    g <- df %>%
        ggplot(aes(E_theta, E_q, ...)) 


    if (heatmap) {
        g <- g + stat_bin2d(bins = 200, geom = "raster") + scale_fill_distiller(palette = "Spectral")
    }

    g <- g +
        geom_smooth() +
        geom_abline(slope = 1, col = "grey", lty = "dashed") +
        scale_x_continuous(n.breaks = 3,
                           expand = expansion()) +
        scale_y_continuous(n.breaks = 3,
                           expand = expansion()) +
        labels +
        coord_cartesian(xlim = c(-0.05, 0.0499), ylim = c(-0.0499, 0.05)) +
        guides(fill = "none")
    
    if (marginal) {
        g.q <- df %>%
            ggplot(aes(x = E_q, ...)) +
            stat_ecdf(n = 200, pad = FALSE, geom = "line") +
            scale_x_continuous(n.breaks = 3,
                               expand = expansion()) +
            scale_y_continuous(labels = NULL,
                               expand = expansion()) +
            coord_flip(xlim = c(-0.0499, 0.05)) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                  axis.text.y = element_blank()) +
            theme(legend.position="none")
        g.theta <- df %>%
            ggplot(aes(x = E_theta, ...)) +
            stat_ecdf(n = 200, pad = FALSE, geom = "line") +
            scale_x_continuous(n.breaks = 3,
                               expand = expansion()) +
            scale_y_continuous(labels = NULL,
                               expand = expansion()) +
            coord_cartesian(xlim = c(-0.0499, 0.05)) +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                  axis.text.x = element_blank()) +
            theme(legend.position="none")
        wrap_plots(list(marginal.theta = g.theta, plot_spacer(), correlation = g, marginal.q = g.q),
                   widths = c(1, 1/3), heights = c(1/3, 1),
                   ncol = 2, nrow = 2, byrow = TRUE)
    } else {
        g
    }
    
}
