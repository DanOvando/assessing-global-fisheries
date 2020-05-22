plot_prior_fit <- function(metric,fit) {

  fit_r2 <- bayes_R2(fit)

  dat <- fit$data

  ppc_plot <- bayesplot::ppc_scatter_avg(y = dat$log_value, yrep = rstanarm::posterior_predict(fit)) +
    labs(
      x = paste0("Mean Posterior Predicted log(",metric,")"),
      y = paste0("Observed log(",metric,")")
    ) +
    geom_smooth(method = "lm", se = TRUE)

  br2_plot <- ggplot() +
    geom_density(data = data_frame(r2 = fit_r2), aes(r2), fill = "lightgrey") +
    labs(x = bquote(R^2), y = "Count")

  ppc_plot + br2_plot + plot_layout(ncol = 2, widths = c(3, 1)) +
    plot_annotation(tag_levels = "A")


}