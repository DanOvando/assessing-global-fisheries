#' assess fit of sraplus
#'
#' compares the fit of sraplus to some known values
#'
#' @param observed the observed values
#' @param estimated the estimated values
#' @param interval the credible interval you want to plot
#' @param window the number of years to calculate summary statistocs
#' @param plots TRUE or FALSE to include plots
#'
#' @return a list of performance fits
#' @export
#'
assess_fit <-
  function(observed,
           estimated,
           interval = 0.9,
           window = 5,
           plots = FALSE) {
    # observed <- ram_fits$data[[1]]
    #
    # estimated <- ram_fits$fit[[1]]
    #
    # ram_fits$fit[[1]] %>% plot_fit()
    if (class(estimated) == "srafit") {
      estimated = list(estimated)
    }


    if (is.null(names(estimated))) {
      names(estimated) <- paste0("fit-", seq_along(estimated))

    }

    fit_names <- names(estimated)

    status_fits <-
      map2_df(estimated, names(estimated), tidy_fits, years = observed$year)

    sraplus_worked <- TRUE

    flat_observed <- observed %>%
      select(year, b_v_bmsy, u_v_umsy, b_rel) %>%
      rename(depletion = b_rel) %>%
      gather(metric, value,-year) %>%
      group_by(metric)

    status_hat <- status_fits %>%
      gather(metric, value,-year,-draw,-fit_name) %>%
      group_by(year, metric, fit_name) %>%
      mutate(percentile = percent_rank(value)) %>%
      filter(percentile >= (1 - interval) / 2, percentile <=  1 - (1 - interval) / 2) %>%
      summarise(
        lower = quantile(value, (1 - interval) / 2, na.rm = TRUE),
        upper = quantile(value, 1 - (1 - interval) / 2, na.rm = TRUE),
        mean_value = mean(value, na.rm = TRUE),
        median_value = median(value, na.rm = TRUE)
      ) %>%
      ungroup()

    if (plots == TRUE) {
      catch_plot <- observed %>%
        ggplot(aes(year, catch)) +
        geom_point() +
        theme(axis.title.x = element_blank()) +
        scale_y_continuous(labels = scales::comma)

      trend_plot <- status_hat %>%
        ggplot() +
        geom_hline(aes(yintercept = 1), linetype = 2, color = "red") +
        geom_ribbon(aes(
          x = year,
          ymin = lower,
          ymax = upper,
          fill = fit_name
        ), alpha = 0.5) +
        geom_line(aes(x = year, y = median_value, color = fit_name), size = 1.5) +
        geom_point(
          data = flat_observed,
          aes(x = year, y = value),
          color = "blue",
          shape = 16,
          alpha = 0.5,
          size = 2
        ) +
        facet_wrap(~ metric, scales = "free_y")

      if (nrow(status_hat) > 0) {
        summary_plot <- catch_plot + {
          trend_plot
        } +
          plot_layout(nrow = 2,
                      ncol = 1,
                      heights = c(1, 2)) +
          plot_annotation(title = observed$scientificname,
                          subtitle = observed$stocklong.x)
      } else{
        summary_plot <- catch_plot +
          geom_text(
            data = observed,
            aes(
              x = 1980,
              y = mean(catch),
              label = "SRAPLUS FAILED"
            ),
            color = "red",
            size = 10
          )

        sraplus_worked <- FALSE

      }


      kobe_plot <-  status_fits %>%
        filter(year > (max(year) - 5)) %>%
        ggplot(aes(b_v_bmsy, u_v_umsy)) +
        geom_bin2d() +
        facet_grid(fit_name ~ year) +
        geom_hline(aes(yintercept = 1), linetype = 2) +
        geom_vline(aes(xintercept = 1), linetype = 2) +
        geom_point(
          data = observed %>% filter(year > (max(year) - 5)),
          aes(b_v_bmsy, u_v_umsy),
          size = 4,
          color = "red",
          alpha = 0.75
        ) +
        scale_fill_viridis(option = "D",
                           guide = guide_colorbar(frame.colour = "black",
                                                  barheight = 10)) +
        labs(title = observed$scientificname,
             subtitle = observed$stocklong.x)

    } # close plots


    fit <-  status_fits %>%
      gather(metric, value,-year,-draw,-fit_name) %>%
      rename(value_hat = value) %>%
      left_join(flat_observed, by = c("year", "metric")) %>%
      filter(year >= (max(year) - window)) %>%
      ungroup() %>%
      mutate(error = value - value_hat)

    accuracy <- fit %>%
      filter(metric == "b_v_bmsy") %>%
      nest(-fit_name, -metric) %>%
      mutate(accuracy = map_dbl(data, ~ check_classification(
        bhat = .x$value_hat, b = .x$value
      ))) %>%
      select(-data) %>%
      unnest()

    fit_summary <- fit  %>%
      group_by(metric, fit_name) %>%
      summarise(
        observed = mean(value),
        predicted = mean(value_hat),
        rmse = sqrt(mean(error ^ 2)),
        rmedse = sqrt(median(error ^ 2)),
        bias = median(error)
      ) %>%
      ungroup() %>%
      left_join(accuracy, by = c("metric", "fit_name"))


    rm(observed, estimated)

    if (plots == FALSE) {
      kobe_plot <- NA

      summary_plot <- NA

    }

    out <-
      list(
        kobe_plot = kobe_plot,
        summary_plot = summary_plot,
        sraplus_worked = sraplus_worked,
        fit_summary = fit_summary
      )

  }