

# setup -------------------------------------------------------------------

library(tidyverse)
library(ggridges)
library(gganimate)
library(rstan)
library(mvtnorm)
library(FishLife)
library(furrr)
library(future)
library(viridis)
library(sraplus)
library(patchwork)
library(scales)
library(hrbrthemes)
library(TMBhelper)
library(TMB)
library(here)
library(tmbstan)
library(tidybayes)
library(rstan)
library(rstanarm)
library(sf)
Sys.unsetenv("PKG_CXXFLAGS")

options(dplyr.summarise.inform = FALSE)

rstan::rstan_options(auto_write = TRUE)

sraplus::get_tmb_model()
# options -----------------------------------------------------------------

min_years_catch <- 20

crazy_b <- 4 # threshold for suspect B/Bmsy value

crazy_u <- 5 # threshold for suspect U/Umsy value

draws <- 2500

min_draws <- 2000 # minimum number of unique SIR draws

n_cores <- 3
# number of cores for parallel processing

future::plan(multiprocess, workers = n_cores)


# options(mc.cores = 1)

results_name <- "v0.5"


results_description <-
  "publication version of results updated VOI from version 0.5 "

run_voi_models <- FALSE
# sub options for run_voi_models
fit_models <- FALSE

write_results <- TRUE

process_fits <- FALSE

# other things

run_continent_examples <- FALSE

run_ei_example <- FALSE

run_sofia_comparison <- TRUE

run_ram_tests <- TRUE

run_ram_comparison <- TRUE



engine = "stan"

catchability = 1e-2

theme_set(theme_classic() + theme(strip.background = element_rect(color = "transparent")))
# data(Return)

# return is from here, WARNING, changes rapidly, things break check and make sure this isn't why
# https://drive.google.com/drive/u/0/folders/1J46tM6PYDdPwhx5zGrlHMdxUyGRrky7X?ogsrc=32
# load(here::here("data","Return.Rdata"))

functions <- list.files(here::here("functions"))

purrr::walk(functions, ~ source(here::here("functions", .x)))

results_path <- here::here("results", results_name)


experiment_path <-
  here::here("results", results_name, "experiments")


if (dir.exists(here::here("results", results_name)) == FALSE) {
  results_path <- here::here("results", results_name)
  
  dir.create(results_path, recursive = TRUE)
  
  write(results_description,
        file = here::here("results", results_name, "description.txt"))
  
  dir.create(experiment_path, recursive = TRUE)
  
  
}


prepare_sofia_data(lookup_fmi_names = FALSE)


has_total_biomass = ram_data %>%
  group_by(stockid) %>%
  summarise(hastb = !any(is.na(total_biomass)),
            hasb = !any(is.na(b_v_bmsy))) %>%
  filter(hastb & hasb)

ram_data <- ram_data %>%
  filter(stockid %in% has_total_biomass$stockid) %>%
  mutate(has_things = !(is.na(u_v_umsy) |
                          is.na(total_biomass) |
                          is.na(catch) | catch == 0)) %>%
  filter(has_things) %>%
  group_by(stockid) %>%
  mutate(delta_year = as.integer(year - lag(year))) %>%
  mutate(delta_year = case_when(year == min(year) ~ as.integer(1),
                                TRUE ~ delta_year)) %>%
  mutate(missing_gaps = any(delta_year > 1)) %>%
  filter(missing_gaps == FALSE) %>%
  group_by(stockid) %>%
  mutate(n = length(catch)) %>%
  filter(n >= min_years_catch) %>%
  ungroup()

# fit voi models --------------------------------------------------------------


if (run_voi_models == TRUE) {
  stress_tests <- data_frame(
    use_index = sample(c(TRUE, FALSE),  draws, replace = TRUE),
    index_type = sample(c("survey"), draws, replace = TRUE),
    initial_state_type = sample(c(
      "unfished", "known", "prior", "heuristic"
    ), draws, replace = TRUE),
    initial_state_cv =  runif(draws, 0.05, 0.2),
    use_terminal_state = sample(c(TRUE, FALSE),  draws, replace = TRUE),
    terminal_state_cv = runif(draws, 0.05, 0.2),
    error_cv = sample(c(0, 1), draws, replace = TRUE),
    stockid = sample(unique(ram_data$stockid), draws, replace = TRUE),
    use_u_priors =  sample(c(TRUE, FALSE),  draws, replace = TRUE),
    u_cv = runif(draws, 0.05, 0.2),
    index_window = sample(c(1, 0.25, .5), draws, replace = TRUE),
    index_freq = sample(c(1, 2), draws, replace = TRUE),
    u_window = sample(c("snapshot", "recent", "complete"), draws, replace = TRUE),
    b_ref_type = sample(c("k", "b"),  draws, replace = TRUE),
    f_ref_type = sample(c("f", "fmsy"),  draws, replace = TRUE),
    estimate_shape = sample(c(FALSE, TRUE), draws, replace = TRUE),
    estimate_proc_error = sample(c(FALSE, TRUE), draws, replace = TRUE)
  ) %>%
    mutate(error_cv = map_dbl(error_cv, ~ ifelse(.x == 0, 1e-6, runif(1, 0.05, 0.2)))) %>%
    mutate(b_ref_type = ifelse(initial_state_type == "prior",
                               "b",
                               b_ref_type))
  
  # stress_tests <- data_frame(
  #   use_index = sample(c(FALSE),  draws, replace = TRUE),
  #   index_type = sample(c("survey"), draws, replace = TRUE),
  #   initial_state_type = sample(c("unfished", "known", "prior","heuristic"), draws, replace = TRUE),
  #   initial_state_cv =  runif(draws, 0.05, 0.2),
  #   use_terminal_state = sample(c(TRUE),  draws, replace = TRUE),
  #   terminal_state_cv = runif(draws, 0.05, 0.2),
  #   error_cv = sample(c(0, 1), draws, replace = TRUE),
  #   stockid = sample(unique(ram_data$stockid), draws, replace = TRUE),
  #   use_u_priors =  sample(c(FALSE),  draws, replace = TRUE),
  #   u_cv = runif(draws, 0.05, 0.2),
  #   index_window = sample(c(1, 0.25, .5), draws, replace = TRUE),
  #   index_freq = sample(c(1, 2), draws, replace = TRUE),
  #   u_window = sample(c("snapshot", "recent", "complete"), draws, replace = TRUE),
  #   b_ref_type = sample(c("b"),  draws, replace = TRUE),
  #   f_ref_type = sample(c("f", "fmsy"),  draws, replace = TRUE),
  #   estimate_shape = sample(c(FALSE, TRUE), draws, replace = TRUE),
  #   estimate_proc_error = sample(c(FALSE, TRUE), draws, replace = TRUE)
  # ) %>%
  #   mutate(error_cv = map_dbl(error_cv, ~ ifelse(.x == 0, 1e-6, runif(1, 0.05, 0.2)))) %>%
  #   mutate(
  #     b_ref_type = ifelse(
  #       initial_state_type == "prior",
  #       "b",
  #       b_ref_type
  #     )
  #   )
  
  
  if (fit_models == TRUE) {
    ram_fits <- stress_tests %>%
      left_join(ram_data %>% group_by(stockid) %>%  nest() %>% ungroup(), by = "stockid") %>%
      # filter(f_ref_type == "f",
      #        use_u_priors == TRUE,
      #        u_window == 1,
      #        use_terminal_u == 0) %>%
      mutate(priors = pmap(
        list(
          dat = data,
          use_index = use_index,
          index_type = index_type,
          initial_state_type = initial_state_type,
          initial_state_cv =  initial_state_cv,
          use_terminal_state = use_terminal_state,
          terminal_state_cv = terminal_state_cv,
          error_cv = error_cv,
          use_u_priors = use_u_priors,
          u_cv = u_cv,
          f_ref_type = f_ref_type,
          index_window = index_window,
          index_freq = index_freq,
          u_window = u_window
        ),
        generate_priors,
        q = catchability
      )) %>%
      ungroup() %>%
      mutate(experiment = 1:nrow(.))
    
    
    saveRDS(ram_fits,
            file = here::here("results", results_name, "ram_experiments.rds"))
    
    # sra_model <-
    #   rstan::stan_model(file = here::here("src", "biosra.stan"))
    #
    ram_fits <- ram_fits %>%
      # filter(initial_state_type == "unfished",
      #        u_window == "snapshot",
      #        use_index == TRUE,
      #        use_u_priors == TRUE,
      #        estimate_shape == FALSE) %>%
      # filter(use_u_priors == FALSE,
      #        use_index == TRUE,
      #        index_window == 1) %>%
      # mutate(sciname = map_chr(data,~unique(.x$"scientificname"))) %>%
      # filter(sciname == "Ophiodon elongatus") %>%
      # sample_n(20) %>%
    # slice(6) %>%
    mutate(
      fits = future_pmap(
        list(
          dat = data,
          priors = priors,
          experiment = experiment,
          estimate_shape = estimate_shape,
          estimate_proc_error = estimate_proc_error,
          estimate_initial_state = !initial_state_type == "unfished"
        ),
        safely(run_assessments),
        assessments = c("cmsy", "sraplus"),
        model = "sraplus_tmb",
        results_name = results_name,
        write_results = write_results,
        engine = engine,
        keep = 2000,
        .progress = TRUE
      )
    )
    
    saveRDS(ram_fits,
            file = here::here("results", results_name, "ram_fits.Rds"))
  } else {
    # ram_fits <-
    #   readRDS(file = here::here("results", results_name, "ram_fits.Rds"))
    #
    ram_fits <-
      readRDS(file = here::here("results", results_name, "ram_experiments.rds"))
    
  }
  
  
  # wtf <- map(ram_fits$fits,'error') %>% map_lgl(is.null)
  #
  # a <- ram_fits %>%
  #   filter(!wtf)
  
  # process fits ------------------------------------------------------------
  
  if (process_fits == TRUE) {
    fits <-
      list.files(here::here("results", results_name, "experiments"))
    
    load_foo <-
      function(experiment,
               data,
               results_name,
               produce = "summary",
               plots = FALSE) {
        if (produce == "fit") {
          out <- readRDS(here::here(
            "results",
            results_name,
            "experiments",
            paste0("experiment-", experiment, ".rds")
          ))
        }
        if (produce == "summary") {
          fit <- readRDS(here::here(
            "results",
            results_name,
            "experiments",
            paste0("experiment-", experiment, ".rds")
          ))
          
          out <-
            assess_fit(
              data,
              fit,
              interval = 0.9,
              window = 5,
              plots = plots
            )
          
          
        }
        
        return(out)
      }
    
    # future::plan(multiprocess, workers = n_cores)
    
    ram_fits <- ram_fits %>%
      mutate(
        fits = map2(
          experiment,
          data,
          safely(load_foo),
          results_name = results_name,
          produce = "summary",
          plots = FALSE
        )
      )
    
    
    saveRDS(ram_fits,
            file = here::here("results", results_name, "ram_results.Rds"))
  } else {
    ram_fits <-
      readRDS(here::here("results", results_name, "ram_results.Rds"))
  }
  
  # assess performance ------------------------------------------------------
  
  
  fit_failed <-
    !(ram_fits$fits %>% map("error") %>% map_lgl(is.null))
  
  ram_fits <- ram_fits %>%
    ungroup() %>%
    filter(map(ram_fits$fits, "error") %>% map_lgl(is.null)) %>%
    mutate(fit = map(fits, "result")) %>%
    select(-fits)
  
  ram_error <- ram_fits %>%
    select(experiment, data, priors, index_type) %>%
    mutate(error = pmap(
      list(
        data = data,
        priors = priors,
        index_type = index_type
      ),
      calc_error
    )) %>%
    select(-data,-priors,-index_type) %>%
    unnest(cols = "error")
  
  ram_fits <- ram_fits %>%
    mutate(
      summary_plot = map(fit, "summary_plot"),
      summary = map(fit, "fit_summary")
    )
  # test <- test %>%
  #   mutate(
  #     summary_plot = map_plot(fit_performance, "summary_plot"),
  #     summary = map(fit_performance, "fit_summary")
  #   )
  
  
  
  # trelliscopejs::trelliscope(
  #   ram_fits %>%
  #     filter(use_index == TRUE),
  #   name = "blah",
  #   panel_col = "summary_plot"
  # )
  
  
  fits <- ram_fits %>%
    select(-data,-summary_plot,-fit,-summary_plot,-priors) %>%
    unnest(cols = summary) %>%
    filter(!is.na(rmse), is.finite(rmse)) %>%
    left_join(ram_error, by = "experiment") %>%
    group_by(stockid, metric) %>%
    mutate(baseline_rmse = mean(rmse)) %>%
    mutate(delta_rmse = rmse - baseline_rmse) %>%
    ungroup()
  
  # why infinite RMSE
  
  
  fits$initial_state_type <-
    forcats::fct_relevel(fits$initial_state_type, "heuristic")
  
  fits$u_window[fits$use_u_priors == FALSE] <- "none"
  
  fits$u_window <- forcats::fct_relevel(fits$u_window, "none")
  
  
  wtf <- fits %>%
    filter(
      use_terminal_state == TRUE,
      metric == "b_v_bmsy",
      b_ref_type == "b",
      error_cv == min(error_cv),
      use_index == FALSE
    ) %>%
    arrange(desc(rmse))
  
  View(wtf)
  
  i <- 2
  
  huh <-
    read_rds(file.path(
      results_path,
      "experiments",
      paste0("experiment-", wtf$experiment[i], ".rds")
    ))
  
  driors <-
    ram_fits$priors[ram_fits$experiment == wtf$experiment[i]][[1]]
  
  data <-
    ram_fits$data[ram_fits$experiment == wtf$experiment[i]][[1]]
  
  plot_sraplus(huh$sirplus)
  
  plot_prior_posterior(huh$sraplus, driors)
  
  
  # plot diagnostics --------------------------------------------------------
  
  noerror_data <- fits %>%
    filter(error_cv == min(error_cv),
           fit_name != 'cmsy')
  
  
  voi_data <- fits %>%
    filter(fit_name != "cmsy")
  
  
  voi_data %>%
    filter(use_index == FALSE,
           metric == "b_v_bmsy",
           b_ref_type == "b") %>%
    ggplot(aes(rmse, fill = use_terminal_state)) +
    geom_density(position = "dodge")
  
  
  
  noerror_data %>%
    filter(use_index == FALSE,
           metric == "b_v_bmsy",
           b_ref_type == "b") %>%
    ggplot(aes(observed, predicted, color = use_terminal_state)) +
    geom_point()
  
  obs_v_pred_plot <- noerror_data %>%
    filter(metric == "b_v_bmsy") %>%
    ggplot(aes(observed, predicted, color = fit_name)) +
    geom_point(alpha = 0.25) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    geom_smooth(method = 'lm', aes(fill = fit_name)) +
    facet_wrap(use_index ~ metric)
  
  obs_v_pred_plot <- noerror_data %>%
    filter(metric == "depletion") %>%
    ggplot(aes(observed, predicted, color = fit_name)) +
    geom_point(alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    geom_smooth(method = 'lm', aes(fill = fit_name)) +
    facet_wrap(use_index ~ metric)
  
  
  
  overall_rmse_plot <- noerror_data %>%
    filter(is.finite(rmse)) %>%
    ggplot(aes(rmse, metric, fill = fit_name)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 2))
  
  
  overall_rmse_plot <- noerror_data %>%
    filter(is.finite(rmse)) %>%
    ggplot(aes(rmse, metric, fill = u_window)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 2)) +
    facet_wrap( ~ use_index)
  
  
  overall_bias_plot <- noerror_data %>%
    ggplot(aes(bias, metric, fill = fit_name)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(-1, 1))
  
  
  u_effect_plot <- noerror_data %>%
    ggplot(aes(rmse, use_u_priors)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0, 1)) +
    facet_wrap( ~ metric)
  
  relative_u_effect_plot <- noerror_data %>%
    ggplot(aes(delta_rmse, use_u_priors)) +
    geom_vline(aes(xintercept = 0), color = "red", linetype = 2) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(-1, 1)) +
    facet_wrap( ~ metric)
  
  
  rmse_plot <- noerror_data %>%
    ggplot(aes(rmse, fit_name, fill = 0.5 - abs(0.5 - ..ecdf..))) +
    stat_density_ridges(alpha = 0.75,
                        geom = "density_ridges_gradient",
                        calc_ecdf = TRUE) +
    scale_fill_viridis(name = "Tail probability", direction = -1) +
    scale_x_continuous(limit = c(0, 1)) +
    facet_wrap(~ metric)
  
  
  # conditional on the same data, how well does it perform?
  
  # cmsy_b_voi_fit <-
  #   rstanarm::stan_glm(
  #     log(rmse) ~  use_u_priors + use_initial_state + stockid,
  #     data = noerror_data %>% filter(metric == "b_v_bmsy", fit_name == "cmsy"),
  #     cores = 4
  #   )
  #
  # #
  # cmsy_b_voi_plot <- bayesplot::mcmc_areas(
  #   as.array(cmsy_b_voi_fit),
  #   regex_pars = c("index", "effort", "terminal", "initial", "fit_"),
  #   prob = 0.8,
  #   # 80% intervals
  #   prob_outer = 0.9,
  #   # 99%
  #   point_est = "mean"
  # ) +
  #   scale_x_percent(name = "Change in RMSE")
  
  # value of information analyses
  
  # how does the comparison stack up when you have more data?
  rmse_bias_plot <- noerror_data %>%
    ggplot(aes(rmse, bias)) +
    geom_hex(aes(fill = ..density..), binwidth = c(0.1, 0.1)) +
    facet_grid(metric ~ fit_name) +
    scale_x_continuous(limits = c(0, .5)) +
    scale_y_continuous(limits = c(-.5, .5))  +
    scale_fill_viridis(option = "A")
  
  b_rmse_fit <-
    rstanarm::stan_glmer(
      log(rmse) ~ fit_name + (1 | stockid),
      data = noerror_data %>% filter(metric == "b_v_bmsy"),
      cores = 4
    )
  
  b_rmse_fit_posterior <- as.array(b_rmse_fit)
  
  b_rmse_plot <- bayesplot::mcmc_areas(
    b_rmse_fit_posterior,
    regex_pars = "fit_name",
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.99,
    # 99%
    point_est = "mean"
  ) +
    scale_x_percent(name = "Percent change in B/Bmsy RMSE")
  
  
  
  # value of information
  
  b_voi_fit <-
    rstanarm::stan_glmer(
      rmse ~  initial_state_type  + use_index * u_window + estimate_shape + estimate_proc_error + use_terminal_state + (1 |
                                                                                                                          stockid),
      data = voi_data %>% filter(metric == "b_v_bmsy"),
      cores = 4,
      family = Gamma(link = "log")
    )
  
  
  b_voi_plot <- bayesplot::mcmc_areas(
    as.array(b_voi_fit),
    regex_pars = c("use", "state", "estimate", "index", "u_"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x)
      exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in RMSE")
  
  u_voi_fit <-
    rstanarm::stan_glmer(
      rmse ~  initial_state_type  + use_index * u_window + estimate_shape + estimate_proc_error + use_terminal_state + (1 |
                                                                                                                          stockid),
      data = voi_data %>% filter(metric == "u_v_umsy", is.finite(rmse)),
      cores = 4,
      family = Gamma(link = "log")
    )
  
  u_voi_plot <- bayesplot::mcmc_areas(
    as.array(u_voi_fit),
    regex_pars = c("use", "state", "estimate", "index", "u_"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x)
      exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in RMSE")
  
  
  
  index_voi_fit <-
    rstanarm::stan_glmer(
      rmse ~  factor(index_freq) + factor(index_window) + (1 | stockid),
      data = voi_data %>% filter(metric == "b_v_bmsy", use_index == TRUE),
      family = Gamma(link = "log"),
      cores = 4
    )
  
  
  index_voi_plot <- bayesplot::mcmc_areas(
    as.array(index_voi_fit),
    regex_pars = c("freq", "window"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x)
      exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in RMSE")
  
  
  
  # effect on accuracy
  
  acc_data <- voi_data %>%
    filter(metric == "b_v_bmsy") %>%
    mutate(bin_acc = accuracy > 0.5)
  
  
  acc_voi_fit <-
    rstanarm::stan_glmer(
      bin_acc ~  initial_state_type  + use_index * u_window + estimate_shape + estimate_proc_error + use_terminal_state + (1 |
                                                                                                                             stockid),
      data = acc_data,
      cores = 4,
      family = binomial()
    )
  
  
  acc_voi_plot <- bayesplot::mcmc_areas(
    as.array(acc_voi_fit),
    regex_pars = c("use", "state", "estimate", "index", "u_"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x)
      exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in Accuracy")
  
  acc_index_voi_fit <-
    rstanarm::stan_glmer(
      bin_acc ~  factor(index_freq) + factor(index_window) + (1 | stockid),
      data = acc_data %>% filter(use_index == TRUE),
      cores = 4,
      family = binomial()
    )
  
  
  acc_index_voi_plot <- bayesplot::mcmc_areas(
    as.array(acc_index_voi_fit),
    regex_pars = c("use", "state", "estimate", "index"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x)
      exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in Accuracy")
  
  # how wrong can you be?
  
  quality_data <- fits %>%
    filter(use_index == TRUE,
           fit_name == "sraplus",
           error_cv > min(error_cv))
  
  quality_obs_v_pred_plot <- quality_data %>%
    filter(metric == "b_v_bmsy") %>%
    ggplot(aes(observed, predicted, color = pmin(20, index_rmse))) +
    geom_point(alpha = 0.75, size = 4) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    geom_smooth(method = 'lm') +
    facet_wrap( ~ metric) +
    scale_color_viridis()
  
  quality_effect_plot <- quality_data %>%
    filter(metric == "b_v_bmsy") %>%
    ggplot(aes(index_rmse, rmse)) +
    geom_point() +
    geom_abline(aes(slope = 0, intercept = mean(rmse))) +
    geom_smooth(method = "lm") +
    scale_x_log10() +
    scale_y_log10()
  
  flist <- ls()[str_detect(ls(), "_fit")]
  
  save(list = flist,
       file = file.path(results_path, "voi_fits.RData"))
  
  # voi_data %>%
  #   ggplot(aes((index_rmse), (rmse))) +
  #   geom_point() +
  #   scale_y_continuous(limits = c(0, 1)) +
  #   scale_x_continuous(limits = c(0, 1)) +
  #   geom_smooth(method = "lm") +
  #   facet_grid(index_window ~ index_freq)
  #
  # voi_data %>%
  #   ggplot(aes(log(u_rmse), log(rmse))) +
  #   geom_point() +
  #   scale_y_continuous(limits = c(0, 1)) +
  #   geom_smooth(method = "lm")
  #
  # voi_data %>%
  #   ggplot(aes(log(terminal_b_rmse), log(rmse))) +
  #   geom_point() +
  #   scale_y_continuous(limits = c(0, 1)) +
  #   geom_smooth(method = "lm")
  #
  #
  # voi_data %>%
  #   ggplot(aes(rmse, y = index_window, group = index_window)) +
  #   geom_density_ridges() +
  #   scale_x_continuous(limits = c(0, 1))
  #
  #
  # voi_data %>%
  #   ggplot(aes(rmse, y = index_freq, group = index_freq)) +
  #   geom_density_ridges() +
  #   scale_x_continuous(limits = c(0, 1))
  #
  #
  # sofiaplus_plots <- ls()[str_detect(ls(), "_plot")]
  #
  # sofiaplus_fits <- ls()[str_detect(ls(), "_fit")]
  #
  # saveRDS(fits,
  #         file = here::here("results", results_name, "sofiaplus_fits.rds"))
  #
  # saveRDS(ram_fits,
  #         file = here::here("results", results_name, "sofiaplus_rawfits.rds"))
  #
  #
  # save(
  #   list = sofiaplus_plots,
  #   file = here::here("results", results_name, "sofiaplus_plots.Rdata")
  # )
  #
  # save(
  #   list = sofiaplus_fits,
  #   file = here::here("results", results_name, "sofiaplus_fits.Rdata")
  # )
  
  
} # close fit voi models

# run asia-eu comparison ------------------------------------------------------


support_data <- list(effort_data = effort_data,
                     effort_region_to_country = effort_region_to_country)



total_nominal_effort <- rous_data %>%
  left_join(effort_region_to_country, by = "country") %>%
  group_by(region, year) %>%
  summarise(total_effort = sum(effort_cell_reported_nom, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(region = tolower(region))


# total_effecive_effort <- effort_data %>%
#   filter(effort_type == "effective") %>%
#   group_by(region,year) %>%
#   summarise(total_effort = sum(effort)) %>%
#   ungroup() %>%
#   mutate(region = tolower(region))


total_nominal_effort %>%
  ggplot(aes(year, total_effort, color = region)) +
  geom_line()

sfe <- safely(fit_continent_examples)

continent_countries <- effort_region_to_country %>%
  filter(region %in% c("Northeast Asia", "Europe"))

total_stocks <- fao %>%
  filter(country %in% continent_countries$country,
         continent %in% c("Asia", "Europe")) %>%
  mutate(continent = case_when(continent == "Asia" ~ "Northeast Asia",
                               TRUE ~ continent)) %>%
  group_by(year, scientific_name, continent) %>%
  summarise(total_catch = sum(capture, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(scientific_name, continent) %>%
  mutate(
    lifetime_catch = sum(total_catch),
    min_catch = min(total_catch),
    catch_length = length(total_catch)
  ) %>%
  ungroup() %>%
  filter(lifetime_catch > 1e5, min_catch > 1e3, catch_length >= 20) %>%
  mutate(continent = tolower(continent)) %>%
  left_join(total_nominal_effort, by = c("year", "continent" = "region"))



total_stocks %>%
  group_by(continent, year) %>%
  summarise(cpue = mean(total_catch / total_effort, na.rm = TRUE)) %>%
  ggplot(aes(year, cpue, color = continent)) +
  geom_line()

total_stocks %>%
  group_by(continent, year, scientific_name) %>%
  summarise(catch = sum(total_catch)) %>%
  group_by(continent, scientific_name) %>%
  mutate(catch = scale(catch)) %>%
  ggplot(aes(year, catch)) +
  geom_smooth(show.legend = FALSE) +
  facet_wrap( ~ continent)


# future::plan(future::multiprocess, workers = 4)

if (run_continent_examples == TRUE) {
  set.seed(42)
  
  continent_fits <- total_stocks %>%
    group_by(scientific_name, continent) %>%
    nest() %>%
    ungroup() #%>%
  # sample_n(10)
  
  sfe <- safely(fit_continent_examples)
  
  continent_fits = continent_fits %>%
    # filter(scientific_name == "Lateolabrax japonicus") %>%
    mutate(
      ex_fits = pmap(
        list(scientific_name = scientific_name,
             data = data),
        sfe,
        model = "sraplus_tmb",
        engine = "stan",
        estimate_proc_error = TRUE,
        estimate_qslope = FALSE,
        estimate_shape = FALSE,
        first_effort_year = 1960
      )
    )
  
  
  # future::plan(future::sequential, workers = 1)
  
  write_rds(continent_fits, path =  file.path(results_path, "continent-fits.rds"))
} else {
  continent_fits <-
    read_rds(path =  file.path(results_path, "continent-fits.rds"))
  
}

fit_worked <-
  continent_fits$ex_fits %>% map("error") %>% map_lgl(is.null)

total_stocks %>%
  group_by(year, continent) %>%
  summarise(`Total Catch` = sum(total_catch)) %>%
  left_join(total_nominal_effort, by = c("continent" = "region", "year")) %>%
  mutate(CPUE = `Total Catch` / total_effort) %>%
  rename(Effort = total_effort) %>%
  gather(metric, value, `Total Catch`:CPUE) %>%
  ggplot(aes(year, value, color = continent)) +
  geom_line(size = 2, alpha = 0.75) +
  facet_wrap( ~ metric, scales = "free_y") +
  labs(x = "Year", y = "") +
  scale_color_discrete(name = "Region") +
  theme_minimal() +
  scale_y_comma()


# total_stocks %>%
#   group_by(year, continent, scientific_name) %>%
#   summarise(`Total Catch` = sum(total_catch)) %>%
#   left_join(total_effecive_effort, by = c("continent" = "region", "year")) %>%
#   mutate(CPUE = `Total Catch` / total_effort) %>%
#   group_by(scientific_name, continent) %>%
#   mutate(CPUE = scale(CPUE),
#          lifetime_catch = sum(`Total Catch`)) %>%
#   rename(Effort = total_effort) %>%
#   gather(metric, value, `Total Catch`:CPUE) %>%
#   filter(metric == "CPUE") %>%
#   ggplot(aes(year, value, color = scientific_name, alpha = lifetime_catch)) +
#   geom_line(size = 1, show.legend = FALSE) +
#   facet_wrap(~continent, scales = "free_y") +
#   labs(x = "Year", y = "Centered and Scaled Effective Abundance Index") +
#   scale_color_discrete(name = "Region") +
#   theme_minimal() +
#   scale_y_comma()



a <- continent_fits %>%
  mutate(lifetime_catch = map_dbl(data,  ~ unique(.x$lifetime_catch))) %>%
  select(scientific_name, ex_fits, continent, lifetime_catch) %>%
  filter(fit_worked) %>%
  mutate(ex_fits = map(ex_fits, "result")) %>%
  unnest(cols = ex_fits) %>%
  filter(year == 2016)


continent_kobe_plot <- a %>%
  filter(variable %in% c("b_div_bmsy", "u_div_umsy")) %>%
  group_by(scientific_name, continent) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  mutate(year = as.factor(year)) %>%
  select(scientific_name,
         year,
         variable,
         mean,
         data,
         continent,
         lifetime_catch) %>%
  spread(variable, mean) %>%
  ggplot(aes(pmin(4, b_div_bmsy), pmin(4, u_div_umsy), color = continent)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  geom_density2d(alpha = 0.5) +
  # geom_hex(binwidth = c(0.25,0.25)) +
  geom_point(aes(size = lifetime_catch),
             show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 4), name = "B/Bmsy") +
  scale_y_continuous(limits = c(0, 4), name = "U/Umsy") +
  facet_grid(continent ~ data) +
  ggsci::scale_color_d3(name = "Region") +
  guides(size = FALSE)

continent_distribution_plot <- a %>%
  filter(variable %in% c("b_div_bmsy", "u_div_umsy")) %>%
  group_by(scientific_name, continent) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  mutate(year = as.factor(year)) %>%
  select(scientific_name,
         year,
         variable,
         mean,
         data,
         continent,
         lifetime_catch) %>%
  ggplot() +
  geom_density_ridges(
    aes(
      x = pmin(4, mean),
      y = continent,
      fill = continent
    ),
    show.legend = FALSE,
    alpha = 0.75
  ) +
  facet_grid(Region ~ variable) +
  scale_x_continuous(limits = c(0, NA))

continent_kobe_plot



# run east-india analysis -----------------------------------------------------

fao_names <- fao %>%
  select(common_name, scientific_name) %>%
  unique()


ei_data <- ei_capture %>%
  mutate(capture = dplyr::case_when(is.na(capture) ~ 0,
                                    TRUE ~ capture)) %>%
  group_by(species_asfis_species) %>%
  mutate(first_year = year[capture > 0][1],
         total_capture = sum(capture)) %>%
  filter(year >= first_year) %>%
  ungroup() %>%
  mutate(capture_rank = percent_rank(total_capture)) %>%
  filter(capture_rank > 0.25) %>%
  left_join(ei_cpue, by = "year") %>%
  left_join(fao_names, by = c("species_asfis_species" = "common_name")) %>%
  ungroup() %>%
  filter(!is.na(scientific_name)) %>%
  group_by(species_asfis_species, scientific_name) %>%
  nest() %>%
  ungroup()


if (run_ei_example == TRUE) {
  set.seed(42)
  
  # future::plan(future::multiprocess, workers = n_cores)
  
  sfi <- safely(fit_india)
  
  ei_fits <- ei_data %>%
    # sample_n(4) %>%
    # filter(scientific_name == "Sphyraena spp") %>%
    mutate(ex_fits = map2(
      scientific_name,
      data,
      ~ sfi(
        scientific_name = .x,
        catches = .y$capture,
        years = .y$year,
        index = .y$cpue[!is.na(.y$cpue) & .y$year >= 1990],
        index_years = .y$year[!is.na(.y$cpue) & .y$year >= 1990],
        support_data = support_data,
        country = "India",
        engine = "stan",
        sar = 18.7,
        chains = 1,
        cores = 1,
        n_keep = 2500
      )
    ))
  
  # future::plan(future::sequential, workers = 1)
  
  write_rds(ei_fits, path = file.path(results_path, "ei_fits.rds"))
} else {
  ei_fits <-   read_rds(path = file.path(results_path, "ei_fits.rds"))
  
}

fit_worked <- ei_fits$ex_fits %>% map("error") %>% map_lgl(is.null)

a <- ei_fits %>%
  select(scientific_name, ex_fits) %>%
  filter(fit_worked) %>%
  mutate(ex_fits = map(ex_fits, "result")) %>%
  unnest(cols = ex_fits) %>%
  mutate(data = ifelse(data == "cpue", "cpue + sar", data))

# a %>%
#   filter(scientific_name == "Sphyraena spp") %>%
#   filter(variable == "dep", data == "cpue + sar") %>%
#   ggplot(aes(seq_along(mean), mean)) +
#   geom_line()

# ei_data %>%
#   filter(scientific_name == "Sphyraena spp") %>%
#   unnest() %>%
#   ggplot(aes(year, capture)) +
#   geom_line()a

ei_kobe_plot <- a %>%
  filter(year == max(year)) %>%
  filter(variable %in% c("b_div_bmsy", "u_div_umsy")) %>%
  group_by(scientific_name) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  mutate(year = as.factor(year)) %>%
  select(scientific_name, year, variable, mean, data) %>%
  spread(variable, mean) %>%
  ggplot(aes(pmin(2, b_div_bmsy), pmin(4, u_div_umsy))) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  geom_density2d() +
  # geom_hex(binwidth = c(0.25,0.25)) +
  geom_point(
    show.legend = FALSE,
    size = 3,
    fill = "grey",
    alpha = 0.5
  ) +
  scale_x_continuous(limits = c(0, 2), name = "B/Bmsy") +
  scale_y_continuous(limits = c(0, 4), name = "U/Umsy") +
  facet_wrap( ~ data)

ei_kobe_plot



eu_cpue_plot <- ei_cpue %>%
  ggplot(aes(year, cpue)) +
  geom_line(size = 2) +
  labs(x = "Year", y = "CPUE (kg/hour)")



# run sofia-comparison --------------------------------------------------------

fao_country_region_key <- fao %>%
  group_by(country, fao_area_code) %>%
  summarise(tc = sum(capture, na.rm = TRUE)) %>%
  group_by(country) %>%
  mutate(ptc = percent_rank(tc)) %>%
  filter(ptc > 0.25) %>%
  ungroup()

fmi$country <-
  countrycode::countrycode(fmi$country_rfmo, "country.name", "un.name.en")



fmi_country_species <- fmi %>%
  filter(!is.na(country),!is.na(scientificname)) %>%
  select(country, scientificname) %>%
  unique() %>%
  rename(scientific_name = scientificname)


fao_country_species_regions <- fao %>%
  group_by(country, scientific_name, fao_area, fao_area_code) %>%
  summarise(tc = sum(capture, na.rm = TRUE)) %>%
  group_by(country) %>%
  mutate(ptc = percent_rank(tc)) %>%
  filter(ptc > 0.25) %>%
  ungroup()

fmi_country_species_region <- fmi_country_species %>%
  left_join(fao_country_species_regions,
            by = c("country", "scientific_name")) %>%
  filter(!is.na(fao_area))

mean_regional_isscaap_fmi <- fmi %>%
  left_join(
    fmi_country_species_region,
    by = c("country_rfmo" = "country", "scientificname" = "scientific_name")
  ) %>%
  select(
    lookup_code,
    fao_area_code,
    research,
    management,
    socioeconomics,
    enforcement,
    isscaap_group,
    tc
  ) %>%
  na.omit() %>%
  pivot_longer(
    cols = c("research", "management", "socioeconomics", "enforcement"),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(fao_area_code, isscaap_group, metric) %>%
  summarise(mean_value = weighted.mean(value, tc)) %>%
  pivot_wider(names_from = metric, values_from = mean_value)

# pull in fao classification data

sofia_status <-
  read_csv(here("data", "StockSTatus_FAOSofia.csv")) %>%
  janitor::clean_names() %>%
  select(area, sp_group, name, species, x2017) %>%
  rename(
    status = x2017,
    bad_fao_code = area,
    isscaap_number = sp_group,
    common_name = name,
    scientific_name = species
  ) %>%
  mutate(fao_area_code = as.numeric(bad_fao_code)) %>%
  filter(!is.na(status), status != "?",!is.na(fao_area_code)) %>%
  mutate(
    stock_complex = paste(common_name, isscaap_number, fao_area_code, sep = '_'),
    full_status = status
  )

stats <-
  str_trim(str_split(sofia_status$status, pattern = "[[:punct:]]", simplify = TRUE)[, 1])


sofia_status$status <- stats

# there are MANY mispellings of scientific names in SOFIA

correct_sofia <- sofia_status %>% 
  select(common_name, scientific_name) %>% 
  rename(bad_scientific_name = scientific_name) %>% 
  unique() %>% 
  mutate(sciname = map(common_name, ~ taxize::comm2sci(com = .x, db = "worms")[[1]]))

correct_sofia <-  correct_sofia %>% 
  mutate(sci_match = map_lgl(sciname, ~length(.x) > 0)) %>% 
  filter(sci_match == TRUE) %>% 
  mutate(good_scientific_name = map_chr(sciname, ~.x[1])) %>% 
  select(common_name,good_scientific_name )

# correct_sofia$foo <- correct_sofia$bad_scientific_name != correct_sofia$good_scientific_name

sofia_status <- sofia_status %>%
  left_join(correct_sofia, by = "common_name") %>%
  mutate(scientific_name = ifelse(
    is.na(good_scientific_name),
    scientific_name,
    good_scientific_name
  )) %>% 
  select(-good_scientific_name)

# replicate rows for multiple status

# tmp <- sofia_status %>%
#   group_by(stock_complex) %>%
#   nest() %>%
#   ungroup()
#
# repfoo <- function(x){
#
#   statusi <- str_trim(str_split(x$status, pattern = "[[:punct:]]", simplify = TRUE))
#
#   out <- map_df(statusi, ~ x) %>%
#     mutate(status = statusi)
#
# }
# ogss <- sofia_status
#
# sofia_status <- tmp %>%
#   mutate(data = map(data, repfoo)) %>%
#   unnest(cols = data)

# merge with fao data

# fuck <- sofia_status %>% 
#   filter(fao_area_code == 67)
# 
# 
# you <- fao_status <- fao %>%
#   group_by(
#     year,
#     fao_area,
#     fao_area_code,
#     isscaap_group,
#     isscaap_number,
#     common_name,
#     scientific_name
#   ) %>%
#   summarise(catch = sum(capture), .groups = "drop") %>% 
#   filter(fao_area_code == 67) 
# 
# 
# wtf <- you %>% 
#   filter(scientific_name %in% unique(fuck$scientific_name))
# 

fao_status <- fao %>%
  group_by(
    year,
    fao_area,
    fao_area_code,
    isscaap_group,
    isscaap_number,
    # common_name,
    scientific_name
  ) %>%
  summarise(catch = sum(capture), .groups = "drop") %>%
  left_join(
    sofia_status,
    by = c(
      "fao_area_code",
      "isscaap_number",
      # "common_name",
      "scientific_name"
    )
  ) %>%
  filter(!is.na(status)) %>%
  mutate(stockid = paste(common_name, isscaap_number, fao_area_code, sep = '_')) %>%
  group_by(stockid) %>%
  mutate(first_year = year[catch > 0 & !is.na(catch)][1],
         last_year = last(year[catch >= 0 & !is.na(catch)])) %>%
  filter(year >= first_year, year <= last_year) %>%
  mutate(interp_catch = zoo::na.approx(catch, na.rm = FALSE)) %>%
  mutate(catch = dplyr::case_when(is.na(catch) ~ interp_catch,
                                  TRUE ~ catch)) %>%
  filter(!is.na(catch)) %>%
  select(-first_year,-last_year,-interp_catch) %>%
  filter(!is.na(fao_area_code),!is.na(common_name),
         status != 'N') %>%
  mutate(status = forcats::fct_relevel(status, c("O", "F", "U")))


n_status_plot <- fao_status %>%
  group_by(year, status) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  ggplot() +
  geom_area(aes(year, n, fill = status), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao_status$year), max(fao_status$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
# made it here

n_region_status_plot <- fao_status %>%
  group_by(year, fao_area, status) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  group_by(year, fao_area) %>%
  mutate(pn = n / sum(n)) %>%
  ggplot() +
  geom_area(aes(year, pn, fill = status), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao_status$year), max(fao_status$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap( ~ fao_area)

catch_status_plot <- fao_status %>%
  group_by(year, status) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  ggplot() +
  geom_area(aes(year, total_catch, fill = status),
            color = "black",
            alpha = 0.75) +
  labs(x = "Year") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao_status$year), max(fao_status$year)), expand = c(0, 0)) +
  scale_y_comma(name = "Total Capture (mt)")


# load sar data

load(here::here("data", "PNAS-data", "TBPdata.Rdata"))

extra_sar <- readr::read_csv(here("data", "sar.csv")) %>%
  mutate(country = countrycode::countrycode(Country, "country.name", 'un.name.en')) %>%
  left_join(fao_country_region_key, by = "country") %>%
  filter(ptc == 1) %>%
  mutate(AreaKm2_0_200 = tc) %>%
  select(fao_area_code, SAR, AreaKm2_0_200)

regional_sar <- TBPdata@SummaryByRegion

sar_to_fao_region <-
  read_csv(file = here("data", "sar_region_lookup-maually-created.csv"))

approx_sar_by_fao_region <- regional_sar %>%
  right_join(sar_to_fao_region, by = "RegionName") %>%
  select(fao_area_code, SAR, AreaKm2_0_200) %>%
  bind_rows(extra_sar) %>%
  group_by(fao_area_code) %>%
  summarise(sar = weighted.mean(SAR, w = AreaKm2_0_200, na.rm = TRUE))

# load in ram status

load(here("data", "RAMCore.rdata"))

ram_to_region <-
  data_frame(stockid = colnames(UvT), region = Region)

u_v_t <- UvT %>%
  as.data.frame() %>%
  mutate(year = rownames(.)) %>%
  gather(stockid, u_umsy,-year) %>%
  left_join(ram_to_region, by = "stockid") %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup()

clean_u <- function(u) {
  if (any(is.na(u$u_umsy))) {
    first_u <- which(!is.na(u$u_umsy))[1]
    
    u$u_umsy[1:(first_u - 1)] <- 0
    
    last_u <- last(which(!is.na(u$u_umsy)))
    if (last_u < nrow(u)) {
      u$u_umsy[(last_u + 1):(nrow(u))] <- u$u_umsy[last_u]
    }
  }
  
  return(u)
}

tempu <- u_v_t %>%
  mutate(u2 = map(data, safely(clean_u)))

u_worked <- map(tempu$u2, "error") %>% map_lgl(is.null)

u_v_t <- tempu %>%
  filter(u_worked) %>%
  mutate(u2 = map(u2, "result")) %>%
  select(-data) %>%
  unnest()


b_v_t <- BvT %>%
  as.data.frame() %>%
  mutate(year = rownames(.)) %>%
  gather(stockid, b_bmsy,-year)

status_t <- u_v_t %>%
  left_join(b_v_t, by = c("year", "stockid"))

breaks <- c(0.8, 1.2)

breaks <- c(0, breaks, Inf)

labels <- c("over", "fully", "under")

status_t <- status_t %>%
  mutate(bin = cut(b_bmsy, breaks = breaks, labels = labels)) %>%
  mutate(year = as.numeric(year))

# add state space u estimates to FAO

fao_area_to_ram <-
  readxl::read_xlsx(here("data", "FAO SOFIA Stock status table.xlsx"), sheet = "cleaned_table") %>%
  janitor::clean_names() %>%
  mutate(area = as.numeric(fao_area)) %>%
  mutate(tempid = 1:nrow(.)) %>%
  select(-sar, -fmi)


state_space_results <-
  readxl::read_xlsx(here("data", "state_space_results_08_28_2018 by region.xlsx"),
                    sheet = "state_space_results_08_28_2018") %>%
  janitor::clean_names()


fao_status <- fao_status %>%
  left_join(fao_area_to_ram %>% select(-fao_area),
            by = c("fao_area_code" = "area")) %>%
  left_join(
    state_space_results %>% filter(variable == "UvU"),
    by = c("year", "state_space_region" = "region")
  ) %>%
  left_join(approx_sar_by_fao_region,
            by = c("fao_area_code" = "fao_area_code")) %>%
  rename(mean_u_umsy = dlm_geomean)

fao_ram_u_plot <- fao_status %>%
  ggplot(aes(year,
             (mean_u_umsy))) +
  geom_point() +
  facet_wrap( ~ state_space_region, scales = "free_y")


ram_status_plot <- status_t %>%
  group_by(year, bin) %>%
  summarise(n = length(b_bmsy)) %>%
  ggplot() +
  geom_area(aes(year, n, fill = bin), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao_status$year), max(fao_status$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# link up generic effort to fao 2011



rough_fao_region_effort <- effort_data %>%
  filter(effort_type == "nominal") %>%
  right_join(effort_region_to_country, by = "region") %>%
  right_join(fao_country_region_key, by = "country") %>%
  group_by(fao_area_code, year) %>%
  summarise(effort_index = sum(effort)) %>%
  ungroup()



temp_fao <- fao_status %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup()

assign_effort <-
  function(fao_stock,
           data,
           effort,
           fao_catch,
           scalar = 1000) {
    # data <- temp_fao$data[[which(temp_fao$stockid == huh)]]
    #
    # fao_catch <- fao
    #
    # effort <- rous_data
    #
    # fao_stock <- temp_fao$stockid[which(temp_fao$stockid == huh)]
    
    comm_name <-
      str_split(fao_stock, pattern = '_')[[1]][1] %>% str_remove_all("(\\d)|(-)")
    
    isscp_number <-
      as.numeric(str_split(fao_stock, pattern = '_')[[1]][2])
    
    fao_code <-
      as.numeric(str_split(fao_stock, pattern = '_')[[1]][3])
    
    fao_matches <- fao_catch %>% {
      if (any(.$common_name == comm_name, na.rm = TRUE)) {
        filter(., common_name == comm_name & fao_area_code == fao_code)
        # filter(fao, common_name == comm_name & fao_area_code == fao_code)
        #
        
      } else {
        filter(.,
               isscaap_number == isscp_number &
                 fao_area_code == fao_code)
        
      }
    }
    matched_effort <- effort %>%
      filter(area == fao_code) %>%
      filter(country %in% unique(fao_matches$country)) %>%
      group_by(year, area) %>%
      summarise(effort_index = sum(effort_cell_reported_nom) / scalar,
                .groups = "drop") %>%
      ungroup()
    
    data <- data %>%
      left_join(matched_effort, by = c("year", "fao_area_code" =  "area"))
    
    
  }


temp_fao <- temp_fao %>%
  mutate(data = map2(
    stockid,
    data,
    assign_effort,
    fao_catch = fao,
    effort = rous_data
  ))


fao_status <- temp_fao %>%
  unnest(cols = data)


# fao_status %>%
#   group_by(year, area) %>%
#   summarise(effort_index = unique(effort_index)) %>%
#   ggplot(aes(year, effort_index)) +
#   geom_line() +
#   facet_wrap(~area)

# temp <- fao_status %>%
#   # filter(area %in% 67) %>%
#   group_by(stockid) %>%
#   nest() %>%
#   ungroup()
#
# huh <- temp$stockid[4]
#
# wtf <- fao_status %>%
#   filter(stockid == huh)


fao_status <- fao_status %>%
  group_by(stockid) %>%
  mutate(
    lifetime_catch = sum(catch, na.rm = TRUE),
    min_catch = min(catch, na.rm = TRUE),
    nas = sum(is.na(catch))
  ) %>%
  ungroup() %>%
  filter(lifetime_catch > 50000,
         min_catch > 0,
         nas == 0,!is.na(common_name)) %>%
  ungroup() %>%
  mutate(sp_group = isscaap_number) %>%
  # left_join(fao_species %>% select(isscaap_group, isscaap_number) %>% unique(), by = c("sp_group" = "isscaap_number")) %>%
  mutate(species = scientific_name)


support_data <-
  list(mean_regional_isscaap_fmi = mean_regional_isscaap_fmi)

# areas <- c(67, 57, 37,71)

areas <- unique(fao_status$fao_area_code)

# fao_status %>%
#   filter(year == max(year)) %>%
#   summarise(catch = sum(catch))
#
# fao %>%
#   filter(year == max(year)) %>%
#   summarise(catch = sum(capture))

# annnnnd try and run assessments
# browser()
if (run_sofia_comparison == TRUE) {
  # future::plan(future::multiprocess, workers = 4)
  
  set.seed(42)
  
  future::plan(multisession, workers = n_cores)
  
  fao_status_fits <- fao_status %>%
    # filter(fao_area_code %in% 67) %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup() %>%
    # sample_n(3) %>%
    # slice(10) %>% 
    mutate(
      fits = future_map(
        data,
        safely(fit_fao),
        support_data = support_data,
        default_initial_state = NA,
        default_initial_state_cv = NA,
        min_effort_year = 1975,
        engine = "stan",
        cores = 2,
        .progress = TRUE,
        .options = future_options(
          globals = support_data,
          packages = c("tidyverse", "sraplus")
        )
      )
    )
  

  write_rds(fao_status_fits,
            path = file.path(results_path, "fao_status_fits.rds"))
  
} else {
  fao_status_fits <-
    read_rds(path = file.path(results_path, "fao_status_fits.rds"))
  
  
}

reserve <- fao_status_fits

ugh <- map_dbl(fao_status$data, ~ unique(.x$area)) == 67

fao_worked <- map(fao_status_fits$fits, "error")


fao_worked <-
  map(fao_status_fits$fits, "error") %>% map_lgl(is.null)

fao_status_fits <- fao_status_fits %>%
  filter(fao_worked)

fao_status_fits <- fao_status_fits %>%
  mutate(fits = map(fits, "result"))

faofits <- fao_status_fits %>%
  select(-data) %>%
  unnest(cols = fits) %>%
  filter(variable == "b_div_bmsy") %>%
  mutate(bin = cut(mean, breaks = breaks, labels = labels)) %>%
  separate(stockid, c("name", "speciesgroup", "area"), sep = '_') %>%
  mutate(fao_area_code = as.numeric(area))


# wtf <- faofits %>%
#   filter(stockid == stockid[1],
#          str_detect(data, "cpue")) %>%
#   ggplot(aes(year, mean, color = data)) +
#   geom_line() +
#   facet_wrap(~data)

fao_status_data <- fao_status %>%
  mutate(
    bin = dplyr::case_when(
      status == "F" ~ "fully",
      status == "U" ~ "under",
      status == "O" ~ "over",
      TRUE ~ "unknown"
    )
  ) %>%
  mutate(data = "FAO Report 569",
         mean = rep(1, nrow(.))) %>%
  select(fao_area_code, year, bin, data, mean)


fao_regions_key <- read_csv(here("data", "fao_regions_key.csv"))

temp <- faofits %>%
  select(fao_area_code, year, bin, data, mean) %>%
  bind_rows(fao_status_data) %>%
  filter(!is.na(fao_area_code),
         !fao_area_code %in% c(88, 48),!bin == "unknown") %>%
  left_join(fao_regions_key, by = c("fao_area_code" = "fao_region_num")) %>%
  group_by(fao_region_name, fao_area_code, year, bin, data) %>%
  summarise(n = length(mean)) %>%
  group_by(fao_region_name, fao_area_code, data, year) %>%
  mutate(pn = round(n / sum(n), 2)) %>%
  group_by(fao_region_name, fao_area_code, data) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  tidyr::complete(bin,
                  nesting(fao_region_name, fao_area_code, year, data),
                  fill = list(n = 0, pn = 0)) %>%
  filter(bin == "over") %>%
  mutate(data = case_when(data == "u_umsy" ~ "ram_u_umsy", TRUE ~ data)) %>%
  mutate(data = fct_reorder(data, pn))

# temp$pn[temp$fao_region_name == "Pacific, Northeast" & temp$bin == "over" & temp$data == "FAO Report 569"] <- 0.09


fao_status_plot <- temp %>%
  ggplot(aes(data, pn, fill = fao_region_name)) +
  geom_col(
    size = 1,
    alpha = 1,
    size = 1,
    show.legend = FALSE,
    color = "black"
  ) +
  coord_flip() +
  facet_grid(fao_region_name ~ .) +
  theme_minimal() +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = paste0(seq(0, 100, by = 25), "%"),
    name = "Percent of Stocks Overfished"
  ) +
  theme(axis.text.x = element_text(size = 8),
        axis.title.y = element_blank()) +
  scale_fill_discrete(name = "FAO Category")

fao_status_plot


if (!dir.exists(here("data", "FAO_AREAS_NOCOASTLINE"))) {
  download.file(url = "http://www.fao.org/figis/geoserver/area/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=area:FAO_AREAS_NOCOASTLINE&outputFormat=SHAPE-ZIP",
                destfile = here("data", "FAO_AREAS_NOCOASTLINE.zip"))
  
  unzip(
    here("data", "FAO_AREAS_NOCOASTLINE.zip"),
    exdir = here("data", "FAO_AREAS_NOCOASTLINE")
  )
  
}

fao_areas <- sf::st_read(here('data', "FAO_AREAS_NOCOASTLINE")) %>%
  janitor::clean_names()

fao_areas <- fao_areas %>%
  group_by(f_area) %>%
  nest() %>%
  mutate(geometry = map(data, st_union)) %>%
  select(-data)

fao_areas = fao_areas %>%
  unnest(cols = geometry) %>%
  ungroup() %>%
  sf::st_as_sf() %>%
  sf::st_simplify() %>%
  mutate(fao_area_code = as.numeric(f_area))

fao_area_status <- fao_areas %>%
  left_join(temp %>% mutate(f_area = as.character(fao_area_code)), by = "f_area") %>%
  filter(!is.na(data))


world_map <-
  rnaturalearth::ne_countries(returnclass = "sf", scale = "small") #%>%
# sf::st_union() %>%
# sf::st_polygonize()



po_map_plot <- fao_area_status %>%
  filter(!is.na(data)) %>%
  ggplot() +
  geom_sf(aes(fill = pn)) +
  geom_sf(
    data = world_map,
    fill = "lightgrey",
    color = "black",
    size = 0.05
  ) +
  facet_wrap( ~ data) +
  scale_fill_viridis_b(
    option = "magma",
    labels = percent,
    name = "% Overfished",
    guide = guide_coloursteps(
      barwidth = ggplot2::unit(8, "lines"),
      reverse = FALSE,
      axis.colour = "white",
      axis.linewidth = 1,
      show.limits = TRUE
    )
  ) +
  ggthemes::theme_base() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "white"))

po_map_plot

ggsave("percent_overfished_by_region.pdf", po_map_plot)

ram_catches <- ram_data %>%
  select(scientificname, catch, year, primary_country, primary_FAOarea) %>%
  mutate(fao_area_code = as.numeric(primary_FAOarea)) %>%
  group_by(year, scientificname, fao_area_code) %>%
  summarise(total_catch = sum(catch, na.rm = TRUE)) %>%
  rename(scientific_name = scientificname) %>%
  ungroup() %>%
  mutate(combo = paste(scientific_name, fao_area_code, sep = '_'),
         data = "RAM")

fao_catches <- fao_status %>%
  select(year, scientific_name, fao_area_code, catch) %>%
  group_by(year, scientific_name, fao_area_code) %>%
  summarise(total_catch = sum(catch, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(combo = paste(scientific_name, fao_area_code, sep = '_'),
         data = "FAO") %>%
  filter(combo %in% unique(ram_catches$combo))

fao_ram_catches <- ram_catches %>%
  bind_rows(fao_catches)  %>%
  group_by(scientific_name, fao_area_code, data) %>%
  mutate(scatch = as.numeric(scale(total_catch))) %>%
  ungroup()



fao_ram_catches %>%
  ggplot(aes(
    year,
    scatch,
    color = data,
    group = interaction(data, scientific_name)
  )) +
  geom_line(alpha = 0.5) +
  facet_wrap( ~ fao_area_code)


fao_ram_catches %>%
  select(-total_catch) %>%
  pivot_wider(names_from = data, values_from = scatch) %>%
  ggplot(aes(RAM, FAO)) +
  geom_point(alpha = 0.5) +
  geom_abline(slopee = 1, intercept = 0) +
  facet_wrap( ~ fao_area_code)


a = fao_ram_catches %>%
  select(-total_catch) %>%
  pivot_wider(names_from = data, values_from = scatch) %>%
  filter(!is.na(FAO)) %>%
  select(RAM, FAO)







ram_status <- ram_data %>%
  select(scientificname,
         b_v_bmsy,
         year,
         primary_country,
         primary_FAOarea) %>%
  mutate(fao_area_code = as.numeric(primary_FAOarea)) %>%
  group_by(year, scientificname, fao_area_code) %>%
  summarise(mean_b = mean(b_v_bmsy, na.rm = TRUE)) %>%
  left_join(fao_regions_key, by = c("fao_area_code" = "fao_region_num")) %>%
  ungroup() %>%
  filter(!is.na(mean_b)) %>%
  rename(scientific_name = scientificname) %>%
  mutate(combo = paste(scientific_name, fao_area_code, sep = '_'),
         data = "RAM")


rng_status <- ram_status %>%
  mutate(mean_b = runif(n(), 0, 2.5),
         data = "rng")

eo_status <- ram_status %>%
  mutate(mean_b = sample(c(1.5, .5), n(), replace = TRUE),
         data = "either or")



ram_status %>%
  group_by(year, fao_area_code) %>%
  summarise(
    mmb = mean(mean_b),
    po = mean(mean_b < 0.8),
    n = n_distinct(scientific_name)
  ) %>%
  ggplot(aes(year, po)) +
  geom_line() +
  facet_wrap( ~ fao_area_code)



ram_status %>%
  group_by(fao_area_code, scientific_name) %>%
  filter(year == max(year)) %>%
  group_by(fao_area_code) %>%
  summarise(
    mmb = mean(mean_b),
    po = mean(mean_b < 0.8),
    n = n_distinct(scientific_name)
  ) %>%
  ggplot(aes(factor(fao_area_code), po, fill = mmb)) +
  geom_col()



ram_fao_status <- faofits %>%
  left_join(fao_names, by = c("name" = "common_name")) %>%
  mutate(combo = paste(scientific_name, fao_area_code, sep = '_')) %>%
  filter(combo %in% unique(ram_status$combo)) %>%
  rename(mean_b = mean) %>%
  select(scientific_name, fao_area_code, year, data,  mean_b) %>%
  bind_rows(ram_status %>% select(scientific_name, fao_area_code, year, data, mean_b)) %>%
  bind_rows(rng_status %>% select(scientific_name, fao_area_code, year, data, mean_b)) %>%
  bind_rows(eo_status %>% select(scientific_name, fao_area_code, year, data, mean_b)) %>%
  mutate(bin = cut(mean_b, breaks = breaks, labels = labels)) %>%
  group_by(fao_area_code, scientific_name) %>%
  mutate(max_ram_year = max(year[data == "RAM"])) %>%
  ungroup() %>%
  filter(year == max_ram_year) %>% # filter to the most recent year that RAM is available for
  ungroup()



ram_fao_status_bins <- ram_fao_status %>%
  group_by(fao_area_code, data, bin) %>%
  count() %>%
  group_by(data, fao_area_code) %>%
  mutate(pn = n / sum(n)) %>%
  ungroup() %>%
  tidyr::complete(bin, nesting(fao_area_code, data),
                  fill = list(n = 0, pn = 0)) %>%
  mutate(data = case_when(data == "u_umsy" ~ "ram_u_umsy", TRUE ~ data)) %>%
  mutate(data = fct_reorder(data, pn))


ram_v_sraplus_map <- fao_areas %>%
  left_join(ram_fao_status_bins, by = "fao_area_code")

ram_v_sraplus_map_plot <- ram_v_sraplus_map %>%
  filter(bin == "over") %>%
  ggplot() +
  geom_sf(aes(fill = pn)) +
  geom_sf(
    data = world_map,
    fill = "lightgrey",
    color = "black",
    size = 0.05
  ) +
  facet_wrap( ~ data) +
  scale_fill_viridis_b(
    option = "magma",
    labels = percent,
    name = "% Overfished",
    guide = guide_coloursteps(
      barwidth = ggplot2::unit(8, "lines"),
      reverse = FALSE,
      axis.colour = "white",
      axis.linewidth = 1,
      show.limits = TRUE
    )
  ) +
  ggthemes::theme_base() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "white"))

ram_v_sraplus_map_plot

# now do some accuracy and RMSE calculations

faofits %>%
  ggplot(aes(year, pmin(2.5, mean), group = name)) +
  geom_line() +
  facet_grid(area ~ data, scales = "free_y")


ram_status %>%
  ggplot(aes(year, pmin(2.5, mean_b), group = scientific_name)) +
  geom_line(alpha = 0.5) +
  facet_wrap( ~ fao_area_code, scales = "free_y")

ram_v_sraplus_acc <- ram_fao_status %>%
  ungroup() %>%
  select(-year, -max_ram_year, -mean_b) %>%
  pivot_wider(names_from = data, values_from = bin) %>%
  pivot_longer(
    c(-scientific_name, -fao_area_code, -RAM),
    names_to = "data",
    values_to = "bin"
  ) %>%
  mutate(correct = bin == RAM)

ram_v_sraplus_acc_plot <- ram_v_sraplus_acc %>%
  group_by(data) %>%
  summarise(accuracy = mean(correct, na.rm = TRUE)) %>%
  ggplot(aes(reorder(data,-accuracy), accuracy)) +
  geom_col() +
  scale_x_discrete(name = "Data Source") +
  scale_y_continuous(labels = label_percent(accuracy = 1), name = "Accuracy")


ram_v_sraplus_rmse <- ram_fao_status %>%
  ungroup() %>%
  select(-year, -max_ram_year, -bin) %>%
  pivot_wider(names_from = data, values_from = mean_b) %>%
  pivot_longer(
    c(-scientific_name, -fao_area_code, -RAM),
    names_to = "data",
    values_to = "mean_b"
  ) %>%
  mutate(se = (RAM - mean_b) ^ 2)



ram_v_sraplus_rmse %>%
  ggplot(aes(RAM, pmin(2.5, mean_b), color = data)) +
  geom_point() +
  facet_wrap( ~ data, scales = "free_y")

ram_v_sraplus_rmse_plot <- ram_v_sraplus_rmse %>%
  group_by(data) %>%
  summarise(rmse = sqrt(mean(se, na.rm = TRUE)),
            r2 = yardstick::rsq_vec(RAM, mean_b)) %>%
  ggplot(aes(reorder(data, rmse), rmse)) +
  geom_col() +
  scale_x_discrete(name = "Data Source") +
  scale_y_continuous(name = "Root Mean Squared Error")

ram_v_sraplus_r2_plot <- ram_v_sraplus_rmse %>%
  group_by(data) %>%
  summarise(rmse = sqrt(mean(se, na.rm = TRUE)),
            r2 = yardstick::rsq_vec(RAM, mean_b)) %>%
  ggplot(aes(reorder(data,-r2), r2)) +
  geom_col() +
  scale_x_discrete(name = "Data Source") +
  scale_y_continuous(name = "R Squared")


ram_v_sraplus_fao_acc <- ram_v_sraplus_acc %>%
  group_by(fao_area_code, data) %>%
  summarise(acc = mean(correct, na.rm = TRUE)) %>%
  ungroup()

ram_v_sraplus_acc_map <- fao_areas %>%
  left_join(ram_v_sraplus_fao_acc, by = "fao_area_code")

ram_v_sraplus_acc_map_plot <- ram_v_sraplus_acc_map %>%
  ggplot() +
  geom_sf(aes(fill = acc)) +
  geom_sf(
    data = world_map,
    fill = "lightgrey",
    color = "black",
    size = 0.05
  ) +
  facet_wrap( ~ data) +
  ggthemes::theme_base() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "white"))

ram_v_sraplus_acc_map_plot


ram_v_sraplus_fao_rmse <- ram_v_sraplus_rmse %>%
  group_by(fao_area_code, data) %>%
  summarise(rmse = sqrt(mean(se, na.rm = TRUE))) %>%
  ungroup()

ram_v_sraplus_rmse_map <- fao_areas %>%
  left_join(ram_v_sraplus_fao_rmse, by = "fao_area_code")

ram_v_sraplus_rmse_map_plot <- ram_v_sraplus_rmse_map %>%
  ggplot() +
  geom_sf(aes(fill = pmin(2, rmse))) +
  geom_sf(
    data = world_map,
    fill = "lightgrey",
    color = "black",
    size = 0.05
  ) +
  facet_wrap( ~ data) +
  ggthemes::theme_base() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "white"))

ram_v_sraplus_rmse_map_plot



# run RAM tests ------------------------------------------------------

# for now, only including things that have total biomass estimates

if (run_ram_tests) {
  ram_fit_data <- ram_data %>%
    select(
      year,
      stockid,
      scientificname,
      commonname,
      year,
      catch,
      b_v_bmsy,
      u_v_umsy,
      total_biomass
    ) %>%
    group_by(stockid, scientificname) %>%
    nest() %>%
    ungroup()
  
  fit_model <-
    function(sciname,
             data,
             survey_q = 1,
             estimate_proc_error = FALSE,
             #set to TRUE to estimate process error
             estimate_shape = FALSE) {

      index_years <- data$year[!is.na(data$b_v_bmsy)]
      
      driors <- format_driors(taxa = sciname,
                              catch = data$catch,
                              years = data$year,
                              index =  (data$b_v_bmsy * survey_q * 100)[data$year %in% index_years],
                              index_years = index_years,
                              initial_state = NA,
                              initial_state_cv = 0.05)
      
      # browser()
      # sraplus::plot_driors(driors)
      # 
      fit <- fit_sraplus(driors = driors,
                         engine = "stan",
                         model = "sraplus_tmb",
                         estimate_shape = estimate_shape, 
                         estimate_proc_error = estimate_proc_error,
                         estimate_initial_state = TRUE,
                         n_keep = 2000)
      
      
      # b <- fit$results %>%
      #   filter(variable == "b_t")
      
      # plot(b$mean, data$total_biomass)
      # abline(a = 0, b = 1)
      #
      out <- list(fit = fit,
                  driors = driors)
    }
  
  future::plan(multiprocess, workers = 8)
  
  a <- Sys.time()
  set.seed(42)
  ram_fit_data <- ram_fit_data %>%
    # sample_n(10) %>% 
    # filter(stockid == "ARFLOUNDBSAI") %>% 
    mutate(
      fit = future_map2(
        scientificname,
        data,
       safely(fit_model),
        estimate_proc_error = TRUE,
        estimate_shape = FALSE,
        .progress = TRUE
      )
    )
  
  fit_worked <- map_lgl(map(ram_fit_data$fit, "error"), is.null)
  Sys.time() - a
  
  ram_fit_data <- ram_fit_data %>%
    filter(fit_worked)
  
  ram_fit_data$fit <- map(ram_fit_data$fit, "result")
  
  write_rds(ram_fit_data, path = file.path(results_path,"ram_tests.rds"))
  
  
  i <- 8
  
  sraplus::plot_prior_posterior(ram_fit_data$fit[[i]]$fit, ram_fit_data$fit[[i]]$driors)
  

  compare_to_ram <- function(data, fit){
    
    # data <- ram_fit_data$data[[1]]
    # 
    # fit <- ram_fit_data$fit[[1]]
    
    comparison <- tibble(observed = data$b_v_bmsy) %>% 
      bind_cols(fit$fit$results[fit$fit$results$variable == "b_div_bmsy",])
    
  }
  
  ram_comparison <- map2_df(ram_fit_data$data, ram_fit_data$fit, compare_to_ram,.id = "stock")
  
  ram_comparison %>% 
    ggplot(aes(observed, mean, color = stock)) + 
    geom_point(show.legend = FALSE, alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0)
  
  ram_comparison %>% 
    ggplot(aes(observed, mean)) + 
    geom_bin2d(show.legend = TRUE, alpha = 0.5) + 
    geom_abline(slope = 1, intercept = 0)

  
  
}

# run RAM comparison  -----------------------------------------------------------



# re run, but with a pared down focus on area 67

# filter down to just species that have matches between FAO and ram by fao area


if (run_ram_comparison == TRUE) {
  
  ram_comp_data <- ram_data %>%
    select(
      year,
      stockid,
      country,
      primary_FAOarea,
      scientificname,
      commonname,
      year,
      catch,
      b_v_bmsy,
      u_v_umsy,
      total_biomass
    ) %>%
    rename(scientific_name = scientificname) %>% 
    mutate(fao_area_code = as.numeric(primary_FAOarea)) %>% 
    left_join(fao_species, by = "scientific_name")
  
  # join regional U/Umsy and SAR data
  ram_comp_data <- ram_comp_data %>%
    left_join(fao_area_to_ram %>% select(-fao_area),
              by = c("fao_area_code" = "area")) %>%
    left_join(
      state_space_results %>% filter(variable == "UvU"),
      by = c("year", "state_space_region" = "region")
    ) %>%
    left_join(approx_sar_by_fao_region,
              by = c("fao_area_code" = "fao_area_code")) %>%
    rename(mean_u_umsy = dlm_geomean)

  # assign effort data to RAM
  
  temp_ram <- ram_comp_data %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup()
  
  assign_effort <-
    function(ram_stock,
             data,
             effort,
             fao_catch,
             scalar = 1000) {
      # data <- temp_fao$data[[which(temp_fao$stockid == huh)]]
      #
      # fao_catch <- fao
      #
      # effort <- rous_data
      #
      # fao_stock <- temp_fao$stockid[which(temp_fao$stockid == huh)]
      
      comm_name <- unique(data$commonname)
        
      isscp_number <- unique(data$isscaap_number)

      if (is.na(isscp_number)){
       
        isscp_number <- 33 
      }
      
      fao_code <- unique(data$fao_area_code)

      fao_matches <- fao_catch %>% {
        if (any(.$common_name == comm_name & .$fao_area_code == fao_code, na.rm = TRUE)) {
          filter(., common_name == comm_name & fao_area_code == fao_code)
          # filter(fao, common_name == comm_name & fao_area_code == fao_code)
          #
          
        } else {
          filter(.,
                 isscaap_number == isscp_number &
                   fao_area_code == fao_code)
          
        }
      }
      matched_effort <- effort %>%
        filter(area == fao_code) %>%
        filter(country %in% unique(fao_matches$country)) %>%
        group_by(year, area) %>%
        summarise(effort_index = sum(effort_cell_reported_nom) / scalar,
                  .groups = "drop") %>%
        ungroup()
      
      data <- data %>%
        left_join(matched_effort, by = c("year", "fao_area_code" =  "area"))
      
      
    }
  
  
  temp_ram <- temp_ram %>%
    mutate(data = map2(
      stockid,
      data,
      assign_effort,
      fao_catch = fao,
      effort = rous_data
    ))
  
  
  ram_comp_data <- temp_ram %>%
    unnest(cols = data)
  
  rous_data %>% 
    group_by(year, area) %>% 
    summarise(effort = mean(effort_cell_reported_nom, na.rm = TRUE)) %>% 
    ungroup() %>% 
    ggplot(aes(year, effort)) + 
    geom_line() + 
    facet_wrap(~area, scales = "free_y")
  
  
  ram_comp_data %>% 
    group_by(year, fao_area_code) %>% 
    summarise(effort = mean(effort_index, na.rm = TRUE)) %>% 
    ungroup() %>% 
    ggplot(aes(year, effort)) + 
    geom_line() + 
    facet_wrap(~fao_area_code, scales = "free_y")
  
  set.seed(24)
  
  future::plan(multisession, workers = n_cores)
  
  ram_status_fits <- ram_comp_data %>%
    filter(fao_area_code %in% 67) %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup() %>% 
    # sample_n(1) %>%
    mutate(
      fits = future_map(
        data,
       safely( fit_ram),
        support_data = support_data,
        default_initial_state = NA,
        default_initial_state_cv = NA,
        min_effort_year = 1975,
        engine = "stan",
        cores = 2,
        .progress = TRUE,
        .options = future_options(
          globals = support_data,
          packages = c("tidyverse", "sraplus")
        )
      )
    )
  
  
  check <- ram_status_fits$fits[[1]]$result

  check %>%
    filter(variable == "b_div_bmsy")  %>%
    ggplot() +
    geom_line(aes(year, mean, color = data)) +
    geom_point(data = ram_status_fits$data[[1]], aes(year, b_v_bmsy))
  #
  write_rds(ram_status_fits,
            path = file.path(results_path, "ram_status_fits.rds"))
  
} else {
  ram_status_fits <-
    read_rds(path = file.path(results_path, "ram_status_fits.rds"))
  
  
}




# ----save-plots----------------------------------------------------------

plots <- ls()[str_detect(ls(), "_plot")]

save(file = here("results", results_name, "paper_plots.Rdata"),
     list = plots)
