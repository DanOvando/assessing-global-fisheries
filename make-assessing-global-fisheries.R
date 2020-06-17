
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

n_cores <- 2
# number of cores for parallel processing

future::plan(multiprocess, workers = n_cores)


# options(mc.cores = 1)

results_name <- "v0.5"


results_description <- "publication version of results updated VOI from version 0.5 "

run_voi_models <- TRUE
# sub options for run_voi_models
fit_models <- FALSE

write_results <- TRUE

process_fits <- FALSE

# other things
run_case_studies <- TRUE

run_continent_examples <- TRUE

run_ei_example <- TRUE

run_sofia_comparison <- FALSE

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
                          is.na(total_biomass) | is.na(catch) | catch == 0)) %>%
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


if (run_voi_models == TRUE){
  
  
  stress_tests <- data_frame(
    use_index = sample(c(TRUE, FALSE),  draws, replace = TRUE),
    index_type = sample(c("survey"), draws, replace = TRUE),
    initial_state_type = sample(c("unfished", "known", "prior","heuristic"), draws, replace = TRUE),
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
    mutate(error_cv = map_dbl(error_cv, ~ ifelse(.x == 0, 1e-6, runif(1, 0.01, .25)))) %>%
    mutate(
      b_ref_type = ifelse(
        initial_state_type == "prior",
        "b",
        b_ref_type
      )
    )
  
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
  #   b_ref_type = sample(c("k", "b"),  draws, replace = TRUE),
  #   f_ref_type = sample(c("f", "fmsy"),  draws, replace = TRUE),
  #   estimate_shape = sample(c(FALSE, TRUE), draws, replace = TRUE),
  #   estimate_proc_error = sample(c(FALSE, TRUE), draws, replace = TRUE)
  # ) %>%
  #   mutate(error_cv = map_dbl(error_cv, ~ ifelse(.x == 0, 1e-6, runif(1, 0.01, .25)))) %>%
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
          assessments = c("cmsy","sraplus"),
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
    ram_fits <-  readRDS(file = here::here("results", results_name, "ram_experiments.rds"))
    
  }
  
  
  # wtf <- map(ram_fits$fits,'error') %>% map_lgl(is.null)
  #
  # a <- ram_fits %>%
  #   filter(!wtf)
  
  # process fits ------------------------------------------------------------
  
  if (process_fits == TRUE) {
    fits <- list.files(here::here("results", results_name, "experiments"))
    
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
      mutate(fits = map2(
        experiment,
        data,
        safely(load_foo),
        results_name = results_name,
        produce = "summary",
        plots = FALSE))
    
    
    saveRDS(ram_fits,
            file = here::here("results", results_name, "ram_results.Rds"))
  } else {
    ram_fits <-
      readRDS(here::here("results", results_name, "ram_results.Rds"))
  }
  
  # assess performance ------------------------------------------------------
  
  
  fit_failed <- !(ram_fits$fits %>% map("error") %>% map_lgl(is.null))
  
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
    select(-data, -priors, -index_type) %>%
    unnest(cols = "error")
  
  ram_fits <- ram_fits %>%
    mutate(
      summary_plot = map(fit, "summary_plot"),
      summary = map(fit, "fit_summary")
    )
  browser()
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
    select(-data,
           -summary_plot,
           -fit,
           -summary_plot,
           -priors) %>%
    unnest(cols = summary) %>%
    filter(!is.na(rmse), is.finite(rmse)) %>%
    left_join(ram_error, by = "experiment") %>%
    group_by(stockid, metric) %>%
    mutate(baseline_rmse = mean(rmse)) %>%
    mutate(delta_rmse = rmse - baseline_rmse) %>%
    ungroup()
  
  # why infinite RMSE
  
  
  fits$initial_state_type <- forcats::fct_relevel(fits$initial_state_type,"heuristic")
  
  fits$u_window[fits$use_u_priors == FALSE] <- "none"
  
  fits$u_window <- forcats::fct_relevel(fits$u_window,"none")
  
  
  wtf <- fits %>% 
    filter(use_terminal_state == TRUE,
           metric == "b_v_bmsy",
           b_ref_type == "b") %>% 
    arrange(desc(rmse))
  
  View(wtf)
  
  i <- 2
  
  huh <- read_rds(file.path(results_path, "experiments", paste0("experiment-", wtf$experiment[i],".rds")))
  
  driors <- ram_fits$priors[ram_fits$experiment == wtf$experiment[i]][[1]]
  
  data <- ram_fits$data[ram_fits$experiment == wtf$experiment[i]][[1]]
  
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
    facet_wrap(use_index~metric)
  
  obs_v_pred_plot <- noerror_data %>%
    filter(metric == "depletion") %>%
    ggplot(aes(observed, predicted, color = fit_name)) +
    geom_point(alpha = 0.5) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    geom_smooth(method = 'lm', aes(fill = fit_name)) +
    facet_wrap(use_index~metric)
  
  
  
  overall_rmse_plot <- noerror_data %>%
    filter(is.finite(rmse)) %>%
    ggplot(aes(rmse, metric, fill = fit_name)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0,2))
  
  
  overall_rmse_plot <- noerror_data %>%
    filter(is.finite(rmse)) %>%
    ggplot(aes(rmse, metric, fill = u_window)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0,2)) +
    facet_wrap(~use_index)
  
  
  overall_bias_plot <- noerror_data %>%
    ggplot(aes(bias, metric, fill = fit_name)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(-1,1))
  
  
  u_effect_plot <- noerror_data %>%
    ggplot(aes(rmse, use_u_priors)) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(0,1)) +
    facet_wrap(~metric)
  
  relative_u_effect_plot <- noerror_data %>%
    ggplot(aes(delta_rmse, use_u_priors)) +
    geom_vline(aes(xintercept = 0), color = "red", linetype = 2) +
    geom_density_ridges(alpha = 0.5) +
    scale_x_continuous(limits = c(-1,1)) +
    facet_wrap(~metric)
  
  
  rmse_plot <- noerror_data %>%
    ggplot(aes(rmse, fit_name, fill = 0.5 - abs(0.5 - ..ecdf..))) +
    stat_density_ridges(alpha = 0.75,
                        geom = "density_ridges_gradient",
                        calc_ecdf = TRUE) +
    scale_fill_viridis(name = "Tail probability", direction = -1) +
    scale_x_continuous(limit = c(0, 1)) +
    facet_wrap( ~ metric)
  
  
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
  browser()
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
      rmse ~  initial_state_type  + use_index*u_window + estimate_shape + estimate_proc_error + use_terminal_state + (1|stockid),
      data = voi_data %>% filter(metric == "b_v_bmsy"),
      cores = 4,
      family = Gamma(link = "log")
    )
  
  
  b_voi_plot <- bayesplot::mcmc_areas(
    as.array(b_voi_fit),
    regex_pars = c("use", "state", "estimate","index","u_"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x) exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in RMSE")
  
  u_voi_fit <-
    rstanarm::stan_glmer(
      rmse ~  initial_state_type  + use_index  + use_u_priors+ u_window + estimate_shape + use_terminal_state + (1|stockid),
      data = voi_data %>% filter(metric == "u_v_umsy", is.finite(rmse)),
      cores = 4,
      family = Gamma(link = "log")
    )
  
  u_voi_plot <- bayesplot::mcmc_areas(
    as.array(u_voi_fit),
    regex_pars = c("use", "state", "estimate"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x) exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in RMSE")
  
  
  
  index_voi_fit <-
    rstanarm::stan_glmer(
      rmse ~  factor(index_freq) + factor(index_window) + (1|stockid),
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
    transformations = function(x) exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in RMSE")
  
  
  
  # effect on accuracy
  
  acc_data <- voi_data %>%
    filter(metric == "b_v_bmsy") %>%
    mutate(bin_acc = accuracy > 0.5)
  
  
  acc_voi_fit <-
    rstanarm::stan_glmer(
      bin_acc ~  use_initial_state + use_index + use_u_priors+ use_terminal_u + estimate_shape + (1|stockid) ,
      data = acc_data,
      cores = 4,
      family = binomial()
    )
  
  
  acc_voi_plot <- bayesplot::mcmc_areas(
    as.array(acc_voi_fit),
    regex_pars = c("use", "state", "estimate"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x) exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in Accuracy")
  
  acc_index_voi_fit <-
    rstanarm::stan_glmer(
      bin_acc ~  factor(index_freq) + factor(index_window) + use_initial_state + use_heuristics + use_u_priors+ use_terminal_u + estimate_shape + (1|stockid) ,
      data = acc_data %>% filter(use_index == TRUE),
      cores = 4,
      family = binomial()
    )
  
  
  acc_index_voi_plot <- bayesplot::mcmc_areas(
    as.array(acc_index_voi_fit),
    regex_pars = c("use", "state", "estimate","index"),
    prob = 0.8,
    # 80% intervals
    prob_outer = 0.95,
    # 99%
    point_est = "mean",
    transformations = function(x) exp(x) - 1
  ) +
    scale_x_percent(name = "% Change in Accuracy")
  
  # how wrong can you be?
  
  quality_data <- fits %>%
    filter(use_index == TRUE, fit_name == "sraplus",
           error_cv > min(error_cv))
  
  quality_obs_v_pred_plot <- quality_data %>%
    filter(metric == "b_v_bmsy") %>%
    ggplot(aes(observed, predicted, color = pmin(20,index_rmse))) +
    geom_point(alpha = 0.75, size = 4) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    geom_smooth(method = 'lm') +
    facet_wrap(~metric) +
    scale_color_viridis()
  
  quality_effect_plot <- quality_data %>%
    filter(metric == "b_v_bmsy") %>%
    ggplot(aes(index_rmse, rmse)) +
    geom_point() +
    geom_abline(aes(slope = 0, intercept = mean(rmse))) +
    geom_smooth(method = "lm") +
    scale_x_log10() +
    scale_y_log10()
  
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


support_data <- list(
  effort_data = effort_data,
  effort_region_to_country = effort_region_to_country
)



total_nominal_effort <- effort_data %>%
  filter(effort_type == "nominal") %>%
  group_by(region,year) %>%
  summarise(total_effort = sum(effort)) %>%
  ungroup() %>%
  mutate(region = tolower(region))



total_effecive_effort <- effort_data %>%
  filter(effort_type == "effective") %>%
  group_by(region,year) %>%
  summarise(total_effort = sum(effort)) %>%
  ungroup() %>%
  mutate(region = tolower(region))


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
  mutate(lifetime_catch = sum(total_catch),
         min_catch = min(total_catch),
         catch_length = length(total_catch)) %>%
  ungroup() %>%
  filter(lifetime_catch > 1e5, min_catch > 1e3, catch_length >= 20) %>%
  mutate(continent = tolower(continent)) %>%
  left_join(total_nominal_effort, by = c("year", "continent" = "region"))

fao %>% 
  group_by(year) %>% 
  summarise(catch = sum(capture, na.rm = TRUE)) %>% 
  ggplot(aes(year, catch)) + 
  geom_line() + 
  scale_y_continuous(labels = comma)


total_stocks %>%
  ggplot(aes(year, total_catch, color = scientific_name)) +
  geom_line(show.legend = FALSE) +
  facet_wrap(~continent) +
  scale_y_log10()

# future::plan(future::multiprocess, workers = 4)


if (run_continent_examples == TRUE) {
  set.seed(42)
  
  continent_fits <- total_stocks %>%
    group_by(scientific_name, continent) %>%
    nest() %>%
    ungroup()
  
  sfe <- safely(fit_continent_examples)
  
  continent_fits = continent_fits %>%
    # filter(scientific_name == "Lateolabrax japonicus") %>%
    mutate(ex_fits = pmap(
      list(
        scientific_name = scientific_name,
        data = data
      ),
      sfe,
      model = "sraplus_tmb",
      engine = "stan",
      estimate_proc_error = TRUE,
      estimate_qslope = FALSE,
      estimate_shape = FALSE
    ))
  
  
  # future::plan(future::sequential, workers = 1)
  
  write_rds(continent_fits, path =  file.path(results_path,"continent-fits.rds"))
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
  facet_wrap(~metric, scales = "free_y") +
  labs(x = "Year", y = "") +
  scale_color_discrete(name = "Region") +
  theme_minimal() +
  scale_y_comma()


total_stocks %>%
  group_by(year, continent, scientific_name) %>%
  summarise(`Total Catch` = sum(total_catch)) %>%
  left_join(total_effecive_effort, by = c("continent" = "region", "year")) %>%
  mutate(CPUE = `Total Catch` / total_effort) %>%
  group_by(scientific_name, continent) %>%
  mutate(CPUE = scale(CPUE),
         lifetime_catch = sum(`Total Catch`)) %>%
  rename(Effort = total_effort) %>%
  gather(metric, value, `Total Catch`:CPUE) %>%
  filter(metric == "CPUE") %>%
  ggplot(aes(year, value, color = scientific_name, alpha = lifetime_catch)) +
  geom_line(size = 1, show.legend = FALSE) +
  facet_wrap(~continent, scales = "free_y") +
  labs(x = "Year", y = "Centered and Scaled Effective Abundance Index") +
  scale_color_discrete(name = "Region") +
  theme_minimal() +
  scale_y_comma()



a <- continent_fits %>%
  mutate(lifetime_catch = map_dbl(data,~unique(.x$lifetime_catch))) %>%
  select(scientific_name, ex_fits, continent,lifetime_catch) %>%
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
  select(scientific_name, year, variable, mean, data, continent, lifetime_catch) %>%
  spread(variable, mean) %>%
  ggplot(aes(pmin(4, b_div_bmsy), pmin(4, u_div_umsy), color = continent)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  geom_density2d(alpha = 0.5) +
  # geom_hex(binwidth = c(0.25,0.25)) +
  geom_point(
    aes(size = lifetime_catch),
    show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 4), name = "B/Bmsy") +
  scale_y_continuous(limits = c(0, 4), name = "U/Umsy") +
  facet_grid(continent ~data) +
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
  geom_density_ridges(aes(
    x = pmin(4, mean),
    y = continent,
    fill = continent
  ), show.legend = FALSE,
  alpha = 0.75) +
  facet_grid(Region ~ variable) +
  scale_x_continuous(limits = c(0, NA))

continent_kobe_plot



# run east-india analysis -----------------------------------------------------

fao_names <- fao %>%
  select(common_name, scientific_name) %>%
  unique()


ei_data <- ei_capture %>%
  mutate(capture = dplyr::case_when(
    is.na(capture) ~ 0,
    TRUE ~ capture
  )) %>%
  group_by(species_asfis_species) %>%
  mutate(
    first_year = year[capture > 0][1],
    total_capture = sum(capture)
  ) %>%
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
    # s(4) %>%
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
        chains = 4,
        cores = 4,
        n_keep = 5000
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
  facet_wrap(~data)

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

fmi$country <- countrycode::countrycode(fmi$country_rfmo, "country.name","un.name.en")



fmi_country_species <- fmi %>%
  filter(!is.na(country), !is.na(scientificname)) %>%
  select(country, scientificname) %>%
  unique() %>%
  rename(scientific_name = scientificname)


fao_country_species_regions <- fao %>%
  group_by(country, scientific_name,fao_area, fao_area_code) %>%
  summarise(tc = sum(capture, na.rm = TRUE)) %>%
  group_by(country) %>%
  mutate(ptc = percent_rank(tc)) %>%
  filter(ptc > 0.25) %>%
  ungroup()

fmi_country_species_region <- fmi_country_species %>%
  left_join(fao_country_species_regions, by = c("country", "scientific_name")) %>%
  filter(!is.na(fao_area))

mean_regional_isscaap_fmi <- fmi %>%
  left_join(fmi_country_species_region, by = c("country_rfmo" ="country", "scientificname" = "scientific_name")) %>%
  select(lookup_code, fao_area_code, research, management, socioeconomics, enforcement, isscaap_group, tc) %>%
  na.omit() %>%
  pivot_longer(cols = c("research", "management", "socioeconomics", "enforcement"),names_to = "metric", values_to = "value") %>%
  group_by(fao_area_code,isscaap_group, metric) %>%
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
  filter(!is.na(status)) %>% 
  mutate(fao_area_code = as.numeric(bad_fao_code)) %>% 
  mutate(stock_complex = paste(common_name, isscaap_number, fao_area_code, sep = '_'))


# merge with fao data

fao_status <- fao %>% 
  group_by(year, fao_area, fao_area_code, isscaap_group, isscaap_number, common_name, 
           scientific_name) %>% 
  summarise(catch = sum(capture), .groups = "drop") %>% 
  left_join(sofia_status, by = c("fao_area_code", "isscaap_number", "common_name", "scientific_name")) %>% 
  filter(!is.na(status)) %>% 
  mutate(stockid = paste(common_name, isscaap_number, fao_area_code, sep = '_')) %>% 
  group_by(stockid) %>%
  mutate(
    first_year = year[catch > 0 & !is.na(catch)][1],
    last_year = last(year[catch >= 0 & !is.na(catch)])
  ) %>%
  filter(year >= first_year, year <= last_year) %>%
  mutate(interp_catch = zoo::na.approx(catch, na.rm = FALSE)) %>%
  mutate(catch = dplyr::case_when(
    is.na(catch) ~ interp_catch,
    TRUE ~ catch
  )) %>%
  filter(!is.na(catch)) %>%
  select(-first_year, -last_year, -interp_catch) %>%
  filter(!is.na(fao_area_code),
         !is.na(common_name))


n_status_plot <- fao_status %>%
  group_by(year, status) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  ggplot() +
  geom_area(aes(year, n, fill = status), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao2011$year), max(fao2011$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
# made it here

huh <- fao_status %>% 
  filter(fao_a)

fao2011 <-
  readxl::read_xlsx(here("data", "FAO SOFIA Stock status table.xlsx"), sheet = "Status2011") %>%
  janitor::clean_names() %>%
  mutate(stockid = paste(name, area, sp_group, sep = "_")) %>%
  select(stockid, everything()) %>%
  filter(!is.na(short)) %>%
  gather(year, catch, contains("x")) %>%
  mutate(year = as.numeric(str_replace(year, "\\D", ""))) %>%
  arrange(stockid, year) %>%
  group_by(stockid) %>%
  mutate(
    first_year = year[catch > 0 & !is.na(catch)][1],
    last_year = last(year[catch >= 0 & !is.na(catch)])
  ) %>%
  filter(year >= first_year, year <= last_year) %>%
  mutate(interp_catch = zoo::na.approx(catch, na.rm = FALSE)) %>%
  mutate(catch = dplyr::case_when(
    is.na(catch) ~ interp_catch,
    TRUE ~ catch
  )) %>%
  filter(!is.na(catch)) %>%
  select(-first_year, -last_year, -interp_catch) %>%
  mutate(
    year = as.numeric(year),
    area = as.numeric(area)
  ) %>%
  filter(!is.na(area),
         !is.na(name))


n_status_plot <- fao2011 %>%
  group_by(year, short) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  ggplot() +
  geom_area(aes(year, n, fill = short), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao2011$year), max(fao2011$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

n_region_status_plot <- fao2011 %>%
  group_by(year,area, short) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  group_by(year, area) %>%
  mutate(pn = n / sum(n)) %>%
  ggplot() +
  geom_area(aes(year, pn, fill = short), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao2011$year), max(fao2011$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~area)

catch_status_plot <- fao2011 %>%
  group_by(year, short) %>%
  summarise(n = length(catch), total_catch = sum(catch)) %>%
  ggplot() +
  geom_area(aes(year, total_catch, fill = short), color = "black", alpha = 0.75) +
  labs(x = "Year") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao2011$year), max(fao2011$year)), expand = c(0, 0)) +
  scale_y_comma(name = "Total Capture (mt)")


# load sar data

load(here::here("data","PNAS-data","TBPdata.Rdata"))

extra_sar <- readr::read_csv(here("data","sar.csv")) %>%
  mutate(country = countrycode::countrycode(Country, "country.name",'un.name.en')) %>%
  left_join(fao_country_region_key, by = "country") %>%
  filter(ptc > 0.5) %>%
  mutate(AreaKm2_0_200 = tc) %>%
  select(fao_area_code,SAR,AreaKm2_0_200)

regional_sar <- TBPdata@SummaryByRegion

sar_to_fao_region <- read_csv(file = here("data","sar_region_lookup-maually-created.csv"))

approx_sar_by_fao_region <- regional_sar %>%
  right_join(sar_to_fao_region, by = "RegionName") %>%
  select(fao_area_code, SAR,AreaKm2_0_200) %>%
  bind_rows(extra_sar) %>%
  group_by(fao_area_code) %>%
  summarise(sar = weighted.mean(SAR, w = AreaKm2_0_200, na.rm = TRUE))

# load in ram status

load(here("data", "RAMCore.rdata"))

ram_to_region <- data_frame(stockid = colnames(UvT), region = Region)

u_v_t <- UvT %>%
  as.data.frame() %>%
  mutate(year = rownames(.)) %>%
  gather(stockid, u_umsy, -year) %>%
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
  gather(stockid, b_bmsy, -year)

status_t <- u_v_t %>%
  left_join(b_v_t, by = c("year", "stockid"))

breaks <- c(0.8, 1.2)

breaks <- c(0, breaks, Inf)

labels <- c("over", "fully", "under")

status_t <- status_t %>%
  mutate(bin = cut(b_bmsy, breaks = breaks, labels = labels)) %>%
  mutate(year = as.numeric(year))

# add state space u estimates to FAO

fao_area_to_ram <- readxl::read_xlsx(here("data", "FAO SOFIA Stock status table.xlsx"), sheet = "cleaned_table") %>%
  janitor::clean_names() %>%
  mutate(
    area = as.numeric(fao_area)
  ) %>%
  mutate(tempid = 1:nrow(.)) %>%
  select(-sar,-fmi)


state_space_results <- readxl::read_xlsx(here("data", "state_space_results_08_28_2018 by region.xlsx"), sheet = "state_space_results_08_28_2018") %>%
  janitor::clean_names()


fao2011 <- fao2011 %>%
  left_join(fao_area_to_ram, by = "area") %>%
  left_join(state_space_results %>% filter(variable == "UvU"), by = c("year", "state_space_region" = "region")) %>%
  left_join(approx_sar_by_fao_region, by = c("area" = "fao_area_code")) %>%
  rename(mean_u_umsy = dlm_geomean)

fao_ram_u_plot <- fao2011 %>%
  ggplot(aes(
    year,
    (mean_u_umsy)
  )) +
  geom_point() +
  facet_wrap(~state_space_region, scales = "free_y")


ram_status_plot <- status_t %>%
  group_by(year, bin) %>%
  summarise(n = length(b_bmsy)) %>%
  ggplot() +
  geom_area(aes(year, n, fill = bin), color = "black", alpha = 0.75) +
  labs(x = "Year", y = "Number of Stocks") +
  scale_fill_viridis_d() +
  scale_x_continuous(limits = c(min(fao2011$year), max(fao2011$year)), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

# link up generic effort to fao 2011


rous_data <- read.csv(here("data", "MappedFAO.csv")) %>% 
  na.omit() %>% 
  as_tibble() %>% 
  janitor::clean_names() %>% 
  mutate(country = countrycode::countrycode(iso3, "iso3c", "un.name.en")) %>% 
  filter(type2 == "I") %>% 
  select(year,fao, effort_cell_reported_nom, country ) %>% 
  rename(area = fao)

rough_fao_region_effort <- effort_data %>%
  filter(effort_type == "nominal") %>%
  right_join(effort_region_to_country, by = "region" ) %>%
  right_join(fao_country_region_key, by = "country") %>%
  group_by(fao_area_code, year) %>%
  summarise(effort_index = sum(effort)) %>%
  ungroup()



temp_fao <- fao2011 %>%
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

    # fao_catch <- fao

    # effort <- rous_data

    # fao_stock <- temp_fao$stockid[which(temp_fao$stockid == huh)]

    comm_name <-
      str_split(fao_stock, pattern = '_')[[1]][1] %>% str_remove_all("(\\d)|(-)")
    
    isscp_number <-
      as.numeric(str_split(fao_stock, pattern = '_')[[1]][3])
    
    fao_code <- as.numeric(str_split(fao_stock, pattern = '_')[[1]][2])
    
    fao_matches <- fao %>% {
      if (any(.$common_name == comm_name,na.rm = TRUE)) {
        filter(., common_name == comm_name & fao_area_code == fao_code)
        
      } else {
        filter(.,
               isscaap_number == isscp_number & fao_area_code == fao_code)
        
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
      left_join(matched_effort, by = c("year", "area"))
    
    
  }


temp_fao <- temp_fao %>%
  mutate(
    data = map2(
      stockid,
      data,
      assign_effort,
      fao_catch = fao,
      effort = rous_data
    )
  )


fao2011 <- temp_fao %>% 
  unnest(cols = data)


# fao2011 %>% 
#   group_by(year, area) %>% 
#   summarise(effort_index = unique(effort_index)) %>% 
#   ggplot(aes(year, effort_index)) + 
#   geom_line() + 
#   facet_wrap(~area)

# temp <- fao2011 %>%
#   # filter(area %in% 67) %>%
#   group_by(stockid) %>%
#   nest() %>%
#   ungroup()
# 
# huh <- temp$stockid[4]
# 
# wtf <- fao2011 %>% 
#   filter(stockid == huh)


fao2011 <- fao2011 %>%
  group_by(stockid) %>%
  mutate(lifetime_catch = sum(catch, na.rm = TRUE),
         min_catch = min(catch, na.rm = TRUE),
         nas = sum(is.na(catch))) %>%
  ungroup() %>%
  filter(lifetime_catch > 50000,
         min_catch > 0,
         nas == 0,
         !is.na(species)) %>%
  ungroup() %>%
  mutate(sp_group = as.numeric(sp_group)) %>%
  left_join(fao_species %>% select(isscaap_group, isscaap_number) %>% unique(), by = c("sp_group" = "isscaap_number"))


support_data <- list(mean_regional_isscaap_fmi = mean_regional_isscaap_fmi)

areas <- c(67, 57, 37,71)

areas <- unique(fao2011$area)
# annnnnd try and run assessments

if (run_sofia_comparison == TRUE) {
  # future::plan(future::multiprocess, workers = 4)
  
  set.seed(42)
  
  future::plan(multisession,workers = n_cores)
  
  fao2011_fits <- fao2011 %>%
    # filter(area %in% 67) %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup() %>%
    # sample_n(25) %>%
    mutate(
      fits = future_map(
        data,
        safely(fit_fao),
        support_data = support_data,
        min_effort_year = 1960,
        engine = "stan",
        cores = 1,
        .progress = TRUE,
        .options = future_options(globals = support_data, packages = c("tidyverse","sraplus"))
      )
    )
  

  # "Pink(=Humpback)salmon_67_23"
  # 
  # temp <- fao2011 %>%
  #   group_by(stockid) %>%
  #   nest() %>%
  #   ungroup() %>%
  #   # sample_n(1) %>%
  #   # filter(stockid ==  "Pink(=Humpback)salmon_67_23") %>%
  #   slice(141) %>%
  #   mutate(
  #     fits = map(
  #       data,
  #       (fit_fao),
  #       support_data = support_data,
  #       min_effort_year = 1960,
  #       engine = "stan",
  #       cores = 1
  #     )
  #   )

  
  write_rds(fao2011_fits, path = file.path(results_path, "fao2011-fits.rds"))
  
} else {
  fao2011_fits <-   read_rds(path = file.path(results_path, "fao2011-fits.rds"))
  

}
reserve <- fao2011_fits

ugh <- map_dbl(fao2011$data, ~unique(.x$area)) == 67
#
fao_worked <- map(fao2011_fits$fits, "error")


fao_worked <- map(fao2011_fits$fits, "error") %>% map_lgl(is.null)

fao2011_fits <- fao2011_fits %>%
  filter(fao_worked)

fao2011_fits <- fao2011_fits %>%
  mutate(fits = map(fits, "result"))

faofits <- fao2011_fits %>%
  select(-data) %>%
  unnest(cols = fits) %>%
  filter(variable == "b_div_bmsy") %>%
  mutate(bin = cut(mean, breaks = breaks, labels = labels)) %>%
  separate(stockid,c("name","area","speciesgroup"), sep = '_')


# wtf <- faofits %>% 
#   filter(stockid == stockid[1],
#          str_detect(data, "cpue")) %>% 
#   ggplot(aes(year, mean, color = data)) + 
#   geom_line() + 
#   facet_wrap(~data)

fao2011_data <- fao2011 %>%
  mutate(bin = dplyr::case_when(
    short == "F" ~ "fully",
    short == "U" ~ "under",
    short == "O" ~ "over",
    TRUE ~ "unknown"
  )) %>%
  mutate(
    data = "FAO Report 569",
    mean = rep(1, nrow(.))
  ) %>%
  select(area, year, bin, data, mean) %>%
  mutate(area = as.numeric(as.character(area)))


fao_regions_key <- read_csv(here("data", "fao_regions_key.csv"))

temp <- faofits %>%
  mutate(area = as.numeric(area)) %>%
  select(area, year, bin, data, mean) %>%
  bind_rows(fao2011_data) %>%
  filter(
    !is.na(area), !area %in% c(88, 48),
    !bin == "unknown"
  ) %>%
  left_join(fao_regions_key, by = c("area" = "fao_region_num")) %>%
  group_by(fao_region_name,area, year, bin, data) %>%
  summarise(n = length(mean)) %>%
  group_by(fao_region_name,area, data, year) %>%
  mutate(pn = round(n / sum(n), 2)) %>%
  group_by(fao_region_name, area,data) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  tidyr::complete(bin, nesting(fao_region_name, area, year, data),
                  fill = list(n = 0, pn = 0)
  ) %>%
  filter(bin == "over") %>%
  mutate(data = case_when(data == "u_umsy" ~ "ram_u_umsy", TRUE ~ data)) %>%
  mutate(data = fct_reorder(data, pn))

# temp$pn[temp$fao_region_name == "Pacific, Northeast" & temp$bin == "over" & temp$data == "FAO Report 569"] <- 0.09


fao2011_plot <- temp %>%
  ggplot(aes(data, pn, fill = fao_region_name)) +
  geom_col(size = 1, alpha = 1, size = 1, show.legend = FALSE, color = "black") +
  coord_flip() +
  facet_grid(fao_region_name ~ .) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), labels = paste0(seq(0, 100, by = 25), "%"), name = "Percent of Stocks Overfished") +
  theme(
    axis.text.x = element_text(size = 8),
    axis.title.y = element_blank()
  ) +
  scale_fill_discrete(name = "FAO Category")

fao2011_plot


if (!dir.exists(here("data","FAO_AREAS_NOCOASTLINE"))){

download.file(url = "http://www.fao.org/figis/geoserver/area/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=area:FAO_AREAS_NOCOASTLINE&outputFormat=SHAPE-ZIP",
              destfile = here("data","FAO_AREAS_NOCOASTLINE.zip"))

unzip(here("data","FAO_AREAS_NOCOASTLINE.zip"), exdir = here("data","FAO_AREAS_NOCOASTLINE"))

}

fao_areas <- sf::st_read(here('data',"FAO_AREAS_NOCOASTLINE")) %>% 
  janitor::clean_names() 

fao_areas <- fao_areas %>% 
  group_by(f_area) %>% 
  nest() %>% 
  mutate(geometry = map(data,st_union)) %>%
  select(-data)

fao_areas = fao_areas %>% 
  unnest(cols = geometry) %>% 
  ungroup() %>% 
  sf::st_as_sf() %>% 
  sf::st_simplify()

fao_area_status <- fao_areas %>% 
  left_join(temp %>% mutate(f_area = as.character(area)), by = "f_area") %>% 
  filter(!is.na(data))


world_map <- rnaturalearth::ne_countries(returnclass = "sf", scale = "small") #%>% 
  # sf::st_union() %>% 
  # sf::st_polygonize()

ggplot() + 
  geom_sf(data = world_map, size = .01, fill = "lightgrey", color = "lightgrey")


po_map_plot <- fao_area_status %>% 
  filter(!is.na(data)) %>% 
  ggplot() + 
  geom_sf(aes(fill = pn)) +
  geom_sf(data = world_map, fill = "lightgrey", color = "black", size = 0.05) +
  facet_wrap(~data) + 
  scale_fill_viridis_b(option = "magma", labels = percent, name = "% Overfished",
                       guide = guide_coloursteps(barwidth = ggplot2::unit(8, "lines"),
                                  reverse = FALSE,
                                  axis.colour = "white",
                                  axis.linewidth = 1,
                                  show.limits = TRUE)) + 
  ggthemes::theme_base() +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "white"))

po_map_plot

ggsave("percent_overfished_by_region.pdf",po_map_plot)
# ----save-plots----------------------------------------------------------

plots <- ls()[str_detect(ls(), "_plot")]

save(file = here("results", results_name, "paper_plots.Rdata"), list = plots)
