

# setup -------------------------------------------------------------------

library(tidyverse)
library(ggridges)
library(ggtext)
library(gganimate)
library(rstan)
library(mvtnorm)
library(FishLife)
library(furrr)
library(doParallel)
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
library(portedcmsy)
extrafont::loadfonts()

Sys.unsetenv("PKG_CXXFLAGS")

sraplus::get_tmb_model()

options(dplyr.summarise.inform = FALSE)

rstan::rstan_options(auto_write = TRUE)

# options -----------------------------------------------------------------

min_years_catch <- 20

crazy_b <- 5 # threshold for suspect B/Bmsy value

crazy_u <- 8 # threshold for suspect U/Umsy value

draws <- 3000

min_draws <- 2000 # minimum number of unique SIR draws

n_cores <- 2
# number of cores for parallel processing

# options(mc.cores = 1)

results_name <- "v1.0"


results_description <-
  "publication version of results"

run_voi_models <- TRUE
# sub options for run_voi_models
fit_models <- FALSE

write_results <- FALSE

process_fits <- FALSE

run_case_studies <- FALSE

run_sofia_comparison <- FALSE

run_ram_tests <- FALSE

run_ram_comparison <- FALSE

knit_paper <- FALSE

warning("Running full analysis takes upwards of 24 hours on 2 cores. Recommend starting on a Friday night and then having a nice weekend. Given memory constraints of models, more than 2 cores is not recommended (memory routinely runs out on a 16 core 36GB machine when cores > 2)")

engine <-  "stan"

catchability = 1e-2

pub_theme <- theme_ipsum(base_size = 10,
                         axis_text_size = 10,
                         axis_title_size = 12) +
    theme(panel.spacing = unit(.5,"lines"))

theme_set(pub_theme)


# theme_set(theme_classic() + theme(strip.background = element_rect(color = "transparent")))
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


prepare_sofia_data(lookup_fmi_names = TRUE)


has_total_biomass = ram_data %>%
  group_by(stockid) %>%
  summarise(hastb = !any(is.na(total_biomass)),
            hasb = !any(is.na(b_v_bmsy))) %>%
  filter(hastb & hasb)

og_ram_data <- ram_data

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
    error_cv = sample(c(0), draws, replace = TRUE),
    stockid = sample(unique(ram_data$stockid), draws, replace = TRUE),
    use_u_priors =  sample(c(TRUE, FALSE),  draws, replace = TRUE),
    u_cv = runif(draws, 0.05, 0.2),
    index_window = sample(c(1, 0.25, .5), draws, replace = TRUE),
    index_freq = sample(c(1), draws, replace = TRUE),
    u_window = sample(c("snapshot", "recent", "complete"), draws, replace = TRUE),
    b_ref_type = sample(c("k", "b"),  draws, replace = TRUE),
    f_ref_type = sample(c("f", "fmsy"),  draws, replace = TRUE),
    estimate_shape = sample(c(FALSE, TRUE), draws, replace = TRUE),
    shape_prior_source = sample(c("thorson", "fishlife"), draws, replace = TRUE),
    estimate_proc_error = sample(c(FALSE, TRUE), draws, replace = TRUE)
  ) %>%
    mutate(error_cv = map_dbl(error_cv, ~ ifelse(.x == 0, 1e-6, runif(1, 0.05, 0.2)))) %>%
    mutate(b_ref_type = ifelse(initial_state_type == "prior",
                               "b",
                               b_ref_type))

  
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
          u_window = u_window,
          shape_prior_source = shape_prior_source
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
    
    future::plan(multiprocess, workers = n_cores)
    
    
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
      # sample_n(4) %>%
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
        assessments = c("sraplus"),
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
  

  # View(wtf)
  
  i <- 2
  
  # huh <-
  #   read_rds(file.path(
  #     results_path,
  #     "experiments",
  #     paste0("experiment-", wtf$experiment[i], ".rds")
  #   ))
  # 
  # driors <-
  #   ram_fits$priors[ram_fits$experiment == wtf$experiment[i]][[1]]
  # 
  # data <-
  #   ram_fits$data[ram_fits$experiment == wtf$experiment[i]][[1]]
  # 
  # plot_sraplus(huh$sirplus)
  
  # plot_prior_posterior(huh$sraplus, driors)
  
  
  # plot diagnostics --------------------------------------------------------
  noerror_data <- fits %>%
    filter(error_cv == min(error_cv),
           fit_name != 'cmsy') %>%
    mutate(index = index_window * use_index,
           shape = case_when(estimate_shape == FALSE ~  "Fixed", 
                             TRUE ~ shape_prior_source)) %>%
    mutate(
      index = case_when(
        index == 0 ~ "No Abundance Index",
        index == 0.25 ~ "Abundance Index (1/4)",
        index == 0.5 ~ "Abundance Index (1/2)",
        index == 1 ~ "Abundance Index"
      )
    ) %>% 
    mutate(index = fct_relevel(index, "No Abundance Index")) %>% 
    mutate(shape = fct_relevel(shape, "Fixed"))
  
  
  
  voi_data <- fits %>%
    filter(fit_name != "cmsy")
  
  
  # voi_data %>%
  #   filter(use_index == FALSE,
  #          metric == "b_v_bmsy",
  #          b_ref_type == "b") %>%
  #   ggplot(aes(rmse, fill = use_terminal_state)) +
  #   geom_density(position = "dodge")
  
  
  
  # noerror_data %>%
  #   filter(use_index == FALSE,
  #          metric == "b_v_bmsy",
  #          b_ref_type == "b") %>%
  #   ggplot(aes(observed, predicted, color = use_terminal_state)) +
  #   geom_point()
  
  obs_v_pred_plot <- noerror_data %>%
    filter(metric == "b_v_bmsy") %>%
    ggplot(aes(observed, predicted, color = fit_name)) +
    geom_point(alpha = 0.25) +
    geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
    geom_smooth(method = 'lm', aes(fill = fit_name)) +
    facet_wrap(use_index ~ metric)
  
  
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
  
  # value of information analyses
  
  # how does the comparison stack up when you have more data?
  rmse_bias_plot <- noerror_data %>%
    ggplot(aes(rmse, bias)) +
    geom_hex(aes(fill = ..density..), binwidth = c(0.1, 0.1)) +
    facet_grid(metric ~ fit_name) +
    scale_x_continuous(limits = c(0, .5)) +
    scale_y_continuous(limits = c(-.5, .5))  +
    scale_fill_viridis(option = "A")

  # b_rmse_fit <-
  #   rstanarm::stan_glmer(
  #     log(rmse) ~ fit_name + (1 | stockid),
  #     data = noerror_data %>% filter(metric == "b_v_bmsy"),
  #     cores = 4
  #   )
  
  # b_rmse_fit_posterior <- as.array(b_rmse_fit)
  
  # b_rmse_plot <- bayesplot::mcmc_areas(
  #   b_rmse_fit_posterior,
  #   regex_pars = "fit_name",
  #   prob = 0.8,
  #   # 80% intervals
  #   prob_outer = 0.99,
  #   # 99%
  #   point_est = "mean"
  # ) +
  #   scale_x_percent(name = "Percent change in B/Bmsy RMSE")
  # 
  # 
  # 
  # value of information regressions
  
  b_voi_fit <-
    rstanarm::stan_glm(
      rmse ~  initial_state_type  + index +  u_window + shape + use_terminal_state,
      data = noerror_data %>% filter(metric == "b_v_bmsy"),
      cores = 4,
      family = Gamma(link = "log")
    )
  
  a = b_voi_fit %>%
    tidy_draws() %>%
    select(
      -contains("b["),
      -contains("__"),
      -contains("Sigma"),
      -contains("shape"),
      -contains("(Intercept)"),
      -contains("proc")
    ) %>%
    pivot_longer(-contains("."), names_to = "variable", values_to = "value") %>% 
    mutate(value = exp(value) - 1)
  
  b_voi_plot <-  a %>%
    mutate(
      clean_var = case_when(
        variable == "use_terminal_stateTRUE" ~ "Most Recent B/Bmsy",
        variable == "indexAbundance Index" ~ "Abundance Index (Complete)",
        variable == "indexAbundance Index (1/4)" ~ "Abundance Index (1/4)",
        variable == "indexAbundance Index (1/2)" ~ "Abundance Index (1/2)",
        variable == "u_windowsnapshot" ~ "F/Fmsy (Most Recent)",
        variable == "u_windowrecent" ~ "F/Fmsy (Last Five Years)",
        variable == "u_windowcomplete" ~ "F/Fmsy (Complete Series)",
        variable == "initial_state_typeunfished" ~ "Initial State Unfished" ,
        variable == "initial_state_typeprior" ~ "Initial State Informed by Catch",
        variable == "initial_state_typeknown" ~ "Initial State Known",
        TRUE ~ variable
      )
    ) %>%
    ggplot() +
    geom_vline(aes(xintercept = 0)) +
    ggdist::stat_halfeye(aes(value, clean_var)) +
    scale_x_percent(name = "% Change in RMSE") + 
    scale_y_discrete(name = '')
  
  
  # u_voi_fit <-
  #   rstanarm::stan_glmer(
  #     rmse ~  initial_state_type  + use_index + u_window + estimate_shape + estimate_proc_error + use_terminal_state + (1 |
  #                                                                                                                         stockid),
  #     data = voi_data %>% filter(metric == "u_v_umsy", is.finite(rmse)),
  #     cores = 4,
  #     family = Gamma(link = "log")
  #   )
  
  # u_voi_plot <- bayesplot::mcmc_areas(
  #   as.array(u_voi_fit),
  #   regex_pars = c("use", "state", "estimate", "index", "u_"),
  #   prob = 0.8,
  #   # 80% intervals
  #   prob_outer = 0.95,
  #   # 99%
  #   point_est = "mean",
  #   transformations = function(x)
  #     exp(x) - 1
  # ) +
  #   scale_x_percent(name = "% Change in RMSE")
  # 
  
  # 
  # index_voi_fit <-
  #   rstanarm::stan_glmer(
  #     rmse ~  factor(index_freq) + factor(index_window) + index_rmse + (1 |
  #                                                                         stockid),
  #     data = voi_data %>% filter(metric == "b_v_bmsy", use_index == TRUE, error_cv > min(error_cv)) %>% mutate(index_rmse = scale(index_rmse)),
  #     family = Gamma(link = "log"),
  #     cores = 4
  #   )
  # 
  # 
  # index_voi_plot <- bayesplot::mcmc_areas(
  #   as.array(index_voi_fit),
  #   regex_pars = c("freq", "window","index"),
  #   prob = 0.8,
  #   # 80% intervals
  #   prob_outer = 0.95,
  #   # 99%
  #   point_est = "mean",
  #   transformations = function(x)
  #     exp(x) - 1
  # ) +
  #   scale_x_percent(name = "% Change in RMSE")
  # 
  # 
  # effect on accuracy
  
  # acc_data <- voi_data %>%
  #   filter(metric == "b_v_bmsy") %>%
  #   mutate(bin_acc = accuracy > 0.5)
  # 
  # 
  # acc_voi_fit <-
  #   rstanarm::stan_glmer(
  #     bin_acc ~  initial_state_type  + use_index + u_window + estimate_shape + estimate_proc_error + use_terminal_state + (1 |
  #                                                                                                                            stockid),
  #     data = acc_data,
  #     cores = 4,
  #     family = binomial()
  #   )
  # 
  # 
  # acc_voi_plot <- bayesplot::mcmc_areas(
  #   as.array(acc_voi_fit),
  #   regex_pars = c("use", "state", "estimate", "index", "u_"),
  #   prob = 0.8,
  #   # 80% intervals
  #   prob_outer = 0.95,
  #   # 99%
  #   point_est = "mean",
  #   transformations = function(x)
  #     exp(x) - 1
  # ) +
  #   scale_x_percent(name = "% Change in Accuracy")
  # 
  # acc_index_voi_fit <-
  #   rstanarm::stan_glmer(
  #     bin_acc ~  factor(index_freq) + factor(index_window) + (1 | stockid),
  #     data = acc_data %>% filter(use_index == TRUE),
  #     cores = 4,
  #     family = binomial()
  #   )
  # 
  # 
  # acc_index_voi_plot <- bayesplot::mcmc_areas(
  #   as.array(acc_index_voi_fit),
  #   regex_pars = c("use", "state", "estimate", "index"),
  #   prob = 0.8,
  #   # 80% intervals
  #   prob_outer = 0.95,
  #   # 99%
  #   point_est = "mean",
  #   transformations = function(x)
  #     exp(x) - 1
  # ) +
  #   scale_x_percent(name = "% Change in Accuracy")
  # 
  # # how wrong can you be?
  # 
  # quality_data <- fits %>%
  #   filter(use_index == TRUE,
  #          fit_name == "sraplus",
  #          error_cv > min(error_cv))
  # 
  # quality_obs_v_pred_plot <- quality_data %>%
  #   filter(metric == "b_v_bmsy") %>%
  #   ggplot(aes(observed, predicted, color = pmin(20, index_rmse))) +
  #   geom_point(alpha = 0.75, size = 4) +
  #   geom_abline(aes(slope = 1, intercept = 0), linetype = 2) +
  #   geom_smooth(method = 'lm') +
  #   facet_wrap( ~ metric) +
  #   scale_color_viridis()
  # 
  # quality_effect_plot <- quality_data %>%
  #   filter(metric == "b_v_bmsy") %>%
  #   ggplot(aes(index_rmse, rmse)) +
  #   geom_point() +
  #   geom_abline(aes(slope = 0, intercept = mean(rmse))) +
  #   geom_smooth(method = "lm") +
  #   scale_x_log10() +
  #   scale_y_log10()
  # 
  flist <- ls()[str_detect(ls(), "_fit")]
  
  save(list = flist,
       file = file.path(results_path, "voi_fits.RData"))
  
 
  
} # close fit voi models



# run case study ----------------------------------------------------------


exs <- ram_data %>%
  filter(
    stockid %in% unique(sraplus::fmi$stockid) &
      stockid %in% unique(sraplus::sar$stockid)
  ) %>%
  mutate(fao_area_code = as.numeric(primary_FAOarea),
         country = as.character(primary_country)) %>%
  mutate(country =
           countrycode::countrycode(country, "country.name", "un.name.en"))

comp_stocks <- exs %>%
  select(scientificname) %>%
  unique() %>%
  mutate(fishbase_vuln = map_dbl(
    scientificname,
    ~ rfishbase::species(.x, fields = "Vulnerability")$Vulnerability
  )) %>%
  mutate(
    resilience = dplyr::case_when(
      fishbase_vuln < 33 ~ "Low",
      fishbase_vuln >= 66 ~ "High",
      TRUE ~ "Medium"
    )
  ) %>%
  select(scientificname, resilience)


# exs %>%
#   ggplot(aes(year, b_v_bmsy, color = stockid)) +
#   geom_line(show.legend = FALSE)

total_nominal_effort <- rous_data %>%
  left_join(effort_region_to_country, by = "country") %>%
  rename(fao_area_code = area) %>%
  group_by(fao_area_code, year, country) %>%
  summarise(total_effort = sum(effort_cell_reported_nom, na.rm = TRUE)) %>%
  ungroup()

matched_fmi <- sraplus::fmi %>%
  filter(stockid %in% unique(exs$stockid)) %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup() %>%
  mutate(tmp = map(data, ~ map_df(.x[,c("research", "management", "enforcement","socioeconomics")], mean))) %>%
  select(-data) %>%
  unnest(cols = tmp) %>%
  select(stockid, research, management, enforcement, socioeconomics)


matched_sar <- sraplus::sar %>%
  filter(stockid %in% unique(exs$stockid)) %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup() %>%
  mutate(tmp = map(data, ~ purrrlyr::dmap_at(.x, c("sar"), mean))) %>%
  select(-data) %>%
  unnest(cols = tmp) %>%
  select(stockid, sar)


exs <- exs %>%
  left_join(total_nominal_effort, by = c("fao_area_code", "year", "country")) %>%
  left_join(matched_fmi, by = "stockid") %>%
  left_join(matched_sar, by = "stockid") %>%
  left_join(comp_stocks, by = "scientificname")

if (run_case_studies){
  
  exs <- exs %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup() %>% 
    mutate(fit = map(data,
                     safely(fit_case_studies)))
  
  write_rds(exs, path = file.path(results_path, "raw-case-study-fits.rds"))
  
} else {
  
  exs <- read_rds(path = file.path(results_path, "raw-case-study-fits.rds"))
  
}

fit_worked <-
  exs$fit %>% map("error") %>% map_lgl(is.null)


truth <- exs %>% 
  select(stockid, data) %>% 
  unnest(cols = data) %>% 
  select(stockid, year, catch, b_v_bmsy, u_v_umsy) %>% 
  group_by(stockid) %>% 
  filter(year == max(year)) %>% 
  ungroup()

tmp <- exs %>%
  mutate(lifetime_catch = map_dbl(data,  ~ sum(.x$catch))) %>%
  select(stockid, fit, lifetime_catch) %>%
  filter(fit_worked) %>%
  mutate(fit = map(fit, "result")) %>%
  unnest(cols = fit) %>%
  group_by(stockid) %>%
  filter(year == max(year)) %>%
  filter(variable %in% c("b_div_bmsy", "u_div_umsy")) %>%
  ungroup() %>%
  select(stockid,
         year,
         variable,
         mean,
         data,
         lifetime_catch) %>% 
  pivot_wider(names_from = variable, values_from = mean)


case_study_fits <- tmp %>%
  left_join(
    truth %>% rename(ram_b_v_msy = b_v_bmsy, ram_u_v_umsy = u_v_umsy),
    by = c("stockid", "year")
  )

write_rds(case_study_fits, path = file.path(results_path, "case-study-fits.rds"))

ex_kobe_plot <- tmp %>%
  mutate(data = fct_relevel(data, c("CMSY", "SAR & FMI"))) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  geom_hline(aes(yintercept = 0), linetype = 1) +
  # geom_density2d(alpha = 0.5) +
  geom_hex(data = truth,aes(pmin(4, b_v_bmsy), pmin(4, u_v_umsy)), binwidth = c(0.5,.5), show.legend = FALSE, alpha = 0.9) +
  geom_point(aes(pmin(4, b_div_bmsy), pmin(4, u_div_umsy), color = data,size = lifetime_catch),
             alpha = 0.5,
             show.legend = FALSE) +
  scale_size(range = c(2,7)) +
  scale_x_continuous(limits = c(-.1, 4), name = "B/Bmsy") +
  scale_y_continuous(limits = c(-.1, 4), name = "F/Fmsy") +
  facet_wrap(~ data) +
  scale_fill_continuous(low = "gainsboro", high = "black") +
  ggsci::scale_color_d3(name = "Region") +
  guides(size = FALSE)


tmp2 <-  tmp %>%
  pivot_longer(contains("div"), names_to = "variable", values_to = "value")

truth2 <- truth %>%
  select(stockid, year, contains("_v_")) %>%
  pivot_longer(contains("_v_"), names_to = "variable", values_to = "truth")  %>%
  mutate(variable = str_replace_all(variable, "_v_", "_div_"))


tmp2 <- tmp2 %>%
  left_join(truth2, by = c("stockid", "year", "variable")) %>%
  mutate(
    variable = case_when(
      variable == "b_div_bmsy" ~ "B/Bmsy",
      variable == "u_div_umsy" ~ "F/Fmsy",
      TRUE ~ variable
    )
  )

ex_performance <- tmp2 %>% 
  group_by(variable, data) %>% 
  summarise(r2 = round(yardstick::rsq_vec(truth, value),2)) 


ex_scatter_plot <- tmp2 %>% 
  ggplot(aes(truth, value)) + 
  geom_vline(aes(xintercept = 0)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_abline(aes(slope = 1, intercept = 0),linetype = 2) +
  geom_point(alpha = 0.5, size = 2) + 
  ggtext::geom_richtext(data = ex_performance, aes(x = 1, y = 4, label = paste0("R<sup>2</sup> = ",r2))) +
  facet_grid(variable ~ data, scales = "free_x")  + 
  # scale_size(trans = "sqrt", name = "Lifetime Catch") +
  scale_x_continuous(name = "RLSADB Value", expand = expansion(add = c(0, .1))) + 
  scale_y_continuous("Estimated Value", expand = expansion(add = c(0, .1)))


sraplus_v_truth <- tmp %>% 
  left_join(truth, by = c("stockid","year")) %>% 
  filter(data == "cpue") 

# sraplus_v_truth %>% 
#   ggplot(aes(catch, b_div_bmsy)) + 
#   geom_point(size = 2)


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
# 
# 
comp_stocks <- fao_status %>%
  select(scientific_name) %>%
  unique() %>%
  mutate(fishbase_vuln = map_dbl(
    scientific_name,
    ~ rfishbase::species(.x, fields = "Vulnerability")$Vulnerability
  )) %>%
  mutate(
    resilience = dplyr::case_when(
      fishbase_vuln < 33 ~ "Low",
      fishbase_vuln >= 66 ~ "High",
      TRUE ~ "Medium"
    )
  ) %>% 
  select(scientific_name, resilience)

fao_status <- fao_status %>% 
  left_join(comp_stocks, by = "scientific_name") %>% 
  mutate(resilience = ifelse(is.na(resilience), "Medium",resilience))
if (run_sofia_comparison == TRUE) {
  # future::plan(future::multiprocess, workers = 4)
  
  set.seed(42)
  
  future::plan(multisession, workers = n_cores, .cleanup = TRUE)
  
  fao_status_fits <- fao_status %>%
    # filter(fao_area_code %in% 67) %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup() %>%
    # sample_n(10) %>%
    # slice(10) %>% 
    mutate(
      fits = future_map2(
        data,
        stockid,
        safely(fit_fao),
        support_data = support_data,
        default_initial_state = NA,
        default_initial_state_cv = NA,
        min_effort_year = 1975,
        engine = "stan",
        cores = 1,# run in to lots of memory problems if greater than 1
        cmsy_cores = 1, # run in to lots of memory problems if greater than 1
        write_results = FALSE, 
        results_path = results_path,
        .progress = TRUE,
        .options = future_options(
          globals = support_data,
          packages = c("tidyverse", "sraplus", "portedcmsy")
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
  separate(stockid, c("name", "speciesgroup", "area"), sep = '_', remove = FALSE) %>%
  mutate(fao_area_code = as.numeric(area))


wtf <- faofits %>%
  filter(stockid == unique(stockid)[3]) %>%
  ggplot(aes(year, mean, color = data)) +
  geom_line() +
  facet_wrap(~data)

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

fao_status_data2 <- fao_status %>%
  mutate(
    bin = dplyr::case_when(
      status == "F" ~ "fully",
      status == "U" ~ "under",
      status == "O" ~ "over",
      TRUE ~ "unknown"
    )
  ) %>%
  mutate(data = "FAO Report 569",
         mean = rep(1, nrow(.))) 

faofits2 <- fao_status_fits %>%
  select(-data) %>%
  unnest(cols = fits) %>%
  filter(variable == "b_div_bmsy") %>%
  mutate(bin = cut(mean, breaks = breaks, labels = labels)) %>%
  separate(stockid, c("name", "speciesgroup", "area"), sep = '_', remove = FALSE) %>%
  mutate(fao_area_code = as.numeric(area))


fao_sraplus_comp_data <-  faofits2 %>%
  group_by(name, area) %>%
  filter(year == max(year)) %>%
  ungroup() %>%
  left_join(
    fao_status_data2 %>% rename(fao_bin = bin) %>% select(stockid, year, fao_bin),
    by = c("stockid", "year")
  ) %>%
  mutate(correct = bin == fao_bin)



if (!dir.exists(here("data", "FAO_AREAS_NOCOASTLINE"))) {
  download.file(url = "http://www.fao.org/figis/geoserver/area/ows?service=WFS&request=GetFeature&version=1.0.0&typeName=area:FAO_AREAS_NOCOASTLINE&outputFormat=SHAPE-ZIP",
                destfile = here("data", "FAO_AREAS_NOCOASTLINE.zip"))
  
  unzip(
    here("data", "FAO_AREAS_NOCOASTLINE.zip"),
    exdir = here("data", "FAO_AREAS_NOCOASTLINE")
  )
  
}

write_rds(fao_sraplus_comp_data, file.path(results_path, "fao_sraplus_comp_data.rds"))

fao_sraplus_acc <- fao_sraplus_comp_data %>% 
  group_by(data, fao_area_code) %>% 
  summarise(accuracy = mean(correct, na.rm = TRUE))

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
  # sf::st_simplify() %>%
  mutate(fao_area_code = as.numeric(f_area)) #%>% 
  # st_transform(crs = "+proj=moll")

fao_area_accuracy <- fao_areas %>%
  left_join(fao_sraplus_acc %>% mutate(f_area = as.character(fao_area_code)), by = "f_area") %>%
  filter(!is.na(data))


# world_map <-
#   rnaturalearth::ne_countries(returnclass = "sf", scale = "small")

world_map <-
  rnaturalearth::ne_download(
    scale = 50,
    type = 'land',
    category = 'physical',
    returnclass = "sf"
  ) #%>% 
  # st_transform(crs = "+proj=moll")
  

fao_acc_map_plot <- fao_area_accuracy %>%
  filter(!is.na(data)) %>%
  mutate(data = fct_reorder(data, -accuracy, .fun = mean)) %>% 
  ggplot() +
  geom_sf(aes(fill = accuracy), size = .5) +
  geom_sf(
    data = world_map,
    fill = "darkgrey",
    color = "black",
    size = 0.01
  ) +
  facet_wrap(~ data) +
  scale_fill_viridis(
    limits = c(0,1),
    breaks = c(0,.33, .66,1),
    name = "Accuracy",
    labels = percent,
    guide = guide_colorbar(
      barwidth = ggplot2::unit(9, "lines"),
      axis.linewidth = 1,
      ticks.colour = "white",
      frame.colour = "black"
    )) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white"))


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
  
  
 
  
  
} else {
  
  ram_fit_tests <-   read_rds(path = file.path(results_path,"ram_tests.rds"))
  
  
}


i <- 8

sraplus::plot_prior_posterior(ram_fit_tests$fit[[i]]$fit, ram_fit_tests$fit[[i]]$driors)


compare_to_ram <- function(data, fit){
  
  # data <- ram_fit_data$data[[1]]
  # 
  # fit <- ram_fit_data$fit[[1]]
  
  comparison <- tibble(observed = data$b_v_bmsy) %>% 
    bind_cols(fit$fit$results[fit$fit$results$variable == "b_div_bmsy",])
  
}

ram_test_comparison <- map2_df(ram_fit_tests$data, ram_fit_tests$fit, compare_to_ram,.id = "stock")

# ram_test_comparison %>% 
#   ggplot(aes(observed, mean, color = stock)) + 
#   geom_point(show.legend = FALSE, alpha = 0.5) + 
#   geom_abline(slope = 1, intercept = 0)
# 
# ram_test_comparison %>% 
#   ggplot(aes(observed, mean)) + 
#   geom_hex(show.legend = TRUE, alpha = 0.5, binwidth = c(0.25, 0.25)) + 
#   geom_abline(slope = 1, intercept = 0) + 
#   scale_fill_gradient(low = "lightgrey", high = "tomato")


# run RAM comparison  -----------------------------------------------------------

ram_comp_data <- og_ram_data %>%
  mutate(has_things = !(is.na(catch) | catch == 0)) %>%
  filter(has_things) %>%
  group_by(stockid) %>%
  mutate(delta_year = as.integer(year - lag(year))) %>%
  mutate(delta_year = case_when(year == min(year) ~ as.integer(1),
                                TRUE ~ delta_year)) %>%
  mutate(missing_gaps = any(delta_year > 1)) %>%
  filter(missing_gaps == FALSE) %>%
  group_by(stockid) %>%
  mutate(n = length(catch),
         has_ram_index = any(!is.na(b_v_bmsy))) %>%
  filter(n >= min_years_catch,
         has_ram_index == TRUE) %>%
  ungroup()

ram_comp_data <- ram_comp_data %>%
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
    comm_name <- unique(data$commonname)
    
    isscp_number <- unique(data$isscaap_number)
    
    if (is.na(isscp_number)){
      
      isscp_number <- 33 
    }
    
    fao_code <- unique(data$fao_area_code)
    
    fao_matches <- fao_catch %>% {
      if (any(.$common_name == comm_name & .$fao_area_code == fao_code, na.rm = TRUE)) {
        filter(., common_name == comm_name & fao_area_code == fao_code)
 
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

comp_stocks <- ram_comp_data %>%
  select(scientific_name) %>%
  unique() %>%
  mutate(fishbase_vuln = map_dbl(
    scientific_name,
    ~ rfishbase::species(.x, fields = "Vulnerability")$Vulnerability
  )) %>%
  mutate(
    resilience = dplyr::case_when(
      fishbase_vuln < 33 ~ "Low",
      fishbase_vuln >= 66 ~ "High",
      TRUE ~ "Medium"
    )
  ) %>% 
  select(scientific_name, resilience)

ram_comp_data <- ram_comp_data %>% 
  left_join(comp_stocks, by = "scientific_name") %>% 
  mutate(resilience = ifelse(is.na(resilience), "Medium",resilience))
if (run_ram_comparison == TRUE) {
  
 
  
  # rous_data %>% 
  #   group_by(year, area) %>% 
  #   summarise(effort = mean(effort_cell_reported_nom, na.rm = TRUE)) %>% 
  #   ungroup() %>% 
  #   ggplot(aes(year, effort)) + 
  #   geom_line() + 
  #   facet_wrap(~area, scales = "free_y")
  # 
  # 
  # ram_comp_data %>% 
  #   group_by(year, fao_area_code) %>% 
  #   summarise(effort = mean(effort_index, na.rm = TRUE)) %>% 
  #   ungroup() %>% 
  #   ggplot(aes(year, effort)) + 
  #   geom_line() + 
  #   facet_wrap(~fao_area_code, scales = "free_y")
  # 
  future::plan(multisession, workers = n_cores)
  
  ram_status_fits <- ram_comp_data %>%
    # filter(fao_area_code %in% 67) %>%
    group_by(stockid) %>%
    nest() %>%
    ungroup() %>% 
    # sample_n(40) %>%
    mutate(
      fits = future_map(
        data,
       safely(fit_ram),
        support_data = support_data,
        default_initial_state = NA,
        default_initial_state_cv = NA,
        min_effort_year = 1975,
        engine = "stan",
        cores = 1,
       estimate_shape = FALSE,
       cmsy_cores = 1,
        .progress = TRUE,
        .options = future_options(
          globals = FALSE,
          packages = c("tidyverse", "sraplus", "portedcmsy")
        )
      )
    )
  
  
  check <- ram_status_fits$fits[[1]]$result

  write_rds(ram_status_fits,
            path = file.path(results_path, "ram_status_fits.rds"))
  
} else {
  ram_status_fits <-
    read_rds(path = file.path(results_path, "ram_status_fits.rds"))
  
  
}


# make plots --------------------------------------------------------------


ram_fit_worked <- map_lgl(map(ram_status_fits$fits, "error"), is.null)


ram_status_fits <- ram_status_fits %>% 
  filter(ram_fit_worked) %>% 
  mutate(fits = map(fits, "result"))


compare_to_ram <- function(data, fit){
  
  bdat <- data %>% 
    select(year, b_v_bmsy) %>% 
    rename(ram_b_v_bmsy = b_v_bmsy )
  
  fit <- fit %>% 
    filter(variable == "b_div_bmsy")
  
  comparison <- fit %>% 
    left_join(bdat, by = "year")

}

ram_status_fits <- ram_status_fits %>% 
  mutate(performance = map2(data, fits, compare_to_ram))

i = 24
 # ram_status_fits$performance[[i]] %>%
 #  ggplot() +
 #  geom_point(aes(year, ram_b_v_bmsy, fill = "Observed from RAM"), shape = 21) +
 #   geom_line(aes(year, mean, color = data)) +
 #  facet_wrap(~data) +
 #  scale_x_continuous(name = "B/Bmsy") +
 #  labs(title = ram_status_fits$stockid[[i]])

assess_ram_fits <- ram_status_fits %>% 
  select(stockid, performance) %>% 
  unnest(cols = performance) %>% 
  filter(!is.na(ram_b_v_bmsy)) %>% 
  group_by(stockid) %>% 
  filter(year == max(year)) %>% 
  ungroup() %>% 
  mutate(ram_b_v_bmsy = pmin(ram_b_v_bmsy,5),
         mean = pmin(mean, 5)) %>% 
  mutate(resid = ram_b_v_bmsy - mean,
         ae = abs(resid)) %>% 
  filter(!data %in%  c("heuristic", "catch_only")) 

null_model <- assess_ram_fits %>% 
  filter(data == "cmsy") %>% 
  mutate(mean = sample(c(.4,1,1.6), n(), replace = TRUE)) %>% 
  mutate(data = "guess") %>% 
  mutate(resid = ram_b_v_bmsy - mean,
         ae = abs(resid))

assess_ram_fits <- assess_ram_fits %>% 
  bind_rows(null_model)

write_rds(assess_ram_fits, path = file.path(results_path,"assess_ram_fits.rds"))

write_rds(ram_comp_data, path = file.path(results_path,"ram_comp_data.rds"))

# assess_ram_fits %>% 
#   ggplot(aes(pmin(5,ram_b_v_bmsy),pmin(5,mean))) + 
#   geom_hex(binwidth = c(0.33, 0.33), color = "white") + 
#   geom_smooth(method = "lm", aes(color = "fit")) +
#   geom_abline(slope = 1, intercept = 0) +
#   facet_wrap(~data, scales = "free") 
#   

ram_v_sraplus_plot <- assess_ram_fits %>% 
  ggplot(aes(pmin(ram_b_v_bmsy,5), pmin(mean,5))) +
  geom_point(alpha = 0.5) +
  # geom_hex(binwidth = c(0.25, 0.25)) +
  geom_smooth(method = "lm", aes(color = "Observed v. Predicted Fit")) +
  geom_abline(size = 2,
              linetype = 2,
              aes(
                slope = 1,
                intercept = 0,
                color = "1:1 Line"
              )) +
  facet_wrap( ~ data, scales = "free") +
  scale_fill_viridis() +
  scale_color_discrete(name = '') +
  scale_x_continuous(name = "Observed B/Bmsy", limits = c(0, NA)) +
  scale_y_continuous(name = "Predicted B/Bmsy", limits =  c(0, NA))

ram_v_sraplus_resid_plot <- assess_ram_fits %>% 
  ggplot(aes(resid, fill = data)) + 
  geom_vline(aes(xintercept = 0) ,linetype = 2) +
  geom_histogram(bins = 50) +
  facet_wrap(~data)

ram_v_sraplus_ae_plot <- assess_ram_fits %>% 
  ggplot(aes(ae, fill = data)) + 
  geom_vline(aes(xintercept = 0) ,linetype = 2) +
  geom_histogram(bins = 50) +
  facet_wrap(~data)

ram_v_sraplus_area <- assess_ram_fits %>%
  left_join(ram_comp_data %>% select(stockid, primary_FAOarea) %>% unique(),
            by = "stockid") %>%
  rename(primary_fao_area = primary_FAOarea) %>% 
  group_by(primary_fao_area, data) %>%
  mutate(sraplus_bin = cut(mean, breaks = breaks, labels = labels)) %>%
  mutate(ram_bin = cut(ram_b_v_bmsy, breaks = breaks, labels = labels)) %>%
  summarise(median_sraplus = median(mean, na.rm = TRUE),
            median_ram = median(ram_b_v_bmsy, na.rm = TRUE),
            accuracy = mean(sraplus_bin == ram_bin, na.rm = TRUE),
            mape = median(abs((mean -  ram_b_v_bmsy) / ram_b_v_bmsy), na.rm = TRUE),
            mpe = median(((mean -  ram_b_v_bmsy) / ram_b_v_bmsy), na.rm = TRUE),
            rmse = mean(((mean -  ram_b_v_bmsy)^2)), na.rm = TRUE) %>% 
  group_by(primary_fao_area)

write_rds(ram_v_sraplus_area, path = file.path(results_path,"ram_v_sraplus_area.rds"))

  

fao_area_ram_status <- fao_areas %>%
  left_join(ram_v_sraplus_area,
            by = c("f_area" = "primary_fao_area")) %>%
  mutate(pe = (median_sraplus - median_ram) / median_ram) %>% 
  mutate(data = fct_relevel(data, "ram-data")) %>% 
  filter(data != "u_umsy")


ram_status <- fao_area_ram_status %>% 
  filter(!is.na(data)) %>% 
  group_by(f_area) %>% 
  slice(1) %>% 
  select(-data, -median_sraplus) %>% 
  mutate(data = "RLSADB Assessment", median_sraplus = median_ram) 
  
  

fao_area_ram_status %>% 
  group_by(data) %>% 
  summarise(mape = median(abs(pe))) %>% 
  arrange(mape)


ram_labeller <- c(
  "ram-data" = "RLSADB Index",
  "sar" = "SAR",
  "fmi" = "FMI",
  "cpue" = "Effective CPUE",
  "cpue-plus" = "Effective CPUE+",
  "nominal-cpue" = "Nominal CPUE",
  "nominal-cpue-plus" = "Nominal CPUE+",
  "cmsy" = "CMSY",
  "guess" = "Guess",
  "u_umsy" = "RLSADB U/Umsy"
)


collabs <-
  c(paste0(seq(-100, 75, by = 25),"%"),
    expression("" >= "100%"))

ram_mpe_map_plot <- fao_area_ram_status %>%
  filter(!is.na(data)) %>% 
  mutate(data = fct_reorder(data, abs(mpe), .fun = median)) %>% 
  ggplot() +
  geom_sf(aes(fill = pmin(1,mpe)), size = .01) +
  geom_sf(
    data = world_map,
    fill = "darkgrey",
    color = "black",
    size = 0.01
  ) +
  facet_wrap(~ data, labeller = labeller(data = ram_labeller)) +
  scale_fill_gradient2(
    low = "darkblue",
    high = "tomato",
    mid = "white",
    name = "% Bias (MPE)",
    labels = collabs,
    midpoint = 0,
    breaks = seq(-1, 1, by = .25),
    guide = guide_colorbar(
      barwidth = ggplot2::unit(15, "lines"),
      axis.linewidth = 1,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) + 
  theme(legend.position = "top",
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white"),
        legend.text = element_text(size = 8))

collabs <-
  c(paste0(seq(0, 75, by = 25),"%"),
    expression("" >= "100%"))

ram_mape_map_plot <- fao_area_ram_status %>%
  filter(!is.na(data)) %>% 
  mutate(data = fct_reorder(data, mape, .fun = median)) %>% 
  ggplot() +
  geom_sf(aes(fill = pmin(1,mape)), size = .01) +
  geom_sf(
    data = world_map,
    fill = "darkgrey",
    color = "black",
    size = 0.01
  ) +
  facet_wrap(~ data, labeller = labeller(data = ram_labeller)) +
  scale_fill_gradient2(
    low = "blue",
    mid = "tomato",
    high = "yellow",
    name = "% Error (MAPE)",
    labels = collabs,
    breaks = seq(0,1, by = .25),
    midpoint = 0.5,
    limits = c(0,NA),
    # trans = "log10",
    guide = guide_colorbar(
      barwidth = ggplot2::unit(9, "lines"),
      axis.linewidth = 1,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) + 
  theme(legend.position = "top",
        legend.direction = "horizontal",
        panel.background = element_rect(fill = "white"))

  


ram_acc_map_plot <- fao_area_ram_status %>%
  filter(!is.na(data)) %>%
  mutate(data = fct_reorder(data, -accuracy, .fun = mean)) %>% 
  ggplot() +
  geom_sf(aes(fill = accuracy), size = .01) +
  geom_sf(
    data = world_map,
    fill = "darkgrey",
    color = "black",
    size = 0.01
  ) +
  facet_wrap(~ data, labeller = labeller(data = ram_labeller)) +
  scale_fill_viridis(
    limits = c(0,1),
    breaks = c(0,.33, .66,1),
    name = "Accuracy",
    labels = percent,
    guide = guide_colorbar(
      barwidth = ggplot2::unit(9, "lines"),
      axis.linewidth = 1,
      ticks.colour = "white",
      frame.colour = "black"
    )) +
    theme(legend.position = "top",
          legend.direction = "horizontal",
          panel.background = element_rect(fill = "white"))
    
ram_status <- fao_area_ram_status %>% 
  filter(!is.na(data)) %>% 
  group_by(f_area) %>% 
  slice(1) %>% 
  select(-data, -median_sraplus) %>% 
  mutate(data = "RLSADB Assessment", median_sraplus = median_ram) 

ram_status <- ram_status[,colnames(fao_area_ram_status)]

combo <- fao_area_ram_status %>% 
  bind_rows(ram_status)


ram_b_map_plot <- combo %>%
  filter(!is.na(data)) %>%
  mutate(data = fct_relevel(data, "RLSADB Assessment")) %>%
  ggplot() +
  geom_sf(aes(fill = median_sraplus), size = .01) +
  geom_sf(
    data = world_map,
    fill = "grey",
    color = "black",
    size = 0.01
  ) +
  facet_wrap(~ data, labeller = labeller(data = ram_labeller)) +
  scale_fill_viridis(
    limits = c(0, NA),
    name = "Median B/Bmsy",
    guide = guide_colorbar(
      barwidth = ggplot2::unit(12, "lines"),
      axis.linewidth = 1,
      ticks.colour = "black",
      frame.colour = "black"
    )
  ) + 
  theme(legend.position = "top",
        panel.background = element_rect(fill = "white"))



# ----save-plots----------------------------------------------------------

plots <- ls()[str_detect(ls(), "_plot")]

save(file = here("results", results_name, "paper_plots.Rdata"),
     list = plots)



# knit paper --------------------------------------------------------------
if (knit_paper == TRUE){
  
  rmarkdown::render(
    here::here("documents", "ovando-etal-assessing-global-fisheries.Rmd"),
    params = list(results_name = results_name,
                  min_years_catch = 25)
  )
  
  # rmarkdown::render(
  #   here::here("documents", "ovando-etal-assessing-global-fisheries-si.Rmd"),
  #   params = list(results_name = results_name,
  #                 min_years_catch = 25)
  # )
  # 
}
