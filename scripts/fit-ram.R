  library(tidyverse)
  library(here)
  library(sraplus)
  library(furrr)
  
  # options -----------------------------------------------------------------
  
  min_years_catch <- 20
  
  crazy_b <- 4 # threshold for suspect B/Bmsy value
  
  crazy_u <- 5 # threshold for suspect U/Umsy value
  
  
  # load data ---------------------------------------------------------------
  
  
  if(!dir.exists("data-raw")){
    dir.create("data-raw")
  }
  
  if (!file.exists(here("data-raw","ram.zip"))) {
    
    download.file("https://www.dropbox.com/s/jpgz0a5s5of3qev/RAM%20v4.491%20Files%20(1-14-20).zip?dl=1", destfile = here::here("data-raw","ram.zip"))
    
    unzip(here::here("data-raw","ram.zip"), exdir = "data-raw")
    
  }
  
  
  
  # process RAM data --------------------------------------------------------
  
  ram_dirs <- list.files("data-raw")
  
  ram_dirs <- ram_dirs[str_detect(ram_dirs,"RAM v\\d")]
  
  ram_files <- list.files(here("data-raw",ram_dirs), recursive = TRUE)
  
  ram_files <- ram_files[str_detect(ram_files,".RData")]
  
  ram_files <- ram_files[str_detect(ram_files,"Assessment Data")]
  
  load(here("data-raw",ram_dirs,ram_files[1]))
  
  stock <- stock %>%
    left_join(area, by = "areaid")
  # catches
  ram_catches <- tcbest.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, catch,-year)
  
  # B/Bmsy
  ram_b_v_bmsy <- divbpref.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, b_v_bmsy,-year)
  
  # U/Umsy
  ram_u_v_umsy <- divupref.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, u_v_umsy,-year)
  
  # Effort
  ram_effort <- effort.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, effort,-year)
  
  # biomass
  
  
  ram_total_biomass <- tbbest.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, total_biomass,-year)
  
  # ssb
  
  ram_ss_biomass <- ssb.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, ss_biomass,-year)
  
  
  ram_exp_rate <- ram_catches %>%
    left_join(ram_total_biomass, by = c("stockid", "year")) %>%
    mutate(exploitation_rate = catch / total_biomass) %>%
    select(-catch, -total_biomass)
  
  
  # put it together
  
  ram_data <- ram_catches %>%
    left_join(bioparams_values_views, by = "stockid") %>%
    left_join(ram_b_v_bmsy, by = c("stockid", "year")) %>%
    left_join(ram_u_v_umsy, by = c("stockid", "year")) %>%
    left_join(ram_exp_rate, by = c("stockid", "year")) %>%
    left_join(ram_effort, by = c("stockid", "year")) %>%
    left_join(ram_total_biomass, by = c("stockid", "year")) %>%
    left_join(ram_ss_biomass, by = c("stockid", "year")) %>%
    left_join(stock, by = "stockid") %>%
    select(stockid, scientificname, commonname, everything())
  
  
  # create new variables
  
  ram_data <- ram_data %>%
    mutate(tb_v_tb0 = total_biomass / TB0,
           ssb_v_ssb0 = ss_biomass / SSB0)
  
  # filter data
  
  # for now only include stocks taht have total biomass
  ram_data <- ram_data %>%
    filter(is.na(catch) == FALSE & is.na(total_biomass) == FALSE) %>%
    group_by(stockid) %>%
    mutate(delta_year = year - lag(year)) %>%
    mutate(delta_year = case_when(year == min(year) ~ as.integer(1),
                                  TRUE ~ delta_year)) %>%
    mutate(missing_gaps = any(delta_year > 1)) %>%
    filter(missing_gaps == FALSE) %>%
    mutate(n_years = n_distinct(year)) %>%
    filter(n_years >= min_years_catch) %>%
    filter(all(b_v_bmsy < crazy_b, na.rm = TRUE),
           all(u_v_umsy < crazy_u, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(stockid) %>%
    mutate(
      has_tb0 = !all(is.na(TB0)),
      has_tb = !all(is.na(total_biomass)),
      first_catch_year = year[which(catch > 0)[1]]
    ) %>%
    filter(year >= first_catch_year) %>%
    mutate(
      pchange_effort = lead(u_v_umsy) / (u_v_umsy + 1e-6),
      cs_effort = (u_v_umsy - mean(u_v_umsy)) / sd(u_v_umsy),
      b_rel = dplyr::case_when(
        has_tb0 ~ total_biomass / max(TB0),
        has_tb ~ total_biomass / max(total_biomass),
        TRUE ~ b_v_bmsy / 2
      )
    ) %>%
    ungroup()
  
  
  # fit ram stocks ----------------------------------------------------------
  
  
  
  # for now, only including things that hve total biomass estimates
  ram_fit_data <- ram_data %>% 
    select(year, stockid, scientificname, commonname, year, catch,b_v_bmsy, u_v_umsy, total_biomass) %>% 
    group_by(stockid,scientificname) %>% 
    nest() %>% 
    ungroup() 
  
  fit_model <-
    function(sciname,
             data,
             survey_q = 1,
             estimate_proc_error = FALSE, #set to TRUE to estimate process error
             estimate_shape = FALSE) { # set to TRUE to estimate shape parameter of PT model
      driors <- sraplus::format_driors(
        taxa = sciname,
        catch = data$catch,
        years = data$year,
        index = data$total_biomass * survey_q,
        index_years = data$year,
        q_prior = survey_q,
        q_prior_cv = 0.01
      )
      
      # sraplus::plot_driors(driors)
      
      fit <- sraplus::fit_sraplus(driors = driors,
                                  engine = "stan",
                                  estimate_proc_error = estimate_proc_error,
                                  estimate_shape = estimate_shape,
                                  estimate_q = FALSE,
                                  workers = 4,
                                  refresh = 0)
      
      # b <- fit$results %>% 
      #   filter(variable == "b_t")
      
      # plot(b$mean, data$total_biomass)
      # abline(a = 0, b = 1)
      # 
      out <- list(fit = fit,
                  driors = driors)
    }
  
  future::plan(multiprocess, workers = 5)
  
  a <- Sys.time()
  ram_fit_data <- ram_fit_data %>%
    mutate(fit = future_map2(
      scientificname,
      data,
      safely(fit_model),
      estimate_proc_error = TRUE,
      estimate_shape = TRUE,
      .progress= TRUE
    ))
  
  fit_worked <- map_lgl(map(ram_fit_data$fit, "error"), is.null)
  Sys.time() - a

  ram_fit_data <- ram_fit_data %>% 
  filter(fit_worked)

ram_fit_data$fit <- map(ram_fit_data$fit, "result")

write_rds(ram_fit_data, path = "ram_fit_data.rds")

  
i <- 5
sraplus::plot_prior_posterior(ram_fit_data$fit[[i]]$fit, ram_fit_data$fit[[i]]$driors)

ram_fits <- ram_fit_data %>% 
  mutate(fit_summary = map(fit, c("fit","results"))) %>% 
  select(stockid, scientificname, fit_summary) %>% 
  unnest(cols = fit_summary)


write_rds(ram_fits, path = "ram_fits.rds")

write_rds(ram_data, path = "ram_data.rds")

b <- ram_fits %>% 
  select(stockid, year, variable, mean) %>% 
  filter(variable == "b_t")


error <- ram_data %>%
  left_join(b, by = c("stockid", "year")) %>% 
left_join(assessment, by = "stockid")

error %>% 
  filter(year > 1950, !is.na(mean)) %>% 
  ggplot(aes(total_biomass, mean)) + 
  geom_point(aes(color = stockid), show.legend = FALSE) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(name = "RAM total biomass") + 
  scale_y_continuous(name = "sraplus estimated total biomass")

error %>% 
  group_by(stockid) %>% 
  mutate(total_biomass = scale(total_biomass), mean = scale(mean)) %>% 
  ggplot(aes(total_biomass, mean)) + 
  geom_point(aes(color = stockid), show.legend = FALSE) + 
  geom_abline(slope = 1, intercept = 0) +
  scale_x_continuous(name = "Centered and scaled RAM total biomass") + 
  scale_y_continuous(name = "Centered and scaled sraplus estimated total biomass")

error %>% 
  yardstick::mape(truth = total_biomass, estimate = mean)
  
  
  
head(ram_fits)
