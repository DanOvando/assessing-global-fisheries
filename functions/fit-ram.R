fit_ram <-
  function(data,
           support_data,
           model = "sraplus_tmb",
           engine = "tmb",
           include_fit = FALSE,
           min_effort_year = 1980,
           u_cv = 0.2,
           q_slope_prior = 0.025,
           estimate_proc_error = TRUE,
           estimate_qslope = FALSE,
           estimate_shape = FALSE,
           default_initial_state = NA,
           default_initial_state_cv = 0.2,
           chains = 1,
           cores = 1,
           cmsy_cores = 1,
           reg_cv = NA,
           draws = 2e6,
           refresh = 0,
           thin_draws = TRUE) {
    #
    #     data <- test$data[[1]]
    # browser()
    
    sfs <- safely(fit_sraplus)
    
    ogengine <- engine
    
    scientific_name <- unique(data$scientific_name)
    
    catches <- data$catch
    
    years <- data$year
    
    sfc <- purrr::safely(portedcmsy::funct_cmsy)
    # browser()
    ported_cmsy <- sfc(
      catches = catches[years >= 1950],
      catch_years = years[years >= 1950],
      stock = scientific_name,
      common_name = unique(data$common_name),
      scientific_name = unique(data$scientific_name),
      res = unique(data$resilience),
      start.yr = years[years >= 1950][1],
      end.yr = dplyr::last(years[years >= 1950]),
      cores = cmsy_cores
    )
    write_rds(ported_cmsy, "wtf.rds")
    
    if (is.null(ported_cmsy$error)){
      
      cmsy_results <- ported_cmsy$result %>% 
        mutate(c_div_msy  = ct / msy) %>% 
        select(year, bt, b_bmsy, f_fmsy, c_div_msy) %>% 
        rename(depletion = bt, b_div_bmsy = b_bmsy, u_div_umsy = f_fmsy) %>% 
        pivot_longer(-year, names_to = "variable", values_to = "mean") %>% 
        mutate(sd = NA, lower = NA, upper = NA, data = "cmsy")
      
      # tests <- cmsy_results %>% 
      #   bind_rows(basic_fit_results)
      # 
      # tests %>% 
      #   filter(variable == "c_div_msy") %>% 
      #   select(year, mean, data) %>% 
      #   ggplot(aes(year, mean, color = data)) + 
      #   geom_line()
      
      cmsy_worked <- TRUE 
    } else {
      cmsy_worked <- FALSE
    }
    
    u_years <- seq_along(years)[!is.na(data$mean_u_umsy)]
    
    u_umsy <- data$mean_u_umsy[u_years]
    basic_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        use_heuristics = TRUE,
        b_ref_type = "k",
        isscaap_group = unique(data$isscaap_group)
      )
    
    basic_fit <-
      sfs(
        driors = basic_driors,
        include_fit = include_fit,
        estimate_shape = estimate_shape,
        model = model,
        engine = "sir",
        draws = draws,
        thin_draws = thin_draws
      )
    
    if (is.null(basic_fit$error)){
      
      basic_fit_results <- basic_fit$result$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
        mutate(year = rep(years, each = 4)) %>%
        mutate(data = "heuristic")
 
      basic_worked <- TRUE 
    } else {
      
      basic_fit_results <- tibble(year = NA,variable = NA, mean = NA, sd = NA, lower = NA, upper = NA, data = NA)
      
      basic_worked <- FALSE
    }
    
    
 
    
    
    # updated cmsy fit
    
  

    # catch only fit
    
    com_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        use_heuristics = FALSE,
        initial_state = NA,
        initial_state_cv = NA,
        b_ref_type = "k",
        isscaap_group = unique(data$isscaap_group),
        use_catch_prior = TRUE
      )
    com_fit <-
      fit_sraplus(
        driors = com_driors,
        include_fit = include_fit,
        model = model,
        engine = "sir",
        draws = draws,
        estimate_shape = estimate_shape,
        
        thin_draws = thin_draws
      )
    
    com_fit_results <- com_fit$results %>%
      filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
      mutate(year = rep(years, each = 4)) %>%
      mutate(data = "catch_only")
    
    
    # sar fit
    
    if (!is.na(last(data$sar))) {
      used_sar <- TRUE
      
      sar_driors <-
        sraplus::format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(catches),
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state,
          sar = last(data$sar),
          sar_cv = reg_cv,
          isscaap_group = unique(data$isscaap_group)
        )
      
      sar_fit <-
        fit_sraplus(
          driors = sar_driors,
          include_fit = include_fit,
          model = model,
          engine = "sir",
          draws = draws,
          estimate_shape = estimate_shape,
          
          thin_draws = thin_draws
          
        )
      
      sar_fit_results <- sar_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
        mutate(year = rep(years, each = 4)) %>%
        mutate(data = "sar")
      
      rm(sar_fit)
      
    } else {
      used_sar <- FALSE
    }
    
    
    # fmi fit
    
    fmi_dat <- support_data$mean_regional_isscaap_fmi %>%
      filter(fao_area_code == unique(as.numeric(data$fao_area_code)),
             isscaap_group == unique(data$isscaap_group))
    
    
    fmi <-
      c(
        research = fmi_dat$research,
        management = fmi_dat$management,
        socioeconomics = fmi_dat$socioeconomics,
        enforcement = fmi_dat$enforcement
      )
    
    if (nrow(fmi_dat) == 0) {
      fmi_dat <- support_data$mean_regional_isscaap_fmi %>%
        filter(fao_area_code == unique(as.numeric(data$fao_area_code))) %>%
        pivot_longer(
          cols = c("research", "management", "socioeconomics", "enforcement"),
          names_to = "metric",
          values_to = "value"
        ) %>%
        group_by(fao_area_code, metric) %>%
        summarise(mean_value = mean(value)) %>%
        pivot_wider(names_from = metric, values_from = mean_value)
      
      fmi <-
        c(
          research = fmi_dat$research,
          management = fmi_dat$management,
          socioeconomics = fmi_dat$socioeconomics,
          enforcement = fmi_dat$enforcement
        )
      
    }
    
    if (nrow(fmi_dat) == 0) {
      fmi <- NA
      
    }
    
    
    if (!all(is.na(fmi))) {
      # term_b_guess <- pmax(0.05,(2.5 - exp(c(last(data$mean_fmi_ppf)))))
      
      fmi_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(years),
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state_cv,
          fmi = fmi,
          f_ref_type = "fmsy",
          fmi_cv = reg_cv,
          use_b_reg = FALSE,
          isscaap_group = unique(data$isscaap_group)
        )
      fmi_fit <-
        fit_sraplus(
          driors = fmi_driors,
          include_fit = include_fit,
          model = model,
          estimate_shape = estimate_shape,
          
          engine = "sir",
          draws = draws,
          thin_draws = thin_draws
          
        )
      
      fmi_fit_results <- fmi_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
        mutate(year = min(years) + year - 1) %>%
        mutate(data = "fmi")
      
      rm(fmi_fit)
      
      used_fmi <- TRUE
    } else {
      used_fmi <- FALSE
    }
    
    # U/Umsy fit
    
    if (any(is.na(u_umsy) == FALSE)) {
      used_u <- TRUE
      
      u_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(catches),
          u = u_umsy,
          u_cv = u_cv,
          u_years = u_years,
          f_ref_type = "fmsy",
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state_cv,
          isscaap_group = unique(data$isscaap_group)
        )
      
      u_fit <-
        sfs(
          driors = u_driors,
          include_fit = include_fit,
          model = model,
          estimate_shape = estimate_shape,
          
          engine = "sir",
          draws = draws,
          thin_draws = thin_draws
        )
      
      
      if (is.null(u_fit$error)) {
        u_fit <- u_fit$result
        
        u_fit_results <- u_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
          mutate(year = min(years) + year - 1) %>%
          mutate(data = "u_umsy")
        
      } else {
        used_u <- FALSE
        
      }
      
    }
    else {
      used_u <- FALSE
    }
    
    # cpue fit
    
    if (any(is.na(data$effort_index) == FALSE)) {
      used_cpue <- TRUE
      
      effort_years <-
        seq_along(years)[!is.na(data$effort_index) &
                           years >= min_effort_year]
      cpue_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(years),
          effort = data$effort_index[effort_years],
          effort_years = effort_years,
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state_cv,
          isscaap_group = unique(data$isscaap_group),
          q_slope_prior = q_slope_prior
        )
      sfs <- safely(fit_sraplus)
      
      cpue_fit <-
        sfs(
          driors = cpue_driors,
          include_fit = include_fit,
          model = model,
          engine = engine,
          estimate_proc_error = estimate_proc_error,
          estimate_qslope = estimate_qslope,
          workers = cores,
          estimate_shape = estimate_shape,
          
          refresh = refresh,
          thin_draws = thin_draws
        )
      
      nom_cpue_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(years),
          effort = data$effort_index[effort_years],
          effort_years = effort_years,
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state_cv,
          isscaap_group = unique(data$isscaap_group),
          q_slope_prior = 0
        )
      
      nom_cpue_fit <-
        sfs(
          driors = nom_cpue_driors,
          include_fit = include_fit,
          model = model,
          engine = engine,
          estimate_proc_error = estimate_proc_error,
          estimate_qslope = estimate_qslope,
          estimate_shape = estimate_shape,
          
          workers = cores,
          refresh = refresh
        )
      # browser()
      # plot_sraplus(nom_cpue_fit$result)
      #
      #
      # plot_prior_posterior(nom_cpue_fit$result, nom_cpue_driors)
      #
      # plot_prior_posterior(cpue_fit$result, cpue_driors)
      #
      # plot_sraplus(cpue_fit$result, nom_cpue_fit$result)
      #
      # plot_driors(cpue_driors)
      #
      # plot_driors(nom_cpue_driors)
      #
      # nom_cpue_fit$result$fit
      
      
      if (is.null(cpue_fit$error)) {
        cpue_fit <- cpue_fit$result
        
        nom_cpue_fit <- nom_cpue_fit$result
        
      } else {
        # cpue_fit <-
        #   fit_sraplus(driors = cpue_driors,
        #               include_fit = TRUE,
        #               model = model,
        #               engine = "stan",
        #               estimate_proc_error = estimate_proc_error,
        #               estimate_qslope = estimate_qslope,
        #               chains = chains,
        #               cores = cores)
        #
        # engine = "stan"
      }
      
      if (engine == "tmb") {
        cpue_fit_results <- cpue_fit$results %>%
          filter(
            variable %in% c(
              "log_b_div_bmsy",
              "log_u_div_umsy",
              "log_c_div_msy",
              "log_depletion"
            )
          ) %>%
          mutate(year = rep(years, 4)) %>%
          modify_at(c("mean", "lower", "upper"), exp) %>%
          mutate(data = "cpue") %>%
          mutate(variable = str_replace_all(variable, "log_", ""))
        
        rm(cpue_fit)
        
        nom_cpue_fit_results <- nom_cpue_fit$results %>%
          filter(
            variable %in% c(
              "log_b_div_bmsy",
              "log_u_div_umsy",
              "log_c_div_msy",
              "log_depletion"
            )
          ) %>%
          mutate(year = rep(years, 4)) %>%
          modify_at(c("mean", "lower", "upper"), exp) %>%
          mutate(data = "nominal-cpue") %>%
          mutate(variable = str_replace_all(variable, "log_", ""))
        
        rm(nom_cpue_fit)
      } else {
        cpue_fit_results <- cpue_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
          mutate(year = year - 1 + min(years)) %>%
          mutate(data = "cpue")
        
        rm(cpue_fit)
        
        nom_cpue_fit_results <- nom_cpue_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
          mutate(year = year - 1 + min(years)) %>%
          mutate(data = "nominal-cpue")
        
        rm(nom_cpue_fit)
        
      }
      
    } else {
      used_cpue <- FALSE
    }
    
    engine <- ogengine
    
    # cpue plus
    
    if (any(is.na(data$effort_index) == FALSE) &
        !all(is.na(fmi)) | !all(is.na(data$sar))) {
      used_cpue_plus <- TRUE
      
      effort_years <-
        seq_along(years)[!is.na(data$effort_index) &
                           years >= min_effort_year]
      
      
      # dyn.load(TMB::dynlib(here::here("src", "sraplus_v2")))
      cpue_plus_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(years),
          effort = data$effort_index[effort_years],
          effort_years = effort_years,
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state_cv,
          fmi = fmi,
          f_ref_type = "fmsy" ,
          fmi_cv = reg_cv,
          sar = unique(data$sar),
          isscaap_group = unique(data$isscaap_group),
          q_slope_prior = q_slope_prior
        )
      
      sfs <- safely(fit_sraplus)
      
      cpue_plus_fit <-
        sfs(
          driors = cpue_plus_driors,
          include_fit = include_fit,
          model = model,
          engine = engine,
          estimate_proc_error = estimate_proc_error,
          estimate_qslope = estimate_qslope,
          workers = cores,
          refresh = refresh,
          thin_draws = thin_draws,
          estimate_shape = estimate_shape
        )
      
      nom_cpue_plus_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(years),
          effort = data$effort_index[effort_years],
          effort_years = effort_years,
          initial_state = default_initial_state,
          initial_state_cv = default_initial_state_cv,
          fmi = fmi,
          f_ref_type = "fmsy" ,
          fmi_cv = reg_cv,
          sar = unique(data$sar),
          isscaap_group = unique(data$isscaap_group),
          q_slope_prior = 0
        )
      
      nom_cpue_plus_fit <-
        sfs(
          driors = nom_cpue_plus_driors,
          include_fit = include_fit,
          model = model,
          engine = engine,
          estimate_proc_error = estimate_proc_error,
          estimate_qslope = estimate_qslope,
          estimate_shape = estimate_shape,
          workers = cores,
          refresh = refresh,
          thin_draws = thin_draws
        )
      
      
      if (is.null(cpue_plus_fit$error)) {
        cpue_plus_fit <- cpue_plus_fit$result
        
        nom_cpue_plus_fit <- nom_cpue_plus_fit$result
        
      } else {
        # cpue_plus_fit <-
        #   fit_sraplus(driors = cpue_plus_driors,
        #               include_fit = TRUE,
        #               model = model,
        #               fit_catches = fit_catches,
        #               engine = "stan",
        #               estimate_proc_error = estimate_proc_error,
        #               estimate_qslope = estimate_qslope,
        #               chains = chains,
        #               cores = cores)
        #
        # engine = "stan"
      }
      
      if (engine == "tmb") {
        cpue_plus_fit_results <- cpue_plus_fit$results %>%
          filter(
            variable %in% c(
              "log_b_div_bmsy",
              "log_u_div_umsy",
              "log_c_div_msy",
              "log_depletion"
            )
          ) %>%
          mutate(year = rep(years, 4)) %>%
          modify_at(c("mean", "lower", "upper"), exp) %>%
          mutate(data = "cpue-plus") %>%
          mutate(variable = str_replace_all(variable, "log_", ""))
        
        rm(cpue_plus_fit)
        
        nom_cpue_plus_fit_results <- nom_cpue_plus_fit$results %>%
          filter(
            variable %in% c(
              "log_b_div_bmsy",
              "log_u_div_umsy",
              "log_c_div_msy",
              "log_depletion"
            )
          ) %>%
          mutate(year = rep(years, 4)) %>%
          modify_at(c("mean", "lower", "upper"), exp) %>%
          mutate(data = "nominal-cpue-plus") %>%
          mutate(variable = str_replace_all(variable, "log_", ""))
        
        rm(nom_cpue_plus_fit)
        
      } else {
        cpue_plus_fit_results <- cpue_plus_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
          mutate(year = year - 1 + min(years)) %>%
          mutate(data = "cpue-plus")
        
        rm(cpue_plus_fit)
        
        nom_cpue_plus_fit_results <- nom_cpue_plus_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
          mutate(year = year - 1 + min(years)) %>%
          mutate(data = "nominal-cpue-plus")
        
        rm(nom_cpue_plus_fit)
        
      }
      
    } else {
      used_cpue_plus <- FALSE
    }
    
    # fit with real RAM data
    
    index_years <- seq_along(years)[!is.na(data$b_v_bmsy)]
    
    used_ram_data <- TRUE
    
    ram_data_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        index = data$b_v_bmsy[!is.na(data$b_v_bmsy)] * 100,
        index_years = index_years,
        initial_state = default_initial_state,
        initial_state_cv = default_initial_state_cv,
        isscaap_group = unique(data$isscaap_group)
      )
    
    # plot_driors(ram_data_driors)
    sfs <- safely(fit_sraplus)
    
    ram_data_fit <-
      sfs(
        driors = ram_data_driors,
        include_fit = include_fit,
        model = model,
        engine = engine,
        estimate_proc_error = estimate_proc_error,
        estimate_qslope = estimate_qslope,
        estimate_shape = estimate_shape,
        workers = cores,
        refresh = refresh,
        thin_draws = thin_draws
      )
    # browser()
    # plot_prior_posterior(ram_data_fit$result, ram_data_driors)
    if (is.null(ram_data_fit$error)) {
      ram_data_fit <- ram_data_fit$result

    } else {
      used_ram_data <- FALSE
    }
    # 
    if (engine == "tmb") {
      ram_data_fit_results <- ram_data_fit$results %>%
        filter(
          variable %in% c(
            "log_b_div_bmsy",
            "log_u_div_umsy",
            "log_c_div_msy",
            "log_depletion"
          )
        ) %>%
        mutate(year = rep(years, 4)) %>%
        modify_at(c("mean", "lower", "upper"), exp) %>%
        mutate(data = "ram-data") %>%
        mutate(variable = str_replace_all(variable, "log_", ""))

      rm(ram_data_fit)


    } else {
      ram_data_fit_results <- ram_data_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
        mutate(year = year - 1 + min(years)) %>%
        mutate(data = "ram-data")

      rm(ram_data_fit)

    } 

    
    
    # store results
    
    results <- 
      basic_fit_results  %>%
      bind_rows(com_fit_results) %>% {
        if (used_sar == TRUE) {
          bind_rows(., sar_fit_results)
        }
        else{
          .
        }
      } %>% {
        if (used_cpue == TRUE) {
          bind_rows(., cpue_fit_results) %>%
            bind_rows(nom_cpue_fit_results)
        }
        else{
          .
        }
      } %>%
      {
        if (used_u == TRUE) {
          bind_rows(., u_fit_results)
        }
        else{
          .
        }
      } %>% {
        if (used_fmi == TRUE) {
          bind_rows(., fmi_fit_results)
        }
        else{
          .
        }
      } %>% {
        if (used_cpue_plus == TRUE) {
          bind_rows(., cpue_plus_fit_results) %>%
            bind_rows(nom_cpue_plus_fit_results)
        }
        else{
          .
        }
      } %>% {
        if (used_ram_data == TRUE) {
          bind_rows(., ram_data_fit_results)
        } else{
          .
        }
      } %>% {
        if (cmsy_worked == TRUE) {
          bind_rows(., cmsy_results)
        } else {
          .
        }
      } %>% 
      filter(!is.na(data))
    
    
    # write(wtf,file = "wtf.txt", append = TRUE)
    
    # results %>%
    #   filter(data == "sar") %>%
    #   ggplot(aes(year, mean, color = data)) +
    #   geom_line() +
    #   facet_wrap(~variable, scales = "free_y")
    
    rm(list = ls()[!str_detect(ls(), "results")])
    gc()
    return(results)
    
  } 