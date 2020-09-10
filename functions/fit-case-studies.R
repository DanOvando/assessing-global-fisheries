fit_case_studies <-
  function(
           data,
           model = "sraplus_tmb",
           engine = "stan",
           estimate_shape = TRUE,
           estimate_proc_error = TRUE,
           estimate_qslope = FALSE,
           use_baranov = TRUE,
           first_effort_year = 1950,
           chains = 1,
           cores = 1,
           cmsy_cores = 4,
           q_slope_prior = 0.025) {
 
    catches <- data$catch
    
    years <- data$year
    
    scientific_name <- data$scientificname[[1]]
    
    # fit cmsy
    
    sfc <- purrr::safely(portedcmsy::funct_cmsy)
    
    ported_cmsy <- sfc(
      catches = catches[years >= 1950],
      catch_years = years[years >= 1950],
      stock = scientific_name,
      common_name = NA,
      scientific_name = scientific_name,
      res = unique(data$resilience),
      start.yr = years[years >= 1950][1],
      end.yr = dplyr::last(years[years >= 1950]),
      cores = cmsy_cores
    )
    if (is.null(ported_cmsy$error)){
      
      cmsy_results <- ported_cmsy$result %>% 
        mutate(c_div_msy  = ct / msy) %>% 
        select(year, bt, b_bmsy, f_fmsy, c_div_msy) %>% 
        rename(depletion = bt, b_div_bmsy = b_bmsy, u_div_umsy = f_fmsy) %>% 
        pivot_longer(-year, names_to = "variable", values_to = "mean") %>% 
        mutate(sd = NA, lower = NA, upper = NA, data = "CMSY")
      
      cmsy_worked <- TRUE 
    } else {
      cmsy_worked <- FALSE
    }
    
    
    # fit heuristic
    
    sar_fmi_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        fmi = data[,c("research", "management", "enforcement","socioeconomics")] %>% colMeans(),
        sar = mean(data$sar),
        initial_state = NA,
        initial_state_cv = 0.2,
        terminal_state = NA,
        terminal_state_cv =  0.05,
        use_heuristics = FALSE,
        b_ref_type = "k"
      )
    
    
    sar_fmi_fit <-
      fit_sraplus(
        driors = sar_fmi_driors,
        include_fit = TRUE,
        model = model,
        engine = "sir",
        estimate_shape = estimate_shape,
        estimate_qslope = estimate_qslope,
        estimate_proc_error = estimate_proc_error
      )
    
 
    sar_fmi_fit <- sar_fmi_fit$results %>%
      filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
      mutate(year = rep(years, each = 4)) %>%
      mutate(data = "SAR & FMI")
    
    # effort fit
    
    effort_years <- which(years >= first_effort_year & !is.na(data$total_effort))
    
    # dyn.load(TMB::dynlib(here::here("src", "sraplus_v2")))
    cpue_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        initial_state = NA,
        fmi = data[,c("research", "management", "enforcement","socioeconomics")] %>% colMeans(),
        sar = mean(data$sar),
        effort = data$total_effort[effort_years],
        effort_years = effort_years,
        q_slope_prior = q_slope_prior,
        q_slope_prior_cv = 0.2,
        b_ref_type = "k"
      )
    
    sfs <- safely(fit_sraplus)
    
    cpue_fit <-
      sfs(
        driors = cpue_driors,
        include_fit = FALSE,
        model = model,
        engine = engine,
        estimate_shape = estimate_shape,
        estimate_qslope = estimate_qslope,
        estimate_proc_error = estimate_proc_error,
        use_baranov = TRUE,
        workers = 4
      )
    
    # browser()
    # plot(cpue_driors$catch[cpue_driors$effort_years] / cpue_driors$effort)
    # plot_prior_posterior(fit = cpue_fit$result, driors = cpue_driors)

    # plot_sraplus(cpue_fit$result, years = cpue_driors$years)
    #
    # cpue = cpue_driors$catch[cpue_driors$effort_years] / (cpue_driors$effort / 10000)
    # #
    # index = cpue_fit$result$results %>%
    #   filter(variable == "index_hat_t")
    #
    # plot(scale(index$mean))
    # lines(scale(cpue))
    #
    # plot_sraplus(cpue_fit$result, years = years)
    # plot_prior_posterior(cpue_fit$result, cpue_driors)
    
    # browser()
    if (is.null(cpue_fit$error)){
      cpue_fit <- cpue_fit$result
    } 
    if (engine == "tmb") {
      cpue_fit_results <- cpue_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
        mutate(year = rep(years, 4)) %>%
        # modify_at(c("mean", "lower", "upper"), exp) %>%
        mutate(data = "cpue") %>%
        mutate(variable = str_replace_all(variable, "log_", ""))
    } else {
      
      cpue_fit_results <- cpue_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
        mutate(year = year - 1 + min(years)) %>%
        mutate(data = "Effort & SAR & FMI")
      
    }
    
    results <- cpue_fit_results %>%
      bind_rows(sar_fmi_fit) %>% {
        if (cmsy_worked){
          bind_rows(., cmsy_results)
        } else {
          .
        }
      }
    
    return(results)
    
  }