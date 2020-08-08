fit_continent_examples <-
  function(scientific_name,
           data,
           model = "sraplus_tmb",
           engine = "tmb",
           fmi_f_reg = NA,
          estimate_shape = TRUE,
          estimate_proc_error = TRUE,
          estimate_qslope = FALSE,
          use_baranov = TRUE,
          first_effort_year = 1950,
          chains = 1,
          cores = 1,
          cmsy_cores = 4,
          q_slope_prior = 0.025) {
    # heuristic example
    # i = 2
    # engine = "tmb"
    # model = "sraplus_tmb"
    # data <- huh$data[[1]]
    # scientific_name <- huh$scientific_name[[1]]
    # catches <- data$total_catch

    # estimate_proc_error = FALSE
    # estimate_qslope = TRUE
    # estimate_shape = TRUE
    # use_baranov = TRUE
    # first_effort_year = 1950

    
    catches <- data$total_catch
    
    years <- data$year
    
 
  
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
    
    
    # fit heuristic
    
    basic_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        initial_state = NA,
        initial_state_cv = 0.2,
        terminal_state = NA,
        terminal_state_cv =  0.05,
        use_heuristics = TRUE,
        b_ref_type = "k"
      )


    basic_fit <-
      fit_sraplus(
        driors = basic_driors,
        include_fit = TRUE,
        model = model,
        engine = "sir",
        estimate_shape = estimate_shape,
        estimate_qslope = estimate_qslope,
        estimate_proc_error = estimate_proc_error
      )

    # basic_fit_results <- basic_fit$results %>%
    #   filter(variable != "b_t", variable != "index_hat_t") %>%
    #   mutate(variable = case_when(
    #     variable == "b_bmsy_t" ~ "b",
    #     variable == "c_msy_t" ~ "c",
    #     variable == "dep_t" ~ "dep",
    #     variable == "u_umsy_t" ~ "u",
    #     TRUE ~ variable
    #   )) %>%
    #   mutate(year = min(years) + year - 1,
    #          data = "heuristic")

    basic_fit_results <- basic_fit$results %>%
      filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
      mutate(year = rep(years, each = 4)) %>%
      mutate(data = "heuristic")

    # effort fit

    effort_years <- which(years >= first_effort_year & !is.na(data$total_effort))

    # dyn.load(TMB::dynlib(here::here("src", "sraplus_v2")))
    cpue_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        effort = data$total_effort[effort_years],
        effort_years = effort_years,
        q_slope_prior = q_slope_prior,
        q_slope_prior_cv = 0.2,
        b_ref_type = "k"
      )

    sfs <- safely(fit_sraplus)

    # cpue_driors$initial_state <- 1
    #
    # cpue_driors$initial_state_cv <- 0.1
    #
    # cpue_driors$q_slope_prior <- 0.025
    #
    # cpue_driors$q_slope_prior_cv <- 0.5
    #
    # cpue_driors$shape_prior_cv <- 0.5
    cpue_fit <-
      sfs(
        driors = cpue_driors,
        include_fit = FALSE,
        model = model,
        engine = engine,
        estimate_shape = estimate_shape,
        estimate_qslope = estimate_qslope,
        estimate_proc_error = estimate_proc_error,
        use_baranov = use_baranov,
        workers = 4
      )
    #
    # plot(cpue_driors$catch[cpue_driors$effort_years] / cpue_driors$effort)
    #
    # plot_prior_posterior(fit = cpue_fit$result, driors = cpue_driors)
    #
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
    } else {

      # cpue_fit <-
      #   fit_sraplus(
      #     driors = cpue_driors,
      #     include_fit = TRUE,
      #     model = model,
      #     engine = "stan",
      #     estimate_shape = estimate_shape,
      #     estimate_qslope = estimate_qslope,
      #     estimate_proc_error = estimate_proc_error,
      #     chains = chains,
      #     cores = cores
      #   )
      #
      # # plot_sraplus(cpue_fit$result, years = years)
      #
      #
      # engine = "stan"
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
        mutate(data = "cpue")

    }

    results <- cpue_fit_results %>%
      bind_rows(basic_fit_results) %>% {
        if (cmsy_worked){
          bind_rows(., cmsy_results)
        } else {
          .
        }
      }

    return(results)

  }