fit_india <-
  function(scientific_name,
           catches,
           years,
           index,
           index_years,
           support_data,
           model = "sraplus_tmb",
           engine = "tmb",
           include_fit = FALSE,
           country = "India",
           sar = NA,
           min_year = 1990,
           estimate_shape = FALSE,
           estimate_proc_error = TRUE,
           n_keep = 2500,
           draws = 1e6,
           chains = 1,
           refresh = 500,
           cores = 1,
           cmsy_cores = 8) {
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
    if (is.null(ported_cmsy$error)) {
      cmsy_results <- ported_cmsy$result %>%
        mutate(c_div_msy  = ct / msy) %>%
        select(year, bt, b_bmsy, f_fmsy, c_div_msy) %>%
        rename(depletion = bt,
               b_div_bmsy = b_bmsy,
               u_div_umsy = f_fmsy) %>%
        pivot_longer(-year, names_to = "variable", values_to = "mean") %>%
        mutate(
          sd = NA,
          lower = NA,
          upper = NA,
          data = "cmsy"
        )
      
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
    
    
    # heuristic example
    
    basic_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches,
        years = seq_along(years),
        initial_state = 1,
        initial_state_cv = 0.2,
        terminal_state = 0.5,
        terminal_state_cv = 0.05,
        use_heuristics = TRUE,
        b_ref_type = "k"
      )
    
    basic_fit <-
      fit_sraplus(
        driors = basic_driors,
        include_fit = include_fit,
        model = model,
        engine = "sir",
        draws = draws
      )
    
    
    basic_fit_results <- basic_fit$results %>%
      filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
      mutate(year = rep(years, each = 4)) %>%
      mutate(data = "heuristic")
    
    
    if (!is.na(sar)) {
      sar_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(catches),
          sar = sar,
          f_ref_type = "fmsy"
        )
      
      
      sar_fit <-
        fit_sraplus(
          driors = sar_driors,
          include_fit = include_fit,
          model = model,
          engine = "sir",
          estimate_proc_error = estimate_proc_error,
          estimate_shape = estimate_proc_error,
          draws = draws
        )
      
      # plot_sraplus(sar_fit)
      
      sar_fit_results <- sar_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
        mutate(year = rep(years, each = 4)) %>%
        mutate(data = "sar")
      
      sar_fit_worked <- TRUE
      
      
    }
    
    # index fit
    
    index_years <- which(years %in% index_years)
    
    cpue_driors <-
      format_driors(
        taxa = scientific_name,
        catch = catches + 1e-3,
        years = seq_along(years),
        index = index ,
        index_years = index_years,
        sar = sar,
        initial_state = 1,
        growth_rate_prior_cv = .5,
        initial_state_cv = 0.05,
        sigma_ratio_prior = 0.5,
        sigma_ratio_prior_cv = 0.05,
        f_ref_type = "fmsy"
      )
    
    sfs <- safely(fit_sraplus)
    
    cpue_fit <-
      sfs(
        driors = cpue_driors,
        include_fit = include_fit,
        model = model,
        engine = engine,
        estimate_shape = estimate_shape,
        estimate_proc_error = estimate_proc_error,
        estimate_qslope = FALSE,
        n_keep = n_keep,
        chains = chains,
        refresh = refresh,
        cores = cores,
        workers = 4
      )
    
    if (is.null(cpue_fit$error)) {
      cpue_fit <- cpue_fit$result
    } else {
      # sraplus::get_tmb_model(model_name = model)
      #
      # cpue_fit <-
      #   fit_sraplus(driors = cpue_driors,
      #       include_fit = TRUE,
      #       model = model,
      #       engine = "stan",
      #       estimate_shape = estimate_shape,
      #       estimate_proc_error = estimate_proc_error,
      #       estimate_qslope = FALSE,
      #       chains = chains,
      #       cores = cores)
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
    } else {
      cpue_fit_results <- cpue_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy", "depletion")) %>%
        mutate(year = rep(years, 4)) %>%
        mutate(data = "cpue")
    }
    
    
    results <- cpue_fit_results %>%
      bind_rows(basic_fit_results) %>% {
        if (sar_fit_worked == TRUE) {
          bind_rows(., sar_fit_results)
        }
        else{
          .
        }
      } %>% {
        if (cmsy_worked) {
          bind_rows(., cmsy_results)
        } else {
          .
        }
      }
  
return(results)

}