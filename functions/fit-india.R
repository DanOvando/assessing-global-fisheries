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
           cores = 1){
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
      filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
      mutate(year = rep(years, each = 4)) %>%
      mutate(data = "heuristic")


    if (!is.na(sar)) {

      sar_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(catches),
          sar = sar,
          f_ref_type = "fmsy")


sar_fit <-
        fit_sraplus(driors = sar_driors,
                    include_fit = include_fit,
                    model = model,
                    engine = "sir",
                    estimate_proc_error = estimate_proc_error,
                    estimate_shape = estimate_proc_error,
                    draws = draws)

      # plot_sraplus(sar_fit)

    sar_fit_results <- sar_fit$results %>%
      filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
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
      sfs(driors = cpue_driors,
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
          workers = 4)

    if (is.null(cpue_fit$error)){
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
        filter(variable %in% c("log_b_div_bmsy", "log_u_div_umsy", "log_c_div_msy","log_depletion")) %>%
        mutate(year = rep(years, 4)) %>%
        modify_at(c("mean", "lower", "upper"), exp) %>%
        mutate(data = "cpue") %>%
        mutate(variable = str_replace_all(variable, "log_", ""))
    } else {

      cpue_fit_results <- cpue_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
        mutate(year = rep(years,4)) %>%
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
      }


    return(results)

  }