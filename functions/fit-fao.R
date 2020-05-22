fit_fao <-
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
           chains = 1,
           cores = 1,
           reg_cv = NA,
           draws = 2e6,
           refresh = 0) {
#
#     data <- test$data[[1]]

    ogengine <- engine

    scientific_name <- unique(data$species)

    catches <- data$catch

    years <- data$year

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


    # sar fit

    if (!is.na(last(data$sar))) {
      used_sar <- TRUE

      sar_driors <-
        sraplus::format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(catches),
          initial_state = basic_driors$initial_state,
          initial_state_cv = basic_driors$initial_state_cv,
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
          draws = draws
        )

      sar_fit_results <- sar_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
        mutate(year = rep(years, each = 4)) %>%
        mutate(data = "sar")

    } else {
      used_sar <- FALSE
    }


    # fmi fit

    fmi_dat <- support_data$mean_regional_isscaap_fmi %>%
      filter(fao_area_code == unique(as.numeric(data$fao_area)),
             isscaap_group == unique(data$isscaap_group))


    fmi <- c(research = fmi_dat$research, management = fmi_dat$management, socioeconomics = fmi_dat$socioeconomics, enforcement = fmi_dat$enforcement)

    if (nrow(fmi_dat) == 0){

      fmi_dat <- support_data$mean_regional_isscaap_fmi %>%
        filter(fao_area_code == unique(as.numeric(data$fao_area))) %>%
        pivot_longer(cols = c("research", "management", "socioeconomics", "enforcement"),names_to = "metric", values_to = "value") %>%
        group_by(fao_area_code, metric) %>%
        summarise(mean_value = mean(value)) %>%
        pivot_wider(names_from = metric, values_from = mean_value)

      fmi <- c(research = fmi_dat$research, management = fmi_dat$management, socioeconomics = fmi_dat$socioeconomics, enforcement = fmi_dat$enforcement)

    }

     if (nrow(fmi_dat) == 0){

       fmi <- NA

     }


    if (!all(is.na(fmi))) {

      # term_b_guess <- pmax(0.05,(2.5 - exp(c(last(data$mean_fmi_ppf)))))

      fmi_driors <-
        format_driors(
          taxa = scientific_name,
          catch = catches,
          years = seq_along(years),
          initial_state = basic_driors$initial_state,
          initial_state_cv = 0.2,
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
          engine = "sir",
          draws = draws
        )

      fmi_fit_results <- fmi_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
        mutate(year = min(years) + year - 1) %>%
        mutate(data = "fmi")

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
          initial_state = basic_driors$initial_state,
          initial_state_cv = 0.2,
          isscaap_group = unique(data$isscaap_group)
        )
      
      u_fit <-
        fit_sraplus(
          driors = u_driors,
          include_fit = include_fit,
          model = model,
          engine = "sir",
          draws = draws
        )


      u_fit_results <- u_fit$results %>%
        filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
        mutate(year = min(years) + year - 1) %>%
        mutate(data = "u_umsy")

    } else {
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
          initial_state = basic_driors$initial_state,
          initial_state_cv = 0.2,
          isscaap_group = unique(data$isscaap_group),
          q_slope_prior = q_slope_prior
        )

      sfs <- safely(fit_sraplus)

      cpue_fit <-
        sfs(driors = cpue_driors,
            include_fit = include_fit,
            model = model,
            engine = engine,
            estimate_proc_error = estimate_proc_error,
            estimate_qslope = estimate_qslope,
            workers = cores,
            refresh = refresh)

      if (is.null(cpue_fit$error)){
        cpue_fit <- cpue_fit$result
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
          filter(variable %in% c("log_b_div_bmsy", "log_u_div_umsy", "log_c_div_msy","log_depletion")) %>%
          mutate(year = rep(years, 4)) %>%
          modify_at(c("mean", "lower", "upper"), exp) %>%
          mutate(data = "cpue") %>%
          mutate(variable = str_replace_all(variable, "log_", ""))
      } else {

        cpue_fit_results <- cpue_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
          mutate(year = year - 1 + min(years)) %>%
          mutate(data = "cpue")

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
          initial_state = 1,
          initial_state_cv = 0.2,
          fmi = fmi,
          f_ref_type = "fmsy" ,
          fmi_cv = reg_cv,
          sar = unique(data$sar),
          isscaap_group = unique(data$isscaap_group),
          q_slope_prior = q_slope_prior
        )

      sfs <- safely(fit_sraplus)

      cpue_plus_fit <-
        sfs(driors = cpue_plus_driors,
            include_fit = include_fit,
            model = model,
            engine = engine,
            estimate_proc_error = estimate_proc_error,
            estimate_qslope = estimate_qslope,
            workers = cores,
            refresh = refresh)

      if (is.null(cpue_plus_fit$error)){
        cpue_plus_fit <- cpue_plus_fit$result
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
          filter(variable %in% c("log_b_div_bmsy", "log_u_div_umsy", "log_c_div_msy","log_depletion")) %>%
          mutate(year = rep(years, 4)) %>%
          modify_at(c("mean", "lower", "upper"), exp) %>%
          mutate(data = "cpue-plus") %>%
          mutate(variable = str_replace_all(variable, "log_", ""))
      } else {

        cpue_plus_fit_results <- cpue_plus_fit$results %>%
          filter(variable %in% c("b_div_bmsy", "u_div_umsy", "c_div_msy","depletion")) %>%
          mutate(year = year - 1 + min(years)) %>%
          mutate(data = "cpue-plus")

      }

    } else {
      used_cpue_plus <- FALSE
    }

    results <- basic_fit_results %>% {
      if (used_sar == TRUE) {
        bind_rows(., sar_fit_results)
      }
      else{
        .
      }
    } %>% {
      if (used_cpue == TRUE) {
        bind_rows(., cpue_fit_results)
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
          bind_rows(., cpue_plus_fit_results)
        }
        else{
          .
        }
      }

    # write(wtf,file = "wtf.txt", append = TRUE)

    # results %>%
    #   filter(data == "sar") %>%
    #   ggplot(aes(year, mean, color = data)) +
    #   geom_line() +
    #   facet_wrap(~variable, scales = "free_y")

    return(results)

  }