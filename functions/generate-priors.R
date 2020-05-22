generate_priors <-
  function(dat,
           use_index = FALSE,
           index_type = "survey",
           use_initial_state = TRUE,
           initial_state_cv = 0.1,
           use_terminal_state = TRUE,
           terminal_state_cv = 0.1,
           use_terminal_u = FALSE,
           terminal_u = NA,
           terminal_u_cv = 0.1,
           known_carry  = FALSE,
           carry_cv = 0.25,
           use_u_priors = FALSE,
           fit_effort = FALSE,
           error_cv = 0.1,
           q = 1e-6,
           b_ref_type = "k",
           f_ref_type = "fmsy",
           index_window = 1,
           effort_window = 1,
           u_window = 1,
           u_freq = 1,
           u_cv = 0.1,
           index_freq = 1,
           effort_freq = 1,
           use_heuristics = FALSE) {
    if (use_initial_state == TRUE) {
      initial_state = dplyr::case_when(b_ref_type == "k" ~ dat$b_rel[1],
                                       TRUE ~ dat$b_v_bmsy[1])

    } else {
      initial_state <- NA
    }

    if (use_terminal_state == TRUE) {
      terminal_state = dplyr::case_when(b_ref_type == "k" ~ last(dat$b_rel),
                                        TRUE ~ last(dat$b_v_bmsy))

    } else {
      terminal_state <- NA
    }

    if (use_heuristics == TRUE) {
      temp <-
        if (dat$catch[1] / max(dat$catch, na.rm = TRUE) < 0.2) {
          0.7
        } else {
          0.4
        }

      initial_state = dplyr::case_when(b_ref_type == "k" ~ temp,
                                       TRUE ~ temp * 2)


      temp_terminal <-
        ifelse((last(dat$catch) / max(dat$catch)) > 0.5, 0.6, 0.2)

      terminal_state = dplyr::case_when(b_ref_type == "k" ~ temp_terminal,
                                        TRUE ~ temp_terminal * 2)


    } # close if use heuristics


    if (known_carry == TRUE) {
      carry <-
        ifelse(all(is.na(dat$TB0)), max(dat$total_biomass), max(dat$TB0))

    } else {
      carry <- NA
    }

    if (use_index == TRUE) {
      index_years <-
        seq((length(dat$index) - floor(index_window * length(dat$index)) + 1), length(dat$index), by = index_freq)

      if (index_type == "survey") {
        index <-  dat$index[index_years]
      } else {
        index_years <- index_years[!is.na(cs$u_v_umsy)]

        index <- NA

        effort <-  cs$u_v_umsy[index_years] * 1e6

        effort_years = index_years
        # index <- dat$approx_cpue[index_years]


      }
    } else {
      index_years <- seq_along(dat$catch)

      index <-  NA
    }

    if (fit_effort == TRUE) {
      effort_years <-
        seq((length(dat$cs_effort) - floor(
          effort_window * length(dat$cs_effort)
        ) + 1), length(dat$cs_effort), by = effort_freq)


    }  else {
      effort_years <- seq_along(dat$catch)

    }

    if (use_u_priors == TRUE) {
      if (f_ref_type == "fmsy") {
        u_years <-
          seq((length(dat$u_v_umsy) - floor(u_window * length(dat$u_v_umsy)) + 1), length(dat$u_v_umsy), by = u_freq)

      } else {
        u_years <-
          seq((
            length(dat$exploitation_rate) - floor(u_window * length(dat$exploitation_rate)) + 1
          ), length(dat$exploitation_rate), by = u_freq)

      }



      if (use_terminal_u == TRUE) {

        u_years <-  NA

        if (f_ref_type == "fmsy") {
          terminal_u <-  dplyr::last(dat$u_v_umsy)

        } else {
          terminal_u <-   dplyr::last(dat$exploitation_rate)

        }

      }

    }  else {
      u_years <- seq_along(dat$catch)

    }

    if (fit_effort == TRUE) {
      temp <- dat$u_v_umsy[effort_years]

      effort <- (temp - mean(temp)) / sd(temp)

    } else {
      effort <-  NA

    }

    if (use_u_priors == TRUE & use_terminal_u == FALSE) {
      if (f_ref_type == "fmsy") {
        u <-  dat$u_v_umsy[u_years]

      } else {
        u <- dat$exploitation_rate[u_years]

      }
    } else {
      u <-  NA
    }
    priors <-
      list(
        carry = carry,
        carry_cv = sqrt(log(carry_cv ^ 2 + 1)),
        terminal_state = terminal_state,
        terminal_state_cv = terminal_state_cv,
        initial_state = initial_state,
        initial_state_cv = initial_state_cv,
        terminal_u = terminal_u,
        terminal_u_cv = u_cv,
        index = index,
        effort = effort,
        u = u,
        index_years = index_years,
        effort_years = effort_years,
        u_years = u_years,
        u_cv = u_cv
      )

    # apply error to priors
    error_prone <- !str_detect(names(priors), "cv|years|anchors")

    priors <-
      map_if(priors, error_prone, ~ .x * exp(rnorm(length(.x), 0, error_cv) - error_cv ^
                                               2 / 2))

    priors$b_ref_type <- b_ref_type

    priors$f_ref_type <- f_ref_type


    return(priors)

  } # close function