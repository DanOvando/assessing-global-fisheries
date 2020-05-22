generate_old_priors <-
  function(dat,
           initial_dep = NA,
           initial_dep_cv = 0.82,
           terminal_b = TRUE,
           terminal_b_cv = 0.82,
           terminal_f = TRUE,
           terminal_f_cv = 0.82,
           f_dist = 2,
           terminal_u = FALSE,
           terminal_u_cv = 0.82,
           u_dist = 2,
           pchange_effort = TRUE,
           pchange_effort_cv = 0.42,
           pchange_effort_dist = 2,
           carry = TRUE,
           carry_cv = 0.42,
           carry_dist = 2,
           initial_dist = 2,
           b_dist = 2,
           default_initial_dep = 1,
           error_cv = 0
  ) {

    if (is.na(initial_dep)) {
      initial_dep <- default_initial_dep
    }


    if (terminal_b == TRUE) {
      bstatus_mean <-
        ifelse(b_dist == 2, log(dat$b_v_bmsy[dat$year == max(dat$year)]), dat$b_v_bmsy[dat$year == max(dat$year)])

    } else{
      bstatus_mean <-  NA

      b_dist <- NULL

    }

    if (terminal_u == TRUE) {
      ustatus_mean <-
        ifelse(u_dist == 2, log(dat$u_v_umsy[dat$year == max(dat$year)]), dat$u_v_umsy[dat$year == max(dat$year)])

    } else{
      ustatus_mean <-  NA

      u_dist <- NULL

    }

    if (terminal_f == TRUE) {
      fstatus_mean <-
        ifelse(f_dist == 2, log(dat$u_v_umsy[dat$year == max(dat$year)]), dat$u_v_umsy[dat$year == max(dat$year)])

    } else{
      fstatus_mean <-  NA

      f_dist <- NULL

    }

    if (pchange_effort == TRUE) {
      pchange_effort_mean <-
        dplyr::case_when(pchange_effort_dist == 2 ~ log(dat$pchange_effort[dat$year < max(dat$year)]),
                         TRUE ~ dat$pchange_effort[dat$year < max(dat$year)])

    } else{
      pchange_effort_mean <-  NA

      pchange_effort_dist <- NULL

    }


    carry_mean <- NA

    if (carry == TRUE) {
      carry_mean <- case_when(carry_dist == 2 ~ log(unique(dat$SSB0)),
                              TRUE ~ unique(dat$SSB0))

    }

    if (carry == FALSE | is.na(carry_mean) | is.null(carry_mean)) {

      carry_mean <- NA

      carry_dist <- NULL

    }
    pen <-
      list(
        carry.mean = carry_mean,
        carry.sd = sqrt(log(carry_cv ^ 2 + 1)),
        carry.dist = carry_dist,
        bstatus.mean = bstatus_mean,
        bstatus.sd = sqrt(log(terminal_b_cv ^ 2 + 1)),
        bstatus.dist = b_dist,
        fstatus.mean = fstatus_mean,
        fstatus.sd = sqrt(log(terminal_f_cv ^ 2 + 1)),
        fstatus.dist = f_dist,
        ustatus.mean = ustatus_mean,
        ustatus.sd = sqrt(log(terminal_u_cv ^ 2 + 1)),
        ustatus.dist = u_dist,
        initial.mean = log(initial_dep) ,
        initial.sd = sqrt(log(initial_dep_cv ^ 2 + 1)),
        initial.dist = initial_dist,
        pchange_effort_mean = pchange_effort_mean ,
        pchange_effort_sd = sqrt(log(pchange_effort_cv ^ 2 + 1)),
        pchange_effort_dist = pchange_effort_dist
      )

    # apply error to priors
    error_prone <- str_detect(names(pen),"mean")

    pen<- map_if(pen, error_prone, ~ .x + rnorm(length(.x),0, error_cv * mean(.x)))


  } # close function