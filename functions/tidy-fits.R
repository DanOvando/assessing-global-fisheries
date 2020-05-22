tidy_fits <- function(fit, name, years) {

  if (name == "cmsy") {
    b_v_bmsy <- (fit$biomass / fit$bmsy) %>%
      as_data_frame() %>%
      select_at(1:(ncol(.) - 1)) %>%
      set_names(as.character(fit$bbmsy$year)) %>%
      mutate(draw = 1:nrow(.)) %>%
      gather(year, b_v_bmsy, -draw) %>%
      mutate(year = as.numeric(year))

    depletion <- (fit$biomass / fit$theta$k) %>%
      as_data_frame() %>%
      select_at(1:(ncol(.) - 1)) %>%
      set_names(as.character(fit$bbmsy$year)) %>%
      mutate(draw = 1:nrow(.)) %>%
      gather(year, depletion, -draw) %>%
      mutate(year = as.numeric(year))

    msys <- data_frame(msy = fit$msy, draw = 1:length(fit$msy))

    out <-
      data_frame(year = fit$bbmsy$year,
                 catch = fit$bbmsy$catch) %>%
      mutate(msy = list(msys)) %>%
      unnest() %>%
      left_join(b_v_bmsy, by = c("year", "draw")) %>%
      left_join(depletion, by = c("year", "draw")) %>%
      mutate(u_v_umsy = (catch / msy) / b_v_bmsy) %>%
      select(year, draw, b_v_bmsy, u_v_umsy, depletion) %>%
      mutate(fit_name = name)


  } else if (name == "biosraplus") {


    b_v_bmsy <- fit$b_v_bmsy %>%
      select(year, .draw, b_v_bmsy) %>%
      rename(draw = .draw)


    u_v_umsy <- fit$u_v_umsy %>%
      select(year, .draw, u_v_umsy) %>%
      rename(draw = .draw)


    out <- b_v_bmsy %>%
      left_join(u_v_umsy, by = c("year", "draw"))  %>%
      mutate(fit_name = name) %>%
      ungroup() %>%
      mutate(year = year - 1 + min(observed$year))



  } else if (name == "sirplus") {


    b_v_bmsy <- fit$fit %>%
      filter(variable == "b_bmsy_t") %>%
      select(year, draw, value) %>%
        rename(b_v_bmsy = value)


    u_v_umsy <- fit$fit %>%
      filter(variable == "u_umsy_t") %>%
      select(year, draw, value) %>%
      rename(u_v_umsy = value)

    depletion <- fit$fit %>%
      filter(variable == "dep_t") %>%
      select(year, draw, value) %>%
      rename(depletion = value)



    out <- b_v_bmsy %>%
      left_join(u_v_umsy, by = c("year", "draw"))  %>%
      left_join(depletion, by = c("year", "draw"))  %>%
      mutate(fit_name = name) %>%
      ungroup() %>%
      mutate(year = min(years) + year - 1)


  }


  else {

    samps <- 4000
    b_v_bmsy <- fit$results %>%
      filter(variable == "log_b_div_bmsy") %>%
      # filter(variable == "b_v_bmsy") %>%
      mutate(samps = map2(mean, sd, ~data_frame(draw = 1:samps, samp = pmin(5,exp(rnorm(samps,.x,.y))))),
             year = years) %>%
      unnest() %>%
      select(year, draw, samp) %>%
      rename(b_v_bmsy = samp)

    depletion <- fit$results %>%
      filter(variable == "log_depletion") %>%
      # filter(variable == "dep_t") %>%
      mutate(samps = map2(mean, sd, ~data_frame(draw = 1:samps, samp = pmin(5,exp(rnorm(samps,.x,.y))))),
             year = years) %>%
      unnest() %>%
      select(year, draw, samp) %>%
      rename(depletion = samp)


    u_v_umsy <- fit$results %>%
      filter(variable == "log_u_div_umsy") %>%
      mutate(samps = map2(mean, sd, ~data_frame(draw = 1:samps, samp = pmin(5,exp(rnorm(samps,.x,.y))))),
             year = years) %>%
      unnest() %>%
      select(year, draw, samp) %>%
      rename(u_v_umsy = samp)

    out <- b_v_bmsy %>%
      left_join(u_v_umsy, by = c("year", "draw"))  %>%
      left_join(depletion, by = c("year", "draw"))  %>%
      mutate(fit_name = name) %>%
      ungroup()

  }

  return(out)
}
