calc_error <- function(data, priors,index_type, q = 1e-6) {
  #
  # data <- ram_fits$data[[1]]
  #
  # priors <- ram_fits$priors[[1]]

  if (index_type == "survey") {
    true <-  data$index[priors$index_years]
  } else {
    true <- data$approx_cpue[priors$index_years]

  }

  observed <- priors$index

  index_rmse <- sqrt(mean((true - observed) ^ 2))

  true <-  data$u_v_umsy[priors$u_years]

  observed <- priors$u_v_umsy

  u_rmse <- sqrt(mean((true - observed) ^ 2))

  ref_b <- dplyr::case_when(priors$b_ref_type == "k" ~ last(data$b_rel),
                            TRUE ~ last(data$b_v_bmsy))

  terminal_b_rmse = sqrt((ref_b - priors$terminal_state) ^ 2)


  ref_b <-
    dplyr::case_when(priors$b_ref_type == "k" ~ first(data$b_rel),
                     TRUE ~ first(data$b_v_bmsy))

  initial_b_rmse = sqrt((ref_b - priors$initial_state) ^ 2)

  outs <- ls()[str_detect(ls(), "rmse")]

  out <- map_dfc(outs, ~ data_frame(.x = get(.x))) %>%
    set_names(outs)

  return(out)

}