run_biosraplus <- function(iter = 1000,
                           dat,
                           priors,
                           refresh = 0,
                           show_messages = FALSE,
                           sra_model,
                           q = 1e-6) {


  genus_species <-
    unique(dat$scientificname) %>% str_split(" ", simplify = TRUE)

  shh <- quietly(FishLife::Search_species)

  fish_search <-
    shh(Genus = genus_species[1], Species = genus_species[2])

  if (is.null(fish_search$error)) {
    taxon <- fish_search$result$match_taxonomy[1] %>%
      str_split("_") %>%
      unlist()

  } else{
    shh <- quietly(FishLife::Search_species)

    taxon <-
      shhh(Genus = genus_species[1])$result$match_taxonomy[1] %>%
      str_split("_") %>%
      unlist()

  }




  taxon <-
    data_frame(
      Class = taxon[[1]],
      Order = taxon[[2]],
      Family = taxon[[3]],
      Genus = taxon[[4]],
      Species = taxon[[5]]
    )

  sp <- shh(
    Class = taxon["Class"],
    Order = taxon["Order"],
    Family = taxon["Family"],
    Genus = taxon["Genus"],
    Species = taxon["Species"],
    ParentChild_gz = Return$ParentChild_gz
  )$result$match_taxonomy[1]

  taxa_location <- grep(sp, Return$ParentChild_gz[, "ChildName"])[1]

  mean_lh <- Return$beta_gv[taxa_location,]

  cov_lh <- Return$Cov_gvv[taxa_location, ,]
  ##------------ deviates from multivariate normal -------------##
  ## parameters for multvariate normal distribution
  ##  params_mvn <- c("tm", "M", "tmax", "Loo", "K", "Winfinity") #"ln_var", "Loo", "h", "logitbound_h")
  params_mvn <-
    c("intrinsic_growth_rate", "ln_var")
  mean_lh <- mean_lh[which(names(mean_lh) %in% params_mvn)]

  cov_lh <-
    cov_lh[which(rownames(cov_lh) %in% params_mvn), which(colnames(cov_lh) %in% params_mvn)] %>% diag()


  years <- length(dat$catch)

  index_years <- priors$index_years

  effort_years <- priors$effort_years

  if (all(is.na(priors$index))){

    log_index<- rep(0, length(index_years))

  } else {

    log_index <- log(priors$index)
  }

  if (all(is.na(priors$effort))){

    cs_effort <- rep(0, length(effort_years))

  } else {

    cs_effort <- priors$effort
  }

  if (all(is.na(priors$u_v_umsy))){

    u_v_umsy <- rep(0, length(effort_years))

  } else {

    u_v_umsy <- priors$u_v_umsy
  }


  if (length(effort_years) == 1){

    cs_effort = array(cs_effort, dim = 1)

    u_v_umsy = array(u_v_umsy, dim = 1)

    effort_years = array(effort_years, dim = 1)
  }

  if (length(index_years) == 1){

    log_index = array(log_index, dim = 1)

    index_years = array(index_years, dim = 1)

  }
  sra_data <- list(
    nt = years,
    n_index = length(index_years),
    n_effort = length(effort_years),
    index_years = index_years,
    catch_t = dat$catch,
    log_index = log_index,
    growth_prior = mean_lh["intrinsic_growth_rate"],
    growth_cv = sqrt(cov_lh["intrinsic_growth_rate"]),
    sigma_r_prior = exp(mean_lh["ln_var"])/2,
    sigma_r_cv = exp(cov_lh["ln_var"]),
    init_ref = priors$initial_b,
    init_cv = priors$initial_b_cv,
    final_ref = priors$terminal_b,
    final_cv = priors$terminal_b_cv,
    final_u = ifelse(is.na(priors$terminal_u),1,priors$terminal_u),
    final_u_cv = priors$terminal_u_cv,
    fit_index = !all(is.na(priors$index)),
    fit_effort = !all(is.na(priors$effort)),
    fit_final_u = !is.na(priors$terminal_u),
    fit_u = !all(is.na(priors$u_v_umsy)),
    u = u_v_umsy,
    cs_effort = cs_effort,
    ref_type = ifelse(priors$ref_type == "k",0,1),
    effort_years = effort_years
  )

  sra_fit <- rstan::sampling(
    sra_model,
    data = sra_data,
    init = list(list(
      k_mult = 10,
      phi = 0.188,
      log_r = log(sra_data$growth_prior)
    )),
    chains = 1,
    cores = 1,
    iter = iter,
    seed = 42,
    control = list(adapt_delta = 0.8,
                   max_treedepth = 10),
    show_messages = show_messages,
    refresh = refresh
  )


  b <- tidybayes::spread_draws(sra_fit, b_v_bmsy[year])

  u <- tidybayes::spread_draws(sra_fit, u_v_umsy[year])

  msy <- tidybayes::spread_draws(sra_fit, msy)

  pp_index = rstan::extract(sra_fit, "pp_index")[[1]]

  pp_u = rstan::extract(sra_fit, "pp_u")[[1]]

  prior_r = rstan::extract(sra_fit, "prior_r")[[1]]

  prior_sigma_r = rstan::extract(sra_fit, "prior_sigma_r")[[1]]

  divergences <- rstan::get_num_divergent(sra_fit)

  max_treedepth <- rstan::get_num_max_treedepth(sra_fit)

  out <- list(
    b_v_bmsy = b,
    u_v_umsy = u,
    msy = msy,
    divergences = divergences,
    max_treedepth = max_treedepth,
    iter = iter,
    pp_index = pp_index,
    pp_u = pp_u,
    prior_r = prior_r,
    prior_sigma_r = prior_sigma_r
  )


  return(out)

}