fit_cmsy <- function(catches, years, scientific_name, priors) {
  genus_species <-
    unique(scientific_name) %>% str_split(" ", simplify = TRUE)

  sss <- purrr::quietly(FishLife::Search_species)

  fish_search <-
    sss(Genus = genus_species[1], Species = genus_species[2])


  if (is.null(fish_search$error)) {
    taxon <- fish_search$result$match_taxonomy[1] %>%
      str_split("_") %>%
      unlist()

  } else{
    taxon <-
      FishLife::sss(Genus = genus_species[1])$result$match_taxonomy[1] %>%
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

  taxon[which(taxon == "")] <- "predictive"

  sp <- sss(
    Class = taxon["Class"],
    Order = taxon["Order"],
    Family = taxon["Family"],
    Genus = taxon["Genus"],
    Species = taxon["Species"],
    ParentChild_gz = FishLife::FishBase_and_RAM$ParentChild_gz
  )$result$match_taxonomy
  ## find means and covariance matrix


  taxa_location <-
    grep(sp[1], FishLife::FishBase_and_RAM$ParentChild_gz[, "ChildName"])[1]

  traits_mean <- FishLife::FishBase_and_RAM$beta_gv[taxa_location,]
  traits_cov <-
    FishLife::FishBase_and_RAM$Cov_gvv[taxa_location, ,]

  params_mvn <-
    c("r", "ln_var")

  mean_mvn <- traits_mean[which(names(traits_mean) %in% params_mvn)]

  sigma_mvn <-
    traits_cov[which(rownames(traits_cov) %in% params_mvn), which(colnames(traits_cov) %in% params_mvn)]

  traits <-
    mvtnorm::rmvnorm(100, mean = mean_mvn, sigma = sigma_mvn)

  if (is.na(priors$initial_state)) {
    priors$initial_state <-
      ifelse((dplyr::first(catches) / max(catches)) < 0.2, 0.7, 0.4)


  }

  if (is.na(priors$terminal_state)) {
    priors$terminal_state <-
      ifelse((dplyr::last(catches) / max(catches)) > 0.5, 0.6, 0.2)

  }

  initial_depletions <-
    pmax(0.1,
         rnorm(100, priors$initial_state, priors$initial_state_cv))

  final_depletions <-
    pmax(0.01,
         rnorm(100, priors$terminal_state, priors$terminal_state_cv))


  if (is.na(priors$carry)) {
    start_k <- c(max(catches) * 5, 50 * max(catches))

  } else {
    start_k <- c(priors$carry / 5, 50 * priors$carry)

  }
  fit <- datalimited::cmsy(
    yr = years,
    ct = catches,
    sig_r = exp(sample(traits[, "ln_var"], 1)) / 2,
    start_r = c(0.5 * min(traits[, "r"]), 2 * max(traits[, "r"])),
    startbio = c(
      quantile(initial_depletions, 0.05),
      quantile(initial_depletions, 0.95)
    ),
    start_k = start_k,
    finalbio = c(
      quantile(final_depletions, 0.05),
      quantile(final_depletions, 0.95)
    )
  )



}