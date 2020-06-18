




run_assessments <-
  function(dat,
           priors,
           experiment,
           draws = 1e6,
           assessments = c("sraplus"),
           seed = 42,
           results_name = "foo",
           write_results = FALSE,
           inst_f = TRUE,
           model = "sraplus_tmb",
           engine = "sir",
           keep = 2500,
           Kscale = NA,
           estimate_shape = FALSE,
           estimate_proc_error = TRUE,
           estimate_initial_state = TRUE,
           max_time = Inf,
           refresh = 0) {
    # fits <- vector("list", length(assessments))

    # names(fits) <- assessments

    fits <- list()
    
    if ("sraplus" %in% assessments) {
      driors <-
        sraplus::format_driors(
          taxa = dat$scientificname[[1]],
          catch = dat$catch,
          years = seq_along(dat$catch),
          index = priors$index,
          index_years = priors$index_years,
          terminal_state_cv = priors$terminal_state_cv,
          terminal_state = priors$terminal_state,
          initial_state = pmin(1, priors$initial_state),
          initial_state_cv = priors$initial_state_cv,
          terminal_u = priors$terminal_u,
          terminal_u_cv = priors$terminal_u_cv,
          u = priors$u,
          f_ref_type = priors$f_ref_type,
          u_years = priors$u_years,
          b_ref_type = priors$b_ref_type,
          shape_prior = 1.01,
          u_cv = priors$u_cv)

      if (all(is.na(driors$index)) & (length(driors$u) == 1)) {
        engine <- "sir"

      }

      if (engine != "sir") {
        fits$sraplus <- sraplus::fit_sraplus(
          driors = driors,
          include_fit = TRUE,
          model = model,
          engine = engine,
          draws = draws,
          n_keep = keep,
          estimate_proc_error = estimate_proc_error,
          estimate_shape = estimate_shape,
          estimate_initial_state = estimate_initial_state,
          refresh = refresh
        )

      } else
      {
        fits$sirplus <- sraplus::fit_sraplus(
          driors = driors,
          include_fit = TRUE,
          model = model,
          engine = engine,
          draws = draws,
          n_keep = keep,
          estimate_proc_error = estimate_proc_error,
          estimate_shape = estimate_shape,
          estimate_initial_state = estimate_initial_state
        )

      }


    }

    if ("colesra" %in% assessments) {
      genus_species <-
        unique(dat$scientificname) %>% str_split(" ", simplify = TRUE)

      shhh <- purrr::quietly(FishLife::Search_species)

      fish_search <-
        shhh(Genus = genus_species[1], Species = genus_species[2])$result

      taxon <- fish_search$match_taxonomy[1] %>%
        str_split("_") %>%
        unlist()


      taxon <-
        data_frame(
          Class = taxon[[1]],
          Order = taxon[[2]],
          Family = taxon[[3]],
          Genus = taxon[[4]],
          Species = taxon[[5]]
        )

      cole_priors <- priors

      cole_priors$ustatus.mean <- cole_priors$fstatus.mean

      cole_priors$ustatus.sd <- cole_priors$fstatus.sd

      cole_priors$ustatus.dist <- cole_priors$fstatus.dist

      cole_priors$fstatus.dist <- NULL

      fits$colesra <-
        colesraplus::run.SIR(
          nrep = draws,
          Catch = dat$catch,
          Taxon = taxon,
          penalties = cole_priors,
          years = dat$year,
          Kscale = Kscale
        )
    }

    if ("cmsy" %in% assessments) {
      safe_cmsy <- safely(fit_cmsy)
      cmsy_fit <- safe_cmsy(
        catches = dat$catch,
        years = dat$year,
        scientific_name = unique(dat$scientificname),
        priors = priors
      )

      if (is.null(cmsy_fit$error)) {
        fits$cmsy <- cmsy_fit$result

      } else {
        fits <- purrr::compact(fits) # remove null list elements

      }


    }

    if (write_results == TRUE) {
      saveRDS(fits,
              file = here::here(
                "results",
                results_name,
                "experiments",
                paste0('experiment-', experiment, ".rds")
              ))

      fits <- NA

    }

    out <- fits


  }