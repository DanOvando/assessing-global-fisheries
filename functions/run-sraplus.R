run_sraplus <-
  function(dat,
           priors,
           nrep = 2000,
           seed = 42,
           Kscale = 5,
           inst_f  = TRUE) {
    genus_species <-
      unique(dat$scientificname) %>% str_split(" ", simplify = TRUE)

    sss <- quietly(FishLife::Search_species)

    fish_search <-
      sss(Genus = genus_species[1], Species = genus_species[2])

    if (is.null(fish_search$error)) {
      taxon <- fish_search$result$match_taxonomy[1] %>%
        str_split("_") %>%
        unlist()

    } else{
      shhh <- quietly(FishLife::Search_species)

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



    fit <-
      sraplus::run.SIR(
        nrep = nrep,
        Catch = dat$catch,
        Taxon = taxon,
        penalties = priors,
        years = dat$year,
        seed = seed,
        Kscale = Kscale,
        inst_f = inst_f
      )

    # stop <- Sys.time()

    # time_took <- stop - start
    # out <- data_frame(keepers = n_distinct(fit$Keepers), time_took = time_took)

  }