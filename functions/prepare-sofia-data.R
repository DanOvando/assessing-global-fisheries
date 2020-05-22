prepare_sofia_data <- function(min_years_catch = 20,
                               crazy_b = 4,
                               crazy_u = 5,
                               lookup_fmi_names = FALSE,
                               q = 1e-6) {
  # min_years_catch = 20
  # crazy_b = 4
  # crazy_u = 5
  # lookup_fmi_names = FALSE
  
  # load(here::here("data", "Return.Rdata"))
  
  
  if (file.exists(here("data", "ram.zip")) == FALSE) {
    # for now storing in my google drive... will need to put this in a better and public location
    
    download.file("https://www.dropbox.com/s/jpgz0a5s5of3qev/RAM%20v4.491%20Files%20(1-14-20).zip?dl=1", destfile = here("data","ram.zip"))
    
    unzip(here("data","ram.zip"), exdir = here("data"))
  }
  
  
  ram_dirs <- list.files(here("data","ram"))
  
  ram_dirs <- ram_dirs[str_detect(ram_dirs,"RAM")]
  
  ram_files <- list.files(here("data","ram",ram_dirs), recursive = TRUE)
  
  ram_files <- ram_files[str_detect(ram_files,".RData")]
  
  ram_files <- ram_files[str_detect(ram_files,"Model Fit")]
  
  load(here("data","ram",ram_dirs,ram_files[1]))
  
  # process data ------------------------------------------------------------
  
  stock <- stock %>%
    left_join(area, by = "areaid")
  # catches
  ram_catches <- tcbest.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, catch,-year)
  
  # B/Bmsy
  ram_b_v_bmsy <- divbpref.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    tibble() %>%
    gather(stockid, b_v_bmsy,-year)
  
  
  # U/Umsy
  ram_u_v_umsy <- divupref.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, u_v_umsy,-year)
  
  # Effort
  ram_effort <- effort.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, effort,-year)
  
  # biomass
  
  
  ram_total_biomass <- tbbest.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, total_biomass,-year)
  
  # ssb
  
  ram_ss_biomass <- ssb.data %>%
    mutate(year = rownames(.) %>% as.integer()) %>%
    as_data_frame() %>%
    gather(stockid, ss_biomass,-year)
  
  
  ram_exp_rate <- ram_catches %>%
    left_join(ram_total_biomass, by = c("stockid", "year")) %>%
    mutate(exploitation_rate = catch / total_biomass) %>%
    select(-catch, -total_biomass)
  
  # ram_exp_rate <- erbest.data %>%
  #   mutate(year = rownames(.) %>% as.integer()) %>%
  #   as_data_frame() %>%
  #   gather(stockid, exploitation_rate, -year)
  
  # put it together
  
  ram_data <- ram_catches %>%
    left_join(bioparams_values_views, by = "stockid") %>%
    left_join(ram_b_v_bmsy, by = c("stockid", "year")) %>%
    left_join(ram_u_v_umsy, by = c("stockid", "year")) %>%
    left_join(ram_exp_rate, by = c("stockid", "year")) %>%
    left_join(ram_effort, by = c("stockid", "year")) %>%
    left_join(ram_total_biomass, by = c("stockid", "year")) %>%
    left_join(ram_ss_biomass, by = c("stockid", "year")) %>%
    left_join(stock, by = "stockid") %>%
    select(stockid, scientificname, commonname, everything())
  
  
  # create new variables
  
  ram_data <- ram_data %>%
    mutate(tb_v_tb0 = total_biomass / TB0,
           ssb_v_ssb0 = ss_biomass / SSB0)
  
  # filter data -------------------------------------------------------------
  
  
  # for now, only include continuous catch series
  
  ram_data <- ram_data %>%
    filter(is.na(catch) == FALSE) %>%
    group_by(stockid) %>%
    mutate(delta_year = year - lag(year)) %>%
    mutate(delta_year = case_when(year == min(year) ~ as.integer(1),
                                  TRUE ~ delta_year)) %>%
    mutate(missing_gaps = any(delta_year > 1)) %>%
    filter(missing_gaps == FALSE) %>%
    mutate(n_years = n_distinct(year)) %>%
    filter(n_years >= min_years_catch) %>%
    filter(all(b_v_bmsy < crazy_b, na.rm = TRUE),
           all(u_v_umsy < crazy_u, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(stockid) %>%
    mutate(
      has_tb0 = !all(is.na(TB0)),
      has_tb = !all(is.na(total_biomass)),
      first_catch_year = year[which(catch > 0)[1]]
    ) %>%
    filter(year >= first_catch_year) %>%
    mutate(
      pchange_effort = lead(u_v_umsy) / (u_v_umsy + 1e-6),
      cs_effort = (u_v_umsy - mean(u_v_umsy)) / sd(u_v_umsy),
      index = total_biomass * q,
      approx_cpue = catch / (u_v_umsy / q + 1e-3),
      b_rel = dplyr::case_when(
        has_tb0 ~ total_biomass / max(TB0),
        has_tb ~ total_biomass / max(total_biomass),
        TRUE ~ b_v_bmsy / 2
      )
    ) %>%
    mutate(approx_cpue = pmin(quantile(approx_cpue, 0.9, na.rm = TRUE), approx_cpue)) %>%
    ungroup()
  
  ram_b_plot <- ram_data %>%
    ggplot(aes(x = year, y = b_v_bmsy)) +
    geom_bin2d() +
    scale_fill_viridis_c()
  
  kobe_panels <- ram_data %>%
    filter(year >= 1950) %>%
    mutate(year_block = plyr::round_any(year, 10, floor)) %>%
    ggplot(aes(x = b_v_bmsy, y = u_v_umsy)) +
    geom_bin2d(binwidth = c(0.5, 0.5)) +
    facet_wrap(~ year_block) +
    scale_fill_viridis_c()
  
  # kobe_animation <- ram_data %>%
  #   ggplot(aes(x = b_v_bmsy, y = u_v_umsy)) +
  #   geom_bin2d() +
  #   transition_states(factor(year),
  #                     transition_length = 2,
  #                     state_length = 1) +
  #   scale_fill_viridis_c() +
  #   labs(title = "Year: {closest_state}")
  
  
  # load other data ---------------------------------------------------------
  
  # load prices
  
  prices <<-
    readr::read_csv(here::here("data", "Exvessel Price Database.csv")) %>%
    janitor::clean_names() %>%
    rename(scientificname = scientific_name) %>%
    mutate(log_exvessel = log(exvessel)) %>%
    group_by(asfis_species, pooled_commodity, group_for_pairing) %>%
    mutate(lag_exvessel = lag(exvessel)) %>%
    ungroup() %>%
    group_by(scientificname, year) %>%
    summarise(exvessel = mean(exvessel, na.rm = TRUE)) %>%
    group_by(scientificname) %>%
    mutate(lag_exvessel = lag(exvessel))
  
  ram_data <- ram_data %>%
    left_join(prices, by = c("scientificname", "year"))
  
  # fao and effort data
  
  fao_to_effort <-
    read_csv(here::here("data", "fao-to-bell-region.csv")) %>%
    rename(bell_region = region)
  
  country_to_fao <-
    read_csv(here::here("data", "country-to-fao-area.csv")) %>%
    unique() %>%
    janitor::clean_names() %>%
    rename(fao_fishing_area = fishing_area_fao_major_fishing_area_1) %>%
    left_join(fao_to_effort, by = "fao_fishing_area")
  
  
  fishstat_dat <-
    readr::read_csv(
      here::here("data", "fao_capture_1950-2016.csv")
    ) %>%
    janitor::clean_names() %>% # clean up the names to something usable
    tibble::rownames_to_column(var = "id") %>% # add an id column
    dplyr::filter(
      !stringr::str_detect(tolower(country_country), "total"),!stringr::str_detect(country_country, "\\d"),!stringr::str_detect(
        tolower(production_source_detailed_production_source),
        "aquaculture"
      )
    ) # get rid of totals rows, numbered rows, and aquaculture datta
  
  
  fishstat_dat <- fishstat_dat %>% # convert the individual year columns to a tidy format
    tidyr::pivot_longer(
      cols = starts_with("x"), # noticed that all year columns started with x
      names_prefix = "x", # gets rid of that x in front of the years
      # names_ptypes = list(year = integer()), # convert each column years to an integer after stripping off the x
      names_to = "year",# the variable the former column names will go into
      values_to = "capture" # the variable the capture values will go ti
    ) %>%
    mutate(
      year = as.integer(year),
      capture_notes = str_replace_all(capture, "\\d", ''), #pull out the text
      capture_numbers = str_replace_all(capture, "\\D", "") %>% as.numeric(), #pull out the numbers and convert the numerics. Note the use of a nested pipe! R is cool.
      clean_capture = as.numeric(capture) # SERIOUSLY awful step to accounnt for the like 10 stocks entered as 70000000e-7.... as.numeric will keep any clean numbers and convert any numbers + text to NA
    ) %>%
    mutate(capture = ifelse(test = is.na(clean_capture), yes = capture_numbers, no = clean_capture)) %>%  # keep any stock that was a clean number (which catches the 7e-7, and add in anything else that was like 500 F
    select(-capture_numbers,-clean_capture)
  #
  
  #serious hack, very dangerous
  names <-
    c("id",
      "country",
      "economic_class",
      "common_name",
      "family",
      "scientific_name",
      "asfis_species_code",
      "isscaap_group",
      "isscaap_number",
      "order",
      "fao_area",
      "fao_production_region",
      "fao_area_code",
      "bad",
      "production_source",
      "capture_units",
      "year",
      "capture",
      "capture_notes"
    )
  
  
  colnames(fishstat_dat) <- names
  
  # fishstat_dat %>%
  #   group_by(year) %>%
  #   summarise(tc = sum(capture, na.rm = TRUE)) %>%
  #   ungroup() %>%
  #   filter(year == max(year))
  
  fao <- fishstat_dat %>%
    select(-bad) %>%
    select(
      id,
      country,
      common_name,
      scientific_name,
      fao_area,
      fao_area_code,
      order,
      family,
      everything()
    ) %>%
    filter(str_detect(production_source, "Capture")) %>%
    arrange(id) %>%
    group_by(id) %>%
    mutate(first_year = year[capture > 0 & !is.na(capture)][1]) %>%
    filter(year >= first_year,
           capture_units == "Tonnes") %>%
    ungroup() %>%
    group_by(id) %>%
    mutate(missing_catch = sum(is.na(capture))) %>%
    # filter(missing_catch == 0) %>%
    ungroup() %>%
    mutate(fao_area_code = as.numeric(fao_area_code)) %>%
    filter(!str_detect(country, "Totals"),
           isscaap_number < 60)
  
  fao_species <<- fao %>%
    select(scientific_name, common_name, isscaap_group, isscaap_number) %>%
    unique()
  
  fao_genus <-
    str_split(fao_species$scientific_name, ' ', simplify = TRUE)[, 1]
  
  fao_genus = fao_species %>%
    mutate(genus = fao_genus) %>%
    group_by(genus, isscaap_group) %>%
    count() %>%
    group_by(genus) %>%
    filter(n == max(n)) %>%
    select(-n) %>%
    ungroup()
  
  
  fao$fao_country_name <-
    countrycode::countrycode(fao$country, "country.name", "fao.name")
  
  fao <- fao %>%
    mutate(country = case_when(is.na(fao_country_name) ~ country, TRUE ~ fao_country_name))
  
  
  fao$continent <-
    countrycode::countrycode(fao$country, "country.name", "continent")
  
  fao_stock_lookup <<- fao %>%
    select(scientific_name,
           common_name,
           country,
           fao_area,
           fao_area_code) %>%
    unique()
  
  fao <<- fao %>%
    mutate(year = as.numeric(year))
  
  # load fmi data -----------------------------------------------------------
  
  
  ram_fmi_sheets <-
    readxl::excel_sheets(here("data", "RAM FMI stocks and taxonomic groupings.xlsx"))
  
  
  ram_fmi <-
    map(ram_fmi_sheets, ~ readxl::read_xlsx(
      here("data", "RAM FMI stocks and taxonomic groupings.xlsx"),
      sheet = .x
    )) %>%
    set_names(str_replace(tolower(ram_fmi_sheets), " ", "_"))
  
  ram_species <- ram_fmi$ram_species %>%
    janitor::clean_names() %>%
    rename(scientificname = scientificname_ram)
  
  ram_data <<- ram_data %>%
    left_join(ram_species, by = "scientificname") %>%
    ungroup()
  
  
  ram_fmi_linkages <-
    readxl::read_xlsx(
      here::here("data", "RAM-FMI linkages for DO 2019-09-16.xlsx"),
      sheet = "RAM-FMI linkages",
      skip = 1
    ) %>%
    janitor::clean_names() %>%
    select(-contains("primary")) %>%
    select(stockid,
           areaid,
           iso3_pref,
           scientificname,
           fao_scientific_name) %>%
    mutate(iso3_pref = stringr::str_trim(stringr::str_replace_all(iso3_pref, ' \\| ', "\\|")))
  
  
  ram_fmi_fao_regions_linkages <-
    readxl::read_xlsx(
      here::here("data", "RAM-FMI linkages for DO 2019-09-16.xlsx"),
      sheet = "csv",
      skip = 0
    ) %>%
    janitor::clean_names() %>%
    mutate(lookup_code = stringr::str_trim(stringr::str_replace_all(cspref, ' \\| ', "\\|"))) %>%
    select(-contains("_recent"), stockid, region, primary_country, faor)
  
  
  fmi <-
    readxl::read_xlsx(here::here("data", "FMI data extract by stock for Ray 2018-11-05.xlsx"),
                      sheet = "summary") %>%
    janitor::clean_names()
  
  
  fmi$fao_country_name <-
    countrycode::countrycode(fmi$country_rfmo, "country.name", "fao.name")
  
  fmi$region <-
    countrycode::countrycode(fmi$country_rfmo, "country.name", "region")
  
  # fmi <- fmi %>%
  #   mutate(country_rfmo = case_when(is.na(fao_country_name) ~ country_rfmo, TRUE ~ fao_country_name))
  #
  # fmi$country <- fmi$country_rfmo
  #
  # fmi$country_rfmo <-
  #   case_when(
  #     (fmi$region == "Southern Europe"  &
  #        (fmi$basin == "Med" |
  #           is.na(fmi$basin))) |
  #       (fmi$basin == "Med" | is.na(fmi$basin))  ~ "GFCM",
  #     fmi$region %in%  c("Northern Europe", "Western Europe") |
  #       (fmi$region == "Southern Europe" &
  #          !(fmi$basin == "Med" | is.na(fmi$basin))) ~ "ICES",
  #     TRUE ~ fmi$country_rfmo
  #   )
  
  if (lookup_fmi_names == TRUE) {
    fmi_names <- fmi %>%
      select(lookup_code, country_rfmo)
    
    temp <-
      fmi_names$lookup_code %>%
      str_split("\\|", simplify = TRUE) %>%
      as_data_frame() %>%
      set_names("code", "species", "region") %>%
      map_df(str_trim) %>%
      select(species) %>%
      unique() %>%
      mutate(sciname = map(species, ~ taxize::comm2sci(commnames = .x, db = "worms")))
    
    get_names <- function(name) {
      name <- name[[1]]
      if (length(name) == 0) {
        out <- NA
      } else {
        out <- name[1]
      }
      return(out)
    }
    
    
    fmi_scinames <- map_chr(temp$sciname, get_names)
    
    temp <- temp %>%
      mutate(scientificname = fmi_scinames) %>%
      select(-sciname)
    
    saveRDS(temp, file = here::here("data", "fmi_scinames.rds"))
    fmi_scinames <- temp
  } else {
    fmi_scinames <-
      readRDS(file = here::here("data", "fmi_scinames.rds"))
  }
  
  fmi <- fmi %>%
    left_join(fmi_scinames %>% rename(scientificname2 = scientificname),
              by = "species")
  
  fmi$lookup_code <-
    stringr::str_trim(stringr::str_replace_all(fmi$lookup_code, ' \\| ', "\\|"))
  
  fmi <- fmi %>%
    left_join(ram_fmi_linkages, by = c("lookup_code" = "iso3_pref")) %>%
    left_join(ram_fmi_fao_regions_linkages,
              by = c("stockid", "lookup_code")) %>%
    mutate(scientificname = ifelse(is.na(scientificname), scientificname2, scientificname)) %>%
    select(-scientificname2)
  
  
  fmi <-
    fmi %>%
    mutate(genus = map_chr(scientificname,  ~ str_split(.x, " ", simplify = TRUE)[, 1])) %>%
    left_join(fao_species, by = c("scientificname" = "scientific_name")) %>%
    left_join(fao_genus, by = "genus") %>%
    mutate(isscaap_group = ifelse(is.na(isscaap_group.x), isscaap_group.y, isscaap_group.x)) %>%
    select(-isscaap_group.x, -isscaap_group.y)
  
  ram_v_fmi <<- ram_data %>%
    group_by(stockid) %>%
    mutate(
      c_maxc = catch / max(catch, na.rm = TRUE),
      c_meanc = catch / mean(catch, na.rm = TRUE)
    ) %>%
    filter(year > (max(year) - 5)) %>%
    summarise(
      mean_bbmsy = mean(b_v_bmsy, na.rm = TRUE),
      mean_uumsy = mean(u_v_umsy, na.rm = TRUE),
      mean_f = mean(-log(1 - pmin(0.95,exploitation_rate)), na.rm = TRUE),
      c_div_max_c = mean(c_maxc)) %>%
    gather(metric, value,-stockid, -c_div_max_c) %>%
    ungroup() %>%
    left_join(fmi, by = "stockid") %>%
    filter(!is.na(lookup_code)) %>%
    select(-basin, -stock, -scientificname,-species,-contains("fao"),-contains(".x")) %>%
    mutate(log_value = log(value)) %>%
    unique() %>%
    na.omit() %>%
    mutate_at(c("research", "management", "enforcement", "socioeconomics"), ~ .x + 1e-6)
  
  
  # load sar data -----------------------------------------------------------
  
  recent_ram <- ram_data %>%
    group_by(stockid) %>%
    mutate(
      c_maxc = catch / max(catch, na.rm = TRUE),
      c_meanc = catch / mean(catch, na.rm = TRUE)
    ) %>%
    filter(year > (max(year[!is.na(u_v_umsy)]) - 5)) %>%
    summarise(
      mean_bbmsy = mean(b_v_bmsy, na.rm = TRUE),
      mean_uumsy = mean(u_v_umsy, na.rm = TRUE),
      mean_f = mean(-log(1 - pmin(
        0.95, exploitation_rate
      )), na.rm = TRUE),
      c_div_max_c = mean(c_maxc),
      c_div_mean_c = mean(c_meanc)
    ) %>%
    na.omit()
  
  
  
  sar_coverage <-
    readr::read_csv(here::here("data", "OverlapTable2.csv")) %>%
    janitor::clean_names() %>%
    group_by(stockid) %>%
    summarise(
      mean_tbp_in_stock = mean(tbp_in_stock),
      mean_stock_in_tbp = mean(stock_in_tbp)
    ) %>%
    ungroup()
  
  
  sar_to_ram <-
    readr::read_csv(here::here("data", "RamStocksWithID2.csv")) %>%
    janitor::clean_names() %>%
    map_df(stringr::str_trim) %>% # sigh, white spaces in the numerics work on mac but not linux
    modify_at(4:7, as.numeric) %>% # fun fun
    mutate(log_f = log(fstatus + 1e-3)) %>%
    mutate(genus = map_chr(latin_binomial,  ~ str_split(.x, " ", simplify = TRUE)[, 1])) %>%
    left_join(fao_species, by = c("latin_binomial" = "scientific_name")) %>%
    left_join(fao_genus, by = "genus") %>%
    mutate(isscaap_group = ifelse(is.na(isscaap_group.x), isscaap_group.y, isscaap_group.x)) %>%
    select(-isscaap_group.x,-isscaap_group.y) %>%
    left_join(sar_coverage, by = "stockid") %>%
    left_join(recent_ram, by = "stockid")
  
  sar = sar_to_ram
  
  ram_v_sar <<- sar_to_ram %>%
    gather(
      metric,
      value,
      contains("mean_"),
      -c_div_mean_c,
      -mean_stock_in_tbp,
      -mean_tbp_in_stock
    ) %>%
    mutate(log_value = log(value + 1e-3)) %>%
    mutate(sar_2 = sar ^ 2) %>%
    select(
      stockid,
      sar,
      sar_2,
      isscaap_group,
      metric,
      value,
      log_value,
      c_div_mean_c,
      c_div_max_c,
      mean_stock_in_tbp
    ) %>%
    filter(mean_stock_in_tbp > 25 | is.na(mean_stock_in_tbp)) %>%
    select(-mean_stock_in_tbp) %>%
    na.omit()
  
  
  # load effort data --------------------------------------------------------
  
  
  # load effort data
  
  effort_data <<-
    readr::read_csv(here::here("data", "rousseau-2019", "Data_Effort_CPUE_forRH.csv")) %>%
    janitor::clean_names() %>%
    select(year, region, sector, contains("_effort_")) %>%
    tidyr::pivot_longer(contains("_effort_"),
                        names_to = c("effort_type",'effort_units'),
                        names_pattern = "(.*)_effort_(.*)",
                        values_to = "effort")
  
  effort_region_to_country <<-
    readxl::read_xlsx(here::here("data", "rousseau-2019", "pnas.1820344116.sd01.xlsx")) %>%
    janitor::clean_names() %>%
    select(country, sociocult_region_for_results) %>%
    unique() %>%
    rename(region = sociocult_region_for_results) %>%
    na.omit() %>%
    mutate(country =countrycode::countrycode(country, "country.name", "fao.name"))
  
  # ugh. manual refactoring of region codes
  
  refactor_regions <-
    forcats::fct_recode(
      effort_region_to_country$region,
      "Sub-Saharan Africa" = "Sub Saharan",
      "Latin America" = "Latin Am.",
      "North America" = "N America",
      "Northeast Asia" = "NE Asia/Asia",
      "Southeast Asia"  = "SE Asia/Asia",
      "Southeast Asia"  = "SE Asia/ Asia",
      "Sub-Saharan Africa" = "Sub Saharan",
      "Indian Peninsula" = "Indian Pen. / Asia"
    )
  
  effort_region_to_country$region <- as.character(refactor_regions)
  
  effort_region_to_country <<- effort_region_to_country %>%
    unique()
  
  effort_region_to_country$country <- countrycode::countrycode( effort_region_to_country$country, "country.name", "fao.name")
  
  
  interpfoo <- function(data) {
    out <-
      approx(data$x, data$mean, xout = ceiling(min(data$x)):floor(max(data$x)))
    
    out <- data_frame(year = out$x, effort = out$y)
  }
  
  # east india
  
  ei_cpue <-
    read_csv(here::here("data", "east-india-cpue.csv")) %>%
    janitor::clean_names()
  
  ei_cpue <<- interpfoo(ei_cpue) %>%
    rename(cpue = effort)
  
  ei_capture <<-
    readxl::read_xlsx(here::here("data", "India East Only.xlsx")) %>%
    janitor::clean_names() %>%
    gather(year, capture, x1950:x2015) %>%
    mutate(year = as.numeric(str_replace_all(year, "\\D", "")))
  
}