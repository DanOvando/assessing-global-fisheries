run_case_studies <- TRUE

exs <- ram_data %>%
  filter(
    stockid %in% unique(sraplus::fmi$stockid) &
      stockid %in% unique(sraplus::sar$stockid)
  ) %>%
  mutate(fao_area_code = as.numeric(primary_FAOarea),
         country = as.character(primary_country)) %>%
  mutate(country =
           countrycode::countrycode(country, "country.name", "un.name.en"))

comp_stocks <- exs %>%
  select(scientificname) %>%
  unique() %>%
  mutate(fishbase_vuln = map_dbl(
    scientific_name,
    ~ rfishbase::species(.x, fields = "Vulnerability")$Vulnerability
  )) %>%
  mutate(
    resilience = dplyr::case_when(
      fishbase_vuln < 33 ~ "Low",
      fishbase_vuln >= 66 ~ "High",
      TRUE ~ "Medium"
    )
  ) %>%
  select(scientificname, resilience)


exs %>%
  ggplot(aes(year, b_v_bmsy, color = stockid)) +
  geom_line(show.legend = FALSE)

total_nominal_effort <- rous_data %>%
  left_join(effort_region_to_country, by = "country") %>%
  rename(fao_area_code = area) %>%
  group_by(fao_area_code, year, country) %>%
  summarise(total_effort = sum(effort_cell_reported_nom, na.rm = TRUE)) %>%
  ungroup()

matched_fmi <- sraplus::fmi %>%
  filter(stockid %in% unique(exs$stockid)) %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup() %>%
  mutate(tmp = map(data, ~ map_df(.x[,c("research", "management", "enforcement","socioeconomics")], mean))) %>%
  select(-data) %>%
  unnest(cols = tmp) %>%
  select(stockid, research, management, enforcement, socioeconomics)


matched_sar <- sraplus::sar %>%
  filter(stockid %in% unique(exs$stockid)) %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup() %>%
  mutate(tmp = map(data, ~ purrrlyr::dmap_at(.x, c("sar"), mean))) %>%
  select(-data) %>%
  unnest(cols = tmp) %>%
  select(stockid, sar)


exs <- exs %>%
  left_join(total_nominal_effort, by = c("fao_area_code", "year", "country")) %>%
  left_join(matched_fmi, by = "stockid") %>%
  left_join(matched_sar, by = "stockid") %>%
  left_join(comp_stocks, by = "scientificname")


if (run_case_studies){

exs <- exs %>%
  group_by(stockid) %>%
  nest() %>%
  ungroup() %>% 
  mutate(fit = map(data,
                   safely(fit_case_studies)))

write_rds(exs, path = file.path(results_path, "case-study-fits.rds"))

} else {
  
  exs <- read_rds(path = file.path(results_path, "case-study-fits.rds"))
  
}

fit_worked <-
  exs$fit %>% map("error") %>% map_lgl(is.null)


truth <- exs %>% 
  select(stockid, data) %>% 
  unnest(cols = data) %>% 
  select(stockid, year, catch, b_v_bmsy, u_v_umsy) %>% 
  group_by(stockid) %>% 
  filter(year == max(year)) %>% 
  ungroup()

tmp <- exs %>%
  mutate(lifetime_catch = map_dbl(data,  ~ sum(.x$catch))) %>%
  select(stockid, fit, lifetime_catch) %>%
  filter(fit_worked) %>%
  mutate(fit = map(fit, "result")) %>%
  unnest(cols = fit) %>%
  group_by(stockid) %>%
  filter(year == max(year)) %>%
  filter(variable %in% c("b_div_bmsy", "u_div_umsy")) %>%
  ungroup() %>%
  select(stockid,
         year,
         variable,
         mean,
         data,
         lifetime_catch) %>% 
  pivot_wider(names_from = variable, values_from = mean)




ex_kobe_plot <- tmp %>%
  ggplot() +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_vline(aes(xintercept = 1), linetype = 2) +
  geom_hline(aes(yintercept = 0), linetype = 1) +
  # geom_density2d(alpha = 0.5) +
  geom_hex(data = truth,aes(pmin(4, b_v_bmsy), pmin(4, u_v_umsy)), binwidth = c(0.5,.5), show.legend = FALSE, alpha = 0.9) +
  geom_point(aes(pmin(4, b_div_bmsy), pmin(4, u_div_umsy), color = data,size = lifetime_catch),
             alpha = 0.5,
             show.legend = FALSE) +
  scale_size(range = c(2,7)) +
  scale_x_continuous(limits = c(-.1, 4), name = "B/Bmsy") +
  scale_y_continuous(limits = c(-.1, 4), name = "U/Umsy") +
  facet_wrap(~ data) +
  scale_fill_continuous(low = "gainsboro", high = "black") +
  ggsci::scale_color_d3(name = "Region") +
  guides(size = FALSE)

sraplus_v_truth <- tmp %>% 
  left_join(truth, by = c("stockid","year")) %>% 
  filter(data == "cpue") 

sraplus_v_truth %>% 
  ggplot(aes(catch, b_div_bmsy)) + 
  geom_point(size = 2)
  
