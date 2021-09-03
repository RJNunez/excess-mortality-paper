# -- Set up
source("code/pr-init.R")

# -- Years to used in the cross-validation
years <- 1999:2013

# -- Wrangling data
dat <- puerto_rico_counts %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Expected values using all data
all <- compute_expected(counts         = dat,
                        weekday.effect = FALSE,
                        exclude        = exclude_dates) %>%
  as_tibble() %>%
  filter(lubridate::year(date) %in% years) %>%
  select(date, outcome, expected)

# -- Expected values taking out each year
res <- map_df(1:length(years), function(i){
  
  print(years[i])
  
  # -- Excluding the year
  tmp_exclude <- seq(make_date(years[i], 01, 01), make_date(years[i], 12, 31), by = "days")
  
  # -- Fitting mean model
  tmp_mu <- compute_expected(counts  = dat,
                             weekday.effect = FALSE,
                             exclude = c(exclude_dates, tmp_exclude)) %>% 
    as_tibble() %>%
    filter(lubridate::year(date) == years[i]) %>%
    select(date, expected) %>%
    rename(c_expected = expected)
})

# -- Weekly average of daily death counts
smooth <- dat %>%
  filter(lubridate::year(date) %in% years) %>%
  mutate(year = lubridate::year(date),
         week = lubridate::week(date)) %>%
  group_by(year, week) %>%
  summarize(outcome_smooth = mean(outcome),
            date           = mean(date)) %>%
  ungroup() %>%
  select(date, outcome_smooth)

# -- To be used for the viz below
viz_dat <- res %>%
  left_join(all, by = "date") %>%
  pivot_longer(cols = c(c_expected, expected)) %>%
  mutate(year = lubridate::year(date),
         day  = lubridate::yday(date),
         name = factor(name, levels = c("expected", "c_expected"))) %>%
  left_join(smooth, by = "date")

# -- Supplementary figure: Cross validation study
ggplot() +
  coord_cartesian(ylim = c(70, 100)) +
  facet_wrap(~year) +
  labs(x = "Day of the year",
       y = "Daily counts",
       color    = "",
       linetype = "") +
  geom_point(aes(day, outcome_smooth), 
             alpha = 0.20,
             data  = viz_dat) +
  geom_line(aes(day, value, color = name, linetype = name),
            # alpha = 0.80,
            size  = 0.80,
            data  = viz_dat) +
  scale_color_manual(values = c(my_palette[["blue"]], my_palette[["orange"]]),
                     labels = c("expected" = "Excluding each year ", "c_expected" = "Including all years ")) +
  scale_linetype_manual(values = c(1, 2),
                        labels = c("expected" = "Excluding each year ", "c_expected" = "Including all years ")) +
  theme(legend.position  = c(0.88, 0.10),
        legend.direction = "vertical",
        legend.title     = element_blank())
