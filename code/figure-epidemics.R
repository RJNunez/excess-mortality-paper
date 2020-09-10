### -- -------------------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- -------------------------------------- ------------------------------------------------------------------
# -- Set up
source("code/pr-init.R")

# -- Age groups
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)

# -- Number of knots per year
nknots <- 6

# -- To be used as parameters in model fitting below
before <- days(365)
after  <- days(365)

# -- Creating breaks for age groups and collapsing data
the_breaks <- c(0, 5, 20, 40, 60, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks)

# -- Event dates
event_dates <- ymd("2014-08-01")
events      <- "Chikungunya"

# -- Fitting model only to the groups of interest
res <- map_df(levels(all_counts$agegroup), function(x){
  
  print(x)
  
  tmp_counts <- filter(all_counts, agegroup == x)
  tmp_fit    <- suppressMessages(excess_model(counts         = tmp_counts, 
                                              start          = event_dates - before,
                                              end            = event_dates + after,
                                              exclude        = exclude_dates,
                                              control.dates  = control_dates,
                                              weekday.effect = TRUE,
                                              aic            = FALSE, 
                                              order.max      = 7,
                                              verbose        = FALSE,
                                              knots.per.year = nknots,
                                              model          = "correlated",
                                              discontinuity  = FALSE))
  
  
  tibble(date       = tmp_fit$date, 
         fitted     = tmp_fit$fitted, 
         observed   = tmp_fit$observed, 
         expected   = tmp_fit$expected, 
         se         = tmp_fit$se, 
         event      = events, 
         event_date = event_dates,
         agegroup   = x)
})

# -- Figure
res %>%
  group_by(date) %>%
  summarize(fitted     = sum(fitted * expected) / sum(expected),
            se         = sqrt(sum(se^2 * expected^2)) / sum(expected),
            event_date = event_date[1]) %>%
  mutate(day = as.numeric(date - event_date),
         lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  filter(day >= -120, day <= 360) %>%
  ggplot(aes(date, fitted)) +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.20) +
  geom_line() +
  xlab("Date") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(-0.10, 0.25, by = 0.05)) +
  theme(text = element_text(size = 10))

ggsave(filename = "figs/figure-chikungunya-effect.pdf",
       width    = 4,
       height   = 3, 
       dpi      = 300)
### -- -------------------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- -------------------------------------- ------------------------------------------------------------------

### -- ------------------------------------- ------------------------------------------------------------------
### -- Figure: Excess deaths for Chikungunya ------------------------------------------------------------------
### -- ------------------------------------- ------------------------------------------------------------------
# -- Set up
source("code/pr-init.R")

# -- Set up to be used below
ndays  <- 365
knots  <- c(6)
disc   <- c(FALSE)
before <- c(365) 
after  <- c(365)
interval_start <- c(Chikungunya = make_date(2014,08,01))

# -- Outer loop that goes through events in PR
chikungunya <- map_df(seq_along(interval_start), function(i)
{
  cat("\nEvent:", names(interval_start)[i], "\n")
  
  # -- Now fit model to compute cumulative excess deaths
  message("\nFitting model to estimate cumulative excess deaths")
  tmp <- map_df(unique(all_counts$agegroup), function(x){
    cat(".")
    f <- suppressMessages(all_counts %>% 
                            filter(agegroup == x) %>%
                            excess_model(event          = interval_start[i],
                                         start          = interval_start[i] - before[i],
                                         end            = interval_start[i] + after[i], 
                                         exclude        = exclude_dates,
                                         control.dates  = control_dates,
                                         knots.per.year = knots[i],
                                         aic            = FALSE, 
                                         order.max      = 7,
                                         weekday.effect = TRUE,
                                         model          = "correlated",
                                         discontinuity  = disc[i], 
                                         verbose = FALSE))
    
    ndays <- 365*2
    excess_cumulative(f, 
                      start = interval_start[i], 
                      end   = ymd(interval_start[i]) + ndays) %>%
      mutate(agegroup = x, event_day = interval_start[i], event = names(interval_start)[i]) %>%
      as_tibble()
    
  })
  tmp %>% 
    group_by(date) %>% 
    summarize(fitted = sum(fitted),
              observed = sum(observed),
              sd = sqrt(sum(sd^2)),
              se = sqrt(sum(se^2)),
              event_day = event_day[1], 
              event = event[1]) %>%
    ungroup()
})

# -- Excess deaths one year after
chikungunya %>%
  mutate(lwr = fitted - 1.96*se,
         upr = fitted + 1.96*se) %>%
  filter(date == "2015-08-01") %>%
  select(fitted, lwr, upr)

# -- Figure
chikungunya %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  ggplot(aes(date, fitted)) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.20) +
  geom_line() +
  ylab("Cumulative excess deaths") +
  xlab("Date") +
  scale_x_date(date_labels = "%b %Y", date_breaks = "3 months") +
  scale_y_continuous(labels = scales::comma) +
  theme(text = element_text(size = 10))

ggsave(filename = "figs/figure-chikungunya-excess-deaths.pdf",
       width    = 4,
       height   = 3, 
       dpi      = 300)
### -- ------------------------------------- ------------------------------------------------------------------
### -- Figure: Excess deaths for Chikungunya ------------------------------------------------------------------
### -- ------------------------------------- ------------------------------------------------------------------

### -- ----------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimate for USA ------------------------------------------------------------------
### -- ----------------------------- ------------------------------------------------------------------
# -- Loading state mortality data 
source("code/pr-init.R")
data("cdc_state_counts")

# -- Expand state abbrevaition objects 
state.name.2 <- c(state.name, "New York City", "Puerto Rico", "District of Columbia")
state.abb.2  <- c(state.abb, "NYC", "PR", "DC")

# -- Importing covd-19 reported deaths data 
covid_nyc    <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
covid_nyc    <- filter(covid_nyc, county =="New York City")
covid_states <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  mutate(date = ymd(date)) %>%
  filter(!is.na(state)) %>%
  arrange(state)

# -- Subsetting new york data
ny <- filter(covid_states, state == "New York")

# -- Covid 19 data for the rest of new york
ny <- left_join(ny, covid_nyc, by = "date") %>%
  mutate(death = deaths.x - deaths.y, 
         state = "Rest of New York") %>%
  select(date, state, death)

# -- Covid 19 data for states
covid_states <- filter(covid_states, state!="New York") %>%
  rename(death = deaths) %>%
  select(date, state, death)

# -- Covid 19 for NYC
covid_nyc <- covid_nyc %>%
  mutate(state = "New York City", death = deaths)%>%
  select(date, state, death)

# -- Putting data together
covid_states <- bind_rows(covid_states, ny, covid_nyc) %>%
  filter(!is.na(death))
rm(covid_nyc, ny)

# -- Denoting periods of interest
flu_season     <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates  <- c(flu_season, seq(make_date(2020, 1, 1), today(), by = "day"))
max_weighted   <- last(cdc_state_counts$date) - weeks(1)
max_unweighted <- max_weighted - weeks(6)

# -- Remove last dates
weight   <- cdc_state_counts %>% filter(date <= max_weighted)
unweight <- cdc_state_counts %>% filter(date <= max_unweighted)
states   <- unique(weight$state)
states   <- setdiff(states, c("Connecticut", "North Carolina", "Puerto Rico"))

# -- Fitting the model to each state
nknots <- 16
fit <- map_df(states, function(x){
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  print(x)
  w <- weight %>% 
    filter(state == x) %>%
    na.omit() %>%
    excess_model(exclude        = exclude_dates,
                 start          = min(weight$date),
                 end            = max_weighted,
                 aic            = FALSE, 
                 order.max      = 7,
                 knots.per.year = nknots,
                 weekday.effect = FALSE,
                 verbose        = FALSE)
  
  w <- with(w,
            tibble(date     = date, 
                   expected = expected, 
                   observed = observed,
                   fitted   = fitted, 
                   se       = se,
                   sd       = sd)) %>%
    mutate(state = x, 
           type  = "weighted")
  
  
  u <- unweight %>% 
    filter(state == x) %>%
    select(-outcome) %>%
    rename(outcome = outcome_unweighted) %>%
    na.omit() %>%
    excess_model(exclude        = exclude_dates,
                 start          = min(unweight$date),
                 end            = max_unweighted,
                 aic            = FALSE, 
                 order.max      = 7,
                 knots.per.year = nknots,
                 weekday.effect = FALSE,
                 verbose        = FALSE)
  
  u <- with(u,
            tibble(date     = date, 
                   expected = expected, 
                   observed = observed,
                   fitted   = fitted, 
                   se       = se,
                   sd       = sd)) %>%
    mutate(state = x, 
           type  = "unweighted")
  
  bind_rows(w, u)
})

# -- To be used in the viz
tmp <- fit %>%
  filter(!state %in% c("Puerto Rico", "North Carolina", "Connecticut")) %>%
  group_by(date, type) %>% 
  summarize(fitted   = sum(expected * fitted) / sum(expected), 
            se       = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd       = sqrt(sum(expected^2 * sd^2)) / sum(expected),
            expected = sum(expected),
            observed = sum(observed)) %>%
  ungroup() %>%
  mutate(type = ifelse(type == "weighted", "CDC weighted", "CDC unweighted"))

# -- Figure
ggplot() +
  geom_hline(yintercept = 0, lty=2, color="gray") +
  geom_ribbon(aes(x = date, ymin = fitted - 1.96*se, ymax = fitted + 1.96*se, fill=type),
              alpha = 0.2, color=NA, show.legend = F, data = tmp) +
  geom_line(aes(date, fitted, color=type), show.legend = F, data = tmp) +
  geom_point(aes(date, observed / expected - 1, color=type), size=1, alpha=0.20, show.legend = F, data = tmp) + 
  xlab("Date") +
  ylab("Percent increase from expected mortality") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_labels = "%b %Y") +
  theme(text = element_text(size = 10))

ggsave(filename = "figs/figure-usa-covid-effect.pdf",
       width    = 4,
       height   = 3, 
       dpi      = 300)
### -- ----------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimate for USA ------------------------------------------------------------------
### -- ----------------------------- ------------------------------------------------------------------

### -- --------------------------- ------------------------------------------------------------------
### -- Figure: Excess deaths in US ------------------------------------------------------------------
### -- --------------------------- ------------------------------------------------------------------
# -- Loading state mortality data 
source("code/pr-init.R")
data("cdc_state_counts")

# -- Expand state abbrevaition objects 
state.name.2 <- c(state.name, "New York City", "Puerto Rico", "District of Columbia")
state.abb.2  <- c(state.abb, "NYC", "PR", "DC")

# -- Importing covd-19 reported deaths data 
covid_nyc    <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv")
covid_nyc    <- filter(covid_nyc, county =="New York City")
covid_states <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
  mutate(abb = state.abb.2[match(state, state.name.2)]) %>%
  mutate(date = ymd(date)) %>%
  filter(!is.na(state)) %>%
  arrange(state)

# -- Subsetting new york data
ny <- filter(covid_states, state == "New York")

# -- Covid 19 data for the rest of new york
ny <- left_join(ny, covid_nyc, by = "date") %>%
  mutate(death = deaths.x - deaths.y, 
         state = "Rest of New York") %>%
  select(date, state, death)

# -- Covid 19 data for states
covid_states <- filter(covid_states, state!="New York") %>%
  rename(death = deaths) %>%
  select(date, state, death)

# -- Covid 19 for NYC
covid_nyc <- covid_nyc %>%
  mutate(state = "New York City", death = deaths)%>%
  select(date, state, death)

# -- Putting data together
covid_states <- bind_rows(covid_states, ny, covid_nyc) %>%
  filter(!is.na(death))
rm(covid_nyc, ny)

# -- Denoting periods of interest
flu_season     <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates  <- c(flu_season, seq(make_date(2020, 1, 1), max(cdc_state_counts$date, na.rm = TRUE), by = "day"))
max_weighted   <- last(cdc_state_counts$date) - weeks(1)
max_unweighted <- max_weighted - weeks(6)

# -- Remove last dates
weight   <- cdc_state_counts %>% filter(date <= max_weighted)
unweight <- cdc_state_counts %>% filter(date <= max_unweighted)
states   <- unique(weight$state)
states   <- setdiff(states, c("Connecticut", "North Carolina", "Puerto Rico"))

# -- Fitting model to each date. Using 6 knots as before
nknots <- 16
fit <- map_df(states, function(x){
  if(x == "Puerto Rico"){
    exclude_dates <- unique(sort(c(exclude_dates, seq(make_date(2017, 9, 20), make_date(2018, 3, 31), by = "day"))))
  }
  print(x)
  w <- weight %>% 
    filter(state == x) %>%
    na.omit() %>%
    excess_model(exclude        = exclude_dates,
                 start          = min(weight$date),
                 end            = max_weighted,
                 aic            = FALSE, 
                 order.max      = 7,
                 knots.per.year = nknots,
                 weekday.effect = FALSE,
                 verbose        = FALSE)
  
  w <- excess_cumulative(w, start = make_date(2020, 03, 01), end = max_weighted) %>%
    mutate(state = x, type = "weighted")
  
  u <- unweight %>% 
    filter(state == x) %>%
    select(-outcome) %>%
    rename(outcome = outcome_unweighted) %>%
    na.omit() %>%
    excess_model(exclude        = exclude_dates,
                 start          = min(unweight$date),
                 end            = max_unweighted,
                 aic            = FALSE, 
                 knots.per.year = nknots,
                 order.max      = 7,
                 weekday.effect = FALSE,
                 verbose        = FALSE)
  
  u <- excess_cumulative(u, start = make_date(2020, 03, 01), end = max_unweighted) %>%
    mutate(state = x, type = "unweighted")
  
  bind_rows(w, u)
}) %>% as_tibble()

# -- Cumulative deaths by Covid-19 in USA
covid_us <- covid_states %>% 
  filter(date >= min(fit$date)) %>%
  group_by(date) %>%
  summarize(covid = sum(death))

# -- Cumulative deaths in USA
fit <- fit %>%
  group_by(date, type) %>%
  summarize(observed = sum(observed),
            sd       = sqrt(sum(sd^2)),
            fitted   = sum(fitted),
            se       = sqrt(sum(se^2))) %>%
  ungroup() %>%
  mutate(lwr = fitted - 1.96 * se, 
         upr = fitted + 1.96 * se) %>%
  mutate(type = ifelse(type == "weighted", "CDC weighted", "CDC unweighted"))

# -- For floating legend
lab <- tibble(posx = c(make_date(2020,05,01), make_date(2020,05,01)),
              posy = c(33000, 20000),
              type = c("CDC weighted", "CDC unweighted"))

# -- Figure
ggplot() +
  geom_ribbon(aes(date, ymin=lwr, ymax=upr, fill=type), alpha=0.20, show.legend = F, data = fit) +
  geom_line(aes(date, fitted, color=type), show.legend = F, data = fit) +
  geom_point(aes(date, observed, color=type), size=1, show.legend = F, alpha=0.20, data = fit) +
  geom_line(aes(date, covid), data = covid_us) +
  ylab("Cumulative excess deaths") +
  xlab("Date") +
  geom_text(aes(posx, posy, color=type, label=type), hjust=0, size=3, show.legend = F, data = lab) +
  geom_text(aes(make_date(2020,05,01), 7000, label="Reported Covid-19 deaths"), hjust=0, size=3) +
  scale_y_continuous(labels = scales::comma,
                     breaks = seq(0, 200000, by=50000)) +
  scale_x_date(date_labels = "%b %Y", date_breaks = "2 months") +
  theme(text = element_text(size = 10))

ggsave(filename = "figs/figure-usa-covid-excess-deaths.pdf",
       width    = 4,
       height   = 3, 
       dpi      = 300)
### -- --------------------------- ------------------------------------------------------------------
### -- Figure: Excess deaths in US ------------------------------------------------------------------
### -- --------------------------- ------------------------------------------------------------------
