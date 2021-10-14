### -- -------------------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimates for Chikungunya ------------------------------------------------------------------
### -- -------------------------------------- ------------------------------------------------------------------
# -- Set up
source("code/pr-init.R")

# -- Number of knots per year
nknots <- 6

# -- To be used as parameters in model fitting below
before <- days(365)
after  <- days(365)

# -- Event dates
event_dates <- ymd("2014-08-01")
events      <- "Chikungunya"

# -- Collapsing sex counts
all_counts <- all_counts %>%
  group_by(date, agegroup) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Fitting model only to the groups of interest
res <- map_df(levels(all_counts$agegroup), function(x){
  
  # -- Print current group
  print(x)
  
  # -- Fitting model to age group
  tmp_counts <- filter(all_counts, agegroup == x)
  tmp_fit    <- suppressMessages(excess_model(counts         = tmp_counts, 
                                              start          = event_dates - before,
                                              end            = event_dates + after,
                                              exclude        = exclude_dates,
                                              control.dates  = control_dates,
                                              weekday.effect = FALSE,#TRUE,
                                              aic            = FALSE, 
                                              order.max      = 7,
                                              verbose        = FALSE,
                                              knots.per.year = nknots,
                                              model          = "correlated",
                                              discontinuity  = FALSE))
  
  # -- Putting everyting together
  tibble(date            = tmp_fit$date, 
         fitted          = tmp_fit$fitted, 
         observed        = tmp_fit$observed, 
         expected        = tmp_fit$expected, 
         log_expected_se = tmp_fit$log_expected_se,
         population      = tmp_fit$population,
         se              = tmp_fit$se, 
         sd              = tmp_fit$sd, 
         event           = events, 
         event_date      = event_dates,
         agegroup        = x)
})

# -- Figure
fig3a <- res %>%
  mutate(expected_se = expected * log_expected_se) %>%
  group_by(date) %>%
  mutate(sum_age_expected       = sum(expected),
         sum_age_gamma_expected = sum(fitted * expected)) %>%
  ungroup() %>%
  mutate(gamma_derivative = expected / sum_age_expected,
         mu_derivative    = (fitted / sum_age_expected) - (sum_age_gamma_expected / (sum_age_expected^2))) %>%
  mutate(tmp_se = se^2 * gamma_derivative^2 + expected_se^2 * mu_derivative^2) %>%
  group_by(date) %>%
  summarize(fitted      = sum(fitted * expected) / sum(expected),
            observed    = sum(observed),
            expected    = sum(expected),
            population  = sum(population),
            expected_se = sqrt(sum(expected_se^2)),
            se          = sqrt(sum(tmp_se)),
            event_date  = event_date[[1]]) %>%
  ungroup() %>%
  mutate(day = as.numeric(date - event_date),
         lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  filter(day >= -120, day <= 360) %>%
  ggplot(aes(date, fitted)) +
  geom_hline(yintercept = 0, 
             lty        = 2, 
             color      = "gray") +
  geom_ribbon(aes(ymin = lwr, 
                  ymax = upr), 
              alpha = 0.20,
              fill  = my_palette[["black"]]) +
  geom_line(size  = 1,
            color = "white") +
  geom_line(color = my_palette[["black"]]) +
  xlab("Date") +
  ylab("Percent increase from expected mortality") +
  scale_x_date(date_labels = "%b", date_breaks = "3 months") +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(-0.10, 0.25, by = 0.05)) +
  theme(text = element_text(size=10))

# -- Saving figure
ggsave(filename = "figures/figure-3a.pdf",
       plot     = fig3a,
       width    = 4,
       height   = 3)
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

# -- Collapsing by age groups
all_counts <- all_counts %>%
  group_by(date, agegroup) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Outer loop that goes through events in PR
chikungunya <- map_df(seq_along(interval_start), function(i)
{
  cat("\nEvent:", names(interval_start)[i], "\n")
  
  # -- Now fit model to compute cumulative excess deaths
  message("\nFitting model to estimate cumulative excess deaths")
  
  # -- Computing excess deaths for all age groups
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
                                         weekday.effect = FALSE,#TRUE,
                                         model          = "correlated",
                                         discontinuity  = disc[i], 
                                         verbose = FALSE))
    
    # -- Computing excess deaths
    excess_cumulative(f,
                      start = interval_start[i], 
                      end   = ymd(interval_start[i]) + ndays*2) %>%
      mutate(agegroup  = x, 
             event_day = interval_start[i], 
             event     = names(interval_start)[i]) %>%
      as_tibble()
  })
  
  # -- Summing excess deaths across all age groups to get total excess deaths
  tmp %>% 
    group_by(date) %>% 
    summarize(fitted    = sum(fitted),
              observed  = sum(observed),
              sd        = sqrt(sum(sd^2)),
              se        = sqrt(sum(se^2)),
              event_day = event_day[1], 
              event     = event[1]) %>%
    ungroup()
})

# -- Figure
fig3b <- chikungunya %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  ggplot(aes(date, fitted)) +
  ylab("Cumulative excess deaths") +
  xlab("Date") +
  scale_x_date(date_labels = "%b", date_breaks = "3 months") +
  scale_y_continuous(labels = scales::comma) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
              alpha = 0.20,
              fill  = my_palette[["black"]]) + 
  geom_line(size  = 1, 
            color = "white") +
  geom_line(color = my_palette[["black"]]) +
  theme(text = element_text(size=10))

# -- Saving figure
ggsave(filename = "figures/figure-3b.pdf",
       plot     = fig3b,
       width    = 4,
       height   = 3)

### -- ------------------------------------- ------------------------------------------------------------------
### -- Figure: Excess deaths for Chikungunya ------------------------------------------------------------------
### -- ------------------------------------- ------------------------------------------------------------------

### -- ----------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimate for USA ------------------------------------------------------------------
### -- ----------------------------- ------------------------------------------------------------------
# -- Loading state mortality data 
source("code/pr-init.R")
data("cdc_state_counts")

# -- Denoting periods of interest
flu_season     <- seq(make_date(2017, 12, 16), make_date(2018, 1, 16), by = "day")
exclude_dates  <- c(flu_season, seq(make_date(2020, 1, 1), today(), by = "day"))

# -- Wrangling the dataset. Using unweighted counts before the end of Jan 2021
unweighted_counts <- cdc_state_counts %>%
  filter(date <= "2021-01-31") %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome_unweighted),
            population = sum(population)) %>%
  ungroup()

# -- Using 16 knots to estimate the event effects
nknots <- 16

# -- Fitting excessmort model
fit <- excess_model(counts         = unweighted_counts,
                    exclude        = exclude_dates,
                    start          = ymd("2017-01-01"),
                    end            = ymd("2021-01-31"),
                    aic            = FALSE, 
                    order.max      = 7,
                    knots.per.year = nknots,
                    weekday.effect = FALSE,
                    verbose        = FALSE)

# -- Vizs: Changes in mortality in the USA
fig3c <- with(fit, tibble(date, fitted, se)) %>%
  mutate(lwr = fitted - 1.96*se, 
         upr = fitted + 1.96*se) %>%
  ggplot(aes(date, fitted)) +
  geom_hline(yintercept = 0, 
             lty        = 2, 
             color      = "gray") +
  geom_ribbon(aes(x    = date, 
                  ymin = lwr, 
                  ymax = upr),
              alpha       = 0.20, 
              color       = NA,
              fill        = my_palette[["black"]]) +
  geom_line(size  = 1,
            color = "white") +
  geom_line(color = my_palette[["black"]]) +
  xlab("Date") +
  ylab("Percent change from expected mortality") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_labels = "%b %Y") +
  theme(text = element_text(size = 10))

# -- Saving figure
ggsave(filename = "figures/figure-3c.pdf",
       plot     = fig3c,
       width    = 4,
       height   = 3)

# -- Getting excess deaths for USA
ed_usa <- excess_cumulative(fit, start = make_date(2020, 03, 01), end = make_date(2021, 01, 31)) %>%
  as_tibble() %>%
  mutate(lwr = fitted - 1.96*se, 
         upr = fitted + 1.96*se)

# -- Getting reported Covid-19 deaths for USA
covid_us <- read_csv("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-states.csv") %>%
  as_tibble() %>%
  filter(!state %in% c("Puerto Rico", "Virgin Islands", "Guam", "West Virginia", "Northern Mariana Islands"),
         date <= "2021-01-31") %>%
  group_by(date) %>%
  summarize(death = sum(deaths)) %>%
  ungroup()

# -- For floating legend
lab <- tibble(posx = c(make_date(2020,05,01), make_date(2020,05,01)),
              posy = c(33000, 20000),
              type = c("CDC weighted", "CDC unweighted"))
# -- Figure
fig3d <- ggplot() +
  geom_line(aes(date, death),
            color = "white",
            size  = 1,
            data  = filter(covid_us, date >= "2020-03-07")) +
  geom_line(aes(date, death),
            lty   = 2,
            color = "gray",
            data  = filter(covid_us, date >= "2020-03-07")) +
  geom_ribbon(aes(x    = date, 
                  ymin = lwr, 
                  ymax = upr), 
              alpha = 0.20,
              fill  = my_palette[["black"]],
              data  = ed_usa) +
  geom_line(aes(x = date, y = fitted), 
            size  = 1,
            color = "white",
            data  = ed_usa) +
  geom_line(aes(x = date, y = fitted), 
            color = my_palette[["black"]],
            data  = ed_usa) + 
  scale_y_continuous(labels = scales::comma,
                     breaks = seq(0, 600000, by=100000)) +
  scale_x_date(date_labels = "%b") +
  ylab("Cumulative excess deaths") +
  xlab("Date") +
  geom_text(aes(make_date(2020,05,01), 40000, label="Cumulative excess deaths"), 
            hjust    = 0, 
            fontface = "bold",
            color    = my_palette[["black"]],
            size     = 3) +
  geom_text(aes(make_date(2020,05,01), 7000, label="Reported Covid-19 deaths"), 
            hjust    = 0, 
            fontface = "bold",
            color    = "gray",
            size     = 3) +
  theme(text = element_text(size = 10))

# -- Saving figure
ggsave(filename = "figures/figure-3d.pdf",
       plot     = fig3d,
       width    = 4,
       height   = 3)
### -- ----------------------------- ------------------------------------------------------------------
### -- Figure: Fhat estimate for USA ------------------------------------------------------------------
### -- ----------------------------- ------------------------------------------------------------------

# -- Figure 3
library(patchwork)
fig3 <- (fig3a | fig3b) / (fig3c | fig3d) + plot_annotation(tag_levels = 'A')

# -- Saving figure
ggsave(filename = "figures/figure-3.pdf",
       plot     = fig3,
       height   = 6,
       width    = 8)
