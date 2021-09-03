# -- Set up
source("code/pr-init.R")
data("new_jersey_counts")
data("louisiana_counts")
# data("florida_counts")
load("florida-counts.rda")
florida_counts <- counts

# -- Remove the outlier from louisana
louisiana_counts$outcome[which.max(louisiana_counts$outcome)] <- 126

# -- Age groups
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)

# -- List of data
counts <- list(florida_counts, louisiana_counts, new_jersey_counts,
               collapse_counts_by_age(puerto_rico_counts, the_breaks))

# -- Hurricane dates
hurricane_dates <- list(irma    = as.Date(c("2017-09-10")),
                        katrina = as.Date(c("2005-08-29")),
                        sandy   = as.Date(c("2012-10-29")),
                        hugo    = as.Date(c("1989-09-18")),
                        georges = as.Date(c("1998-09-21")),
                        maria   = as.Date("2017-09-20"))

# -- Hurricane names
hurricane_names <-  list(irma    = "FL: Irma",
                         katrina = "LA: Katrina",
                         sandy   = "NJ: Sandy",
                         hugo    = "PR: Hugo", 
                         georges = "PR: Georges",
                         maria   = "PR: Maria")

# -- To be used below
count_index <- c(1, 2, 3, 4, 4, 4)

# -- Control dates for each hurricane
control_dates <- list(irma    = seq(make_date(2015, 01, 01), make_date(2016, 12, 31), by = "day"),
                      katrina = seq(make_date(2003, 01, 01), make_date(2005, 08, 01), by = "day"),
                      sandy   = seq(make_date(2007, 01, 01), make_date(2012, 10, 01), by = "day"),
                      hugo    = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"),
                      georges = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"),
                      maria   = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"))

# -- Period of effect for hurricanes in Puerto Rico
puerto_rico_hurricane_dates       <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
puerto_rico_hurricane_effect_ends <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))

# -- Dates to exclude for hurricanes in Puerto Rico
puerto_rico_out_dates <- exclude_dates

# -- Dates to exclude for each hurricane
exclude_dates  <- list(irma    = hurricane_dates[["irma"]] + 0:180,
                       katrina = hurricane_dates[["katrina"]]  + 0:180,
                       sandy   = hurricane_dates[["sandy"]]  + 0:180,
                       hugo    = puerto_rico_out_dates,
                       georges = puerto_rico_out_dates,
                       maria   = puerto_rico_out_dates)

# -- To be used as parameters in model fitting below
before <- days(365)
after  <- days(365)

# -- Number of knots per year
nknots <- 6

# -- Loop to fit models to hurricanes
fits <- map_df(seq_along(count_index), function(i){
  
  # -- Current hurricane
  print(hurricane_names[[i]])
  
  # -- Current index
  h <- count_index[i]
  
  if(i < 4){ # -- Here we fit model to hurricanes outside of PR
    
    # -- Fitting model
    fit <- suppressMessages(excess_model(counts[[h]],
                                         event          = hurricane_dates[[i]],
                                         start          = hurricane_dates[[i]] - before,
                                         end            = hurricane_dates[[i]] + after, 
                                         exclude        = exclude_dates[[i]],
                                         control.dates  = control_dates[[i]],
                                         knots.per.year = nknots,
                                         weekday.effect = TRUE,
                                         verbose        = FALSE,
                                         aic            = FALSE, 
                                         discontinuity  = T,
                                         order.max      = 7,
                                         model          = "correlated"))
    
    # -- Dataframe with fit results
    ret <- tibble(date            = fit$date, 
                  observed        = fit$observed,
                  expected        = fit$expected,
                  log_expected_se = fit$log_expected_se,
                  population      = fit$population,
                  fitted          = fit$fitted, 
                  se              = fit$se,
                  sd              = fit$sd,
                  hurricane       = hurricane_names[[i]],
                  hurricane_date  = hurricane_dates[[i]]) %>%
      select(date, fitted, observed, expected, se, sd, hurricane, hurricane_date)
    
  } else { # -- Here we fit model to hurricanes in PR
    
    # -- Collapsing sex counts
    tmp_counts <- counts[[h]] %>%
      group_by(date, agegroup) %>%
      summarize(outcome    = sum(outcome),
                population = sum(population)) %>%
      ungroup()
    
    # -- This loops through the demographic groups for PR
    tmp <- map_df(unique(counts[[h]]$agegroup), function(x){
      
      # -- Fitting model to each demographic group
      fit <- suppressMessages(tmp_counts %>% 
                                filter(agegroup == x) %>%
                                excess_model(event          = hurricane_dates[[i]],
                                             start          = hurricane_dates[[i]] - before,
                                             end            = hurricane_dates[[i]] + after, 
                                             exclude        = exclude_dates[[i]],
                                             control.dates  = control_dates[[i]],
                                             knots.per.year = nknots,
                                             weekday.effect = TRUE,
                                             verbose        = FALSE,
                                             aic            = FALSE, 
                                             discontinuity  = T,
                                             order.max      = 7,
                                             model          = "correlated"))
      
      # -- Dataframe with fit results for each demographics
      tibble(date            = fit$date, 
             observed        = fit$observed,
             expected        = fit$expected,
             log_expected_se = fit$log_expected_se,
             population      = fit$population,
             fitted          = fit$fitted, 
             se              = fit$se,
             sd              = fit$sd)
    })
    
    # -- Putting PR results together
    ret <- tmp %>%
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
                sd          = sqrt(sum(sd^2)),
                se          = sqrt(sum(tmp_se))) %>%
      ungroup() %>%
      mutate(hurricane      = hurricane_names[[i]],
             hurricane_date = hurricane_dates[[i]]) %>%
      select(date, fitted, observed, expected, se, sd, hurricane, hurricane_date) %>%
      arrange(date)
  }
  return(ret)
}) %>%
  mutate(hurricane = reorder(hurricane, date, min))

# -- Figure 2A
fig2a <- fits %>%
  mutate(lwr = fitted - 1.96 * se, 
         upr = fitted + 1.96 * se,
         day = as.numeric(date - hurricane_date),
         hurricane = factor(hurricane, levels = c("PR: Maria", "PR: Georges", "PR: Hugo",
                                                  "NJ: Sandy", "LA: Katrina", "FL: Irma"))) %>%
  filter(day >= -120, day <= 300) %>%
  ggplot(aes(day, fitted, color = hurricane, fill = hurricane)) +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 0.80, by = 0.20)) +
  ylab("Percent increase from expected mortality") +
  xlab("Days since the hurricane") +
  coord_cartesian(ylim = c(-0.10, 0.80)) +
  geom_hline(yintercept = 0, 
             lty        = 2, 
             color      = "gray") +
  geom_line(size  = 1,
            color = "white") +
  geom_line() +
  scale_color_manual(values = c(my_palette[["red"]]  , my_palette[["blue"]], 
                                my_palette[["green"]], my_palette[["pink"]], 
                                my_palette[["black"]], my_palette[["yellow"]])) +
  theme(text             = element_text(size = 10),
        legend.position  = c(0.80, 0.60),
        legend.direction = "vertical",
        legend.key.size = unit(0.80, "lines"),
        legend.title     = element_blank())

# -- Saving figure
ggsave(filename = "figures/figure-2a.pdf",
       plot     = fig2a,
       height   = 3,
       width    = 4)

# -- Set up
rm(list = ls()[-which(ls() == "fig1b")])
source("code/pr-init.R")
data("new_jersey_counts")
data("louisiana_counts")
# data("florida_counts")
load("florida-counts.rda")
florida_counts <- counts

# -- Remove the outlier from louisana
louisiana_counts$outcome[which.max(louisiana_counts$outcome)] <- 126

# -- Age groups
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)

# -- Set up to be used below
ndays  <- 365
knots  <- c(6, 6, 6, 6, 6, 6)
disc   <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
before <- c(365, 365, 365, 365, 365, 365) 
after  <- c(365, 365, 365, 365, 365, 365)
interval_start <- c("Katrina" = ymd("2005-08-29"),
                    "Sandy"   = ymd("2012-10-29"),
                    "Irma"    = ymd("2017-09-10"),
                    "Hugo"    = hurricane_dates[[1]], # Hugo
                    "Georges" = hurricane_dates[[2]], # Georges
                    "Maria"   = hurricane_dates[[3]]) # Maria

# -- Control dates for each disaster
control_dates <- list(Katrina = seq(make_date(2003, 01, 01), interval_start[["Katrina"]] - 365, by = "day"),
                      Sandy   = seq(make_date(2007, 01, 01), interval_start[["Sandy"]] - 365, by = "day"),
                      Irma    = seq(make_date(2015, 01, 01), interval_start[["Irma"]] - 365, by = "day"),
                      Hugo    = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"),
                      Georges = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"),
                      Maria   = seq(make_date(2002, 01, 01), make_date(2013, 12, 31), by = "day"))

# -- Period of effect for hurricanes in Puerto Rico
puerto_rico_hurricane_dates <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))

# -- Dates to exclude for hurricanes in Puerto Rico
puerto_rico_out_dates <- exclude_dates

# -- Dates to exclude for each hurricane
exclude_dates  <- list(Katrina = interval_start[["Katrina"]] + 0:180,
                       Sandy   = interval_start[["Sandy"]]   + 0:180,
                       Irma    = interval_start[["Irma"]]    + 0:180,
                       Hugo    = puerto_rico_out_dates,
                       Georges = puerto_rico_out_dates,
                       Maria   = puerto_rico_out_dates)

# -- Computing excess deaths
excess_deaths <- map_df(seq_along(interval_start), function(i)
{
  # -- Message of current hurricane
  cat("\nEvent:", names(interval_start)[[i]], "\n")
  message("Fitting model to estimate period of effect")
  
  if(names(interval_start)[[i]] == "Katrina")
  {
    f <- excess_model(counts         = louisiana_counts,
                      event          = interval_start[[i]],
                      start          = interval_start[[i]] - before[[i]],
                      end            = interval_start[[i]] + after[[i]],
                      exclude        = exclude_dates[[i]],
                      control.dates  = control_dates[[i]],
                      knots.per.year = knots[[i]],
                      weekday.effect = TRUE,
                      aic            = FALSE, 
                      order.max      = 7,
                      model          = "correlated",
                      discontinuity  = disc[i], 
                      verbose = FALSE)
    
    fit <- excess_cumulative(f, 
                             start = interval_start[[i]], 
                             end   = ymd(interval_start[[i]]) + ndays) %>%
      mutate(event_day = interval_start[[i]], event = names(interval_start)[[i]]) %>%
      as_tibble()
    
    return(fit)
    # fit <- tibble(date = f$date, expected = f$expected, fitted = f$fitted, se = f$se)
  }
  
  if(names(interval_start)[[i]] == "Sandy")
  {
    f <- excess_model(counts         = new_jersey_counts,
                      event          = interval_start[[i]],
                      start          = interval_start[[i]] - before[[i]],
                      end            = interval_start[[i]] + after[[i]],
                      exclude        = exclude_dates[[i]],
                      control.dates  = control_dates[[i]],
                      knots.per.year = knots[[i]],
                      weekday.effect = TRUE,
                      aic            = FALSE, 
                      order.max      = 7,
                      model          = "correlated",
                      discontinuity  = disc[i], 
                      verbose = FALSE)
    
    fit <- excess_cumulative(f, 
                             start = interval_start[[i]], 
                             end   = ymd(interval_start[[i]]) + ndays) %>%
      mutate(event_day = interval_start[[i]], event = names(interval_start)[[i]]) %>%
      as_tibble()
    
    return(fit)
    # fit <- tibble(date = f$date, expected = f$expected, fitted = f$fitted, se = f$se)
  }
  
  if(names(interval_start)[[i]] == "Irma")
  {
    
    f <- excess_model(counts         = florida_counts,
                      event          = interval_start[[i]],
                      start          = interval_start[[i]] - before[[i]],
                      end            = interval_start[[i]] + after[[i]],
                      exclude        = exclude_dates[[i]],
                      control.dates  = control_dates[[i]],
                      knots.per.year = knots[[i]],
                      weekday.effect = TRUE,
                      aic            = FALSE, 
                      order.max      = 7,
                      model          = "correlated",
                      discontinuity  = disc[i], 
                      verbose = FALSE)
    
    fit <- excess_cumulative(f, 
                             start = interval_start[[i]], 
                             end   = ymd(interval_start[[i]]) + ndays) %>%
      mutate(event_day = interval_start[[i]], event = names(interval_start)[[i]]) %>%
      as_tibble()
    
    return(fit)
    # fit <- tibble(date = f$date, expected = f$expected, fitted = f$fitted, se = f$se)
  }
  
  if(names(interval_start)[[i]] %in% c("Hugo", "Georges", "Maria"))
  {
    # -- Collapsing sex counts
    tmp_counts <- all_counts %>%
      group_by(date, agegroup) %>%
      summarize(outcome    = sum(outcome),
                population = sum(population)) %>%
      ungroup()
    
    # -- This loops through the demographic groups for PR
    tmp <- map_df(unique(all_counts$agegroup), function(x)
    {
      cat(".")
      f <- suppressMessages(tmp_counts %>% 
                              filter(agegroup == x) %>%
                              excess_model(event          = interval_start[[i]],
                                           start          = interval_start[[i]] - before[[i]],
                                           end            = interval_start[[i]] + after[[i]],
                                           exclude        = exclude_dates[[i]],
                                           control.dates  = control_dates[[i]],
                                           knots.per.year = knots[[i]],
                                           weekday.effect = TRUE,
                                           aic            = FALSE, 
                                           order.max      = 7,
                                           model          = "correlated",
                                           discontinuity  = disc[i], 
                                           verbose        = FALSE))
      
      excess_cumulative(f, 
                        start = interval_start[[i]], 
                        end   = ymd(interval_start[[i]]) + ndays) %>%
        mutate(agegroup = x, event_day = interval_start[[i]], event = names(interval_start)[[i]]) %>%
        as_tibble()
    })
    
    tmp %>% 
      group_by(date) %>% 
      summarize(fitted    = sum(fitted),
                observed  = sum(observed),
                sd        = sqrt(sum(sd^2)),
                se        = sqrt(sum(se^2)),
                event_day = event_day[1], 
                event     = event[1]) %>%
      ungroup()
    
    # -- Computing marginal effect
    fit <- tmp %>% 
      group_by(date) %>% 
      summarize(fitted    = sum(fitted),
                observed  = sum(observed),
                sd        = sqrt(sum(sd^2)),
                se        = sqrt(sum(se^2)),
                event_day = event_day[1], 
                event     = event[1]) %>%
      ungroup()
    
    return(fit)
  }
}) %>%
  mutate(event = case_when(event == "Sandy" ~ "NJ: Sandy",
                           event == "Irma" ~ "FL: Irma",
                           event == "Katrina" ~ "LA: Katrina",
                           event == "Hugo" ~ "PR: Hugo",
                           event == "Georges" ~ "PR: Georges",
                           event == "Maria" ~ "PR: Maria"),
         event = factor(event, levels = c("PR: Maria", "PR: Georges", "PR: Hugo",
                                          "NJ: Sandy", "LA: Katrina", "FL: Irma")),
         lwr = fitted-1.96*se,
         upr = fitted+1.96*se,
         day = as.numeric(date - event_day))

# -- Excess deaths estimate for the hurricaes
excess_deaths %>%
  mutate(period = case_when(event == "PR: Maria"   ~ 197,
                            event == "PR: Georges" ~ 90,
                            event == "PR: Hugo"    ~ 12,
                            event == "NJ: Sandy"   ~ 12,
                            event == "FL: Irma"    ~ 48,
                            event == "LA: Katrina" ~ 109)) %>%
  filter(day == period) %>%
  mutate(fitted = ifelse(event == "LA: Katrina", fitted + 834, fitted),
         lwr    = ifelse(event == "LA: Katrina", lwr + 834, lwr),
         upr    = ifelse(event == "LA: Katrina", upr + 834, upr),
         digits = nchar(as.character(round(se))) - 2,
         fitted = round(fitted, digits = -digits),
         lwr    = round(lwr, digits = -digits),
         upr    = round(upr, digits = -digits)) %>%
  select(event, fitted, lwr, upr) %>%
  mutate(event = factor(event, levels = c("PR: Hugo", "PR: Georges", "PR: Maria", "LA: Katrina", "NJ: Sandy", "FL: Irma"))) %>%
  arrange(event)

# -- Adding excess deaths point estimates. To be used in the viz below
excess_deaths_res <- excess_deaths %>%
  filter(event %in% c("PR: Maria", "PR: Georges", "LA: Katrina")) %>%
  mutate(estimate = case_when(event == "PR: Maria"   ~ 3279,
                              event == "PR: Georges" ~ 1296,
                              event == "LA: Katrina" ~ 834 + 734),
         period = case_when(event == "PR: Maria"   ~ 197,
                            event == "PR: Georges" ~ 90,
                            event == "LA: Katrina" ~ 109)) %>%
  mutate(fitted = ifelse(event == "LA: Katrina", fitted + 834, fitted),
         lwr = ifelse(event == "LA: Katrina", lwr + 834, lwr),
         upr = ifelse(event == "LA: Katrina", upr + 834, upr)) %>%
  mutate(fitted = ifelse(event == "LA: Katrina" & day == 0, 0, fitted))

# -- To be used for floating legend in the viz below
lab <- tibble(posx = rep(100, 3), 
              posy = c(100, 300, 500), 
              labs = c("PR: Georges", "LA: Katrina", "PR: Maria"))

# -- Figure 2B
fig2b <- ggplot() +
  geom_ribbon(aes(x    = day, 
                  ymin = lwr, 
                  ymax = upr, 
                  fill = event), 
              color       = NA, 
              alpha       = 0.20, 
              show.legend = F, 
              data        = excess_deaths_res) +
  geom_line(aes(day, fitted, group = event), 
            color       = "white",
            size        = 1,
            show.legend = FALSE, 
            data        = excess_deaths_res) +
  geom_line(aes(day, fitted, color = event), 
            show.legend = FALSE, 
            data        = excess_deaths_res) +
  geom_point(aes(period, estimate, fill = event),
             shape       = 21,
             color       = "white",
             size        = 2,
             show.legend = F, 
             data        = unique(select(excess_deaths_res, period, estimate, event))) +
  xlab("Days since the event") +
  ylab("Cumulative excess deaths") +
  scale_y_continuous(breaks = seq(0, 6500, by = 1000),
                     labels = scales::comma) +
  scale_color_manual(values = c(my_palette[["red"]], my_palette[["blue"]], my_palette[["black"]])) +
  scale_fill_manual(values = c(my_palette[["red"]], my_palette[["blue"]], my_palette[["black"]])) +
  theme(text = element_text(size=10))

# -- Saving figure
ggsave(filename = "figures/figure-2b.pdf",
       plot     = fig2b,
       height   = 3,
       width    = 4)

# -- Figure 1
library(patchwork)
fig2 <- (fig2a | fig2b) + plot_annotation(tag_levels = 'A')

# -- Saving figure
ggsave(filename = "figures/figure-2.pdf",
       plot     = fig2,
       height   = 3,
       width    = 8)

