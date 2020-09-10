# -- Set up
source("code/pr-init.R")
data("new_jersey_counts")
data("louisiana_counts")
data("florida_counts")

# -- Remove the outlier from louisana
louisiana_counts$outcome[which.max(louisiana_counts$outcome)] <- 126

# -- Set up to be used below
ndays  <- 365
knots  <- c(6, 6, 6, 6, 6, 6)
disc   <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
before <- c(365, 365, 365, 365, 365, 365) 
after  <- c(365, 365, 365, 365, 365, 365)
interval_start <- c("Katrina" = ymd("2005-08-29"),
                    "Sandy"   = ymd("2012-10-29"),
                    "Irma"    = ymd("2017-09-10"),
                    hurricane_dates[1], # Hugo
                    hurricane_dates[2], # Georges
                    hurricane_dates[3]) # Maria

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
exclude_dates  <- list(Katrina = interval_start[["Katrina"]]  + -365:365,
                       Sandy   = interval_start[["Sandy"]]  + -365:365,
                       Irma    = interval_start[["Irma"]] + -365:365,
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
    f <- excess_model(counts    = louisiana_counts,
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
    f <- excess_model(counts    = new_jersey_counts,
                      event          = interval_start[[i]],
                      start          = interval_start[[i]] - before[[i]],
                      end            = interval_start[[i]] + after[[i]],
                      exclude        = exclude_dates[[i]],
                      control.dates  = control_dates[[i]],
                      knots.per.year = knots[[i]],
                      weekday.effect = TRUE,
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
                      knots.per.year = knots[i],
                      weekday.effect = TRUE,
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
    
    tmp <- map_df(unique(all_counts$agegroup), function(x)
    {
      cat(".")
      f <- suppressMessages(all_counts %>% 
                              filter(agegroup == x) %>%
                              excess_model(event          = interval_start[[i]],
                                           start          = interval_start[[i]] - before[[i]],
                                           end            = interval_start[[i]] + after[[i]],
                                           exclude        = exclude_dates[[i]],
                                           control.dates  = control_dates[[i]],
                                           knots.per.year = knots[[i]],
                                           weekday.effect = TRUE,
                                           model          = "correlated",
                                           discontinuity  = disc[[i]], 
                                           verbose = FALSE))
      
      excess_cumulative(f, 
                        start = interval_start[[i]], 
                        end   = ymd(interval_start[[i]]) + ndays) %>%
        mutate(agegroup = x, event_day = interval_start[[i]], event = names(interval_start)[[i]]) %>%
        as_tibble()
    })
    # -- Computing marginal effect
    fit <- tmp %>% 
      group_by(date) %>% 
      summarize(fitted = sum(fitted),
                observed = sum(observed),
                sd = sqrt(sum(sd^2)),
                se = sqrt(sum(se^2)),
                event_day = event_day[1], 
                event = event[1]) %>%
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
                                          "NJ: Sandy", "LA: Katrina", "FL: Irma")))

# -- Things to be use in the viz
tmp <- excess_deaths %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  mutate(day = as.numeric(date - event_day)) %>%
  filter(event %in% c("PR: Maria", "PR: Georges", "LA: Katrina")) %>%
  mutate(estimate = case_when(event == "PR: Maria" ~ 3274,
                              event == "PR: Georges" ~ 1213,
                              event == "LA: Katrina" ~ 1428),
         period = case_when(event == "PR: Maria" ~ 191,
                            event == "PR: Georges" ~ 85,
                            event == "LA: Katrina" ~ 108)) %>%
  mutate(fitted = ifelse(event == "LA: Katrina", fitted + 834, fitted),
         lwr = ifelse(event == "LA: Katrina", lwr + 834, lwr),
         upr = ifelse(event == "LA: Katrina", upr + 834, upr)) %>%
  mutate(fitted = ifelse(event == "LA: Katrina" & day == 0, 0, fitted))

# -- For floating legend
lab <- tibble(posx = rep(100, 3), 
              posy = c(100, 300, 500), 
              labs = c("PR: Georges", "LA: Katrina", "PR: Maria"))

# -- Figure
ggplot() +
  geom_ribbon(aes(day, ymin = lwr, ymax = upr, fill=event), color=NA, alpha=0.20, show.legend = F, data = tmp) +
  geom_line(aes(day, fitted, color=event), show.legend = FALSE, data = tmp) +
  geom_point(aes(period, estimate), color="white", size=2.5, show.legend = F, data = unique(select(tmp, period, estimate, event))) +
  geom_point(aes(period, estimate, color=event), size=2, show.legend = F, data = unique(select(tmp, period, estimate, event))) +
  geom_point(aes(period, estimate), color="black", size=2, show.legend = F, pch = 1, data = unique(select(tmp, period, estimate, event))) +
  scale_y_continuous(breaks = seq(0, 6500, by = 1000),
                     labels = scales::comma) +
  geom_text(aes(posx, posy, label=labs, color=labs), hjust=0, size=3, show.legend=F, data=lab) +
  xlab("Days since the event") +
  ylab("Cumulative excess deaths") +
  theme(text = element_text(size=10))

ggsave(filename = "figs/figure-hurricane-excess-deaths.pdf",
       width    = 4,
       height   = 3, 
       dpi      = 300)
### -- ------------------------------ ------------------------------------------------------------------
### -- Figure 2B: Excess deaths in PR ------------------------------------------------------------------
### -- ------------------------------ ------------------------------------------------------------------
