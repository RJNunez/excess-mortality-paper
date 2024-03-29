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

# -- Hurricanes
vec1 <- c("PR: Maria", "PR: Georges", "PR: Hugo", "NJ: Sandy", "LA: Katrina", "FL: Irma")
vec2 <- c("maria", "georges", "hugo", "sandy", "katrina", "irma")

# -- Saving figures
for(i in seq_along(vec1)){
  
  # -- Figure
  p <- fits %>%
    mutate(lwr = fitted - 1.96 * se, 
           upr = fitted + 1.96 * se,
           day = as.numeric(date - hurricane_date),
           hurricane = factor(hurricane, levels = c("PR: Maria", "PR: Georges", "PR: Hugo",
                                                    "NJ: Sandy", "LA: Katrina", "FL: Irma"))) %>%
    mutate(year  = lubridate::year(date),
           month = month(date),
           day   = day(date),
           year  = ifelse(year %in% c(1989, 1998, 2005, 2012, 2017), 2017, year),
           year  = ifelse(year %in% c(1990, 1999, 2006, 2013, 2018), 2018, year),
           date  = make_date(year, month, day)) %>%
    filter(date >= "2017-06-20", date <= "2018-05-20") %>%
    filter(hurricane == vec1[i]) %>%
    ggplot(aes(date, fitted)) +
    scale_y_continuous(labels = scales::percent,
                       breaks = seq(0, 0.80, by = 0.20)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b") +
    ylab("Percent increase from expected mortality") +
    xlab("Date") +
    coord_cartesian(ylim = c(-0.10, 0.80)) +
    geom_hline(yintercept = 0, 
               lty        = 2, 
               color      = "gray") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), 
                alpha = 0.20, 
                fill  = my_palette[["blue"]]) +
    geom_line(size  = 1,
              color = "white") +
    geom_line(color = my_palette[["blue"]]) +
    theme(text = element_text(size = 10))
  
  print(p)
}

