# -- Loading data
source("code/pr-init.R")
data("new_jersey_counts")
data("louisiana_counts")
data("florida_counts")

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
effect <- map_df(seq_along(count_index), function(i){
  
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
    ret <- tibble(date           = fit$date, 
                  observed       = fit$observed,
                  expected       = fit$expected,
                  fitted         = fit$fitted, 
                  se             = fit$se, 
                  hurricane      = hurricane_names[[i]],
                  hurricane_date = hurricane_dates[[i]])
    
  } else { # -- Here we fit model to hurricanes in PR
    
    # -- This loops through the demographic groups for PR
    tmp <- map_df(unique(counts[[h]]$agegroup), function(x){
      
      cat(".")
      # -- Fitting model to each demographic group
      fit <- suppressMessages(counts[[h]] %>% 
                                filter(agegroup == x) %>%
                                excess_model(event          = hurricane_dates[[i]],
                                             start          = hurricane_dates[[i]] - before,
                                             end            = hurricane_dates[[i]] + after, 
                                             exclude        = exclude_dates[[i]],
                                             weekday.effect = TRUE,
                                             control.dates  = control_dates[[i]],
                                             knots.per.year = nknots, 
                                             verbose        = FALSE,
                                             model          = "correlated"))
      
      # -- Dataframe with fit results for each demographics
      tibble(date     = fit$date, 
             observed = fit$observed,
             expected = fit$expected, 
             fitted   = fit$fit, 
             se       = fit$se)
    })
    
    # -- Putting PR results together
    ret <- tmp %>% 
      dplyr::group_by(date) %>% 
      dplyr::summarize(fitted    = sum(expected * fitted) / sum(expected), 
                       se        = sqrt(sum(expected^2 * se^2)) / sum(expected),
                       observed  = sum(observed),
                       expected  = sum(expected)) %>%
      dplyr::mutate(hurricane      = hurricane_names[[i]],
                    hurricane_date = hurricane_dates[[i]])
  }
  return(ret)
}) %>%
  mutate(hurricane = reorder(hurricane, date, min),
         day       = as.numeric(date - hurricane_date),
         hurricane = factor(hurricane, levels=rev(unlist(hurricane_names))),
         lwr       = fitted - 1.96 * se,
         upr       = fitted + 1.96 * se)

# -- Table 1 (more below)
table1 <- effect %>% 
  dplyr::mutate(fitted   = round(100 * fitted), 
                lwr      = round(100 * lwr), 
                upr      = round(100 * upr), 
                increase = paste0(fitted, "% (",lwr,"%, ", upr,"%)")) %>%
  dplyr::filter(date == hurricane_date) %>%
  dplyr::rename(landfall = date) %>%
  dplyr::select(hurricane, landfall, increase) %>%
  dplyr::arrange(hurricane, landfall)

# -- Computing indirect effect duration (Katrina's indirect effect duration is a special case)
period <- sapply(table1$hurricane, function(x){
  
  # -- Period
  idx <- effect %>%
    filter(hurricane == x) %>%
    filter(date >= hurricane_date) %>%
    filter(lwr >= 0) %>%
    mutate(d = c(NA, diff(date))) %>%
    select(date, d)
  
  if(all(idx$d == 1, na.rm=T)) {
    id <- last(idx$date)
  } else {
    id <- idx$date[which(idx$d != 1) - 1][1]
  }
  # id
  return(id - first(idx$date))
})
table1 <- cbind(table1, period)
table1[which(table1$hurricane == "LA: Katrina"), 4] <- 108

# -- Set up to be used below for excess deaths
ndays  <- 365
knots  <- c(6, 6, 6, 6, 6, 6)
disc   <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
before <- c(365, 365, 365, 365, 365, 365) 
after  <- c(365, 365, 365, 365, 365, 365)

# -- Computing excess deaths
excess_deaths <- map_df(seq_along(hurricane_dates), function(i)
{
  # -- Message of current hurricane
  cat("\nEvent:", names(hurricane_dates)[[i]], "\n")
  message("Fitting model to estimate period of effect")
  
  if(names(hurricane_dates)[[i]] == "katrina")
  {
    f <- excess_model(counts    = louisiana_counts,
                      event          = hurricane_dates[[i]],
                      start          = hurricane_dates[[i]] - before[[i]],
                      end            = hurricane_dates[[i]] + after[[i]],
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
                             start = hurricane_dates[[i]], 
                             end   = ymd(hurricane_dates[[i]]) + ndays) %>%
      mutate(event_day = hurricane_dates[[i]], event = names(hurricane_dates)[[i]]) %>%
      as_tibble()
    
    return(fit)
  }
  
  if(names(hurricane_dates)[[i]] == "sandy")
  {
    f <- excess_model(counts    = new_jersey_counts,
                      event          = hurricane_dates[[i]],
                      start          = hurricane_dates[[i]] - before[[i]],
                      end            = hurricane_dates[[i]] + after[[i]],
                      exclude        = exclude_dates[[i]],
                      control.dates  = control_dates[[i]],
                      knots.per.year = knots[[i]],
                      weekday.effect = TRUE,
                      model          = "correlated",
                      discontinuity  = disc[i], 
                      verbose = FALSE)
    
    fit <- excess_cumulative(f, 
                             start = hurricane_dates[[i]], 
                             end   = ymd(hurricane_dates[[i]]) + ndays) %>%
      mutate(event_day = hurricane_dates[[i]], event = names(hurricane_dates)[[i]]) %>%
      as_tibble()
    
    return(fit)
  }
  
  if(names(hurricane_dates)[[i]] == "irma")
  {
    
    f <- excess_model(counts         = florida_counts,
                      event          = hurricane_dates[[i]],
                      start          = hurricane_dates[[i]] - before[[i]],
                      end            = hurricane_dates[[i]] + after[[i]],
                      exclude        = exclude_dates[[i]],
                      control.dates  = control_dates[[i]],
                      knots.per.year = knots[i],
                      weekday.effect = TRUE,
                      model          = "correlated",
                      discontinuity  = disc[i],
                      verbose = FALSE)
    
    fit <- excess_cumulative(f, 
                             start = hurricane_dates[[i]], 
                             end   = ymd(hurricane_dates[[i]]) + ndays) %>%
      mutate(event_day = hurricane_dates[[i]], event = names(hurricane_dates)[[i]]) %>%
      as_tibble()
    
    return(fit)
  }
  
  if(names(hurricane_dates)[[i]] %in% c("hugo", "georges", "maria"))
  {
    
    tmp <- map_df(unique(all_counts$agegroup), function(x)
    {
      cat(".")
      f <- suppressMessages(all_counts %>% 
                              filter(agegroup == x) %>%
                              excess_model(event          = hurricane_dates[[i]],
                                           start          = hurricane_dates[[i]] - before[[i]],
                                           end            = hurricane_dates[[i]] + after[[i]],
                                           exclude        = exclude_dates[[i]],
                                           control.dates  = control_dates[[i]],
                                           knots.per.year = knots[[i]],
                                           weekday.effect = TRUE,
                                           model          = "correlated",
                                           discontinuity  = disc[[i]], 
                                           verbose = FALSE))
      
      excess_cumulative(f, 
                        start = hurricane_dates[[i]], 
                        end   = ymd(hurricane_dates[[i]]) + ndays) %>%
        mutate(agegroup = x, event_day = hurricane_dates[[i]], event = names(hurricane_dates)[[i]]) %>%
        as_tibble()
    })
    # -- Computing marginal effect
    fit <- tmp %>% 
      dplyr::group_by(date) %>% 
      dplyr::summarize(fitted    = sum(fitted),
                       observed  = sum(observed),
                       sd        = sqrt(sum(sd^2)),
                       se        = sqrt(sum(se^2)),
                       event_day = event_day[1], 
                       event     = event[1]) %>%
      ungroup()
    return(fit)
  }
}) %>%
  dplyr::mutate(event = case_when(event == "sandy" ~ "NJ: Sandy",
                                  event == "irma" ~ "FL: Irma",
                                  event == "katrina" ~ "LA: Katrina",
                                  event == "hugo" ~ "PR: Hugo",
                                  event == "georges" ~ "PR: Georges",
                                  event == "maria" ~ "PR: Maria"),
         event = factor(event, levels = c("PR: Maria", "PR: Georges", "PR: Hugo",
                                          "NJ: Sandy", "LA: Katrina", "FL: Irma")),
         day = as.numeric(date - event_day),
         lwr = fitted-1.96*se,
         upr = fitted+1.96*se)

# -- Table 1
table1 <- excess_deaths %>%
  dplyr::left_join(table1, by = c("event" = "hurricane")) %>%
  dplyr::group_by(event) %>%
  dplyr::filter(day == period) %>%
  dplyr::mutate(tmp           = 10^(nchar(as.character(round(se/10)))-1),
                fitted        = round_any(fitted, tmp),
                lwr           = round_any(lwr, tmp),
                upr           = round_any(upr, tmp),
                excess_deaths = paste0(round(fitted), " (",round(lwr),", ", round(upr),")")) %>%
  dplyr::select(event, landfall, increase, period, excess_deaths) %>%
  dplyr::arrange(event)
table1[which(table1$event == "LA: Katrina"), 5] <- "1520 (1220, 1820)"
table1
