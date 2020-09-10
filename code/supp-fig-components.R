# -- Set up
source("code/pr-init.R")

# -- Daily mortality for 75 and over
counts <- filter(all_counts, agegroup == "75-Inf")

# -- Simulation function
simulation_components <- function(counts, num_sims)
{
  # -- Dates
  dates <- counts$date
  
  # -- True effect
  f <- rep(0, nrow(counts))
  # if(!normal) { f[dates  %in% fit$date] <- fit$fitted }
  
  # -- Computing expected mortality
  res <- suppressMessages(compute_expected(counts = counts, exclude = exclude_dates, weekday.effect = TRUE, keep.components = TRUE))
  
  # -- Expected counts
  mu <- res$counts$expected
  
  # -- True mortality
  change <- mu * (1 + f)
  
  # -- Random sample of 10
  idx <- sample(1:10, size = 10)
  
  # -- Getting ar component from data with no event
  autores <- ar(counts$outcome[dates %in% control_dates], aic = FALSE, order.max = 7)
  md      <- list(order = c(autores$order, 0, 0), ar = autores$ar)
  
  # -- Simulating data
  message("Simulating data")
  sims <- map_df(1:num_sims, function(i){
    
    # -- Simulated dataset
    sim <- counts
    
    # -- Simulated values from an ARIMA
    sim_ar <- arima.sim(model = md,
                        sd    = 0.05,
                        n     = nrow(counts))
    
    # -- Correlated errors (Centered at more or less 1)
    e <- pmax(0, 1 - sim_ar)
    
    # -- Generating Poisson correlated data
    y <- rpois(nrow(counts), change * e)
    # y <- rpois(nrow(counts), change)
    
    # -- Adding y to the dataset
    sim$expected        <- mu
    sim$expected_change <- change
    sim$deaths          <- y
    sim$sim             <- i
    return(sim)
  })
  
  message("Fitting mean model")
  fits <- lapply(1:num_sims, function(x){
    
    cat(x, "\n")
    
    # -- Temporarily changing name of columns
    sim <- filter(sims, sim == x)
    sim <- sim %>%
      rename(tmp     = outcome,
             outcome = deaths)
    
    fit <- suppressWarnings(suppressMessages(compute_expected(counts = sim, keep.components = TRUE, weekday.effect = TRUE)))
  })

  # -- Trend  
  trend <- map_df(1:num_sims, function(x){ tibble(date = dates, trend = fits[[x]]$trend, sim = x) })
  tmp <- trend %>%
    group_by(date) %>%
    summarize(avg_trend = mean(trend)) %>%
    ungroup() %>%
    mutate(trend = res$trend)
  fig_trend <- ggplot() +
    geom_line(aes(date, trend, group = sim), color="gray", alpha=0.50, data = filter(trend, sim %in% idx)) +
    geom_line(aes(date, trend), color="red", size=0.70, data = tmp) +
    geom_line(aes(date, avg_trend), lty=2, size=0.70, data = tmp) +
    ylab("Average yearly rate of mortality") +
    xlab("Date") +
    theme(text = element_text(size=10))
    
  # -- Seasonal
  seasonal <- map_df(1:num_sims, function(x){ mutate(fits[[x]]$seasonal, sim = x) })
  tmp <- seasonal %>%
    group_by(day) %>%
    summarize(avg_seasonal = mean(s)) %>%
    ungroup() %>%
    mutate(seasonal = res$seasonal$s)
  fig_seasonal <- ggplot() +
    geom_line(aes(day, s, group = sim), color="gray", alpha=0.50, data = filter(seasonal, sim %in% idx)) +
    geom_line(aes(day, seasonal), color="red", data = tmp) +
    geom_line(aes(day, avg_seasonal), lty=2, data = tmp) +
    ylab("Seasonal trend") +
    xlab("Day of the year") +
    theme(text = element_text(size=10))
  
  # -- Weekday
  weekday <- map_df(1:num_sims, function(x){ 
    mutate(fits[[x]]$weekday, sim = x) %>%
      mutate(wnames = factor(c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"), levels = c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")))
    })
  tmp <- weekday %>%
    group_by(weekday, wnames) %>%
    summarize(avg_effect = mean(effect)) %>%
    ungroup() %>%
    left_join(res$weekday, by = "weekday")
  fig_weekday <- ggplot() +
    geom_point(aes(wnames, effect, group = sim), color="gray", alpha=0.50, data = filter(weekday, sim %in% idx)) +
    geom_point(aes(wnames, effect), color="red", size=3, data = tmp) +
    geom_point(aes(wnames, effect), size=3, pch=1, data = tmp) +
    geom_point(aes(wnames, avg_effect), data = tmp) +
    ylab("Weekday effect") +
    xlab("Weekday") +
    theme(text = element_text(size=10))
    
   return(list("fig_trend" = fig_trend, "fig_seasonal" = fig_seasonal, "fig_weekday" = fig_weekday))
}

# -- Running simulations
res <- simulation_components(counts, num_sims = 2000)

# -- Long term trend
ggsave("figs/supp-figure-trend", 
       plot   = res$fig_trend,
       dpi    = 300, 
       width  = 4,
       height = 3)

# -- Seasonal trend
ggsave("figs/supp-figure-seasonal",
       plot   = res$fig_seasonal,
       dpi    = 300, 
       width  = 4,
       height = 3)

# -- Weekday effects
ggsave("figs/supp-figure-weekday", 
       plot   = res$fig_weekday,
       dpi    = 300, 
       width  = 4,
       height = 3)
