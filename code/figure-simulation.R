# -- Set up
source("code/pr-init.R")

# -- Daily mortality for 75 and over
counts <- filter(all_counts, agegroup == "75-Inf")

# -- Function for simulations
simulation <- function(event_day, discontinuity, num_sims, normal = FALSE, model = "correlated")
{
  # -- Fitting model to obtain the true fit
  fit <- excess_model(counts         = counts,
                      event          = event_day,
                      start          = event_day - months(2),
                      end            = event_day + months(7),
                      knots.per.year = 6,
                      control.dates  = control_dates,
                      exclude        = exclude_dates,
                      model          = model,
                      aic            = FALSE,
                      order.max      = 1,
                      discontinuity  = discontinuity)
  
  # -- Dates
  dates <- counts$date
  
  # -- True effect
  f <- rep(0, nrow(counts))
  if(!normal) { f[dates  %in% fit$date] <- fit$fitted }
  
  # -- Computing expected mortality
  res <- suppressMessages(compute_expected(counts = counts, exclude = exclude_dates))
  
  # -- Expected counts
  mu <- res$expected
  
  # -- True mortality
  change <- mu * (1 + f)

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
  
  # -- Fitting models
  message("Fitting models")
  fits <- lapply(1:num_sims, function(i){
    
    cat(i,"\n")
    
    # -- Temporarily changing name of columns
    sim <- filter(sims, sim == i)
    sim <- sim %>%
      rename(tmp     = outcome,
             outcome = deaths)
    
    # -- Fitting model
    tmp <- suppressMessages(excess_model(counts         = sim,
                                         event          = event_day,
                                         start          = event_day - months(2),
                                         end            = event_day + months(7),
                                         model          = model,
                                         control.dates  = control_dates,
                                         exclude        = exclude_dates,
                                         discontinuity  = discontinuity,
                                         aic            = FALSE,
                                         order.max      = 7,
                                         knots.per.year = 6))
    return(tmp)
  })
  
  # -- Finding correct dates
  index <- match(fits[[1]]$date, dates)
  
  # -- Putting everything together
  fitted <- map_df(1:num_sims, function(i){
    tibble(date       = fits[[i]]$date,
           expected   = fits[[i]]$expected,
           fitted     = fits[[i]]$fitted, 
           se         = fits[[i]]$se, 
           sd         = fits[[i]]$sd,
           dispersion = fits[[i]]$dispersion,
           sigma      = fits[[i]]$ar$sigma,
           # rho        = fits[[i]]$ar$ar,
           sim        = i)
  }) %>%
    left_join(tibble(date = dates[index], true_f = f[index]), by="date")
  
  # -- Computing avg and sd of fitted values 
  avg_fitted <- fitted %>% 
    group_by(date) %>% 
    summarize(se     = sd(fitted), 
              fitted = mean(fitted))
  
  # -- Viz effect size
  tmp <- unique(select(fitted, date, true_f))
  idx <- sample(1:10, size = 10)
  fig_effect <- ggplot() +
    geom_hline(yintercept = 0, color="gray", lty=2) +
    geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
    geom_line(aes(date, true_f), size=0.70, color="red", data = tmp) +
    geom_line(aes(date, fitted), size=0.70, lty=2, data = avg_fitted) +
    scale_y_continuous(labels = scales::percent,
                       breaks = seq(0.0, 0.80, by=0.10)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    ylab("Percent increase from expected mortality") +
    coord_cartesian(ylim = c(0, 0.80)) +
    xlab("Date") +
    theme(text = element_text(size = 10))
  
  # -- Variability estimates
  fig_standard_errors <- fitted %>%
    group_by(date) %>%
    summarize(sd_fit = sd(fitted),
              avg_se = mean(se)) %>%
    ungroup() %>%
    ggplot(aes(date, sd_fit)) +
    geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
    geom_line(size=0.70, color="red") +
    geom_line(aes(date, avg_se), size=0.70, lty=2) +
    scale_y_continuous(labels = scales::percent,
                       breaks = seq(0.02, 0.11, by=0.01)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    ylab("Standard errors") +
    coord_cartesian(ylim = c(0.018, 0.11)) +
    xlab("Date") +
    theme(text = element_text(size = 10))
  
  return(list("fitted" = fitted, "avg_fitted" = avg_fitted, "fig_effect" = fig_effect, "fig_standard_errors" = fig_standard_errors))
}

# -- Hurricane like event
set.seed(4321)
h <- simulation(event_day = ymd("2017-09-20"), discontinuity = TRUE, num_sims = 2000)
ggsave("figs/figure-hurricane-sim-effect.pdf",
       plot   = h$fig_effect,
       dpi    = 300, 
       width  = 4,
       height = 3)
ggsave("figs/figure-hurricane-sim-se.pdf",
       plot   = h$fig_standard_errors,
       dpi    = 300, 
       width  = 4,
       height = 3)

# -- Epidemic/outbreak like event
set.seed(4321)
d <- simulation(event_day = ymd("2014-07-14"), discontinuity = FALSE, num_sims = 2000)
ggsave("figs/figure-epidemic-sim-effect.pdf",
       plot   = d$fig_effect,
       dpi    = 300, 
       width  = 4,
       height = 3)
ggsave("figs/figure-epidemic-sim-se.pdf",
       plot   = d$fig_standard_errors,
       dpi    = 300, 
       width  = 4,
       height = 3)

# -- Normal period
set.seed(4321)
n <- simulation(event_day = ymd("2014-07-14"), discontinuity = F, normal = TRUE, num_sims = 2000)
ggsave("figs/figure-normal-sim-effect.pdf",
       plot   = n$fig_effect,
       dpi    = 300, 
       width  = 4,
       height = 3)
ggsave("figs/figure-normal-sim-se.pdf",
       plot   = n$fig_standard_errors,
       dpi    = 300, 
       width  = 4,
       height = 3)
