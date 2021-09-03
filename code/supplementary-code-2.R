###################################################################################################
##################### ------------------------ SET UP ----------------------- #####################
###################################################################################################
# -- Set up
source("code/pr-init.R")

###################################################################################################
################### ------------------------ END SET UP ----------------------- ###################
###################################################################################################

###################################################################################################
##################### -- SIMULATIONS: ASSESING MODEL FIT FOR SMALL COUNTS  -- #####################
###################################################################################################
# -- Simulation function: Assesing model fit for small counts
simulation_low_mu <- function(counts, 
                              event_day, 
                              discontinuity, 
                              num_sims, 
                              string,
                              min_mu    = 0.001,
                              normal    = FALSE,
                              order.max = 2,
                              model     = "correlated")
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
                      order.max      = order.max,
                      discontinuity  = discontinuity)
  
  # -- Dates
  dates <- counts$date
  
  # -- True effect
  f <- rep(0, nrow(counts))
  if(!normal) { f[dates %in% fit$date] <- fit$fitted }
  
  # -- Computing expected mortality
  # res <- suppressMessages(compute_expected(counts = counts, exclude = exclude_dates))
  
  # -- Expected counts
  mu <- res$expected
  
  # -- True mortality
  change <- mu * (1 + f)
  
  # -- Getting ar component from data with no event
  autores <- ar(counts$outcome[dates %in% control_dates], aic = FALSE, order.max = order.max)
  md      <- list(order = c(autores$order, 0, 0), ar = autores$ar)
  
  # -- Finding correct dates
  index <- match(seq(event_day - months(2), event_day + months(7), by = "days"), dates)
  
  # -- Simulating data and fitting models
  message("Simulating data and fitting models")
  fitted <- foreach(i = 1:num_sims) %dopar% {
    
    # -- Simulated dataset
    set.seed(i)
    sim <- counts
    
    # -- Simulated values from an ARIMA
    sim_ar <- arima.sim(model = md,
                        sd    = 0.05,
                        n     = nrow(counts))
    
    # -- Correlated errors (Centered at more or less 1)
    e <- pmax(0, 1 - sim_ar)
    
    # -- Generating Poisson correlated data
    y <- rpois(nrow(counts), (min_mu * change / mean(change)) * e)
    
    # -- Adding y to the dataset
    sim$expected        <- (min_mu * mu / mean(change))
    sim$expected_change <- sim$expected * (1 + f)# change
    sim$deaths          <- y
    sim$sim             <- i
    
    # -- Temporarily changing name of columns
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
    
    # -- Extracting from model object
    with(tmp, 
         tibble(date, expected, fitted, 
                se, sd, dispersion)) %>%
      mutate(sim = i) %>%
      left_join(tibble(date = dates[index], true_f = f[index]), by="date")
  }
  fitted <- bind_rows(fitted)
  
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
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    ylab("Percent increase from expected mortality") +
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
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    ylab("Standard errors") +
    xlab("Date") +
    theme(text = element_text(size = 10))
  
  return(list("fitted" = fitted, "avg_fitted" = avg_fitted, "fig_effect" = fig_effect, "fig_standard_errors" = fig_standard_errors))
}

# -- Daily mortality for 75 and over. We are using these data as the basis for the simulations
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Simulation parameters
num_sims  <- 100000
event_day <- ymd("2014-08-01")

# -- Set up for parallelization
library(doParallel)
no_cores <- detectCores() - 1
cl       <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

# -- Simulation: Assuming 1 expected deaths per day
start_time <- Sys.time()
res_mu_1 <- simulation_low_mu(counts        = counts,
                              event_day     = event_day,
                              min_mu        = 1,
                              discontinuity = FALSE,
                              string        = "Average = 1",
                              num_sims      = num_sims)
end_time <- Sys.time()
end_time - start_time

# -- Saving the results that assumes 1 deaths on average per day
save(res_mu_1, file = "simulation-results/simulation-small-counts/res-mu-1.rda", compress = "xz")

# -- Simulation: Assuming 0.50 expected deaths per day
start_time <- Sys.time()
res_mu_0.50 <- simulation_low_mu(counts        = counts,
                                 event_day     = event_day,
                                 min_mu        = 0.50,
                                 discontinuity = FALSE,
                                 string        = "Average = 0.50",
                                 num_sims      = num_sims)
end_time <- Sys.time()
end_time - start_time

# -- Saving the results that assumes 0.50 deaths on average per day
save(res_mu_0.50, file = "simulation-results/simulation-small-counts/res-mu-0.50.rda", compress = "xz")

# -- Simulation: Assuming 0.10 expected deaths per day
start_time <- Sys.time()
res_mu_0.10 <- simulation_low_mu(counts        = counts,
                                 event_day     = event_day,
                                 min_mu        = 0.10,
                                 discontinuity = FALSE,
                                 string        = "Average = 0.10",
                                 num_sims      = num_sims)
end_time <- Sys.time()
end_time - start_time

# -- Saving the results that assumes 0.10 deaths on average per day
save(res_mu_0.10, file = "simulation-results/simulation-small-counts/res-mu-0.50.rda", compress = "xz")

# -- Simulation: Assuming 0.05 expected deaths per day
start_time <- Sys.time()
res_mu_0.05 <- simulation_low_mu(counts        = counts,
                                 event_day     = event_day,
                                 min_mu        = 0.05,
                                 discontinuity = FALSE,
                                 string        = "Average = 0.05",
                                 num_sims      = num_sims)
end_time <- Sys.time()
end_time - start_time

# -- Saving the results that assumes 0.05 deaths on average per day
save(res_mu_0.05, file = "simulation-results/simulation-small-counts/res-mu-0.05.rda", compress = "xz")
stopCluster(cl)

###################################################################################################
################### -- END SIMULATIONS: ASSESING MODEL FIT FOR SMALL COUNTS  -- ###################
###################################################################################################

###################################################################################################
################ -- SIMULATIONS: ASSESING MODEL FIT FOR YEARS OF TRAINING DATA  -- ################
###################################################################################################
# -- Simulation function: Assesing model fit for different number of years of training data
simulation_training_years <- function(counts, 
                                      training_years, 
                                      weekday_effect,
                                      event_day, 
                                      discontinuity, 
                                      num_sims, 
                                      string,
                                      normal    = FALSE,
                                      order.max = 2, 
                                      model     = "correlated")
{
  
  # -- Fitting model to obtain the true fit
  fit <- excess_model(counts         = counts,
                      event          = event_day,
                      start          = event_day - months(2),
                      end            = event_day + months(7),
                      knots.per.year = 6,
                      exclude        = exclude_dates,
                      control.dates  = control_dates,
                      model          = model,
                      aic            = FALSE,
                      order.max      = order.max,
                      discontinuity  = discontinuity)
  
  # -- Subsetting data
  tmp_counts <- filter(counts, lubridate::year(date) >= 2014 - training_years) %>%
    arrange(date)
  
  # -- Dates
  dates <- tmp_counts$date
  
  # -- True effect
  f <- rep(0, nrow(tmp_counts))
  if(!normal) { f[dates  %in% fit$date] <- fit$fitted }
  
  # -- Computing expected mortality
  res <- suppressMessages(compute_expected(counts         = tmp_counts, 
                                           exclude        = exclude_dates,
                                           weekday.effect = weekday_effect))
  
  # -- Expected counts
  mu <- res$expected
  
  # -- True mortality
  change <- mu * (1 + f)
  
  # -- Getting ar component from data with no event
  autores <- ar(counts$outcome[dates %in% control_dates], aic = FALSE, order.max = order.max)
  md      <- list(order = c(autores$order, 0, 0), ar = autores$ar)
  
  # -- Finding correct dates
  index <- match(seq(event_day - months(2), event_day + months(7), by = "days"), dates)
  
  # -- Simulating data and fitting models
  message("Simulating data and fitting models")
  fitted <- foreach(i = 1:num_sims) %dopar% {
    
    # -- Simulated dataset
    set.seed(i)
    sim <- tmp_counts
    
    # -- Simulated values from an ARIMA
    sim_ar <- arima.sim(model = md,
                        sd    = 0.05,
                        n     = nrow(tmp_counts))
    
    # -- Correlated errors (Centered at more or less 1)
    e <- pmax(0, 1 - sim_ar)
    
    # -- Generating Poisson correlated data
    y <- rpois(nrow(tmp_counts), change * e)
    
    # -- Adding y to the dataset
    sim$expected        <- mu
    sim$expected_change <- change
    sim$deaths          <- y
    sim$sim             <- i
    # -- Temporarily changing name of columns
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
    
    # -- Extracting from model object
    with(tmp, 
         tibble(date, expected, fitted, 
                se, sd, dispersion)) %>%
      mutate(sim = i) %>%
      left_join(tibble(date = dates[index], true_f = f[index]), by="date")
  }
  fitted <- bind_rows(fitted)
  
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
    summarize(sd_fit = sqrt(((n() - 1) / n()) * var(fitted)),
              # sd_fit = sd(fitted),
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

# -- Daily mortality for 75 and over. We are using these data as the basis for the simulations
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Simulation parameters
num_sims  <- 100000
event_day <- ymd("2014-08-01")

# -- Set up for parallelization
library(doParallel)
no_cores <- detectCores() - 2
cl       <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

# -- Simulation: Using 8 years of training data
start_time <- Sys.time()
res_8_years <- simulation_training_years(counts         = counts,
                                         training_years = 8,
                                         weekday_effect = T,
                                         event_day      = event_day,
                                         discontinuity  = FALSE,
                                         string         = "Training years = 8",
                                         num_sims       = num_sims)
Sys.time() - start_time

# -- Saving the results that uses 8 years of training data
save(res_8_years, file = "simulation-results/simulation-training-years/res-8-years.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation: Using 6 years of training data
start_time <- Sys.time()
res_6_years <- simulation_training_years(counts         = counts,
                                         training_years = 6,
                                         weekday_effect = T,
                                         event_day      = event_day,
                                         discontinuity  = FALSE,
                                         string         = "Training years = 6",
                                         num_sims       = num_sims)
Sys.time() - start_time

# -- Saving the results that uses 6 years of training data
save(res_6_years, file = "simulation-results/simulation-training-years/res-6-years.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation: Using 4 years of training data
start_time <- Sys.time()
res_4_years <- simulation_training_years(counts         = counts,
                                         training_years = 4,
                                         weekday_effect = T,
                                         event_day      = event_day,
                                         discontinuity  = FALSE,
                                         string         = "Training years = 4",
                                         num_sims       = num_sims)
Sys.time() - start_time

# -- Saving the results that uses 4 years of training data
save(res_4_years, file = "simulation-results/simulation-training-years/res-4-years.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation: Using 2 years of training data
start_time <- Sys.time()
res_2_years <- simulation_training_years(counts         = counts,
                                         training_years = 2,
                                         weekday_effect = T,
                                         event_day      = event_day,
                                         discontinuity  = FALSE,
                                         string         = "Training years = 2",
                                         num_sims       = num_sims)
Sys.time() - start_time

# -- Saving the results that uses 2 years of training data
save(res_2_years, file = "simulation-results/simulation-training-years/res-2-years.rda", compress = "xz")
Sys.time() - start_time
###################################################################################################
############## -- END SIMULATIONS: ASSESING MODEL FIT FOR YEARS OF TRAINING DATA  -- ##############
###################################################################################################

###################################################################################################
#################### -- VIZ: MEDIAN ABSOLUTE DEVIATION V AVG DEATHS PER DAY  -- ###################
###################################################################################################
# -- Loading simulation results (if these are not present, run the simulations above)
load("simulation-results/simulation-small-counts/res-mu-1.rda")
load("simulation-results/simulation-small-counts/res-mu-0.50.rda")
load("simulation-results/simulation-small-counts/res-mu-0.10.rda")
load("simulation-results/simulation-small-counts/res-mu-0.05.rda")

# -- Fits and mu values
sim_fits <- list(res_mu_1$fitted, res_mu_0.50$fitted, res_mu_0.10$fitted, res_mu_0.05$fitted)
sim_vals <- c(1, 0.50, 0.10, 0.05)

# -- Computing Median Absolute Deviation for simulation fits
res_mad_mu <- map_df(seq_along(sim_fits), function(i){
  
  cat(i, "\n")
  
  # -- SD over all fits
  sd_fit <- sim_fits[[i]] %>%
    group_by(date) %>%
    summarize(sd_fit = sd(fitted)) %>%
    ungroup() %>%
    pull(sd_fit)
  
  # -- B x T matrix with fits per simulation
  A <- sim_fits[[i]] %>%
    select(date, sim, se) %>%
    pivot_wider(names_from = date, values_from = se) %>%
    select(-sim)
  
  # -- B x T matrix with each row equal to SD over all fits
  B <- matrix(sd_fit, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)
  
  # -- RMSE over B and median over T
  rmse_se <- median(sqrt(colMeans((A - B)^2)))
  
  # -- Getting MAD value for SE estimates
  mad_se <- sim_fits[[i]] %>%
    group_by(date) %>%
    summarize(sd_fit = sd(fitted),
              avg_se = mean(se)) %>%
    ungroup() %>%
    summarize(mab_se  = median(abs(avg_se - sd_fit)),
              rmse_se = rmse_se) %>%
    ungroup()
  
  # -- Getting MAD value for fit estimates
  sim_fits[[i]] %>%
    group_by(date) %>%
    mutate(avg_f = mean(fitted)) %>%
    ungroup() %>%
    filter(sim == 1) %>%
    select(date, true_f, avg_f) %>%
    summarize(mab_fit = median(abs(avg_f - true_f))) %>%
    ungroup() %>%
    bind_cols(mad_se)
}) %>%
  mutate(expected = sim_vals) 

# -- Supplementary Figure 2 (A-D)
set.seed(4321)
idx <- sample(1:10000, size = 10)
p1 <- with(res_mu_1,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

p2 <- with(res_mu_0.50,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

p3 <- with(res_mu_0.10,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

p4 <- with(res_mu_0.05,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

# -- Supplementary Figure 3 (A-D)
p1 <- with(res_mu_1,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

p2 <- with(res_mu_0.50,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

p3 <- with(res_mu_0.10,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

p4 <- with(res_mu_0.05,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

###################################################################################################
################## -- END VIZ: MEDIAN ABSOLUTE DEVIATION V AVG DEATHS PER DAY  -- #################
###################################################################################################

###################################################################################################
################## -- VIZ: MEDIAN ABSOLUTE DEVIATION V YEARS OF TRAINING DATA  -- #################
###################################################################################################
# -- Loading simulation results (if these are not present, run the simulations above)
load("simulation-results/simulation-training-years/res-8-years.rda")
load("simulation-results/simulation-training-years/res-6-years.rda")
load("simulation-results/simulation-training-years/res-4-years.rda")
load("simulation-results/simulation-training-years/res-2-years.rda")

# -- Fits and training values
sim_fits <- list(res_8_years$fitted, res_6_years$fitted, res_4_years$fitted, res_2_years$fitted)
sim_vals <- seq(8, 2, by = -2)

# -- Computing Median Absolute Deviation for simulation fits
res_mad_train_years <- map_df(seq_along(sim_fits), function(i){
  
  cat(i, "\n")
  
  # -- SD over all fits
  sd_fit <- sim_fits[[i]] %>%
    group_by(date) %>%
    summarize(sd_fit = sd(fitted)) %>%
    ungroup() %>%
    pull(sd_fit)
  
  # -- B x T matrix with fits per simulation
  A <- sim_fits[[i]] %>%
    select(date, sim, se) %>%
    pivot_wider(names_from = date, values_from = se) %>%
    select(-sim)
  
  # -- B x T matrix with each row equal to SD over all fits
  B <- matrix(sd_fit, nrow = nrow(A), ncol = ncol(A), byrow = TRUE)
  
  # -- RMSE over B and median over T
  rmse_se <- median(sqrt(colMeans((A - B)^2)))
  
  # -- Getting MAD value for SE estimates
  mad_se <- sim_fits[[i]] %>%
    group_by(date) %>%
    summarize(sd_fit = sd(fitted),
              avg_se = mean(se)) %>%
    ungroup() %>%
    summarize(mab_se  = median(abs(avg_se - sd_fit)),
              rmse_se = rmse_se) %>%
    ungroup()
  
  # -- Getting MAD value for fit estimates
  sim_fits[[i]] %>%
    group_by(date) %>%
    mutate(avg_f = mean(fitted)) %>%
    ungroup() %>%
    filter(sim == 1) %>%
    select(date, true_f, avg_f) %>%
    summarize(mab_fit = median(abs(avg_f - true_f))) %>%
    ungroup() %>%
    bind_cols(mad_se)
  
}) %>%
  mutate(years = sim_vals)

# -- Supplementary Figure 4 (A-D)
set.seed(4321)
idx <- sample(1:10000, size = 10)
p1 <- with(res_8_years,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

p2 <- with(res_6_years,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

p3 <- with(res_4_years,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

p4 <- with(res_2_years,
           ggplot() +
             coord_cartesian(ylim = c(-0.05, 0.25)) +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Percent change from expected mortality") +
             theme(text = element_text(size = 10)))

# -- Supplementary Figure 5 (A-D)
p1 <- with(res_8_years,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

p2 <- with(res_6_years,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

p3 <- with(res_4_years,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))

p4 <- with(res_2_years,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             labs(x = "Date",
                  y = "Standard errors") +
             theme(text = element_text(size = 10)))
###################################################################################################
################ -- END VIZ: MEDIAN ABSOLUTE DEVIATION V YEARS OF TRAINING DATA  -- ###############
###################################################################################################