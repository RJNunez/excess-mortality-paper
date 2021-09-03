###################################################################################################
##################### ------------------------ SET UP ----------------------- #####################
###################################################################################################
# -- Set up
source("code/pr-init.R")

###################################################################################################
################### ------------------------ END SET UP ----------------------- ###################
###################################################################################################

###################################################################################################
################## ------------ SIMULATION: ASSESING TYPE 1 ERROR RATE ---------- #################
###################################################################################################
# -- Function for simulations
simulation_type1_error_rate <- function(counts,
                                        event_day,
                                        num_sims,
                                        nknots,
                                        order_max = 2,
                                        model     = "correlated")
{
  # -- Dates
  dates <- counts$date
  
  # -- True effect
  f <- rep(0, nrow(counts))
  
  # -- Computing expected mortality
  res <- compute_expected(counts = counts, exclude = exclude_dates)
  
  # -- Expected counts
  mu <- res$expected
  
  # -- True mortality
  change <- mu * (1 + f)
  
  # -- Getting ar component from data with no event
  autores <- ar(counts$outcome[dates %in% control_dates], aic = FALSE, order.max = order_max)
  md      <- list(order = c(autores$order, 0, 0), ar = autores$ar)
  
  # -- Finding correct dates
  index <- match(seq(make_date(2014, 01, 01), make_date(2014, 12, 31), by = "days"), dates)
  
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
    y <- rpois(nrow(counts), change * e)
    
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
                                         start          = make_date(2014, 01, 01),
                                         end            = make_date(2014, 12, 31),
                                         model          = model,
                                         control.dates  = control_dates,
                                         exclude        = exclude_dates,
                                         discontinuity  = FALSE,
                                         aic            = FALSE,
                                         order.max      = 7,
                                         knots.per.year = nknots))
    
    # -- Extracting from model object
    with(tmp, 
         tibble(date, observed, expected, log_expected_se,
                fitted, se, sd, dispersion)) %>%
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
  idx <- sample(1:num_sims, size = 10)
  fig_effect <- ggplot() +
    geom_hline(yintercept = 0, color="gray", lty=2) +
    geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
    geom_line(aes(date, true_f), size=0.70, color="red", data = tmp) +
    geom_line(aes(date, fitted), size=0.70, lty=2, data = avg_fitted) +
    scale_y_continuous(labels = scales::percent,
                       breaks = seq(0.0, 0.80, by=0.10)) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
    ylab("Percent increase from expected mortality") +
    # coord_cartesian(ylim = c(0, 0.80)) +
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
    # coord_cartesian(ylim = c(0.018, 0.11)) +
    xlab("Date") +
    theme(text = element_text(size = 10))
  
  return(list("fitted" = fitted, "avg_fitted" = avg_fitted, "fig_effect" = fig_effect, "fig_standard_errors" = fig_standard_errors))
}

# -- Daily death counts for those 75 and older in PR
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Simulation parameters
num_sims  <- 100000
event_day <- make_date(2014, 08, 01)

# -- Set up for parallelization
library(doParallel)
no_cores <- detectCores() - 2
cl       <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

# -- Simulation:
start_time <- Sys.time()
res_6_knots_type1_error <- simulation_type1_error_rate(counts     = counts,
                                                       event_day = event_day,
                                                       nknots    = 6,
                                                       num_sims  = num_sims)
Sys.time() - start_time

save(res_6_knots_type1_error, file = "simulation-results/simulation-type1-error/res-knots-6.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation:
start_time <- Sys.time()
res_12_knots_type1_error <- simulation_type1_error_rate(counts    = counts,
                                                        event_day = event_day,
                                                        nknots    = 12,
                                                        num_sims  = num_sims)
Sys.time() - start_time

save(res_12_knots_type1_error, file = "simulation-results/simulation-type1-error/res-knots-12.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation:
start_time <- Sys.time()
res_363_knots_type1_error <- simulation_type1_error_rate(counts    = counts,
                                                         event_day = event_day,
                                                         nknots    = 363,
                                                         num_sims  = num_sims)
Sys.time() - start_time

save(res_363_knots_type1_error, file = "simulation-results/simulation-type1-error/res-knots-363.rda", compress = "xz")
Sys.time() - start_time
###################################################################################################
################ ------------ END SIMULATION: ASSESING TYPE 1 ERROR RATE ---------- ###############
###################################################################################################

###################################################################################################
######################## ------------ SIMULATION: ASSESING POWER ---------- ######################
###################################################################################################
# -- Function for simulations
simulation_power <- function(counts,
                             event_day,
                             num_sims,
                             nknots,
                             height    = 0.20,
                             width     =  45,
                             order_max = 2,
                             model     = "correlated")
{
  
  # -- Tukey's tri weight function to generate the true f
  tri_weight <- function(u){
    return((1 - abs(u)^3)^3 * (abs(u) <= 1))
  }
  
  # -- Dates
  dates <- counts$date
  
  # -- True effect
  f <- height * tri_weight((as.numeric(dates - event_day)) / width)
  
  # -- Computing expected mortality
  res <- suppressMessages(compute_expected(counts = counts, exclude = exclude_dates))
  
  # -- Expected counts
  mu <- res$expected
  
  # -- True mortality
  change <- mu * (1 + f)
  
  # -- Getting ar component from data with no event
  autores <- ar(counts$outcome[dates %in% control_dates], aic = FALSE, order.max = order_max)
  md      <- list(order = c(autores$order, 0, 0), ar = autores$ar)
  
  # -- Finding correct dates
  index <- match(seq(make_date(2014, 01, 01), make_date(2014, 12, 31), by = "days"), dates)
  
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
    y <- rpois(nrow(counts), change * e)
    
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
                                         start          = make_date(2014, 01, 01),
                                         end            = make_date(2014, 12, 31),
                                         model          = model,
                                         control.dates  = control_dates,
                                         exclude        = exclude_dates,
                                         discontinuity  = FALSE,
                                         aic            = FALSE,
                                         order.max      = 7,
                                         knots.per.year = nknots))
    
    # -- Extracting from model object
    with(tmp, 
         tibble(date, observed, expected, log_expected_se,
                fitted, se, sd, dispersion)) %>%
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

# -- Daily death counts for those 75 and older in PR
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Simulation parameters
num_sims  <- 100000
event_day <- make_date(2014, 08, 01)

# -- Set up for parallelization
library(doParallel)
no_cores <- detectCores() - 2
cl       <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

# -- Simulation:
start_time <- Sys.time()
res_6_knots_power <- simulation_power(counts    = counts,
                                      event_day = event_day,
                                      nknots    = 6,
                                      num_sims  = num_sims)
Sys.time() - start_time

save(res_6_knots_power, file = "simulation-results/simulation-power/res-knots-power-6.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation:
start_time <- Sys.time()
res_12_knots_power <- simulation_power(counts    = counts,
                                       event_day = event_day,
                                       nknots    = 12,
                                       num_sims  = num_sims)
Sys.time() - start_time

save(res_12_knots_power, file = "simulation-results/simulation-power/res-knots-power-12.rda", compress = "xz")
Sys.time() - start_time

# -- Simulation:
start_time <- Sys.time()
res_363_knots_power <- simulation_power(counts    = counts,
                                        event_day = event_day,
                                        nknots    = 363,
                                        num_sims  = num_sims)
Sys.time() - start_time

save(res_363_knots_power, file = "simulation-results/simulation-power/res-knots-power-363.rda", compress = "xz")
Sys.time() - start_time
stopCluster(cl)
###################################################################################################
####################### ------------ END SIMULATION: ASSESING POWER ---------- ####################
###################################################################################################

###################################################################################################
################## ----------- VIZ: TYPE 1 ERROR AR A FUNCTION OF KNOTS --------- #################
###################################################################################################
# -- Loading data
start_time <- Sys.time()
load("simulation-results/simulation-type1-error/res-knots-363.rda")
Sys.time() - start_time

# -- Loading data
start_time <- Sys.time()
load("simulation-results/simulation-type1-error/res-knots-12.rda")
Sys.time() - start_time

# -- Loading data
start_time <- Sys.time()
load("simulation-results/simulation-type1-error/res-knots-6.rda")
Sys.time() - start_time

# -- Putting all the simulations results together
sim_fits <- mutate(res_363_knots_type1_error$fitted, knots = Inf) %>%
  bind_rows(mutate(res_12_knots_type1_error$fitted, knots = 12)) %>%
  bind_rows(mutate(res_6_knots_type1_error$fitted, knots = 6)) %>%
  mutate(lwr  = fitted - 1.96 * se,
         upr  = fitted + 1.96 * se,
         flag = as.numeric(!lwr <= 0 & upr >= 0))

# -- Subsetting to the simulations with excess periods (for computational efficiency)
sim_fits_trimmed <- sim_fits %>%
  group_by(knots, sim) %>%
  filter(sum(flag) != 0) %>%
  ungroup()

# -- Index of the simulations with excess periods (for computational efficiency)
idx_dat <- sim_fits_trimmed %>%
  select(knots, sim) %>%
  unique()

# -- Parallelizing code 
library(doParallel)
no_cores <- detectCores() - 2
cl       <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)

# -- Getting number of excess death periods and length of the periods for knots == 6 fit
start_time <- Sys.time()
res_type1_error_knots_6 <- foreach(i = (1:nrow(filter(idx_dat, knots == 6)))) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- filter(idx_dat, knots == 6)$sim[[i]]
  tmp_knot <- filter(idx_dat, knots == 6)$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(sim_fits_trimmed, sim == tmp_sim, knots == tmp_knot)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
  
}
res_type1_error_knots_6 <- bind_rows(res_type1_error_knots_6)
Sys.time() - start_time
save(res_type1_error_knots_6, file = "simulation-results/simulation-type1-error/res-type1-error-knots-6.rda", compress = "xz")
Sys.time() - start_time

# -- Getting number of excess death periods and length of the periods for knots == 12 fit
start_time <- Sys.time()
res_type1_error_knots_12 <- foreach(i = (1:nrow(filter(idx_dat, knots == 12)))) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- filter(idx_dat, knots == 12)$sim[[i]]
  tmp_knot <- filter(idx_dat, knots == 12)$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(sim_fits_trimmed, sim == tmp_sim, knots == tmp_knot)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
  
}
res_type1_error_knots_12 <- bind_rows(res_type1_error_knots_12)
Sys.time() - start_time
save(res_type1_error_knots_12, file = "simulation-results/simulation-type1-error/res-type1-error-knots-12.rda", compress = "xz")
Sys.time() - start_time

# -- Getting number of excess death periods and length of the periods for knots == INF fit
start_time <- Sys.time()
res_type1_error_knots_Inf <- foreach(i = (1:nrow(filter(idx_dat, knots == Inf)))) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- idx_dat$sim[[i]]
  tmp_knot <- idx_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(sim_fits_trimmed, sim == tmp_sim, knots == tmp_knot)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
  
}
res_type1_error_knots_Inf <- bind_rows(res_type1_error_knots_Inf)
Sys.time() - start_time
save(res_type1_error_knots_Inf, file = "simulation-results/simulation-type1-error/res-type1-error-knots-Inf.rda", compress = "xz")
Sys.time() - start_time
stopCluster(cl)

# -- Supp table
load("simulation-results/simulation-type1-error/res-type1-error-knots-6.rda")
load("simulation-results/simulation-type1-error/res-type1-error-knots-12.rda")
load("simulation-results/simulation-type1-error/res-type1-error-knots-Inf.rda")
res_type1_error <- bind_rows(res_type1_error_knots_6, res_type1_error_knots_12, res_type1_error_knots_Inf)
num_sims <- 10e4

tab <- res_type1_error %>%
  group_by(knots) %>%
  summarize(`>= 1 days`   = sum(length_ed_period >= 1) / num_sims,
            `>= 3 days`   = sum(length_ed_period >= 3) / num_sims,
            `>= 5 days`   = sum(length_ed_period >= 5) / num_sims,
            `>= 10 days`  = sum(length_ed_period >= 10) / num_sims,
            `>= 1 month`  = sum(length_ed_period >= 30) / num_sims,
            `>= 2 months` = sum(length_ed_period >= 60) / num_sims) %>%
  mutate(knots = ifelse(knots == Inf, "Saturated", knots))

###################################################################################################
################ ----------- END VIZ: TYPE 1 ERROR AR A FUNCTION OF KNOTS --------- ###############
###################################################################################################

###################################################################################################
###################### ----------- VIZ: POWER AS A FUNCTION OF KNOTS --------- ####################
###################################################################################################
# -- Loading data
start_time <- Sys.time()
load("simulation-results/simulation-power/res-knots-power-363.rda")
Sys.time() - start_time

# -- Loading data
start_time <- Sys.time()
load("simulation-results/simulation-power/res-knots-power-12.rda")
Sys.time() - start_time

# -- Loading data
start_time <- Sys.time()
load("simulation-results/simulation-power/res-knots-power-6.rda")
Sys.time() - start_time

# -- Putting all the simulations results together
sim_fits <- mutate(res_363_knots_power$fitted, knots = Inf) %>%
  bind_rows(mutate(res_12_knots_power$fitted, knots = 12)) %>%
  bind_rows(mutate(res_6_knots_power$fitted, knots = 6))  %>%
  mutate(lwr  = fitted - 1.96 * se,
         upr  = fitted + 1.96 * se,
         flag = as.numeric(!lwr <= 0 & upr >= 0))

# -- Subsetting to the simulations with excess periods (for computational efficiency)
sim_fits_trimmed <- sim_fits %>%
  group_by(knots, sim) %>%
  filter(sum(flag) != 0) %>%
  ungroup()

# -- Index of the simulations with excess periods (for computational efficiency)
idx_dat <- sim_fits_trimmed %>%
  select(knots, sim) %>%
  unique()

# -- Getting period of true excess mortality
true_dates <- select(res_12_knots_power$fitted, date, true_f, sim) %>%
  filter(true_f != 0, sim == 1) %>%
  pull(date)

# -- Parallelizing code 
library(doParallel)
no_cores <- detectCores() - 2
cl       <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)

# -- Getting number of excess death periods and length of the periods for knots == 6 fit

K         <- nrow(filter(idx_dat, knots == 12))
k_cuts    <- c(seq(1, K, by = K / 10) - 1, K)
k_cuts[1] <- 1

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice(k_cuts[1]:k_cuts[2])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part1 <- foreach(i = ( k_cuts[1]:k_cuts[2] )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part1 <- bind_rows(res_power_knots_6_part1)
Sys.time() - start_time
save(res_power_knots_6_part1, file = "simulation-results/simulation-power/res-power-knots-6-part1.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[2]+1):k_cuts[3])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part2 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part2 <- bind_rows(res_power_knots_6_part2)
Sys.time() - start_time
save(res_power_knots_6_part2, file = "simulation-results/simulation-power/res-power-knots-6-part2.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[3]+1):k_cuts[4])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part3 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part3 <- bind_rows(res_power_knots_6_part3)
Sys.time() - start_time
save(res_power_knots_6_part3, file = "simulation-results/simulation-power/res-power-knots-6-part3.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[4]+1):k_cuts[5])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part4 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part4 <- bind_rows(res_power_knots_6_part4)
Sys.time() - start_time
save(res_power_knots_6_part4, file = "simulation-results/simulation-power/res-power-knots-6-part4.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[5]+1):k_cuts[6])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part5 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part5 <- bind_rows(res_power_knots_6_part5)
Sys.time() - start_time
save(res_power_knots_6_part5, file = "simulation-results/simulation-power/res-power-knots-6-part5.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[6]+1):k_cuts[7])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part6 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part6 <- bind_rows(res_power_knots_6_part6)
Sys.time() - start_time
save(res_power_knots_6_part6, file = "simulation-results/simulation-power/res-power-knots-6-part6.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[7]+1):k_cuts[8])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part7 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part7 <- bind_rows(res_power_knots_6_part7)
Sys.time() - start_time
save(res_power_knots_6_part7, file = "simulation-results/simulation-power/res-power-knots-6-part7.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[8]+1):k_cuts[9])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part8 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part8 <- bind_rows(res_power_knots_6_part8)
Sys.time() - start_time
save(res_power_knots_6_part8, file = "simulation-results/simulation-power/res-power-knots-6-part8.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[9]+1):k_cuts[10])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part9 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part9 <- bind_rows(res_power_knots_6_part9)
Sys.time() - start_time
save(res_power_knots_6_part9, file = "simulation-results/simulation-power/res-power-knots-6-part9.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 6) %>% slice((k_cuts[10]+1):k_cuts[11])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 6, sim %in% tmp_dat$sim)
res_power_knots_6_part10 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_6_part10 <- bind_rows(res_power_knots_6_part10)
Sys.time() - start_time
save(res_power_knots_6_part10, file = "simulation-results/simulation-power/res-power-knots-6-part10.rda", compress = "xz")
Sys.time() - start_time

res_power_knots_6 <- bind_rows(res_power_knots_6_part1,
                               res_power_knots_6_part2,
                               res_power_knots_6_part3,
                               res_power_knots_6_part4,
                               res_power_knots_6_part5,
                               res_power_knots_6_part6,
                               res_power_knots_6_part7,
                               res_power_knots_6_part8,
                               res_power_knots_6_part9,
                               res_power_knots_6_part10)
save(res_power_knots_6, file = "simulation-results/simulation-power/res-power-knots-6.rda", compress = "xz")

# -- Getting number of excess death periods and length of the periods for knots == 12 fit
K         <- nrow(filter(idx_dat, knots == 12))
k_cuts    <- c(seq(1, K, by = K / 10) - 1, K)
k_cuts[1] <- 1

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice(k_cuts[1]:k_cuts[2])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part1 <- foreach(i = ( k_cuts[1]:k_cuts[2] )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part1 <- bind_rows(res_power_knots_12_part1)
Sys.time() - start_time
save(res_power_knots_12_part1, file = "simulation-results/simulation-power/res-power-knots-12-part1.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[2]+1):k_cuts[3])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part2 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part2 <- bind_rows(res_power_knots_12_part2)
Sys.time() - start_time
save(res_power_knots_12_part2, file = "simulation-results/simulation-power/res-power-knots-12-part2.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[3]+1):k_cuts[4])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part3 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part3 <- bind_rows(res_power_knots_12_part3)
Sys.time() - start_time
save(res_power_knots_12_part3, file = "simulation-results/simulation-power/res-power-knots-12-part3.rda", compress = "xz")
Sys.time() - start_time


# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[4]+1):k_cuts[5])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part4 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part4 <- bind_rows(res_power_knots_12_part4)
Sys.time() - start_time
save(res_power_knots_12_part4, file = "simulation-results/simulation-power/res-power-knots-12-part4.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[5]+1):k_cuts[6])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part5 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part5 <- bind_rows(res_power_knots_12_part5)
Sys.time() - start_time
save(res_power_knots_12_part5, file = "simulation-results/simulation-power/res-power-knots-12-part5.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[6]+1):k_cuts[7])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part6 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part6 <- bind_rows(res_power_knots_12_part6)
Sys.time() - start_time
save(res_power_knots_12_part6, file = "simulation-results/simulation-power/res-power-knots-12-part6.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[7]+1):k_cuts[8])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part7 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part7 <- bind_rows(res_power_knots_12_part7)
Sys.time() - start_time
save(res_power_knots_12_part7, file = "simulation-results/simulation-power/res-power-knots-12-part7.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[8]+1):k_cuts[9])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part8 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part8 <- bind_rows(res_power_knots_12_part8)
Sys.time() - start_time
save(res_power_knots_12_part8, file = "simulation-results/simulation-power/res-power-knots-12-part8.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[9]+1):k_cuts[10])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part9 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part9 <- bind_rows(res_power_knots_12_part9)
Sys.time() - start_time
save(res_power_knots_12_part9, file = "simulation-results/simulation-power/res-power-knots-12-part9.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == 12) %>% slice((k_cuts[10]+1):k_cuts[11])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == 12, sim %in% tmp_dat$sim)
res_power_knots_12_part10 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_12_part10 <- bind_rows(res_power_knots_12_part10)
Sys.time() - start_time
save(res_power_knots_12_part10, file = "simulation-results/simulation-power/res-power-knots-12-part10.rda", compress = "xz")
Sys.time() - start_time

res_power_knots_12 <- bind_rows(res_power_knots_12_part1,
                                res_power_knots_12_part2,
                                res_power_knots_12_part3,
                                res_power_knots_12_part4,
                                res_power_knots_12_part5,
                                res_power_knots_12_part6,
                                res_power_knots_12_part7,
                                res_power_knots_12_part8,
                                res_power_knots_12_part9,
                                res_power_knots_12_part10)
save(res_power_knots_12, file = "simulation-results/simulation-power/res-power-knots-12.rda", compress = "xz")

# -- Getting number of excess death periods and length of the periods for knots == Inf fit
K      <- nrow(filter(idx_dat, knots == Inf))
k_cuts <- c(seq(1, K, by = floor(K / 10)))
k_cuts[length(k_cuts)] <- K

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice(k_cuts[1]:k_cuts[2])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part1 <- foreach(i = ( k_cuts[1]:k_cuts[2] )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part1 <- bind_rows(res_power_knots_Inf_part1)
Sys.time() - start_time
save(res_power_knots_Inf_part1, file = "simulation-results/simulation-power/res-power-knots-Inf-part1.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[2]+1):k_cuts[3])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part2 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part2 <- bind_rows(res_power_knots_Inf_part2)
Sys.time() - start_time
save(res_power_knots_Inf_part2, file = "simulation-results/simulation-power/res-power-knots-Inf-part2.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[3]+1):k_cuts[4])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part3 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part3 <- bind_rows(res_power_knots_Inf_part3)
Sys.time() - start_time
save(res_power_knots_Inf_part3, file = "simulation-results/simulation-power/res-power-knots-Inf-part3.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[4]+1):k_cuts[5])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part4 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part4 <- bind_rows(res_power_knots_Inf_part4)
Sys.time() - start_time
save(res_power_knots_Inf_part4, file = "simulation-results/simulation-power/res-power-knots-Inf-part4.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[5]+1):k_cuts[6])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part5 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part5 <- bind_rows(res_power_knots_Inf_part5)
Sys.time() - start_time
save(res_power_knots_Inf_part5, file = "simulation-results/simulation-power/res-power-knots-Inf-part5.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[6]+1):k_cuts[7])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part6 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part6 <- bind_rows(res_power_knots_Inf_part6)
Sys.time() - start_time
save(res_power_knots_Inf_part6, file = "simulation-results/simulation-power/res-power-knots-Inf-part6.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[7]+1):k_cuts[8])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part7 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part7 <- bind_rows(res_power_knots_Inf_part7)
Sys.time() - start_time
save(res_power_knots_Inf_part7, file = "simulation-results/simulation-power/res-power-knots-Inf-part7.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[8]+1):k_cuts[9])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part8 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part8 <- bind_rows(res_power_knots_Inf_part8)
Sys.time() - start_time
save(res_power_knots_Inf_part8, file = "simulation-results/simulation-power/res-power-knots-Inf-part8.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[9]+1):k_cuts[10])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part9 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part9 <- bind_rows(res_power_knots_Inf_part9)
Sys.time() - start_time
save(res_power_knots_Inf_part9, file = "simulation-results/simulation-power/res-power-knots-Inf-part9.rda", compress = "xz")
Sys.time() - start_time

# --
start_time <- Sys.time()
tmp_dat   <- filter(idx_dat, knots == Inf) %>% slice((k_cuts[10]+1):k_cuts[11])
tmp_dat_trimmed <- filter(sim_fits_trimmed, knots == Inf, sim %in% tmp_dat$sim)
res_power_knots_Inf_part10 <- foreach(i = ( 1:nrow(tmp_dat) )) %dopar% {
  
  # -- Getting simulation and knot value
  tmp_sim  <- tmp_dat$sim[[i]]
  tmp_knot <- tmp_dat$knots[[i]]
  
  # -- Excess deaths period info
  excess_death_periods <- rle(filter(tmp_dat_trimmed, sim == tmp_sim, date %in% true_dates)$flag)
  
  # -- Number of periods with excess deaths
  num_periods <- length(which(excess_death_periods$values == 1))
  
  # -- Number of consecutive days flagged as excess deaths
  length_periods <- excess_death_periods$lengths[which(excess_death_periods$values == 1)]
  
  # -- Gathering resutls
  tibble(sim = tmp_sim, knots = tmp_knot, num_ed_periods = num_periods, length_ed_period = length_periods)
}
res_power_knots_Inf_part10 <- bind_rows(res_power_knots_Inf_part10)
Sys.time() - start_time
save(res_power_knots_Inf_part10, file = "simulation-results/simulation-power/res-power-knots-Inf-part10.rda", compress = "xz")
Sys.time() - start_time

res_power_knots_Inf <- bind_rows(res_power_knots_Inf_part1,
                                 res_power_knots_Inf_part2,
                                 res_power_knots_Inf_part3,
                                 res_power_knots_Inf_part4,
                                 res_power_knots_Inf_part5,
                                 res_power_knots_Inf_part6,
                                 res_power_knots_Inf_part7,
                                 res_power_knots_Inf_part8,
                                 res_power_knots_Inf_part9,
                                 res_power_knots_Inf_part10)
save(res_power_knots_Inf, file = "simulation-results/simulation-power/res-power-knots-Inf.rda", compress = "xz")
stopCluster(cl)

# -- Possible supplemental figure table (THIS IS THE TABLE)
load("simulation-results/simulation-power/res-power-knots-6.rda")
load("simulation-results/simulation-power/res-power-knots-12.rda")
load("simulation-results/simulation-power/res-power-knots-Inf.rda")
res_power <- bind_rows(res_power_knots_6,
                       res_power_knots_12,
                       res_power_knots_Inf)
tab <- lapply(list(1,3,5,10,30,60), function(x){
  
  if(x <= 5)
  {
    tmp_tab <- res_power %>%
      filter(length_ed_period >= x) %>%
      select(sim, knots) %>%
      unique() %>%
      group_by(knots) %>%
      summarize(tmp = n()  / num_sims) %>%
      ungroup() %>%
      select(tmp)
  } else{
    tmp_tab <- res_power %>%
      filter(length_ed_period >= x) %>%
      select(sim, knots) %>%
      unique() %>%
      group_by(knots) %>%
      summarize(tmp = n()  / num_sims) %>%
      ungroup() %>%
      select(tmp)
    tmp_tab <- rbind(tmp_tab, 0)
  }
})
tab <- do.call(cbind, tab)
tab <- cbind(c("6", "12", "Saturated"), tab)
colnames(tab) <- c("knots", ">= 1 day", ">= 3 days", ">= 5 days", ">= 10 days", ">= 1 month", ">= 2 month")
###################################################################################################
###################### ----------- END VIZ: POWER AS A FUNCTION OF KNOTS --------- ################
###################################################################################################