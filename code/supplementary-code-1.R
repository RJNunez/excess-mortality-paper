# -- Set up
source("code/pr-init.R")

# -- Daily mortality for 75 and over
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()


# -- Function for simulations
simulation <- function(event_day, 
                       discontinuity, 
                       num_sims, 
                       order_max = 2,
                       normal    = FALSE, 
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
                      order.max      = order_max,
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
  autores <- ar(counts$outcome[dates %in% control_dates], aic = FALSE, order.max = order_max)
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

# -- Set up for parallelization
library(doParallel)
no_cores <- detectCores() - 2
cl       <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)

start_time <- Sys.time()
end_time <- Sys.time()
end_time - start_time
num_sims <- 100000

# -- Hurricane like event
start_time <- Sys.time()
res_hurricane <- simulation(event_day     = ymd("2017-09-20"),
                            discontinuity = TRUE,
                            num_sims      = 100000)
end_time <- Sys.time()
end_time - start_time
save(res_hurricane, file = "simulation-results/simulation-model-assessment/res-hurricane.rda", compress = "xz")

# -- Epidemic/outbreak like event
start_time <- Sys.time()
res_epidemic <- simulation(event_day     = ymd("2014-07-14"), 
                           discontinuity = FALSE,
                           num_sims      = 100000)
end_time <- Sys.time()
end_time - start_time
save(res_epidemic, file = "simulation-results/simulation-model-assessment/res-epidemic.rda", compress = "xz")

# -- Normal period
start_time <- Sys.time()
res_normal <- simulation(event_day     = ymd("2014-07-14"), 
                         discontinuity = FALSE, 
                         normal        = TRUE, 
                         num_sims      = 100000)
end_time <- Sys.time()
end_time - start_time
save(res_normal, file = "simulation-results/simulation-model-assessment/res-normal.rda", compress = "xz")
stopCluster(cl)

# -- Loading data
load("simulation-results/simulation-model-assessment/res-hurricane.rda")
load("simulation-results/simulation-model-assessment/res-epidemic.rda")
load("simulation-results/simulation-model-assessment/res-normal.rda")

# -- Random sample of ten fits
set.seed(4321)
idx <- sample(1:10000, size = 10)

# -- Hurricane scenario: effect
f1 <- with(res_hurricane,
           ggplot() +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             ylab("Percent increase from expected mortality") +
             coord_cartesian(ylim = c(0, 0.8)) +
             xlab("Date") +
             theme(text = element_text(size = 10)))

# -- Epidemic scenario: effect
f2 <- with(res_epidemic,
           ggplot() +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             ylab("Percent increase from expected mortality") +
             coord_cartesian(ylim = c(0, 0.8)) +
             xlab("Date") +
             theme(text = element_text(size = 10)))

# -- Normal scenario: effect
f3 <- with(res_normal,
           ggplot() +
             geom_hline(yintercept = 0, color="gray", lty=2) +
             geom_line(aes(date, fitted, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(aes(date, true_f), size=0.70, color = my_palette[["red"]], data = unique(select(fitted, date, true_f))) +
             geom_line(aes(date, fitted), size=0.70, color = my_palette[["black"]], lty=2, data = avg_fitted) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.0, 0.8, by=0.10)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             ylab("Percent increase from expected mortality") +
             coord_cartesian(ylim = c(0, 0.8)) +
             xlab("Date") +
             theme(text = element_text(size = 10)))


# -- Hurricane scenario: standard errors
s1 <- with(res_hurricane,
           fitted %>%
             group_by(date) %>%
             summarize(sd_fit = sd(fitted),
                       avg_se = mean(se)) %>%
             ungroup() %>%
             ggplot(aes(date, sd_fit)) +
             geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
             geom_line(size = 0.70, color = my_palette[["red"]]) +
             geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
             scale_y_continuous(labels = scales::percent,
                                breaks = seq(0.02, 0.11, by=0.01)) +
             scale_x_date(date_breaks = "2 months", date_labels = "%b") +
             ylab("Standard errors") +
             coord_cartesian(ylim = c(0.02, 0.10)) +
             xlab("Date") +
             theme(text = element_text(size = 10)))

# -- Epidemic scenario: standard errors
s2 <-with(res_epidemic,
          fitted %>%
            group_by(date) %>%
            summarize(sd_fit = sd(fitted),
                      avg_se = mean(se)) %>%
            ungroup() %>%
            ggplot(aes(date, sd_fit)) +
            geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
            geom_line(size = 0.70, color = my_palette[["red"]]) +
            geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
            scale_y_continuous(labels = scales::percent,
                               breaks = seq(0.02, 0.11, by=0.01)) +
            scale_x_date(date_breaks = "2 months", date_labels = "%b") +
            ylab("Standard errors") +
            coord_cartesian(ylim = c(0.02, 0.10)) +
            xlab("Date") +
            theme(text = element_text(size = 10)))

# -- Normal scenario: standard errors
s3 <-with(res_normal,
          fitted %>%
            group_by(date) %>%
            summarize(sd_fit = sd(fitted),
                      avg_se = mean(se)) %>%
            ungroup() %>%
            ggplot(aes(date, sd_fit)) +
            geom_line(aes(date, se, group = sim), alpha=0.40, color="gray", data = filter(fitted, sim %in% idx)) +
            geom_line(size = 0.70, color = my_palette[["red"]]) +
            geom_line(aes(date, avg_se), size = 0.70, color = my_palette[["black"]], lty = 2) +
            scale_y_continuous(labels = scales::percent,
                               breaks = seq(0.02, 0.11, by=0.01)) +
            scale_x_date(date_breaks = "2 months", date_labels = "%b") +
            ylab("Standard errors") +
            coord_cartesian(ylim = c(0.02, 0.10)) +
            xlab("Date") +
            theme(text = element_text(size = 10)))

library(patchwork)
(f1 | f2 | f3) / (s1 | s2 | s3)

