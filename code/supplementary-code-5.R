# -- Set up
source("code/pr-init.R")
library(surveillance)

# -- Daily death counts for Puerto Rico
counts <- puerto_rico_counts %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Getting census population data
pop <- read_csv("pr-population.csv")

# -- Interpolating population values and adding it to the data
counts <- mutate(counts, population = round(approx(x = pop$date, y = pop$population, xout = unique(counts$date), rule = 2)$y))

# -- Creating time series object for daily Farrington fit
farrington <- sts(observed       = matrix(counts$outcome),
                  populationFrac = matrix(counts$population),
                  state          = rep(0, nrow(counts)),
                  freq           = 365,
                  start          = c(1985, 1))

# -- Dates for Maria examples
maria_start_date <- "2017-07-01"
maria_end_date   <- "2018-07-01"

# -- Dates for Georges examples
georges_start_date <- "1998-07-01"
georges_end_date   <- "1999-07-01"

# -- Control list to be use for Farrington daily fit
control_farrington_maria <- list(range                = which(counts$date >= maria_start_date & counts$date <= maria_end_date),
                                 populationOffset     = TRUE,
                                 noPeriods            = 10, 
                                 b                    = 10,
                                 w                    = 3*7,
                                 pastWeeksNotIncluded = 26*7, 
                                 weightsThreshold     = 2.58, 
                                 pThresholdTrend      = 1,
                                 thresholdMethod      = "nbPlugin")

# -- Farrington daily fit
farrington_maria <- farringtonFlexible(farrington, control = control_farrington_maria)

# -- Results
res <- tibble(date                  = counts$date[which(counts$date >= maria_start_date & counts$date <= maria_end_date)],
              observed              = as.numeric(farrington_maria@observed),
              farrington_expected   = as.numeric(farrington_maria@control$expected),
              farrington_upperbound = as.numeric(farrington_maria@upperbound))

# -- Computing period of indirect effect for Farrington
indirect_period_maria_farrington <- res %>%
  filter(date     >= hurricane_dates[3], 
         observed >= farrington_upperbound) %>%
  mutate(date_diff = c(NA, diff(date))) %>%
  filter(date <= "2017-10-08") %>%
  summarize(mini = min(date),
            maxi = max(date))

# -- Viz
res %>%
  ggplot(aes(date, observed)) +
  coord_cartesian(ylim = c(65, 140)) +
  labs(x = "Date",
       y = "Daily counts") +
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b") +
  geom_point(size  = 0.60,
             alpha = 0.10) +
  geom_line(aes(y = farrington_expected), 
            linetype = 1,
            size     = 0.80,
            color    = "white") +
  geom_line(aes(y = farrington_expected), 
            linetype = 1,
            color    = my_palette[["black"]]) +
  geom_line(aes(y = farrington_upperbound), 
            linetype = 1,
            size     = 0.80,
            color    = "white") +
  geom_line(aes(y = farrington_upperbound), 
            linetype = 1,
            color    = my_palette[["orange"]]) +
  geom_vline(xintercept = c(indirect_period_maria_farrington$mini, indirect_period_maria_farrington$maxi),
             color      = my_palette[["red"]],
             linetype   = 2) +
  annotate("rect",
           xmin = indirect_period_maria_farrington$mini, xmax = indirect_period_maria_farrington$maxi,
           ymin = 145, ymax = 142,
           fill = my_palette[["red"]]) +
  annotate("rect",
           xmin  = indirect_period_maria_farrington$mini, xmax = indirect_period_maria_farrington$maxi,
           ymin  = 145, ymax = 142,
           fill  = NA,
           color = my_palette[["black"]])

# -- Control list to be use for Farrington daily fit for Hurricane George
control_farrington_georges <- list(range                = which(counts$date >= georges_start_date & counts$date <= georges_end_date),
                                   populationOffset     = TRUE,
                                   noPeriods            = 10,
                                   b                    = 8,
                                   w                    = 3*7,
                                   pastWeeksNotIncluded = 26*7,
                                   weightsThreshold     = 2.58,
                                   pThresholdTrend      = 1,
                                   thresholdMethod      = "nbPlugin")

# -- Farrington daily fit for Hurricane George
farrington_georges <- farringtonFlexible(farrington, control = control_farrington_georges)

# -- Acosta & Irizarry daily fit for Hurricane George
acosta_georges <- excess_model(counts         = counts,
                               start          = ymd(georges_start_date),
                               end            = ymd(georges_end_date),
                               event          = ymd(hurricane_dates[[2]]),
                               discontinuity  = TRUE,
                               knots.per.year = 6,
                               model          = "correlated",
                               control.dates  = control_dates,
                               exclude        = exclude_dates)

# -- Putting results together
res_georges <- tibble(date                  = acosta_georges$date, 
                      observed              = acosta_georges$observed, 
                      expected              = acosta_georges$expected,
                      log_expected_se       = acosta_georges$log_expected_se,
                      fitted                = acosta_georges$fitted,
                      se                    = acosta_georges$se,
                      farrington_expected   = as.numeric(farrington_georges@control$expected),
                      farrington_upperbound = as.numeric(farrington_georges@upperbound)) %>%
  mutate(lwr_expected = exp(log(expected) - 2 * log_expected_se),
         upr_expected = exp(log(expected) + 2 * log_expected_se), 
         se_smooth    = sqrt((se^2 * expected^2 * log_expected_se^2) + (fitted^2 * expected^2 * log_expected_se^2) + (expected^2 * se^2)),
         smooth       = expected * (fitted + 1),
         lwr_smooth   = smooth - 2 * se_smooth,
         upr_smooth   = smooth + 2 * se_smooth,
         farrington_period = observed >= farrington_upperbound,
         acosta_period     = ifelse(date >= hurricane_dates[3] & fitted - 2 * se >= 0, TRUE, FALSE))

# -- Computing period of indirect effect for Farrington
indirect_period_georges_farrington <- res_georges %>%
  filter(date     >= hurricane_dates[2], 
         observed >= farrington_upperbound) %>%
  mutate(date_diff = c(NA, diff(date))) %>%
  filter(date <= "1998-09-26") %>%
  summarize(mini = min(date),
            maxi = max(date))

# -- Computing period of indirect effect for our method
indirect_period_georges_acosta <- res_georges %>%
  filter(date >= hurricane_dates[2], 
         fitted - 2 * se >= 0) %>%
  summarize(mini = min(date),
            maxi = max(date))

# -- Supplementary Figure: Farrington daily fit for George example
res_georges %>%
  ggplot(aes(date, observed)) +
  coord_cartesian(ylim = c(60, 130)) +
  labs(x = "Date",
       y = "Daily counts") +
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b") +
  geom_point(size  = 0.60,
             alpha = 0.10) +
  geom_line(aes(y = farrington_expected), 
            linetype = 1,
            size     = 0.80,
            color    = "white") +
  geom_line(aes(y = farrington_expected), 
            linetype = 1,
            color    = my_palette[["black"]]) +
  geom_line(aes(y = farrington_upperbound), 
            linetype = 1,
            size     = 0.80,
            color    = "white") +
  geom_line(aes(y = farrington_upperbound), 
            linetype = 1,
            color    = my_palette[["orange"]]) +
  geom_vline(xintercept = c(indirect_period_georges_farrington$mini, indirect_period_georges_farrington$maxi),
             color      = my_palette[["red"]],
             linetype   = 2) +
  annotate("rect",
           xmin = indirect_period_georges_farrington$mini, xmax = indirect_period_georges_farrington$maxi,
           ymin = 135, ymax = 132,
           fill = my_palette[["red"]]) +
  annotate("rect",
           xmin  = indirect_period_georges_farrington$mini, xmax = indirect_period_georges_farrington$maxi,
           ymin  = 135, ymax = 132,
           fill  = NA,
           color = my_palette[["black"]])

# -- Supplementary Figure: Acosta & Irizarry daily fit for George example
res_georges %>%
  ggplot(aes(date, observed)) +
  coord_cartesian(ylim = c(60, 130)) +
  labs(x = "Date",
       y = "Daily counts") +
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b") +
  geom_point(size  = 0.60,
             alpha = 0.10) +
  geom_ribbon(aes(ymin = lwr_expected,
                  ymax = upr_expected),
              fill  = my_palette[["black"]],
              alpha = 0.50) +
  geom_line(aes(y = expected),
            color    = "white",
            size     = 0.80,
            linetype = 1) +
  geom_line(aes(y = expected),
            color    = my_palette[["black"]],
            linetype = 1) +
  geom_ribbon(aes(ymin = lwr_smooth,
                  ymax = upr_smooth),
              fill  = my_palette[["blue"]],
              alpha = 0.50) + 
  geom_line(aes(y = smooth),
            color    = "white",
            size     = 0.80,
            linetype = 1) +
  geom_line(aes(y = smooth),
            color    = my_palette[["blue"]],
            linetype = 1) +
  geom_vline(xintercept = c(indirect_period_georges_acosta$mini, indirect_period_georges_acosta$maxi),
             color      = my_palette[["red"]],
             linetype   = 2) +
  annotate("rect",
           xmin = indirect_period_georges_acosta$mini, xmax = indirect_period_georges_acosta$maxi,
           ymin = 135, ymax = 132,
           fill = my_palette[["red"]]) +
  annotate("rect",
           xmin  = indirect_period_georges_acosta$mini, xmax = indirect_period_georges_acosta$maxi,
           ymin  = 135, ymax = 132,
           fill  = NA,
           color = my_palette[["black"]])

# -- Weekly death counts for Puerto Rico
counts_weekly <- counts %>%
  mutate(date = round_date(date, unit = "week")) %>%
  group_by(date) %>%
  summarize(outcome    = sum(outcome), 
            population = round(mean(population))) %>%
  ungroup() %>%
  filter(outcome >= 250)

# -- Creating time series object for weekly Farrington fit for Georges
farrington_weekly <- sts(observed       = matrix(counts_weekly$outcome),
                         populationFrac = matrix(counts_weekly$population),
                         state          = rep(0, nrow(counts_weekly)),
                         freq           = 52,
                         start          = c(1985, 1))

# -- Control list to be use for Farrington weekly fit for Hurricane George
control_farrington_georges_weekly <- list(range                = which(counts_weekly$date >= georges_start_date & counts_weekly$date <= georges_end_date),
                                          populationOffset     = TRUE,
                                          noPeriods            = 10, 
                                          b                    = 8,
                                          w                    = 3,
                                          pastWeeksNotIncluded = 26, 
                                          weightsThreshold     = 2.58, 
                                          pThresholdTrend      = 1,
                                          thresholdMethod      = "nbPlugin")

# -- Farrington weekly fit for Hurricane George
farrington_georges_weekly <- farringtonFlexible(farrington_weekly, control = control_farrington_georges_weekly)
dres <- tibble(date       = ymd(counts_weekly$date[which(counts_weekly$date >= georges_start_date & counts_weekly$date <= georges_end_date)]),
               observed   = as.numeric(farrington_georges_weekly@observed),
               expected   = as.numeric(farrington_georges_weekly@control$expected),
               upperbound = as.numeric(farrington_georges_weekly@upperbound))

# -- Supplemental Figure: Farrington fit for Hurricane Georges in Puerto Rico using weekly data
dres %>%  
  ggplot(aes(date, observed)) +
  labs(x = "Date",
       y = "Weekly counts") +
  coord_cartesian(ylim = c(500, 750)) +
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b") +
  geom_point(size  = 0.60,
             alpha = 0.10) +
  geom_line(aes(y = expected),
            size  = 0.80,
            color = "white") +
  geom_line(aes(y = expected),
            color = my_palette[["black"]]) +
  geom_line(aes(y = upperbound),
            size  = 0.80,
            color = "white") +
  geom_line(aes(y = upperbound),
            color = my_palette[["orange"]]) +
  geom_vline(xintercept = c(make_date(1998, 09, 21), make_date(1998, 11, 22)),
             color      = my_palette[["red"]],
             linetype   = 2) +
  annotate("rect",
           xmin = make_date(1998, 09, 21),
           xmax = make_date(1998, 11, 22),
           ymin = 758, 
           ymax = 765,
           fill = my_palette[["red"]]) +
  annotate("rect",
           xmin  = make_date(1998, 09, 21),
           xmax  = make_date(1998, 11, 22),
           ymin = 758, 
           ymax = 765,
           fill  = NA,
           color = my_palette[["black"]])

