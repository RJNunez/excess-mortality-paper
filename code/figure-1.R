# -- Set up
source("code/pr-init.R")

# -- Loading PR data from CDC
cdc_pr_dat <- read_csv("Excess_Deaths_Associated_with_COVID-19.csv") %>%
  filter(State == "Puerto Rico", Outcome == "All causes", Type == "Unweighted") %>%
  select(`Week Ending Date`, `Observed Number`, `Upper Bound Threshold`, `Average Expected Count`) %>%
  setNames(c("date", "outcome", "threshold", "expected"))

# -- Figure 2A: Farrington model-based estimates for Puerto Rico
fig1a <- cdc_pr_dat %>%
  filter(date >= "2017-07-01", date <= "2018-07-01") %>%
  ggplot(aes(date, outcome)) +
  coord_cartesian(ylim = c(500, 900)) +
  labs(x = "Date",
       y = "Weekly counts") +
  scale_x_date(date_breaks = "3 months", 
               date_labels = "%b") +
  geom_point(size  = 0.60,
             alpha = 0.10) +
  geom_line(aes(y = expected),
            size  = 0.80,
            color = "white") +
  geom_line(aes(y = expected),
            color = my_palette[["black"]]) +
  geom_line(aes(y = threshold),
            size  = 0.80,
            color = "white") +
  geom_line(aes(y = threshold),
            color = my_palette[["orange"]]) +
  geom_vline(xintercept = c(make_date(2017, 09, 20), make_date(2017, 10, 14)),
             color      = my_palette[["red"]],
             linetype   = 2) +
  annotate("rect",
           xmin = make_date(2017, 09, 20),
           xmax = make_date(2017, 10, 14),
           ymin = 920, 
           ymax = 908,
           fill = my_palette[["red"]]) +
  annotate("rect",
           xmin  = make_date(2017, 09, 20),
           xmax  = make_date(2017, 10, 14),
           ymin  = 920, 
           ymax  = 908,
           fill  = NA,
           color = my_palette[["black"]])

# -- Saving figure
ggsave(filename = "figures/figure-1a.pdf",
       plot     = fig1a,
       height   = 3,
       width    = 4)

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

# -- Dates for Maria examples
maria_start_date <- "2017-07-01"
maria_end_date   <- "2018-07-01"

# -- Acosta & Irizarry daily fit
acosta_maria <- excess_model(counts         = counts,
                             start          = ymd(maria_start_date),
                             end            = ymd(maria_end_date),
                             event          = ymd(hurricane_dates[[3]]),
                             discontinuity  = TRUE,
                             knots.per.year = 6,
                             model          = "correlated",
                             control.dates  = control_dates,
                             exclude        = exclude_dates)

# -- Results from Acosta & Irizarry daily fit
acosta_maria_res <- with(acosta_maria, tibble(date, observed, expected, log_expected_se, fitted, se)) %>%
  mutate(lwr_expected = exp(log(expected) - 2 * log_expected_se),
         upr_expected = exp(log(expected) + 2 * log_expected_se), 
         se_smooth    = sqrt((se^2 * expected^2 * log_expected_se^2) + (fitted^2 * expected^2 * log_expected_se^2) + (expected^2 * se^2)),
         smooth       = expected * (fitted + 1),
         lwr_smooth   = smooth - 2 * se_smooth,
         upr_smooth   = smooth + 2 * se_smooth)

# -- Computing period of indirect effect for our method
indirect_period_maria_acosta <- acosta_maria_res %>%
  filter(date >= hurricane_dates[3], 
         fitted - 2 * se >= 0) %>%
  summarize(mini = min(date),
            maxi = max(date))

# -- Figure 1B: Acosta & Irizarry Maria fit
fig1b <- acosta_maria_res %>%
  ggplot(aes(date, observed)) +
  coord_cartesian(ylim = c(65, 140)) +
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
  geom_vline(xintercept = c(indirect_period_maria_acosta$mini, indirect_period_maria_acosta$maxi),
             color      = my_palette[["red"]],
             linetype   = 2) +
  annotate("rect",
           xmin = indirect_period_maria_acosta$mini, xmax = indirect_period_maria_acosta$maxi,
           ymin = 145, ymax = 142,
           fill = my_palette[["red"]]) +
  annotate("rect",
           xmin  = indirect_period_maria_acosta$mini, xmax = indirect_period_maria_acosta$maxi,
           ymin  = 145, ymax = 142,
           fill  = NA,
           color = my_palette[["black"]])

# -- Saving figure
ggsave(filename = "figures/figure-1b.pdf",
       plot     = fig1b,
       height   = 3,
       width    = 4)

# -- Figure 1
library(patchwork)
fig1 <- (fig1a | fig1b) + plot_annotation(tag_levels = 'A')

# -- Saving figure
ggsave(filename = "figures/figure-1.pdf",
       plot     = fig1,
       height   = 3,
       width    = 8)
