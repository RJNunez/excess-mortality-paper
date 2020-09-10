# -- Libraries
library(scales)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(excessmort)

# -- Daily data
data("puerto_rico_counts")
the_breaks <- c(0, 5, 20, 40, 60, 75, Inf)
all_counts <- collapse_counts_by_age(puerto_rico_counts, the_breaks) %>% filter(date <= "2020-05-19")

# -- Hurricanes information
hurricane_dates        <- as.Date(c("1989-09-18","1998-09-21","2017-09-20"))
hurricane_effect_ends  <- as.Date(c("1990-03-18","1999-03-21","2018-03-20"))
names(hurricane_dates) <- c("Hugo", "Georges", "Maria")

# -- Control & exclude periods
control_dates <- seq(as.Date("2002-01-01"), as.Date("2013-12-31"), by = "day")
exclude_dates <- c(seq(hurricane_dates[1], hurricane_effect_ends[1], by = "day"),
                   seq(hurricane_dates[2], hurricane_effect_ends[2], by = "day"),
                   seq(hurricane_dates[3], hurricane_effect_ends[3], by = "day"),
                   seq(as.Date("2014-09-01"), as.Date("2015-03-21"), by = "day"),
                   seq(as.Date("2001-01-01"), as.Date("2001-01-15"), by = "day"),
                   seq(as.Date("2020-01-01"), lubridate::today(), by = "day"))

# -- Loading CDC data
dat <- read_csv("~/Documents/Excess_Deaths_Associated_with_COVID-19.csv")
pr  <- filter(dat, State == "Puerto Rico")

# -- Viz Hurricane Maria in PR in 2017
p1 <- pr %>%
  filter(Type == "Predicted (weighted)",
         Outcome == "All causes",
         year(`Week Ending Date`) %in% 2017:2018) %>%
  ggplot(aes(`Week Ending Date`, `Observed Number`)) +
  geom_point(size=1, alpha = 0.30, color="#2171b5") +
  geom_line(aes(`Week Ending Date`, `Average Expected Count`)) +
  geom_line(aes(`Week Ending Date`, `Upper Bound Threshold`), color="orange2", lty=1) +
  geom_point(aes(`Week Ending Date`, `Observed Number` + 10), pch=3, size=2, color="red3", data = filter(pr, `Exceeds Threshold` == TRUE, year(`Week Ending Date`) %in% 2017:2018)) +
  ylab("Weekly counts") +
  xlab("Date") +
  coord_cartesian(ylim = c(450, 900)) +
  scale_x_date(date_labels = "%b %Y") +
  theme_bw() +
  theme(text = element_text(size = 10))

ggsave(filename = "figs/figure-cdc-maria.pdf",
       plot     = p1, 
       width    = 4,
       height   = 3, 
       dpi      = 300)

# -- Fitting our model
res <- map_df(levels(all_counts$agegroup), function(x){
  
  fit <- all_counts %>%
    filter(agegroup == x) %>%
    excess_model(counts = .,
                 start         = make_date(2017, 01, 01),
                 event         = ymd("2017-09-20"),
                 discontinuity = TRUE,
                 end           = make_date(2018, 12, 31),
                 knots.per.year = 6,
                 model         = "correlated",
                 control.dates = control_dates,
                 exclude       = exclude_dates)
  
  tibble(date = fit$date, observed = fit$observed, expected = fit$expected, fitted = fit$fitted, se = fit$se, sd = fit$sd, agegroup = x)
})

# -- Computing statistics to be used below
fit <- res %>% 
  group_by(date) %>% 
  summarize(fitted   = sum(expected * fitted) / sum(expected), 
            se       = sqrt(sum(expected^2 * se^2)) / sum(expected),
            sd       = sqrt(sum(expected^2 * sd^2)) / sum(expected),
            observed = sum(observed),
            expected = sum(expected)) %>%
  ungroup() %>%
  mutate(upr     = fitted + 1.96 * se, 
         lwr     = fitted - 1.96 * se,
         obs_upr = expected + 1.96 * expected * sd, 
         obs_lwr = expected - 1.96 * expected * sd,
         gamma   = fitted * expected + expected)

# -- Computing period of indirect effect
indirect_period <- fit %>%
  filter(date >= hurricane_dates[3], 
         lwr >= 0) %>%
  summarize(mini = min(date),
            maxi = max(date))

# -- Some wrangling for viz
tmp <- pr %>%
  select(`Week Ending Date`, `Average Expected Count`, `Upper Bound Threshold`) %>%
  setNames(c("date", "cdc_expected", "cdc_threshold")) %>%
  mutate(cdc_expected  = cdc_expected / 7,
         cdc_threshold = cdc_threshold / 7) %>%
  filter(year(date) <= 2018)

# -- Figure
p2 <- ggplot() +
  geom_point(aes(date, observed), alpha=0.20, size=1, color="#2171b5", data = fit) +
  geom_line(aes(date, cdc_expected), data = tmp) +
  geom_line(aes(date, cdc_threshold), data = tmp, color="orange2") +
  geom_line(aes(date, expected), lty=2, data = fit) +
  geom_ribbon(aes(date, 
                  ymin = gamma - 1.96 * expected * se,
                  ymax = gamma + 1.96 * expected * se), 
              fill  = "#2171b5", 
              alpha = 0.40,
              data = filter(fit, date < ymd("2017-09-20"))) +
  geom_line(aes(date, gamma), 
            color = "#2171b5",
            size  = 0.70,
            data  = filter(fit, date < ymd("2017-09-20"))) +
  geom_ribbon(aes(date, 
                  ymin = gamma - 1.96 * expected * se,
                  ymax = gamma + 1.96 * expected * se), 
              fill  = "#2171b5", 
              alpha = 0.40,
              data = filter(fit, date >= ymd("2017-09-20"))) +
  geom_line(aes(date, gamma), 
            color = "#2171b5",
            size  = 0.70,
            data  = filter(fit, date >= ymd("2017-09-20"))) +
  ylab("Daily counts") +
  xlab("Date") +
  coord_cartesian(ylim = c(55, 140)) +
  scale_y_continuous(breaks = seq(55, 135, by=20)) +
  scale_x_date(date_labels = "%b %Y") +
  geom_vline(xintercept = indirect_period$mini, lty = 2, color="red3", alpha=0.70) +
  geom_vline(xintercept = indirect_period$maxi, lty = 2, color="red3", alpha=0.70) +
  theme_bw() +
  theme(text = element_text(size = 10)) +
  annotate("rect",
           xmin = indirect_period$mini, xmax = indirect_period$maxi,
           ymin = 143, ymax = 145,
           alpha=0.70, fill="red3")

ggsave(filename = "figs/figure-our-approach-maria.pdf",
       plot     = p2, 
       width    = 4,
       height   = 3, 
       dpi      = 300)
