# -- Set up
source("code/pr-init.R")

# -- Number of knots per year
nknots <- 6

# -- Determining period of interest
event  <- ymd("2014-07-14")
before <- days(365)
after  <- days(365)

# -- Fitting model only to the groups of interest
res <- map_df(levels(all_counts$agegroup), function(x){
  
  tmp_counts <- filter(all_counts, agegroup == x) %>%
    group_by(date, agegroup) %>%
    summarize(outcome    = sum(outcome),
              population = sum(population)) %>%
    ungroup()
  
  tmp_fit    <- excess_model(counts         = tmp_counts, 
                             start          = event - before,
                             end            = event + after,
                             exclude        = exclude_dates,
                             control.dates  = control_dates,
                             knots.per.year = nknots,
                             weekday.effect = TRUE,
                             model          = "correlated",
                             discontinuity  = FALSE)
  
  tibble(date = tmp_fit$date, observed = tmp_fit$observed, expected = tmp_fit$expected,fitted = tmp_fit$fitted, se = tmp_fit$se, agegroup = x)
})

# -- Supplemental figure
res %>%
  filter(date >= "2014-05-01", date <= "2015-05-01") %>%
  mutate(lwr = fitted-1.96*se,
         upr = fitted+1.96*se) %>%
  mutate(agegroup = factor(agegroup, levels = c("0-4", "5-19", "20-39", "40-59", "60-74", "75-Inf"))) %>%
  ggplot(aes(date, fitted)) +
  geom_hline(yintercept = 0, lty=2, color="red2") +
  geom_point(aes(date, (observed/expected - 1)), size=1, alpha=0.10) +
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.30, color=NA, fill="steelblue") +
  geom_line(color="steelblue") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b") +
  ylab("Percent increase from expected mortality") +
  xlab("Date") +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap(~agegroup) +
  theme(text  = element_text(size=10))
