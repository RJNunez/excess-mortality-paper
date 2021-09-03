### -- ------------------------------------------------------ -----------------------------------------------------
### -- Supp Figure: QQ-plot without adjusting for correlation -----------------------------------------------------
### -- ------------------------------------------------------ -----------------------------------------------------
# -- Set up 
source("code/pr-init.R")

# -- Wrangling mortality data
counts <- filter(all_counts, agegroup == "75-Inf") %>%
  group_by(date, agegroup) %>%
  summarize(outcome = sum(outcome),
            population = sum(population)) %>%
  ungroup()

# -- Computing expected mortality counts
counts <- compute_expected(counts, exclude = exclude_dates, weekday.effect = TRUE)

# -- Dates for model check
example_dates <- control_dates[lubridate::year(control_dates) > 2005]
first(example_dates)
last(example_dates)

# -- Computing z scores for example dates (H0 true)
r <- tibble(date = counts$date, observed = counts$outcome, expected = counts$expected) %>%
  filter(date %in% example_dates) %>%
  mutate(r = (observed - expected)/sqrt(expected)) %>%
  pull(r)

# -- Supp Figure: Poisson model doesn't fit the tails
tibble(r=r) %>%
  ggplot(aes(sample=r)) +
  stat_qq(alpha = 0.50,
          color = my_palette[["black"]]) + 
  geom_abline(intercept = 0, 
              slope     = 1,
              color     = my_palette[["red"]], 
              lty       = 2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  scale_y_continuous(limits = c(-4, 5),
                     breaks = seq(-4, 5, by=2)) +
  scale_x_continuous(limits = c(-4, 4),
                     breaks = seq(-4, 4, by=2)) +
  theme(text = element_text(size=10))
### -- ------------------------------------------------------ -----------------------------------------------------
### -- Supp Figure: QQ-plot without adjusting for correlation -----------------------------------------------------
### -- ------------------------------------------------------ -----------------------------------------------------

### -- ----------------------------------------- -----------------------------------------------------
### -- Supp Figure: Show correlation in the data -----------------------------------------------------
### -- ----------------------------------------- -----------------------------------------------------
# -- Supp Figure
auto_cor <- acf(r, plot=FALSE)
tibble(acf = auto_cor$acf, lag = auto_cor$lag) %>%
  ggplot(aes(lag, acf)) +
  geom_col(color = "black", 
           fill  = my_palette[["black"]], 
           width = 0.5) +
  ylab("ACF") +
  xlab("Lag") +
  geom_hline(yintercept = qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             color    = my_palette[["red"]],
             linetype = 2) +
  geom_hline(yintercept = -qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             color    = my_palette[["red"]],
             linetype = 2) +
  theme(text  = element_text(size=10))

# -- The sd is not 1
sd(r)
### -- ----------------------------------------- -----------------------------------------------------
### -- Supp Figure: Show correlation in the data -----------------------------------------------------
### -- ----------------------------------------- -----------------------------------------------------

### -- ---------------------------------------------- -----------------------------------------------------
### -- Supp Figure: QQ-plot adjusting for correlation -----------------------------------------------------
### -- ---------------------------------------------- -----------------------------------------------------
# -- Now show that we if adjust for correlation things get better
tmp <- counts %>% 
  as_tibble() %>%
  filter(date %in% example_dates) %>%
  mutate(r = (outcome - expected)/expected)

# -- Extracting percent change
r <- tmp$r

# -- Expected value
mu <- tmp$expected

# -- Number of data points
n <- length(r)

# -- Estimating AR(p) with control dates
arfit <- excessmort:::fit_ar(tmp,  control.dates = control_dates)

# -- Estimated rhos
rhos <- ARMAacf(ar = arfit$ar, ma = 0, lag.max = n)

# -- Estimated sd
s <- arfit$sigma

# -- Variance of log(mu)
log_mu_vari <- tmp$log_expected_se^2

# -- Computing covariance matrix
Sigma <- apply(abs(outer(1:n, 1:n, "-")) + 1, 1, function(i) rhos[i]) * outer(sqrt(s^2 + 1/mu + log_mu_vari), sqrt(s^2 + 1/mu + log_mu_vari))

# -- Cholesky decomposition of Sigma & computing inverse
V     <- chol(Sigma)
V_inv <- backsolve(r = V, x = diag(ncol(V))) ## V is upper triangular so backsolve faster

# -- Supp Figure
tibble(r=V_inv %*% r) %>%
  ggplot(aes(sample=r)) +
  stat_qq(alpha = 0.50,
          color = my_palette[["black"]]) + 
  geom_abline(intercept = 0, 
              slope     = 1, 
              color     = my_palette[["red"]], 
              lty       = 2) +
  ylab("Sample quantiles") +
  xlab("Theoretical quantiles") +
  scale_y_continuous(limits = c(-4, 5),
                     breaks = seq(-4, 5, by=2)) +
  scale_x_continuous(limits = c(-4, 4),
                     breaks = seq(-4, 4, by=2)) +
  theme(text = element_text(size=10))

### -- ---------------------------------------------- -----------------------------------------------------
### -- Supp Figure: QQ-plot adjusting for correlation -----------------------------------------------------
### -- ---------------------------------------------- -----------------------------------------------------

### -- ------------------------------------------------- -----------------------------------------------------
### -- Supp Figure: No correlation in adjusted residuals -----------------------------------------------------
### -- ------------------------------------------------- -----------------------------------------------------
# -- Supp Figure
auto_cor   <- acf(V_inv %*% r, plot = FALSE)
tibble(acf = auto_cor$acf, lag = auto_cor$lag) %>%
  ggplot(aes(lag, acf)) +
  geom_col(color = "black", 
           fill  = my_palette[["black"]], 
           width = 0.5) +
  ylab("ACF") +
  xlab("Lag") +
  geom_hline(yintercept = qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             color      = my_palette[["red"]],
             linetype   = 2) +
  geom_hline(yintercept = -qnorm((1+0.95)/2)/sqrt(auto_cor$n.used), 
             color      = my_palette[["red"]],
             linetype   = 2) +
  theme(text = element_text(size=10))

# -- Estimated sd
sd(V_inv %*% r)
### -- ------------------------------------------------- -----------------------------------------------------
### -- Supp Figure: No correlation in adjusted residuals -----------------------------------------------------
### -- ------------------------------------------------- -----------------------------------------------------