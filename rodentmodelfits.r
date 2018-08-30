# Fitting density functions to rodent data.

logbin_setedges <- function(x, y = NULL, bin_edges) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- log10(bin_edges)
  n <- length(bin_edges) - 1
  logxbin <- rep(NA, length(logx))                           # create data structure to assign trees to bins
  b <- bin_edges                                             # add a little to the biggest bin temporarily
  b[length(b)] <- b[length(b)] + 1                           # (so that the biggest single tree is put in a bin)
  for (i in 1:length(logx)) {
    logxbin[i] <- sum(logx[i] >= b)                          # assign each tree to a bin
  }
  bin_midpoints <- 10^(bin_edges[-1] - diff(bin_edges)/2)
  bin_widths <- diff(10^bin_edges)                              # get linear width of each bin
  bin_factor <- factor(logxbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
  bin_counts <- table(bin_factor)                            # find number of trees in each bin
  if (!is.null(y)) {
    rawy <- tapply(y, bin_factor, sum)                       # sum y value (production) in each bin
    rawy[is.na(rawy)] <- 0                                   # add zeroes back in if present
    bin_values <- as.numeric(rawy/bin_widths)                # divide production by width for each bin 
  }
  else {
    bin_values <- as.numeric(bin_counts/bin_widths)          # 1-dimensional case.
  }
  
  return(data.frame(bin_midpoint = bin_midpoints,            # return result!
                    bin_value = bin_values,                  # also add bin min and max for bar plot purposes
                    bin_count = as.numeric(bin_counts),
                    bin_min = 10^bin_edges[1:n],
                    bin_max = 10^bin_edges[2:(n+1)]))
  
}


mammal_data <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/mammal_reduced.csv', stringsAsFactors = FALSE)

library(dplyr)
library(lubridate)

mammal_data <- mammal_data %>%
  filter(!recapture %in% c('U','Y')) %>%
  #filter(!lifeStage %in% c('juvenile')) %>%
  mutate(month = month(collectDate), year = year(collectDate)) %>%
  filter(year == 2016, !is.na(weight)) %>%
  select(domainID, siteID, plotID, collectDate, month, year, taxonID, lifeStage, weight)

mammal_data <- mammal_data %>%
  mutate(relativeEnergy = weight ^ 0.75)

range(mammal_data$weight)
# One is zero so let's get rid of it.
mammal_data <- filter(mammal_data, weight > 0)
(weight_range <- range(mammal_data$weight))

n_bins <- 8

# If we did 8 bins these would be the bin bounds.
(bin_bounds <- exp(seq(log(weight_range[1]), log(weight_range[2]), length.out = n_bins + 1)))

mass_bins <- mammal_data %>%
  group_by(siteID) %>%
  do(logbin_setedges(x = .$weight, bin_edges = bin_bounds))

energy_bins <- mammal_data %>%
  group_by(siteID) %>%
  do(logbin_setedges(x = .$weight, y = .$relativeEnergy, bin_edges = bin_bounds))

n_individuals <- mass_bins %>%
  group_by(siteID) %>%
  summarize(n = sum(bin_count)) %>%
  filter(n >= 100)


#### 
# Do all mammals' mass together in one for now.

x <- mammal_data$weight
# Bin 20 bins
allsite_bin <- logbin_setedges(x = x, bin_edges = exp(seq(log(weight_range[1]), log(weight_range[2]), length.out = 21)))

# Plot mass.
library(ggplot2)
ggplot(data.frame(x), aes(x)) + geom_histogram(bins=20) + scale_x_log10() + scale_y_log10()

# Plot logbinmass
ggplot(subset(allsite_bin,bin_value>0), aes(x=bin_midpoint,y=bin_value)) + geom_point() + scale_x_log10() + scale_y_log10()

# Specify stan model 1. (energy equivalence only; Pareto)

model1 <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
    real<lower=0> x_min;
  }
  
  parameters {
    real<lower=0, upper=5> alpha;
  }
  
  model {
    // Priors (truncated lognormal)
    alpha ~ lognormal(1, 1) T[0, 5];
    // Likelihood (power law)
    x ~ pareto(x_min, alpha);
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = pareto_lpdf(x[i] | x_min, alpha);
    }
  }
'

# Model 1, estimating xmin
model1_xmin <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
  }
  
  parameters {
    real<lower=0, upper=5> alpha;
    real<lower=0, upper=100> x_min;
  }
  
  model {
    // Priors (truncated lognormal)
    alpha ~ lognormal(1, 1) T[0, 5];
    x_min ~ lognormal(1, 1) T[0, 100];
    // Likelihood (power law)
    for (i in 1:N) {
      if (x[i] >= x_min) x[i] ~ pareto(x_min, alpha);
    }
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      if (x[i] >= x_min) log_lik[i] = pareto_lpdf(x[i] | x_min, alpha);
    }
  }
'

# Lomax distribution

model_lomax <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
  }
  
  parameters {
    real<lower=0, upper=50> alpha;
    real<lower=0> lambda;
    real<lower=0> mu;
  }
  
  model {
    // Priors (truncated lognormal)
    alpha ~ lognormal(1, 1) T[0, 50];
    lambda ~ lognormal(1, 1);
    // Likelihood (power law)
    x ~ pareto_type_2(mu, lambda, alpha);
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = pareto_type_2_lpdf(x[i] | mu, lambda, alpha);
    }
  }
'

# Lognormal distribution
model_lognorm <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
  }
  
  parameters {
    real<lower=0> mu;
    real<lower=0> sigma;
  }
  
  model {
    // Priors (lognormal)
    mu ~ lognormal(1, 1);
    sigma ~ lognormal(1, 1);
    // Likelihood (power law)
    x ~ lognormal(mu, sigma);
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      log_lik[i] = lognormal_lpdf(x[i] | mu, sigma);
    }
  }
'

# Specify Stan model 2 (optimal size)

# Edited 11 June: stan model 2 with corrected breakpoint
# Edit so that prior truncations on optimum can be input as data

model2 <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
    real<lower=0> x_min;
    real<lower=0> x_max;
    real<lower=0> x_opt_lim[2];
  }
  
  parameters {
    real<lower=0, upper=5> alpha_low;
    real<lower=0, upper=5> alpha_high;
    // Set x_opt to be uniform within reasonable boundaries
    real<lower=x_opt_lim[1], upper=x_opt_lim[2]> x_opt;
  }
  
  transformed parameters {
    real<lower=0> x_min_high;
    x_min_high = ((alpha_low/alpha_high) * (x_min^alpha_low) * (x_opt^(alpha_high - alpha_low)))^(-alpha_high);
  }
  
  model {
    // Priors (truncated lognormal)
    // Do not set prior on x_opt for now.
    alpha_low ~ lognormal(1, 1) T[0, 5];
    alpha_high ~ lognormal(1, 1) T[0, 5];
    // Likelihood (power law)
    for (i in 1:N) {
      if (x[i] <= x_opt) {
        x[i] ~ pareto(x_min, alpha_low);
      } else {
        x[i] ~ pareto(x_min_high, alpha_high);
      }
    }
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      if (x[i] <= x_opt) {
        log_lik[i] = pareto_lpdf(x[i] | x_min, alpha_low);
      } else {
        log_lik[i] = pareto_lpdf(x[i] | x_min_high, alpha_high);
      }
    }
  }
'

# Specify Stan model 3 (plateau or truncated energy equivalence)

model3 <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
    real<lower=0> x_min;
    real<lower=0> x_max;
    real<lower=0> x_opt_lim[3];
  }
  
  parameters {
    real<lower=0, upper=5> alpha_low;
    real<lower=0, upper=5> alpha_mid;
    real<lower=0, upper=5> alpha_high;
    // Set x_opt1 and x_opt2 to be uniform within reasonable boundaries
    real<lower=x_opt_lim[1], upper=x_opt_lim[2]> x_opt1;
    real<lower=x_opt_lim[2], upper=x_opt_lim[3]> x_opt2;
  }
  
  transformed parameters {
    real<lower=0> x_min_mid;
    real<lower=0> x_min_high;
    x_min_mid = ((alpha_low/alpha_mid) * (x_min^alpha_low) * (x_opt1^(alpha_mid - alpha_low)))^(-alpha_mid);
    x_min_high = ((alpha_mid/alpha_high) * (x_min_mid^alpha_mid) * (x_opt2^(alpha_high - alpha_mid)))^(-alpha_high);
  }


  model {
    // Priors (truncated lognormal)
    // Do not set prior on x_opt for now.
    alpha_low ~ lognormal(1, 1) T[0, 5];
    alpha_mid ~ lognormal(1, 1) T[0, 5];
    alpha_high ~ lognormal(1, 1) T[0, 5];
    // Likelihood (power law)
    for (i in 1:N) {
      if (x[i] <= x_opt1) {
        x[i] ~ pareto(x_min, alpha_low);
      } else if (x[i] > x_opt1 && x[i] <= x_opt2) {
        x[i] ~ pareto(x_min_mid, alpha_mid);
      } else {
        x[i] ~ pareto(x_min_high, alpha_high);
      }
    }
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      if (x[i] <= x_opt1) {
        log_lik[i] = pareto_lpdf(x[i] | x_min, alpha_low);
      } else if (x[i] > x_opt1 && x[i] <= x_opt2) {
        log_lik[i] = pareto_lpdf(x[i] | x_min_mid, alpha_mid);
      } else {
        log_lik[i] = pareto_lpdf(x[i] | x_min_high, alpha_high);
      }
    }
  }
'

library(rstan)
library(bayesplot)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

standata <- list(x = x, N = length(x), x_min = min(x), x_max = max(x))

######## Model 1

stanmodel1 <- stan_model(model_code = model1)

stanfit1 <- sampling(stanmodel1, data = standata, pars = c('alpha'), chains = 3, iter = 2000, warmup = 1000, seed = 919)

summary(stanfit1)$summary
stanfit1_pars <- extract(stanfit1)
mcmc_trace(as.array(stanfit1))

######## Model 1, estimating xmin

stanmodel1xmin <- stan_model(model_code = model1_xmin)

stanfit1xmin <- sampling(stanmodel1xmin, data = list(x=x, N=length(x)), pars = c('alpha', 'x_min'), chains = 3, iter = 2000, warmup = 1000, seed = 919)

summary(stanfit1xmin)$summary
mcmc_trace(as.array(stanfit1xmin))

######## Model 1b: Lomax

stan_lomax <- stan_model(model_code = model_lomax)

stanfit_lomax <- sampling(stan_lomax, data = list(x=x, N=length(x)), pars = c('alpha', 'lambda', 'mu'), chains = 3, iter = 2000, warmup = 1000, seed = 22)

summary(stanfit_lomax)$summary
mcmc_trace(as.array(stanfit_lomax))

######## Model 1c: Lognorm

stan_ln <- stan_model(model_code = model_lognorm)

stanfit_ln <- sampling(stan_ln, data = list(x=x, N=length(x)), pars = c('mu', 'sigma'), chains = 3, iter = 2000, warmup = 1000, seed = 22)

summary(stanfit_ln)$summary
mcmc_trace(as.array(stanfit_ln))

######## Model 2

stanmodel2 <- stan_model(model_code = model2)

stanfit2 <- sampling(stanmodel2, data = c(standata, list(x_opt_lim = c(3, 50))), pars = c('alpha_low', 'alpha_high', 'x_min_high', 'x_opt'), chains = 3, iter = 2000, warmup = 1000, seed = 99)

summary(stanfit2)$summary
stanfit2_pars <- extract(stanfit2)
mcmc_trace(as.array(stanfit2))

######## Model 3

stanmodel3 <- stan_model(model_code = model3)

stanfit3 <- sampling(stanmodel3, data = standata, pars = c('alpha_low', 'alpha_mid', 'alpha_high', 'x_opt1', 'x_opt2'), chains = 3, iter = 7000, warmup = 5000, seed = 313)

summary(stanfit3)
stanfit3_pars <- extract(stanfit3)
mcmc_trace(as.array(stanfit3))

######## Model 2 with fixed optimum

stanmodel2fixed <- stan_model(model_code = model2_fixedoptimum)

stanfit2fixed <- sampling(stanmodel2fixed, data = c(standata, x_opt = 10), pars = c('alpha_low', 'alpha_high'), chains = 3, iter = 7000, warmup = 5000, seed = 1010)

summary(stanfit2fixed)$summary
mcmc_trace(as.array(stanfit2fixed))

######## Maximum likelihood estimation

stanopt1 <- optimizing(stanmodel1, data = standata, seed = 1)
stanopt2 <- optimizing(stanmodel2, data = c(standata, list(x_opt_lim = c(3, 30))), seed = 11)
stanopt3 <- optimizing(stanmodel3, data = c(standata, list(x_opt_lim = c(3, 20, 50))), seed = 1111)
