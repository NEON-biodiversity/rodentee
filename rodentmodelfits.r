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


mammal_data <- read.csv('C:/Users/Q/google_drive/NEON_EAGER/Manuscript5_RodentEE/mammal_reduced.csv', stringsAsFactors = FALSE)

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
# Do all mammals' energy together in one for now.

x <- mammal_data$relativeEnergy

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

# Specify Stan model 2 (optimal size)

model2 <- '
  data {
    int<lower=0> N;
    vector<lower=0>[N] x;
    real<lower=0> x_min;
    real<lower=0> x_max;
  }
  
  parameters {
    real<lower=0, upper=5> alpha_low;
    real<lower=0, upper=5> alpha_high;
    // Set x_opt to be uniform within reasonable boundaries
    real<lower=5, upper=30> x_opt;
  }
  
  model {
    // Priors (truncated lognormal)
    // Do not set prior on x_opt for now.
    alpha_low ~ lognormal(1, 1) T[0, 5];
    alpha_high ~ lognormal(1, 1) T[0, 5];
    // Likelihood (power law)
    for (i in 1:N) {
      if (x[N] <= x_opt) {
        x[i] ~ pareto(x_min, alpha_low);
      } else {
        x[i] ~ pareto(x_opt, alpha_high);
      }
    }
  }

  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      if (x[N] <= x_opt) {
        log_lik[i] = pareto_lpdf(x[i] | x_min, alpha_low);
      } else {
        log_lik[i] = pareto_lpdf(x[i] | x_opt, alpha_high);
      }
    }
  }
'

# Stan model 2 with fixed optimum

model2_fixedoptimum <- '
  data {
int<lower=0> N;
vector<lower=0>[N] x;
real<lower=0> x_min;
real<lower=0> x_max;
real<lower=0> x_opt;
}

parameters {
real<lower=0, upper=5> alpha_low;
real<lower=0, upper=5> alpha_high;

}

model {
// Priors (truncated lognormal)
// Do not set prior on x_opt for now.
alpha_low ~ lognormal(1, 1) T[0, 5];
alpha_high ~ lognormal(1, 1) T[0, 5];
// Likelihood (power law)
for (i in 1:N) {
if (x[N] <= x_opt) {
x[i] ~ pareto(x_min, alpha_low);
} else {
x[i] ~ pareto(x_opt, alpha_high);
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
  }
  
  parameters {
    real<lower=0, upper=5> alpha_low;
    real<lower=0, upper=5> alpha_mid;
    real<lower=0, upper=5> alpha_high;
    // Set x_opt1 and x_opt2 to be uniform within reasonable boundaries
    real<lower=5, upper=10> x_opt1;
    real<lower=10, upper=30> x_opt2;
  }
  
  model {
    // Priors (truncated lognormal)
    // Do not set prior on x_opt for now.
    alpha_low ~ lognormal(1, 1) T[0, 5];
    alpha_mid ~ lognormal(1, 1) T[0, 5];
    alpha_high ~ lognormal(1, 1) T[0, 5];
    // Likelihood (power law)
    for (i in 1:N) {
      if (x[N] <= x_opt1) {
        x[i] ~ pareto(x_min, alpha_low);
      } else if (x[N] > x_opt1 && x[N] <= x_opt2) {
        x[i] ~ pareto(x_opt1, alpha_mid);
      } else {
        x[i] ~ pareto(x_opt2, alpha_high);
      }
    }
  }
  
  generated quantities {
    // Log-likelihood (needed for calculating info criteria)
    vector[N] log_lik;
    for (i in 1:N) {
      if (x[N] <= x_opt1) {
        log_lik[i] = pareto_lpdf(x[i] | x_min, alpha_low);
      } else if (x[N] > x_opt1 && x[N] <= x_opt2) {
        log_lik[i] = pareto_lpdf(x[i] | x_opt1, alpha_mid);
      } else {
        log_lik[i] = pareto_lpdf(x[i] | x_opt2, alpha_high);
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

summary(stanfit1)
stanfit1_pars <- extract(stanfit1)
mcmc_trace(as.array(stanfit1))

######## Model 2

stanmodel2 <- stan_model(model_code = model2)

stanfit2 <- sampling(stanmodel2, data = standata, pars = c('alpha_low', 'alpha_high', 'x_opt'), chains = 3, iter = 7000, warmup = 5000, seed = 574)

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
