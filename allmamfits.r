# Using all mammals from 2016, fit models corresponding to the three hypotheses.
# Use information criteria to differentiate between them.

library(rstan)
library(bayesplot)
library(loo)
library(purrr)
library(dplyr)
library(lubridate)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

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


mammal_data <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/mammal_reduced.csv', stringsAsFactors = FALSE)

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
# Test this on a single site, Konza.

x <- subset(mammal_data, siteID == 'KONZ')$weight

standata <- list(x = x, N = length(x), x_min = min(x), x_max = max(x))

stanmod1 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h1.stan') # Hypothesis 1: power law (global ee)
stanmod2 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h2.stan') # Hypothesis 2: peaked
stanmod3 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h3.stan') # Hypothesis 3: tail+EE+tail

n_iter <- 5000
n_warm <- 4000
n_chain <- 3

# Initial value function for fit 3 to ensure that the parameters have realistic values
init_h3 <- function() {
  require(truncdist)
  tau_low <- runif(1, standata$x_min, standata$x_max/4)
  tau_high <- 2*tau_low
  alphas <- rtrunc(3, 'lnorm', 0, 2)
  list(tau_low=tau_low, tau_high=tau_high, alpha_low=alphas[1], alpha_mid=alphas[2], alpha_high=alphas[3])
  
}

fit1 <- sampling(stanmod1, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 11111)
fit2 <- sampling(stanmod2, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 22222)
fit3 <- sampling(stanmod3, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 303, control = list(max_treedepth = 25, adapt_delta = 0.9), init = init_h3)


summary(fit1)$summary[1,]
summary(fit2)$summary[1:3,]
summary(fit3)$summary[1:5,]

parnames1 <- c('alpha')
parnames2 <- c('alpha_low', 'alpha_high', 'tau')
parnames3 <- c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high')

mcmc_trace(as.array(fit1), pars = parnames1)
mcmc_trace(as.array(fit2), pars = parnames2)
mcmc_trace(as.array(fit3), pars = parnames3)

# Information criteria.
ll1 <- extract_log_lik(fit1)
ll2 <- extract_log_lik(fit2)
ll3 <- extract_log_lik(fit3)

loo1 <- loo(ll1)
loo2 <- loo(ll2)
loo3 <- loo(ll3)

# Get fitted values with credible intervals using all parameter values
x_pred <- exp(seq(log(min(x)), log(max(x)), length.out = 51)) # Vector of x values to get fitted values

# Functions to get fitted value
pdf_mod1 <- function(x, xmin, alpha) (alpha * xmin^alpha) / (x ^ (alpha+1))
pdf_mod2 <- function(x, xmin, alpha_low, alpha_high, tau) {
  C_con = tau ^ -(alpha_high + alpha_low)
  C_norm = ( (C_con / alpha_low) * (tau ^ alpha_low - xmin ^ alpha_low) + ( tau ^ (-alpha_high) ) / alpha_high ) ^ -1
  
  ifelse(x < tau, 
         C_con * C_norm * ( x ^ (alpha_low - 1) ),
         C_norm * ( x ^ - (alpha_high + 1) ) )
}
pdf_mod3 <- function(x, xmin, alpha_low, alpha_mid, alpha_high, tau_low, tau_high) {
  C_con_low = tau_low ^ -(alpha_mid + alpha_low)
  C_con_high = tau_high ^ (alpha_high - alpha_mid)
  C_norm = ( (C_con_low / alpha_low) * (tau_low ^ alpha_low - xmin ^ alpha_low) + (1 / alpha_mid) * (tau_low ^ -alpha_mid - tau_high ^ -alpha_mid) + (C_con_high / alpha_high) * (tau_high ^ -alpha_high) ) ^ -1
  
  ifelse(x < tau_low,
         C_con_low * C_norm * ( x ^ (alpha_low - 1) ),
         ifelse(x > tau_high,
                C_con_high * C_norm * ( x ^ - (alpha_high + 1) ),
                C_norm * ( x ^ - (alpha_mid + 1) ) ))
}

# Extract parameters
pars1 <- as.data.frame(do.call('cbind', extract(fit1, parnames1)))
pars2 <- as.data.frame(do.call('cbind', extract(fit2, parnames2)))
pars3 <- as.data.frame(do.call('cbind', extract(fit3, parnames3)))

# Calculate fitted values
fitted1 <- length(x) * pmap_dfc(pars1, pdf_mod1, x = x_pred, xmin = standata$x_min)
fitted2 <- length(x) * pmap_dfc(pars2, pdf_mod2, x = x_pred, xmin = standata$x_min)
fitted3 <- length(x) * pmap_dfc(pars3, pdf_mod3, x = x_pred, xmin = standata$x_min)

# Get quantiles of fitted values
qprobs <- c(0.025, 0.5, 0.975)
fittedquant1 <- t(apply(fitted1, 1, quantile, probs = qprobs))
fittedquant2 <- t(apply(fitted2, 1, quantile, probs = qprobs))
fittedquant3 <- t(apply(fitted3, 1, quantile, probs = qprobs))

# Make fitted value df
fittedquant_all <- data.frame(x = x_pred,
                              model = rep(c('H1', 'H2', 'H3'), each = length(x_pred)),
                              rbind(fittedquant1, fittedquant2, fittedquant3)) %>%
  setNames(c('x','model','q025','q50','q975'))

# Make a plot

library(ggplot2)

ggplot(mass_bins %>% filter(siteID == 'KONZ', bin_value > 0)) +
  geom_ribbon(aes(x = x, ymin = q025, ymax = q975, fill = model, group = model), data = fittedquant_all, alpha = 0.6) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  geom_line(aes(x = x, y = q50, color = model, group = model), data = fittedquant_all) +
  scale_x_log10(name = 'Mass (g)') + scale_y_log10(name = 'Density (individuals/g)') +
  theme_bw()

# Make a plot showing the information criterion results for the site
# Mean +/- SE

loo_all <- data.frame(model = c('H1','H2','H3'),
                      LOOIC = c(loo1$estimates['looic','Estimate'], loo2$estimates['looic','Estimate'], loo3$estimates['looic','Estimate']),
                      se_LOOIC = c(loo1$estimates['looic','SE'], loo2$estimates['looic','SE'], loo3$estimates['looic','SE']))

ggplot(loo_all, aes(x=model, y=LOOIC, ymin=LOOIC-se_LOOIC, ymax=LOOIC+se_LOOIC)) +
  geom_pointrange() + theme_bw()
