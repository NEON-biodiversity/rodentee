# Fit 3 models to each site and year, and save the results.
# Does not need to be done in parallel for the moment because it runs fairly quickly.
# QDR/NEON Rodents/31 Aug 2018

# Modified 13 Sep: get rid of the non target taxa using NEON taxonomy

# Load and process data ---------------------------------------------------


library(rstan)
library(bayesplot)
library(loo)
library(purrr)
library(dplyr)
library(lubridate)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

# Standard logbinning algorithm to give the same boundaries to every site
logbin_setedges <- function(x, y = NULL, bin_edges) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- log10(bin_edges)
  n <- length(bin_edges) - 1
  
  logxbin <- findInterval(logx, bin_edges, rightmost.closed = TRUE) # assign each individual to a bin
  
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

# New logbinning algorithm to give different ones to every site
logbin_split <- function(x, y = NULL, bin_edges) {
  logx <- log10(x)                                           # log transform x value (biomass)
  bin_edges <- log10(bin_edges)
  
  logxbin <- findInterval(logx, bin_edges, rightmost.closed = TRUE) # assign each individual to a bin
  
  n <- length(bin_edges) - 1 # Desired number of nonzero bins
  
  # While loop to keep splitting the biggest bin in half until we get n nonzero bins
  while(length(unique(logxbin)) < n) {
    uniquebin <- unique(logxbin)
    n_mode <- uniquebin[which.max(tabulate(match(logxbin, uniquebin)))] # Which logbin is the mode (most individuals)?
    # Add a midpoint in the middle of that bin.
    bin_edges <- c(bin_edges[1:n_mode], mean(bin_edges[n_mode:(n_mode+1)]), bin_edges[(n_mode+1):length(bin_edges)])
    logxbin <- findInterval(logx, bin_edges, rightmost.closed = TRUE)
  }
  
  n <- length(bin_edges) - 1 # New n
  
  bin_midpoints <- 10^(bin_edges[-1] - diff(bin_edges)/2)
  bin_widths <- diff(10^bin_edges)                              # get linear width of each bin
  bin_factor <- factor(logxbin, levels=1:n)                  # convert bin to factor (required to deal with zeroes if present)
  bin_counts <- table(bin_factor)                            # find number of individuals in each bin
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


source('load_mammal_data_withimputed.r')




# Look at where there are too few individuals. Get rid of anything less than 50 for now and with less than 4 species.
n_individuals <- mammal_data %>%
  group_by(siteID, year) %>%
  summarize(n=n(), site_richness = length(unique(taxonID)), xmin = min(weight), xmax = max(weight))

mammal2016 <- mammal_data %>%
  left_join(n_individuals) %>%
  filter(n >= 50, site_richness >= 4, year == 2016)

# Sum up total relative energy by site so the fits can be normalized
total_energy <- mammal2016 %>%
  group_by(siteID, year) %>%
  summarize(total = sum(metabolicRate))

# bin up the mass and energy. 
# Modification 21 Sep: Use 10 bins and use the new binning algorithm.
n_bins <- 10
(mass_range <- range(mammal2016$weight))

(bin_bounds <- exp(seq(log(mass_range[1]), log(mass_range[2]), length.out = n_bins + 1)))

mass_bins <- mammal2016 %>%
  group_by(siteID, year) %>%
  do(logbin_split(x = .$weight, bin_edges = bin_bounds))

energy_bins <- mammal2016 %>%
  group_by(siteID, year) %>%
  do(logbin_split(x = .$weight, y = .$metabolicRate, bin_edges = bin_bounds))


# Load stan models --------------------------------------------------------


stanmod1 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h1.stan') # Hypothesis 1: power law (global ee)
stanmod2 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h2.stan') # Hypothesis 2: peaked
stanmod3 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_h3.stan') # Hypothesis 3: tail+EE+tail

n_iter <- 5000
n_warm <- 4000
n_chain <- 3

parnames1 <- c('alpha')
parnames2 <- c('alpha_low', 'alpha_high', 'tau')
parnames3 <- c('alpha_low', 'alpha_mid', 'alpha_high', 'tau_low', 'tau_high')

# Initial value function for fit 3 to ensure that the parameters have realistic values
init_h3 <- function() {
  require(truncdist)
  tau_low <- runif(1, standata$x_min, standata$x_max/4)
  tau_high <- 2*tau_low
  alphas <- rtrunc(3, 'lnorm', 0, 2)
  list(tau_low=tau_low, tau_high=tau_high, alpha_low=alphas[1], alpha_mid=alphas[2], alpha_high=alphas[3])
  
}


# Fit the models! ---------------------------------------------------------

# Make a big loop to do each site separately
# It isn't that long to run so the whole thing can probably be run on a laptop.
# Just do 2016 for now.

fit_model <- function(dat, mod, ctrl = NULL, init = 'random', random_seed = NULL) {
  x <- dat$weight
  standata <<- list(x = x, N = length(x), x_min = min(x), x_max = max(x))
  if (is.null(random_seed)) random_seed <- strtoi(dat$siteID[1], 36) + dat$year[1] + 101
  sampling(mod, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = random_seed, control = ctrl, init = init)
}

all_fits_mod1 <- mammal2016 %>%
  group_by(siteID, year) %>%
  do(fit = fit_model(., mod = stanmod1))

all_fits_mod2 <- mammal2016 %>%
  group_by(siteID, year) %>%
  do(fit = fit_model(., mod = stanmod2))

all_fits_mod3 <- mammal2016 %>%
  group_by(siteID, year) %>%
  do(fit = fit_model(., mod = stanmod3, ctrl = list(max_treedepth = 25, adapt_delta = 0.9), init = init_h3))

# Refit some of the models that did not fit well the first time
bad_mod2 <- c(3,6,12)
bad_mod3 <- c(1,2,16,21,8,15,17)

n_iter <- 10000
n_warm <- 8000

bad_mod2 <- c(3)
bad_mod3 <- c(17)

n_iter <- 20000
n_warm <- 15000

for (i in bad_mod2) {
  print(i)
  all_fits_mod2$fit[[i]] <- fit_model(dat = mammal2016 %>% filter(siteID == all_fits_mod2$siteID[i], year == all_fits_mod2$year[i]),
                                      mod = stanmod2, random_seed = 444)
}


for (i in bad_mod3) {
  print(i)
  all_fits_mod3$fit[[i]] <- fit_model(dat = mammal2016 %>% filter(siteID == all_fits_mod2$siteID[i], year == all_fits_mod2$year[i]),
                                      mod = stanmod3, ctrl = list(max_treedepth = 25, adapt_delta = 0.9), init = init_h3, random_seed = 333)
}

summary(all_fits_mod2$fit[[3]])$summary[1:3,]
mcmc_trace(as.array(all_fits_mod2$fit[[3]]), pars = c('alpha_low', 'alpha_high', 'tau'))

# Get parameter values and CIs --------------------------------------------

get_pars <- function(fit, parnames) {
  summ <- try(summary(fit)$summary[parnames,,drop = FALSE], TRUE)
  if (is.null(summ) | inherits(summ, 'try-error')) {
    data.frame(parameter = parnames, mean = NA, se_mean = NA, sd = NA, q025 = NA, q25 = NA, q50 = NA, q75 = NA, q975 = NA, n_eff = NA, Rhat = NA)
  } else {
    data.frame(parameter = parnames, summ) %>%
      setNames(c('parameter', 'mean', 'se_mean', 'sd', 'q025', 'q25', 'q50', 'q75', 'q975', 'n_eff', 'Rhat'))
  }
}

all_pars1 <- all_fits_mod1 %>%
  group_by(siteID, year) %>%
  do(get_pars(.$fit[[1]], parnames1))

all_pars2 <- all_fits_mod2 %>%
  group_by(siteID, year) %>%
  do(get_pars(.$fit[[1]], parnames2))

all_pars3 <- all_fits_mod3 %>%
  group_by(siteID, year) %>%
  do(get_pars(.$fit[[1]], parnames3))

all_pars <- rbind(data.frame(model = 'H1', all_pars1),
                  data.frame(model = 'H2', all_pars2),
                  data.frame(model = 'H3', all_pars3))

# Get information criteria ------------------------------------------------

get_looic <- function(fit) {
  ll <- extract_log_lik(fit)
  looic <- loo(ll)
  data.frame(LOOIC = looic$estimates['looic', 'Estimate'],
             se_LOOIC = looic$estimates['looic', 'SE'])
}

all_loo1 <- all_fits_mod1 %>%
  group_by(siteID, year) %>%
  do(get_looic(.$fit[[1]]))

all_loo2 <- all_fits_mod2 %>%
  group_by(siteID, year) %>%
  do(get_looic(.$fit[[1]]))

all_loo3 <- all_fits_mod3 %>%
  group_by(siteID, year) %>%
  do(get_looic(.$fit[[1]]))

all_loo <- rbind(data.frame(model = 'H1', all_loo1),
                  data.frame(model = 'H2', all_loo2),
                  data.frame(model = 'H3', all_loo3))

# Get fitted values from functions and their CIs --------------------------

x_pred <- exp(seq(log(min(mammal2016$weight)), log(max(mammal2016$weight)), length.out = 51)) # Vector of x values to get fitted values

# Functions to return fitted value given function parameters
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

# Function to return fitted value quantiles from function fit
get_fit_quant <- function(fit, parnames, pdf, n, xmin) {
  # Extract parameters
  pars <- as.data.frame(do.call('cbind', extract(fit, parnames)))
  # Get the number of individuals and the minimum value from the original data
  # Calculate fitted values
  fitted <- n * pmap_dfc(pars, pdf, x = x_pred, xmin =xmin)
  # Get quantiles of fitted values
  qprobs <- c(0.025, 0.5, 0.975)
  fittedquant <- t(apply(fitted, 1, quantile, probs = qprobs))
  # Make fitted value df
  data.frame(x = x_pred, fittedquant) %>%
    setNames(c('x','q025','q50','q975'))
}

all_fitted1 <- all_fits_mod1 %>%
  group_by(siteID, year) %>%
  do(get_fit_quant(.$fit[[1]], parnames1, pdf_mod1, 
                   n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]]))

all_fitted2 <- all_fits_mod2 %>%
  group_by(siteID, year) %>%
  do(get_fit_quant(.$fit[[1]], parnames2, pdf_mod2, 
                   n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]]))

all_fitted3 <- all_fits_mod3 %>%
  group_by(siteID, year) %>%
  do(get_fit_quant(.$fit[[1]], parnames3, pdf_mod3, 
                   n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]]))

all_fitted <- rbind(data.frame(model = 'H1', all_fitted1),
                    data.frame(model = 'H2', all_fitted2),
                    data.frame(model = 'H3', all_fitted3))

# Function to return *energy* quantiles
# Must normalize
# Edit this 13 Sep 2018: now use the more accurate values for energy.
get_fitted_energy_quant <- function(fit, parnames, pdf, n, xmin, total_e) {
  require(pracma)
  # Extract parameters
  pars <- as.data.frame(do.call('cbind', extract(fit, parnames)))
  # Get the number of individuals and the minimum value from the original data
  # Calculate fitted values
  fitted <- n * 0.02446148 * pmap_dfc(pars, pdf, x = x_pred, xmin =xmin) ^ 0.6777703
  
  # Integrate fitted total energy and multiply total energy fitted values by the 
  # ratio of total observed energy and integral of fitted energy
  # (Use trapezoidal integration)
  fitted <- sapply(1:ncol(fitted), function(i) {
    fitted_integral <- trapz(x = x_pred, y = fitted[,i])
    fitted[,i] * total_e / fitted_integral
  })
  
  
  # Get quantiles of fitted values
  qprobs <- c(0.025, 0.5, 0.975)
  fittedquant <- t(apply(fitted, 1, quantile, probs = qprobs))
  # Make fitted value df
  data.frame(x = x_pred, fittedquant) %>%
    setNames(c('x','q025','q50','q975'))
}

all_fittedenergy1 <- all_fits_mod1 %>%
  group_by(siteID, year) %>%
  do(get_fitted_energy_quant(.$fit[[1]], parnames1, pdf_mod1, 
                   n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   total_e = total_energy$total[total_energy$siteID==.$siteID[1] & total_energy$year==.$year[1]]))

all_fittedenergy2 <- all_fits_mod2 %>%
  group_by(siteID, year) %>%
  do(get_fitted_energy_quant(.$fit[[1]], parnames2, pdf_mod2, 
                   n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   total_e = total_energy$total[total_energy$siteID==.$siteID[1] & total_energy$year==.$year[1]]))

all_fittedenergy3 <- all_fits_mod3 %>%
  group_by(siteID, year) %>%
  do(get_fitted_energy_quant(.$fit[[1]], parnames3, pdf_mod3, 
                   n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                   total_e = total_energy$total[total_energy$siteID==.$siteID[1] & total_energy$year==.$year[1]]))

all_fittedenergy <- rbind(data.frame(model = 'H1', all_fittedenergy1),
                    data.frame(model = 'H2', all_fittedenergy2),
                    data.frame(model = 'H3', all_fittedenergy3))


# Get fitted slopes and their CIs -----------------------------------------

get_fitted_slope_quant <- function(fit, parnames, pdf, n, xmin) {
  # Extract parameters
  pars <- as.data.frame(do.call('cbind', extract(fit, parnames)))
  # Get the number of individuals and the minimum value from the original data
  # Calculate fitted values
  fitted <- n * pmap_dfc(pars, pdf, x = x_pred, xmin = xmin)
  
  # Function to get log slope
  log_slope <- function(x, y) (x[-1]/y[-1]) * diff(y)/diff(x)
  
  fitted_slopes <- map_dfr(as.data.frame(fitted), ~ log_slope(x_pred, .))
  
  # Get quantiles of fitted values
  qprobs <- c(0.025, 0.5, 0.975)
  fittedslopequant <- t(apply(fitted_slopes, 1, quantile, probs = qprobs))
  # Make fitted value df
  data.frame(x = x_pred[-length(x_pred)] + diff(x_pred)/2, fittedslopequant) %>%
    setNames(c('x','q025','q50','q975'))
}

all_fittedslope1 <- all_fits_mod1 %>%
  group_by(siteID, year) %>%
  do(get_fitted_slope_quant(.$fit[[1]], parnames1, pdf_mod1, 
                             n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                             xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]]))

all_fittedslope2 <- all_fits_mod2 %>%
  group_by(siteID, year) %>%
  do(get_fitted_slope_quant(.$fit[[1]], parnames2, pdf_mod2, 
                             n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                             xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]]))

all_fittedslope3 <- all_fits_mod3 %>%
  group_by(siteID, year) %>%
  do(get_fitted_slope_quant(.$fit[[1]], parnames3, pdf_mod3, 
                             n = n_individuals$n[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]],
                             xmin = n_individuals$xmin[n_individuals$siteID==.$siteID[1] & n_individuals$year==.$year[1]]))

all_fittedslope <- rbind(data.frame(model = 'H1', all_fittedslope1),
                          data.frame(model = 'H2', all_fittedslope2),
                          data.frame(model = 'H3', all_fittedslope3))

# Save results ------------------------------------------------------------

save(all_pars, all_fitted, all_fittedenergy, all_fittedslope, all_loo, file = '~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/modelfitstats2016.RData')


# Determine best model at each site ---------------------------------------

load('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/modelfitstats2016.RData')

library(ggplot2)
library(dplyr)

# Decision rule for best model

# 1. Get rid of any models that have "significantly" worse information criteria
# 2. Out of the remaining models pick the simplest model (i.e. tiebreaker goes to less number of parameters)


# Select best model for each site by finding simplest model from the group of non overlapping models

get_best_model <- function(dat) {
  # remove a model from list of best models if its interval is entirely worse than any other model
  dat <- dat %>%
    arrange(model) %>%
    mutate(min = LOOIC - se_LOOIC, max = LOOIC + se_LOOIC)
  best_models <- 1:3
  for (i in 1:3) {
    if (any(dat$min[i] > dat$max[-i])) {
      best_models <- best_models[best_models != i]
    }
  }
  # return the lowest number from remaining list
  return(data.frame(best = min(best_models)))
}

all_best <- all_loo %>%
  group_by(siteID, year) %>%
  do(get_best_model(.))


# Plot of function fits ---------------------------------------------------


# Limit the plots to the range of the data at each site
all_fitted_toplot <- all_fitted %>%
  left_join(n_individuals) %>%
  filter(x >= (0.9*xmin) & x <= (xmax*1.1))

all_fittedenergy_toplot <- all_fittedenergy %>%
  left_join(n_individuals) %>%
  filter(x >= (0.9*xmin) & x <= (xmax*1.1))

best_fitted_toplot <- all_fitted_toplot %>%
  left_join(all_best) %>%
  filter(model == paste0('H',best))

best_fittedenergy_toplot <- all_fittedenergy_toplot %>%
  left_join(all_best) %>%
  filter(model == paste0('H',best))

model_names <- c('EE','Optimal size', 'Truncated')

th <- theme_bw() +
  theme(legend.position = 'none', strip.background = element_blank())

p_massbin <- ggplot(mass_bins %>% 
         filter(bin_value > 0) %>%
         left_join(all_best)) +
  facet_wrap(~ siteID) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  geom_line(aes(x = x, y = q025, color = model), data = best_fitted_toplot, linetype = 'dotted') +
  geom_line(aes(x = x, y = q975, color = model), data = best_fitted_toplot, linetype = 'dotted') +
  geom_line(aes(x = x, y = q50, color = model, group = model), data = best_fitted_toplot, size = 1) +
  geom_text(aes(x = 30, y = 1e-3, label = model_names[best])) +
  scale_x_log10(name = 'Mass (g)') + scale_y_log10(name = 'Density (individuals/g)') +
  th

p_energybin <- ggplot(energy_bins %>% 
         filter(bin_value > 0) %>%
         left_join(all_best)) +
  facet_wrap(~ siteID) +
  geom_point(aes(x = bin_midpoint, y = bin_value)) +
  geom_line(aes(x = x, y = q025, color = model), data = best_fittedenergy_toplot, linetype = 'dotted') +
  geom_line(aes(x = x, y = q975, color = model), data = best_fittedenergy_toplot, linetype = 'dotted') +
  geom_line(aes(x = x, y = q50, color = model, group = model), data = best_fittedenergy_toplot, size =  1) +
  geom_text(aes(x = 30, y = 0.001, label = model_names[best])) +
  scale_x_log10(name = 'Mass (g)') + scale_y_log10(name = 'Metabolic rate (W/g)') +
  th

# Plot of information criteria --------------------------------------------

all_loo <- all_loo %>%
  group_by(siteID, year) %>%
  mutate(deltaLOOIC = LOOIC - min(LOOIC))

p_looic <- ggplot(all_loo, aes(x=model, y=deltaLOOIC, ymin=deltaLOOIC-se_LOOIC, ymax=deltaLOOIC+se_LOOIC)) +
  facet_wrap(~ siteID) +
  geom_errorbar(width = 0.3) +
  geom_point() +
  scale_y_reverse(name = expression(Delta*LOOIC)) +
  scale_x_discrete(labels = c('EE','Optimal','Truncated EE')) +
  th


# Plot of slopes, given best model ----------------------------------------

# Fitted energy slopes

p_slopes <- all_fittedslope %>%
  left_join(all_best) %>%
  filter(paste0('H',best)==model) %>%
  mutate_at(vars(starts_with('q')), funs(.+beta1)) %>%
  ggplot(aes(x = x, y = q50, ymin = q025, ymax = q975)) +
    facet_wrap(~ siteID) +
    geom_hline(yintercept = 0, color = 'red', size = 0.7) +
    geom_line(aes(y = q025), linetype = 'dotted') +
    geom_line(aes(y = q975), linetype = 'dotted') +
    geom_line(size = 1) +
    scale_x_log10(name = 'Mass (g)') + scale_y_continuous(name = 'Fitted energy slope') + 
    th


# Save pdfs ---------------------------------------------------------------

pdf('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Analyses/threemodelfits.pdf', height = 9, width = 9)
  p_looic + ggtitle('Information criteria')
  p_massbin + ggtitle('Fits to mass')
  p_energybin + ggtitle('Fits to energy')
  p_slopes + ggtitle('Fitted slopes')
dev.off()
