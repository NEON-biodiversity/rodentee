# Fit distribution to average species mass and total species energy

# Run bin_by_species.r to get mammal_sp_sums object

# Fit peaked density, linear production

library(rstan)
library(bayesplot)
library(loo)
library(purrr)
options(mc.cores = 3)
rstan_options(auto_write = TRUE)

model_sp_h1 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_bysp_h1.stan')
model_sp_h2 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_bysp_h2.stan')
model_sp_h3 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/model_bysp_h3.stan')

# Define function to fit this for each siteID.
fit_by_species <- function(dat, mod) {
  standata <- with(dat, list(N = length(global_sp_avg_mass),
                             x = global_sp_avg_mass,
                             y = site_sp_total_metabolicRate,
                             x_min = min(global_sp_avg_mass),
                             x_max = max(global_sp_avg_mass)))
  sampling(mod, data = standata, iter = 5000, warmup = 4000, chains = 3)
  
}

# filter out any site with less than 4 species or 50 indiv, and only do 2016.
# This only takes a few minutes to run for all sites.
n_individuals <- mammal_data %>%
  group_by(siteID, year) %>%
  summarize(n=n(), site_richness = length(unique(taxonID)), xmin = min(weight), xmax = max(weight))

all_fits_m1 <- mammal_sp_sums %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 50, site_richness >=4) %>%
  do(fit = fit_by_species(., model_sp_h1))

all_fits_m2 <- mammal_sp_sums %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 50, site_richness >=4) %>%
  do(fit = fit_by_species(., model_sp_h2))

all_fits_m3 <- mammal_sp_sums %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 50, site_richness >=4) %>%
  do(fit = fit_by_species(., model_sp_h3))

# Extract the fitted values and their credible intervals, for plotting.
source('~/Documents/GitHub/forestlight/stan/piecewise_workflow/extraction_functions_piecewise.r')

parnames <- list(c('alpha', 'beta0', 'beta1'),
                 c('alpha_low', 'alpha_high', 'beta0', 'beta1', 'tau'),
                 c('alpha_low', 'alpha_mid', 'alpha_high', 'beta0', 'beta1', 'tau_low', 'tau_high'))

all_pars_m1 <- all_fits_m1 %>%
  group_by(siteID) %>%
  do(param_values(.$fit[[1]], parnames[[1]]))


all_pars_m2 <- all_fits_m2 %>%
  group_by(siteID) %>%
  do(param_values(.$fit[[1]], parnames[[2]]))


all_pars_m3 <- all_fits_m3 %>%
  group_by(siteID) %>%
  do(param_values(.$fit[[1]], parnames[[3]]))

# Extract log likelihood values and use these to get information criteria

get_ic <- function(fit) {
  ll_dens <- extract_log_lik(fit, 'log_lik_dens')
  waic_dens <- waic(ll_dens)
  loo_dens <- loo(ll_dens)
  data.frame(WAIC = waic_dens$estimates['waic','Estimate'],
             se_WAIC = waic_dens$estimates['waic','SE'],
             LOOIC = loo_dens$estimates['looic','Estimate'],
             se_LOOIC = loo_dens$estimates['looic','SE'])
}

all_ic_m1 <- all_fits_m1 %>%
  group_by(siteID) %>%
  do(get_ic(.$fit[[1]]))

all_ic_m2 <- all_fits_m2 %>%
  group_by(siteID) %>%
  do(get_ic(.$fit[[1]]))

all_ic_m3 <- all_fits_m3 %>%
  group_by(siteID) %>%
  do(get_ic(.$fit[[1]]))

all_ic <- cbind(model = rep(1:3, each = nrow(all_fits_m1)), rbind(all_ic_m1, all_ic_m2, all_ic_m3))

# Get best model
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

all_best <- all_ic %>%
  group_by(siteID) %>%
  do(get_best_model(.))

# Get mass values where we want to find the fitted values
mass_limits <- range(mammal_sp_sums$global_sp_avg_mass)
x_pred <- exp(seq(log(min(mass_limits)), log(max(mass_limits)), length.out = 51))

# Get minimum mass value and total energy value from each site for normalization of total energy curve
min_n_2016 <- mammal_sp_sums %>%
  filter(year == 2016) %>%
  group_by(siteID) %>%
  summarize(x_min = min(global_sp_avg_mass), total_energy = sum(site_sp_total_metabolicRate), site_richness = n())
  
all_fitted_m1 <- all_fits_m1 %>%
  left_join(min_n_2016) %>%
  group_by(siteID) %>%
  do(fitted_predicted_values(.$fit[[1]], x_pred, dens_form = 1, prod_form = 1, x_min = .$x_min[1], n_indiv = .$site_richness[1], total_prod = .$total_energy[1], pars_to_get = parnames[[1]]))

all_fitted_m2 <- all_fits_m2 %>%
  left_join(min_n_2016) %>%
  group_by(siteID) %>%
  do(fitted_predicted_values(.$fit[[1]], x_pred, dens_form = 2, prod_form = 1, x_min = .$x_min[1], n_indiv = .$site_richness[1], total_prod = .$total_energy[1], pars_to_get = parnames[[2]]))

all_fitted_m3 <- all_fits_m3 %>%
  left_join(min_n_2016) %>%
  group_by(siteID) %>%
  do(fitted_predicted_values(.$fit[[1]], x_pred, dens_form = 3, prod_form = 1, x_min = .$x_min[1], n_indiv = .$site_richness[1], total_prod = .$total_energy[1], pars_to_get = parnames[[3]]))


# Edit fitted value data frame to change the names from the BCI names to mammal names
names(all_fitted_m1)[2] <- 'mass'
names(all_fitted_m2)[2] <- 'mass'
names(all_fitted_m3)[2] <- 'mass'

all_fitted_m1 <- all_fitted_m1 %>%
  mutate(variable = factor(variable, labels = c('richness per size', 'energy per species', 'energy per species fitted', 'energy per size', 'energy per size fitted'))) %>%
  filter(variable %in% c('richness per size', 'energy per species', 'energy per size'))

all_fitted_m2 <- all_fitted_m2 %>%
  mutate(variable = factor(variable, labels = c('richness per size', 'energy per species', 'energy per species fitted', 'energy per size', 'energy per size fitted'))) %>%
  filter(variable %in% c('richness per size', 'energy per species', 'energy per size'))

all_fitted_m3 <- all_fitted_m3 %>%
  mutate(variable = factor(variable, labels = c('richness per size', 'energy per species', 'energy per species fitted', 'energy per size', 'energy per size fitted'))) %>%
  filter(variable %in% c('richness per size', 'energy per species', 'energy per size'))


# Create spiffy plot.

# Reorder factors
all_fitted_m1 <- all_fitted_m1 %>%
  mutate(variable = factor(variable, levels = c('energy per species', 'richness per size', 'energy per size')))
all_fitted_m2 <- all_fitted_m2 %>%
  mutate(variable = factor(variable, levels = c('energy per species', 'richness per size', 'energy per size')))
all_fitted_m3 <- all_fitted_m3 %>%
  mutate(variable = factor(variable, levels = c('energy per species', 'richness per size', 'energy per size')))

# Just show first three sites
ggplot(all_fitted_m1 %>% filter(siteID %in% c('BART','BLAN','CLBJ')), aes(x = mass)) +
  facet_grid(siteID ~ variable) +
  geom_line(aes(y = q50), size = 1) +
  geom_line(aes(y = q025), linetype = 'dotted') +
  geom_line(aes(y = q975), linetype = 'dotted') +
  scale_x_log10(name = 'mass (g)') + scale_y_log10(name = 'fitted value') +
  theme_bw() +
  theme(strip.background = element_blank())
  
# Create logbins of richness per mass and energy per mass

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

n_bins <- 8
(bin_bounds <- exp(seq(log(mass_limits[1]), log(mass_limits[2]), length.out = n_bins + 1)))

energy_logbins <- mammal_sp_sums %>%
  do(logbin_setedges(x = .$global_sp_avg_mass, y = .$global_sp_avg_energy, bin_edges = bin_bounds))

richness_logbins <- mammal_sp_sums %>%
  do(logbin_setedges(x = .$global_sp_avg_mass, y = rep(1, nrow(.)), bin_edges = bin_bounds))

ggplot(all_fitted_m3 %>% filter(siteID %in% c('BART','BLAN','CLBJ')), aes(x = mass)) +
  facet_grid(siteID ~ variable) +
  geom_point(data = energy_logbins %>% 
               mutate(variable = 'energy per size') %>%
               filter(year == 2016, bin_value > 0, siteID %in% c('BART','BLAN','CLBJ')), 
             aes(x = bin_midpoint, y = bin_value)) +
  geom_point(data = mammal_sp_sums %>%
               mutate(variable = 'energy per species') %>%
               filter(year == 2016, siteID %in% c('BART','BLAN','CLBJ')),
             aes(x = global_sp_avg_mass, y = site_sp_total_metabolicRate)) +
  geom_point(data = richness_logbins %>%
               mutate(variable = 'richness per size') %>%
               filter(year == 2016, bin_value > 0, siteID %in% c('BART','BLAN','CLBJ')),
             aes(x = bin_midpoint, y = bin_value)) +
  geom_line(aes(y = q50), size = 1, color = 'slateblue3') +
  geom_line(aes(y = q025), linetype = 'dotted', color = 'slateblue3') +
  geom_line(aes(y = q975), linetype = 'dotted', color = 'slateblue3') +
  scale_x_log10(name = 'mass (g)') + scale_y_log10(name = 'fitted value') +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank())

# Plot everything across multiple pages

library(ggforce)
n_pages <- ceiling(nrow(all_fits) / 3)

plot_pages <- lapply(1:n_pages, function(i)
ggplot(all_fitted, aes(x = mass)) +
  facet_grid_paginate(siteID ~ variable, nrow = 3, ncol = 3, page = i, scales = 'free_y') +
  geom_point(data = energy_logbins %>% 
               mutate(variable = factor('energy per size')) %>%
               filter(year == 2016, bin_value > 0, siteID %in% all_fits$siteID), 
             aes(x = bin_midpoint, y = bin_value)) +
  geom_point(data = mammal_sp_sums %>%
               mutate(variable = factor('energy per species')) %>%
               filter(year == 2016, siteID %in% all_fits$siteID),
             aes(x = global_sp_avg_mass, y = site_sp_total_metabolicRate)) +
  geom_point(data = richness_logbins %>%
               mutate(variable = factor('richness per size')) %>%
               filter(year == 2016, bin_value > 0, siteID %in% all_fits$siteID),
             aes(x = bin_midpoint, y = bin_value)) +
  geom_line(aes(y = q50), size = 1, color = 'slateblue3') +
  geom_line(aes(y = q025), linetype = 'dotted', color = 'slateblue3') +
  geom_line(aes(y = q975), linetype = 'dotted', color = 'slateblue3') +
  scale_x_log10(name = 'mass (g)') + scale_y_log10(name = 'fitted value') +
  theme_bw() +
  theme(strip.background = element_blank(), panel.grid = element_blank())
)

pdf('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Analyses/allfitsbyspecies.pdf', height = 8, width = 8)
  for (i in plot_pages) print(i)
dev.off()
