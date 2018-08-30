# Using all mammals from 2016, fit models corresponding to the three hypotheses.
# Use information criteria to differentiate between them.

library(rstan)
library(bayesplot)
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
# Test this on a single site, Konza.

x <- subset(mammal_data, siteID == 'KONZ')$weight

standata <- list(x = x, N = length(x), x_min = min(x), x_max = max(x))

stanmod1 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/pareto.stan')
stanmod2 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/triang_continuous.stan') 
stanmod3 <- stan_model('~/Documents/GitHub/NEON_repos/rodentee/threepart_continuous.stan') 

n_iter <- 5000
n_warm <- 4000
n_chain <- 3

fit1 <- sampling(stanmod1, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 11111)
fit2 <- sampling(stanmod2, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 22222)
fit3 <- sampling(stanmod3, data = standata, chains = n_chain, iter = n_iter, warmup = n_warm, seed = 33333)

summary(fit1)$summary[1,]
summary(fit2)$summary[1:5,]
summary(fit3)$summary[1:8,]

mcmc_trace(as.array(fit1))
mcmc_trace(as.array(fit2))
mcmc_trace(as.array(fit3))
