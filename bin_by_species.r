# Bin mass and energy by species
# QDR / NEON Rodents / 04 Sep 2018

library(dplyr)
library(lubridate)

# Load data
mammal_data <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/mammal_reduced.csv', stringsAsFactors = FALSE)

# Get rid of recaptures, mammals where no mass was recorded, extract month and year, and get rid of unnecessary columns
mammal_data <- mammal_data %>%
  filter(!recapture %in% c('U','Y')) %>% # 
  mutate(month = month(collectDate), year = year(collectDate)) %>%
  filter(!is.na(weight), weight > 0) %>%
  select(domainID, siteID, plotID, collectDate, month, year, taxonID, lifeStage, weight)

# Energy defined as proportional to 3/4 power of mass
mammal_data <- mammal_data %>%
  mutate(relativeEnergy = weight ^ 0.75)

# Get average mass of each species, across all sites/years
# (Do geometric mean, though this can be changed)

sp_avgs <- mammal_data %>%
  group_by(taxonID) %>%
  summarize(global_sp_avg_mass = exp(mean(log(weight))), global_sp_avg_energy = global_sp_avg_mass ^ 0.75)

# Get the range of average masses
(avg_mass_range <- range(sp_avgs$global_sp_avg_mass))

# Define 8 bins
n_bins <- 8

# Bin boundaries (9 values so there are 8 intervals)
(bin_bounds <- exp(seq(log(avg_mass_range[1]), log(avg_mass_range[2]), length.out = n_bins + 1)))

# Assign each species to a bin, globally defined.
# So regardless of what individuals appear of that species at a site-year combo, it will always be in the same bin

sp_avgs <- sp_avgs %>%
  mutate(bin = findInterval(global_sp_avg_mass, bin_bounds[-length(bin_bounds)]))

table(sp_avgs$bin) # Richness per log bin

# Get individual per-species mass and energy for each site-year combo
# Join it with the global numbers.
mammal_sp_sums <- mammal_data %>%
  group_by(siteID, year, taxonID) %>%
  summarize(site_sp_total_mass = sum(weight), site_sp_total_relativeEnergy = sum(relativeEnergy)) %>%
  left_join(sp_avgs)

# Get per bin averages
# Again use geometric means
mammal_bin_avgs <- mammal_sp_sums %>%
  ungroup %>%
  group_by(siteID, year, bin) %>%
  summarize(bin_sp_avg_mass = exp(mean(log(site_sp_total_mass))),
            bin_sp_avg_relativeEnergy = exp(mean(log(site_sp_total_relativeEnergy))))

# Join the species averages with the actual bin midpoints so they can be plotted
(bin_midpoints <- bin_bounds[-length(bin_bounds)] + diff(bin_bounds)/2)

mammal_bin_avgs <- mammal_bin_avgs %>%
  mutate(bin_midpoint = bin_midpoints[bin])


# Create a plot of the results

library(ggplot2)

# Figure out which ones have 100 or more individuals and plot only those
n_individuals <- mammal_data %>%
  group_by(siteID, year) %>%
  summarize(n = n())

mammal_bin_avgs %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 100) %>%
ggplot(aes(x = bin_midpoint, y = bin_sp_avg_relativeEnergy)) +
  facet_wrap(~ siteID) +
  geom_point() + 
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  ggtitle('Average energy of a species in each log bin')

mammal_sp_sums %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 100) %>%
ggplot(aes(x = global_sp_avg_mass, y = site_sp_total_relativeEnergy)) +
  facet_wrap(~ siteID) +
  geom_point() + 
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  ggtitle('Total energy of each species, disregarding bins')
