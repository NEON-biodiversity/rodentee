# Bin mass and energy by species
# QDR / NEON Rodents / 04 Sep 2018

# Modified 13 Sep: use new metabolic rate

library(dplyr)
library(lubridate)

source('load_mammal_data_withimputed.r')

# Get average mass of each species, across all sites/years
# (Do geometric mean, though this can be changed)

sp_avgs <- mammal_data %>%
  group_by(taxonID) %>%
  summarize(global_sp_avg_mass = exp(mean(log(weight))), global_sp_avg_metabolicRate = exp(mean(log(metabolicRate))))

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
  summarize(site_sp_total_mass = sum(weight), site_sp_total_metabolicRate = sum(metabolicRate)) %>%
  left_join(sp_avgs)

# Get per bin averages
# Again use geometric means
mammal_bin_avgs <- mammal_sp_sums %>%
  ungroup %>%
  group_by(siteID, year, bin) %>%
  summarize(bin_sp_avg_mass = exp(mean(log(site_sp_total_mass))),
            bin_sp_avg_metabolicRate = exp(mean(log(site_sp_total_metabolicRate))),
            bin_richness = length(unique(taxonID)))

# Join the species averages with the actual bin midpoints so they can be plotted
(bin_midpoints <- bin_bounds[-length(bin_bounds)] + diff(bin_bounds)/2)

mammal_bin_avgs <- mammal_bin_avgs %>%
  mutate(bin_midpoint = bin_midpoints[bin])


# Create a plot of the results

library(ggplot2)

# Figure out which ones have 50 or more individuals, and 4 or more species and plot only those
n_individuals <- mammal_data %>%
  group_by(siteID, year) %>%
  summarize(n = n(), site_richness = length(unique(taxonID)))

p_by_bin <- mammal_bin_avgs %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 50, site_richness >= 4) %>%
ggplot(aes(x = bin_midpoint, y = bin_sp_avg_metabolicRate)) +
  facet_wrap(~ siteID) +
  geom_point() + 
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(strip.background = element_blank())
  ggtitle('Average total energy of a species in each log bin')

p_by_sp <- mammal_sp_sums %>%
  left_join(n_individuals) %>%
  filter(year == 2016, n >= 50, site_richness >= 4) %>%
ggplot(aes(x = global_sp_avg_mass, y = site_sp_total_metabolicRate)) +
  facet_wrap(~ siteID) +
  geom_point() + 
  scale_x_log10() + scale_y_log10() +
  theme_bw() + theme(strip.background = element_blank())
ggtitle('Total energy of each species, disregarding bins')

ggsave('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Analyses/energy_by_species_by_bin.pdf', p_by_bin, height = 8, width = 8)
ggsave('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Analyses/energy_by_species.pdf', p_by_sp, height = 8, width = 8)
