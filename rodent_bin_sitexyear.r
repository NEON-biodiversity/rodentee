# Binning by mass and by energy
# By site and by year


# Define binning function -------------------------------------------------

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


# Load and process data ---------------------------------------------------------------


# This CSV was pulled from NEON in April 2018. I deleted some unneeded columns and rows denoting empty traps so that the file is smaller but other than that it isn't modified at all from what I got directly from NEON.

mammal_data <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/mammal_reduced.csv', stringsAsFactors = FALSE)

# Get rid of all recaptures. 
# Other than that do no other quality control for now. Probably a good idea to look back through later and get rid of odd outliers. 

library(dplyr)
library(lubridate)

mammal_data <- mammal_data %>%
  filter(!recapture %in% c('U','Y'), weight > 0) %>%
  mutate(month = month(collectDate), year = year(collectDate)) %>%
  select(domainID, siteID, plotID, collectDate, month, year, taxonID, lifeStage, weight)

# Create a rough guess for relative metabolic rate (energy flux) by taking mass^.75 for each individual.

mammal_data <- mammal_data %>%
  mutate(relativeEnergy = weight ^ 0.75)


# Do binning --------------------------------------------------------------

# Set global bin bounds based on range

(weight_range <- range(mammal_data$weight))

n_bins <- 8

# If we did 8 bins these would be the bin bounds.
(bin_bounds <- exp(seq(log(weight_range[1]), log(weight_range[2]), length.out = n_bins + 1)))

# Use the log-binning algorithm to bin each site's body mass measurements based on the bin bounds we just determined. 
# We don't care about the date of sampling or the plot within site that's sampled. They are all lumped.

mass_bins <- mammal_data %>%
  group_by(siteID, year) %>%
  do(logbin_setedges(x = .$weight, bin_edges = bin_bounds))

# Use the same bins to get the binned energy flux. 

energy_bins <- mammal_data %>%
  group_by(siteID, year) %>%
  do(logbin_setedges(x = .$weight, y = .$relativeEnergy, bin_edges = bin_bounds))

# Make a plot, by site, of binned body mass. I only plotted the ones with at least 100 individuals for now.
# Color code the year to get darker closer to the present. Can be easily modified.

library(ggplot2)

n_individuals <- mass_bins %>%
  group_by(siteID, year) %>%
  summarize(n = sum(bin_count)) %>%
  filter(n >= 100)

mass_bins %>% 
  ungroup %>%
  filter(bin_value > 0, siteID %in% n_individuals$siteID) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = bin_midpoint, y = bin_value, color = year, group = year)) + 
  facet_wrap(~ siteID) +
  scale_color_manual(values = rev(gray.colors(6))) +
  theme_bw() + geom_line() + geom_point() +
  scale_x_log10(name = 'Mass (g)') + scale_y_log10('Density (individuals/g)')

# Make a plot, by site, of binned (relative) total energy flux (assuming each individual has a metabolic rate exactly proportional to 3/4 power of its mass). 
# Again only sites with at least 100 individuals are shown. 

energy_bins %>% 
  ungroup %>%
  filter(bin_value > 0, siteID %in% n_individuals$siteID) %>%
  mutate(year = factor(year)) %>%
  ggplot(aes(x = bin_midpoint, y = bin_value, color = year, group = year)) + 
  facet_wrap(~ siteID) +
  scale_color_manual(values = rev(gray.colors(6))) +
  theme_bw() + geom_line() + geom_point() +
  scale_x_log10(name = 'Mass (g)') + scale_y_log10(name = 'Total relative energy')
