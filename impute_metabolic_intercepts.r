# Imputation of rodent metabolic rates

library(dplyr)

source('load_mammal_data_tax.r')

# Clean up data by getting rid of commas in numbers
met_rate <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/rodent_metabolic_rates.csv', stringsAsFactors = FALSE) %>%
  mutate(mass_g = as.numeric(gsub(',', '', mass_g)))

# Find slope and intercept from data
lm_met <- lm(log(met_rate_W) ~ log(mass_g), data = met_rate)

global_slope <- lm_met$coefficients[2]
global_intercept <- exp(lm_met$coefficients[1])

# Assuming beta1 = 0.75, impute intercept beta0 for each species
# metabolic rate = beta0 * mass ^ beta1
# beta0 = metabolic rate * mass ^ -0.737
met_rate <- met_rate %>%
  mutate(beta0 = met_rate_W * mass_g ^ -global_slope)


# See how many are in the mammal data for the taxa we actually need
neon_mammal <- unique(mammal_data[,c('taxonID','dwc.scientificName')]) # 73 spp

matches <- neon_mammal$dwc.scientificName %in% met_rate$scientificName
table(matches) # 37 of 73 match

# Load phylogeny for imputation
library(Rphylopars)
mam_tree <- read.tree('~/google_drive/NEON_EAGER/Manuscript1_RodentOverlap/R_data/mammaltol.nwk')

neon_mammal <- neon_mammal %>%
  rename(scientificName = dwc.scientificName) %>%
  left_join(met_rate)

# Get masses for the ones that don't have a mass in the dataset
# Use median of observed adult mass
neon_masses <- mammal_data %>%
  filter(lifeStage == 'adult') %>%
  group_by(taxonID) %>%
  summarize(observed_mass = quantile(weight, 0.5))

neon_mammal <- neon_mammal %>%
  left_join(neon_masses) %>%
  mutate(mass_g = if_else(is.na(mass_g), observed_mass, mass_g))

# Use imputation to get the missing intercepts
# All taxa only identifiable at the genus level have already been added at the root of each genus
neon_mammal$taxonID %in% mam_tree$tip.label

imputation_traits <- neon_mammal %>%
  filter(taxonID %in% mam_tree$tip.label) %>%
  select(taxonID, mass_g, met_rate_W, beta0) %>%
  rename(species = taxonID) %>%
  mutate(mass_g = log(mass_g)) # Log transform mass vars

imputation_IDs <- neon_mammal %>%
  filter(taxonID %in% mam_tree$tip.label) %>%
  select(taxonID, scientificName)

set.seed(555)
phyimp <- phylopars(trait_data = imputation_traits,
                    tree = mam_tree,
                    model = 'OU')

# Put imputed and non imputed traits back together, and log-transform
# Add back the rows from spp not in the phylogeny
traits_imputed <- imputation_IDs %>%
  cbind(phyimp$anc_recon[1:nrow(imputation_traits), ]) %>%
  mutate(mass_g = exp(mass_g)) %>%
  bind_rows(neon_mammal %>% filter(!taxonID %in% mam_tree$tip.label) %>% select(-Order, -Family, -observed_mass))

hist(traits_imputed$beta0)

# There are three remaining species without an intercept, and then three taxa not identified to genus. Give them the overall number.
traits_imputed$beta0[is.na(traits_imputed$beta0)] <- global_intercept

write.csv(traits_imputed, '~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/imputed_intercepts.csv', row.names = FALSE)
