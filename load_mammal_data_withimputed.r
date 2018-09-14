# Load mammal data and taxonomy
# This version also loads the imputed traits
data_path <- '~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data'

mammal_data <- read.csv(file.path(data_path, 'mammal_reduced.csv'), stringsAsFactors = FALSE)
mammal_tax <- read.csv(file.path(data_path, 'small_mammal_taxonomy.csv'), stringsAsFactors = FALSE) %>%
  select(acceptedTaxonID, taxonProtocolCategory, dwc.scientificName) %>%
  rename(taxonID = acceptedTaxonID)
met_intercepts <- read.csv(file.path(data_path, 'imputed_intercepts.csv'), stringsAsFactors = FALSE)

# Join mammal data with taxonomy and get rid of the unneeded rows
mammal_data <- mammal_data %>%
  filter(!recapture %in% c('U','Y')) %>% # 
  mutate(month = month(collectDate), year = year(collectDate)) %>%
  filter(!is.na(weight)) %>%
  select(domainID, siteID, plotID, collectDate, month, year, taxonID, lifeStage, weight) %>%
  left_join(mammal_tax) %>%
  filter(taxonProtocolCategory %in% 'target')

# Join with metabolic rate intercepts and get the metabolic rate for each individual
beta1 <- 0.6777703

mammal_data <- mammal_data %>%
  left_join(met_intercepts[,c('taxonID', 'beta0')]) %>%
  mutate(metabolicRate = beta0 * weight ^ beta1)

range(mammal_data$weight)
# One is zero so let's get rid of it.
mammal_data <- filter(mammal_data, weight > 0)
(weight_range <- range(mammal_data$weight))
