# Load mammal data and taxonomy
mammal_data <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/mammal_reduced.csv', stringsAsFactors = FALSE)
mammal_tax <- read.csv('~/google_drive/NEON_EAGER/Manuscript5_RodentEE/Data/small_mammal_taxonomy.csv', stringsAsFactors = FALSE) %>%
  select(acceptedTaxonID, taxonProtocolCategory, dwc.scientificName) %>%
  rename(taxonID = acceptedTaxonID)

mammal_data <- mammal_data %>%
  filter(!recapture %in% c('U','Y')) %>% # 
  mutate(month = month(collectDate), year = year(collectDate)) %>%
  filter(!is.na(weight)) %>%
  select(domainID, siteID, plotID, collectDate, month, year, taxonID, lifeStage, weight) %>%
  left_join(mammal_tax) %>%
  filter(taxonProtocolCategory %in% 'target')

mammal_data <- mammal_data %>%
  mutate(relativeEnergy = weight ^ 0.75)

range(mammal_data$weight)
# One is zero so let's get rid of it.
mammal_data <- filter(mammal_data, weight > 0)
(weight_range <- range(mammal_data$weight))