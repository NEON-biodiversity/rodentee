# Code to reduce size of small mammal data file to be used for new NEON rodent analysis
# QDR 11 May 2018

smt <- read.csv('/mnt/research/neon/raw_data/organismal_data_apr2018/small_mammal_trapping.csv')

library(dplyr)

# Get rid of all the rows denoting empty traps
# Get rid of a lot of ID columns we don't care about

keep_cols <- c("uid", "nightuid", "namedLocation", "domainID", "siteID", "plotID",
"trapCoordinate", "plotType", "nlcdClass", "decimalLatitude",
"decimalLongitude", "elevation", "trapStatus", "trapType", "collectDate",
"tagID", "taxonID", "scientificName", "taxonRank", "nativeStatusCode",
"sex", "recapture", "fate", "replacedTag", "lifeStage", "testes",
"nipples", "pregnancyStatus", "vagina", "hindfootLength", "earLength",
"tailLength", "totalLength", "weight", "larvalTicksAttached",
"nymphalTicksAttached", "adultTicksAttached")


smt <- smt %>%
	filter(grepl('capture', trapStatus)) %>%
	select(keep_cols)
	
write.csv(smt, file = '/mnt/research/neon/MS5_RodentEE/data/mammal_reduced.csv', row.names = FALSE)