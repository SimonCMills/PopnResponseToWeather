## Script to create final df of bird abundance data to use in analyses. 
## 
## Merging in contextual information (lat&lon, site climateID, and species names)
## from datakeys & add new cols (relative growth rate and log-count in year prior)
##
## Formatting/cleaning steps comprise: 
## (1) Remove canary islands 
## (2) Remove obs w/o prior obs or where both prior and current obs ==0. After 
## this remove sites with 5 or fewer obs
## (3) remove feral or ambiguous species
## (4) Finally, drop sites w/o climate data
## (5) Remove redundant columns & reorder columns to produce final cleaned df

# housekeeping ----
rm(list=ls()); graphics.off()
require(lubridate); require(ggplot2); require(dplyr)

# helper functions
# function to get count in year prior (checks that year interval ==1)
source("scripts/functions/count_prior.R")

# files ----
# read in climate key 
climateKey <- readRDS("./files/climateKey_siteID-climateID.rds") %>% unique

# Read in bird dataset, species data key, and site data key
bird <- readRDS("../../PECBMS/cleanedFiles/countData_allSchemes_extraPECBMSfiles.rds") 
birdKey <- readRDS("../../PECBMS/cleanedFiles/dataKey_Euring-species.rds")
siteKey <- readRDS("../../PECBMS/cleanedFiles/dataKey_siteID-WGS84coord.rds") %>%
    as_data_frame %>%
    unique

# join in species names and lats and lons
bird_full <- left_join(bird, birdKey) %>%
    left_join(., siteKey)

# add climate ID
bird_clim <- left_join(bird_full, climateKey)

# check: #observations should be unchanged 
nrow(bird_clim) == nrow(bird) # TRUE

## formatting steps ----
# first create column with count in year prior
bird_formatted <- bird_clim %>%
    group_by(species, siteID) %>% 
    arrange(species, siteID, year) %>%
    mutate(diff_priorYr = year - c(NA, year[-length(year)]), 
           count_pre = count_prior(count_t, diff_priorYr, 1)) %>%
    ungroup

# (1) Remove canary islands 
bird_clean1 <- filter(bird_formatted, lon_WGS84 > -10)

# (2) Remove obs w/o prior obs or where both prior and current obs ==0. After 
# this remove sites with 5 or fewer obs (i.e. screening subsequent to first step)
bird_clean2 <- bird_clean1 %>%
    ungroup %>%
    # removing NA and non-1 year differences
    filter(!is.na(diff_priorYr)) %>%
    filter(diff_priorYr==1) %>%
    # removing observations where both counts are == 0
    filter(!(count_pre == 0 & count_t == 0)) %>%
    # finally, remove sites that have fewer than 5 observations for a species &
    # add relative growth rate column
    group_by(species, siteID) %>%
    mutate(n_site = n()) %>%
    ungroup %>%
    filter(!n_site <= 5)

# (3) remove feral or ambiguous species
sp_toRemove <- c("Feral Pigeon", 
                 "Rock Dove", 
                 "Egyptian Goose", 
                 "Ring-necked Parakeet", 
                 "Canada Goose", 
                 "unidentified crossbill", 
                 "Carrion x Hooded Crow")

bird_clean3 <- bird_clean2 %>% 
    ungroup %>%
    filter(!species %in% sp_toRemove)

# (4) Finally, drop sites w/o climate data
# read in monthly climate data (note, cells that are not present here are those 
# that, by definition, have no climate data)
climate_summ <- readRDS("files/climateDF_monthlySummary.rds")
hasClimate <- climate_summ %>%
    .$climateID %>%
    unique

bird_clean4 <- bird_clean3 %>% 
    filter(climateID %in% hasClimate)

# (5) Remove redundant columns & reorder columns to produce final cleaned df
bird_final <- bird_clean4 %>%
    mutate(rgr = log(count_t+1) - log(count_pre+1), 
           logCount_pre = log(count_pre + 1)) %>%
    select(species, species_LB, Euring, schemeID, climateID, siteID, lat_WGS84, 
           lon_WGS84, n_site, year, count_t, count_pre, rgr, logCount_pre)

# summarise # datapoints by species ----
n_bySpecies <- bird_final %>% 
    select(species, siteID, n_site) %>%
    unique %>%
    group_by(species) %>%
    summarise(n_species = sum(n_site)) %>%
    arrange(desc(n_species)) %>%
    mutate(species = paste0(1:n(), ". ", species)) %>%
    arrange(n_species) %>%
    mutate(species = factor(species, levels=species))

n_bySpecies %>%
    filter(n_species > 20000) %>%
    ggplot(aes(x=n_species, y=species)) + geom_point() +
    labs(x = "Number of datapoints", y = "Species", 
         caption = bquote(bold("Figure: ")*"Number of datapoints following data-cleaning steps"))

ggsave("figures/nDP_bySpecies_cleanedDataset.png", width=210, height=297, dpi = 150, units="mm")

## save cleaned datasets----
saveRDS(bird_final, "files/birdDF.rds")
saveRDS(n_bySpecies, "files/nDP in final dataset by species.rds")
