# Script to extract climate data & match it with site ID
#
# There are three main components: 
# 
# First match bird site latitude and longitudes with their corresponding climate
# cell. Achieve this by converting to an equal area coordinate system to flatten
# Europe and then finding closest climate cell centre with pythagoras. Generates
# a key that links each siteID to a corresponding climateID and array position 
# (lon_ID, lat_ID)
#
# extract all climate info for each climate ID for which there are bird data
# and discard rest. Check these patterns make sense (do this for temp for which
# there is more obvious expectation) & save

## Housekeeping ----
# rm(list=ls()); graphics.off()
library(ncdf4); library(dplyr); library(ggplot2); library(sp); library(lubridate)
CRS_WGS84 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
CRS_LAEA <- "+proj=laea +lat_0=90 +lon_0=10 +ellps=WGS84 +datum=WGS84 +units=m"

# Read in NetCDF files 
## Create key for climate array ----
# overall dimensions of these are the same, just contain different climate info 
ncdf_rain <- nc_open("files/data-raw/climate-netCDF/rr_0.22deg_rot_v12.0.nc") 
ncdf_temp <- nc_open("files/data-raw/climate-netCDF/tg_0.22deg_rot_v12.0.nc")

# get lats and lons and assign a climate ID
lat <- ncvar_get(ncdf_rain,"Actual_latitude")
lon <- ncvar_get(ncdf_rain,"Actual_longitude")
climateKey <- data_frame(lat_WGS84 = c(lat), lon_WGS84=c(lon)) %>%
    mutate(climateID = 1:n())

## The explicit approach to getting array indices for each climate ID (using which()
## to find the lats and lons that correspond to a particular array index) is slow, 
## so replaced with code that just replicates the logic. 
## Can confirm correctness of result using the identical(). 
# climID <- climateKey %>% 
#     group_by(climateID) %>%
#     mutate(rowID = which(lat==lat_WGS84 & lon == lon_WGS84, arr.ind = T)[1],
#            colID = which(lat==lat_WGS84 & lon == lon_WGS84, arr.ind = T)[2])
# climIDtest <- climateKey %>%
#     slice(subset) %>%
#     mutate(rowID = rep(1:272, 214),
#            colID = rep(1:214, each=272))
# identical(ungroup(climID), ungroup(climIDtest)) # ==TRUE

climID <- climateKey %>%
    mutate(lon_ID = rep(1:272, 214),
           lat_ID = rep(1:214, each=272))

## Sense check 1: ----
# do the extracted tiles have the correct configuration (particularly
# do the integer 'tile' values give the correct outline) 
rain <- ncvar_get(ncdf_temp, "tg", start = c(x=1,y=1,t=23922-1000), count = c(272, 214, 1000))
hasClim <- apply(rain, c(1, 2), function(x) sum(!is.na(x)))

climateCheck <- climID %>% ungroup %>%
    mutate(hasClim = c(hasClim))

p1 <- ggplot(climateCheck, aes(lon_WGS84, lat_WGS84, col=hasClim)) + geom_point() +
    labs(caption = bquote(bold("Sense check:") ~ "extracted lats and lons with climate data correspond to Europe outline"))
p2 <- ggplot(climateCheck, aes(lon_ID, lat_ID, fill=hasClim)) + geom_tile() +
    labs(caption = bquote(bold("Sense check:") ~ "extracted lats and lons with climate data correspond to Europe outline"))
pBoth <- egg::ggarrange(p1,p2)
ggsave("figures/check extracted climate/extractedClimate.png",plot=pBoth, height=200, width=160, units="mm")

## Match sites with climate cells ----
# next step is to place on a flattened grid (LAEA)
climateCells <- SpatialPoints(climateKey[c("lon_WGS84", "lat_WGS84")])
proj4string(climateCells) <- CRS_WGS84

climateCells_transformed <- spTransform(climateCells, CRS_LAEA) %>%
    coordinates(.) %>%
    as_data_frame %>%
    rename(lon_LAEA = lon_WGS84, lat_LAEA=lat_WGS84) %>%
    bind_cols(climateKey, .)

## Sense check 2: ----
# does transformation look sensible? 
climateCells_transformed %>%
    select(climateID, lon_LAEA, lat_LAEA) %>%
    left_join(., select(climateCheck, climateID, hasClim)) %>%
    ggplot(aes(lon_LAEA, lat_LAEA, col=hasClim)) + geom_point()
ggsave("figures/check extracted climate/extractedClimate_flattened_LAEA.png", height=100, width=160, units="mm")

# get bird sites (to extract climate for them)
bird_sites <- readRDS("../../PECBMS/cleanedFiles/dataKey_siteID-WGS84coord.rds") %>%
    as_data_frame
bird_sp <- SpatialPoints(bird_sites[c("lon_WGS84", "lat_WGS84")])
proj4string(bird_sp) <- CRS_WGS84

# flatten
birdSites_transformed <- spTransform(bird_sp, CRS_LAEA) %>%
    coordinates(.) %>%
    as_data_frame %>%
    rename(lon_LAEA = lon_WGS84, lat_LAEA=lat_WGS84) %>%
    bind_cols(bird_sites, .)

# Sense check 3: ----
# is spatial object sensible?
plot(bird_sp)
plot(birdSites_transformed[c("lon_LAEA", "lat_LAEA")])
# no point saving these

## function to calculate distances between climate cell centre and site
getClimateCell <- function(x) {
    with(x,
         data_frame(dist = ((lon_LAEA-climateCells_transformed$lon_LAEA)^2 +
                        (lat_LAEA-climateCells_transformed$lat_LAEA)^2)^.5,
         siteID, climateID = climateCells_transformed$climateID)) %>%
        slice(which.min(dist))
}

# match sites: note this takes a minute or two
matchedSites <- split(birdSites_transformed, f=1:nrow(birdSites_transformed)) %>%
    lapply(., getClimateCell) %>%
    bind_rows()

matchedSites_arrInd <- left_join(matchedSites, climID, by="climateID") %>%
    select(climateID, lon_ID, lat_ID) %>% unique

# Sense check 4: ----
# do coords correspond, i.e. is configuration sensible?
full <- left_join(bird_sites, matchedSites) %>%
    full_join(., climateCells_transformed, by="climateID", suffix=c("_sites", "_climateCell"))

p3 <- full %>%
    ggplot(aes(lon_WGS84_climateCell, lon_WGS84_sites)) + geom_point() +
    geom_abline(col="red") + coord_equal()

p4 <- full %>%
    ggplot(aes(lat_WGS84_climateCell, lat_WGS84_sites)) + geom_point() +
    geom_abline(col="red") + coord_equal()
pBoth2 <- egg::ggarrange(p3,p4)
ggsave("figures/check extracted climate/siteIDs vs climateIDs.png",plot=pBoth2, height=200, width=160, units="mm")
ggplot(matchedSites_arrInd, aes(lon_ID, lat_ID)) + geom_tile() +
    coord_equal()
ggsave("figures/check extracted climate/climate cells with bird site.png", height=150, width=150, units="mm")

## Extract climate data ----
# there seems to be a large overhead on the extracting of data, i.e. extracting
# a single cells worth of values takes a number of seconds, 10s even (which I then
# need to repeat 3494 times), while just extracting the lot takes less than a minute
# (although eats almost all of 32GB of RAM; would not work on smaller machine)
# I open sequentially and close and remove array after for this reason

array_rain <- ncvar_get(ncdf_rain, "rr", start = c(1,1,1))
timeLength <- dim(array_rain)[3]

# initialise full dataframe that will get filled with rain and temp values
fullClimateDF <- data_frame(climateID = rep(matchedSites_arrInd$climateID, each = timeLength),
                            lon_ID = rep(matchedSites_arrInd$lon_ID, each = timeLength),
                            lat_ID = rep(matchedSites_arrInd$lat_ID, each = timeLength)) %>%
    group_by(climateID) %>%
    mutate(time = 1:n(),
           rain = rep(NA, n()),
           temp = rep(NA, n()))

# fill rain column
fullClimateDF <- fullClimateDF %>%
    group_by(climateID) %>%
    mutate(rain = array_rain[lon_ID[1], lat_ID[1],])

## close to make space for temperature
nc_close(ncdf_rain)
rm(array_rain)

## repeat for temp
array_temp <- ncvar_get(ncdf_temp, "tg", start = c(1,1,1))

# update temp col
fullClimateDF <- fullClimateDF %>%
    group_by(climateID) %>%
    mutate(temp = array_temp[lon_ID[1], lat_ID[1],])

## close to make space
nc_close(ncdf_temp)
rm(array_temp)


## Sense check5: ----
# do spatial & temporal patterns correspond to expectation?
fullClimateDF_check <- fullClimateDF %>%
    mutate(date = as.Date(time, origin ="1-1-1950", format="%d-%m-%Y"),
           year=year(date),
           day=yday(date),
           month=month(date))

summerWinter <- fullClimateDF_check %>%
    filter(month %in% c(6, 1)) %>%
    group_by(month, climateID, lon_ID, lat_ID) %>%
    summarise(meanT = mean(temp, na.rm=T)) 

# winter temp
filter(summerWinter, month==1) %>%
    ggplot(aes(lon_ID, lat_ID, fill=meanT)) + geom_tile() +
    scale_fill_gradient2(low="lightblue", high="indianred") +
    coord_equal() +
    labs(title="Mean January temperature by climate cell")
ggsave("figures/check extracted climate/spatial_meanJanTemp.png", height=150, width=160, units="mm")

# summer temp
filter(summerWinter, month==6) %>%
    ggplot(aes(lon_ID, lat_ID, fill=meanT)) + geom_tile() +
    scale_fill_gradient2(low="lightblue", high="indianred", midpoint=15) +
    coord_equal() +
    labs(title="Mean June temperature by climate cell")
ggsave("figures/check extracted climate/spatial_meanJuneTemp.png", height=150, width=160, units="mm")

## temporal pattern
set.seed(101)
sample_climID <- sample(unique(fullClimateDF_check$climateID), 10)
cols <- colorRampPalette(c("lightblue", "indianred"))(66)
# monthly pattern
fullClimateDF_check %>%
    filter(climateID %in% sample_climID) %>%
    group_by(climateID, month, year) %>%
    summarise(meanT = mean(temp, na.rm=T), start_date = date[1]) %>%
    ggplot(aes(month, meanT, col=factor(year))) + geom_line() + facet_wrap(~climateID) +
    scale_color_manual(values = cols) +
    scale_x_continuous(breaks=1:12) +
    guides(col=F) +
    labs(title="Monthly mean temps, by year, for subset of climate cells", 
         caption=bquote(bold("Notes: ")*"blue -> red from 1950:2015"))
ggsave("figures/check extracted climate/temporal_monthlyTemp.png", height=150, width=160, units="mm")

# fine-scale daily pattern
fullClimateDF_check %>%
    filter(climateID %in% sample_climID[1:5], year %in% c(1990, 2000, 2010)) %>%
    ggplot(aes(yday(date), temp)) + geom_line() +
    facet_grid(climateID~year) +
    labs(title="Daily temps, by year, for subset of climate cells")
ggsave("figures/check extracted climate/temporal_dailyTemp.png", height=150, width=160, units="mm")
# yearly pattern
fullClimateDF_check %>%
    filter(climateID %in% sample_climID) %>%
    group_by(climateID, year) %>%
    summarise(meanT = mean(temp, na.rm=T), start_date = date[1]) %>%
    ggplot(aes(year, meanT, col=factor(year))) + geom_point() + facet_wrap(~climateID) +
    scale_color_manual(values = cols) +
    # scale_x_continuous(breaks=1:12) +
    guides(col=F) +
    labs(title="Annual mean temp for subset of climate cells", 
         caption=bquote(bold("Notes: ")*"blue -> red from 1950:2015"))
ggsave("figures/check extracted climate/temporal_annualTemp.png", height=150, width=160, units="mm")

## Save outputs ----
# save climate key (that relates siteID to climate ID), and climate datasets 
saveRDS(fullClimateDF, "files/fullClimateDF.rds")

# get monthly summaries (note: slow)
summClimateDF <- fullClimateDF_check %>% 
    group_by(climateID, year, month) %>%
    summarise(temp_monthlyMean = mean(temp, na.rm=T), 
              precip_monthlySum = sum(rain, na.rm=T))
# remove cells with no data
summClimateDF %>%
    filter(!is.nan(temp_monthlyMean)) %>%
saveRDS(., "files/monthlySummaryClimateDF.rds")

# climateKey

full %>%
    select(schemeID, siteID, climateID, 
           climLat_WGS84 = lat_WGS84_climateCell, 
           climLon_WGS84 = lon_WGS84_climateCell) %>%
    left_join(., select(climID, climateID, lon_ID, lat_ID)) %>%
    saveRDS(., "files/climateKey_siteID-climateID.rds")
