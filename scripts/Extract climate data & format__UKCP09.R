## Extract climate data for bird monitoring sites in the UK, using the 5km 
## resolution UKCP09 dataset. 

## housekeeping ----
library(ncdf4); library(ggplot2); library(raster); library(dplyr); library(lubridate)

BNG <- CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

## read site data ----
## create spatial object of sites to use to extract data from climate raster
# get sites from bird dataset
PECBMS_site <- readRDS("files/birdDF.rds") %>%
    dplyr::select(schemeID, siteID, lon_WGS84, lat_WGS84) %>%
    unique

UKBMS <- filter(PECBMS_site,schemeID=="UKI.0")
UKBMS_sites <- UKBMS$siteID

UKBMS_spatial <- SpatialPoints(dplyr::select(UKBMS , lon_WGS84, lat_WGS84))
proj4string(UKBMS_spatial) <- WGS84
UKBMS_BNG <- spTransform(UKBMS_spatial, BNG)

# check spatial configuration -looks fine, and this check is replicated in climate
# plots below
# plot(UKBMS_BNG)

## extract climatic data ----
# get UKCP09 filenames, for temperature and rainfall respectively
fnames_temp <- c(list.files("../UKCP09/download_1990/", full.names = TRUE),
            list.files("../UKCP09/download_2000s/", full.names = TRUE))
fnames_rain <- c(list.files("../UKCP09/rainfall_1990s/", full.names = TRUE, pattern=".nc"),
            list.files("../UKCP09/rainfall_2000s/", full.names = TRUE, pattern=".nc"))
years <- gsub(".*rainfall_(....).*", "\\1", fnames_rain)
# note years in temp and rainfall are identical

# iterate across years, extracting temperature data
extracted_temp <- list()
for(i in seq(fnames_temp)) {
    print(years[i])
    r <- brick(fnames_temp[i])
    proj4string(r) <- BNG
    extracted <- raster::extract(r, UKBMS_BNG)
    extracted_temp[[i]] <- as_data_frame(extracted) %>%
        mutate(siteID = 1:n()) %>%
        reshape2::melt(., id.vars="siteID") %>%
        as_data_frame %>%
        group_by(siteID) %>%
        mutate(day = 1:length(siteID), year=years[i])
}

# iterate across years, extracting rainfall data
extracted_rain <- list()
for(i in seq(fnames_rain)) {
    print(years[i])
    r <- brick(fnames_rain[i])
    proj4string(r) <- BNG
    extracted <- raster::extract(r, UKBMS_BNG)
    extracted_rain[[i]] <- as_data_frame(extracted) %>%
        mutate(siteID = 1:n()) %>%
        reshape2::melt(., id.vars="siteID") %>%
        as_data_frame %>%
        group_by(siteID) %>%
        mutate(day = 1:length(siteID), year=years[i])
    
}

## format extracted data ----
UKCP_temp <- bind_rows(extracted_temp) %>%
    ungroup %>% 
    mutate(siteID = UKBMS_sites[siteID],
           date = as.Date(paste(year, day), format="%Y %j"),
           month = month(date)) %>%
    dplyr::select(siteID, year, month, day, date, temp=value)

UKCP_rain <- bind_rows(extracted_rain) %>%
    ungroup %>% 
    mutate(siteID = UKBMS_sites[siteID],
           date = as.Date(paste(year, day), format="%Y %j"),
           month = month(date)) %>%
    dplyr::select(siteID, year, month, day, date, rain=value)

# join these to form single df (note: slowish)
UKCP_full <- full_join(UKCP_temp, dplyr::select(UKCP_rain, siteID, date, rain))

## sensibility check ----
## Are spatial and temporal patterns sensible?

# temporal plot
ggplot(filter(UKCP_full, siteID %in% UKBMS_sites[c(1,5,50,500,3000)]), aes(date, temp)) + 
    geom_line() +
    scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
    theme(axis.text.x = element_text(angle = 90, vjust=0.5)) +
    facet_wrap(~siteID, ncol=1) +
    labs(y="Daily temperature", x="Time", 
         caption=bquote(bold("Figure: ")*"daily temperatures for 5 sites in UK, using UKCP09"))
ggsave("figures/check extracted climate/UKCP09_temporalPattern.png", height=290, width=210, units="mm")

# spatial plots
coord_BNG <- cbind(UKBMS_sites, coordinates(UKBMS_BNG)) %>%
    as_data_frame %>%
    dplyr::select(siteID=1, lon_BNG=2, lat_BNG=3)

temp_monthlySumm <- UKCP_full %>%
    group_by(siteID, year, month) %>%
    summarise(meanT = mean(temp)) %>%
    left_join(., coord_BNG) %>%
    mutate(lon_BNG = as.numeric(lon_BNG), 
           lat_BNG = as.numeric(lat_BNG))

theme_spatial <- theme_classic() +
    theme(legend.position = "bottom", 
          panel.background = element_rect(fill="grey40"), 
          strip.background = element_rect(colour=NA), 
          strip.text = element_text(hjust=0, face="bold"), 
          axis.text = element_blank(), axis.ticks=element_blank()) 

# June temperature
temp_June <- temp_monthlySumm %>% 
    filter(month == 6) %>%
    filter(year %in% c(1990, 2000, 2010, 2014))

ggplot(temp_June, aes(lon_BNG, lat_BNG, col=meanT)) + 
    geom_point() +
    facet_wrap(~year, nrow=1) +
    scale_colour_gradient2(midpoint=median(temp_June$meanT, na.rm=T), 
                           low = "lightblue3", high="indianred") +
    coord_fixed() +
    theme_classic() +
    theme_spatial +
    labs(caption=bquote(bold("Figure:")~"June temperatures for three time-slices, using UKCP09 dataset"))
ggsave("figures/check extracted climate/UKCP09_spatialPattern_June.png", height=210, width=290, units="mm")


# Jan temperature
temp_Jan <- temp_monthlySumm %>% 
    filter(month == 1) %>%
    filter(year %in% c(1990, 2000, 2010, 2014))

ggplot(temp_Jan, aes(lon_BNG, lat_BNG, col=meanT)) + 
    geom_point() +
    facet_wrap(~year, nrow=1) +
    scale_colour_gradient2(midpoint=median(temp_Jan$meanT, na.rm=T), 
                           low = "lightblue3", high="indianred") +
    coord_fixed() +
    theme_classic() +
    theme_spatial +
    labs(caption=bquote(bold("Figure:")~"January temperatures for three time-slices, using UKCP09 dataset"))
ggsave("figures/check extracted climate/UKCP09_spatialPattern_Jan.png", height=210, width=290, units="mm")

## save files ----
# full df
saveRDS(UKCP_full, "files/UKCP09_full.rds")

# monthly averages 
UKCP_summ <- UKCP_full %>% 
    group_by(siteID, year, month) %>%
    summarise(temp_monthlyMean = mean(temp, na.rm=T), 
              precip_monthlySum = sum(rain, na.rm=T))
saveRDS(UKCP_summ, "files/UKCP09_monthlySumm_UK monitoring sites.rds")
