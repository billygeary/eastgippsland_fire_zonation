
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Prepare Species Occurrence Data for East Gippsland Project ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Three data sets included in this analysis
### a) Southern Ark 2016-2017 Monitoring
### b) Long-footed Potoroo Monitoring at Tulloch Ard
### c) Spot-tailed Quoll Monitoring at W Tree

#### Setup ####
library(sf)
library(raster)
library(fasterize)
library(ggplot2)
library(tidyverse)

std.crs = st_crs(28355) # Standard CRS will be EPSG:28355

# Read in Site Info
camera_data = read.csv(here::here("Data_Processing", "sark_lfp_stq_1617_presenceabsence_allspp.csv"))
camera_coords = camera_data %>% 
  distinct(SiteID,Program, DateDeployed,Lure, Easting, Northing) %>% 
  st_as_sf(coords = c('Easting', 'Northing'), crs = std.crs)
dates = unique(camera_data$DateDeployed)

# Read in spatial data to create Study Area Footprint
sark_footprint = read_sf(here::here("Data_Raw", "spatial", "SouthernArkFootprint.shp"))
sark_footprint = sark_footprint %>% st_transform(std.crs)
victoria = read_sf("~/Library/CloudStorage/Dropbox/Billy/_research/data/spatial/outline_shapefiles/Vic Outline/STE_2021_AUST_GDA2020.shp") %>%
  filter(STE_NAME21=="Victoria") %>% st_transform(std.crs)

# Visualise sites to check it all looks good
site_plot <- ggplot() +
  geom_sf(data = sark_footprint, fill = "grey", inherit.aes = FALSE, alpha = 0.2, colour = "black") +
  geom_sf(mapping = aes(color=Program), data = camera_coords, size = 2, shape = 1, inherit.aes = FALSE) +
  scale_color_manual(values = c(`Southern Ark` = 'green', `W Tree` = 'black', `Quoll` = 'blue')) +
  ylab("Latitude") +
  xlab("Longitude") +
  theme(axis.title = element_text(size = 20), 
        axis.text = element_text(size = 15), 
        axis.text.x = element_text(hjust = -1),
        legend.title = element_text(size=24), 
        legend.text = element_text(size=18),
        legend.key.size = unit(7, 'cm')) +
  theme_bw()
site_plot

#### 1: Build base raster ####
study_extent = camera_coords %>% st_buffer(dist=5000) %>% st_union(sark_footprint)

base_raster = raster(ext = extent(study_extent), resolution = 50, crs = std.crs$input)
base_raster = fasterize::fasterize(victoria, base_raster, background = 0) # Create base raster, all cells in study region = 1
names(base_raster) <- "base_raster"
base_raster_mask = base_raster
base_raster_mask[base_raster_mask==0] <- NA

plot(base_raster_mask)
plot(camera_coords, add=TRUE)
plot(sark_footprint, add=TRUE)

#### 2: Fox Baiting Density ####
baiting_data = read.csv(here::here("Data_Raw","NEP_SouthernArk", "baiting_data_20170228.csv"))

bait_coords = baiting_data %>% 
  distinct(uniqueid_t, easting, northing) %>% 
  st_as_sf(coords = c('easting', 'northing'), crs = std.crs) %>%
  mutate(Dummy =1)

bait_station_raster = rasterize(bait_coords, base_raster, field = 'Dummy', background = 0)
window_2500m = raster::focalWeight(base_raster, d = 2500, type='circle')
window_2500m[window_2500m > 0] <- 1 #Number of baits or proportion of cells?
bait_density_2.5km = focal(bait_station_raster, w=window_2500m, fun ='sum')

names(bait_density_2.5km) <- "bait_intensity"

plot(bait_density_2.5km)
rm(baiting_data)

#### 3: Time Since Fire Pre Survey ####
firehistory = read_sf(here::here("Data_Raw", "spatial", "FIRE_HISTORY.shp"))

# First need to apply burn cover model to clean up fire history data
burncover = raster(here::here("Data_Raw","spatial", "BurnCovModYN.tif"))

base_raster_transformed = projectRaster(base_raster, crs = crs(burncover))
burncover = resample(burncover, base_raster_transformed, "ngb")

# Make nv mask
evcs = sf::read_sf(here::here("Data_Raw", "spatial", "EVCs_gippsland.shp")) %>% st_transform(st_crs(burncover))
evcs$Dummy = 1
nv = fasterize::fasterize(evcs, base_raster_transformed, field = "Dummy")

rm(evcs)

firehistoryT = 
  firehistory %>% 
  st_transform(st_crs(burncover))

### Loop through each fire SEASON to create raster of fires
seasons = sort(unique(firehistory$SEASON))
seasons = subset(seasons, seasons < 2017)

fire.stack = list()
burn.stack = list()
for (i in 1:length(seasons)){
  s = seasons[i]
  fires = firehistoryT %>% filter(SEASON == s)
  fires.certain = fires %>% filter(FIRE_SVRTY %in% c("BURNT_1", "BURNT_2", "BURNT_2F", "BURNT_2P", "BURNT_3", "BURNT_NONFOREST") |
                                     ACCURACY %in% c("High - 25m or less", "Medium - 26m to 100m"))
  fires.certain.ras = base_raster_transformed
  values(fires.certain.ras) <- NA
  if(length(fires.certain$SEASON) > 0){fires.certain.ras = fasterize::fasterize(fires.certain, base_raster_transformed, field = "SEASON", background = NA)}
  
  fires.uncertain = fires %>% filter(FIRE_SVRTY == "BURNT_UNKNOWN" | ACCURACY %in% c("Low - greater than 100m", "Unknown"))
  fires.uncertain.ras = base_raster_transformed
  if(length(fires.uncertain$SEASON) > 0){fires.uncertain.ras = fasterize::fasterize(fires.uncertain, base_raster_transformed, field = "SEASON", background = NA)}
  fires.uncertain.ras.burncover = fires.uncertain.ras * burncover
  fires.uncertain.ras.burncover[fires.uncertain.ras.burncover == 0] <- NA
  fire.ras = mosaic(fires.certain.ras, fires.uncertain.ras.burncover, fun = "max", na.rm=TRUE)
  names(fire.ras) <- s
  fire.ras[is.na(fire.ras)] <-0
  fire.stack[[i]] <- fire.ras
  burn.ras = fire.ras; burn.ras[burn.ras>0] <-1
  burn.stack[[i]] <- burn.ras
}

rm(firehistory)
rm(evcs)
#Add other arguments for mosaic function to the list.
fire.stack$fun <- max
fire.stack$na.rm <- TRUE

#Use do.call to apply the resulting "list of arguments" to the mosaic 
#function.
fireage <- do.call(mosaic, fire.stack)
# Patch in nv with no fire history
nv[nv ==1] = 1900
fireage = mosaic(fireage, nv, fun = 'max')
fireage[fireage<1900] <- NA

tsf_startsurvey_burn = lubridate::year(min(as.Date(dates, format="%d/%m/%Y"))) - fireage
names(tsf_startsurvey_burn) <- "tsf_startsurvey_burn"
tsf_startsurvey_burn <- projectRaster(tsf_startsurvey_burn, crs = crs(base_raster), method = 'ngb')
tsf_startsurvey_burn = resample(tsf_startsurvey_burn, base_raster, method = 'ngb')
tsf_endsurvey_burn = lubridate::year(min(as.Date(dates, format="%d/%m/%Y"))) - fireage
names(tsf_endsurvey_burn) <- "tsf_endsurvey_burn"
tsf_endsurvey_burn <- projectRaster(tsf_endsurvey_burn, crs = crs(base_raster), method = 'ngb')
tsf_endsurvey_burn = resample(tsf_endsurvey_burn, base_raster, method = 'ngb')

# Make a prop recent landscape 
recent_burn = tsf_startsurvey_burn
recent_burn[recent_burn < 11] <- 1
recent_burn[recent_burn > 10] <- 0 
recent_burn[is.na(recent_burn)] <- 0
window_1000m = raster::focalWeight(base_raster, d = 1000, type='circle')
prop_recent_burn = focal(recent_burn, w=window_1000m, fun = 'sum')
names(prop_recent_burn) <- "prop_recentburn"
# Make a prop long unburnt landscape
longunburnt_burn = tsf_startsurvey_burn
longunburnt_burn[longunburnt_burn < 40] <- 0
longunburnt_burn[longunburnt_burn > 39] <- 1
longunburnt_burn[is.na(longunburnt_burn)] <- 0
window_1000m = raster::focalWeight(base_raster, d = 1000, type='circle')
prop_longunburnt_burn = focal(longunburnt_burn, w=window_1000m, fun = 'sum')
names(prop_longunburnt_burn) <- "prop_longunburnt"


#### 4: Fire frequency Pre Survey ####
# Can be static because there were no fires in the study period
# burn.stack.sub = burn.stack[21:60]
# burn.stack.sub$fun <- sum
# burn.stack.sub$na.rm <- TRUE
# burnt_since_1978 <- do.call(mosaic, burn.stack.sub)

burn.stack.sub = burn.stack[41:60]
burn.stack.sub$fun <- sum
burn.stack.sub$na.rm <- TRUE
burnt_since_1998 <- do.call(mosaic, burn.stack.sub)

rm(fire.stack)
rm(burn.stack)
rm(burn.stack.sub)
gc()

burnt_twice = burnt_since_1998
burnt_twice[burnt_twice == 1] <- 0
burnt_twice[burnt_twice >1] <- 1
plot(burnt_twice)
prop_burnt_twice = focal(burnt_twice, w = window_1000m, fun='sum')
prop_burnt_twice <- projectRaster(prop_burnt_twice, crs = crs(base_raster))
prop_burnt_twice = resample(prop_burnt_twice , base_raster, method = 'ngb')
names(prop_burnt_twice) = "prop_burnt_twice"

window_1000m[window_1000m > 0] <- 1 #Number of baits or proportion of cells?
mean_repeat_burn = focal(burnt_since_1998, w = window_1000m, fun='mean')
mean_repeat_burn <- projectRaster(mean_repeat_burn, crs = crs(base_raster))
mean_repeat_burn = resample(mean_repeat_burn, base_raster)
names(mean_repeat_burn) = "mean_repeat_burn"

#### 3: Time Since Fire Post 2019/2020 ####
## Repeat for Post 2019/2020 Fires
### Loop through each fire SEASON to create raster of fires
seasons = sort(unique(firehistoryT$SEASON))
seasons = subset(seasons, seasons < 2021)

# Read in the severity raster for the 19/20 fires and resample 
# to the covariate raster format
sev = raster(here::here("Data_Raw", "spatial", "FireSeverityFinal_20200414.tif"))
sev = projectRaster(sev, crs=crs(fireage))
sev.burnt = sev
sev.burnt[sev.burnt < 3] <- NA 
sev.burnt[sev.burnt > 2] <- 2020 
sev.re = resample(sev.burnt, fireage, "ngb")

fire.stack = list()
burn.stack = list()
for (i in 1:length(seasons)){
  s = seasons[i]
  fires = firehistoryT %>% filter(SEASON == s)
  fires.certain = fires %>% filter(FIRE_SVRTY %in% c("BURNT_1", "BURNT_2", "BURNT_2F", "BURNT_2P", "BURNT_3", "BURNT_NONFOREST") |
                                     ACCURACY %in% c("High - 25m or less", "Medium - 26m to 100m"))
  fires.certain.ras = base_raster_transformed
  values(fires.certain.ras) <- NA
  if(length(fires.certain$SEASON) > 0){fires.certain.ras = fasterize::fasterize(fires.certain, base_raster_transformed, field = "SEASON", background = NA)}
  
  fires.uncertain = fires %>% filter(FIRE_SVRTY == "BURNT_UNKNOWN" | ACCURACY %in% c("Low - greater than 100m", "Unknown"))
  fires.uncertain.ras = base_raster_transformed
  if(length(fires.uncertain$SEASON) > 0){fires.uncertain.ras = fasterize::fasterize(fires.uncertain, base_raster_transformed, field = "SEASON", background = NA)}
  fires.uncertain.ras.burncover = fires.uncertain.ras * burncover
  fires.uncertain.ras.burncover[fires.uncertain.ras.burncover == 0] <- NA
  fire.ras = mosaic(fires.certain.ras, fires.uncertain.ras.burncover, fun = "max", na.rm=TRUE)
  if(s == 2020){fire.ras = mosaic(fire.ras, sev.re, fun="max", na.rm=TRUE)}
  names(fire.ras) <- s
  fire.ras[is.na(fire.ras)] <-0
  fire.stack[[i]] <- fire.ras
  burn.ras = fire.ras; burn.ras[burn.ras>0] <-1
  burn.stack[[i]] <- burn.ras
}

#Add other arguments for mosaic function to the list.
fire.stack$fun <- max
fire.stack$na.rm <- TRUE

#Use do.call to apply the resulting "list of arguments" to the mosaic 
#function.
fireage_post1920 <- do.call(mosaic, fire.stack)
# Patch in nv with no fire history
nv[nv ==1] = 1900
fireage_post1920  = mosaic(fireage_post1920, nv, fun = 'max')
fireage_post1920[fireage_post1920 <1900] <- NA

rm(fire.stack)
gc()
new.dates = "01/01/2022"
tsf_startsurvey_burn_post1920  = lubridate::year(min(as.Date(new.dates, format="%d/%m/%Y"))) - fireage_post1920 
names(tsf_startsurvey_burn_post1920) <- "tsf_startsurvey_burn_post1920"
tsf_startsurvey_burn_post1920  <- projectRaster(tsf_startsurvey_burn_post1920 , crs = crs(base_raster), method = 'ngb')
tsf_startsurvey_burn_post1920  = resample(tsf_startsurvey_burn_post1920 , base_raster, method = 'ngb')
tsf_endsurvey_burn_post1920  = lubridate::year(min(as.Date(new.dates, format="%d/%m/%Y"))) - fireage_post1920 
names(tsf_endsurvey_burn_post1920) <- "tsf_endsurvey_burn_post1920"
tsf_endsurvey_burn_post1920 <- projectRaster(tsf_endsurvey_burn_post1920, crs = crs(base_raster), method = 'ngb')
tsf_endsurvey_burn_post1920 = resample(tsf_endsurvey_burn_post1920, base_raster, method = 'ngb')

# Make a prop recent landscape 
recent_burn_post1920 = tsf_startsurvey_burn_post1920
recent_burn[recent_burn_post1920 < 11] <- 1
recent_burn_post1920[recent_burn_post1920 > 10] <- 0 
recent_burn_post1920[is.na(recent_burn_post1920)] <- 0
window_1000m = raster::focalWeight(base_raster, d = 1000, type='circle')
prop_recent_burn_post1920 = focal(recent_burn_post1920, w=window_1000m, fun = 'sum')
names(prop_recent_burn_post1920) <- "prop_recentburn_post1920"

# Make a prop long unburnt landscape
longunburnt_burn_post1920 = tsf_startsurvey_burn_post1920
longunburnt_burn_post1920[longunburnt_burn_post1920 < 40] <- 0
longunburnt_burn_post1920[longunburnt_burn_post1920 > 39] <- 1
longunburnt_burn_post1920[is.na(longunburnt_burn_post1920)] <- 0
window_1000m = raster::focalWeight(base_raster, d = 1000, type='circle')
prop_longunburnt_burn_post1920 = focal(longunburnt_burn_post1920, w=window_1000m, fun = 'sum')
names(prop_longunburnt_burn_post1920) <- "prop_longunburnt_post1920"


#### 4: Fire frequency Post 2019/2020 ####
# Can be static because there were no fires in the study period
# burn.stack.sub_post1920 = burn.stack[21:60]
# burn.stack.sub_post1920$fun <- sum
# burn.stack.sub_post1920$na.rm <- TRUE
# burnt_since_last40_post1920 <- do.call(mosaic, burn.stack.sub_post1920)

burn.stack.sub_post1920 = burn.stack[67:86]
burn.stack.sub_post1920$fun <- sum
burn.stack.sub_post1920$na.rm <- TRUE
burnt_since_last20_post1920 <- do.call(mosaic, burn.stack.sub_post1920)

burnt_twice_post1920 = burnt_since_last20_post1920
burnt_twice_post1920[burnt_twice_post1920 == 1] <- 0
burnt_twice_post1920[burnt_twice_post1920 >1] <- 1
plot(burnt_twice_post1920)
prop_burnt_twice_post1920 = focal(burnt_twice_post1920, w = window_1000m, fun='sum')
prop_burnt_twice_post1920 <- projectRaster(prop_burnt_twice_post1920, crs = crs(base_raster))
prop_burnt_twice_post1920 = resample(prop_burnt_twice_post1920, base_raster, method = 'ngb')
names(prop_burnt_twice_post1920) = "prop_burnt_twice_post1920"

window_1000m[window_1000m > 0] <- 1 #Number of baits or proportion of cells?
mean_repeat_burn_post1920 = focal(burnt_since_last20_post1920, w = window_1000m, fun='mean')
mean_repeat_burn_post1920 <- projectRaster(mean_repeat_burn_post1920, crs = crs(base_raster))
mean_repeat_burn_post1920 = resample(mean_repeat_burn_post1920, base_raster)
names(mean_repeat_burn_post1920) = "mean_repeat_burn_post1920"

### Ten years post 19-20 fires
# Time since fire
tsf_burn_10yrspost1920 = tsf_endsurvey_burn_post1920 + 10
names(tsf_burn_10yrspost1920) = "tsf_burn_10yrspost1920"

# Mean repeat burn
burn.stack.sub_10yrspost1920 = burn.stack[77:86] # last ten years of dataset
burn.stack.sub_10yrspost1920$fun <- sum
burn.stack.sub_10yrspost1920$na.rm <- TRUE
burnt_since_last20_10yrspost1920 <- do.call(mosaic, burn.stack.sub_10yrspost1920)

window_1000m[window_1000m > 0] <- 1 #Number of baits or proportion of cells?
mean_repeat_burn_10yrspost1920 = focal(burnt_since_last20_10yrspost1920, w = window_1000m, fun='mean')
mean_repeat_burn_10yrspost1920 <- projectRaster(mean_repeat_burn_10yrspost1920, crs = crs(base_raster))
mean_repeat_burn_10yrspost1920 = resample(mean_repeat_burn_10yrspost1920, base_raster)
names(mean_repeat_burn_10yrspost1920) = "mean_repeat_burn_10yrspost1920"

rm(burn.stack)
rm(burn.stack.sub_post1920)
gc()

#### 5: Timber Harvesting ####
harvest_history = read_sf(here::here("Data_Raw", "spatial", "LOG_SEASON.shp"))
harvest_history = harvest_history %>% 
  st_transform(std.crs) %>% 
  mutate(ENDDATE = as.Date(ENDDATE, format = "%Y-%m-%d"))
# 
# dates = unique(camera_data$DateDeployed)
# 
# # Extract survey specific tsf for each site
# harvest.list = list()
# for (i in 1:length(dates)) {
#   start.date = dates[i]
#   harvest_history$TSH = as.numeric(difftime(as.Date(start.date, format="%d/%m/%Y"),harvest_history$ENDDATE, unit="weeks"))/52.25
#   #time since harvest
#   harvest_history_temp = subset(harvest_history, TSH > 0)
#   tsh_raster = fasterize(harvest_history_temp, base_raster, field="TSH", fun='min')
#   extract_traps = camera_coords %>% dplyr::filter(DateDeployed == start.date)
#   tsh.point = raster::extract(tsh_raster, extract_traps)
#   harvest.data = cbind(extract_traps, tsh.point)
#   harvest.list[[as.character(start.date)]] = harvest.data
# }
# tsh.out = do.call(rbind, harvest.list)
# tsh.out = tsh.out %>% st_drop_geometry()
# 

harvest_history$TSH = as.numeric(difftime(max(as.Date(dates, format="%d/%m/%Y")),harvest_history$ENDDATE, unit="weeks"))/52.25

tsh_endsurvey = fasterize::fasterize(harvest_history, base_raster, field = "TSH", fun="min")
plot(tsh_endsurvey)

harvest_history_pa = harvest_history %>%
  filter(X_SILVSYS %in% c("Clearfelling", "Seed Tree (includes retained overwood)","Clearfelling Salvage",
                          "Clearfellilng Salvage", "Unknown")) %>%
  filter(TSH < 40 & TSH>0)
harvest_pa = fasterize::fasterize(harvest_history_pa, base_raster, field = "TSH", fun = "min", background = 0)
harvest_pa[harvest_pa > 0] <- 1
window_1000m = raster::focalWeight(base_raster, d = 1000, type='circle')
harvest_prop_40 <- focal(harvest_pa, w = window_1000m, fun = 'sum')
names(harvest_prop_40) <- "harvest_prop_40"

#### 7: Rainfall deciles ####
rain_deciles_2015 = raster(here::here("Data_Raw", "spatial", "precip_percentile_r005_20150101_20151231.nc"))
rain_deciles_2015 = projectRaster(rain_deciles_2015, crs=std.crs$input)
rain_deciles_2015 = resample(rain_deciles_2015, base_raster)
names(rain_deciles_2015) <- "rain_deciles_2015"

#### 8: Long-term rainfall ####
mean_annual_rainfall = raster("~/Library/CloudStorage/Dropbox/Billy/_Research/_PhD/07_UpperWarren_WA/data/Climate Data/rainan 4.txt")
crs(mean_annual_rainfall) <- "+proj=longlat +datum=WGS84"
mean_annual_rainfall = projectRaster(mean_annual_rainfall, crs = "EPSG:28355")
mean_annual_rainfall = resample(mean_annual_rainfall, base_raster)
names(mean_annual_rainfall) <- 'mean_annual_rainfall'

#### 9: Topographic Wetness Index ####
twi = raster("~/Library/CloudStorage/Dropbox/Billy/_Research/Data/spatial/twi/twi_3s.tif")
clippoly = study_extent %>% st_buffer(dist=5000) %>% st_transform(st_crs(twi)) 
twi = raster::crop(twi, extent(clippoly))
twi = projectRaster(twi, crs="EPSG:28355")
twi = resample(twi, base_raster)
names(twi) <- "topographic_wetness_index"

#### 10: Burnability Raw ####
burnability = raster(here::here("Data_Raw","spatial", "BurnCovModInt.tif"))
burnability = projectRaster(burnability, crs="EPSG:28355")
burnability = resample(burnability, twi, method = 'bilinear')
names(burnability) <- "burnability"

#### 12: Soil Grids ####
clay = raster(here::here("Data_Raw", "spatial", "CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif"))
clay = projectRaster(clay, crs="EPSG:28355")
clay = resample(clay, twi, method = 'bilinear')
names(clay) <- "clay"

carbon = raster(here::here("Data_Raw", "spatial", "SOC_000_005_EV_N_P_AU_TRN_N_20220727.tif"))
carbon = projectRaster(carbon, crs="EPSG:28355")
carbon = resample(carbon, twi, method = 'bilinear')
names(carbon) <- "carbon"

sand = raster(here::here("Data_Raw", "spatial", "SND_000_005_EV_N_P_AU_TRN_N_20210902.tif"))
sand = projectRaster(sand, crs="EPSG:28355")
sand = resample(sand, twi, method = 'bilinear')
names(sand) <- "sand"


#### 11: Easting and Northing ####
xm<-matrix(xFromCell(base_raster_mask,c(1:ncell(base_raster))),nrow=nrow(base_raster),byrow=TRUE)
easting = raster(xm, xmn=592500, xmx=762900, ymn=5807531, ymx=5930631)
crs(easting) <- crs(base_raster)
names(easting)<-"easting"
ym<-matrix(yFromCell(base_raster_mask,c(1:ncell(base_raster))),nrow=nrow(base_raster),byrow=TRUE)
northing = raster(ym, xmn=592500, xmx=762900, ymn=5807531, ymx=5930631)
crs(northing) <- crs(base_raster_mask)
names(northing)<-"northing"


#### 12: Bring Together ####
covariate_stack = raster::stack(base_raster_mask, easting, northing,
                        # Fox Baiting
                        bait_density_2.5km,
                        # Fire
                        tsf_startsurvey_burn, tsf_endsurvey_burn,
                        mean_repeat_burn, prop_recent_burn, 
                        prop_longunburnt_burn, prop_burnt_twice,
                        tsf_startsurvey_burn_post1920, tsf_endsurvey_burn_post1920,
                        mean_repeat_burn_post1920, prop_recent_burn_post1920, 
                        prop_longunburnt_burn_post1920, prop_burnt_twice_post1920,
                        tsf_burn_10yrspost1920, mean_repeat_burn_10yrspost1920,
                        
                        # Harvesting
                        tsh_endsurvey, harvest_prop_40,
                        # Rainfall
                        rain_deciles_2015, 
                        mean_annual_rainfall,
                        # Topographic Wetness Index
                        twi,
                        # Soil
                        clay, carbon, sand,
                        burnability)

covariate_stack_masked = mask(covariate_stack, base_raster_mask)

saveRDS(covariate_stack, "Data_Clean/covariate_stack_updated.RData")
saveRDS(covariate_stack_masked, "Data_Clean/covariate_stack_masked.RData")

#### 13: Ecological Vegetation Classes ####
evcs = sf::read_sf(here::here("Data_Raw", "spatial", "EVCs_gippsland.shp")) %>% st_transform(std.crs)
evcs$Dummy = 1
nv = fasterize::fasterize(evcs, base_raster, field = "Dummy")

evcs_sites = evcs %>% dplyr::select(XGROUPNAME) %>% st_intersection(camera_coords) %>%
  mutate(ForestType = ifelse(XGROUPNAME %in% c("Dry Forests", "Rocky Outcrop or Escarpment Scrubs"), "DryForest",
                             ifelse(XGROUPNAME %in% c("Lowland Forests", "Heathlands"), "LowlandForest",
                                    ifelse(XGROUPNAME %in% c("Wet or Damp Forests", "Rainforests"), "WetForest",
                                           ifelse(XGROUPNAME %in% c("Coastal Scrubs Grasslands and Woodlands", "Lower Slopes or Hills Woodlands", 
                                                                    "Montane Grasslands, Shrublands or Woodlands","Riparian Scrubs or Swampy Scrubs and Woodlands",
                                                                    "Sub-alpine Grasslands, Shrublands or Woodlands"), 
                                                  "Woodlands", XGROUPNAME))))) %>%
  st_drop_geometry()


#### 14: Extract for Sites ####

#Extract static data
camera_modelling = camera_data %>%
  mutate(bait_intensity = raster::extract(bait_density_2.5km, camera_coords),
         rain_deciles_2015 = raster::extract(rain_deciles_2015, camera_coords),
         mean_rainfall = raster::extract(mean_annual_rainfall, camera_coords),
         tsf_burn = raster::extract(tsf_startsurvey_burn, camera_coords),
         twi = raster::extract(twi, camera_coords),
         burnability = raster::extract(burnability, camera_coords),
         mean_repeat_burn = raster::extract(mean_repeat_burn, camera_coords),
         harvest_prop_40 = raster::extract(harvest_prop_40, camera_coords),
         prop_recent_burn = raster::extract(prop_recent_burn, camera_coords),
         prop_burnt_twice = raster::extract( prop_burnt_twice, camera_coords),
         prop_longunburnt_burn = raster::extract(prop_longunburnt_burn, camera_coords),
         clay = raster::extract(clay, camera_coords),
         sand = raster::extract(sand, camera_coords),
         carbon = raster::extract(carbon, camera_coords))




# Add date specific data
#camera_modelling = left_join(camera_modelling, tsf.out, by=c("Program", "Lure", "SiteID", "DateDeployed"))
#camera_modelling = left_join(camera_modelling, tsh.out, by=c("Program", "Lure", "SiteID", "DateDeployed"))

# Add evc data
camera_modelling = left_join(camera_modelling, evcs_sites, by=c("Program", "Lure", "SiteID", "DateDeployed"))

# Fill NAs
#camera_modelling$tsf.point[is.na(camera_modelling$tsf.point)] <- 80 # Make NAs the max tsf in study region
#camera_modelling$tsf_burn[is.na(camera_modelling$tsf_burn)] <- max(camera_modelling$tsf_burn, na.rm=TRUE) # Make NAs the max tsf in study region
#camera_modelling$tsh.point[is.na(camera_modelling$tsh.point)] <- 60 # Make NAs the max tsh in study region
#camera_modelling$fire_frequency_1980[is.na(camera_modelling$fire_frequency_1980)] <- 0 # make fire frequency since 1980 zero as these sites are unburnt
camera_modelling$burnability[is.na(camera_modelling$burnability)] <- 0 # make burnability zero as likely area with no native veg

summary(camera_modelling)

# Check correlations
M <- cor(camera_modelling%>%
           dplyr::select(48:(ncol(camera_modelling))), use='complete.obs', method='pearson')           

# Find pairs that shouldnt be in the same model
cors = as.data.frame(as.table(M))

corplot = ggplot(cors) + geom_tile(aes(x=Var1, y=Var2, fill=Freq)) + scale_fill_distiller(palette = "Spectral") +
  theme(axis.text.x = element_text(angle = 90))
corplot

high.cors = cors %>% filter(Freq < -0.7 | Freq >0.7) %>% filter(Freq !=1) %>%
  group_by(Var1) %>% summarise(Variables=paste(Var2, collapse=", "))

#### 13: Save data frame for modelling
write.csv(camera_modelling, here::here("Data_Clean", "pa_covariates_sark_lfp_stq_allspp.csv"))


