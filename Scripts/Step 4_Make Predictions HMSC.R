# Step 4 - Make the model predictions (HMSC)
library(dplyr)
library(Hmsc)
library(tidyr)
library(AHMbook)
library(parallel)
library(foreach)
library(doParallel)
library(terra)
library(raster)

# Read in the model
m.spatial = readRDS(here::here("Data_Processing", "jsdm_hmsc_gippsland.RDS"))

# Read in the predictors
covariate_stack = readRDS("Data_Clean/covariate_stack_updated.RData")
base_raster_mask = covariate_stack$base_raster
base_raster_mask[base_raster_mask==0] <- NA
covariate_stack_masked = raster::mask(covariate_stack, base_raster_mask)
data = read.csv(here::here("Data_Clean", "pa_covariates_sark_lfp_stq_allspp.csv"))
data = data %>% 
  filter(Program %in% c("W Tree", "Southern Ark", "Quoll")) %>%
  drop_na()

grad = constructGradient(m.spatial, focalVariable = "bait_intensity")
pred = predict(m.spatial, Gradient = grad, expected = TRUE)
plotGradient(m.spatial, grad, predY = pred, measure = "Y", index=9, q=c(0.05,0.5, 0.95))

plotGradient(m.spatial, grad, predY = pred, measure = "S", q=c(0.05,0.5, 0.95))

##########################
#### PRE FIRE RASTERS ####
##########################

cov.stack = subset(covariate_stack_masked, c('easting', 'northing', "tsf_startsurvey_burn", "rain_deciles_2015", "mean_repeat_burn",
                                             "bait_intensity", "harvest_prop_40", "mean_annual_rainfall", "burnability",
                                             "prop_burnt_twice", "topographic_wetness_index", "prop_recentburn", "clay", "carbon"))

covras.df = data.frame(rasterToPoints(cov.stack))

names(covras.df) <- c('x', 'y', "easting", "northing", "tsf_burn", "rain_deciles_2015", "mean_repeat_burn",
                      "bait_intensity", "harvest_prop_40", "mean_rainfall", "burnability", "prop_burnt_twice", "twi",
                      "prop_recent_burn", "clay", "carbon")

covras.df$easting = standardize2match(covras.df$easting, data$Easting)
covras.df$northing = standardize2match(covras.df$northing, data$Northing)
covras.df$tsf_burn = standardize2match(covras.df$tsf_burn, data$tsf_burn)
covras.df$rain_deciles_2015 = standardize2match(covras.df$rain_deciles_2015, data$rain_deciles_2015)
covras.df$mean_repeat_burn = standardize2match(covras.df$mean_repeat_burn, data$mean_repeat_burn)
covras.df$bait_intensity = standardize2match(covras.df$bait_intensity, data$bait_intensity)
covras.df$harvest_prop_40 = standardize2match(covras.df$harvest_prop_40, data$harvest_prop_40)
covras.df$mean_rainfall = standardize2match(covras.df$mean_rainfall, data$mean_rainfall)
covras.df$burnability = standardize2match(covras.df$burnability, data$burnability)
covras.df$twi = standardize2match(covras.df$twi, data$twi)
covras.df$prop_burnt_twice = standardize2match(covras.df$prop_burnt_twice, data$prop_burnt_twice)
covras.df$prop_recent_burn = standardize2match(covras.df$prop_recent_burn, data$prop_recent_burn)
covras.df$clay = standardize2match(covras.df$clay, data$clay)
covras.df$carbon = standardize2match(covras.df$carbon, data$carbon)

nd = covras.df %>% drop_na()
SP = data.frame(x = nd$easting,
                y = nd$northing)

SP.lookup = data.frame(x = nd$easting, y = nd$northing,
                       easting = nd$x, northing = nd$y)

# Make prediction grids
xy.grid = data.frame(SP)
XData.grid = data.frame(nd)

pred.list = list()
# Split the data into chunks
chunk_size = 5000
chunks <- split(seq_len(nrow(xy.grid)), ceiling(seq_len(nrow(xy.grid))/chunk_size))

# Function to perform the prediction for a chunk
# Function to calculate quantiles for each species across all slots
calculate_quantiles <- function(predictions_list, lower_quantile, upper_quantile) {
  # Extract the number of species and sites from the first matrix
  num_species <- ncol(predictions_list[[1]])
  num_sites <- nrow(predictions_list[[1]])
  # Create an array to store quantiles
  quantiles_array <- array(NA, dim = c(num_species, num_sites, 2))
  # Loop over species and sites
  for (species in 1:num_species) {
    for (site in 1:num_sites) {
      # Extract values for the specific species and site across all slots
      values <- sapply(predictions_list, function(matrix) matrix[site, species])
      # Calculate quantiles
      quantiles <- quantile(values, c(lower_quantile, upper_quantile))
      # Store quantiles in the array
      quantiles_array[species, site, ] <- quantiles
    }
  }
  
  # Return the quantiles array
  return(quantiles_array)
}


predict_mean_occurrence <- function(coords, covs, model.object) {
  predict.gradient =Hmsc::prepareGradient(model.object, XDataNew = covs, 
                                          sDataNew = list(sample=coords))
  predY = predict(model.object, Gradient = predict.gradient, 
                  expected = TRUE, predictEtaMean = TRUE)
  # Calculate mean and upper and lower 95% quantiles
  EpredY = Reduce("+", predY) / length(predY)
  aP = array(unlist(predY), dim = c(nrow(coords), 40, length(predY)))
  quantiles = apply(aP, c(1,2), quantile, probs = c(0.025,0.975))
  colnames(quantiles[1,,]) = paste0(colnames(EpredY),"_LCI")
  colnames(quantiles[2,,]) = paste0(colnames(EpredY),"_UCI")
  prediction.out = cbind(xy.grid_chunk, EpredY, quantiles[1,,], quantiles[2,,])
}

predictions = data.frame()
for(i in 1:length(chunks)){
  chunk = chunks[[i]]
  xy.grid_chunk = xy.grid[chunk, ]
  XData.grid_chunk = XData.grid[chunk, ]
  pred.out = predict_mean_occurrence(coords = xy.grid_chunk, covs = XData.grid_chunk, model.object = m.spatial)
  predictions = rbind(predictions, pred.out)
  print(paste0("Chunk number ", i, " done at ", Sys.time()))
} 

# pred.list <- foreach(i = chunks, .options.snow = opts, .combine = rbind) %dopar% {
#   xy.grid_chunk = xy.grid[i, ]
#   XData.grid_chunk = XData.grid[i, ]
#   predictions = predict_mean_occurrence(coords = xy.grid_chunk, covs = XData.grid_chunk, model.object = m.spatial)
#   predictions
# }
# 
# # Close the parallel cluster
# stopCluster(cl)
names.cols = c("x", "y", m.spatial$spNames, paste0(m.spatial$spNames,"_LCI"), paste0(m.spatial$spNames,"_UCI"))
names(predictions) <- names.cols
saveRDS(predictions, file = "Data_Clean/hmsc_prefire_spatialpredictions_uncertainty_gippsland.RDS")

# Check the plot
ggplot(data = predictions) +
  geom_tile(aes(x=x, y=y, fill=Black.Wallaby)) + 
  ggtitle(expression(italic("Pr(Occurrence)")))+ 
  scale_fill_viridis_c() + coord_equal()


###########################
#### POST FIRE RASTERS ####
###########################
cov.stack = subset(covariate_stack_masked, c('easting', 'northing', 
                                             "tsf_startsurvey_burn_post1920", "rain_deciles_2015", 
                                             "mean_repeat_burn_post1920",
                                             "bait_intensity", "harvest_prop_40", "mean_annual_rainfall", "burnability",
                                             "prop_burnt_twice_post1920", "topographic_wetness_index", 
                                             "prop_recentburn_post1920", 
                                             "clay", "carbon"))

covras.df = data.frame(rasterToPoints(cov.stack))

names(covras.df) <- c('x', 'y', "easting", "northing", "tsf_burn", "rain_deciles_2015", "mean_repeat_burn",
                      "bait_intensity", "harvest_prop_40", "mean_rainfall", "burnability", "prop_burnt_twice", "twi",
                      "prop_recent_burn", "clay", "carbon")

covras.df$easting = standardize2match(covras.df$easting, data$Easting)
covras.df$northing = standardize2match(covras.df$northing, data$Northing)
covras.df$tsf_burn = standardize2match(covras.df$tsf_burn, data$tsf_burn)
covras.df$rain_deciles_2015 = standardize2match(covras.df$rain_deciles_2015, data$rain_deciles_2015)
covras.df$mean_repeat_burn = standardize2match(covras.df$mean_repeat_burn, data$mean_repeat_burn)
covras.df$bait_intensity = standardize2match(covras.df$bait_intensity, data$bait_intensity)
covras.df$harvest_prop_40 = standardize2match(covras.df$harvest_prop_40, data$harvest_prop_40)
covras.df$mean_rainfall = standardize2match(covras.df$mean_rainfall, data$mean_rainfall)
covras.df$burnability = standardize2match(covras.df$burnability, data$burnability)
covras.df$twi = standardize2match(covras.df$twi, data$twi)
covras.df$prop_burnt_twice = standardize2match(covras.df$prop_burnt_twice, data$prop_burnt_twice)
covras.df$prop_recent_burn = standardize2match(covras.df$prop_recent_burn, data$prop_recent_burn)
covras.df$clay = standardize2match(covras.df$clay, data$clay)
covras.df$carbon = standardize2match(covras.df$carbon, data$carbon)

nd = covras.df %>% drop_na()
SP = data.frame(x = nd$easting,
                y = nd$northing)

# Make prediction grids
xy.grid = data.frame(SP)
XData.grid = data.frame(nd)

pred.list = list()
# Split the data into chunks
chunk_size = 5000
chunks <- split(seq_len(nrow(xy.grid)), ceiling(seq_len(nrow(xy.grid))/chunk_size))

# Function to calculate quantiles for each species across all slots
calculate_quantiles <- function(predictions_list, lower_quantile, upper_quantile) {
  # Extract the number of species and sites from the first matrix
  num_species <- ncol(predictions_list[[1]])
  num_sites <- nrow(predictions_list[[1]])
  # Create an array to store quantiles
  quantiles_array <- array(NA, dim = c(num_species, num_sites, 2))
  # Loop over species and sites
  for (species in 1:num_species) {
    for (site in 1:num_sites) {
      # Extract values for the specific species and site across all slots
      values <- sapply(predictions_list, function(matrix) matrix[site, species])
      # Calculate quantiles
      quantiles <- quantile(values, c(lower_quantile, upper_quantile))
      # Store quantiles in the array
      quantiles_array[species, site, ] <- quantiles
    }
  }
  
  # Return the quantiles array
  return(quantiles_array)
}


predict_mean_occurrence <- function(coords, covs, model.object) {
  predict.gradient =Hmsc::prepareGradient(model.object, XDataNew = covs, 
                                          sDataNew = list(sample=coords))
  predY = predict(model.object, Gradient = predict.gradient, 
                  expected = TRUE, predictEtaMean = TRUE)
  # Calculate mean and upper and lower 95% quantiles
  EpredY = Reduce("+", predY) / length(predY)
  aP = array(unlist(predY), dim = c(nrow(coords), 40, length(predY)))
  quantiles = apply(aP, c(1,2), quantile, probs = c(0.025,0.975))
  colnames(quantiles[1,,]) = paste0(colnames(EpredY),"_LCI")
  colnames(quantiles[2,,]) = paste0(colnames(EpredY),"_UCI")
  prediction.out = cbind(xy.grid_chunk, EpredY, quantiles[1,,], quantiles[2,,])
}


# Apply the function in parallel
predictions = data.frame()
for(i in 1:length(chunks)){
  chunk = chunks[[i]]
  xy.grid_chunk = xy.grid[chunk, ]
  XData.grid_chunk = XData.grid[chunk, ]
  pred.out = predict_mean_occurrence(coords = xy.grid_chunk, covs = XData.grid_chunk, model.object = m.spatial)
  predictions = rbind(predictions, pred.out)
  print(paste0("Chunk number ", i, " done at ", Sys.time()))
} # started at 9.25am on 26/01/2024 - 

# Save
names.cols = c("x", "y", m.spatial$spNames, paste0(m.spatial$spNames,"_LCI"), paste0(m.spatial$spNames,"_UCI"))
names(predictions) <- names.cols
saveRDS(predictions, file = "Data_Clean/hmsc_postfire_spatialpredictions_gippsland.RDS")

####################################
#### 10 YEARS POST FIRE RASTERS ####
####################################
cov.stack = subset(covariate_stack_masked, c('easting', 'northing', 
                                             "tsf_burn_10yrspost1920", "rain_deciles_2015", 
                                             "mean_repeat_burn_10yrspost1920",
                                             "bait_intensity", "harvest_prop_40", "mean_annual_rainfall", "burnability",
                                             "prop_burnt_twice_post1920", "topographic_wetness_index", 
                                             "prop_recentburn_post1920", 
                                             "clay", "carbon"))

covras.df = data.frame(rasterToPoints(cov.stack))

names(covras.df) <- c('x', 'y', "easting", "northing", "tsf_burn", "rain_deciles_2015", "mean_repeat_burn",
                      "bait_intensity", "harvest_prop_40", "mean_rainfall", "burnability", "prop_burnt_twice", "twi",
                      "prop_recent_burn", "clay", "carbon")

covras.df$easting = standardize2match(covras.df$easting, data$Easting)
covras.df$northing = standardize2match(covras.df$northing, data$Northing)
covras.df$tsf_burn = standardize2match(covras.df$tsf_burn, data$tsf_burn)
covras.df$rain_deciles_2015 = standardize2match(covras.df$rain_deciles_2015, data$rain_deciles_2015)
covras.df$mean_repeat_burn = standardize2match(covras.df$mean_repeat_burn, data$mean_repeat_burn)
covras.df$bait_intensity = standardize2match(covras.df$bait_intensity, data$bait_intensity)
covras.df$harvest_prop_40 = standardize2match(covras.df$harvest_prop_40, data$harvest_prop_40)
covras.df$mean_rainfall = standardize2match(covras.df$mean_rainfall, data$mean_rainfall)
covras.df$burnability = standardize2match(covras.df$burnability, data$burnability)
covras.df$twi = standardize2match(covras.df$twi, data$twi)
covras.df$prop_burnt_twice = standardize2match(covras.df$prop_burnt_twice, data$prop_burnt_twice)
covras.df$prop_recent_burn = standardize2match(covras.df$prop_recent_burn, data$prop_recent_burn)
covras.df$clay = standardize2match(covras.df$clay, data$clay)
covras.df$carbon = standardize2match(covras.df$carbon, data$carbon)

nd = covras.df %>% drop_na()
SP = data.frame(x = nd$easting,
                y = nd$northing)

# Make prediction grids
xy.grid = data.frame(SP)
XData.grid = data.frame(nd)

pred.list = list()
# Split the data into chunks
chunk_size = 5000
chunks <- split(seq_len(nrow(xy.grid)), ceiling(seq_len(nrow(xy.grid))/chunk_size))

predict_mean_occurrence <- function(coords, covs, model.object) {
  predict.gradient =Hmsc::prepareGradient(model.object, XDataNew = covs, 
                                          sDataNew = list(sample=coords))
  predY = predict(model.object, Gradient = predict.gradient, 
                  expected = TRUE, predictEtaMean = TRUE)
  # Calculate mean and upper and lower 95% quantiles
  EpredY = Reduce("+", predY) / length(predY)
  aP = array(unlist(predY), dim = c(nrow(coords), 40, length(predY)))
  quantiles = apply(aP, c(1,2), quantile, probs = c(0.025,0.975))
  colnames(quantiles[1,,]) = paste0(colnames(EpredY),"_LCI")
  colnames(quantiles[2,,]) = paste0(colnames(EpredY),"_UCI")
  prediction.out = cbind(xy.grid_chunk, EpredY, quantiles[1,,], quantiles[2,,])
}

# Apply the function in parallel
predictions = data.frame()
for(i in 1:length(chunks)){
  chunk = chunks[[i]]
  xy.grid_chunk = xy.grid[chunk, ]
  XData.grid_chunk = XData.grid[chunk, ]
  pred.out = predict_mean_occurrence(coords = xy.grid_chunk, covs = XData.grid_chunk, model.object = m.spatial)
  predictions = rbind(predictions, pred.out)
  print(paste0("Chunk number ", i, " done at ", Sys.time()))
} # started at 9.25am on 26/01/2024 - 


names.cols = c("x", "y", m.spatial$spNames, paste0(m.spatial$spNames,"_LCI"), paste0(m.spatial$spNames,"_UCI"))
names(predictions) <- names.cols
saveRDS(predictions, file = "Data_Clean/hmsc_10yrspostfire_spatialpredictions_gippsland.RDS")

#########################
#### Save as rasters ####
#########################

preds = list.files(path="Data_Clean", pattern = "spatialpredictions_gippsland.RDS", full=T)
species = m.spatial$spNames

SP.lookup = data.frame(x = nd$easting, y = nd$northing,
                       easting = nd$x, northing = nd$y)

# Pre Fire
pred.pre = readRDS(preds[3])
pred.pre = left_join(SP.lookup, pred.pre, by = c("x"="x", "y"="y"))
ras.out.pre = list()
for (i in 1:length(species)){
  spp = species[i]
  sdm = rasterFromXYZ(dplyr::select(pred.pre, c(easting,northing,spp)), res=res(cov.stack), crs=crs(cov.stack))
  names(sdm) <- spp
  ras.out.pre[[i]] <- sdm
  writeRaster(sdm, paste0("Data_Clean/SDMs/prefire_", spp, ".tif"), overwrite = TRUE)
}
ras.stack.pre = stack(ras.out.pre)
plot(ras.stack.pre)

# Post Fire
pred.post = readRDS(preds[2])
pred.post = left_join(SP.lookup, pred.post, by = c("x"="x", "y"="y"))
ras.out.post = list()
for (i in 1:length(species)){
  spp = species[i]
  sdm = rasterFromXYZ(dplyr::select(pred.post, c(easting,northing,spp)), res=res(cov.stack), crs=crs(cov.stack))
  names(sdm) <- spp
  ras.out.post[[i]] <- sdm
  writeRaster(sdm, paste0("Data_Clean/SDMs/postfire_", spp, ".tif"), overwrite = TRUE)
}
ras.stack.post = stack(ras.out.post)
plot(ras.stack.post)

# Ten Years Post Fire
pred.10post = readRDS(preds[1])
pred.10post = left_join(SP.lookup, pred.10post, by = c("x"="x", "y"="y"))
ras.out.10post = list()
for (i in 1:length(species)){
  spp = species[i]
  sdm = rasterFromXYZ(dplyr::select(pred.10post, c(easting,northing,spp)), res=res(cov.stack), crs=crs(cov.stack))
  names(sdm) <- spp
  ras.out.10post[[i]] <- sdm
  writeRaster(sdm, paste0("Data_Clean/SDMs/10yearspostfire_", spp, ".tif"), overwrite = TRUE)
}
ras.stack.10post = stack(ras.out.10post)
plot(ras.stack.10post)
