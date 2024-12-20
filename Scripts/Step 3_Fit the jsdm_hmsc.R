library(dplyr)
library(Hmsc)
library(tidyr)
data = read.csv(here::here("Data_Clean", "pa_covariates_sark_lfp_stq_allspp.csv"))

spp.names = data.frame(species = c("Agile.Antechinus", "Australian.Magpie", "Australian.Owlet.nightjar",  "Australian.Raven", "Bassian.Thrush",              
                                   "Beautiful.Firetail", "Black.Rat", "Black.Wallaby","Blotched.Blue.tongued.Lizard", "Brown.Goshawk",               
                                   "Brush.Bronzewing", "Bush.Rat", "Cat","Cattle",                      
                                   "Common.Blackbird", "Common.Bronzewing", "Common.Brushtail.Possum", "Common.Ringtail.Possum", "Common.Wombat",               
                                   "Crescent.Honeyeater", "Crimson.Rosella", "Deer", "Dog", "Dusky.Antechinus",            
                                   "Dusky.Woodswallow", "Eastern.Brown.Snake", "Eastern.Grey.Kangaroo", "Eastern.Pygmy.possum", "Eastern.Whipbird",            
                                   "Eastern.Yellow.Robin", "Emu", "European.Hare", "European.Rabbit", "Fallow.Deer",                 
                                   "Flame.Robin","Grey.Butcherbird","Grey.Currawong","Grey.Fantail","Grey.Shrike.thrush",          
                                   "Hooded.Robin","Horse","Koala","Lace.Monitor","Laughing.Kookaburra",         
                                   "Long.footed.Potoroo","Long.nosed.Bandicoot","Long.nosed.Potoroo","Masked.Owl","Mountain.Brushtail.Possum",   
                                   "Mountain.Dragon","New.Holland.Honeyeater","Olive.Whistler","Pied.Currawong","Pig",                         
                                   "Pilotbird","Red.bellied.Black.Snake","Red.browed.Finch","Red.browed.Treecreeper","Red.necked.Wallaby",          
                                   "Red.Fox","Red.Wattlebird","Rufous.Fantail","Sambar","Satin.Bowerbird",             
                                   "Scarlet.Robin","Short.beaked.Echidna","Southern.Brown.Bandicoot","Spot.tailed.Quoll","Spotted.Quail.thrush",        
                                   "Sugar.Glider","Superb.Fairy.wren","Superb.Lyrebird","Swamp.Rat","Tawny.Frogmouth",             
                                   "Tiger.Snake","Water.Rat","White.browed.Scrubwren","White.eared.Honeyeater","White.footed.Dunnart",        
                                   "White.throated.Nightjar","White.throated.Treecreeper","White.winged.Chough","Wonga.Pigeon","Yellow.faced.Honeyeater",     
                                   "Yellow.tailed.Black.Cockatoo","Yellow.tufted.Honeyeater"), 
                       group = c("Mammal", "Bird","Bird","Bird","Bird",
                                 "Bird", "Invasive", "Mammal", "Reptile", "Bird",
                                 "Bird", "Mammal", "Invasive", "Invasive",
                                 "Bird","Bird", "Mammal", "Mammal", "Mammal",
                                 "Bird","Bird", "Invasive", "Mammal", "Mammal",
                                 "Bird", "Reptile", "Mammal", "Mammal", "Bird",
                                 "Bird","Bird", "Invasive", "Invasive", "Invasive",
                                 "Bird","Bird","Bird","Bird","Bird",
                                 "Bird", "Invasive", "Mammal", "Reptile", "Bird",
                                 "Mammal", "Mammal", "Mammal", "Bird", "Mammal",
                                 "Reptile", "Bird", "Bird", "Bird", "Invasive",
                                 "Bird", "Reptile", "Bird","Bird", "Mammal", 
                                 "Invasive", "Bird","Bird", "Invasive", "Bird",
                                 "Bird", "Mammal", "Mammal", "Mammal", "Bird",
                                 "Mammal", "Bird","Bird","Mammal", "Bird",
                                 "Reptile", "Mammal", "Bird","Bird", "Mammal",
                                 "Bird","Bird","Bird","Bird","Bird",
                                 "Bird", "Bird"))



data = data %>% 
  filter(Program %in% c("W Tree", "Southern Ark", "Quoll")) %>%
  drop_na() %>%
  group_by(Program, Easting, Northing, bait_intensity, rain_deciles_2015, 
           mean_rainfall, tsf_burn, twi, burnability, 
           mean_repeat_burn, harvest_prop_40, prop_recent_burn, 
           prop_burnt_twice, prop_longunburnt_burn, clay, sand, carbon) %>%
  summarise(across(spp.names$species, max)) %>% ungroup()

#### Check species being included
species = data %>% 
  dplyr::select(all_of(spp.names$species)) %>% 
  as.matrix() %>% colSums() %>% sort()
species
species.det = names(species[species>20])

spp.names = filter(spp.names, species %in% species.det)

## Prepare and scale the coveriates
covs = data %>% dplyr::select(bait_intensity, mean_rainfall, twi, rain_deciles_2015,
                              burnability, harvest_prop_40, mean_repeat_burn, 
                              tsf_burn, prop_recent_burn, prop_burnt_twice,
                              prop_longunburnt_burn, sand, clay, carbon) 
covs_scaled = as.matrix(data.frame(sapply(covs, function(x) {scale(x)})))

spatial = data %>% dplyr::select(Easting, Northing) 
spatial_scaled = as.matrix(data.frame(sapply(spatial, function(x){scale(x)})))

## Prepare species input data
species = data %>% 
  dplyr::select(all_of(spp.names$species)) %>% 
  as.matrix()

#####
model.formula = "~ tsf_burn + I(tsf_burn^2) + mean_repeat_burn + bait_intensity + harvest_prop_40 + mean_rainfall + twi + clay + carbon"

n = nrow(data)
xycoords = spatial_scaled
colnames(xycoords) = c("x-coordinate","y-coordinate")
rownames(xycoords) = 1:n
studyDesign = data.frame(sample = as.factor(1:n))
rL.spatial = HmscRandomLevel(sData = xycoords, sMethod = 'NNGP', nNeighbours = 10)
rL.spatial = setPriors(rL.spatial,nfMin=1,nfMax=1) #We limit the model to one latent variables for visum.spatial = Hmsc(Y=yprob, XData=XData, XFormula=~x1,
m.spatial = Hmsc(Y=species, 
                 XData = data.frame(covs_scaled), 
                 XFormula=as.formula(paste(model.formula)),
                 studyDesign=studyDesign, 
                 ranLevels=list("sample"=rL.spatial), 
                 distr="probit")

nChains = 4
thin = 100
samples = 250 # we want 1000 samples at the end. But the total iteraations is thin*samples
transient = 50*thin
verbose = 25*thin
start = Sys.time()
m.spatial = sampleMcmc(m.spatial, 
                       thin = thin, 
                       samples = samples, 
                       transient = transient,
                       nChains = nChains, 
                       verbose = verbose,
                       nParallel = nChains,
                       updater=list(GammaEta=FALSE))
end = Sys.time()
end-start

mpost = convertToCodaObject(m.spatial)

# Model convergance checks 
# Rhats
hist(gelman.diag(mpost$Beta, multivariate=FALSE)$psrf, main="psrf(beta)")
hist(gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf, main="psrf(gamma)")
hist(gelman.diag(mpost$V, multivariate=FALSE)$psrf, main="psrf(V)")
hist(gelman.diag(mpost$Sigma, multivariate=FALSE)$psrf, main="psrf(sigma)")
eta = as.mcmc(mpost$Eta)
hist(gelman.diag(eta[[1]], multivariate=FALSE)$psrf, main="psrf(Eta)")
lambda = as.mcmc(mpost$Lambda)
hist(gelman.diag(lambda[[1]], multivariate=FALSE)$psrf, main="psrf(Lambda)")
omega = as.mcmc(mpost$Omega)
hist(gelman.diag(omega[[1]], multivariate=FALSE)$psrf, main="psrf(Omega)")
alpha = as.mcmc(mpost$Alpha)
hist(gelman.diag(alpha[[1]], multivariate=FALSE)$psrf, main="psrf(Alpha)")
psi = as.mcmc(mpost$Psi)
hist(gelman.diag(psi[[1]], multivariate=FALSE)$psrf, main="psrf(Psi)")
delta = as.mcmc(mpost$Delta)
hist(gelman.diag(delta[[1]], multivariate=FALSE)$psrf, main="psrf(Delta)")

# ESS
hist(effectiveSize(mpost$Beta), main="ess(beta)")
hist(effectiveSize(mpost$Gamma), main="ess(gamma)")
hist(effectiveSize(mpost$V), main="ess(V)")
hist(effectiveSize(mpost$Sigma), main="ess(sigma)")
hist(effectiveSize(mpost$Eta), main="ess(Eta)")

## Model looks good, so lets save it. 

saveRDS(m.spatial, file = "Data_Processing/jsdm_hmsc_gippsland.RDS")

m.spatial = readRDS("Data_Processing/jsdm_hmsc_gippsland.RDS")

# Look at predictive capacity - RMSE and AUC for this 
# Not cross validated
preds = computePredictedValues(m.spatial)
model.fit = data.frame(evaluateModelFit(hM = m.spatial, predY = preds))
model.fit$Species = m.spatial$spNames

write.csv(model.fit, "Outputs/jsdm_hmsc_gippsland_modelfit.csv")

# Quick look at the betas and which ones have suppport 
postBeta = getPostEstimate(m.spatial, parName = "Beta", q = c(0.05, 0.95))
plotBeta(m.spatial, post = postBeta, param = "Support", supportLevel = 0.90)

## Cross Validation of model performance
## Allows us to filter which species have reasonable models or not. 
n.folds = 4
partition = createPartition(m.spatial, nfolds = n.folds, column = "sample") 

cvpreds.spatial = computePredictedValues(m.spatial, 
                                         partition=partition,
                                         nParallel = nChains, 
                                         updater=list(GammaEta=FALSE))

cvMF.spatial = data.frame(evaluateModelFit(hM=m.spatial, predY=cvpreds.spatial)) 
cvMF.spatial$Species = m.spatial$spNames
cvMF.spatial

## Save the cf outputs
write.csv(cvMF.spatial, "Outputs/jsdm_hmsc_gippsland_CV_modelfit.csv")



