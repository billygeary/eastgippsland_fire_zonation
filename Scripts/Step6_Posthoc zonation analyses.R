####################################### 
#### POST HOC ANALYSIS OF ZONATION ####
#######################################

# Create the dataset
z_post = terra::rast(here::here("Data_Clean", "z_postfire", "rankmap.tif"))
z_post_top = z_post
z_post_top=ifel(z_post_top>0.8,1,0)
names(z_post_top) <- "rankcat"

covariate_stack = readRDS("Data_Clean/covariate_stack_updated.RData")
base_raster_mask = covariate_stack$base_raster
base_raster_mask[base_raster_mask==0] <- NA
covariate_stack_masked = raster::mask(covariate_stack, base_raster_mask)

cat_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Cat.tif"))
sambar_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Sambar.tif"))
fox_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Red.Fox.tif"))

stack = c(z_post, z_post_top, cat_post, sambar_post, fox_post)
tsf = rast(covariate_stack_masked$tsf_startsurvey_burn_post1920); tsf = resample(tsf, stack, method="near")
bait = rast(covariate_stack_masked$bait_intensity); bait = resample(bait, stack, method="near")
repeat_burn = rast(covariate_stack_masked$mean_repeat_burn); repeat_burn = resample(repeat_burn,stack, method="near")

stack = c(stack, tsf, bait, repeat_burn)
names(stack) = c("rankmap", "rankcat","Cat","Sambar", "Red.Fox", "tsf", "bait_intensity", "repeat_burn")

sample = spatSample(stack, size = 10000, method="random", na.rm=TRUE, xy=TRUE)

# logit transform the z value
sample$rankmap_logit = log((sample$rankmap-0.0000001)/(1-(sample$rankmap-0.0000001)))

# Create autocovariate 
SP = data.frame(x = sample$x,
                y = sample$y)
SP = as.matrix(SP)

sample$autocov = spdep::autocov_dist(z= sample$rankmap, 
                              xy = SP, 
                              nbs = 10000, 
                              type = "inverse", 
                              zero.policy = NULL,
                              style = "B", 
                              longlat=NULL)


# Check correlations
cor.covs = cor(sample, use='complete.obs', method='pearson')
corrplot::corrplot(cor.covs)

# Model 
library(lme4)
library(ggeffects)
library(MuMIn)
library(AHMbook)
library(brms)
AHMbook::standardize()
sample$Cat_s = scale(sample$Cat)
sample$Red.Fox_s = scale(sample$Red.Fox)
sample$Sambar_s = scale(sample$Sambar)
sample$tsf_s = scale(sample$tsf)
sample$bait_s = scale(sample$bait_intensity)
sample$repeatburn_s = scale(sample$repeat_burn)
sample$autocov_s = scale(sample$autocov)

m.cubic = brm(rankmap ~ poly(Cat_s,3) + poly(Sambar_s,3) + poly(Red.Fox_s,3) + 
                poly(tsf_s,3) + poly(bait_s,3) + poly(repeatburn_s,3) + poly(autocov_s, 3),
              family = "beta",
              data = sample)


summary(m.cubic)

AIC(m,m.cubic)
summary(m.cubic)

brms::bayes_R2(m.cubic)
out.coefs = fixef(m.cubic)
write.csv(out.coefs, "Outputs/zonation_model_coefs.csv")

catpred = data.frame(ggpredict(m.cubic, "Cat_s"))
sambarpred = data.frame(ggpredict(m.cubic, "Sambar_s"))
foxpred = data.frame(ggpredict(m.cubic, "Red.Fox_s"))
tsfpred = data.frame(ggpredict(m.cubic, "tsf_s"))
baitpred = data.frame(ggpredict(m.cubic, "bait_s"))
repeatpred = data.frame(ggpredict(m.cubic, "repeatburn_s"))

pred.plot = function(pred.data,unscaled, xname, samples.name){
  pred.data[,"unscaled"] <- pred.data[,1] * sd(unscaled) + mean(unscaled)
  
  plot.out = ggplot() +  
    geom_ribbon(aes(x = pred.data[,"unscaled"], 
                    ymin = pred.data$conf.low, 
                    ymax = pred.data$conf.high), fill="#21908CFF", alpha=0.4) + 
    geom_line(aes(x=pred.data[,"unscaled"], y = pred.data$predicted)) + theme_cowplot() + xlab(paste(xname)) + ylab("Conservation Value")
              return(plot.out)
}


a=pred.plot(catpred,sample$Cat, "Pr(Feral Cat Occurrence)", "Cat")
b=pred.plot(sambarpred,sample$Sambar, "Pr(Sambar Deer Occurrence)", "Sambar")
c=pred.plot(foxpred,sample$Red.Fox, "Pr(Red Fox Occurrence)", "Red.Fox")
d=pred.plot(tsfpred,sample$tsf, "Time Since Fire", "tsf")
e=pred.plot(baitpred,sample$bait_intensity, "Baiting Intensity", "bait_intensity")
f=pred.plot(repeatpred, sample$repeat_burn, "Repeat Burns", "repeat_burn")

final.plot = plot_grid(d,e,f,a,b,c, labels = c("a)","b)","c)","d)", "e)","f)"))
ggsave(final.plot, filename = "Outputs/zonation_associations.pdf", width = 7, height = 4,scale=1.5, dpi=300)
          