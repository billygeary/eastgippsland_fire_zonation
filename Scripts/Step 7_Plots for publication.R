## Make plots for paper
library(terra)
library(tmap)
library(tidyverse)
library(sf)
library(Hmsc)
library(ggplot2)
library(cowplot)
library(tidybayes)

######################
#### JSDM OUTPUTS ####
######################
m.spatial = readRDS("Data_Processing/jsdm_hmsc_gippsland.RDS")
mpost = convertToCodaObject(m.spatial)

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
                                 "Bird", "Bird"),
                      SpeciesLab = c("Agile Antechinus", "Australian Magpie", "Australian Owlet-nightjar", "Australian Raven",
                         "Bassian Thrush", "Beautiful Firetail", "Black Rat", "Swamp Wallaby",
                         "Blotched Blue tongued Lizard", "Brown Goshawk", "Brush Bronzewing", "Bush Rat",
                         "Feral Cat", "Cattle", "Common Blackbird", "Common Bronzewing",
                         "Common Brushtail Possum", "Common Ringtail Possum", "Common Wombat", "Crescent Honeyeater",
                         "Crimson Rosella", "Deer", "Dingo", "Dusky Antechinus",
                         "Dusky Woodswallow", "Eastern Brown Snake", "Eastern Grey Kangaroo", "Eastern Pygmy possum",
                         "Eastern Whipbird", "Eastern Yellow Robin", "Emu", "European Hare",
                         "European Rabbit", "Fallow Deer", "Flame Robin", "Grey Butcherbird",
                         "Grey Currawong", "Grey Fantail", "Grey Shrike thrush", "Hooded Robin",
                         "Feral Horse", "Koala", "Lace Monitor", "Laughing Kookaburra",
                         "Long-footed Potoroo", "Long-nosed Bandicoot", "Long-nosed Potoroo", "Masked Owl",
                         "Mountain Brushtail Possum", "Mountain Dragon", "New Holland Honeyeater", "Olive Whistler",
                         "Pied Currawong", "Feral Pig", "Pilotbird", "Red-bellied Black Snake",
                         "Red-browed Finch", "Red-browed Treecreeper", "Red-necked Wallaby", "Red Fox",
                         "Red Wattlebird", "Rufous Fantail", "Sambar Deer", "Satin Bowerbird",
                         "Scarlet Robin", "Short-beaked Echidna", "Southern Brown Bandicoot", "Spot-tailed Quoll",
                         "Spotted Quail-thrush", "Sugar Glider", "Superb Fairy wren", "Superb Lyrebird",
                         "Swamp Rat", "Tawny Frogmouth", "Tiger Snake", "Water Rat",
                         "White-browed Scrubwren", "White-eared Honeyeater", "White-footed Dunnart", "White-throated Nightjar",
                         "White-throated Treecreeper", "White-winged Chough", "Wonga Pigeon", "Yellow-faced Honeyeater",
                         "Yellow-tailed Black Cockatoo", "Yellow-tufted Honeyeater"))

covariates = data.frame(Covariate = c("(Intercept)","I(tsf_burn^2)","bait_intensity","carbon","clay", "harvest_prop_40",
                                      "mean_rainfall", "mean_repeat_burn","tsf_burn", "twi"),
                        CovLab = c("Intercept", "Time Since Fire^2", "Baiting Intensity", "Carbon", "Clay", "Timber Harvesting", 
                                   "Mean Annual Rainfall", "Mean Repeat Burn", "Time Since Fire", "Topographic Wetness Index"))


betas = mpost$Beta %>% as.mcmc.list() %>% tidy_draws() %>% pivot_longer(cols=4:403, names_to = "Variable", values_to = "Value") %>%
  separate(col=Variable, into =c("Covariate", "Species"), sep=",") %>%
  mutate(Covariate = str_extract(Covariate, pattern = "B\\[([^\\s]+)")) %>%
  mutate(Covariate = gsub("B\\[","",Covariate)) %>%
  mutate(Species = sub(" ", "", Species)) %>% mutate(Species = sub("\\s.*", "",Species)) %>%
  group_by(Species, Covariate) %>%
  summarise(Mean = mean(Value), 
            Median = median(Value), 
            LCI = quantile(Value, probs = 0.05),
            UCI = quantile(Value, probs = 0.95)) %>%
  left_join(spp.names, by=c("Species" = "species")) %>% left_join(covariates, by = "Covariate")

# Create a new column indicating whether CIs cross 0
betas$CrossZero <- ifelse(betas$LCI > 0 | betas$UCI < 0, "No", "Yes")

# Alphabetical order
betas$SpeciesLab <- factor(betas$SpeciesLab, levels = rev(unique(betas$SpeciesLab)))

mammals = betas %>% filter(group %in% c("Mammal", "Reptile")) %>% filter(Covariate != "(Intercept)") %>%
  ggplot(aes(x=Mean, y = SpeciesLab)) + 
  geom_point(aes(fill = CrossZero, shape = CrossZero), size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("Yes" = "white", "No" = "gray")) +
  scale_shape_manual(values = c("Yes" = 1, "No" = 16)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Coefficient + CIs",y="") +
  theme_cowplot() + facet_wrap(~CovLab, scales="free_x", ncol=3)

ggsave(plot = mammals, filename = "Outputs/mammals_coefs.pdf", width = 6, height=8, dpi=300, scale=1.6)

birds = betas %>% filter(group=="Bird") %>% filter(Covariate != "(Intercept)") %>%
  ggplot(aes(x=Mean, y = SpeciesLab)) + 
  geom_point(aes(fill = CrossZero, shape = CrossZero), size = 2, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("Yes" = "white", "No" = "gray")) +
  scale_shape_manual(values = c("Yes" = 1, "No" = 16)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Coefficient + CIs",y="") +
  theme_cowplot() + facet_wrap(~CovLab, scales="free_x", ncol=3)

ggsave(plot = birds, filename = "Outputs/birds_coefs.pdf", width = 6, height=8, dpi=300, scale=1.6)


invasive = betas %>% filter(group=="Invasive") %>% filter(Covariate != "(Intercept)") %>%
  ggplot(aes(x=Mean, y = SpeciesLab)) + 
  geom_point(aes(fill = CrossZero, shape = CrossZero), size = 2, color = "black",  show.legend = FALSE) +
  scale_fill_manual(values = c("Yes" = "white", "No" = "gray")) +
  scale_shape_manual(values = c("Yes" = 1, "No" = 16)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Coefficient + CIs",y="") +
  theme_cowplot() + facet_wrap(~CovLab, scales="free_x", ncol=3)

ggsave(plot = invasive, filename = "Outputs/invasive_coefs.pdf", width = 6, height=5, dpi=300, scale=1.6)

#########################
#### SPATIAL FIGURES ####
#########################

col_palette <- viridis::viridis(n=100)
change_pal = RColorBrewer::brewer.pal(n=11, "RdBu")

#### Figure 4,5: SDM changes #### 
sl_pre = terra::rast(here::here("Data_Clean", "SDMs", "prefire_Superb.Lyrebird.tif"))
sl_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Superb.Lyrebird.tif"))
sl_tenpost = terra::rast(here::here("Data_Clean", "SDMs", "10yearspostfire_Superb.Lyrebird.tif"))
sl_change = sl_post-sl_pre

lnb_pre = terra::rast(here::here("Data_Clean", "SDMs", "prefire_Long.nosed.Bandicoot.tif"))
lnb_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Long.nosed.Bandicoot.tif"))
lnb_change = lnb_post - lnb_pre

lfb_pre = terra::rast(here::here("Data_Clean", "SDMs", "prefire_Long.footed.Potoroo.tif"))
lfb_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Long.footed.Potoroo.tif"))
lfb_change = lfb_post - lfb_pre

bt_pre = terra::rast(here::here("Data_Clean", "SDMs", "prefire_Bassian.Thrush.tif"))
bt_post = terra::rast(here::here("Data_Clean", "SDMs", "postfire_Bassian.Thrush.tif"))
bt_change = bt_post - bt_pre

fire_poly = read_sf(here::here("Data_Raw", "spatial", "FH_2020_burn.shp")) %>% group_by() %>% summarise() %>% st_transform(crs=crs(bt_change))
fire_poly = fire_poly %>% terra::vect()
fire_poly_crop = terra::crop(fire_poly,bt_change) %>% st_as_sf()

# make some bbox magic
bbox_new <- st_bbox(bt_change) # current bounding box
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values
bbox_new[3] <- bbox_new[3] + (0.2 * xrange) # xmax - right
bbox_new[4] <- bbox_new[4] + (0.2 * yrange) # ymax - top
bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon


occurrence_map = function(sdm, label){
  sdm_map = tm_shape(sdm, bbox = bbox_new, raster.downsample = FALSE) +
    tm_raster(palette = col_palette, style='cont', 
              title = "Pr(Occurrence)", 
              legend.is.portrait = FALSE, 
              breaks=c(min(values(sdm), na.rm=TRUE), max(values(sdm), na.rm=TRUE)))+ 
    tm_layout(title = paste0(label),
              title.position = c('left', 'top'),
              legend.position=c("right","top"),
              legend.outside=FALSE, 
              #legend.width = 0.5, 
              frame = FALSE)
  return(sdm_map)
}

change_map = function(change_sdm, label){
  change_out = tm_shape(change_sdm, bbox = bbox_new, raster.downsample = FALSE) +
    tm_raster(palette = change_pal, style='cont', 
              title = "Change in Pr(Occurrence)", 
              legend.is.portrait = FALSE, 
              breaks=c(min(values(change_sdm), na.rm=TRUE), max(values(change_sdm), na.rm=TRUE)))+ 
    tm_shape(fire_poly_crop) + tm_borders() +
    tm_layout(legend.position=c("right","top"),
              title = paste0(label),
              title.position = c('left', 'top'),
              legend.outside=FALSE, 
              legend.width = 0.6, 
              frame = FALSE)
  return(change_out)
}

sl_pre_map = occurrence_map(sl_pre, "a)")
sl_change_map = change_map(sl_change, "b)")
lnb_pre_map = occurrence_map(lnb_pre, "c)")
lnb_change_map = change_map(lnb_change, "d)")
bt_pre_map = occurrence_map(bt_pre, "e)")
bt_change_map = change_map(bt_change, "f)")

sdm_plots = tmap_arrange(sl_pre_map, sl_change_map,
                        lnb_pre_map, lnb_change_map,
                        bt_pre_map, bt_change_map,
                        ncol=2)

tmap_save(sdm_plots, filename = here::here("Outputs", "sdmchanges.pdf"), dpi=300, width=6, height=9)

#### Figure 6: Zonation outputs ####
z_pre = terra::rast(here::here("Data_Clean", "z_prefire", "rankmap.tif"))
z_post = terra::rast(here::here("Data_Clean", "z_postfire", "rankmap.tif"))
z_tenpost= terra::rast(here::here("Data_Clean", "z_10yearspostfire", "rankmap.tif")) 

z_map = function(z_raster, label, year.lab){
  z_map_out = tm_shape(z_raster, bbox = bbox_new, raster.downsample = FALSE) +
    tm_raster(palette = col_palette, style='pretty', n=5, 
              title = paste0(year.lab, " Rank"), 
              legend.is.portrait = FALSE, 
              breaks=c(min(values(z_raster), na.rm=TRUE), max(values(z_raster), na.rm=TRUE)))+ 
    tm_layout(title = paste0(label),
              title.position = c('left', 'top'),
              legend.position=c("right","top"),
              legend.outside=FALSE, 
              legend.width = 0.5, 
              frame = FALSE)
  return(z_map_out)
}

z_plots = tmap_arrange(z_map(z_pre, "a)", "2017"),
                       z_map(z_post, "b)", "2022"),
                       z_map(z_tenpost, "c)", "2030"),
                       ncol=3)

z_plots

tmap_save(z_plots, filename = here::here("Outputs", "zonation_maps.pdf"), dpi=300, width=12, height=4)

z_plot_2017 = z_map(z_pre, " ", "2017")
z_plot_2017
tmap_save(z_plot_2017, filename = here::here("Outputs", "zonation_map_2017.pdf"), dpi=300, width=6, height=4)



#### Figure 7: Zonation Change in Rank ####
library(terra)
z_pre_top = z_pre 
z_pre_top=ifel(z_pre_top>0.8,1,0)
z_post_top = z_post
z_post_top=ifel(z_post_top>0.8,2,0)
z_tenpost_top = z_tenpost
z_tenpost_top=ifel(z_tenpost_top>0.8,3,0)

# 3 = top rank in both, 2 = top rank post fire only, 1 = top rank in pre fire only
z_change_top = as.factor(z_pre_top + z_post_top)
z_tenchange_top = as.factor(z_pre_top + z_tenpost_top)

bbox_new <- st_bbox(bt_change) # current bounding box
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values
bbox_new[3] <- bbox_new[3] + (0.1 * xrange) # xmax - right
bbox_new[4] <- bbox_new[4] + (0.1 * yrange) # ymax - top
bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon


z_change_map = function(change_z, label, year.lab){
  change_out = tm_shape(change_z, bbox = bbox_new, raster.downsample = FALSE) +
    tm_raster(palette = (pal=c("lightgrey", "#443A83FF", "#21908CFF", "#FDE725FF")), style='cat', legend.show=FALSE)+ 
    tm_shape(fire_poly_crop) + tm_borders() +
    tm_add_legend(type = "symbol", shape = 15,
                  labels = c("Never High Rank", "No Longer High Rank", "Became High Rank", "High Rank in Both"),
                  col = c("lightgrey","#443A83FF", "#21908CFF", "#FDE725FF"),
                  border.lwd = 0.5,
                  title = paste0(year.lab, "Change"), is.portrait = TRUE) + 
    tm_layout(title = paste0(label),
              title.position = c('left', 'top'),
              legend.position=c("right", "top"),legend.outside=FALSE, 
              #legend.width = 0.6, 
              frame=FALSE)
  return(change_out)
}

z_change = z_change_map(z_change_top, "a)", "2017-2022 ")
z_tenchange = z_change_map(z_tenchange_top, "b)", "2017-2030 ")

z_change_plots = tmap_arrange(z_change, z_tenchange)

tmap_save(z_change_plots, filename = here::here("Outputs", "zonation_changes.pdf"), dpi=300, width=12, height=6)

z_plot_2017 = z_map(z_pre, "a)", "2017")
z_change = z_change_map(z_change_top, "b)", "2017-2022 ")
z_tenchange = z_change_map(z_tenchange_top, "c)", "2017-2030 ")

z_change_plots_2017 = tmap_arrange(z_plot_2017, z_change, z_tenchange)
tmap_save(z_change_plots_2017, filename = here::here("Outputs", "zonation_2017_andchanges.pdf"), 
          dpi=300, width=12, height=4)


########################
#### Sankey Diagram ####
########################
z_pre_top = z_pre 
z_pre_top=ifel(z_pre_top>0.8,1,0)
z_post_top = z_post
z_post_top=ifel(z_post_top>0.8,10,0)
z_tenpost_top = z_tenpost
z_tenpost_top=ifel(z_tenpost_top>0.8,100,0)

sankey.data = data.frame("y2017"= values(z_pre_top),
                         "y2020" = values(z_post_top),
                         "y2030" = values(z_tenpost_top))

names(sankey.data) <- c("y2017", "y2020", "y2030")

sankey.data = sankey.data %>% drop_na() 

ncells = length(sankey.data$y2017)

sankey.data.t1 = sankey.data %>%
  group_by(y2017, y2020) %>% summarise(cells = n()/ncells)

sankey.data.t2 = sankey.data %>%
  group_by(y2020, y2030) %>% summarise(cells = n()/ncells)



library(networkD3)
nodes = data.frame("name" = 
                     c("Node A", # Node 0
                       "Node B", # Node 1
                       "Node C", # Node 2
                       "Node D"))# Node 3

links = as.data.frame(matrix(c(
  0, 1, 10, # Each row represents a link. The first number
  0, 2, 20, # represents the node being conntected from. 
  1, 3, 30, # the second number represents the node connected to.
  2, 3, 40),# The third number is the value of the node
  byrow = TRUE, ncol = 3))
names(links) = c("source", "target", "value")


sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30)



nodes = data.frame("name" = c("top2017", "low2017", "top2020", "low2020", "top2030", "low2030"))

links = data.frame(matrix(c(
  0, 2, 0.157, # Top to Top
  0, 3, 0.0432, # Top to Low
  1, 2, 0.0430, # Low to Top
  1, 3, 0.757, # Low to Low
  2, 4, 0.178,  # Top to Top
  2, 5, 0.0222,
  3, 4, 0.0222,
  3, 5, 0.778),# The third number is the value of the node
  byrow = TRUE, ncol = 3))
names(links) = c("source", "target", "value")

sankeyNetwork(Links = links, Nodes = nodes,
              Source = "source", Target = "target",
              Value = "value", NodeID = "name",
              fontSize= 12, nodeWidth = 30)



URL <- paste0('https://cdn.rawgit.com/christophergandrud/networkD3/',
              'master/JSONdata/energy.json')
energy <- jsonlite::fromJSON(URL)

a = energy$nodes
b = energy$links
############################################
#### Figure S2: Covariates used in JSDM ####
############################################

# Read in the covariates 
library(terra)
covariate_stack = readRDS("Data_Clean/covariate_stack_updated.RData")
base_raster_mask = covariate_stack$base_raster
base_raster_mask[base_raster_mask==0] <- NA
covariate_stack_masked = raster::mask(covariate_stack, base_raster_mask)

plot_stack = subset(covariate_stack_masked, c("easting","northing", "bait_intensity", "tsf_endsurvey_burn",
                                              "mean_repeat_burn", "harvest_prop_40", "mean_annual_rainfall", 
                                              "topographic_wetness_index", "clay", "carbon","sand"))

col_palette <- viridis::viridis(n=100)

bbox_new <- st_bbox(base_raster_mask) # current bounding box
xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values
bbox_new[3] <- bbox_new[3] + (0.2 * xrange) # xmax - right
bbox_new[4] <- bbox_new[4] + (0.2 * yrange) # ymax - top
bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon


# Create a map object with the raster stack

plot_covariate = function(covname, covlabel){
  plot.ras = plot_stack[[paste0(covname)]]
  plot.out = tm_shape(plot.ras, bbox = bbox_new) +
    tm_raster(palette = col_palette, style='cont', 
              title = paste0(covlabel), 
              legend.is.portrait = FALSE, 
              breaks=c(min(values(plot.ras), na.rm=TRUE), max(values(plot.ras), na.rm=TRUE))) + 
    tm_layout(legend.position=c("right","top"),legend.outside=FALSE, frame=FALSE)
  
  return(plot.out)
}

east.plot = plot_covariate("easting", "Easting")
north.plot = plot_covariate("northing", "Northing")
bait.plot = plot_covariate("bait_intensity", "Bait Intensity")
tsf.plot = plot_covariate("tsf_endsurvey_burn", "Time Since Fire")
repeatburn.plot = plot_covariate("mean_repeat_burn", "Fire Frequency")
harvest.plot = plot_covariate("harvest_prop_40", "Timber Harvesting")
rain.plot = plot_covariate("mean_annual_rainfall", "Mean Annual Rainfall")
twi.plot = plot_covariate("topographic_wetness_index", "TWI")
clay.plot = plot_covariate("clay", "Clay")
carbon.plot = plot_covariate("carbon", "Carbon")

cov_map = tmap_arrange(east.plot, north.plot, bait.plot, tsf.plot, repeatburn.plot, 
             harvest.plot, rain.plot, twi.plot, clay.plot, carbon.plot,
             ncol = 2)

tmap_save(cov_map, filename = here::here("Outputs", "covariate_maps.pdf"), dpi=300, width=6, height=12)



#############################################
#### Change in suitable habitat approach ####
#############################################

preds = list.files(path="Data_Clean", pattern = "spatialpredictions_uncertainty_gippsland.RDS", full=T)
preds
species = names(readRDS(preds[1]))
pre = readRDS(preds[3]) %>% colSums() %>% as.vector()
post = readRDS(preds[2]) %>% colSums() %>% as.vector()
tenpost = readRDS(preds[1]) %>% colSums() %>% as.vector()

sum_occ = data.frame(Species = species,
                     pre_occ = pre,
                     post_occ = post, 
                     tenpost_occ = tenpost) %>% 
  filter(!Species %in% c("x","y")) %>%
  pivot_longer(cols = 2:4, names_to = "TimeStep", values_to = "Change_Occ") %>%
  mutate(TimeStep = factor(TimeStep, levels = c("pre_occ", "post_occ", "tenpost_occ")))



#sum_occ$Change_Occ_Ha = sum_occ$Change_Occ*50^2 / 1000

sum_occ = sum_occ %>% separate(Species, sep="_", into = c("Species", "Measure"))
sum_occ$Measure = ifelse(is.na(sum_occ$Measure)==TRUE, "Mean", sum_occ$Measure)

sum_occ = left_join(sum_occ, spp.names, by = c("Species"="species"))

sum_occ = sum_occ %>% pivot_wider(id_cols = c(SpeciesLab, TimeStep, group), 
                                  names_from=Measure, 
                                  names_prefix = "Change_Occ_", 
                                  values_from = Change_Occ)
label_lookup = data.frame(TimeStep = c("pre_occ", "post_occ", "tenpost_occ"),
                          TimeLab = c(2017, 2022, 2030))
sum_occ = left_join(sum_occ, label_lookup)

occ_change = 
  filter(sum_occ, group=="Mammal") %>% 
  ggplot() + 
  geom_pointrange(aes(x = TimeLab, 
                      y = Change_Occ_Mean/1000, 
                      ymin = Change_Occ_LCI/1000, 
                      ymax = Change_Occ_UCI/1000, 
                      color = SpeciesLab)) + 
  geom_path(aes(x = TimeLab, 
                y = Change_Occ_Mean/1000)) + 
  facet_wrap(~SpeciesLab, scale="free_y") + 
  scale_x_continuous(breaks = c(2018, 2020, 2030)) +
  theme_cowplot() + theme(legend.position = "none") + 
  ylab("Occurence-weighted Habitat ('000s Ha)") + xlab("Year")

occ_change

occ_change = 
  filter(sum_occ, group=="Bird") %>% 
  ggplot() + 
  geom_pointrange(aes(x = TimeLab, 
                      y = Change_Occ_Mean/1000, 
                      ymin = Change_Occ_LCI/1000, 
                      ymax = Change_Occ_UCI/1000, 
                      color = SpeciesLab)) + 
  geom_path(aes(x = TimeLab, 
                y = Change_Occ_Mean/1000)) + 
  facet_wrap(~SpeciesLab, scale="free_y") + 
  scale_x_continuous(breaks = c(2018, 2020, 2030)) +
  theme_cowplot() + theme(legend.position = "none") + 
  ylab("Occurence-weighted Habitat ('000s Ha)") + xlab("Year")
occ_change

occ_change = 
  filter(sum_occ, group=="Invasive") %>% 
  ggplot() + 
  geom_pointrange(aes(x = TimeLab, 
                      y = Change_Occ_Mean/1000, 
                      ymin = Change_Occ_LCI/1000, 
                      ymax = Change_Occ_UCI/1000, 
                      color = SpeciesLab)) + 
  geom_path(aes(x = TimeLab, 
                y = Change_Occ_Mean/1000)) + 
  facet_wrap(~SpeciesLab, scale="free_y") + 
  scale_x_continuous(breaks = c(2018, 2020, 2030)) +
  theme_cowplot() + theme(legend.position = "none") + 
  ylab("Occurence-weighted Habitat ('000s Ha)") + xlab("Year")
occ_change

occ_change = 
  filter(sum_occ,SpeciesLab %in% c("Superb Lyrebird", "Long-nosed Bandicoot", "Bassian Thrush")) %>% 
  ggplot() + 
  geom_pointrange(aes(x = TimeLab, 
                      y = Change_Occ_Mean/1000, 
                      ymin = Change_Occ_LCI/1000, 
                      ymax = Change_Occ_UCI/1000, 
                      color = SpeciesLab)) + 
  geom_path(aes(x = TimeLab, 
                y = Change_Occ_Mean/1000)) + 
  facet_wrap(~SpeciesLab, scale="free_y") + 
  scale_x_continuous(breaks = c(2018, 2020, 2030)) +
  theme_cowplot() + theme(legend.position = "none") + 
  ylab("Occurence-weighted Habitat ('000s Ha)") + xlab("Year")
occ_change

###########################################

