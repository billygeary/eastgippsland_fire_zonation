## Create files for zonation
## Feature list file - pre fire
wd = getwd()
rasters = intersect(list.files("Data_Clean/SDMs", pattern = "\\.tif$", full=T),
                    list.files("Data_Clean/SDMs", pattern = "prefire", full=T))
name = gsub("Data_Clean/SDMs/prefire_","",rasters)
name = gsub(".tif","",name)

features = data.frame(filename = rasters, weight = 1, name = name)

features = features %>% 
  filter(!name %in% c("Black.Rat", "Cat", "European.Rabbit", "Fallow.Deer", "Horse", "Pig", "Red.Fox", "Sambar")) 

features$filename = paste0(wd,"/", features$filename)

write.table(features, "Scripts/Z_files/features_prefire.txt", row.names=FALSE)

## Feature list file - post fire
rasters = intersect(list.files("Data_Clean/SDMs", pattern = "\\.tif$", full=T),
                    list.files("Data_Clean/SDMs", pattern = "postfire", full=T))
name = gsub("Data_Clean/SDMs/postfire_","",rasters)
name = gsub(".tif","",name)

features = data.frame(filename = rasters, weight = 1, name = name)

features = features %>% filter(!grepl("10years",name)) %>%
  filter(!name %in% c("Black.Rat", "Cat", "European.Rabbit", "Fallow.Deer", "Horse", "Pig", "Red.Fox", "Sambar")) 
features$filename = paste0(wd,"/", features$filename)

write.table(features, "Scripts/Z_files/features_postfire.txt", row.names=FALSE)

## Feature list file - ten years post-fire
rasters = intersect(list.files("Data_Clean/SDMs", pattern = "\\.tif$", full=T),
                    list.files("Data_Clean/SDMs", pattern = "10yearspostfire", full=T))
name = gsub("Data_Clean/SDMs/10yearspostfire_","",rasters)
name = gsub(".tif","",name)

features = data.frame(filename = rasters, weight = 1, name = name)

features = features %>%
  filter(!name %in% c("Black.Rat", "Cat", "European.Rabbit", "Fallow.Deer", "Horse", "Pig", "Red.Fox", "Sambar")) 
features$filename = paste0(wd,"/", features$filename)
write.table(features, "Scripts/Z_files/features_10yearspostfire.txt", row.names=FALSE)

## Settings file

?write.csv
