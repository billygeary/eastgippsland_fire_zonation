
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### East Gippsland Project ####
#### Prepare Species Occurrence Data and Detection Histories ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Three data sets included in this analysis
### a) Southern Ark 2016-2017 Monitoring
### b) Long-footed Potoroo Monitoring at Tulloch Ard
### c) Spot-tailed Quoll Monitoring at W Tree

#### Setup ####
library(tidyverse)
library(camtrapR)

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1: Prepare presence-absence dataset from existing data block ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Read in presence-absence data provided by NEP Gippsland
camera_data = read.csv(here::here("Data_Raw", "NEP_SouthernArk", "all_cam_dataset_20170830.csv"))
head(camera_data)

#### Convert species column to wide format
camera_data_wide = camera_data %>% 
  filter(!Species %in% c("Antechinus", "Rats", "Ravens and Crows", "Unidentified Antechinus", "Unidentified Bandicoot",
                         "Unidentified Corvus","Unidentified Dasyurid", "Unidentified Potoroo","Unidentified Rat",
                         "Unidentified Rodent", "Unidentified Skink")) %>%
  # filter(Species %in% c("Agile Antechinus", "Black Rat", "Black Wallaby", "Bush Rat", "Cat", "Cattle", "Common Brushtail Possum",
  #          "Common Ringtail Possum", "Common Wombat","Deer", "Dog", "Dusky Antechinus","Eastern Grey Kangaroo", 
  #          "Eastern Pygmy-possum", "European Hare", "European Rabbit", "Fallow Deer", "Horse", "Koala",
  #          "Lace Monitor", "Long-footed Potoroo", "Long-nosed Bandicoot", "Long-nosed Potoroo","Mountain Brushtail Possum", 
  #          "Pig", "Red-necked Wallaby", "Red Fox","Sambar", "Short-beaked Echidna","Southern Brown Bandicoot",  
  #          "Spot-tailed Quoll", "Sugar Glider", "Swamp Rat", "Water Rat", "White-footed Dunnart")) %>%
  # 
  pivot_wider(id_cols = c(FolderID, Program, Lure, CameraModel, CameraNights, SiteID, DateDeployed, 
                          CameraID, CardID, Easting, Northing), 
              names_from = Species, values_from = present)

write.csv(camera_data_wide, here::here("Data_Processing", "sark_lfp_stq_1617_presenceabsence_allspp.csv"))

number_presences = camera_data %>% group_by(Species) %>% summarise(n = sum(present))

#### 
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2: Prepare detection histories for occupancy modelling ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# The following script is a modified version of the script developed by Lucas Bluff, NEP Gippsland

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Step 1: Read in exif data #####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# create detection histories
project_name<-"Southern Ark remote camera monitoring"
VBA_project_ID<-4879

rawimagedata<-read.table("Data_Raw/NEP_SouthernArk/SthArk_1617_camera_monitoring_DATABLOCK_20171130/EXIFoutput_SARK1617_20171124.txt",header=F,sep="\t")

# for deployment data, check dates in dd/mm/yyyy
# also person columns must be VBA user IDs

deploymentdata<-read.table("Data_Raw/NEP_SouthernArk/SthArk_1617_camera_monitoring_DATABLOCK_20171130/SARK1617_camera_data_MASTER_20171128.csv",header=T,sep=",")
deploymentdata<-deploymentdata[1:17] # remove this if the whole set of columns is desired

VBAtaxonIDlookup<-read.csv("Data_Raw/NEP_SouthernArk/SthArk_1617_camera_monitoring_DATABLOCK_20171130/camera_tag_lookup_VBA_taxonID_20170808.csv",header=T,sep=",")

VBAuploadheader<-read.csv("Data_Raw/NEP_SouthernArk/SthArk_1617_camera_monitoring_DATABLOCK_20171130/VBA_batch_header_20160706.csv",header=T,sep=",")

### setup other important precursors ###

#valid camera models (first four digits of EXIF serial number field)

valid_cams<-c('H500','H600','P850','P900')

# minimum length (in camera nights) to count as a valid deployment (otherwise filtered out)

valid_deployment_length<-5 

#enter the total number of characters in each FolderID (including spaces or underscores)
folder_name_length<-37

#enter the initial battery voltage threshold between NiMH and Alkaline
battery_threshold<-8.5

#enter the battery voltage denoting a non-functional camera
dud_voltage<-6.1

#list all non-species (i.e. method) tags to exclude
method_tags<-c('AAA','IDX','NIL','NTS','SPX','XXX','ZZZ')

#list all species tags for detailed reporting and occupancy matrices
reporting_tags<-c('STQ','CAT','FOX','HOR','COW','PIG','SAM')

#function to extract right side of text strings
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#function to extract left side of text strings
substrLeft <- function(x, n){
  substr(x, 1, nchar(x)-n)
}


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Step 2: clean and format input data #####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### clean and format the EXIF data ###

rawimagedata<-subset(rawimagedata,(substr(V6,1,4) %in% valid_cams))

rawimagedata<-rawimagedata[which(rawimagedata$V8=='1 of 5'),]

rawimagedata<-rawimagedata %>% separate_rows(V12,sep=", ")   # creates new data rows for any image with >1 tag

imagedata<-NULL

imagedata$FolderID<-substrRight(as.character(rawimagedata$V1),folder_name_length)

imagedata<-as.data.frame(imagedata)

imagedata$ImageName<-rawimagedata$V2

imagedata$DateTime<-as.POSIXct(rawimagedata$V7,format="%Y:%m:%d %H:%M:%S")

imagedata$ImageDate<-as.POSIXct(rawimagedata$V7,format="%Y:%m:%d")

#imagedata$ImageDate<-format(as.POSIXct(rawimagedata$V7,format="%Y:%m:%d %H:%M:%S"),"%Y/%m/%d")

#imagedata$ImageTime<-format(as.POSIXct(rawimagedata$V7,format="%Y:%m:%d %H:%M:%S"),"%H:%M")

#imagedata$Night<-format(as.Date(imagedata$DateTime), "%Y/%m/%d")  #assigns images taken prior to midday as the preceding date

imagedata$Night<-round(imagedata$DateTime,"days")-(24*60*60) 

imagedata$Tag<-rawimagedata$V12

imagedata$CamModel<-substr(rawimagedata$V6,1,4)

imagedata$TriggerID<-rawimagedata$V9

imagedata$Temp<-as.numeric(substrLeft(as.character(rawimagedata$V10),2))

imagedata$Voltage<-as.numeric(substrLeft(as.character(rawimagedata$V11),2))

### subset images with fauna data only ie excluding method tags ###

faunaimagedata<-subset(imagedata,!(Tag %in% method_tags))


### clean and format the deployment data including cross reference to camera specs from images ###

deploymentdata$DateDeployed<-as.POSIXct(deploymentdata$DateDeployed,format="%d/%m/%Y")
deploymentdata$DateRetrieved<-as.POSIXct(deploymentdata$DateRetrieved,format="%d/%m/%Y")
deploymentdata$DateNonFunctional<-as.POSIXct(deploymentdata$DateNonFunctional,format="%d/%m/%Y")
deploymentdata$NightsActive<-as.numeric(difftime(deploymentdata$DateNonFunctional,deploymentdata$DateDeployed,units="days"))

cam_specs<-as.data.frame(summarise(group_by(imagedata,FolderID,CamModel),
                                   TotalTriggers=n(),
                                   MeanVolt=mean(Voltage),
                                   MinVolt=min(Voltage)))

cam_specs$BattType[cam_specs$MeanVolt<battery_threshold]<-"Dummy_Lo_Volt"

cam_specs$BattType[cam_specs$MeanVolt>=battery_threshold]<-"Dummy_Hi_Volt"

startdatetemp<-as.data.frame(summarise(group_by(filter(imagedata,Tag=='AAA'),FolderID),
                                       DateDeployImage=max(ImageDate),
                                       VoltDeployImage=mean(Voltage)))

startdatetemp$DateDeployImage<-as.POSIXct(startdatetemp$DateDeployImage,format="%Y-%m-%d")

enddatetemp<-as.data.frame(summarise(group_by(filter(imagedata,Tag=='ZZZ'),FolderID),
                                     DateRetrieveImage=min(ImageDate),
                                     VoltRetrieveImage=mean(Voltage)))

enddatetemp$DateRetrieveImage<-as.POSIXct(enddatetemp$DateRetrieveImage,format="%Y-%m-%d")

faunatriggers<-as.data.frame(summarise(group_by(faunaimagedata,FolderID),FaunaTriggers=n()))

cam_specs<-merge(cam_specs,startdatetemp,by="FolderID",all=TRUE)

cam_specs<-merge(cam_specs,enddatetemp,by="FolderID",all=TRUE)

cam_specs<-merge(cam_specs,faunatriggers,by="FolderID",all=TRUE)

deploymentdata<-merge(deploymentdata,cam_specs,by="FolderID",all=TRUE)


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Step 3: ERROR CHECK #####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### if any issues detected, apply fix then rerun before moving to data export ###

# check for remaining uncertain (XXX) tags; these must be resolved before exporting EXIFoutput.txt
if(nrow(filter(imagedata,Tag=='XXX'))==0){print("No XXX tags")}else{print("WARNING RESOLVE XXX tags")}

filter(imagedata,Tag=='XXX')

# check for untagged images; these must be resolved before exporting EXIFoutput.txt
filter(faunaimagedata,is.na(Tag))
filter(imagedata,Tag=='-')

# check for images with tags not in tag-to-VBA lookup table

filter(faunaimagedata,!Tag %in% VBAtaxonIDlookup$Tag)

# check the length of FolderIDs in the metadata
print(paste0("Number of rows in deployment data: ",nrow(deploymentdata)))

# check the number of folderID/deployments matches between the EXIF data and the deployment metadata
print(paste0("Number of deployments with fauna images: ",nrow(filter(deploymentdata,FaunaTriggers!='NA'))))

filter(deploymentdata,is.na(FaunaTriggers))

# code alternatives for excluding selected FolderIDs from deploymentdata are below:
# deploymentdata<-filter(deploymentdata,!is.na(FaunaTriggers))
# deploymentdata<-filter(deploymentdata,!FolderID=="STQ113_CAM85014_CARD0926")
# deploymentdata<-filter(deploymentdata,CameraNights > valid_deployment_length)


# check the number of folders with matching FolderID between deployment and image datasets
print(paste0("Number of FolderIDs matching between deployment and image data: ",
             nrow(merge(distinct(deploymentdata,FolderID),
                        distinct(imagedata,FolderID),
                        by="FolderID"))))


# check the dates match between the EXIF data and the deployment metadata
deploymentdata$CheckDeployDate<-as.numeric(difftime(deploymentdata$DateDeployed,deploymentdata$DateDeployImage,units="days"))
filter(deploymentdata,CheckDeployDate!=0)


deploymentdata$CheckRetrieveDate<-as.numeric(difftime(deploymentdata$DateRetrieved,deploymentdata$DateRetrieveImage,units="days"))
filter(deploymentdata,CheckRetrieveDate!=0)

print(paste0("Number of deployments with date mismatch on deployment AND retrieval: ",
             nrow(filter(deploymentdata,CheckDeployDate!=0,CheckRetrieveDate!=0))))


# if incorrect camera settings found, can correct by adding/subtracting the relevant number of seconds
# as.POSIXct("12/6/1981",format="%d/%m/%Y")+24*60*60

# check minimum battery voltages
# not yet implemented; examine deployment data table voltages instead

### END ERROR CHECK ###

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Step 4: Create detection histories for occupancy modelling #####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

species_data = faunaimagedata
species_data = left_join(faunaimagedata, VBAtaxonIDlookup, by = 'Tag')

sppsites = unique(species_data$SiteID)
camopsites = unique(deploymentdata$SiteID)
intersect(sppsites, camopsites)

number_presences = camera_data %>% group_by(Species, CameraID, Lure) %>% summarise(n = max(present)) %>%
  group_by(Species) %>% summarise(n = sum(n)) %>%
  filter(n>10) %>% filter(!Species %in% c("Unidentified Rat", "Unidentified Bandicoot", "Unidentified Potoroo",
                                          "Unidentified Antechinus", "Brushtail Possums", "Rats"))

species_data = filter(species_data, COMM_NAME %in% number_presences$Species)

species_data$DateTime2 =as.character(species_data$DateTime)

# Create camera operation table
camop = cameraOperation(CTtable = deploymentdata, 
                        stationCol = "FolderID", 
                        setupCol = "DateDeployed",
                        retrievalCol = "DateRetrieved",
                        dateFormat = '%Y-%m-%d')

det.covariates = deploymentdata %>% 
  dplyr::select(FolderID, round, SiteID, CameraID, Easting, Northing,DateDeployed, NightsActive, Lure) %>%
  mutate(Lure = ifelse(Lure == "PNB-Pistachio", "Herbivore", "Predator")) %>%
  mutate(Year = lubridate::year(DateDeployed),
         JulianDate = as.numeric(format(DateDeployed, "%j")),
         SiteID_Lure = paste0(SiteID, "_", Lure))

species = unique(species_data$COMM_NAME)
species_data = species_data[-29849,] # Remove image that's not cooperating, doesnt affect Det Hist in any way. 

occ_data = list()

for (s in species){
  dethist = detectionHistory(recordTable          = species_data,
                             camOp                = camop,
                             stationCol           = "FolderID",
                             speciesCol           = "COMM_NAME",
                             recordDateTimeCol    = "DateTime",
                             species              = s,
                             occasionLength       = 1,
                             day1                 = "station",
                             timeZone             = "Australia/Victoria" ,
                             includeEffort        = FALSE)
  dethist.out = data.frame(dethist$detection_history)
  dethist.out$FolderID <- rownames(dethist.out)
  covs = det.covariates
  covs$Species = s
  dethist.out <- left_join(covs, dethist.out, by="FolderID")
  occ_data[[s]] <- dethist.out
  print(s)
}

save(occ_data, file=here::here("Data_Processing", "sark_1617_detectionhistories.RData"))

#### Step 5: Create Reporting Rate Table #### 
load(here::here("Data_Processing", "sark_1617_detectionhistories.RData"))
reprate = sapply(occ_data, function(x) {rowSums(x[,14:ncol(x)]/round(x$NightsActive), na.rm=TRUE)})
sitedata = occ_data$`Agile Antechinus`[,1:12]
reprate = cbind(sitedata, reprate)
write.csv(reprate, here::here("Data_Processing", "sark_1617_reprate.csv"))
