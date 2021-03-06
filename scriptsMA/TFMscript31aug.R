##########Load Teresa's custom functions########################
source('scriptsMA/basicFunTEG.R')

##########Load data########################
library(googledrive)
library(tidyverse)
library(ggpubr)
library(broom)
library(tidyr)


drive_download("https://docs.google.com/spreadsheets/d/1Z7HathxScqACg5vBxWm5dBTbiYMsRKGFDWu2okI2JLE/edit?usp=sharing",
               type = 'csv', path = 'dataMA/meta_sample.csv', overwrite = T)
meta <- read.csv('dataMA/meta_sample.csv', sep = ",")

drive_download("https://docs.google.com/spreadsheets/d/1oRfLFvQthr_ZxMYjOSYpnmdbclQgxbz2waraHTafJ5U/edit#gid=0",
               type = 'csv', path = 'dataMA/source_sample.csv', overwrite = T)
source <- read.csv('dataMA/source_sample.csv', sep = ",")

drive_download("https://docs.google.com/spreadsheets/d/1v2zOeuB-b_BdiEQUMSHTz9nmm9pgRseUygD7SDmDgpo/edit?usp=sharing",
               type = 'csv', path = 'dataMA/plant_sample.csv', overwrite = T)
plant <- read.csv('dataMA/plant_sample.csv', sep = ",")


#########WORLDCLIM############

#First we have to fix the data of meta by filling the gaos of map and mat with the worldclim database

library(sp)
library(raster)

WCvars <- getData("worldclim",var="bio",res=10)

# ignore glasshouse studies
meta$lat <- ifelse(meta$lat > 91 | meta$lat < -91, NA, meta$lat)
meta$log <- ifelse(meta$log > 181 | meta$log < -181, NA, meta$log)
# convert 99999 in map and mat into NA
meta[which(meta$map > 10000 | meta$mat > 10000), c('map', 'mat')] <- NA
coords <- data.frame(x=meta$log,y=meta$lat)
# filter for NA's
coords <- coords[which(!is.na(coords$x)), ]

coords<- unique(coords) # remove duplicate values

points <- SpatialPoints(coords, proj4string = WCvars@crs)

values <- raster::extract(WCvars,points)
WCvars <- cbind.data.frame(coordinates(points),values)

metaWC<- WCvars %>%
  dplyr::select(x,y,bio1,bio12)

colnames(metaWC)<- c("log","lat","matWC","mapWC")
#MAT data must be the same as meta database

metaWC$matWC<- 0.1*metaWC$matWC
# make sure you don't lose the glasshouse studies
meta<-merge(meta, metaWC, by=c('log','lat'), all.x = T)
# this line does not work now, but I don't think you need it
# metacompare<- meta %>% select(log,lat,map,mapWC,mat,matWC) # to see the diferences
source('scriptsMA/compare_with_WC.R')

detach("package:raster", unload = TRUE)
rm(WCvars, coords, points, values)
#################SWL##############################

# this is the first thing you should do
str(source)  ###check the type of variable
# these lines are needed if you use the right file

# source$d2H_permil_source<-as.numeric(as.character(source$d2H_permil_source))   #for some reason those are factors...
# source$d18O_permil_source<-as.numeric(as.character(source$d18O_permil_source))
# source<-subset(source,!d2H_permil_source=="NA")
# source<-subset(source,!d18O_permil_source=="NA") ##clean NA
# source$year<-as.factor(source$year)
source$authorYear <- paste0(source$author, '_', source$year)
source$authorYearDate <- paste0(source$author, '_', source$year, '_', source$date)
source$authorYearPlot <- paste0(source$author, '_', source$year, '_', source$plotR)
source$campaign <- paste0(source$authorYearDate, '_', source$plotR)
n_distinct(source$campaign)

source<-subset(source,label_pool=="bulk") ##only bulk soil
multiple<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-campaign) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source, na.action='na.omit',data = .x)),
    tidied = map(fit, tidy)
  ) %>%
  unnest(tidied)

multiple<-multiple[,c('campaign','term','estimate','std.error',"statistic",'p.value')]  ###select only relevant columns

###more parameters from the regression, including r squared
rsquared<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-campaign) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source,data = .x)),
    tidied = map(fit, glance)
  ) %>%
  unnest(tidied)

rsquared<-rsquared[,c('campaign','r.squared')]  ###select only relevant columns

###count the number of soil samples, in case we want to cut off using this
length<- source %>% count(campaign)

###divide into intercept and slope
intercept<-subset(multiple,term=='(Intercept)')
slope<-subset(multiple,term=='d18O_permil_source')

hist(slope$estimate)##the slopes should be positive, it will be fixed in following lines


swl<-merge(intercept,slope,by='campaign')   ###create a new table with separate columns for intercept and slope


###give proper names
colnames(swl)<-c("campaign","term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
                 "statistic.slope","p.value.slope")
##paste rsquareds
swl<-merge(swl,rsquared,by='campaign')

#past N's
swl<-merge(swl,length,by='campaign')

hist(swl$r.squared)

swl<-subset(swl,p.value.slope<0.05&n>2&estimate.slope>0)   #cutoff non-significant and poorly fit regressions
hist(swl$estimate.slope)

meta$authorYearPlot <- paste0(meta$author, '_', meta$year, '_', meta$plotR)
meta_lmwl <- meta[, c('authorYearPlot', 'slope_LMWL', 'intercept_LMWL')]
meta_lmwl <- rmDup(meta_lmwl, 'authorYearPlot')

swl[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(swl$campaign, '_', 4)
swl$authorYearPlot <- paste0(swl$author, '_', swl$year, '_', swl$plotR)
swl <- left_join(swl, meta_lmwl, by = 'authorYearPlot')
swl <- swl[, c("campaign","term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
               "statistic.slope","p.value.slope", 'slope_LMWL', 'intercept_LMWL')]

####lets screen the plots

swlplot <- left_join(source, swl, by = 'campaign')
swlplot <- subset(swlplot, p.value.slope < 0.05 & n > 2 & estimate.slope > 0)
# split in groups to see the plots more clearly
campNames <- data.frame(row.names = 1:length(unique(swlplot$campaign)))
campNames$campaign <- unique(swlplot$campaign)
campNames$crapNumber <- c(1:nrow(campNames))
swlplot <- left_join(swlplot, campNames, by = 'campaign')
swlplotL <- list()
for(i in 1:ceiling((nrow(campNames)/20))){
  swlplotL[[i]] <- swlplot[which(swlplot$crapNumber >= i*20-19 & swlplot$crapNumber <= i*20), ]
}

# not needed if you use the right file
# swlplot<-subset(swlplot,!item_source=="")  #I detected an strange problem. Data with empty item source has the same outliers (around d2H=45, d180=7).

# run this in batches so that we can see somehting
# ggplot(data=swlplot,aes(x=d18O_permil_source,y=d2H_permil_source))+
#   geom_point()+
#   geom_smooth(method=lm,se=F)+
#   facet_wrap(~campaign)+
#   stat_cor()

windows(12, 8)
# enter numbers from 1 to 9 where it says "i" to see batches of 20 plots
ggplot(data=swlplotL[[i]],aes(x=d18O_permil_source,y=d2H_permil_source))+
  geom_point()+
  geom_smooth(method=lm,se=F)+
  facet_wrap(~campaign)+
  stat_cor()

# try keeping a clean environment
rm(campNames, swlplotL, multiple, rsquared, length, slope, intercept)

######SWL-LMWL##########################
#meta$authorYear <- paste0(meta$author, '_', meta$year)

# meta$authorYearPlot <- paste0(meta$author, '_', meta$year, '_',meta$plotR)
# authorYearPlot<- source %>% select (campaign, authorYearPlot)
# authorYearPlot<- unique(authorYearPlot)
# meta <- merge(authorYearPlot, meta, by = 'authorYearPlot', all.x = T, all.y = F)
# # on this step the nrow of meta increases because there are studies with multiple sampling campaigns within a site (plot)
# 
# metasource <- left_join(swl, meta, by= 'campaign')
# # here the nrow of swl increases because there are campaigns where multiple species were measured for the same swl
# 
# #we eliminate studies without LMWL (experimental)
# # turn non-sense values of slope and intercept of lmwl into NA'
# metasource[which(metasource$slope_LMWL >= 100), c('slope_LMWL', 'intercept_LMWL')] <- NA
# 
# metasource$slopediff<- metasource$slope_LMWL-metasource$estimate.slope
# 
# slopedifference<-metasource %>%
#   group_by(campaign)%>%
#   summarise(diff=mean(slopediff,na.rm=T),count=n())
# 
# hist(slopedifference$diff)#looks nice
# # remove repetead code
# rm(authorYearPlot)
# if yo want to run a model with slopedifference, make sure repeated values are removed to avoid pseudo-replication

###########offset############################

hist(plant$d2H_permil_plant)
hist(plant$d18O_permil_plant)

plant$authorYear <- paste0(plant$author, '_', plant$year)
plant$campaign <- paste0(plant$author, '_', plant$year, '_', plant$date, '_', plant$plotR)
meta$authorYear <- paste0(meta$author, '_', meta$year)
crap <- meta %>%
  select(authorYear, pool_plant) %>%
  unique
plant <- left_join(plant, crap, by = 'authorYear')
rm(crap)
# drop values measured on leaves and exclude those that do not represent plant pool water
plant <- subset(plant, !plant_tissue == 'leaf' & pool_plant == 'yes')

offset <- left_join(plant[, c('campaign', 'd2H_permil_plant', 'd18O_permil_plant', 'species_plant', 'season', 'natural')],
                     swl, by = 'campaign')
# the nrow of offset should be same as in plant

###let's calculate the offset!
# equation (1) in Barbeta et al. 2019 HESS
offset$offset <- offset$d2H_permil_plant - offset$estimate.slope*offset$d18O_permil_plant - offset$estimate
# calcualte the d-excess with the slope of the GMWL (8)
offset$dexcess_gmwl <- offset$d2H_permil_plant - 8 * offset$d18O_permil_plant
# calculate d-excess with the slope of the corresponding LMWL
offset$dexcess <- offset$d2H_permil_plant - offset$slope_LMWL * offset$d18O_permil_plant
offset[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(offset$campaign, '_', 4)

hist(offset$offset) ###weird things, lets fix them

rareoffset<-subset(offset, offset > 25 | offset < -75)
#rareoffset<- unique(rareoffset)
hist(rareoffset$offset)
wierdCamps <- unique(rareoffset$campaign)
wierdCampsL <- list()
for(i in 1:length(wierdCamps)){
  wierdCampsL[[i]] <- subset(rareoffset, campaign == wierdCamps[i])
}
windows(12,8)
par(mfrow=c(2, 4))
for(i in 1:length(wierdCampsL)){
  plot(wierdCampsL[[i]][, 'd2H_permil_plant'] ~ wierdCampsL[[i]][, 'd18O_permil_plant'],
       main = wierdCampsL[[i]][1, 'campaign'], ylim = c(-75, -25), xlim = c(-25, 0),
       pch = 19, ylab = 'd2H (permil)', xlab = 'd18O (permil)')
  abline(a = wierdCampsL[[i]][1, 'estimate'], b = wierdCampsL[[i]][1, 'estimate.slope'])
  abline(10, 8, lty = 2)
  legend('bottomright', legend = c('plant', 'swl', 'GMWL'), pch = c(19, NA, NA), lty = c(NA, 1, 2))
}

str (rareoffset)
rareoffset<- rareoffset %>%
  select(campaign, d18O_permil_plant, d2H_permil_plant, estimate.slope, estimate, offset)

view(rareoffset)
offset <- subset(offset, author != 'Bertrand' & author != 'Eggemeyer')
rm(wierdCamps, wierdCampsL, rareoffset)

lengthplant<-offset %>% count(campaign,species_plant)

means_offset<-offset %>%
  group_by(campaign,species_plant)%>%
  summarise(mean_offset=mean(offset,na.rm=T), se_offset = s.err.na(offset),
            mean_dex_gmwl=mean(dexcess_gmwl, na.rm=T), se_dex_gmwl = s.err.na(dexcess_gmwl),
            mean_dex_lmwl=mean(dexcess, na.rm=T), se_dex_lmwl = s.err.na(dexcess),
            count_offset=n(), natural=natural[1], season=season[1])
hist(means_offset$mean_offset)


######database######
modeldata <- left_join(means_offset, swl, by = 'campaign')
modeldata <- modeldata[which(!is.na(modeldata$term.slope)), ]
modeldata[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(modeldata$campaign, '_', 4)
modeldata$authorYearPlot <- paste0(modeldata$author, '_', modeldata$year, '_', modeldata$plotR)
meta_clim_short <- meta[, c('authorYearPlot', 'log', 'lat', 'elevation', 'mapWC', 'matWC', 'climate_class')]
meta_clim_short <- rmDup(meta_clim_short, 'authorYearPlot')
modeldata <- left_join(modeldata, meta_clim_short, by = 'authorYearPlot')

modeldata$species_metaR <- modeldata$species_plant
# check that species_metaR from meta_data and species_plant from plant_data are actually the same

meta_spp_short <- meta[, c('species_metaR', 'pft', 'leaf_habit', 'leaf_shape', 'plant_group', 'growth_form')]
meta_spp_short <- rmDup(meta_spp_short, 'species_metaR')
meta_spp_short <- meta_spp_short[which(!is.na(meta_spp_short$species_metaR)),]
# create variable woody or non-woody
meta_spp_short$woodiness <- ifelse(meta_spp_short$growth_form == 'tree' | meta_spp_short$growth_form == 'shrub',
                                   meta_spp_short$growth_form, 'non-woody')
# turn non-applicable into NA's
meta_spp_short[which(meta_spp_short$leaf_habit == "not applicable"), 'leaf_habit'] <- NA
meta_spp_short[which(meta_spp_short$leaf_shape == "not applicable"), 'leaf_shape'] <- NA
meta_spp_short[which(meta_spp_short$plant_group == "not applicable"), 'plant_group'] <- NA
# get rid of class 'liana' because there is only one observation
meta_spp_short[which(meta_spp_short$growth_form == 'liana'), c('woodiness', 'pft')] <- NA

modeldata <- left_join(modeldata, meta_spp_short, by = 'species_metaR')
modeldata[, c('author', 'year', 'date', 'plot')] <- str_split_fixed(modeldata$campaign, '_', 4)
# this is your random term for the model
modeldata$authorYear <- paste0(modeldata$author, '_', modeldata$year)


# not checked any further

source$authorYearJournal<- paste0(source$authorYear, '_', source$journal)

sourcemerge<-source %>%
  select(campaign,natural,authorYearJournal)


plantmerge<-plant%>%
  select(authorYearPlot,campaign,species_plant,season)


mergeSP<- merge(sourcemerge,plantmerge, by='campaign')
mergeSP<- unique(mergeSP)

metamergePLANT<-meta %>%
  select(authorYearPlot, species_metaR,leaf_habit,leaf_shape,plant_group,growth_form)

metamergePLOT<-meta %>%
  select(authorYearPlot,climate_class,log,lat,elevation,mapWC,matWC)

mergeSPMp<- merge(mergeSP,metamergePLOT, by='authorYearPlot')

mergeSPMp$speciesMP<-paste0(mergeSPMp$authorYearPlot, '_', mergeSPMp$species_plant)

metamergePLANT$speciesMP <-paste0(metamergePLANT$authorYearPlot, '_', metamergePLANT$species_metaR)

mergeSPMpp<- merge(mergeSPMp,metamergePLANT, by='speciesMP')

mergeSPMppOF<- merge(mergeSPMpp,means_offset, by=c('campaign','species_plant'))

mergeSPMofSWL<- merge(mergeSPMppOF,swl, by='campaign')

modeldata<- merge(mergeSPMofSWL, slopedifference, by='campaign')

modeldata<- unique(modeldata) #erase duplicates


#select the useful ones

modeldata<-modeldata %>%
  select(authorYearJournal,campaign,natural ,species_plant,natural,leaf_habit,leaf_shape,plant_group,growth_form,season,
         climate_class,log,lat,elevation,mapWC,matWC,mean_offset,count.x,
         estimate.slope,std.error.slope,p.value.slope,estimate,std.error,p.value,
         r.squared,n,diff,count.y)

#give proper names

colnames(modeldata)<- c("study","campaign","natural","species_plant","leaf_habit","leaf_shape","plant_group","growth_form",
                        "season","climate_class","log","lat","elevation","map","mat","mean_offset",
                        "n_offset","SWLslope","SWLslope.std.error","SWLslope.pvalue","SWLintercept",
                        "SWLintercept.std.error","SWLintercept.pvalue","SWLrsquared","n_SWL","slopediff","n_slopediff")



#lang Index
modeldata$lang<- modeldata$map/modeldata$mat

# #woodness (not working)
# modeldata$woodness<- "non woody"
# modeldata$woodness[modeldata$growth_form == "shrub"] <- "shrub"
# modeldata$woodness[modeldata$growth_form == "tree"] <- "tree"


write.csv(modeldata,"C:\\Users\\Jabier\\Desktop\\modeldataoutliers.csv")


modeldata<- subset(modeldata, modeldata$mean_offset>-50& modeldata$mean_offset<50)
hist(modeldata$mean_offset) #Outliers from Eggemeyer and Geisler are removed
str(modeldata)

write.csv(modeldata,"C:\\Users\\Jabier\\Desktop\\modeldatas.csv")
#For the climate model we'll only use natural conditions
modeldataCLIMATE<- subset(modeldata, natural== 'natural')
write.csv(modeldataCLIMATE,"C:\\Users\\Jabier\\Desktop\\modeldataclimate.csv")

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
theme_set(theme_bw())

#cargo mapa del mundo con los paises:
world <- ne_countries(scale = "medium", returnclass = "sf")

#crear data frame con las coordenadas de mis puntos:
sites <- data.frame(modeldata$log, modeldata$lat)
names(sites)[1] <- "longitude"
names(sites)[2] <- "latitude"
sites<- unique(sites)
sites<- subset(sites, !longitude== "NA")
#para que cuadren las proyecciones: here WGS84, which is the CRS code #4326
sites<- st_as_sf(sites, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")

#plot:
gg1 <- ggplot(data = world) +
  geom_sf() +
  geom_sf(data = sites, size = 3.7, shape = 21, 
          fill = "dodgerblue2", color="dodgerblue3", alpha=0.9) +
  coord_sf(xlim = c(-170, 170), ylim = c(-90, 90), expand = FALSE)
