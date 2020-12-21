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

detach("package:raster", unload = TRUE)
rm(WCvars, coords, points, values, metaWC)

#################SWL##############################
source$authorYear <- paste0(source$author, '-', source$year)
source$authorYearDate <- paste0(source$author, '-', source$year, '-', source$date)
source$authorYearPlot <- paste0(source$author, '-', source$year, '-', source$plotR)
source$campaign <- paste0(source$authorYearDate, '-', source$plotR)


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

swl<-merge(intercept,slope,by='campaign')   ###create a new table with separate columns for intercept and slope

###give proper names
colnames(swl)<-c("campaign","term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
                 "statistic.slope","p.value.slope")
##paste rsquareds
swl<-merge(swl,rsquared,by='campaign')

#past N's
swl<-merge(swl,length,by='campaign')

swl<-subset(swl,p.value.slope<0.055&n>2&estimate.slope>0)   #cutoff non-significant and poorly fit regressions
swl[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(swl$campaign, '-', 4)
swl$authorYearPlot <- paste0(swl$author, '-', swl$year, '-', swl$plotR)

meta$authorYearPlot <- paste0(meta$author, '-', meta$year, '-', meta$plotR)
meta_lmwl <- meta[, c('authorYearPlot', 'slope_LMWL', 'intercept_LMWL')]
meta_lmwl[which(meta_lmwl$slope_LMWL > 1000), c('slope_LMWL', 'intercept_LMWL')] <- NA
meta_lmwl<- rmDup(meta_lmwl, 'authorYearPlot')

swl <- left_join(swl, meta_lmwl, by = 'authorYearPlot')

rm(multiple, rsquared, length, slope, intercept, meta_lmwl)

####lets screen the plots (put back the # after screening plots to use this script faster
#when using genratemodeldata in the model script)

swlplot <- left_join(source, swl, by = 'campaign')
swlplot <- subset(swlplot, p.value.slope < 0.055 & n > 2 & estimate.slope > 0)

#split in groups to see the plots more clearly
campNames <- data.frame(row.names = 1:length(unique(swlplot$campaign)))
campNames$campaign <- unique(swlplot$campaign)
campNames$crapNumber <- c(1:nrow(campNames))
swlplot <- left_join(swlplot, campNames, by = 'campaign')
swlplotL <- list()
for(i in 1:ceiling((nrow(campNames)/20))){
  swlplotL[[i]] <- swlplot[which(swlplot$crapNumber >= i*20-19 & swlplot$crapNumber <= i*20), ]
}

windows(12, 8)
#enter numbers from 1 to 9 where it says "i" to see batches of 20 plots
ggplot(data=swlplotL[[21]],aes(x=d18O_permil_source,y=d2H_permil_source))+
  geom_point()+
  geom_smooth(method=lm,se=F)+
  facet_wrap(~campaign)+
  stat_cor()


###########offset############################
plant$authorYear <- paste0(plant$author, '-', plant$year)
plant$authorYearPlot <- paste0(plant$author, '-', plant$year, '-', plant$plotR)
plant$campaign <- paste0(plant$author, '-', plant$year, '-', plant$date, '-', plant$plotR)
meta$authorYear <- paste0(meta$author, '-', meta$year)
crap <- meta %>%
  select(authorYear, pool_plant) %>%
  unique
plant <- left_join(plant, crap, by = 'authorYear')
rm(crap)
# drop values measured on leaves and exclude those that do not represent plant pool water
plant <- subset(plant, !plant_tissue == 'leaf' & pool_plant == 'yes')

offset <- inner_join(plant[, c('campaign', 'd2H_permil_plant', 'd18O_permil_plant', 'species_plant_complete', 'season', 'natural')],
                     swl, by = 'campaign') #inner join because we dont want plant campaigns matching issin source ones (non significative)
# the nrow of offset should be same as in plant

###let's calculate the offset!
# equation (1) in Barbeta et al. 2019 HESS
offset$offset <- offset$d2H_permil_plant - offset$estimate.slope*offset$d18O_permil_plant - offset$estimate
# calcualte the d-excess with the slope of the GMWL (8)
offset$dexcess<- offset$d2H_permil_plant - 8 * offset$d18O_permil_plant
# calculate lc-excess with the slope of the corresponding LMWL
offset$lcexcess <- offset$d2H_permil_plant - offset$slope_LMWL * offset$d18O_permil_plant -offset$intercept_LMWL


means_offset<-offset %>%
  group_by(campaign,species_plant_complete)%>%
  summarise(mean_offset=mean(offset,na.rm=T), se_offset = s.err.na(offset),
            mean_dexcess=mean(dexcess, na.rm=T), se_dexcess = s.err.na(dexcess),
            mean_lcexcess=mean(lcexcess, na.rm=T), se_lcexcess = s.err.na(lcexcess),
            count_offset=n(), natural=natural[1], season=season[1])

######database######
modeldata <- inner_join(means_offset, swl, by = 'campaign')
modeldata <- modeldata[which(!is.na(modeldata$term.slope)), ]
modeldata[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(modeldata$campaign, '-', 4)
modeldata$authorYearPlot <- paste0(modeldata$author, '-', modeldata$year, '-', modeldata$plotR)
# this is your random term
modeldata$authorYear <- paste0(modeldata$author, '-', modeldata$year)
meta_clim_short <- meta[, c('authorYearPlot', 'log', 'lat', 'elevation', 'mapWC', 'matWC', 
                            'climate_class','extraction_method_soil','extraction_method_plant',
                            'soil_measurement_method','xylem_measurement_method')]
meta_clim_short <- rmDup(meta_clim_short, 'authorYearPlot')
modeldata <- inner_join(modeldata, meta_clim_short, by = 'authorYearPlot')

modeldata$species_meta_complete <- modeldata$species_plant_complete
# check that species_metaR from meta_data and species_plant from plant_data are actually the same

meta_spp_short <- meta[, c('species_meta_complete', 'leaf_habit', 'leaf_shape', 'plant_group', 'growth_form')]
meta_spp_short <- rmDup(meta_spp_short, 'species_meta_complete')
meta_spp_short <- meta_spp_short[which(!is.na(meta_spp_short$species_meta_complete)),]
# create variable woody or non-woody
meta_spp_short$woodiness <- ifelse(meta_spp_short$growth_form == 'tree' | meta_spp_short$growth_form == 'shrub',
                                   meta_spp_short$growth_form, 'non-woody')
# turn non-applicable into NA's
meta_spp_short[which(meta_spp_short$leaf_habit == "not applicable"), 'leaf_habit'] <- NA
meta_spp_short[which(meta_spp_short$leaf_shape == "not applicable"), 'leaf_shape'] <- NA
meta_spp_short[which(meta_spp_short$plant_group == "not applicable"), 'plant_group'] <- NA
# get rid of class 'liana' because there is only one observation
meta_spp_short[which(meta_spp_short$growth_form == 'liana'), c('woodiness', 'pft')] <- NA
# here do left_join because otherwise you lose those for which the species is not defined in the plant_data file
modeldata <- left_join(modeldata, meta_spp_short, by = 'species_meta_complete')

rm(means_offset, meta, meta_clim_short, meta_spp_short,
   offset, plant, source, swl)

modeldata<-modeldata %>%
  select(authorYear,campaign, species_plant_complete,natural,leaf_habit,leaf_shape,plant_group,growth_form,woodiness,
         climate_class,log,lat,elevation,date,mapWC,matWC, extraction_method_soil,extraction_method_plant,
         soil_measurement_method,xylem_measurement_method, mean_offset,se_offset, count_offset, 
         mean_dexcess, se_dexcess, mean_lcexcess ,se_lcexcess,
         estimate.slope,std.error.slope,p.value.slope,estimate,std.error,p.value,r.squared,n)

#give proper names
colnames(modeldata) <- c("study", "campaign", "species", "natural", "leaf_habit", "leaf_shape", "plant_group",
                         "growth_form", "woodiness", "climate_class", "log", "lat", "elevation","date",
                         "map", "mat","extraction_method_soil","extraction_method_plant",
                         "soil_measurement_method","xylem_measurement_method","mean_offset", "se_offset", "n_offset", "mean_dexcess", "se_dexcess",
                         "mean_lcexcess", "se_lcexcess", "SWLslope", "SWLslope.std.error", "SWLslope.pvalue", "SWLintercept",
                         "SWLintercept.std.error", "SWLintercept.pvalue", "SWLrsquared", "n_SWL")



##########RAP, WOOD DENSITY, MYCO#######

#RAP

rap <- read.csv('dataMA/rap.csv', sep = ";")

rap<- rap[, c('species', 'RP', 'AP', 'RAP')]

means_rap<-rap %>%
  group_by(species)%>%
  summarise(RAP=mean(RAP,na.rm=T),
            RP=mean(RP, na.rm=T),
            AP=mean(AP, na.rm=T))

modeldata<- merge(modeldata,means_rap, by='species', all.x=TRUE)

#WOOD DENSITY

wd<- read.csv('dataMA/wood_density.csv', sep = ";")

wd<- wd[, c('species','wd')]

wd<-wd %>%
  group_by(species)%>%
  summarise(wd=mean(wd,na.rm=T))

modeldata <- merge(modeldata,wd, by='species', all.x=TRUE)

##myco

# myco <-read.csv('dataMA/myco.csv', sep = ",")
# 
# myco<- myco[, c('species_plant','mico')]
# colnames(myco)<-c('species','myco')
# 
# modeldata <- merge(modeldata,myco, by='species', all.x=TRUE)

##########ERA5##########
library(sf)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(tidync)
library(raster)
library(doBy)



#abres cada netcdf con brik

t2mERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                  varname="t2m")
t2mERA5_2<- brick(x= "ERA5_data/20012009ERA5wateriso-004.nc",
                  varname="t2m")
t2mERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                  varname="t2m")

#agrupas los rasterBricks en stacks
t2mERA5list<- c(t2mERA5_1,t2mERA5_2,t2mERA5_3)

t2mERA5<-stack(t2mERA5list)

modeldata[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(modeldata$campaign, '-', 4)
modeldata$studyPlot <- paste0(modeldata$study, '-', modeldata$plotR)
modeldataxy1 <- modeldata[,c('studyPlot', "log","lat")]
modeldataxy1 <- subset(modeldataxy1, !log== "NA") # por si acaso, pero me dio lo mismo
modeldataxy1 <- rmDup(modeldataxy1, 'studyPlot')
modeldataxy1$log_orig <- modeldataxy1$log
modeldataxy1$log <- modeldataxy1$log_orig + 180
modeldataxy1 <- doBy::orderBy(~studyPlot, modeldataxy1)
modeldataxy <- as.data.frame(modeldataxy1[, c('log', 'lat')])
coordinates(modeldataxy) <- ~ log + lat

t2mERA5data_1<-raster::extract(t2mERA5,modeldataxy)
t2mERA5data_1 <- cbind(modeldataxy1, as.data.frame(t2mERA5data_1))

temp <-  t2mERA5data_1 %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "temp_2m"
  )

temp$temp_C <- temp$temp_2m - 273.15
temp$date <- lubridate::ymd(str_sub(temp$dateOdd, 2, 11))
temp <- temp[, c('studyPlot', 'log', 'lat','date', 'temp_C')]
temp$month <- lubridate::month(temp$date, label = TRUE)
temp$year <- lubridate::year(temp$date)
temp$fyear <- as.factor(temp$year)


swvl1ERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                    varname="swvl1")
swvl1ERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                    varname="swvl1")
swvl1ERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                    varname="swvl1")

swvl1ERA5list<- c(swvl1ERA5_1,swvl1ERA5_2,swvl1ERA5_3)
swvl1ERA5<-stack(swvl1ERA5list)

smwl1ERA_data <-raster::extract(swvl1ERA5, modeldataxy)
smwl1ERA_data <- cbind(modeldataxy1, as.data.frame(smwl1ERA_data))

smwl1 <-  smwl1ERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "smwl1"
  )

smwl1$date <- lubridate::ymd(str_sub(smwl1$dateOdd, 2, 11))

ERA5 <- left_join(temp, smwl1[, c('studyPlot', 'date', 'smwl1')], by = c('studyPlot', 'date'))

swvl2ERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                    varname="swvl2")
swvl2ERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                    varname="swvl2")
swvl2ERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                    varname="swvl2")

swvl2ERA5list<- c(swvl2ERA5_1,swvl2ERA5_2,swvl2ERA5_3)
swvl2ERA5<-stack(swvl2ERA5list)

smwl2ERA_data <-raster::extract(swvl2ERA5, modeldataxy)
smwl2ERA_data <- cbind(modeldataxy1, as.data.frame(smwl2ERA_data))

smwl2 <-  smwl2ERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "smwl2"
  )

smwl2$date <- lubridate::ymd(str_sub(smwl2$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, smwl2[, c('studyPlot', 'date', 'smwl2')], by = c('studyPlot', 'date'))

swvl3ERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                    varname="swvl3")
swvl3ERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                    varname="swvl3")
swvl3ERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                    varname="swvl3")

swvl3ERA5list<- c(swvl3ERA5_1,swvl3ERA5_2,swvl3ERA5_3)
swvl3ERA5<-stack(swvl3ERA5list)

smwl3ERA_data <-raster::extract(swvl3ERA5, modeldataxy)
smwl3ERA_data <- cbind(modeldataxy1, as.data.frame(smwl3ERA_data))

smwl3 <-  smwl3ERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "smwl3"
  )

smwl3$date <- lubridate::ymd(str_sub(smwl3$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, smwl3[, c('studyPlot', 'date', 'smwl3')], by = c('studyPlot', 'date'))

swvl4ERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                    varname="swvl4")
swvl4ERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                    varname="swvl4")
swvl4ERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                    varname="swvl4")

swvl4ERA5list<- c(swvl4ERA5_1,swvl4ERA5_2,swvl4ERA5_3)
swvl4ERA5<-stack(swvl4ERA5list)

smwl4ERA_data <-raster::extract(swvl4ERA5, modeldataxy)
smwl4ERA_data <- cbind(modeldataxy1, as.data.frame(smwl4ERA_data))

smwl4 <-  smwl4ERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "smwl4"
  )

smwl4$date <- lubridate::ymd(str_sub(smwl4$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, smwl4[, c('studyPlot', 'date', 'smwl4')], by = c('studyPlot', 'date'))

#REVISAR EL CÄLCULO INTEGRADO

for (i in 1:4){
  ERA5[which(ERA5[, paste0('smwl', i)] < 0), paste0('smwl', i)] <- 0
}
ERA5$int_smwlA <- 0.07*1000*ERA5$smwl1 + 0.21*1000*ERA5$smwl2 + 0.72*1000*ERA5$smwl3 + 1.89*1000*ERA5$smwl4
ERA5$int_smwlB <- 0.07*1000*ERA5$smwl1 + 0.21*1000*ERA5$smwl2 + 0.72*1000*ERA5$smwl3
k <- subset(ERA5, date >= as.Date('2013-01-01') & date <= as.Date('2013-12-31'))

ggplot(data=k,aes(x=lat,y=int_smwlB))+
  geom_point()+
  facet_wrap(~month)


lai_lvERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                     varname="lai_lv")
lai_lvERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                     varname="lai_lv")
lai_lvERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                     varname="lai_lv")

lai_lvERA5list<- c(lai_lvERA5_1,lai_lvERA5_2,lai_lvERA5_3)
lai_lvERA5<-stack(lai_lvERA5list)

lai_lvERA_data <-raster::extract(lai_lvERA5, modeldataxy)
lai_lvERA_data <- cbind(modeldataxy1, as.data.frame(lai_lvERA_data))

lai_lv <-  lai_lvERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "lai_lv"
  )

lai_lv$date <- lubridate::ymd(str_sub(lai_lv$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, lai_lv[, c('studyPlot', 'date', 'lai_lv')], by = c('studyPlot', 'date'))

lai_hvERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                     varname="lai_hv")
lai_hvERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                     varname="lai_hv")
lai_hvERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                     varname="lai_hv")

lai_hvERA5list<- c(lai_hvERA5_1,lai_hvERA5_2,lai_hvERA5_3)
lai_hvERA5<-stack(lai_hvERA5list)

lai_hvERA_data <-raster::extract(lai_hvERA5, modeldataxy)
lai_hvERA_data <- cbind(modeldataxy1, as.data.frame(lai_hvERA_data))

lai_hv <-  lai_hvERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "lai_hv"
  )

lai_hv$date <- lubridate::ymd(str_sub(lai_hv$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, lai_hv[, c('studyPlot', 'date', 'lai_hv')], by = c('studyPlot', 'date'))

mperERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                   varname="mper")
mperERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                   varname="mper")
mperERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                   varname="mper")

mperERA5list<- c(mperERA5_1,mperERA5_2,mperERA5_3)
mperERA5<-stack(mperERA5list)

mperERA_data <-raster::extract(mperERA5, modeldataxy)
mperERA_data <- cbind(modeldataxy1, as.data.frame(mperERA_data))

mper <-  mperERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "mper"
  )

mper$date <- lubridate::ymd(str_sub(mper$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, mper[, c('studyPlot', 'date', 'mper')], by = c('studyPlot', 'date'))

mtprERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                   varname="mtpr")
mtprERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                   varname="mtpr")
mtprERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                   varname="mtpr")

mtprERA5list<- c(mtprERA5_1,mtprERA5_2,mtprERA5_3)
mtprERA5<-stack(mtprERA5list)

mtprERA_data <-raster::extract(mtprERA5, modeldataxy)
mtprERA_data <- cbind(modeldataxy1, as.data.frame(mtprERA_data))

mtpr <-  mtprERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "mtpr"
  )

mtpr$date <- lubridate::ymd(str_sub(mtpr$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, mtpr[, c('studyPlot', 'date', 'mtpr')], by = c('studyPlot', 'date'))

sltERA5_1<- brick(x = 'ERA5_data/19892000ERA5wateriso-002.nc',
                  varname="slt")
sltERA5_2<- brick(x = 'ERA5_data/20012009ERA5wateriso-004.nc',
                  varname="slt")
sltERA5_3<- brick(x = 'ERA5_data/20102019ERA5wateriso-003.nc',
                  varname="slt")

sltERA5list<- c(sltERA5_1,sltERA5_2,sltERA5_3)
sltERA5<-stack(sltERA5list)

sltERA_data <-raster::extract(sltERA5, modeldataxy)
sltERA_data <- cbind(modeldataxy1, as.data.frame(sltERA_data))

slt <-  sltERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "slt"
  )

slt$date <- lubridate::ymd(str_sub(slt$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, slt[, c('studyPlot', 'date', 'slt')], by = c('studyPlot', 'date'))


ERA5$month <- as.character(ERA5$month)

ERA5$month[ERA5$month=="ene"]<-"01"
ERA5$month[ERA5$month=="feb"]<-"02"
ERA5$month[ERA5$month=="mar"]<-"03"
ERA5$month[ERA5$month=="abr"]<-"04"
ERA5$month[ERA5$month=="may"]<-"05"
ERA5$month[ERA5$month=="jun"]<-"06"
ERA5$month[ERA5$month=="jul"]<-"07"
ERA5$month[ERA5$month=="ago"]<-"08"
ERA5$month[ERA5$month=="sep"]<-"09"
ERA5$month[ERA5$month=="oct"]<-"10"
ERA5$month[ERA5$month=="nov"]<-"11"
ERA5$month[ERA5$month=="dic"]<-"12"

ERA5$month <- as.factor(ERA5$month)

ERA5$truedate <- paste0(ERA5$year, '', ERA5$month)
modeldata$truedate<-modeldata$date

ERA5[, c('au', 'ye', 'pl')] <- str_split_fixed(ERA5$studyPlot, '-', 3)

ERA5$campaign <- paste0(ERA5$au, '-',ERA5$ye, '-',ERA5$truedate, '-', ERA5$pl)

modeldata <- merge(modeldata,ERA5[, c('campaign', 'smwl1', 'smwl2', 'smwl3', 'smwl4',
                                      'int_smwl', 'temp_C', 'lai_lv', 'lai_hv', 'mper', 'mtpr', 'slt')]
                   , by='campaign', all.x=TRUE)

modeldata_ultimate<- modeldata[, c("study", "campaign", "species", "natural", "leaf_habit", "leaf_shape", "plant_group",
                                   "growth_form", "woodiness", "climate_class", "log", "lat", "elevation","date",
                                   "map", "mat","extraction_method_soil","extraction_method_plant",
                                   "soil_measurement_method","xylem_measurement_method","mean_offset", "se_offset", "n_offset", "mean_dexcess", "se_dexcess",
                                   "mean_lcexcess", "se_lcexcess", "SWLslope", "SWLslope.std.error", "SWLslope.pvalue", "SWLintercept",
                                   "SWLintercept.std.error", "SWLintercept.pvalue", "SWLrsquared", "n_SWL", 'smwl1', 'smwl2', 'smwl3', 'smwl4',
                                   'int_smwl', 'temp_C', 'lai_lv', 'lai_hv', 'mper', 'mtpr', 'slt', 'RAP', "RP", "AP", "wd", "myco")]

modeldata_ultimate$AP[modeldata_ultimate$AP=="NaN"]<-NA
modeldata_ultimate$montRP[modeldata_ultimate$RP=="NaN"]<-NA
modeldata_ultimate$montRAP[modeldata_ultimate$RAP=="NaN"]<-NA

write.csv(modeldata_ultimate, "C:\\Users\\JAVI\\Desktop\\pRojects\\waterIsotopesMetaAnalysis\\dataMA\\modeldata_ultimate.csv")  
