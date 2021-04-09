##########Load Teresa's custom functions########################
source('scriptsMA/basicFunTEG.R')

##########Load data########################
library(googledrive)
library(tidyverse)
library(ggpubr)
library(broom)
library(tidyr)

# double check that all variables that should be NUMERIC or ITEGRER are correct
drive_download("https://docs.google.com/spreadsheets/d/1Z7HathxScqACg5vBxWm5dBTbiYMsRKGFDWu2okI2JLE/edit?usp=sharing",
               type = 'csv', path = 'dataMA/meta_sample.csv', overwrite = T)
meta <- read.csv('dataMA/meta_sample.csv', sep = ",")
str(meta)

drive_download("https://docs.google.com/spreadsheets/d/1oRfLFvQthr_ZxMYjOSYpnmdbclQgxbz2waraHTafJ5U/edit#gid=0",
               type = 'csv', path = 'dataMA/source_sample.csv', overwrite = T)
source <- read.csv('dataMA/source_sample.csv', sep = ",")
str(source)

drive_download("https://docs.google.com/spreadsheets/d/1v2zOeuB-b_BdiEQUMSHTz9nmm9pgRseUygD7SDmDgpo/edit?usp=sharing",
               type = 'csv', path = 'dataMA/plant_sample.csv', overwrite = T)
plant <- read.csv('dataMA/plant_sample.csv', sep = ",")
str(plant)

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
source$dateInt <- source$date
crap <- source %>%
  select(campaign, dateInt) %>%
  unique

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
swl <- left_join(swl, crap, by = 'campaign')

rm(multiple, rsquared, length, slope, intercept, meta_lmwl, crap)

####lets screen the plots (put back the # after screening plots to use this script faster
#when using genratemodeldata in the model script)
# 
#  swlplot <- left_join(source, swl, by = 'campaign')
#  swlplot <- subset(swlplot, p.value.slope < 0.055 & n > 2 & estimate.slope > 0)
# 
# #split in groups to see the plots more clearly
#  campNames <- data.frame(row.names = 1:length(unique(swlplot$campaign)))
#  campNames$campaign <- unique(swlplot$campaign)
#  campNames$crapNumber <- c(1:nrow(campNames))
#  swlplot <- left_join(swlplot, campNames, by = 'campaign')
#  swlplotL <- list()
# for(i in 1:ceiling((nrow(campNames)/20))){
#    swlplotL[[i]] <- swlplot[which(swlplot$crapNumber >= i*20-19 & swlplot$crapNumber <= i*20), ]
#  }
#  
#  windows(12, 8)
#  #enter numbers from 1 to 9 where it says "i" to see batches of 20 plots
#  ggplot(data=swlplotL[[17]],aes(x=d18O_permil_source,y=d2H_permil_source))+
#    geom_point()+
#    geom_smooth(method=lm,se=F)+
#    facet_wrap(~campaign)+
#    stat_cor()


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

offset <- inner_join(plant[, c('campaign', 'd2H_permil_plant', 'd18O_permil_plant', 'species_plant_complete', 'natural','meanvalue_plant')],
                     swl, by = 'campaign') #inner join because we dont want plant campaigns matching issin source ones (non significative)
# the nrow of offset should be same as in plant

# Optional
source('scriptsMA/plot_LMWL_SWL.R')

###lets calculate the variance of 2H and 18O


###let's calculate the offset!
# equation (1) in Barbeta et al. 2019 HESS
offset$offset <- offset$d2H_permil_plant - offset$estimate.slope*offset$d18O_permil_plant - offset$estimate
# calcualte the d-excess with the slope of the GMWL (8)
offset$dexcess<- offset$d2H_permil_plant - 8 * offset$d18O_permil_plant
# calculate lc-excess with the slope of the corresponding LMWL
offset$lcexcess <- offset$d2H_permil_plant - offset$slope_LMWL * offset$d18O_permil_plant -offset$intercept_LMWL


means_offset<-offset %>%
  group_by(campaign,species_plant_complete)%>%
  summarise(mean_offset=mean(offset,na.rm=T),
            mean_dexcess=mean(dexcess, na.rm=T),
            mean_lcexcess=mean(lcexcess, na.rm=T),
            #se_lcexcess = s.err.na(lcexcess),
            mean_d2Hplant=mean(d2H_permil_plant, na.rm=T),
            se_d2Hplant = s.err.na(d2H_permil_plant),
            var_d2Hplant = var(d2H_permil_plant),
            mean_d18Oplant=mean(d18O_permil_plant, na.rm=T),
            se_d18Oplant = s.err.na(d18O_permil_plant),
           var_d18Oplant = var(d18O_permil_plant),
           covar_plant = cov(d2H_permil_plant,d18O_permil_plant),
            mean_lcexcess=mean(lcexcess, na.rm=T), 
            count_offset=n(), natural=natural[1])#, meanvalue_plant=meanvalue_plant[1])

se_means <- means_offset %>%
  group_by(campaign, species_plant_complete) %>%
  summarise(mean_se_d2H = mean(se_d2Hplant, na.rm = T),
            mean_se_d18O = mean(se_d18Oplant, na.rm = T),
            mean_var_d2H = mean(var_d2Hplant, na.rm = T),
            mean_var_d18O = mean(var_d18Oplant, na.rm = T),
            mean_cov = mean(covar_plant, na.rm = T))
            #mean_se_lc = mean(se_lcexcess, na.rm = T))

means_offset <- left_join(means_offset, se_means, by = c('campaign', 'species_plant_complete'))
means_offset[which(means_offset$count_offset < 2), c('se_d2Hplant', 'se_d18Oplant', 'covar_plant',
                                                     'var_d2Hplant', 'var_d18Oplant')] <-
  means_offset[which(means_offset$count_offset < 2), c('mean_se_d2H', 'mean_se_d18O', 'mean_cov',
                                                       'mean_var_d2H', 'mean_var_d18O')]

# meansoffsetraw<- subset(means_offset,!means_offset$meanvalue_plant=='yes') 
# 
# meansoffsetraw$sed2H <- s.err.na(meansoffsetraw$mean_d2Hplant)
# meansoffsetraw$sed18O <- s.err.na(meansoffsetraw$mean_d18Oplant)
# 
# 
# means_offset$se_d2Hplant<- ifelse(means_offset$meanvalue_plant == 'no' , 
#                                   means_offset$se_d2Hplant, meansoffsetraw$sed2H)
# 
# means_offset$se_d2Hplant<- ifelse( is.na(means_offset$se_d2Hplant), '0',
#                                   means_offset$se_d2Hplant)
# 
# means_offset$se_d18Oplant<- ifelse(means_offset$meanvalue_plant == 'no' , 
#                                    means_offset$se_d18Oplant, meansoffsetraw$sed18O)
# 
# means_offset$se_d18Oplant<- ifelse( is.na(means_offset$se_d18Oplant), '0',
#                                     means_offset$se_d18Oplant)
# rm(meansoffsetraw)
######database######
modeldata <- inner_join(means_offset, swl, by = 'campaign')
modeldata[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(modeldata$campaign, '-', 4)
modeldata$authorYearPlot <- paste0(modeldata$author, '-', modeldata$year, '-', modeldata$plotR)
# this is your random term
modeldata$authorYear <- paste0(modeldata$author, '-', modeldata$year)
meta_clim_short <- meta[, c('authorYearPlot', 'log', 'lat', 'elevation', 'mapWC', 'matWC', 
                            'slope_LMWL','intercept_LMWL',
                            'climate_class','extraction_method_soil','extraction_method_plant',
                            'soil_measurement_method','xylem_measurement_method')]
meta_clim_short <- rmDup(meta_clim_short, 'authorYearPlot')
modeldata <- left_join(modeldata, meta_clim_short, by = 'authorYearPlot')

modeldata$species_meta_complete <- modeldata$species_plant_complete
# check that species_metaR from meta_data and species_plant from plant_data are actually the same

meta_spp_short <- meta[, c('species_meta_complete', 'leaf_habit', 'leaf_shape', 'plant_group', 'growth_form')]
meta_spp_short <- rmDup(meta_spp_short, 'species_meta_complete')
meta_spp_short <- meta_spp_short[which(!is.na(meta_spp_short$species_meta_complete)),]
# create variable woody or non-woody
meta_spp_short$woodiness <- ifelse(meta_spp_short$growth_form == 'tree' | meta_spp_short$growth_form == 'shrub',
                                   meta_spp_short$growth_form, 'non-woody')
# turn non-applicable into NA's
meta_spp_short[which(meta_spp_short$growth_form == "not applicable"), c('growth_form', 'woodiness')] <- NA
meta_spp_short[which(meta_spp_short$leaf_habit == "not applicable"), 'leaf_habit'] <- NA
meta_spp_short[which(meta_spp_short$leaf_shape == "not applicable"), 'leaf_shape'] <- NA
meta_spp_short[which(meta_spp_short$plant_group == "not applicable"), 'plant_group'] <- NA
# get rid of class 'liana' because there is only one observation
meta_spp_short[which(meta_spp_short$growth_form == 'liana'), c('woodiness', 'pft')] <- NA
# here do left_join because otherwise you lose those for which the species is not defined in the plant_data file
modeldata <- left_join(modeldata, meta_spp_short, by = 'species_meta_complete')

rm(means_offset, meta_clim_short, meta_spp_short,
   offset, swl)

modeldata<-modeldata %>%
  select(authorYear,campaign, species_plant_complete,natural,leaf_habit,leaf_shape,plant_group,growth_form,woodiness,
         climate_class,log,lat,elevation,date,mapWC,matWC, extraction_method_soil,extraction_method_plant,
         soil_measurement_method,xylem_measurement_method, mean_offset, count_offset, mean_d2Hplant, se_d2Hplant,  
         mean_d18Oplant,  se_d18Oplant,var_d2Hplant, var_d18Oplant, covar_plant, slope_LMWL.x,intercept_LMWL.x, mean_lcexcess, se_lcexcess, 
         estimate.slope,std.error.slope,p.value.slope,estimate,std.error,p.value,r.squared,n,dateInt)

#give proper names
colnames(modeldata) <- c("study", "campaign", "species", "natural", "leaf_habit", "leaf_shape", "plant_group",
                         "growth_form", "woodiness", "climate_class", "log", "lat", "elevation","date",
                         "map", "mat","extraction_method_soil","extraction_method_plant",
                         "soil_measurement_method","xylem_measurement_method","mean_offset","n_offset"
                         , "mean_d2Hplant", "se_d2Hplant","mean_d18Oplant",  "se_d18Oplant",
                         'var_d2Hplant', 'var_d18Oplant', 'covar_plant', 'slope_LMWL','intercept_LMWL',
                         "mean_lcexcess", "se_lcexcess" , "SWLslope", "SWLslope.std.error", "SWLslope.pvalue", "SWLintercept",
                         "SWLintercept.std.error", "SWLintercept.pvalue", "SWLrsquared", "n_SWL", "dateInt")


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

myco <-read.csv('dataMA/myco.csv', sep = ",")

myco<- myco[, c('species_plant','mico')]
colnames(myco)<-c('species','myco')

modeldata <- merge(modeldata,myco, by='species', all.x=TRUE)

rm(means_rap, myco, wd, rap)

##########ERA5##########
library(sf)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(tidync)
library(raster)
library(doBy)



#abres cada netcdf con brik

t2mERA5_1<- brick(x = 'dataMA/era5_99.nc',
                  varname="t2m")
t2mERA5_2<- brick(x= "dataMA/era5_09.nc",
                  varname="t2m")
t2mERA5_3<- brick(x = 'dataMA/era5_19.nc',
                  varname="t2m")
t2mERA5_4<- brick(x = 'dataMA/era5_20.nc',
                  varname="t2m")

#agrupas los rasterBricks en stacks
t2mERA5list<- c(t2mERA5_1,t2mERA5_2,t2mERA5_3,t2mERA5_4)

t2mERA5<-stack(t2mERA5list)

modeldata[, c('author', 'yearPublication', 'dateChr', 'plotR')] <- str_split_fixed(modeldata$campaign, '-', 4)
modeldata$studyPlot <- paste0(modeldata$study, '-', modeldata$plotR)
modeldataxy1 <- modeldata[,c('studyPlot', "log","lat")]
modeldataxy1 <- subset(modeldataxy1, !log== "NA") # por si acaso, pero me dio lo mismo
modeldataxy1 <- rmDup(modeldataxy1, 'studyPlot')
modeldataxy1$log_orig <- modeldataxy1$log
modeldataxy1$log <- modeldataxy1$log_orig %% 360
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


swvl1ERA5_1<- brick(x = 'dataMA/era5_99.nc',
                    varname="swvl1")
swvl1ERA5_2<- brick(x = 'dataMA/era5_09.nc',
                    varname="swvl1")
swvl1ERA5_3<- brick(x = 'dataMA/era5_19.nc',
                    varname="swvl1")
swvl1ERA5_4<- brick(x = 'dataMA/era5_20.nc',
                    varname="swvl1")


swvl1ERA5list<- c(swvl1ERA5_1,swvl1ERA5_2,swvl1ERA5_3,swvl1ERA5_4)
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

swvl2ERA5_1<- brick(x = 'dataMA/era5_99.nc',
                    varname="swvl2")
swvl2ERA5_2<- brick(x = 'dataMA/era5_09.nc',
                    varname="swvl2")
swvl2ERA5_3<- brick(x = 'dataMA/era5_19.nc',
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

swvl3ERA5_1<- brick(x = 'dataMA/era5_99.nc',
                    varname="swvl3")
swvl3ERA5_2<- brick(x = 'dataMA/era5_09.nc',
                    varname="swvl3")
swvl3ERA5_3<- brick(x = 'dataMA/era5_19.nc',
                    varname="swvl3")
swvl3ERA5_4<- brick(x = 'dataMA/era5_20.nc',
                    varname="swvl3")

swvl3ERA5list<- c(swvl3ERA5_1,swvl3ERA5_2,swvl3ERA5_3,swvl3ERA5_4)
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

swvl4ERA5_1<- brick(x = 'dataMA/era5_99.nc',
                    varname="swvl4")
swvl4ERA5_2<- brick(x = 'dataMA/era5_09.nc',
                    varname="swvl4")
swvl4ERA5_3<- brick(x = 'dataMA/era5_19.nc',
                    varname="swvl4")
swvl4ERA5_4<- brick(x = 'dataMA/era5_20.nc',
                    varname="swvl4")

swvl4ERA5list<- c(swvl4ERA5_1,swvl4ERA5_2,swvl4ERA5_3,swvl4ERA5_4)
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


for (i in 1:4){
  ERA5[which(ERA5[, paste0('smwl', i)] < 0), paste0('smwl', i)] <- 0
}
ERA5$int_smwlA <- 0.07*1000*ERA5$smwl1 + 0.21*1000*ERA5$smwl2 + 0.72*1000*ERA5$smwl3 + 1.89*1000*ERA5$smwl4
ERA5$int_smwlB <- 0.07*1000*ERA5$smwl1 + 0.21*1000*ERA5$smwl2 + 0.72*1000*ERA5$smwl3

lai_lvERA5_1<- brick(x = 'dataMA/era5_99.nc',
                     varname="lai_lv")
lai_lvERA5_2<- brick(x = 'dataMA/era5_09.nc',
                     varname="lai_lv")
lai_lvERA5_3<- brick(x = 'dataMA/era5_19.nc',
                     varname="lai_lv")
lai_lvERA5_4<- brick(x = 'dataMA/era5_20.nc',
                     varname="lai_lv")

lai_lvERA5list<- c(lai_lvERA5_1,lai_lvERA5_2,lai_lvERA5_3,lai_lvERA5_4)
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

lai_hvERA5_1<- brick(x = 'dataMA/era5_99.nc',
                     varname="lai_hv")
lai_hvERA5_2<- brick(x = 'dataMA/era5_09.nc',
                     varname="lai_hv")
lai_hvERA5_3<- brick(x = 'dataMA/era5_19.nc',
                     varname="lai_hv")
lai_hvERA5_4<- brick(x = 'dataMA/era5_20.nc',
                     varname="lai_hv")

lai_hvERA5list<- c(lai_hvERA5_1,lai_hvERA5_2,lai_hvERA5_3,lai_hvERA5_4)
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

pevERA5_1<- brick(x = 'dataMA/era5_99.nc',
                   varname="pev")
pevERA5_2<- brick(x = 'dataMA/era5_09.nc',
                   varname="pev")
pevERA5_3<- brick(x = 'dataMA/era5_19.nc',
                   varname="pev")
pevERA5_4<- brick(x = 'dataMA/era5_20.nc',
                  varname="pev")

pevERA5list<- c(pevERA5_1,pevERA5_2,pevERA5_3,pevERA5_4)
pevERA5<-stack(pevERA5list)

pevERA_data <-raster::extract(pevERA5, modeldataxy)
pevERA_data <- cbind(modeldataxy1, as.data.frame(pevERA_data))

pev <-  pevERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "pev"
  )

pev$date <- lubridate::ymd(str_sub(pev$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, pev[, c('studyPlot', 'date', 'pev')], by = c('studyPlot', 'date'))

ERA5$pev<- 1000 * ERA5$pev

eERA5_1<- brick(x = 'dataMA/era5_99.nc',
                  varname="e")
eERA5_2<- brick(x = 'dataMA/era5_09.nc',
                  varname="e")
eERA5_3<- brick(x = 'dataMA/era5_19.nc',
                  varname="e")
eERA5_4<- brick(x = 'dataMA/era5_20.nc',
                varname="e")

eERA5list<- c(eERA5_1,eERA5_2,eERA5_3,eERA5_4)
eERA5<-stack(eERA5list)

eERA_data <-raster::extract(eERA5, modeldataxy)
eERA_data <- cbind(modeldataxy1, as.data.frame(eERA_data))

e <-  eERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "e"
  )

e$date <- lubridate::ymd(str_sub(e$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, e[, c('studyPlot', 'date', 'e')], by = c('studyPlot', 'date'))

ERA5$e<- 1000 * ERA5$e

tpERA5_1<- brick(x = 'dataMA/era5_99.nc',
                   varname="tp")
tpERA5_2<- brick(x = 'dataMA/era5_09.nc',
                   varname="tp")
tpERA5_3<- brick(x = 'dataMA/era5_19.nc',
                   varname="tp")
tpERA5_4<- brick(x = 'dataMA/era5_20.nc',
                 varname="tp")

tpERA5list<- c(tpERA5_1,tpERA5_2,tpERA5_3,tpERA5_4)
tpERA5<-stack(tpERA5list)

tpERA_data <-raster::extract(tpERA5, modeldataxy)
tpERA_data <- cbind(modeldataxy1, as.data.frame(tpERA_data))

tp <-  tpERA_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "tp"
  )

tp$date <- lubridate::ymd(str_sub(tp$dateOdd, 2, 11))

ERA5 <- left_join(ERA5, tp[, c('studyPlot', 'date', 'tp')], by = c('studyPlot', 'date'))

ERA5$tp<- 1000 * ERA5$tp

sltERA5_1<- brick(x = 'dataMA/era5_99.nc',
                  varname="slt")
sltERA5_2<- brick(x = 'dataMA/era5_09.nc',
                  varname="slt")
sltERA5_3<- brick(x = 'dataMA/era5_19.nc',
                  varname="slt")
sltERA5_4<- brick(x = 'dataMA/era5_20.nc',
                  varname="slt")

sltERA5list<- c(sltERA5_1,sltERA5_2,sltERA5_3,sltERA5_4)
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



ERA5year <- ERA5 %>%
  group_by(studyPlot, year) %>%
  summarise(sm1_annual = mean.na(smwl1), sm2_annual = mean.na(smwl2), sm3_annual = mean.na(smwl3),
            sm4_annual = mean.na(smwl4), smIntA_annual = mean.na(int_smwlA), smIntB_annual = mean.na(int_smwlB),
            temp_annual = mean.na(temp_C), laiLV_annual = mean.na(lai_lv), laiHV_annual = mean.na(lai_hv),
            pev_sum = sum(pev, na.rm = T),  e_sum = sum(e, na.rm = T), tp_sum = sum(tp, na.rm = T),pev_annual = mean.na(pev), 
            e_annual = mean.na(e), tp_annual = mean.na(tp))
          



ERA5_all <- left_join(ERA5, ERA5year, by = c('studyPlot', 'year'))

modeldata[which(modeldata$dateInt == 99999), 'dateInt'] <- 1900
modeldata$date <- ifelse(modeldata$dateInt <=  9999, modeldata$dateInt * 10000 + 101, modeldata$dateInt * 100 + 1)
modeldata$date <- lubridate::ymd(modeldata$date)
modeldata$year <- lubridate::year(modeldata$date)
modeldata$month <- lubridate::month(modeldata$date, label = T)
modeldata <- modeldata[ , -which(names(modeldata) %in% c("date", "log", "lat"))]
modeldata <- left_join(modeldata, ERA5_all, by = c('studyPlot', 'year', 'month'))
# for those studies where we cannot define the moth of sample collection only keep annual averages
modeldata[which(modeldata$dateInt <= 9999), c("temp_C", "smwl1", "smwl2", "smwl3", "smwl4",
                                              "int_smwlA", "int_smwlB", "lai_lv", "lai_hv", "pev","e" ,"tp")] <- NA
modeldata$slt<- round(modeldata$slt,0)

modeldata$slt[modeldata$slt=="0"]<-NA
modeldata$slt[modeldata$slt=="1"]<-"coarse"
modeldata$slt[modeldata$slt=="2"]<-"medium"
modeldata$slt[modeldata$slt=="3"]<-"medium fine"
modeldata$slt[modeldata$slt=="4"]<-"fine"
modeldata$slt[modeldata$slt=="5"]<-"very fine"
modeldata$slt[modeldata$slt=="6"]<-"organic"
modeldata$slt[modeldata$slt=="7"]<-"tropical organic"

modeldata$aridUNEP = modeldata$tp/ modeldata$pev
modeldata$aridUNEP_annual =modeldata$tp_sum/modeldata$pev_sum
modeldata$aridUNEP_dif = modeldata$aridUNEP_annual- modeldata$aridUNEP
modeldata$sm1_dif = modeldata$sm1_annual - modeldata$smwl1 
modeldata$sm2_dif = modeldata$sm2_annual - modeldata$smwl2 
modeldata$sm3_dif = modeldata$sm3_annual - modeldata$smwl3 
modeldata$sm4_dif = modeldata$sm4_annual - modeldata$smwl4 
modeldata$int_smwlA_dif = modeldata$smIntA_annual - modeldata$int_smwlA 
modeldata$int_smwlB_dif = modeldata$smIntB_annual - modeldata$int_smwlB 
modeldata$temp_dif = modeldata$temp_annual - modeldata$temp_C
modeldata$laiLV_dif = modeldata$laiLV_annual - modeldata$lai_lv
modeldata$laiHV_dif = modeldata$laiHV_annual - modeldata$lai_hv 
modeldata$e_dif = modeldata$e_annual - modeldata$e
modeldata$pev_dif = modeldata$pev_annual - modeldata$pev
modeldata$tp_dif = modeldata$tp_annual - modeldata$tp


rm (eERA5,eERA5_1,eERA5_2,eERA5_3,eERA5_4, eERA5list,
    pevERA5,pevERA5_1,pevERA5_2,pevERA5_3,pevERA5_4, pevERA5list,
    tpERA5,tpERA5_1,tpERA5_2,tpERA5_3,tpERA5_4, tpERA5list,
    lai_lvERA5,lai_lvERA5_1,lai_lvERA5_2,lai_lvERA5_3,lai_lvERA5_4, lai_lvERA5list,
    lai_hvERA5,lai_hvERA5_1,lai_hvERA5_2,lai_hvERA5_3,lai_hvERA5_4, lai_hvERA5list,
    t2mERA5,t2mERA5_1,t2mERA5_2,t2mERA5_3, t2mERA5_4, t2mERA5list,
    sltERA5,sltERA5_1,sltERA5_2,sltERA5_3, sltERA5_4, sltERA5list,
    swvl4ERA5,swvl4ERA5_1,swvl4ERA5_2,swvl4ERA5_3,swvl4ERA5_4, swvl4ERA5list,
    swvl3ERA5,swvl3ERA5_1,swvl3ERA5_2,swvl3ERA5_3,swvl3ERA5_4, swvl3ERA5list,
    swvl2ERA5,swvl2ERA5_1,swvl2ERA5_2,swvl2ERA5_3,swvl2ERA5_4, swvl2ERA5list,
    swvl1ERA5,swvl1ERA5_1,swvl1ERA5_2,swvl1ERA5_3,swvl1ERA5_4, swvl1ERA5list,
    campNames,ERA5,lai_hv,lai_hvERA_data,lai_lv,lai_lvERA_data, means_rap,
   e, eERA_data,modeldataxy,modeldataxy1,pev, pevERA_data,
    myco,rap,slt,sltERA_data,smwl1,smwl1ERA_data,smwl2, smwl2ERA_data,smwl3,
    smwl3ERA_data,smwl4, smwl4ERA_data, swlplotL, swlplot, t2mERA5data_1, temp,
    wd, mtpr, mtprERA_data, ERA5_all, ERA5year, meta, plant, source, tpERA_data, tp, meansoffsetraw)



modeldata$AP[modeldata$AP=="NaN"]<-NA
modeldata$RP[modeldata$RP=="NaN"]<-NA
modeldata$RAP[modeldata$RAP=="NaN"]<-NA

# modeldata$se_d18Oplant<- as.numeric(modeldata$se_d18Oplant)
# modeldata$se_d2Hplant<- as.numeric(modeldata$se_d2Hplant)

#######weights#####

modeldata$seSW <-sqrt(((modeldata$se_d2Hplant)^2) + 
                            ((modeldata$SWLslope * modeldata$se_d18Oplant)^2) +
                            ((modeldata$SWLslope.std.error * modeldata$mean_d18Oplant)^2)+
                            ((modeldata$SWLintercept.std.error)^2)
                            )

# revisar esta fórmula!!!
modeldata$seLC <-sqrt(((modeldata$se_d2Hplant)^2) + 
                            ((modeldata$SWLslope * modeldata$se_d18Oplant)^2)
                            )

modeldata$varSW <- modeldata$var_d2Hplant + modeldata$SWLslope^2 * modeldata$var_d18Oplant + 2*modeldata$SWLslope* modeldata$covar_plant


modeldata$varLC <- modeldata$var_d2Hplant + modeldata$slope_LMWL^2 * modeldata$var_d18Oplant + 2*modeldata$slope_LMWL* modeldata$covar_plant

modeldata$varSWL<- modeldata$SWLslope.std.error^2 * modeldata$n_SWL

