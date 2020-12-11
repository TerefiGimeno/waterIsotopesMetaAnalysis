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
meta <- merge(meta, metaWC, by=c('log','lat'), all.x = T)

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

swl<-subset(swl,p.value.slope<0.05&n>2&estimate.slope>0)   #cutoff non-significant and poorly fit regressions
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

# swlplot <- left_join(source, swl, by = 'campaign')
# swlplot <- subset(swlplot, p.value.slope < 0.05 & n > 2 & estimate.slope > 0 & label_class == 'soil')
# 
# # split in groups to see the plots more clearly
# campNames <- data.frame(row.names = 1:length(unique(swlplot$campaign)))
# campNames$campaign <- unique(swlplot$campaign)
# campNames$crapNumber <- c(1:nrow(campNames))
# swlplot <- left_join(swlplot, campNames, by = 'campaign')
# swlplotL <- list()
# for(i in 1:ceiling((nrow(campNames)/20))){
#   swlplotL[[i]] <- swlplot[which(swlplot$crapNumber >= i*20-19 & swlplot$crapNumber <= i*20), ]
# }
# 
# windows(12, 8)
# # enter numbers from 1 to 9 where it says "i" to see batches of 20 plots
# ggplot(data=swlplotL[[18]],aes(x=d18O_permil_source,y=d2H_permil_source))+
#   geom_point()+
#   geom_smooth(method=lm,se=F)+
#   facet_wrap(~campaign)+
#   stat_cor()


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

offset <- inner_join(plant[, c('campaign', 'd2H_permil_plant', 'd18O_permil_plant', 'species_plant_complete', 'natural')],
                     swl, by = 'campaign') #inner join because we dont want plant campaigns matching issin source ones (non significative)

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
            count_offset=n(), natural=natural[1])

######database######
modeldata <- inner_join(means_offset, swl, by = 'campaign')
modeldata[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(modeldata$campaign, '-', 4)
modeldata$authorYearPlot <- paste0(modeldata$author, '-', modeldata$year, '-', modeldata$plotR)
# this is your random term
modeldata$authorYear <- paste0(modeldata$author, '-', modeldata$year)
meta_clim_short <- meta[, c('authorYearPlot', 'log', 'lat', 'elevation', 'mapWC', 'matWC')]
meta_clim_short <- rmDup(meta_clim_short, 'authorYearPlot')
# do not use inner_join here or you will lose glasshouse studies
# modeldata <- inner_join(modeldata, meta_clim_short, by = 'authorYearPlot')
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
         log,lat,elevation,date,mapWC,matWC,mean_offset, se_offset, count_offset, mean_dexcess, se_dexcess, mean_lcexcess
         ,se_lcexcess, estimate.slope,std.error.slope,p.value.slope,
         r.squared,n)

#give proper names
colnames(modeldata) <- c("study", "campaign", "species_plant", "natural", "leaf_habit", "leaf_shape", "plant_group",
                         "growth_form", "woodiness", "log", "lat", "elevation","date",
                         "map", "mat", "mean_offset", "se_offset", "n_offset", "mean_dexcess", "se_dexcess",
                         "mean_lcexcess", "se_lcexcess", "SWLslope", "SWLslope.std.error", "SWLslope.pvalue", "SWLrsquared", "n_SWL")

#write.csv(modeldata, "C:\\Users\\JAVI\\Desktop\\pRojects\\waterIsotopesMetaAnalysis\\dataMA\\new_modeldata.csv")  
