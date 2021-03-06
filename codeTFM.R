##########Load data########################
library(googledrive)
library(tidyverse)
library(ggpubr)
library(broom)
library(tidyr)

drive_download("https://docs.google.com/spreadsheets/d/1Z7HathxScqACg5vBxWm5dBTbiYMsRKGFDWu2okI2JLE/edit?usp=sharing",
               type = 'csv', path = 'dataMA/meta_sample.csv', overwrite = T)
meta <- read.csv('dataMA/meta_sample.csv', sep = ",")

drive_download("https://docs.google.com/spreadsheets/d/1udiv5w7YXCZXP2Pjnpshkzxsx5dxV1slG0OX7eztVq4/edit?usp=sharing",
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
coords <- data.frame(x=meta$log,y=meta$lat)
# filter for NA's
coords <- coords[which(!is.na(coords$x)), ]
# remove duplicate values
coords<- unique(coords)

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)
WCvars <- cbind.data.frame(coordinates(points),values)

metaWC<-WCvars %>% 
  select(x,y,bio1,bio12)

colnames(metaWC)<- c("log","lat","matWC","mapWC")
#MAT data must be the same as meta database

metaWC$matWC<- 0.1*metaWC$matWC

meta<-merge(meta,metaWC, by=c('log','lat'))

metacompare<- meta %>% select(log,lat,map,mapWC,mat,matWC)

#MAP and MAt data from the review are a bit different than the one for worldclim
#we have to decide if we fill the gaps with WC at if we use directly WC
#i'm going to use WC at this moment

#################SWL##############################

source$d2H_permil_source<-as.numeric(as.character(source$d2H_permil_source))   #for some reason those are factors...
source$d18O_permil_source<-as.numeric(as.character(source$d18O_permil_source))
source$year<-as.factor(source$year)
source$authorYear <- paste0(source$author, '_', source$year)
source$authorYearDate <- paste0(source$author, '_', source$year, '_', source$date)
source$authorYearPlot <- paste0(source$author, '_', source$year, '_', source$plotR)
source$campaign <- paste0(source$authorYearDate, '_', source$plotR)
n_distinct(source$campaign)

str(source)  ###check the type of variable

source<-subset(source,!d2H_permil_source=="NA")
source<-subset(source,!d18O_permil_source=="NA") ##clean NA

source<-subset(source,label_pool=="bulk") ##only bulk soil

multiple<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-campaign, -species_source) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source, na.action='na.omit',data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

multiple<-multiple[,c('campaign', 'species_source', 'term','estimate','std.error',"statistic",'p.value')]  ###select only relevant columns

###more parameters from the regression, including r squared
rsquared<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-campaign, -species_source) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source,data = .x)),
    tidied = map(fit, glance)
  ) %>% 
  unnest(tidied)

rsquared<-rsquared[,c('campaign', 'species_source','r.squared')]  ###select only relevant columns

###count the number of soil samples, in case we want to cut off using this
length<- source %>% count(campaign, species_source)

###divide into intercept and slope
intercept<-subset(multiple,term=='(Intercept)')
slope<-subset(multiple,term=='d18O_permil_source')

hist(slope$estimate)

weird<-subset(slope,estimate<0)  ##the slope should be positive...anyway those got cut when using the p.value, except Twining

swl<-merge(intercept,slope,by=c('campaign', 'species_source'))   ###create a new table with separate columns for intercept and slope

###give proper names
colnames(swl)<-c("campaign","species_source","term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
                 "statistic.slope","p.value.slope")
##paste rsquareds
swl<-merge(swl,rsquared,by=c('campaign', 'species_source'))
source<-merge(source,rsquared,by=c('campaign', 'species_source'))
#past N's
swl<-merge(swl,length,by=c('campaign', 'species_source'))

hist(swl$r.squared)

swl<-subset(swl,p.value.slope<0.05&n>2)   #cutoff non-significant and poorly fit regressions
swlftt<-subset(swl,p.value.slope<0.05&n>2&r.squared>0.5)
swlftt$campaignspp<- paste0(swlftt$campaign, '_', swlftt$species_source)
  
####lets plot them
souswlftt<-inner_join(swlftt,source,by=c('campaign', 'species_source'))
souswlftt$campaignspp<-paste0(souswlftt$campaign, '_', souswlftt$species_source)
souswlftt<-subset(souswlftt,!item_source=="")  #I detected a strange problem. Data with empty item source has the same outliers (around d2H=45, d180=7). 

ggplot(data=souswlftt,aes(x=d18O_permil_source,y=d2H_permil_source))+
  geom_point()+
  geom_smooth(method=lm,se=F)+
  facet_wrap(~campaignspp)+
  stat_cor()

#the problem must be related with the merging. The number of swl is the same beforme merging and after eliminating the weird ones

######SWL-LMWL##########################
meta$authorYear <- paste0(meta$author, '_', meta$year)

meta$authorYearPlot <- paste0(meta$author, '_', meta$year, '_',meta$plotR)

metasource<- inner_join(souswlftt,meta, by='authorYearPlot')

metasource<- metasource[ ,c('authorYear','campaignspp', 'term','estimate','std.error',"statistic",'p.value','estimate.slope','n','d18O_permil_source','d2H_permil_source','authorYearPlot','slope_LMWL','intercept_LMWL','climate_class','log','lat','elevation','map','mat' )]

#we eliminate studies without LMWL (experimental)
metasource<-subset(metasource,!natural=='experimental')
metasource$slopediff<- metasource$slope_LMWL-metasource$estimate.slope

slopedif<-metasource %>%
  group_by(campaignspp)%>%
  summarise(diff=mean(slopediff,na.rm=T),count=n())

hist(slopedif$diff)#looks nice
###########offset############################

plant$d2H_permil_plant<-as.numeric(as.character(plant$d2H_permil_plant))   #for some reason those are factors...
plant$d18O_permil_plant<-as.numeric(as.character(plant$d18O_permil_plant))

hist(plant$d2H_permil_plant)
hist(plant$d18O_permil_plant)

plant<-subset(plant,!plant_tissue=='leaf')  #we don't want leaves
plant$authorYear <- paste0(plant$author, '_', plant$year)
plant$campaign <- paste0(plant$authorYear, '_', plant$date, '_', plant$plotR)
plant$authorYearPlot <- paste0(plant$authorYear, '_', plant$plotR)

offst<-inner_join(plant,metasource,by=c('campaign')) 

offst<- offst[ ,c('authorYear.y','campaignsppSMP','campaign','species_plant','xylem', 'term','estimate','std.error',"statistic",'p.value','estimate.slope','n','d18O_permil_plant','d2H_permil_plant','d18O_permil_source','d2H_permil_source','season.y')]
###let's calculate the offset!

offst$offset<-offst$d2H_permil_plant-offst$estimate.slope*offst$d18O_permil_plant-offst$estimate

hist(offst$offset) ###weird things

rareoffset<-subset(offst,offset< -100)
hist(rareoffset$offset)

lengthplant<-offst %>% count(campaignspp,species_plant)

means_offset<-offst %>%
  group_by(campaignspp,species_plant)%>%
  summarise(mean_offset=mean(offset,na.rm=T),se_offset=sd(offset)/sqrt(length(offset)),count=n())
#strange that in count some has a lot of n

#######set up the model database##########

source$campaignspp<- paste0(source$campaign, '_', source$species_source)

source$authorYearJournal<- paste0(source$authorYear, '_', source$journal)

sourcemerge<-source %>%
  select(campaign, campaignspp,natural,authorYearJournal)

plantmerge<-plant%>%
  select(authorYearPlot,campaign,species_plant,season)

mergeSP<- merge(sourcemerge,plantmerge, by='campaign')

metamergePLANT<-meta %>%
  select(authorYearPlot, species_metaR,leaf_habit,leaf_shape,plant_group,growth_form)

metamergePLOT<-meta %>%
  select(authorYearPlot,climate_class,log,lat,elevation,mapWC,matWC)

mergeSPMp<- merge(mergeSP,metamergePLOT, by='authorYearPlot')

mergeSPMp$speciesMP<-paste0(mergeSPMp$authorYearPlot, '_', mergeSPMp$species_plant)

metamergePLANT$speciesMP <-paste0(metamergePLANT$authorYearPlot, '_', metamergePLANT$species_metaR) 

mergeSPMpp<- merge(mergeSPMp,metamergePLANT, by='speciesMP')

mergeSPMppOF<- merge(mergeSPMpp,means_offset, by=c('campaignspp','species_plant'))

mergeSPMofSWL<- merge(mergeSPMppOF,swlftt, by='campaignspp')

modeldata<- merge(mergeSPMofSWL, slopedif, by='campaignspp')

modeldata<- unique(modeldata)

#select the useful ones

modeldata<-modeldata %>%
  select(authorYearJournal,campaignspp,natural,species_plant,natural,leaf_habit,leaf_shape,plant_group,growth_form,season,
         climate_class,log,lat,elevation,mapWC,matWC,mean_offset,count.x,
         estimate.slope,std.error.slope,p.value.slope,estimate,std.error,p.value,
         r.squared,n,diff,count.y)
         
#give proper names

colnames(modeldata)<- c("study","campaign","natural","species_plant","leaf_habit","leaf_shape","plant_group","growth_form",
                        "season","climate_class","log","lat","elevation","map","mat","mean_offset",
                        "n_offset","SWLslope","SWLslope.std.error","SWLslope.pvalue","SWLintercept",
                        "SWLintercept.std.error","SWLintercept.pvalue","SWLrsquared","n_SWL","LMWL-SWLslopediff","n_slopediff")


write.csv(modeldata,"C:\\Users\\JAVI\\Desktop\\pRojects\\waterIsotopesMetaAnalysis\\dataMA\\modeldata.csv",row.names=T)

#######looking the data########

#That will be our dataset for the model. but first we have to fill the gaps in MAP ans MAT. 
#there are only agricultural, urban and natural studies

modeldata %>% count(climate_class) #see the distribution of classes

modeldata %>% count(plant_group)

modeldata %>% count(leaf_habit)

modeldata %>% count(leaf_shape)

modeldata %>% count(growth_form) %>% group_by(authorYearJournal)


metamergePLOT %>% select (authorYearPlot,climate_class)
metamergePLOT<- unique(metamergePLOT)
metamerge %>% count (climate_class)

meta$authorYearJournal<- paste0(meta$authorYear,'_',meta$journal) 
metaCOUNT<-meta %>% select (authorYearJournal,climate_class)
metaCOUNT<- unique(metaCOUNT)
metaCOUNT %>% count(climate_class)

climatecountTROPICAL<- subset(modeldata,climate_class=='tropical')
climatecountTROPICAL %>% count(season)
climatecountTROPICAL %>% count(species_plant) #####NO
climatecountTROPICAL %>% mean(SWLslope)

modeldataCLIMATE<-modeldata %>%
  group_by(climate_class)%>%
  summarise(mean_SWL=mean(SWLslope,na.rm=T),mean_offset=mean(mean_offset,na.rm=T),
            mean_LMWLdiff=mean(`LMWL-SWLslopediff`,na.rm=T), count=n())


######Aplicating the models########



