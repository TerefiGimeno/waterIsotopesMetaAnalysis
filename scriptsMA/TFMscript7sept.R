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
source$authorYear <- paste0(source$author, '-', source$year)
source$authorYearDate <- paste0(source$author, '-', source$year, '-', source$date)
source$authorYearPlot <- paste0(source$author, '-', source$year, '-', source$plotR)
source$campaign <- paste0(source$authorYearDate, '-', source$plotR)
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
#source0<- source %>% select(campaign,authorYearPlot)
#multiple<- left_join(multiple, source0, by='campaign')
#multiple<-unique(multiple)
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

swl[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(swl$campaign, '-', 4)
swl$authorYearPlot <- paste0(swl$author, '-', swl$year, '-', swl$plotR)

meta$authorYearPlot <- paste0(meta$author, '-', meta$year, '-', meta$plotR)
meta_lmwl <- meta[, c('authorYearPlot', 'slope_LMWL', 'intercept_LMWL')]
meta_lmwl[which(meta_lmwl$slope_LMWL > 1000), c('slope_LMWL', 'intercept_LMWL')] <- NA
meta_lmwl<- rmDup(meta_lmwl, 'authorYearPlot')
n_distinct(swl$authorYearPlot)
 

#swl[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(swl$campaign, '-', 4)
#swl$authorYearPlot <- paste0(swl$author, '-', swl$year, '-', swl$plotR)
swl <- left_join(swl, meta_lmwl, by = 'authorYearPlot')
# swl <- swl[, c("campaign","authorYearPlot","term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
#                "statistic.slope","p.value.slope",'r.squared','n','slope_LMWL','intercept_LMWL')]

#meta$authorYearPlot <- paste0(meta$author, '-', meta$year, '-', meta$plotR)
#source$authorYearPlot <- paste0(source$author, '-', source$year, '-', source$plotR)
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
#meta$authorYear <- paste0(meta$author, '-', meta$year)

# meta$authorYearPlot <- paste0(meta$author, '-', meta$year, '-',meta$plotR)
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

offset <- inner_join(plant[, c('campaign', 'd2H_permil_plant', 'd18O_permil_plant', 'species_plant', 'season', 'natural')],
                    swl, by = 'campaign') #inner join because we dont want plant campaigns matching issin source ones (non significative)
# offset0<- rmDup(offset, 'campaign')
# the nrow of offset should be same as in plant

###let's calculate the offset!
# equation (1) in Barbeta et al. 2019 HESS
offset$offset <- offset$d2H_permil_plant - offset$estimate.slope*offset$d18O_permil_plant - offset$estimate
# calcualte the d-excess with the slope of the GMWL (8)
offset$dexcess<- offset$d2H_permil_plant - 8 * offset$d18O_permil_plant
# calculate lc-excess with the slope of the corresponding LMWL
offset$lcexcess <- offset$d2H_permil_plant - offset$slope_LMWL * offset$d18O_permil_plant -offset$intercept_LMWL
#offset[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(offset$campaign, '-', 4)

hist(offset$offset) ###weird things, lets fix them


#rareoffset<- unique(rareoffset)
#hist(rareoffset$offset)
#wierdCamps <- unique(rareoffset$campaign)
#wierdCampsL <- list()
#for(i in 1:length(wierdCamps)){
  #wierdCampsL[[i]] <- subset(rareoffset, campaign == wierdCamps[i])
#}
#windows(12,8)
#par(mfrow=c(2, 4))
#for(i in 1:length(wierdCampsL)){
  #plot(wierdCampsL[[i]][, 'd2H_permil_plant'] ~ wierdCampsL[[i]][, 'd18O_permil_plant'],
       #main = wierdCampsL[[i]][1, 'campaign'], ylim = c(-75, -25), xlim = c(-25, 0),
       #pch = 19, ylab = 'd2H (permil)', xlab = 'd18O (permil)')
  #abline(a = wierdCampsL[[i]][1, 'estimate'], b = wierdCampsL[[i]][1, 'estimate.slope'])
  #abline(10, 8, lty = 2)
  #legend('bottomright', legend = c('plant', 'swl', 'GMWL'), pch = c(19, NA, NA), lty = c(NA, 1, 2))
#}

#str (rareoffset)
#rareoffset<- rareoffset %>%
  #select(campaign, d18O_permil_plant, d2H_permil_plant, estimate.slope, estimate, offset)

#view(rareoffset)
offset <- subset(offset, author != 'Bertrand' & author != 'Eggemeyer')
#rm(wierdCamps, wierdCampsL, rareoffset)

lengthplant<-offset %>% count(campaign,species_plant)

means_offset<-offset %>%
  group_by(campaign,species_plant)%>%
  summarise(mean_offset=mean(offset,na.rm=T), se_offset = s.err.na(offset),
            mean_dexcess=mean(dexcess, na.rm=T), se_dexcess = s.err.na(dexcess),
            mean_lcexcess=mean(lcexcess, na.rm=T), se_lcexcess = s.err.na(lcexcess),
            count_offset=n(), natural=natural[1], season=season[1])

hist(means_offset$mean_offset)
hist(means_offset$mean_lcexcess)
hist(means_offset$mean_dexcess)

#lets fix the outliers for the offset

# rareoffset<-subset(means_offset, mean_offset > 25 | mean_offset < -50)
# means_offset <- subset(means_offset,mean_offset < 50) 
#data over 50 belonged to eggmeyer and neef further review

######database######
modeldata <- inner_join(means_offset, swl, by = 'campaign')
modeldata <- modeldata[which(!is.na(modeldata$term.slope)), ]
modeldata[, c('author', 'year', 'date', 'plotR')] <- str_split_fixed(modeldata$campaign, '-', 4)
modeldata$authorYearPlot <- paste0(modeldata$author, '-', modeldata$year, '-', modeldata$plotR)
meta_clim_short <- meta[, c('authorYearPlot', 'log', 'lat', 'elevation', 'mapWC', 'matWC', 'climate_class')]
meta_clim_short <- rmDup(meta_clim_short, 'authorYearPlot')
modeldata <- inner_join(modeldata, meta_clim_short, by = 'authorYearPlot')

modeldata$species_metaR <- modeldata$species_plant
# check that species_metaR from meta_data and species_plant from plant_data are actually the same

meta_spp_short <- meta[, c('author','year','species_metaR', 'pft', 'leaf_habit', 'leaf_shape', 'plant_group', 'growth_form')]
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
meta_spp_short$authorYear<- paste0(meta_spp_short$author,'-',meta_spp_short$year) # this is your random term for the model
modeldata <- inner_join(modeldata, meta_spp_short, by = 'species_metaR')
#modeldata[, c('author', 'year', 'date', 'plot')] <- str_split_fixed(modeldata$campaign, '-', 4)

#modeldata$authorYear <- paste0(modeldata$author, '-', modeldata$year)

rm(lengthplant,means_offset,meta,meta_clim_short,meta_lmwl,meta_spp_short, metaWC,
   offset, plant, source, swl, swlplot)

modeldata<-modeldata %>%
  select(authorYear,campaign, species_plant,natural,leaf_habit,leaf_shape,plant_group,growth_form,woodiness,season,
         climate_class,log,lat,elevation,mapWC,matWC,mean_offset, se_offset, count_offset, mean_dexcess, se_dexcess, mean_lcexcess
         ,se_lcexcess, estimate.slope,std.error.slope,p.value.slope,estimate,std.error,p.value,
         r.squared,n)

#give proper names

colnames(modeldata)<- c("study","campaign","species_plant","natural","leaf_habit","leaf_shape","plant_group","woodiness","growth_form",
                        "season","climate_class","log","lat","elevation","map","mat","mean_offset","se_offset",
                        "n_offset","mean_dexcess","se_dexcess","mean_lcexcess","se_lcexcess","SWLslope","SWLslope.std.error","SWLslope.pvalue","SWLintercept",
                        "SWLintercept.std.error","SWLintercept.pvalue","SWLrsquared","n_SWL") 

#lang Index
modeldata$lang<- modeldata$mapWC/modeldata$matWC

#absolute offset

modeldata$absoluteOffset<-abs(modeldata$mean_offset)

hist(modeldata$absoluteOffset)
mean(modeldata$absoluteOffset)

modeldata$woodiness <- ifelse(modeldata$growth_form == 'tree' | modeldata$growth_form == 'shrub',
                                   modeldata$growth_form, 'non-woody')

write.csv(modeldata, "C:\\Users\\JAVI\\Desktop\\pRojects\\waterIsotopesMetaAnalysis\\dataMA\\modeldata.csv")  

rm(swl0, rareoffset)
#####Lets see some results before the models

####map##########

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
theme_set(theme_bw())

#world map load
world <- ne_countries(scale = "medium", returnclass = "sf")

#dataframe of coords and climate:
sites <- data.frame(modeldata$log, modeldata$lat)
names(sites)[1] <- "longitude"
names(sites)[2] <- "latitude"
sites<- unique(sites)
sites<- subset(sites, !longitude== "NA")

#same projections (for WGS84 CRS code is #4326)
sites<- st_as_sf(sites, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")

#plot:
map<-ggplot(data = world) +
  geom_sf() +
  geom_sf(data = sites, size = 3.7, shape = 21, 
          fill = "dodgerblue2", color="dodgerblue3", alpha=0.9) +
  coord_sf(xlim = c(-170, 170), ylim = c(-90, 90), expand = FALSE)

map

ggsave("map.jpg", map, width = 20, height = 10, dpi = 900)

rm(map, sites, world)
#######database exploration#########

n_distinct(modeldata$campaign)

n_distinct(modeldata$study)

n_distinct(modeldata$species_plant)

ggplot(data=modeldata, aes(x=mean_offset, y=mean_lcexcess)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  
 
#######models (SWL)#######

myresplot <- function(x,df) { #Home-made residual plot
  e <- resid(x, type="pearson")
  f <- fitted(x)
  obs <- x@frame[,1]
  par(mfrow=c(1,2), pty="m",mar=c(3,3,1,1),mgp = c(2, 1, 0), ps=12, bty="l")
  plot(f,e, pch=1, cex=0.8, col=alpha("black",0.3),
       xlab="Fitted values", ylab="Pearson residuals")
  abline(h = 0, lty = 2, col="black", lw=1)
  plot(f ~ obs, xlab="Observed values", ylab="Fitted values", type="n")
  points(obs,f, pch=1, cex=0.8, col=alpha("black",0.3))
  abline(a = 0, b = 1, lty = 2, col="Black", lw=1)
}

library(lme4)
library(tidyverse)
library(GGally)
library(MuMIn)
library(broom.mixed)
library(visreg)

#Lets begin with SWL slope 

natdata<- subset(modeldata, natural== 'natural') #reject hydro-managed plots

#standarize
natdata0<-natdata
natdata0$SWLslope <- (natdata0$SWLslope-mean(natdata0$SWLslope))/sd(natdata0$SWLslope)

natdata0$lang <- (natdata0$lang - mean(natdata0$lang))/sd(natdata0$lang)
natdata0$map <- (natdata0$map -mean(natdata0$map))/sd(natdata0$map)
natdata0$mat <- (natdata0$mat -mean(natdata0$mat))/sd(natdata0$mat)

#we dont want not applicable values
natdata[which(natdata$season == "not applicable"), 'season'] <- NA
natdata0[which(natdata0$season == "not applicable"), 'season'] <- NA
#general model

#season may be related with map and mat but their behaivior can be different in both seasons

#mmyvar<-c("climate_class", "map", "mat", "season")
#pairs(natdata0[,mmyvar],upper.panel=panel.smooth,diag.panel=panel.hist,lower.panel= panel.cor ,cex.cor=4)

swls1<- lmer(SWLslope ~ mat*season + map*season + 
                  (1|study) -1, data=natdata0, na.action = "na.omit", REML = FALSE)
swls2<- lmer(SWLslope ~ lang*season + 
                    (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls3<- lmer(SWLslope ~ climate_class*season +
                    (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls4<- lmer(SWLslope ~ season +
             (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls5<- lmer(SWLslope ~ mat +
             (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls6<- lmer(SWLslope ~ map +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls7<- lmer(SWLslope ~ lang +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls8<- lmer(SWLslope ~ climate_class +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)


AIC(swls1,swls2,swls3,swls4,swls5,swls6,swls7,swls8) 
AICc(swls1,swls2,swls3,swls4,swls5,swls6,swls7,swls8)

#swls3 is the best model

summary(swls3)
anova(swls3) 


#dredge won't work as season has NA's
#result1 <- dredge(swls3, rank = AIC, beta = "sd")

swls3.1<- lmer(SWLslope ~ climate_class*season +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls3.2<- lmer(SWLslope ~ climate_class +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls3.3<- lmer(SWLslope ~ season +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls3.4<- lmer(SWLslope ~ climate_class:season + climate_class +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
swls3.5<- lmer(SWLslope ~ climate_class:season + season +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)

AIC(swls3.1,swls3.2,swls3.3,swls3.4,swls3.5) 
AICc(swls3.1,swls3.2,swls3.3,swls3.4,swls3.5)

#swls3.1,4 and 5 have the same AIC and df. I'll choose swls3.5

swlbest <-lmer(SWLslope ~ climate_class:season + season +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = TRUE) #has to be true

summary(swlbest)
myresplot(swlbest, natdata) #looks good
hist(natdata0$SWLslope) 

visreg(swlbest, xvar = "season",
       scale = "response", gg=TRUE) 
visreg(swlbest, xvar = "climate_class",
       scale = "response", gg=TRUE) #both makes sense


#Unstandarized plots:

ggplot(data=natdata, aes(x=climate_class, y=SWLslope)) +
  geom_point(alpha=0.2, aes()) +
  geom_boxplot()

ggplot(data=natdata, aes(x=season, y=SWLslope)) +
  geom_point(alpha=0.2, aes()) +
  geom_boxplot()

ggplot(data=natdata, aes(x=climate_class, y=SWLslope)) +
  geom_point(alpha=0.2, aes()) +
  geom_boxplot()+
  facet_grid(.~season)  #i dont know how to unplot the NA's

rm(swls1,swls2,swls3,swls4,swls5,swls6,swls7,swls8,swls3.1,swls3.2,swls3.3,swls3.4,swls3.5)
##We are going to see what happends inside each climate class for the SWL

natdata_arid <- natdata %>% 
  filter(climate_class=="arid")

natdata_cold <- natdata %>% 
  filter(climate_class=="cold")

natdata_trop <- natdata %>% 
  filter(climate_class=="tropical")

natdata_warm <- natdata %>% 
  filter(climate_class=="warm")

#standarized ones

natdata_arid0 <- natdata0 %>% 
  filter(climate_class=="arid")

natdata_cold0 <- natdata0 %>% 
  filter(climate_class=="cold")

natdata_trop0 <- natdata0 %>% 
  filter(climate_class=="tropical")

natdata_warm0 <- natdata0 %>% 
  filter(climate_class=="warm")


######SWL_warm########

swlw1<- lmer(SWLslope ~ mat*season + map*season + 
                  (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw2<- lmer(SWLslope ~ lang*season + 
                   (1|study) -1, data=natdata_warm0, na.action = "na.omit",REML = FALSE)
swlw3<- lmer(SWLslope ~ mat + 
               (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw4<- lmer(SWLslope ~ map + 
               (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw5<- lmer(SWLslope ~ map + mat +
               (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw6<- lmer(SWLslope ~ lang + 
               (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)

AIC(swlw1, swlw2, swlw3, swlw4,swlw5, swlw6) 
# swlw1 is better


swlw1.1<- lmer(SWLslope ~ mat + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.2<- lmer(SWLslope ~ map + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.3<- lmer(SWLslope ~ season + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.4<- lmer(SWLslope ~ mat:season + map*season + mat +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.5<- lmer(SWLslope ~ mat:season + map*season + season +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.6<- lmer(SWLslope ~ mat*season + map:season + season +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.7<- lmer(SWLslope ~ mat*season + map:season + map +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.8<- lmer(SWLslope ~ mat:season + map:season + season +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.9<- lmer(SWLslope ~ mat:season + map:season + map +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.10<- lmer(SWLslope ~ mat:season + map:season + mat +
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.11<- lmer(SWLslope ~ mat*season + map*season + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.12<- lmer(SWLslope ~ mat:season + mat + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.13<- lmer(SWLslope ~ mat:season + season + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.14<- lmer(SWLslope ~ season + map:season + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.15<- lmer(SWLslope ~ mat + map:season + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.16<- lmer(SWLslope ~ mat:season + 
                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
swlw1.17<- lmer(SWLslope ~  map:season + 
                  (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)

AIC(swlw1.1,swlw1.2,swlw1.3,swlw1.4,swlw1.5,swlw1.6,swlw1.7,swlw1.8,swlw1.9,swlw1.10,
    swlw1.11,swlw1.12,swlw1.13,swlw1.14,swlw1.15,swlw1.16,swlw1.17)
#1.84 to 1.8 are the same values, i choose swlw1.8


swlbest_warm <- lmer(SWLslope ~ mat:season + map:season + season +
                                 (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = TRUE)

rm(swlw1.1,swlw1.2,swlw1.3,swlw1.4,swlw1.5,swlw1.6,swlw1.7,swlw1.8,swlw1.9,swlw1.10,
    swlw1.11,swlw1.12,swlw1.13,swlw1.14,swlw1.15,swlw1.16,swlw1.17,swlw1, swlw2, swlw3, swlw4,swlw5, swlw6)

summary(swlbest_warm)
anova(swlbest_warm)
myresplot(swlbest_warm, natdata_warm0)
hist(natdata_warm$SWLslope)


visreg(swlbest_warm, xvar = "map", by = "season", 
       scale = "response", gg=TRUE)
visreg(swlbest_warm, xvar = "mat", by = "season", 
       scale = "response", gg=TRUE)

#interesting. SWL slope increases with map (gets closer to lmwl) in dry and wet plots. 
#but its more drastic in dry season plots. For mat, it increases in dry seasons but decreases on wet ones
#(I thought that it'll also decrease un dry seasons)

#unestandarized plots (rawdata)
ggplot(data=natdata_warm, aes(x=map, y=SWLslope, colour=season)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F) +
  facet_grid(season ~ .)

ggplot(data=natdata_warm, aes(x=mat, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .)
#not the same results

######SWL_arid#######
swla1<- lmer(SWLslope ~ mat*season + map*season + 
               (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla2<- lmer(SWLslope ~ lang*season + 
               (1|study) -1, data=natdata_arid0, na.action = "na.omit",REML = FALSE)
swla3<- lmer(SWLslope ~ mat + 
               (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla4<- lmer(SWLslope ~ map + 
               (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla5<- lmer(SWLslope ~ map + mat +
               (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla6<- lmer(SWLslope ~ lang + 
               (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)

AIC(swla1, swla2, swla3, swla4,swla5, swla6) 

# swla1 is better


swla1.1<- lmer(SWLslope ~ mat + 
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.2<- lmer(SWLslope ~ map + 
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.3<- lmer(SWLslope ~ season + 
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.4<- lmer(SWLslope ~ mat:season + map*season + mat +
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.5<- lmer(SWLslope ~ mat:season + map*season + season +
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.6<- lmer(SWLslope ~ mat*season + map:season + season +
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.7<- lmer(SWLslope ~ mat*season + map:season + map +
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.8<- lmer(SWLslope ~ mat:season + map:season + season +
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.9<- lmer(SWLslope ~ mat:season + map:season + map +
                 (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.10<- lmer(SWLslope ~ mat:season + map:season + mat +
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.11<- lmer(SWLslope ~ mat*season + map*season + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.12<- lmer(SWLslope ~ mat:season + mat + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.13<- lmer(SWLslope ~ mat:season + season + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.14<- lmer(SWLslope ~ season + map:season + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.15<- lmer(SWLslope ~ mat + map:season + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.16<- lmer(SWLslope ~ mat:season + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
swla1.17<- lmer(SWLslope ~  map:season + 
                  (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)

AIC(swla1.1,swla1.2,swla1.3,swla1.4,swla1.5,swla1.6,swla1.7,swla1.8,swla1.9,swla1.10,
    swla1.11,swla1.12,swla1.13,swla1.14,swla1.15,swla1.16,swla1.17)
#looks like the best is swla1.17

swlbest_arid<- lmer(SWLslope ~  map:season + 
                      (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = TRUE)

rm(swla1.1,swla1.2,swla1.3,swla1.4,swla1.5,swla1.6,swla1.7,swla1.8,swla1.9,swla1.10,
   swla1.11,swla1.12,swla1.13,swla1.14,swla1.15,swla1.16,swla1.17,swla1, swla2, swla3, swla4,swla5, swla6)

summary(swlbest_arid)
anova(swlbest_arid)
myresplot(swlbest_arid, natdata_arid0)
hist(natdata_warm$SWLslope)


visreg(swlbest_warm, xvar = "season", 
       scale = "response", gg=TRUE) #slopes are stepper in wet seasons (makes sense)

#unestandarized plots (rawdata)
ggplot(data=natdata_arid, aes(x=season, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_boxplot() #different results as before

######SWL_tropical#######

swlt1<- lmer(SWLslope ~ mat*season + map*season + 
               (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt2<- lmer(SWLslope ~ lang*season + 
               (1|study) -1, data=natdata_trop0, na.action = "na.omit",REML = FALSE)
swlt3<- lmer(SWLslope ~ mat + 
               (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt4<- lmer(SWLslope ~ map + 
               (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt5<- lmer(SWLslope ~ map + mat +
               (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt6<- lmer(SWLslope ~ lang + 
               (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)

AIC(swlt1, swlt2, swlt3, swlt4,swlt5, swlt6) 

# swla1 is better


swlt1.1<- lmer(SWLslope ~ mat + 
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.2<- lmer(SWLslope ~ map + 
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.3<- lmer(SWLslope ~ season + 
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.4<- lmer(SWLslope ~ mat:season + map*season + mat +
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.5<- lmer(SWLslope ~ mat:season + map*season + season +
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.6<- lmer(SWLslope ~ mat*season + map:season + season +
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.7<- lmer(SWLslope ~ mat*season + map:season + map +
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.8<- lmer(SWLslope ~ mat:season + map:season + season +
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.9<- lmer(SWLslope ~ mat:season + map:season + map +
                 (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.10<- lmer(SWLslope ~ mat:season + map:season + mat +
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.11<- lmer(SWLslope ~ mat*season + map*season + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.12<- lmer(SWLslope ~ mat:season + mat + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.13<- lmer(SWLslope ~ mat:season + season + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.14<- lmer(SWLslope ~ season + map:season + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.15<- lmer(SWLslope ~ mat + map:season + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.16<- lmer(SWLslope ~ mat:season + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
swlt1.17<- lmer(SWLslope ~  map:season + 
                  (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)

AIC(swlt1.1,swlt1.2,swlt1.3,swlt1.4,swlt1.5,swlt1.6,swlt1.7,swlt1.8,swlt1.9,swlt1.10,
    swlt1.11,swlt1.12,swlt1.13,swlt1.14,swlt1.15,swlt1.16,swlt1.17)
#swlt1.13 is the best model

swlbest_tropical<- lmer(SWLslope ~ mat:season + season + 
                          (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = TRUE)

rm(swlt1.1,swlt1.2,swlt1.3,swlt1.4,swlt1.5,swlt1.6,swlt1.7,swlt1.8,swlt1.9,swlt1.10,
    swlt1.11,swlt1.12,swlt1.13,swlt1.14,swlt1.15,swlt1.16,swlt1.17,swlt1, swlt2, swlt3, swlt4,swlt5, swlt6)

summary(swlbest_tropical)
anova(swlbest_tropical)
myresplot(swlbest_tropical, natdata_trop0) #strange...
hist(natdata_trop0$SWLslope) #...


visreg(swlbest_warm, xvar = "mat", by= "season",
       scale = "response", gg=TRUE) #similar as found in SWLwarm model

#unestandarized plots (rawdata)
ggplot(data=natdata_trop, aes(x=mat, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .) #opposite as the visreg

######SWL_cold#######

swlc1<- lmer(SWLslope ~ mat*season + map*season + 
               (1|study) -1, data=natdata_cold0, na.action = "na.omit", REML = FALSE)
swlc2<- lmer(SWLslope ~ lang*season + 
               (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
swlc3<- lmer(SWLslope ~ mat + 
               (1|study) -1, data=natdata_cold0, na.action = "na.omit", REML = FALSE)
swlc4<- lmer(SWLslope ~ map + 
               (1|study) -1, data=natdata_cold0, na.action = "na.omit", REML = FALSE)
swlc5<- lmer(SWLslope ~ map + mat +
               (1|study) -1, data=natdata_cold0, na.action = "na.omit", REML = FALSE)
swlc6<- lmer(SWLslope ~ lang + 
               (1|study) -1, data=natdata_cold0, na.action = "na.omit", REML = FALSE)

AIC(swlc1, swlc2, swlc3, swlc4,swlc5, swlc6) 

#swlc2 is the best

swlc2.1<- lmer(SWLslope ~ lang*season + 
                 (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
swlc2.2<- lmer(SWLslope ~ lang:season + season +
               (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
swlc2.3<- lmer(SWLslope ~ lang:season + lang +
                 (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
swlc2.4<- lmer(SWLslope ~ lang + season + 
                 (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
swlc2.5<- lmer(SWLslope ~ lang + 
                 (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
swlc2.6<- lmer(SWLslope ~ season + 
                 (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)

AIC(swlc2.1, swlc2.2, swlc2.3, swlc2.4,swlc2.5, swlc2.6) 
#swlc2.3 is the best
swlbest_cold<- lmer(SWLslope ~ lang:season + lang +
                      (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = TRUE)

rm(swlc1, swlc2, swlc3, swlc4,swlc5, swlc6,swlc2.1, swlc2.2, swlc2.3, swlc2.4,swlc2.5, swlc2.6)

summary(swlbest_cold)
anova(swlbest_cold)
myresplot(swlbest_cold, natdata_cold0) #strange...
hist(natdata_cold0$SWLslope) #...


visreg(swlbest_cold, xvar = "lang", by= "season",
       scale = "response", gg=TRUE) #makes sense with the warm and tropical findings

#unestandarized plots (rawdata)
ggplot(data=natdata_cold, aes(x=lang, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .) #similar as the visreg

#####offset model (pft)######
#standarize continuous, fot pft we'll use the whole database
modeldata0<- modeldata
modeldata0$mean_offset <- (modeldata0$SWLslope-mean(modeldata0$SWLslope))/sd(modeldata0$SWLslope)

#pant group
offGroup<-lmer(mean_offset ~plant_group  +
       (1|study) -1, data=modeldata0, na.action = "na.omit",REML = TRUE)

summary(offGroup)
myresplot(offGroup, modeldata0) 
hist(modeldata0$mean_offset) 


visreg(offGroup, xvar = "plant_group", 
       scale = "response", gg=TRUE) #low effect

#unestandarized
ggplot(data=modeldata, aes(x=plant_group, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() #similar results 


#woodiness
offwood<-lmer(mean_offset ~woodiness  +
                 (1|study) -1, data=modeldata0, na.action = "na.omit",REML = TRUE)

summary(offwood)
myresplot(offwood, modeldata0) 
hist(modeldata0$mean_offset) 


visreg(offwood, xvar = "woodiness", 
       scale = "response", gg=TRUE) #low effect

#unestandarized
ggplot(data=modeldata, aes(x=woodiness, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() #similar results 

#leaf habit
offhabit<-lmer(mean_offset ~leaf_habit  +
                (1|study) -1, data=modeldata0, na.action = "na.omit",REML = TRUE)

summary(offhabit)
myresplot(offhabit, modeldata0) 
hist(modeldata0$mean_offset) 


visreg(offhabit, xvar = "leaf_habit", 
       scale = "response", gg=TRUE) 
#unestandarized
ggplot(data=modeldata, aes(x=leaf_habit, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() 

#leaf shape
offshape<-lmer(mean_offset ~leaf_shape  +
                 (1|study) -1, data=modeldata0, na.action = "na.omit",REML = TRUE)

summary(offshape)
myresplot(offhabit, modeldata0) 
hist(modeldata0$mean_offset) 


visreg(offshape, xvar = "leaf_shape", 
       scale = "response", gg=TRUE) 
#unestandarized
ggplot(data=modeldata, aes(x=leaf_shape , y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() 

AIC(offGroup,offhabit,offshape,offwood) #offgroup is the best

#####offset model (climate)######
#now we use natural database
offc1<- lmer(mean_offset ~ mat*season + map*season + 
               (1|study) -1, data=natdata0, na.action = "na.omit", REML = FALSE)
offc2<- lmer(mean_offset ~ lang*season + 
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc3<- lmer(mean_offset ~ climate_class*season +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc4<- lmer(mean_offset ~ season +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc5<- lmer(mean_offset ~ mat +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc6<- lmer(mean_offset ~ map +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc7<- lmer(mean_offset ~ lang +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc8<- lmer(mean_offset ~ climate_class +
               (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)

AIC(offc1,offc2,offc3,offc4,offc5,offc6,offc7, offc8)
#offc2 is the best model

offc2.1<- lmer(mean_offset ~ lang*season + 
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc2.2<- lmer(mean_offset ~ lang:season + season +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc2.3<- lmer(mean_offset ~ lang:season + lang +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc2.4<- lmer(mean_offset ~ lang + season + 
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc2.5<- lmer(mean_offset ~ lang + 
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offc2.6<- lmer(mean_offset ~ season + 
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)


AIC(offc2.1,offc2.2,offc2.3,offc2.4,offc2.5,offc2.6)

#offc2.3 is the best model

offbest_climate<- lmer(mean_offset ~ lang:season + lang +
                 (1|study) -1, data=natdata0, na.action = "na.omit",REML = TRUE)

rm(offc1,offc2,offc3,offc4,offc5,offc6,offc7,offc8,offc2.1,offc2.2,offc2.3,offc2.4,offc2.5,offc2.6)

summary(offbest_climate)
anova(offbest_climate)
myresplot(offbest_climate, natdata0) 
hist(natdata0$mean_offset) 


visreg(offbest_climate, xvar = "lang", by= "season",
       scale = "response", gg=TRUE) #looks similar for dry and wet seasons. less humidity, more offset (makes sense)

#unestandarized plot (rawdata)
ggplot(data=natdata, aes(x=lang, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .) #similar as the visreg

##### offset's general model########

offgeneral1<- lmer(mean_offset ~ lang:season + lang + plant_group +
                    (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
#this should be the best model, but im going to test more similar


offgeneral2<- lmer(mean_offset ~ lang*season + plant_group +
                    (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offgeneral3<- lmer(mean_offset ~ lang:season + season + plant_group +
                    (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offgeneral4<- lmer(mean_offset ~ lang*season*plant_group +
                     (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offgeneral5<- lmer(mean_offset ~  plant_group +
                     (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)
offgeneral6<- lmer(mean_offset ~ lang:season + lang + 
                     (1|study) -1, data=natdata0, na.action = "na.omit",REML = FALSE)


AIC(offgeneral1,offgeneral2, offgeneral3, offgeneral4, offgeneral5, offgeneral6)
#yes, it is

offgeneral<- lmer(mean_offset ~ lang:season + lang + plant_group +
                    (1|study) -1, data=natdata0, na.action = "na.omit",REML = TRUE)

rm(offgeneral1,offgeneral2, offgeneral3, offgeneral4, offgeneral5, offgeneral6)

summary(offgeneral)
anova(offgeneral)
myresplot(offgeneral, natdata0) 
hist(natdata0$mean_offset) 


visreg(offgeneral, xvar = "lang", by= "season",
       scale = "response", gg=TRUE) #looks similar for dry and wet seasons. less humidity, more offset (makes sense)
visreg(offgeneral, xvar = "lang", by= "plant_group",
       scale = "response", gg=TRUE) #same patterns
visreg(offgeneral, xvar = "lang", 
       scale = "response", gg=TRUE) #same pattern
#unestandarized plot (rawdata)
#ggplot(data=natdata, aes(x=lang, y=mean_offset)) +
  #geom_point(alpha=0.2) +
  #geom_smooth(method=lm,se=F)  +
  #facet_grid(season ~ .) #plotear plant_group, lang y season a la vez???

#let's see what happens in each climate class.
#also lets see whats work better. climatic best model + plant_ group or offgeneral applied for class

######offset warm######

offset_warm1 <- lmer(mean_offset ~ mat:season + map:season + season + plant_group +
                       (1|study) -1, data=natdata_warm0, na.action = "na.omit", REML = FALSE)
offset_warm2<- lmer(mean_offset ~ lang:season + lang + plant_group +
                    (1|study) -1, data=natdata_warm0, na.action = "na.omit",REML = FALSE)

AIC(offset_warm1, offset_warm2) #offset_warm1 works better

offset_warm<- lmer(mean_offset ~ lang:season + lang + plant_group +
                      (1|study) -1, data=natdata_warm0, na.action = "na.omit",REML = TRUE)

##### offset arid#######
offset_arid1<- lmer(SWLslope ~  map:season + plant_group +
                      (1|study) -1, data=natdata_arid0, na.action = "na.omit", REML = FALSE)
offset_arid2<- lmer(mean_offset ~ lang:season + lang + plant_group +
                      (1|study) -1, data=natdata_arid0, na.action = "na.omit",REML = FALSE)

AIC(offset_arid1, offset_arid2) #first one much better

offset_arid<- lmer(mean_offset ~ lang:season + lang + plant_group +
                     (1|study) -1, data=natdata_arid0, na.action = "na.omit",REML = TRUE)

#####offset tropical#####
offset_tropical1<- lmer(SWLslope ~ mat:season + season + plant_group +
                          (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = FALSE)
offset_tropical2<- lmer(mean_offset ~ lang:season + lang + plant_group +
                      (1|study) -1, data=natdata_trop0, na.action = "na.omit",REML = FALSE)

AIC(offset_tropical1, offset_tropical2) #first one much better

offset_tropical<- lmer(SWLslope ~ mat:season + season + plant_group +
                          (1|study) -1, data=natdata_trop0, na.action = "na.omit", REML = TRUE)

#####offset cold#######

offset_cold1<- lmer(SWLslope ~ lang:season + lang + plant_group +
                      (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)
offset_cold2<- lmer(mean_offset ~ lang:season + lang + plant_group +
                      (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = FALSE)

AIC(offset_cold1, offset_cold2) #first one much better

offset_cold<- lmer(SWLslope ~ lang:season + lang + plant_group +
                      (1|study) -1, data=natdata_cold0, na.action = "na.omit",REML = TRUE)

par(mfrow=c(2,2))
x11()

#no se muy bien como plotear esto último
