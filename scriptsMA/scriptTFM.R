##########Load data########################
library(googledrive)

drive_download("https://docs.google.com/spreadsheets/d/1Z7HathxScqACg5vBxWm5dBTbiYMsRKGFDWu2okI2JLE/edit?usp=sharing",
               type = 'csv', path = 'dataMA/meta_sample.csv', overwrite = T)
meta <- read.csv('dataMA/meta_sample.csv', sep = ",")

drive_download("https://docs.google.com/spreadsheets/d/1oRfLFvQthr_ZxMYjOSYpnmdbclQgxbz2waraHTafJ5U/edit#gid=0",
               type = 'csv', path = 'dataMA/source_sample.csv', overwrite = T)
source <- read.csv('dataMA/source_sample.csv', sep = ",")

drive_download("https://docs.google.com/spreadsheets/d/1v2zOeuB-b_BdiEQUMSHTz9nmm9pgRseUygD7SDmDgpo/edit?usp=sharing",
               type = 'csv', path = 'dataMA/plant_sample.csv', overwrite = T)
plant <- read.csv('dataMA/plant_sample.csv', sep = ",")

###############################################
library(tidyverse)
library(ggpubr)
library(broom)
library(tidyr)

# I fixed this problem, no need to run these lines any more
# source$d2H_permil_source<-as.numeric(as.character(source$d2H_permil_source))   #for some reason those are factors...
# source$d18O_permil_source<-as.numeric(as.character(source$d18O_permil_source))
source$year<-as.factor(source$year)
source$authorYear <- paste0(source$author, '_', source$year)
source$authorYearDate <- paste0(source$author, '_', source$year, '_', source$date)
# create variable id.campgain
source$id.campaign <- paste0(source$authorYear, '_', source$plotR)
source$campaignDate <- paste0(source$id.campaign, '_', source$date)
hist(source$d18O_permil_source)
# there are wierd values for the paper by Yonggang Chi in Ecohydrology CHECK

str(source)  ###check the type of variable

# cleaned already
# source<-subset(source,!d2H_permil_source=="NA")
# source<-subset(source,!d18O_permil_source=="NA") ##clean NA

source<-subset(source,label_pool=="bulk") ##only bulk soil

multiple<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-campaignDate, -species_source) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source, na.action='na.omit',data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

multiple<-multiple[,c('campaignDate', 'species_source', 'term','estimate','std.error',"statistic",'p.value')]  ###select only relevant columns

###more parameters from the regression, including r squared
rsquared<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-campaignDate, -species_source) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source,data = .x)),
    tidied = map(fit, glance)
  ) %>% 
  unnest(tidied)

rsquared<-rsquared[,c('campaignDate', 'species_source','r.squared')]  ###select only relevant columns

###count the number of soil samples, in case we want to cut off using this
length<- source %>% count(campaignDate, species_source)

###divide into intercept and slope
intercept<-subset(multiple,term=='(Intercept)')
slope<-subset(multiple,term=='d18O_permil_source')

hist(slope$estimate)

weird<-subset(slope,estimate<0)  ##the slope should be positive...anyway those got cut when using the p.value, except Twining

swl<-merge(intercept,slope,by=c('campaignDate', 'species_source'))   ###create a new table with separate columns for intercept and slope

###give proper names
colnames(swl)<-c("campaignDate","species_source","term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
                 "statistic.slope","p.value.slope")
##paste rsquareds
swl<-merge(swl,rsquared,by=c('campaignDate', 'species_source'))
source<-merge(source,rsquared,by=c('campaignDate', 'species_source'))
#past N's
swl<-merge(swl,length,by=c('campaignDate', 'species_source'))

hist(swl$r.squared)

swl <- subset(swl, p.value.slope < 0.05 & n > 2)   #cutoff non-significant and poorly fit regressions

library(ggpubr)

# add code to incorporate r.squared values
# something is off, we are retaining regressions with n < 2 and non significant slopes
sourceswl<- merge(swl,source,by=c('campaignDate', 'species_source'))
names(sourceswl)
ggplot(data=subset(sourceswl,r.squared<0.65&p.value.slope<0.05&n>2),aes(x=d18O_permil_source,y=d2H_permil_source))+
  geom_point()+
  geom_smooth(method=lm,se=F)+
  facet_wrap(~campaignDate)

ggplot(data=subset(sourceswl,rsquared<0.65),aes(x=d18O_permil_source,y=d2H_permil_source))+
  geom_point()+
  geom_smooth(method=lm,se=F)+
  facet_wrap(~campaignDate)+
  stat_cor()


###metadata filled information (MAP and MAT from NA's)

library(raster)
library(sp)

r <- getData("worldclim",var="bio",res=10)
r <- r[[c(1,12)]]
names(r) <- c("MATwc","MAPwc")
# ignore glasshouse studies
meta$lat <- ifelse(meta$lat > 91 | meta$lat < -91, NA, meta$lat)
meta$log <- ifelse(meta$log > 181 | meta$log < -181, NA, meta$log)
coords <- data.frame(x=meta$log,y=meta$lat)
# filter for NA's
coords <- coords[which(!is.na(coords$x)), ]
# remove duplicate values
source('scriptsMA/basicFunTEG.R')
coords <- rmDup(coords, 'x')

points <- SpatialPoints(coords, proj4string = r@crs)

values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)
points <- spsample(as(r@extent, 'SpatialPolygons'),n=100, type="random")   
values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)
# this works but does not fill in all the values and I don't know why!