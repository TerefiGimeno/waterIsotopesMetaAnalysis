source('scriptsMA/new_generatemodeldata.R')

#netCDF descargado de  https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview

library(sf)
library(ncdf4)
library(tidyverse)
library(lubridate)
library(tidync)
library(raster)
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
temp <- temp[, c('studyPlot', 'log', 'lat', 'date', 'temp_C')]
temp$month <- lubridate::month(temp$date, label = TRUE)
temp$year <- lubridate::year(temp$date)
temp$fyear <- as.factor(clean$year)

# explore correlation of ERA5 data with mat from meta-data
tempSumm <- summarise(group_by(temp, studyPlot), temp = mean(temp_C, na.rm = T))

meta$studyPlot <- paste0(meta$author, '-', meta$year, '-', meta$plotR)
metaShort <- rmDup(meta, 'studyPlot')
metaShort <- merge(metaShort, tempSumm, by = 'studyPlot', all.x = F, all.y = F)
windows(12,8)
par(mfrow=c(1,2))
with(metaShort, plot(mat ~ temp, pch = 19))
abline(lm(metaShort$mat ~ metaShort$temp))
legend('topleft', legend = c('R2 = 0.33', 'P < 0.001', 'slope: 0.48 (0.08)'), bty = 'n')
with(metaShort, plot(matWC ~ temp, pch = 19))
abline(lm(metaShort$matWC ~ metaShort$temp))
legend('topleft', legend = c('R2 = 0.49', 'P < 0.001', 'slope: 0.6 (0.05)'), bty = 'n')

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
  