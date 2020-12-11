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
modeldataxy1 <- doBy::orderBy(~studyPlot, modeldataxy1)
modeldataxy <- as.data.frame(modeldataxy1[, c('log', 'lat')])
coordinates(modeldataxy) <- ~ log + lat

t2mERA5data_1<-raster::extract(t2mERA5,modeldataxy)
write.csv(t2mERA5data_1, file ='crap.csv', row.names = F)
write.csv(modeldataxy1, file ='crap.names.csv', row.names = F)
