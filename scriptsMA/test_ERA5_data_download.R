t2mERA5_1<- brick(x = 'ERA5_data/adaptor.mars.internal-1607942638.2488568-25061-33-689a8c97-ec8f-48cd-8e18-804ebfc38d4e.nc',
                  varname="t2m")
t2mERA5data_1<-raster::extract(t2mERA5_1,modeldataxy)
t2mERA5data_1 <- cbind(modeldataxy1, as.data.frame(t2mERA5data_1))

test <-  t2mERA5data_1 %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "dateOdd",
    values_to = "temp_2m"
  )
test$temp_C <- test$temp_2m - 273.15
test[, c('crap', 'date')] <- str_split_fixed(test$dateOdd, 'X', 2)
test$date <- lubridate::ymd(test$date)
clean <- test[, c('studyPlot', 'log', 'lat', 'date', 'temp_C')]
clean$month <- lubridate::month(test$date, label = TRUE)

ggplot(data=clean,aes(x=lat,y=temp_C))+
  geom_point()+
  facet_wrap(~month)