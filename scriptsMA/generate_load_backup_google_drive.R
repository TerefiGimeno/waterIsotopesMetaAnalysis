library(googledrive)
library(googlesheets4)
library(lubridate)
meta_data <- read_sheet("https://docs.google.com/spreadsheets/d/1Z7HathxScqACg5vBxWm5dBTbiYMsRKGFDWu2okI2JLE/edit#gid=0")
metaToday <- paste0(today(), '_meta_data.csv')
write.csv(meta_data, file = paste0('meta_data/', metaToday), row.names = F)
drive_upload(paste0('meta_data/', metaToday),
             path = paste0('backup_meta_analysis/meta_backup/', metaToday))

plant_data <- read_sheet("https://docs.google.com/spreadsheets/d/1v2zOeuB-b_BdiEQUMSHTz9nmm9pgRseUygD7SDmDgpo/edit#gid=0")
plantToday <- paste0(today(), '_plant_data.csv')
write.csv(plant_data, file = paste0('plant_water/', plantToday), row.names = F)
drive_upload(paste0('plant_water/', plantToday),
             path = paste0('backup_meta_analysis/plant_backup/', plantToday))

other_vars_data <- read_sheet("https://docs.google.com/spreadsheets/d/1-2QE95kF_G4uugl2BQylzqFfPVpOmTO0_m57Ke6Engw/edit#gid=0")
otherVarsToday <- paste0(today(), '_other_vars_data.csv')
write.csv(other_vars_data, file = paste0('other_vars/', otherVarsToday), row.names = F)
drive_upload(paste0('other_vars/', otherVarsToday),
             path = paste0('backup_meta_analysis/other_vars_backup/', otherVarsToday))

source_data <- read_sheet("https://docs.google.com/spreadsheets/d/1udiv5w7YXCZXP2Pjnpshkzxsx5dxV1slG0OX7eztVq4/edit#gid=0")
sourceToday <- paste0(today(), '_source_data.csv')
write.csv(sourceToday, file = paste0('source_water/', sourceToday), row.names = F)
drive_upload(paste0('source_water/', sourceToday),
             path = paste0('backup_meta_analysis/source_meta/', sourceToday))
