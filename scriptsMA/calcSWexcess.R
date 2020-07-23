source("scriptsMA/calcSWL.R")
drive_download("https://drive.google.com/file/d/1cAIXdvp2W7jFAMwwHZJm_YYKLVPGKuhp/view?usp=sharing",
               type = 'csv', path = 'dataMA/plant_sample.csv', overwrite = T)
# if this step does not work in your local computer you have a problem with ; and ,
plant <- read.csv('dataMA/plant_sample.csv', sep = ",")
plant$id_date_plot <- paste0(plant$author, '_', plant$year, 
                             '_', plant$date, '_', plant$plot)
plant <- dplyr::left_join(plant, source_results, by = 'id_date_plot')
# calculate soil water excess eq. 1 in Barbeta et al. 2019 HESS
plant$SWexcess <- plant$d2H_permil_plant - plant$slope * plant$d18O_permil_plant - plant$intercept
