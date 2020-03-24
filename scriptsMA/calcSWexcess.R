source("scriptsMA/calcSWL.R")
plantWI <- read.csv("dataMA/data_plant_water.csv")
plantWI$id <- paste0(plantWI$author, '_', plantWI$year, '_', plantWI$journal)
plantWI <- dplyr::left_join(plantWI, swl_results[, c('id', 'intercept', 'slope')], by = 'id')
# calculate soil water excess eq. 1 in Barbeta et al. 2019 HESS
plantWI$SWexcess <- plantWI$d2H_permil - plantWI$slope * plantWI$d18O_permil - plantWI$intercept
