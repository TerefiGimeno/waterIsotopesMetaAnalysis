library(googledrive)
drive_download("https://drive.google.com/file/d/1hP1V9eLKqUyRq9Zlfr6gzkpzSXVLuYtY/view?usp=sharing",
               type = 'csv', path = 'dataMA/source_sample.csv', overwrite = T)
# if this step does not work in your local computer you have a problem with ; and ,
source <- read.csv('dataMA/source_sample.csv', sep = ",")
source$id_species <- paste0(source$author, '_', source$year, '_', source$species)
source$id_species_date <- paste0(source$author, '_', source$year, '_', source$species, '_', source$date)
source$id_species_date_plot <- paste0(source$author, '_', source$year, '_', source$species,
                              '_', source$date, '_', source$plot)
source$id_date <- paste0(source$author, '_', source$year, '_', source$date)
source$id_date_plot <- paste0(source$author, '_', source$year,
                                      '_', source$date, '_', source$plot)

soil <- subset(source, label_class == 'soil')
source('scriptsMA/clean_stuff.R')
source_results <- data.frame(row.names = 1:length(unique(soil$id_date_plot)))
source_results$id_date_plot <- unique(soil$id_date_plot)


sourceL <- list()
for(i in 1:nrow(source_results)){
  sourceL[[i]] <- subset(soil, soil$id_date_plot == source_results$id_date_plot[i])
}

sourceFits <- list()
for(i in 1:length(sourceL)){
  sourceFits[[i]] <- summary(lm(d2H_permil_source ~ d18O_permil_source, data = sourceL[[i]]))
}

for (i in 1:nrow(source_results)){
  source_results$intercept[i] <- sourceFits[[i]]$coeff[1]
  source_results$int_se[i] <- sourceFits[[i]]$coeff[3]
  source_results$int_t[i] <- sourceFits[[i]]$coeff[5]
  source_results$int_p[i] <- sourceFits[[i]]$coeff[7]
  source_results$slope[i] <- sourceFits[[i]]$coeff[2]
  source_results$slp_se[i] <- sourceFits[[i]]$coeff[4]
  source_results$slp_t[i] <- sourceFits[[i]]$coeff[6]
  source_results$slp_p[i] <- sourceFits[[i]]$coeff[8]
  source_results$r2[i] <- sourceFits[[i]]$r.squared
  source_results$N[i] <- nrow(sourceL[[i]])
}

par(mfrow = c(3, 4))
for (i in 1: length(sourceL)){
  plot(d2H_permil_source ~ d18O_permil_source, data = sourceL[[i]],
       ylab ='d2H permil', xlab ='d18O permil', pch = 19, col = as.factor(myData$label_pool),
       main = sourceL[[i]][1, 'id_date_plot'])
  abline(sourceFits[[i]])
  }
