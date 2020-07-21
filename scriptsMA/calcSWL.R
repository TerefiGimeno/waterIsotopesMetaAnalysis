library(googledrive)
drive_download("https://drive.google.com/file/d/1_bzWTiOK2JXatuYI4ovK7kXoSXEsRCKy/view?usp=sharing",
               type = 'csv', path = 'dataMA/fuenteprueba.csv')
source <- read.csv('dataMA/fuenteprueba.csv', sep = ";")[, 1:15]
names(source)[1] <- 'author'
source$id <- paste0(source$author, '_', source$year, '_', source$journal)
source$id_date <- paste0(source$author, '_', source$year, '_', source$journal, '_', source$date)
source <- doBy::orderBy(~id_date, source)

source_results <- data.frame(row.names = 1:length(unique(source$id_date)))
source_results$id_date <- unique(source$id_date)

sourceL <- list()
for(i in 1:nrow(source_results)){
  sourceL[[i]] <- subset(source, source$id_date == source_results$id_date[i])
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

