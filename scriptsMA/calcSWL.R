swl <- read.csv("dataMA/data_SWLN.csv")
swl$id <- paste0(swl$author, '_', swl$year, '_', swl$journal)
swl <- doBy::orderBy(~id, swl)

swl_results <- data.frame(row.names = 1:length(unique(swl$id)))
swl_results$id <- unique(swl$id)

swlL <- list()
for(i in 1:nrow(swl_results)){
  swlL[[i]] <- subset(swl, swl$id == swl_results$id[i])
}

swlFits <- list()
for(i in 1:length(swlL)){
  swlFits[[i]] <- summary(lm(d2H_permil ~ d18O_permil, data = swlL[[i]]))
}

for (i in 1:nrow(swl_results)){
  swl_results$intercept <- swlFits[[i]]$coeff[1]
  swl_results$int_se <- swlFits[[i]]$coeff[3]
  swl_results$int_t <- swlFits[[i]]$coeff[5]
  swl_results$int_p <- swlFits[[i]]$coeff[7]
  swl_results$slope <- swlFits[[i]]$coeff[2]
  swl_results$slp_se <- swlFits[[i]]$coeff[4]
  swl_results$slp_t <- swlFits[[i]]$coeff[6]
  swl_results$slp_p <- swlFits[[i]]$coeff[8]
  swl_results$r2 <- swlFits[[i]]$r.squared
  swl_results$N <- nrow(swlL[[i]])
}

