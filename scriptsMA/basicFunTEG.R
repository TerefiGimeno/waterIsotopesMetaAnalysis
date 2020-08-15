countNA <- function(x){
  l <- length(x)-length(which(is.na(x)))
  return(l)
}

lengthWithoutNA <- function(x){
  l <- length(which(!is.na(x)))
  return(l)
}

rmDup <- function(dfr, whichvar){
  dfr <- dfr[!duplicated(dfr[,whichvar]),]
  return(dfr)
}
s.err <- function(x){sd(x)/sqrt(length(x))}

s.err.na <- function(x){
  se <- sd(x, na.rm = TRUE)/sqrt(length(which(!is.na(x))))
  return(se)
}

coef.var <- function(x){sd(x)/mean(x)}

calcVPD <- function(temp, RH){0.61365 * exp(17.502 * temp/(240.97 + temp)) * (1 -(RH/100))}

sum.na <- function(x){sum(x, na.rm = TRUE)}

mean.na <- function(x){mean(x, na.rm = T)}