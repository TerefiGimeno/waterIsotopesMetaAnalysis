#split in groups to see the plots more clearly
offset <- doBy::order_by(~campaign, data = offset)
campNames <- data.frame(row.names = 1:length(unique(offset$campaign)))
campNames$campaign <- unique(offset$campaign)
campNames <- doBy::orderBy(~campaign, data = campNames)
campNames$crapNumber <- c(1:nrow(campNames))
lmwlswlplot <- left_join(offset, campNames, by = 'campaign')
lmwlswl <- lmwlswlplot %>%
  select(campaign, crapNumber, estimate, estimate.slope, intercept_LMWL, slope_LMWL) %>%
  unique

lmwlswlplotL <- list()
for(i in 1:ceiling((nrow(campNames)/20))){
  lmwlswlplotL[[i]] <- lmwlswlplot[which(lmwlswlplot$crapNumber >= i*20-19 & lmwlswlplot$crapNumber <= i*20), ]
}
lmwlswlL <- list()
for(i in 1:ceiling((nrow(campNames)/20))){
  lmwlswlL[[i]] <- lmwlswl[which(lmwlswl$crapNumber >= i*20-19 & lmwlswl$crapNumber <= i*20), ]
}

windows(12, 8)
#enter numbers from 1 to 9 where it says "i" to see batches of 20 plots
ggplot(data=lmwlswlplotL[[20]], aes(x = d18O_permil_plant, y = d2H_permil_plant))+
  geom_point()+
  facet_wrap(~campaign)+
  geom_abline(aes(slope = slope_LMWL, intercept = intercept_LMWL), lmwlswlL[[20]], col = 'blue')+
  geom_abline(aes(slope = estimate.slope, intercept = estimate), lmwlswlL[[20]], col = 'red')
  
