
###test of natural$position

natural$position[natural$mean_lcexcess<0&natural$mean_offset<0]<-"neg_neg"
natural$position[natural$mean_lcexcess>=0&natural$mean_offset>=0]<-"pos_pos"
natural$position[natural$mean_lcexcess>=0&natural$mean_offset<0]<-"pos_neg"
natural$position[natural$mean_lcexcess<0&natural$mean_offset>=0]<-"neg_pos"

###below-below

summary(model<-lmer(SWLslope~smwl1+(1|study),subset(natural,position=="neg_neg")))

summary(model<-lmer(mean_offset~smwl1+(1|study),subset(natural,position=="neg_neg")))
summary(model<-lmer(mean_offset~temp_C+(1|study),subset(natural,position=="neg_neg")))

ggplot(data=subset(natural,position=="neg_neg"),aes(x=temp_C,y=mean_offset))+
  geom_point(size=4,shape=21,aes(fill=smwl1))+
  stat_cor()+
  theme_bw()+
  geom_smooth(method=lm)+
  scale_fill_gradient2(high="blue",midpoint=0.2,low="red")

summary(model<-lmer(mean_offset~temp_C+(1|study),subset(natural,position=="neg_neg")))

###is the temprature effect mediated by its effect on SWL slope?

summary(model<-lmer(mean_offset~SWLslope+temp_C+(1|study),subset(natural,position=="neg_neg")))

summary(model<-lmer(mean_offset~SWLintercept+temp_C+smwl1+(1|study),subset(natural,position=="neg_neg")))

summary(model<-lmer(SWLintercept~temp_C+smwl1+(1|study),subset(natural,position=="neg_neg")))

summary(model<-lmer(SWLslope~temp_C+smwl1+(1|study),subset(natural,position=="neg_neg")))

performance::check_model(model)

ggplot(natural,aes(x=climate_class,y=SWLintercept))+
  geom_boxplot(aes(fill=climate_class))+
  geom_jitter(shape=21,aes(fill=climate_class))+
  theme_bw()

###it seems that the temperature effect on the mean_offset is created by the negative T effect on the SWL intercept
###it seems that the smwl1 effect on the mean_offset is created by the positive effect on the SWL intercept and slope!

###if the intercept is higher in wet and cold sites, then the SW-excess will be lower there?



###above-above

hist(subset(natural,position=="pos_pos")$SWLslope)

ggplot(natural,aes(x=position,y=mean_lcexcess))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

summary(model<-lmer(SWLslope~smwl1+(1|study),subset(natural,position=="neg_neg")))

###Above LMWL and below SWL

hist(subset(natural,position=="pos_neg")$SWLslope)

ggplot(natural,aes(x=position,y=mean_lcexcess))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

###Below LMWL and above SWL

summary(model<-lmer(smwl1~position+(1|study),natural))

ggplot(natural,aes(x=position,y=smwl1))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

summary(model<-lmer(temp_C~position+(1|study),natural))

ggplot(natural,aes(x=position,y=temp_C))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

summary(model<-lmer(mper~position+(1|study),natural))

ggplot(natural,aes(x=position,y=mper))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()

###

summary(model<-lmer(mean_offset~smwl1+temp_C+(1|study),subset(natural,position=="neg_neg")))

summary(model<-lmer(mean_offset~smwl1+temp_C+(1|study),subset(natural,position=="pos_neg")))

summary(model<-lmer(mean_offset~smwl1+temp_C+(1|study),subset(natural,position=="neg_pos")))

summary(model<-lmer(mean_offset~smwl1+temp_C+(1|study),subset(natural,position=="pos_pos")))


###LET'S SEE THE INTERDEPENDENCE OF SW-EXCESS-SWL SLOPE-SLW INTERCEPT-TO CLIMATIC VARIABLES
library(lavaan)

natural.num<- select_if(natural, is.numeric)

model.Lavaan <- 'mean_offset~smwl1+temp_C
SWLslope~smwl1+temp_C
SWLintercept~smwl1
mean_offset~~SWLslope
mean_offset~~SWLintercept'

fit1 <- lavaan:::cfa(model.Lavaan, data=natural.num,std.lv=TRUE)
summary(fit1,rsq=T,fit.measures=T)
varTable(fit1)
# Plot path diagram:
layout(t(1:1))
labels<-c("SW-excess","SWL slope","SWL intercept","Soil VWC","Air T")
semPlot::semPaths(fit1,title=T,curvePivot =TRUE,layout='tree',what="stand",label.cex=2,nodeLabels = labels)
title("SEM", line=3)
AIC(fit1)