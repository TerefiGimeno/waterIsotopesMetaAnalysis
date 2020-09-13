source('scriptsMA/generate_modeldata.R')

myresplot <- function(x,df) { #Home-made residual plot
  e <- resid(x, type="pearson")
  f <- fitted(x)
  obs <- x@frame[,1]
  par(mfrow=c(1,2), pty="m",mar=c(3,3,1,1),mgp = c(2, 1, 0), ps=12, bty="l")
  plot(f,e, pch=1, cex=0.8, col=alpha("black",0.3),
       xlab="Fitted values", ylab="Pearson residuals")
  abline(h = 0, lty = 2, col="black", lw=1)
  plot(f ~ obs, xlab="Observed values", ylab="Fitted values", type="n")
  points(obs,f, pch=1, cex=0.8, col=alpha("black",0.3))
  abline(a = 0, b = 1, lty = 2, col="Black", lw=1)
}

library(lme4)
library(tidyverse)
library(GGally)
library(MuMIn)
library(broom.mixed)
library(visreg)
library(performance)
# load this package to get p-values
library(lmerTest)
library(emmeans)


natdata <- subset(modeldata, natural== 'natural') #for climatic 
natdata0<-natdata
natdata0[which(natdata0$season == "not applicable"), 'season'] <- NA
################ A) models with study alone as random factor#################

#SWL

swlA0 <-lmer(SWLslope ~ (1|study) , data=modeldata,
                na.action = "na.omit",REML = FALSE)

swlA1<-lmer(SWLslope ~ climate_class*season + (1|study) , data=natdata,
                na.action = "na.omit",REML = FALSE)

swlA2 <-lmer(SWLslope ~ climate_class + (1|study) , data=natdata,
                na.action = "na.omit",REML = FALSE)

swlA3 <-lmer(SWLslope ~ season + (1|study) , data=natdata,
                na.action = "na.omit",REML = FALSE)

swlA4 <-lmer(SWLslope ~ season*lang + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA5 <-lmer(SWLslope ~ season*map + season*mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA6 <-lmer(SWLslope ~ season*map + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA7 <-lmer(SWLslope ~ season*mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA8 <-lmer(SWLslope ~ map*mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA9 <-lmer(SWLslope ~ map + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA10 <-lmer(SWLslope ~ mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlA11 <-lmer(SWLslope ~ lang + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

AIC(swlA0,swlA1,swlA2,swlA3,swlA4,swlA5,swlA6,swlA7,swlA8,swlA9,swlA10,swlA11)

swlAbest1<-lmer(SWLslope ~ season*lang + (1|study) , data=natdata,
                  na.action = "na.omit",REML = TRUE)

swlAbest2<-lmer(SWLslope ~ map + (1|study) , data=natdata,
                 na.action = "na.omit",REML = TRUE)
  
swlAbest3<-lmer(SWLslope ~ season*lang + (1|study) , data=natdata,
                 na.action = "na.omit",REML = TRUE)
  
AIC(swlAbest1,swlAbest2,swlAbest3)  

swlAbest<-lmer(SWLslope ~ map + (1|study) , data=natdata,
               na.action = "na.omit",REML = TRUE)

summary(swlAbest)

#offset climate

offclimA0 <-lmer(mean_offset ~ (1|study) , data=modeldata,
             na.action = "na.omit",REML = FALSE)

offclimA1<-lmer(mean_offset ~ climate_class*season + (1|study) , data=natdata,
            na.action = "na.omit",REML = FALSE)

offclimA2 <-lmer(mean_offset ~ climate_class + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA3 <-lmer(mean_offset ~ season + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA4 <-lmer(mean_offset ~ season*lang + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA5 <-lmer(mean_offset ~ season*map + season*mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA6 <-lmer(mean_offset ~ season*map + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA7 <-lmer(mean_offset ~ season*mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA8 <-lmer(mean_offset ~ map*mat + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA9 <-lmer(mean_offset ~ map + (1|study) , data=natdata,
             na.action = "na.omit",REML = FALSE)

offclimA10 <-lmer(mean_offset ~ mat + (1|study) , data=natdata,
              na.action = "na.omit",REML = FALSE)

offclimA11 <-lmer(mean_offset ~ lang + (1|study) , data=natdata,
              na.action = "na.omit",REML = FALSE)

offclimA12<-lmer(mean_offset ~ SWLslope + (1|study) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclimA13<-lmer(mean_offset ~ SWLslope*season + (1|study) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

AIC(offclimA0,offclimA1,offclimA2,offclimA3,offclimA4,offclimA5,offclimA6,
    offclimA7,offclimA8,offclimA9,offclimA10,offclimA11, offclimA12, offclimA13)

offclimAbest<-lmer(mean_offset ~ SWLslope*season + (1|study) , data=natdata,
                    na.action = "na.omit",REML = TRUE)

summary(offclimAbest)
emmeans(offclimAbest, list(pairwise ~ SWLslope*season), adjust = "tukey") # i dont really know how to interpret this

#offset pft

offplantgroupA<-lmer(mean_offset ~plant_group  +
                     (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

offleafhabitA<-lmer(mean_offset ~leaf_habit  + (1|study),
                    data= subset(modeldata[which(modeldata$woodiness!='non-woody'),], leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
                    na.action = "na.omit",REML = TRUE)

offwoodA<-lmer(mean_offset ~woodiness  +
                 (1|study), data=modeldata, na.action = "na.omit",REML = TRUE)


offleafshapeA<-lmer(mean_offset ~leaf_shape  +
                      (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)


AIC(offplantgroupA,offleafhabitA,offleafshapeA,offwoodA) #shape and group. Group has more ecological logic

offpftA<-lmer(mean_offset ~plant_group  +
                (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

#ofset general

offA<-lmer(mean_offset ~plant_group  + SWLslope*season +
             (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

############# B) models with study and season (with not applicble values) as random factor##########

#SWL

swlB0 <-lmer(SWLslope ~ (1|study/season) , data=modeldata,
             na.action = "na.omit",REML = FALSE)

swlB1 <-lmer(SWLslope ~ climate_class + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlB2 <-lmer(SWLslope ~ map*mat + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlB3 <-lmer(SWLslope ~ map + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swlB4 <-lmer(SWLslope ~ mat + (1|study/season) , data=natdata,
              na.action = "na.omit",REML = FALSE)

swlB5 <-lmer(SWLslope ~ lang + (1|study/season) , data=natdata,
              na.action = "na.omit",REML = FALSE)

AIC(swlB0,swlB1,swlB2,swlB3,swlB4,swlB5)

swlBbest<-lmer(SWLslope ~ map + (1|study/season) , data=natdata,
               na.action = "na.omit",REML = TRUE)

summary(swlBbest)

#offset climate

offclimB0 <-lmer(mean_offset ~ (1|study/season) , data=modeldata,
                 na.action = "na.omit",REML = FALSE)

offclimB1 <-lmer(mean_offset ~ climate_class + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclimB2 <-lmer(mean_offset ~ map*mat + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclimB3 <-lmer(mean_offset ~ map + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclimB4 <-lmer(mean_offset ~ mat + (1|study/season) , data=natdata,
                  na.action = "na.omit",REML = FALSE)

offclimB5 <-lmer(mean_offset ~ lang + (1|study/season) , data=natdata,
                  na.action = "na.omit",REML = FALSE)

offclimB6<-lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)


AIC(offclimB0,offclimB1,offclimB2,offclimB3,offclimB4,offclimB5,offclimB6)

offclimBbest<-lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata,
                    na.action = "na.omit",REML = FALSE)

summary(offclimBbest)

#offset pft

offplantgroupB<-lmer(mean_offset ~plant_group  +
                       (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

offleafhabitB<-lmer(mean_offset ~leaf_habit  + (1|study/season),
                    data= subset(modeldata[which(modeldata$woodiness!='non-woody'),], leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
                    na.action = "na.omit",REML = TRUE)

offwoodB<-lmer(mean_offset ~woodiness  +
                 (1|study/season), data=modeldata, na.action = "na.omit",REML = TRUE)


offleafshapeB<-lmer(mean_offset ~leaf_shape  +
                      (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)


AIC(offplantgroupB,offleafhabitB,offleafshapeB,offwoodB) #shape and group. Group has more ecological logic

offpftB<-lmer(mean_offset ~plant_group  +
                (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)


#ofset general

offB<- lmer(mean_offset ~plant_group  + SWLslope +
              (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)



############# C) models with study and season (without not applicble values) as random factor##########

#SWL

swlC0 <-lmer(SWLslope ~ (1|study/season) , data=subset(modeldata, season != 'not applicable'),
             na.action = "na.omit",REML = FALSE)

swlC1 <-lmer(SWLslope ~ climate_class + (1|study/season) , data=subset(natdata, season != 'not applicable'),
             na.action = "na.omit",REML = FALSE)

swlC2 <-lmer(SWLslope ~ map*mat + (1|study/season) , data=subset(natdata, season != 'not applicable'),
             na.action = "na.omit",REML = FALSE)

swlC3 <-lmer(SWLslope ~ map + (1|study/season) , data=subset(natdata, season != 'not applicable'),
             na.action = "na.omit",REML = FALSE)

swlC4 <-lmer(SWLslope ~ mat + (1|study/season) , data=subset(natdata, season != 'not applicable'),
             na.action = "na.omit",REML = FALSE)

swlC5 <-lmer(SWLslope ~ lang + (1|study/season) , data=subset(natdata, season != 'not applicable'),
             na.action = "na.omit",REML = FALSE)

AIC(swlC0,swlC1,swlC2,swlC3,swlC4,swlC5)

swlCbest1<-lmer(SWLslope ~ map + (1|study/season) , data=subset(natdata, season != 'not applicable'),
               na.action = "na.omit",REML = TRUE)
swlCbest2<-lmer(SWLslope ~ lang+ (1|study/season) , data=subset(natdata, season != 'not applicable'),
                na.action = "na.omit",REML = TRUE)

AIC(swlCbest1, swlCbest2) #both are similar



#offset climate

offclimC0 <-lmer(mean_offset ~ (1|study/season) , data=modeldata,
                 na.action = "na.omit",REML = FALSE)

offclimC1 <-lmer(mean_offset ~ climate_class + (1|study/season) , data=subset(natdata, season != 'not applicable'),
                 na.action = "na.omit",REML = FALSE)

offclimC2 <-lmer(mean_offset ~ map*mat + (1|study/season) , data=subset(natdata, season != 'not applicable'),
                 na.action = "na.omit",REML = FALSE)

offclimC3 <-lmer(mean_offset ~ map + (1|study/season) , data=subset(natdata, season != 'not applicable'),
                 na.action = "na.omit",REML = FALSE)

offclimC4 <-lmer(mean_offset ~ mat + (1|study/season) , data=subset(natdata, season != 'not applicable'),
                 na.action = "na.omit",REML = FALSE)

offclimC5 <-lmer(mean_offset ~ lang + (1|study/season) , data=subset(natdata, season != 'not applicable'),
                 na.action = "na.omit",REML = FALSE)

offclimC6<-lmer(mean_offset ~ SWLslope + (1|study/season) , data=subset(natdata, season != 'not applicable'),
                na.action = "na.omit",REML = FALSE)


AIC(offclimC0,offclimC1,offclimC2,offclimC3,offclimC4,offclimC5,offclimC6)

offclimCbest1<-lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata,
                   na.action = "na.omit",REML = TRUE)

offclimCbest2<-lmer(mean_offset ~ lang + (1|study/season) , data=natdata,
                    na.action = "na.omit",REML = TRUE)

AIC(offclimCbest1, offclimCbest2)

offclimCbest<- lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata,
                    na.action = "na.omit",REML = TRUE)

summary(offclimCbest)

#offset pft

offplantgroupC<-lmer(mean_offset ~plant_group  +
                       (1|study/season), data=subset(modeldata, !woodiness=="non-woody", season != 'not applicable'), na.action = "na.omit",REML = TRUE)

offleafhabitC<-lmer(mean_offset ~leaf_habit  + (1|study/season),
                    data= subset(modeldata[which(modeldata$woodiness!='non-woody'),],leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
                    na.action = "na.omit",REML = TRUE) #struggling with this one


offwoodC<-lmer(mean_offset ~woodiness  +
                 (1|study/season), data=subset(modeldata, season != 'not applicable'), na.action = "na.omit",REML = TRUE)


offleafshapeC<-lmer(mean_offset ~leaf_shape  +
                      (1|study/season), data=subset(modeldata, !woodiness=="non-woody", season != 'not applicable'), na.action = "na.omit",REML = TRUE)


AIC(offplantgroupC,offleafhabitC,offleafshapeC,offwoodC) #wood and by far 8strange how it changed...suspicious)

offpftC<-lmer(mean_offset ~woodiness  +
                (1|study/season), data=subset(modeldata, season != 'not applicable'), na.action = "na.omit",REML = TRUE)

#ofset general

offC<- lmer(mean_offset ~woodiness  + SWLslope +
              (1|study/season), data=subset(natdata, season != 'not applicable'), na.action = "na.omit",REML = TRUE)


#####D) models as A but not applicable(season)= NA############

#SWL

swlD0 <-lmer(SWLslope ~ (1|study) , data=modeldata,
             na.action = "na.omit",REML = FALSE)

swlD1<-lmer(SWLslope ~ climate_class*season + (1|study) , data=natdata0,
            na.action = "na.omit",REML = FALSE)

swlD2 <-lmer(SWLslope ~ climate_class + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD3 <-lmer(SWLslope ~ season + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD4 <-lmer(SWLslope ~ season*lang + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD5 <-lmer(SWLslope ~ season*map + season*mat + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD6 <-lmer(SWLslope ~ season*map + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD7 <-lmer(SWLslope ~ season*mat + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD8 <-lmer(SWLslope ~ map*mat + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD9 <-lmer(SWLslope ~ map + (1|study) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlD10 <-lmer(SWLslope ~ mat + (1|study) , data=natdata0,
              na.action = "na.omit",REML = FALSE)

swlD11 <-lmer(SWLslope ~ lang + (1|study) , data=natdata0,
              na.action = "na.omit",REML = FALSE)

AIC(swlD0,swlD1,swlD2,swlD3,swlD4,swlD5,swlD6,swlD7,swlD8,swlD9,swlD10,swlD11)

swlAbest1<-lmer(SWLslope ~ season*lang + (1|study) , data=natdata0,
                na.action = "na.omit",REML = TRUE)

swlAbest2<-lmer(SWLslope ~ map + (1|study) , data=natdata0,
                na.action = "na.omit",REML = TRUE)

swlAbest3<-lmer(SWLslope ~ season*lang + (1|study) , data=natdata0,
                na.action = "na.omit",REML = TRUE)

AIC(swlAbest1,swlAbest2,swlAbest3)  

swlDbest<-lmer(SWLslope ~ season*map + (1|study) , data=natdata0,
               na.action = "na.omit",REML = FALSE)


summary(swlDbest)

#offset climate

offclimD0 <-lmer(mean_offset ~ (1|study) , data=modeldata,
                 na.action = "na.omit",REML = FALSE)

offclimD1<-lmer(mean_offset ~ climate_class*season + (1|study) , data=natdata0,
                na.action = "na.omit",REML = FALSE)

offclimD2 <-lmer(mean_offset ~ climate_class + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD3 <-lmer(mean_offset ~ season + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD4 <-lmer(mean_offset ~ season*lang + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD5 <-lmer(mean_offset ~ season*map + season*mat + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD6 <-lmer(mean_offset ~ season*map + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD7 <-lmer(mean_offset ~ season*mat + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD8 <-lmer(mean_offset ~ map*mat + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD9 <-lmer(mean_offset ~ map + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD10 <-lmer(mean_offset ~ mat + (1|study) , data=natdata0,
                  na.action = "na.omit",REML = FALSE)

offclimD11 <-lmer(mean_offset ~ lang + (1|study) , data=natdata0,
                  na.action = "na.omit",REML = FALSE)

offclimD12<-lmer(mean_offset ~ SWLslope + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimD13<-lmer(mean_offset ~ SWLslope*season + (1|study) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

AIC(offclimD0,offclimD1,offclimD2,offclimD3,offclimD4,offclimD5,offclimD6,
    offclimD7,offclimD8,offclimD9,offclimD10,offclimD11, offclimD12, offclimD13)

offclimDbest1<-lmer(mean_offset ~ SWLslope*season + (1|study) , data=natdata0,
                   na.action = "na.omit",REML = TRUE)
offclimDbest2<-lmer(mean_offset ~ season*map + (1|study) , data=natdata0,
                    na.action = "na.omit",REML = TRUE)

AIC(offclimDbest1,offclimDbest2)

offclimDbest<-lmer(mean_offset ~ SWLslope*season + (1|study) , data=natdata0,
                   na.action = "na.omit",REML = TRUE)

summary(offclimAbest)
emmeans(offclimAbest, list(pairwise ~ SWLslope*season), adjust = "tukey") # i dont really know how to interpret this

#offset pft

offplantgroupD<-lmer(mean_offset ~plant_group  +
                       (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

offleafhabitD<-lmer(mean_offset ~leaf_habit  + (1|study),
                    data= subset(modeldata[which(modeldata$woodiness!='non-woody'),], leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
                    na.action = "na.omit",REML = TRUE)

offwoodD<-lmer(mean_offset ~woodiness  +
                 (1|study), data=modeldata, na.action = "na.omit",REML = TRUE)


offleafshapeD<-lmer(mean_offset ~leaf_shape  +
                      (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)


AIC(offplantgroupD,offleafhabitD,offleafshapeD,offwoodD) #shape and group. Group has more ecological logic

offpftD<-lmer(mean_offset ~plant_group  +
                (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

#ofset general

offD<-lmer(mean_offset ~plant_group  + map*season +
             (1|study), data=subset(natdata0, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

##### E) models as C but not applicable(season)= NA############

#SWL

swlE0 <-lmer(SWLslope ~ (1|study/season) ,data=modeldata,
             na.action = "na.omit",REML = FALSE)

swlE1 <-lmer(SWLslope ~ climate_class + (1|study/season) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlE2 <-lmer(SWLslope ~ map*mat + (1|study/season) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlE3 <-lmer(SWLslope ~ map + (1|study/season) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlE4 <-lmer(SWLslope ~ mat + (1|study/season) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

swlE5 <-lmer(SWLslope ~ lang + (1|study/season) , data=natdata0,
             na.action = "na.omit",REML = FALSE)

AIC(swlE0,swlE1,swlE2,swlE3,swlE4,swlE5)

swlEbest1<-lmer(SWLslope ~ map + (1|study/season) , data=natdata0,
                na.action = "na.omit",REML = TRUE)
swlEbest2<-lmer(SWLslope ~ lang + (1|study/season) , data=natdata0,
                na.action = "na.omit",REML = TRUE)

AIC(swlEbest1, swlEbest2) #both are similar


#offset climate

offclimE0 <-lmer(mean_offset ~ (1|study/season) , data=modeldata,
                 na.action = "na.omit",REML = FALSE)

offclimE1 <-lmer(mean_offset ~ climate_class + (1|study/season) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimE2 <-lmer(mean_offset ~ map*mat + (1|study/season) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimE3 <-lmer(mean_offset ~ map + (1|study/season) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimE4 <-lmer(mean_offset ~ mat + (1|study/season) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimE5 <-lmer(mean_offset ~ lang + (1|study/season) , data=natdata0,
                 na.action = "na.omit",REML = FALSE)

offclimE6<-lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata0,
                na.action = "na.omit",REML = FALSE) 


AIC(offclimE0,offclimE1,offclimE2,offclimE3,offclimE4,offclimE5,offclimE6)

offclimEbest1<-lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata0,
                    na.action = "na.omit",REML = TRUE)

offclimEbest2<-lmer(mean_offset ~ lang + (1|study/season) , data=natdata0,
                    na.action = "na.omit",REML = TRUE)

AIC(offclimEbest1, offclimEbest2)

offclimEbest<- lmer(mean_offset ~ SWLslope + (1|study/season) , data=natdata0,
                    na.action = "na.omit",REML = TRUE)

summary(offclimCbest)

#offset pft

offplantgroupE<-lmer(mean_offset ~plant_group  +
                       (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

offleafhabitE<-lmer(mean_offset ~leaf_habit  + (1|study/season),
                    data= subset(modeldata[which(modeldata$woodiness!='non-woody'),],leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
                    na.action = "na.omit",REML = TRUE) 


offwoodE<-lmer(mean_offset ~woodiness  +
                 (1|study/season), data=subset(modeldata, season != 'not applicable'), na.action = "na.omit",REML = TRUE)


offleafshapeE<-lmer(mean_offset ~leaf_shape  +
                      (1|study/season), data=subset(modeldata, !woodiness=="non-woody", season != 'not applicable'), na.action = "na.omit",REML = TRUE)


AIC(offplantgroupC,offleafhabitC,offleafshapeC,offwoodC) #group and shape. again we chose group

offpftE<-lmer(mean_offset ~plant_group  +
                (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)
#ofset general

offE<- lmer(mean_offset ~plant_group  + SWLslope +
              (1|study/season), data=natdata0, na.action = "na.omit",REML = TRUE)

####### SUMMARY compare A, B, C, D and E##########

#SWL

AIC(swlAbest,swlBbest,swlCbest1,swlCbest2,swlDbest,swlEbest1,swlEbest2) #D > C(1,2),E(1,2) > B > A

summary(swlAbest)
summary(swlBbest)
summary(swlCbest1)
summary(swlCbest2)
summary(swlDbest)
summary(swlEbest1)
summary(swlEbest2)

#offset climate

AIC (offclimAbest, offclimBbest, offclimCbest, offclimDbest, offclimEbest) #(E>D)>A>(C,B)

summary(offclimAbest)
summary(offclimBbest)
summary(offclimCbest)
summary(offclimDbest)
summary(offclimEbest)

#offset pft

AIC(offpftA, offpftB, offpftC, offpftD, offpftE) #C > (A,B,D,E)

summary(offpftA) 
summary(offpftB)
summary(offpftC)
summary(offpftD)
summary(offpftE)

#ofset general 

AIC(offA,offB,offC,offD,offE) #D>>E>A > (B,C)

summary(offA)
summary(offB)
summary(offC)
summary(offD)
summary(offE)
