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


natdata <- subset(modeldata, natural== 'natural') ¡

######SWL#######

#null model
swlnull <-lmer(SWLslope ~ (1|study/season) , data=modeldata,
             na.action = "na.omit",REML = FALSE)



#best model for swl
swl1 <-lmer(SWLslope ~ climate_class + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swl2 <-lmer(SWLslope ~ map*mat + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swl3 <-lmer(SWLslope ~ map + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swl4 <-lmer(SWLslope ~ mat + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

swl5 <-lmer(SWLslope ~ lang + (1|study/season) , data=natdata,
             na.action = "na.omit",REML = FALSE)

AIC(swlnull,swl1,swl2,swl3,swl4,swl5)

swlbest<-lmer(SWLslope ~ map + (1|study/season) , data=natdata,
              na.action = "na.omit",REML = TRUE)

myresplot(swlbest, natdata) 
check_model(swlbest)
summary(swlbest)
visreg(swlbest, xvar = "map",
       scale = "response", gg=TRUE)

#offset climate


offclim1 <-lmer(mean_offset ~ climate_class + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclim2 <-lmer(mean_offset ~ map*mat + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclim3 <-lmer(mean_offset ~ map + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclim4 <-lmer(mean_offset ~ mat + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)

offclim5 <-lmer(mean_offset ~ lang + (1|study/season) , data=natdata,
                 na.action = "na.omit",REML = FALSE)


AIC(offclim1,offclim2,offclim3,offclim4,offclim5)

offclimbest<-lmer(mean_offset ~ map + (1|study/season) , data=natdata,
                   na.action = "na.omit",REML = TRUE)

check_model(offclimbest)
summary(offclimbest)
visreg(offclimbest, xvar = "map",
       scale = "response", gg=TRUE)

#offset pft

offplantgroup<-lmer(mean_offset ~plant_group  +
                       (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)
check_model(offplantgroup)
summary(offplantgroup)
visreg(offplantgroup, xvar = "plant_group",
       scale = "response", gg=TRUE)

offleafhabit<-lmer(mean_offset ~leaf_habit  + (1|study/season),
                    data= subset(modeldata[which(modeldata$woodiness!='non-woody'),], leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
                    na.action = "na.omit",REML = TRUE)

check_model(offleafhabit)
summary(offleafhabit)
visreg(offleafhabit, xvar = "leaf_habit",
       scale = "response", gg=TRUE)


offwood<-lmer(mean_offset ~woodiness  +
                 (1|study/season), data=modeldata, na.action = "na.omit",REML = TRUE)


check_model(offwood)
summary(offwood)
visreg(offwood, xvar = "woodiness",
       scale = "response", gg=TRUE)

offleafshape<-lmer(mean_offset ~leaf_shape  +
                      (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

check_model(offleafshape)
summary(offleafshape)
visreg(offleafshape)

AIC(offplantgroup,offleafhabit,offleafshape,offwood) # (REML=FALSE) shape and group. Group has more ecological logic


#ofset general and null

#null
offnull<-lmer(mean_offset ~ (1|study/season) , data=natdata,
              na.action = "na.omit",REML = TRUE)


summary(offnull)

#general

offgeneral1<- lmer(mean_offset ~plant_group  + map +
              (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral2<- lmer(mean_offset ~leaf_shape  + map +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral3<- lmer(mean_offset ~leaf_habit  + map +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral4<- lmer(mean_offset ~woodiness  + map +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral5<- lmer(mean_offset ~plant_group  + lang +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral6<- lmer(mean_offset ~leaf_shape  + lang +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral7<- lmer(mean_offset ~leaf_habit  + lang +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)
offgeneral8<- lmer(mean_offset ~woodiness  + lang +
                    (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = F)


AIC(offgeneral1,offgeneral2,offgeneral3,offgeneral4,offgeneral5,offgeneral6,offgeneral7,offgeneral8)

#Ä1 are 2 are the best optios, 1 is more logical

offgeneralbest<-lmer(mean_offset ~plant_group  + map +
                     (1|study), data=subset(natdata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)


check_model(offgeneralbest)
summary(offgeneralbest)

#best models off all

AIC(offnull,offgeneral1,offclim3,offplantgroup)

offbest<-offgeneralbest

summary(offbest)
visreg(offbest)

