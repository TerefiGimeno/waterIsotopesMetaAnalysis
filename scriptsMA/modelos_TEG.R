source('scriptsMA/generate_modeldata.R')

#######models (SWL)#######

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
# load this package to get p-values
library(lmerTest)
library(emmeans)

#Lets begin with SWL slope 

natdata <- subset(modeldata, natural== 'natural') #reject hydro-managed plots

#we dont want not applicable values
# natdata[which(natdata$season == "not applicable"), 'season'] <- NA

#general model
swlbest2 <-lmer(SWLslope ~ climate_class +
                 (1|study) , data=natdata,
               na.action = "na.omit",REML = TRUE)
summary(swlbest2)
anova(swlbest2)
emmeans(swlbest2, list(pairwise ~ climate_class), adjust = "tukey")
visreg(swlbest2, xvar = "climate_class",
       scale = "response", gg=TRUE)
# no significant differences in swl among climate classes

swlbest2 <-lmer(SWLslope ~ lang +
                  (1|study) , data=natdata,
                na.action = "na.omit",REML = TRUE)
myresplot(swlbest2, natdata)
summary(swlbest2)
plot(natdata$SWLslope ~ natdata$lang, pch = 19, col = as.factor(natdata$climate_class),
     ylab = 'log (slope SWL)', xlab = 'Lang Index (mm/c)')
legend('bottomright', legend = levels(as.factor(natdata$climate_class)), pch = 19, col = 1:4, bty = 'n')
plot(natdata$SWLslope ~ natdata$lang, pch = 19, col = as.factor(natdata$climate_class),
     ylab = 'log (slope SWL)', xlab = 'Lang Index (mm/c)', xlim = c(0, 450))

# important disequilibrium in sample size across climate_classes
doBy::summaryBy(SWLslope ~ season + climate_class, FUN = lengthWithoutNA, data = subset(natdata, season != 'not applicable'))
doBy::summaryBy(SWLslope ~ climate_class, FUN = lengthWithoutNA, data = natdata)

swlbest <-lmer(SWLslope ~ climate_class*season +
                 (1|study) , data=subset(natdata, season != 'not applicable'),
               na.action = "na.omit",REML = TRUE) #has to be true

summary(swlbest)
myresplot(swlbest, natdata) #looks good
anova(swlbest)
emmeans(swlbest, list(pairwise ~ climate_class*season), adjust = "tukey")
visreg(swlbest, xvar = "season", by = 'climate_class',
       scale = "response", gg=TRUE)
# no significant differences neither among climate classes, nor among seasons

swlLang <- lmer(SWLslope ~ lang*season + (1|study), data = subset(natdata, season != 'not applicable'),
                na.action = 'na.omit', REML = TRUE)
summary(swlLang)
# no significant effects of the Lang index or season
windows(8, 8)
par(mfrow=c(1, 1))
dry <- subset(natdata, season == 'dry')
plot(dry$SWLslope ~ dry$lang, pch = 19, col = as.factor(dry$climate_class),
     xlim = c(0, 600), ylim= c(0, 12), ylab = 'SWL slope', xlab = 'Lang Index (mm/C)')
wet <- subset(natdata, season == 'wet')
points(wet$SWLslope ~ wet$lang, pch = 17, col = as.factor(wet$climate_class))
legend('bottomright', legend=c(levels(as.factor(dry$climate_class)), 'dry', 'wet'), pch = c(rep(15, 4), 1, 2),
       col =c(1:4, 'black', 'black'), bty = 'n')
rm(wet, dry)

offsetnull <-lmer(mean_offset ~ (1|study/season) , data=modeldata,
                na.action = "na.omit",REML = TRUE)
lcexnull <-lmer(mean_lcexcess ~ (1|study) , data=modeldata,
                  na.action = "na.omit",REML = TRUE)
summary(offsetnull)
summary(lcexnull)
# the overall mean offset is significantly negative :-)
offsetClim <- lmer(mean_offset ~ climate_class + (1|study) , data=natdata,
                      na.action = "na.omit",REML = TRUE)

anova(offsetClim)
visreg(offsetClim, xvar = "climate_class",
       scale = "response", gg=TRUE)
# there are no significant overall differences in offset among climate classes
offsetClimSeason <- lmer(mean_offset ~ climate_class*season + (1|study) ,
                      data=subset(natdata, season != 'not applicable'),
                      na.action = "na.omit",REML = TRUE)
summary(offsetClimSeason)
# no significant differenes neither among climate classes, nor between seasons
offsetLang <- lmer(mean_offset ~ lang + (1|study) , data=natdata,
                      na.action = "na.omit",REML = TRUE)
summary(offsetLang)
# no significant overall effect of the Lang index on the offset

offsetLangSeason <- lmer(mean_offset ~ lang*season + (1|study) ,
                         data=subset(natdata, season != 'not applicable'), na.action = "na.omit",REML = TRUE)
summary(offsetLangSeason)
anova(offsetLangSeason)
# the Lang Index appears to have a negative effect on the magnitude of the offset: the larger the value of 
#  the Lang index (i.e. the more humid), the more negative the offset becomes
windows(8, 8)
par(mfrow=c(1, 1))
dry <- subset(natdata, season == 'dry')
plot(dry$mean_offset ~ dry$lang, pch = 19, col = as.factor(dry$climate_class),
     xlim = c(0, 600), ylim= c(-32, 40), ylab = 'offset (permil)', xlab = 'Lang Index (mm/C)')
wet <- subset(natdata, season == 'wet')
points(wet$mean_offset ~ wet$lang, pch = 17, col = as.factor(wet$climate_class))
legend('bottomright', legend=c(levels(as.factor(dry$climate_class)), 'dry', 'wet'), pch = c(rep(15, 4), 1, 2),
       col =c(1:4, 'black', 'black'), bty = 'n')
rm(wet, dry)


#####offset model (pft)######

offsetnull <- lmer(mean_offset ~ (1|study/season), data=modeldata, na.action = "na.omit",REML = TRUE)
summary(offsetnull)

#plant group
offGroup<-lmer(mean_offset ~plant_group  +
                 (1|study/season), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

summary(offGroup)
anova(offGroup)
myresplot(offGroup, modeldata) 
hist(modeldata$mean_offset) 
visreg(offGroup, xvar = "plant_group", 
       scale = "response", gg=TRUE) #low effect

ggplot(data=subset(modeldata,!plant_group=="NA"), aes(x=plant_group, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() 

#This is the 

#woodiness
offwood<-lmer(mean_offset ~woodiness  +
                (1|study/season), data=modeldata, na.action = "na.omit",REML = TRUE)

summary(offwood)
myresplot(offwood, modeldata) 
hist(modeldata0$mean_offset) 


visreg(offwood, xvar = "woodiness", 
       scale = "response", gg=TRUE) #low effect

ggplot(data=modeldata, aes(x=woodiness, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() #similar results 

#leaf habit
doBy::summaryBy(mean_offset ~ leaf_habit, FUN = lengthWithoutNA, data = modeldata)
# discard semi-deciduous (only 3 observations)
offhabit<-lmer(mean_offset ~leaf_habit  + (1|study/season),
               data= subset(modeldata[which(modeldata$woodiness!='non-woody'),], leaf_habit == 'deciduous' | leaf_habit == 'evergreen'),
               na.action = "na.omit",REML = TRUE)
summary(offhabit)
myresplot(offhabit, modeldata) 
hist(modeldata$mean_offset)
visreg(offhabit, xvar = "leaf_habit", 
       scale = "response", gg=TRUE)
# significant difference: deciduous more negative offset than evergreen

visreg(offhabit, xvar = "leaf_habit", 
       scale = "response", gg=TRUE) 

ggplot(data=modeldata, aes(x=leaf_habit, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() 

