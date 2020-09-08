write.csv(modeldata, "C:\\Users\\JAVI\\Desktop\\pRojects\\waterIsotopesMetaAnalysis\\dataMA\\modeldata.csv")  
source('modeldata/dataMA')
rm(swl0, rareoffset)

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
library(visreg)libray(lmerTest)

#Lets begin with SWL slope 

natdata<- subset(modeldata, natural== 'natural') #reject hydro-managed plots

#we dont want not applicable values
natdata[which(natdata$season == "not applicable"), 'season'] <- NA

#general model

swlbest <-lmer(SWLslope ~ climate_class*season +
                 (1|study) , data=natdata, na.action = "na.omit",REML = TRUE) #has to be true

summary(swlbest)
myresplot(swlbest, natdata) #looks good
hist(natdata$SWLslope) 

visreg(swlbest, xvar = "season",
       scale = "response", gg=TRUE) 
visreg(swlbest, xvar = "climate_class",
       scale = "response", gg=TRUE) #both makes sense

visreg(swlbest, xvar = "climate_class", by="season",
       scale = "response", gg=TRUE)

ggplot(data=subset(natdata,!season=="NA"), aes(x=climate_class, y=SWLslope)) +
  geom_point(alpha=0.2, aes()) +
  geom_boxplot()+
  facet_grid(.~season)  #i dont know how to unplot the NA's


##We are going to see what happends inside each climate class for the SWL

natdata_arid <- natdata %>% 
  filter(climate_class=="arid")

natdata_cold <- natdata %>% 
  filter(climate_class=="cold")

natdata_trop <- natdata %>% 
  filter(climate_class=="tropical")

natdata_warm <- natdata %>% 
  filter(climate_class=="warm")

######SWL_warm########


swlbest_warm <- lmer(SWLslope ~ mat:season + map:season + season +
                       (1|study), data=natdata_warm, na.action = "na.omit", REML = TRUE)

#####este no lo he modificado de : a *, tratadlo como veais

summary(swlbest_warm)
anova(swlbest_warm)
myresplot(swlbest_warm, natdata_warm0)
hist(natdata_warm$SWLslope)


visreg(swlbest_warm, xvar = "map", by = "season", 
       scale = "response", gg=TRUE)
visreg(swlbest_warm, xvar = "mat", by = "season", 
       scale = "response", gg=TRUE)

#interesting. SWL slope increases with map (gets closer to lmwl) in dry and wet plots. 
#but its more drastic in dry season plots. For mat, it increases in dry seasons but decreases on wet ones
#(I thought that it'll also decrease un dry seasons)

ggplot(data=subset(natdata_warm,!season=="NA"), aes(x=map, y=SWLslope, colour=season)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F) +
  facet_grid(season ~ .)

ggplot(data=subset(natdata_warm,!season=="NA"), aes(x=mat, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .)
#not the same results

######SWL_arid#######

swlbest_arid<- lmer(SWLslope ~  map*season + 
                      (1|study), data=natdata_arid, na.action = "na.omit", REML = TRUE)


summary(swlbest_arid)
anova(swlbest_arid)
myresplot(swlbest_arid, natdata_arid)
hist(natdata_warm$SWLslope)


visreg(swlbest_warm, xvar = "season", 
       scale = "response", gg=TRUE) #slopes are stepper in wet seasons (makes sense)

ggplot(data=subset(natdata_arid,!season=="NA"), aes(x=season, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_boxplot() #different results as before

######SWL_tropical#######

swlbest_tropical<- lmer(SWLslope ~ mat*season  + 
                          (1|study), data=natdata_trop, na.action = "na.omit", REML = TRUE)


summary(swlbest_tropical)
anova(swlbest_tropical)
myresplot(swlbest_tropical, natdata_trop) #strange...
hist(natdata_trop$SWLslope) #...


visreg(swlbest_warm, xvar = "mat", by= "season",
       scale = "response", gg=TRUE) #similar as found in SWLwarm model

ggplot(data=subset(natdata_trop,!season=="NA"), aes(x=mat, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .) #opposite as the visreg

######SWL_cold#######

swlbest_cold<- lmer(SWLslope ~ lang*season +
                      (1|study), data=natdata_cold, na.action = "na.omit",REML = TRUE)


summary(swlbest_cold)
anova(swlbest_cold)
myresplot(swlbest_cold, natdata_cold) 
hist(natdata_cold$SWLslope) #...


visreg(swlbest_cold, xvar = "lang", by= "season",
       scale = "response", gg=TRUE) #makes sense with the warm and tropical findings


ggplot(data=subset(natdata_cold,!season=="NA"), aes(x=lang, y=SWLslope)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .) #similar as the visreg

#####offset model (pft)######

#plant group
offGroup<-lmer(mean_offset ~plant_group  +
                 (1|study), data=subset(modeldata, !woodiness=="non-woody"), na.action = "na.omit",REML = TRUE)

summary(offGroup)
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
                (1|study), data=modeldata, na.action = "na.omit",REML = TRUE)

summary(offwood)
myresplot(offwood, modeldata) 
hist(modeldata0$mean_offset) 


visreg(offwood, xvar = "woodiness", 
       scale = "response", gg=TRUE) #low effect

ggplot(data=modeldata, aes(x=woodiness, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() #similar results 

#leaf habit
offhabit<-lmer(mean_offset ~leaf_habit  +
                 (1|study), data=modeldata, na.action = "na.omit",REML = TRUE)

summary(offhabit)
myresplot(offhabit, modeldata) 
hist(modeldata$mean_offset) 


visreg(offhabit, xvar = "leaf_habit", 
       scale = "response", gg=TRUE) 

ggplot(data=modeldata, aes(x=leaf_habit, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() 

#leaf shape
offshape<-lmer(mean_offset ~leaf_shape  +
                 (1|study), data=modeldata, na.action = "na.omit",REML = TRUE)

summary(offshape)
myresplot(offhabit, modeldata) 
hist(modeldata$mean_offset) 


visreg(offshape, xvar = "leaf_shape", 
       scale = "response", gg=TRUE) 


ggplot(data=modeldata, aes(x=leaf_shape , y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_boxplot() 

AIC(offGroup,offhabit,offshape,offwood) #offgroup is the best but shape also

#####offset model (climate)######

offbest_climate<- lmer(mean_offset ~ lang*season + 
                         (1|study), data=natdata, na.action = "na.omit",REML = TRUE)


summary(offbest_climate)
anova(offbest_climate)
myresplot(offbest_climate, natdata) 
hist(natdata$mean_offset) 


visreg(offbest_climate, xvar = "lang", by= "season",
       scale = "response", gg=TRUE) #looks similar for dry and wet seasons. less humidity, more offset (makes sense)

ggplot(data=subset(natdata,!season=="NA"),aes(x=lang, y=mean_offset)) +
  geom_point(alpha=0.2) +
  geom_smooth(method=lm,se=F)  +
  facet_grid(season ~ .) #similar as the visreg

##### offset's general model########

offgeneral<- lmer(mean_offset ~ lang*season + plant_group +
                    (1|study), data=natdata, na.action = "na.omit",REML = TRUE)

summary(offgeneral)
anova(offgeneral)
myresplot(offgeneral, natdata) 
hist(natdata$mean_offset) 


visreg(offgeneral, xvar = "lang", by= "season",
       scale = "response", gg=TRUE) #looks similar for dry and wet seasons. less humidity, more offset (makes sense)
visreg(offgeneral, xvar = "lang", by= "plant_group",
       scale = "response", gg=TRUE) #same patterns
visreg(offgeneral, xvar = "lang", 
       scale = "response", gg=TRUE) #same pattern

