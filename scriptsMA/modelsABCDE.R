source('scriptsMA/generate_modeldata.R')
n_distinct(modeldata$species_plant)
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

natdata <- subset(modeldata, natural== 'natural') #reject hydro-managed plots, modeldata for pft models

########SWL######

swlbest1 <-lmer(SWLslope ~ 
                  (1|study/season), data=natdata,
                na.action = "na.omit",REML = FALSE)

swlbest1.1<- lmer(SWLslope ~  (1|study), data = subset(natdata, season != 'not applicable'),
               na.action = 'na.omit', REML = FALSE)

swlbest1.2<- lmer(SWLslope ~  (1|study/season), data = subset(natdata, season != 'not applicable'),
                na.action = 'na.omit', REML = FALSE)

swlbest2 <-lmer(SWLslope ~ climate_class +
                  (1|study/season), data=natdata,
                na.action = "na.omit",REML = FALSE)

swlbest3 <-lmer(SWLslope ~ map*mat +
                  (1|study/season), data=natdata,
                na.action = "na.omit",REML = FALSE)

swlbest4 <-lmer(SWLslope ~ map +
                (1|study/season), data=natdata,
                na.action = "na.omit",REML = FALSE)

swlbest4.1 <-lmer(SWLslope ~ map +
                  (1|study/season), subset(natdata, season != 'not applicable'),
                  na.action = 'na.omit', REML = FALSE)

swlbest4.2 <-lmer(SWLslope ~ map +
                  (1|study/season), data=natdata,
                na.action = "na.omit",REML = FALSE)

swlbest5 <-lmer(SWLslope ~ mat +
                (1|study), data=natdata,
                na.action = "na.omit",REML = FALSE)

swlbest6 <-lmer(SWLslope ~ lang +
                (1|study/season), data=natdata,
                na.action = "na.omit",REML = FALSE)


swlbest0 <-lmer(SWLslope ~ climate_class*season +
                 (1|study) , data=subset(natdata, season != 'not applicable'),
               na.action = "na.omit",REML = FALSE) 

swlbest00 <-lmer(SWLslope ~ map*season +
                  (1|study) , data=subset(natdata, season != 'not applicable'),
                na.action = "na.omit",REML = FALSE) 
 
AIC(swlbest1,swlbest1.1,swlbest1.2,swlbest2,swlbest3,swlbest4,swlbest4.1,swlbest4.2,swlbest5,swlbest6,swlbest0,swlbest00)

summary(swlbest4.1)
summary(swlbest00)

emmeans(swlbest00, list(pairwise ~ season), adjust = "tukey")
#we dont want not applicable values
# natdata[which(natdata$season == "not applicable"), 'season'] <- NA

#general model
swlbest11 <-lmer(SWLslope ~ climate_class +
                  (1|study), data=natdata,
                na.action = "na.omit",REML = TRUE)
summary(swlbest1)
anova(swlbest1)
emmeans(swlbest1, list(pairwise ~ climate_class), adjust = "tukey")
visreg(swlbest1, xvar = "climate_class",
       scale = "response", gg=TRUE)
check_model(swlbest1)

# no significant differences in swl among climate classes

swlbest2 <-lmer(SWLslope ~ lang +
                  (1|study) , data=natdata,
                na.action = "na.omit",REML = TRUE)
myresplot(swlbest2, natdata)
summary(swlbest2)
visreg(swlbest2, xvar = "lang",
       scale = "response", gg=TRUE)
plot(natdata$SWLslope ~ natdata$lang, pch = 19, col = as.factor(natdata$climate_class),
     ylab = 'log (slope SWL)', xlab = 'Lang Index (mm/c)')
legend('bottomright', legend = levels(as.factor(natdata$climate_class)), pch = 19, col = 1:4, bty = 'n')
plot(natdata$SWLslope ~ natdata$lang, pch = 19, col = as.factor(natdata$climate_class),
     ylab = 'log (slope SWL)', xlab = 'Lang Index (mm/c)', xlim = c(0, 450))

# important disequilibrium in sample size across climate_classes
#doBy::summaryBy(SWLslope ~ season + climate_class, FUN = lengthWithoutNA, data = subset(natdata, season != 'not applicable'))
#doBy::summaryBy(SWLslope ~ climate_class, FUN = lengthWithoutNA, data = natdata)
summarise(group_by(modeldata, climate_class), count = lengthWithoutNA(mean_offset))
summarise(group_by(natdata, climate_class), count = lengthWithoutNA(mean_offset))

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

swlLang <- lmer(SWLslope ~ lang*season + (1|study)-1, data = subset(natdata, season != 'not applicable'),
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

swlnull<- lmer(SWLslope ~  (1|study), data = subset(natdata, season != 'not applicable'),
               na.action = 'na.omit', REML = TRUE)

summary(swlnull)
anova(swlnull)
emmeans(swlnull, list(pairwise ~ climate_class), adjust = "tukey")
visreg(swlnull, xvar = "climate_class",
       scale = "response", gg=TRUE)
check_model(swlnull)
AIC(swlnull,swlbest,swlbest1,swlbest2,swlLang)


####map##########

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(rgeos)
theme_set(theme_bw())

#world map load
world <- ne_countries(scale = "medium", returnclass = "sf")

#dataframe of coords and climate:
sites <- data.frame(modeldata$log, modeldata$lat)
names(sites)[1] <- "longitude"
names(sites)[2] <- "latitude"
sites<- unique(sites)
sites<- subset(sites, !longitude== "NA")

#same projections (for WGS84 CRS code is #4326)
sites<- st_as_sf(sites, coords = c("longitude", "latitude"), 
                 crs = 4326, agr = "constant")

#plot:
map<-ggplot(data = world) +
  geom_sf() +
  geom_sf(data = sites, size = 3.7, shape = 21, 
          fill = "dodgerblue2", color="dodgerblue3", alpha=0.9) +
  coord_sf(xlim = c(-170, 170), ylim = c(-90, 90), expand = FALSE)
x11()
map

ggsave("map.jpg", map, width = 20, height = 10, dpi = 900)

rm(map, sites, world)




summarise(group_by(modeldata, season), count = lengthWithoutNA(season))
n_distinct(modeldata$species_plant)
arid <- natdata %>% 
  filter(climate_class=="arid")

n_distinct(arid$campaign)
n_distinct(arid$study)
n_distinct(arid$species_plant)
summarise(group_by(arid, season), count = lengthWithoutNA(season))
mean(arid$SWLslope)
max(arid$SWLslope)
min(arid$SWLslope)
mean(arid$mean_lcexcess)
max(arid$mean_lcexcess)
min(arid$mean_lcexcess)
mean(arid$mean_offset)
max(arid$mean_offset)
min(arid$lcexcess)


tropical <- natdata %>% 
  filter(climate_class=="tropical")

n_distinct(tropical$campaign)
n_distinct(tropical$study)
n_distinct(tropical$species_plant)
summarise(group_by(tropical, season), count = lengthWithoutNA(season))
mean(tropical$SWLslope)
max(tropical$SWLslope)
min(tropical$SWLslope)
mean(tropical$mean_lcexcess)
max(tropical$mean_lcexcess)
min(tropical$mean_lcexcess)
mean(tropical$mean_offset)
max(tropical$mean_offset)
min(tropical$lcexcess)

warm <- natdata %>% 
  filter(climate_class=="warm")

n_distinct(warm$campaign)
n_distinct(warm$study)
n_distinct(warm$species_plant)
summarise(group_by(warm, season), count = lengthWithoutNA(season))
mean(warm$SWLslope)
max(warm$SWLslope)
min(warm$SWLslope)
mean(warm$mean_lcexcess)
max(warm$mean_lcexcess)
min(warm$mean_lcexcess)
mean(warm$mean_offset)
max(warm$mean_offset)
min(warm$lcexcess)



cold <- modeldata %>% 
  filter(climate_class=="cold")

n_distinct(cold$campaign)
n_distinct(cold$study)
n_distinct(cold$species_plant)
summarise(group_by(cold, season), count = lengthWithoutNA(season))
mean(cold$SWLslope)
max(cold$SWLslope)
min(cold$SWLslope)
mean(cold$mean_lcexcess)
max(cold$mean_lcexcess)
min(cold$mean_lcexcess)
mean(cold$mean_offset)
max(cold$mean_offset)
min(cold$lcexcess)


angio<- modeldata %>% filter(plant_group=="angiosperm")
gymno<- modeldata %>% filter(plant_group=="gymnosperm")

n_distinct(angio$species_plant)
n_distinct(gymno$species_plant)
n_distinct(natdata$campaign, natdata$season)

dry<- natdata %>% filter(season=="dry")
n_distinct(dry$campaign)
wet<- natdata %>% filter(season=="wet")
n_distinct(wet$campaign)
not<- natdata %>% filter(season=="not applicable")          
n_distinct(not$campaign)
