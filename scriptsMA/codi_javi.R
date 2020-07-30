setwd('C:/Users/Adrià/Documents/Teresa Gimeno') ##modificar!
library(tidyverse)
library(ggpubr)
library(broom)
library(tidyr)

###I removed some weird numbers (not even numbers sometimes)
source<-read.csv2('data_water_sources_adrià.csv')

source$d2H_permil_source<-as.numeric(as.character(source$d2H_permil_source))   #for some reason those are factors...
source$d18O_permil_source<-as.numeric(as.character(source$d18O_permil_source))
source$year<-as.factor(source$year)

hist(source$d18O_permil_source)

str(source)  ###check the type of variable

source<-subset(source,!d2H_permil_source=="NA")
source<-subset(source,!d18O_permil_source=="NA") ##clean NA

source<-subset(source,label_pool=="bulk") ##only bulk soil

multiple<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-author,-date,-plot) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source, na.action='na.omit',data = .x)),
    tidied = map(fit, tidy)
  ) %>% 
  unnest(tidied)

multiple<-multiple[,c('author','date','plot','term','estimate','std.error',"statistic",'p.value')]  ###select only relevant columns

###more parameters from the regression, including r squared
rsquared<-subset(source,label_class=='soil') %>%    ###this function computes the linear regressions
  nest(-author,-date,-plot) %>%                       ###here you can choose which factors should be used
  mutate(
    fit = map(data, ~ lm(d2H_permil_source ~ d18O_permil_source,data = .x)),
    tidied = map(fit, glance)
  ) %>% 
  unnest(tidied)

rsquared<-rsquared[,c('author','date','plot','r.squared')]  ###select only relevant columns

###count the number of soil samples, in case we want to cut off using this
length<- source %>% count(author,date,plot)

###divide into intercept and slope
intercept<-subset(multiple,term=='(Intercept)')
slope<-subset(multiple,term=='d18O_permil_source')

hist(slope$estimate)

weird<-subset(slope,estimate<0)  ##the slope should be positive...anyway those got cut when using the p.value, except Twining

swl<-merge(intercept,slope,by=c('author','date','plot'))   ###create a new table with separate columns for intercept and slope

###give proper names
colnames(swl)<-c("author",'date','plot',"term","estimate","std.error","statistic","p.value","term.slope","estimate.slope","std.error.slope",
                 "statistic.slope","p.value.slope")
##paste rsquareds
swl<-merge(swl,rsquared,by=c('author','date','plot'))

#past N's
swl<-merge(swl,length,by=c('author','date','plot'))

swl<-subset(swl,p.value.slope<0.05&n>2)   #cutoff non-significant and poorly fit regressions

swl<-subset(swl,p.value.slope<0.05&n>2&!date=="99999") ###remove undated data?

plant<-read.csv2('data_plant_water_adria.csv')

plant$d2H_permil_plant<-as.numeric(as.character(plant$d2H_permil_plant))   #for some reason those are factors...
plant$d18O_permil_plant<-as.numeric(as.character(plant$d18O_permil_plant))

hist(plant$d2H_permil_plant)
hist(plant$d18O_permil_plant)

plant<-subset(plant,!plant_tissue=='leaf')  #we don't want leaves

plant<-merge(plant,swl,by=c('author','date','plot'))   

###let's calculate the offset!

plant$offset<-plant$d2H_permil_plant-plant$estimate.slope*plant$d18O_permil_plant-plant$estimate

hist(plant$offset) ###it looks OK

###in case they are not charged yet
#library(lme4)
#library(lmerTest)

###for the metadata file, I normalised the 'precipitation season' and the method of analysis
meta<-read.csv2('meta_data_adria.csv')

plant<-merge(plant,meta,by=c('author','species'))  ##attach metadata, I think plot not specified

###it seems that some papers have stem evaporation, if this stated in the paper, we should remove them from the analysis!

###here you can start testing the effect of any factor, some factors need revision though! 
###They are not consistent or have too many levels

plant$elevation<-as.numeric(as.character(plant$elevation))   #most of variables appear as factor, need to transform to numeric
plant$Bar_GW_use<-as.numeric(as.character(plant$Bar_GW_use)) 
plant$mat<-as.numeric(as.character(plant$mat)) 
plant$map<-as.numeric(as.character(plant$map)) 

##create a more simple factor for the analysis type
plant$laser_mass[plant$method_simple=='CRDS']<-'Laser'
plant$laser_mass[plant$method_simple=='IRIS']<-'Laser'
plant$laser_mass[plant$method_simple=='IRMS']<-'MassSpec'

###I think that author should be a random factor, perhaps also date and species. However species is not 
###nested with author because the same species may appear in different studies

summary(model<-lmer(offset~(1|author:date)+(1|species),data=plant)) 
plot(allEffects(model))

summary(model<-lmer(offset~Bar_Precip_Season+(1|author:date)+(1|species),
                    data=subset(plant,!(Bar_Precip_Season=="")))) 

summary(model<-lmer(offset~season+(1|author:date)+(1|species),
                    data=subset(plant,!(season=="not applicable")))) 

summary(model<-lmer(offset~KG+(1|author:date)+(1|species),
                    data=subset(plant,!(KG=="")))) 

summary(model<-lmer(offset~laser_mass+(1|author:date),
                    data=subset(plant,!(laser_mass=="")))) 

summary(model<-lmer(offset~leaf_habit+(1|author:date),
                    data=subset(plant,!(leaf_habit=="")))) 

###example of full model
plant_full<-subset(plant,!(mat=="99999"|map=="99999"|elevation=='99999'|laser_mass==""|growth_form==""|leaf_habit==""|Bar_Precip_Season==""|KG=="")) 

summary(model<-lmer(offset~leaf_habit+growth_form+(1|author:date),
   data=plant_full))

library(sjstats)

##calculates standarised coeffcients
std_beta(model)

###visualise effects
#library(effects)
plot(allEffects(model))

ggplot(data=subset(plant,!laser_mass==''),aes(x=laser_mass,y=offset))+
  geom_boxplot()+
  stat_compare_means()+
  theme_bw()


###means offset and relationship with Morris' and Zanne's datasets

means_offset<-plant %>%
  group_by(author,date,plot,laser_mass)%>%
  summarise(mean_offset=mean(offset,na.rm=T),se_offset=sd(offset)/sqrt(length(offset)))

means_offset$single<-paste(means_offset$author,means_offset$plot,means_offset$date)

ggplot(data=subset(means_offset,!laser_mass==""))+
  geom_bar(aes(x=single,y=mean_offset,fill=laser_mass),stat="identity")+
  theme_bw()+
  theme(panel.grid=element_blank())

means_offset<-merge(means_offset,meta,by=c('author','species'))

ggplot(data=subset(means_offset,!mat=='99999'),aes(x=as.numeric(as.character(mat)),y=mean_offset))+
  geom_point()+
  geom_smooth(method=lm)+
  stat_cor()

ggplot(data=subset(means_offset,!mat=='99999'),aes(x=laser_mass,y=mean_offset))+
  geom_boxplot()+
  stat_compare_means()

wood_density<-read.csv2('wood_density.csv')

colnames(wood_density)

wood_density$wd<-as.numeric(as.character(wood_density$wd))

means_wd<-wood_density %>%
  group_by(species)%>%
  summarise(mean_wd=mean(wd,na.rm=T))

plant2<-merge(plant,means_wd,by='species')

ggplot(plant2,aes(x=mean_wd,y=offset))+
  geom_point()+
  stat_cor()+
  geom_smooth(method=loess)

summary(model<-lmer(offset~mean_wd+(1|author:date),plant2))

plant3<-plant2 %>%
  group_by(species,date)%>%
  summarise(mean_wd=mean(mean_wd),mean_offset=mean(offset,na.rm=T),se_offset=sd(offset,na.rm=T/sqrt(length(offset))))

ggplot(data=subset(plant3,!se_offset=="NA"),aes(x=mean_wd,y=mean_offset))+
  geom_point()+
  stat_cor()+
  geom_smooth(method=lm)+
  geom_errorbar(aes(x=mean_wd,ymin=mean_offset-se_offset,ymax=mean_offset+se_offset))+
  theme_bw()

rap<-read.csv2('rap.csv')

colnames(rap)

rap$RAP<-as.numeric(as.character(rap$RAP))
rap$AP<-as.numeric(as.character(rap$AP))
rap$RP<-as.numeric(as.character(rap$RP))

means_rap<-rap %>%
  group_by(species)%>%
  summarise(mean_rap=mean(RP,na.rm=T))

plant2<-merge(plant,means_rap,by='species')

summary(model<-lmer(offset~mean_rap+laser_mass+(1|author),plant2))

ggplot(plant2,aes(x=mean_rap,y=offset))+
  geom_point()+
  stat_cor()+
  geom_smooth(method=lm)

###altogether

plant4<-merge(plant,means_wd,by='species')
plant5<-merge(plant4,means_rap,by='species')

summary(model<-lmer(offset~mean_rap+mean_wd+(1|author),plant5))

library(effects)

plot(allEffects(model))

ggplot(plant5,aes(x=mean_rap,y=offset,col=mean_wd))+
  geom_point(size=3)+
  stat_cor()+
  geom_smooth(method=lm)+
  theme_bw()

means2<-merge(means_offset,means_rap,by='species')
means3<-merge(means2,means_wd,by='species')

ggplot(data=subset(means3,!mean_offset<(-50)),aes(x=mean_wd,y=mean_offset,col=mean_rap))+
  geom_point(size=3)+
  stat_cor()+
  geom_smooth(method=lm)+
  theme_bw()
