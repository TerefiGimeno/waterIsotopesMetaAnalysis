#We have two databases. One with oultiers (modeldata) and other without
# (modeldata0). We have to repeat the same process in both and see the differences

########modeldata#############
modeldata<-modeldata
#####weights(modeldata)#######

library(lme4)
library(tidyverse)
library(GGally)
library(MuMIn)
library(broom.mixed)
library(visreg)
# load this package to get p-values
library(lmerTest)
library(emmeans)

modeldata$weightSW <-sqrt(((modeldata$se_d2Hplant)^2) + 
                            ((modeldata$SWLslope * modeldata$se_d18Oplant)^2) +
                            ((modeldata$SWLslope.std.error * modeldata$mean_d18Oplant)^2)+
                            ((modeldata$SWLintercept.std.error)^2)
)

modeldata$weightLC <-sqrt(((modeldata$se_d2Hplant)^2) + 
                            ((modeldata$SWLslope * modeldata$se_d18Oplant)^2)
)

# Random-effects weights: recogen la heterogeneidad del estudio en si
Q <- sum(modeldata$mean_offset-weighted.mean(modeldata$mean_offset,w=modeldata$weightSW))^2
#Q: estadastico de heterogeneidad, que toma un valor unico para todo el meta-analisis.
#a mayor Q, mayor heterogeneidad
c <- sum(modeldata$weightSW) - (sum(modeldata$weightSW^2)/sum(modeldata$weightSW))
#C. es como el valor Q pero sin inc?uir la variable dependiente
tDL <- max((Q-(nrow(modeldata)-1))/c,0) #o nos quedamos con el resultado de la ecuaci?n o si no el valor cero
tDL#parametro de heterogeneidad, que toma un valor unico para todo el meta-analisis.
#Es un valor de varianza que solo puede ser positivo.
vi <- 1/modeldata$weightSW #vi=inversa de los pesos
modeldata$weightSWr <- 1/(vi+tDL) # Random-effects weights are defined
i2 <- (Q-(nrow(modeldata)-1))/Q

plot(modeldata$wFE, modeldata$weightSWr)
ggplot(data=modeldata, aes(x=map, y=weightSWr)) +
  geom_boxplot()
ggplot(data=modeldata, aes(x=map, y=weightSW)) +
  geom_boxplot()

#to test which one is better we see wich provides better fits

mf <- lmer(mean_offset ~ (1|study) , data=modeldata,  weights = weightSW, na.action = "na.omit",REML = TRUE)


mr<- lmer(mean_offset ~ (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)

AIC(mf,mr) #mr is better, then we will use forward weights=weightSWr (also for modeldata0)
rm(mr,mf)

######null model(modeldata)##########

null <- lmer(mean_offset ~ (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(null) ###overall existence of the offset (negative effect) WHAT WE WERE LOOKING

#######climate model(modeldata)#########
# this symbol(#) will appear in the ones with significant effects close to summary
mmap <- lmer(mean_offset ~ map + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmap)

mmat <- lmer(mean_offset ~ mat + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmat)#

mtemp <- lmer(mean_offset ~ temp_C + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtemp)###

mtempA <- lmer(mean_offset ~ temp_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtempA)

mtempD <- lmer(mean_offset ~ temp_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtempD)  ###

mSWCa <- lmer(mean_offset ~ int_smwlA + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCa)#?

mSWCb <- lmer(mean_offset ~ int_smwlB + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCb) ##?

mSWCaA <- lmer(mean_offset ~ smIntA_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCaA) 

mSWCbA <- lmer(mean_offset ~ smIntB_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCbA)

mSWCaD <- lmer(mean_offset ~ int_smwlA_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCaA) 

mSWCbD <- lmer(mean_offset ~ int_smwlB_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCbA)

msmwl1 <- lmer(mean_offset ~ smwl1 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl1) ###?

msmwl2 <- lmer(mean_offset ~ smwl2 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl2) ###?

msmwl3 <- lmer(mean_offset ~ smwl3 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl3) ##?

msmwl4 <- lmer(mean_offset ~ smwl4 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl4)

msmwl1A <- lmer(mean_offset ~ sm1_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl1A)

msmwl2A <- lmer(mean_offset ~ sm2_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl2A)

msmwl3A <- lmer(mean_offset ~ sm3_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl3A)

msmwl4A <- lmer(mean_offset ~ sm4_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl4A)

msmwl1D <- lmer(mean_offset ~ sm1_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl1D) ###

msmwl2D <- lmer(mean_offset ~ sm2_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl2D) ###

msmwl3D <- lmer(mean_offset ~ sm3_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl3D) ###

msmwl4D <- lmer(mean_offset ~ sm4_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(msmwl4D) ###

me <- lmer(mean_offset ~ e + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(me) #

meA <- lmer(mean_offset ~ e_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(meA) 

meD <- lmer(mean_offset ~ e_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(meD)

meS <- lmer(mean_offset ~ e_sum + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(meS)

mpev <- lmer(mean_offset ~ pev + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mpev)###

mpevA <- lmer(mean_offset ~ pev_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mpevA) 

mpevD <- lmer(mean_offset ~ pev_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mpevD) #.

mpevS <- lmer(mean_offset ~ pev_sum + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mpevS) 

mtp <- lmer(mean_offset ~ tp + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtp) #

mtpA <- lmer(mean_offset ~ tp_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtpA) #

mtpD <- lmer(mean_offset ~ tp_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtpD) #

mtpS <- lmer(mean_offset ~ tp_sum + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtpS) #


marid <- lmer(mean_offset ~ aridUNEP + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(marid) ###

maridA <- lmer(mean_offset ~ aridUNEP_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(maridA) ###

maridD <- lmer(mean_offset ~ aridUNEP_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(maridD) ###


#climate effects combined


#######biological model(modeldata)########

mpg <- lmer(mean_offset ~ plant_group + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mpg) ###

mwood <- lmer(mean_offset ~ woodiness + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mwood) #

mlh <- lmer(mean_offset ~ leaf_habit + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlh) ###

mls <- lmer(mean_offset ~ leaf_shape + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mls) ###

mwd <- lmer(mean_offset ~ wd + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mwd) #.

mRAP <- lmer(mean_offset ~ RAP + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mRAP) ##

mmyco <- lmer(mean_offset ~ myco + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmyco) #


#biological combined

#######other (modeldata)######

mslt <- lmer(mean_offset ~ slt + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mslt) ###

mmesxyl <- lmer(mean_offset ~ xylem_measurement_method + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmesxyl) ##

mmessoil <- lmer(mean_offset ~ soil_measurement_method + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmessoil) ##

mlailv <- lmer(mean_offset ~ lai_lv + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlailv) 

mlailvA <- lmer(mean_offset ~ laiLV_annual + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlailvA) #?

mlailvD <- lmer(mean_offset ~ laiLV_dif + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlailvD) ###

mlaihv <- lmer(mean_offset ~ lai_hv + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlaihv) #.

mlaihvA <- lmer(mean_offset ~ laiHV_annual + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlaihvA) #

mlaihvD <- lmer(mean_offset ~ laiHV_dif + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlaihvD) ##

#####combined model(modeldata)########



##########SWL#########

smmap <- lmer(SWLslope ~ map + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smmap) ###

smmat <- lmer(SWLslope ~ mat + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smmat) ###

smtemp <- lmer(SWLslope ~ temp_C + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtemp) ###

smtempA <- lmer(SWLslope ~ temp_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtempA) ###

smtempD <- lmer(SWLslope ~ temp_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtempD) ### 

smSWCa <- lmer(SWLslope ~ int_smwlA + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smSWCa) ###

smSWCb <- lmer(SWLslope ~ int_smwlB + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smSWCb) ###

smSWCaA <- lmer(SWLslope ~ smIntA_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smSWCaA) ###

smSWCbA <- lmer(SWLslope ~ smIntB_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smSWCbA) ###

smSWCaD <- lmer(SWLslope ~ int_smwlA_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smSWCaA) ###

smSWCbD <- lmer(SWLslope ~ int_smwlB_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smSWCbA) ###

smsmwl1 <- lmer(SWLslope ~ smwl1 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl1) ###

smsmwl2 <- lmer(SWLslope ~ smwl2 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl2) ###

smsmwl3 <- lmer(SWLslope ~ smwl3 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl3) ###

smsmwl4 <- lmer(SWLslope ~ smwl4 + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl4) ###

smsmwl1A <- lmer(SWLslope ~ sm1_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl1A) ###

smsmwl2A <- lmer(SWLslope ~ sm2_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl2A) ###

smsmwl3A <- lmer(SWLslope ~ sm3_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl3A) ###

smsmwl4A <- lmer(SWLslope ~ sm4_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl4A) ###

smsmwl1D <- lmer(SWLslope ~ sm1_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl1D)  ###

smsmwl2D <- lmer(SWLslope ~ sm2_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl2D)  ###

smsmwl3D <- lmer(SWLslope ~ sm3_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl3D)  ###

smsmwl4D <- lmer(SWLslope ~ sm4_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smsmwl4D) ###

sme <- lmer(SWLslope ~ e + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(sme) ###

meA <- lmer(SWLslope ~ e_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(meA) ###

smeD <- lmer(SWLslope ~ e_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smeD)###

smeS <- lmer(SWLslope ~ e_sum + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smeS)###

smpev <- lmer(SWLslope ~ pev + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smpev)###

smpevA <- lmer(SWLslope ~ pev_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smpevA) ###

smpevD <- lmer(SWLslope ~ pev_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smpevD) ###

smpevS <- lmer(SWLslope ~ pev_sum + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smpevS) ###

smtp <- lmer(SWLslope ~ tp + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtp) ###

smtpA <- lmer(SWLslope ~ tp_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtpA) ###

smtpD <- lmer(SWLslope ~ tp_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtpD) ###

smtpS <- lmer(SWLslope ~ tp_sum + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smtpS) ###


smarid <- lmer(SWLslope ~ aridUNEP + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smarid) ###

smaridA <- lmer(SWLslope ~ aridUNEP_annual + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smaridA) ###

smaridD <- lmer(SWLslope ~ aridUNEP_dif + (1|study) , data=subset(modeldata, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smaridD) ###

smslt <- lmer(SWLslope ~ slt + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smslt) ###

smmessoil <- lmer(SWLslope ~ soil_measurement_method + (1|study) , data=modeldata,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(smmessoil)###

########modeldata0#############

modeldata0 <-subset(modeldata,mean_offset<20)
modeldata0 <-subset(modeldata0,mean_offset>-25)
modeldata0 <-subset(modeldata0,mean_lcexcess<20)
modeldata0 <-subset(modeldata0,mean_lcexcess>-25)

#####weights(modeldata0)#######

modeldata0$weightSW <-sqrt(((modeldata0$se_d2Hplant)^2) + 
                            ((modeldata0$SWLslope * modeldata0$se_d18Oplant)^2) +
                            ((modeldata0$SWLslope.std.error * modeldata0$mean_d18Oplant)^2)+
                            ((modeldata0$SWLintercept.std.error)^2)
)

modeldata0$weightLC <-sqrt(((modeldata0$se_d2Hplant)^2) + 
                            ((modeldata0$SWLslope * modeldata0$se_d18Oplant)^2)
)


Q <- sum(modeldata0$mean_offset-weighted.mean(modeldata0$mean_offset,w=modeldata0$weightSW))^2
c <- sum(modeldata0$weightSW) - (sum(modeldata0$weightSW^2)/sum(modeldata0$weightSW))
tDL <- max((Q-(nrow(modeldata0)-1))/c,0) 
vi <- 1/modeldata0$weightSW 
modeldata0$weightSWr <- 1/(vi+tDL) 

######null model(modeldata0)##########

null0 <- lmer(mean_offset ~ (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(null0) #the effect is the same but lower

#######climate model(modeldata)#########
# this symbol(#) will appear in the ones with significant effects close to summary
mmap0 <- lmer(mean_offset ~ map + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmap0)

mmat0 <- lmer(mean_offset ~ mat + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmat0)#no longer effect

mtemp0 <- lmer(mean_offset ~ temp_C + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtemp0)### stay the same

mtempA0 <- lmer(mean_offset ~ temp_annual + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mtempA0)

mSWCa0 <- lmer(mean_offset ~ int_smwlA + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCa0)#. less effect

mSWCb0 <- lmer(mean_offset ~ int_smwlB + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCb0) # less effect

mSWCaA0 <- lmer(mean_offset ~ smIntA_annual + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCaA0) 

mSWCbA0 <- lmer(mean_offset ~ smIntB_annual + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mSWCbA0)

mmer0 <- lmer(mean_offset ~ mer + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmer0) #

mmerA0 <- lmer(mean_offset ~ mer_annual + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmerA0) 

mmper0 <- lmer(mean_offset ~ mper + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmper0)###

mmperA0 <- lmer(mean_offset ~ pet_annual + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmperA0) 

mmtpr0 <- lmer(mean_offset ~ mtpr + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmap0)

mmtprA0 <- lmer(mean_offset ~ p_annual + (1|study) , data=subset(modeldata0, natural == 'natural'),  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmap0)

#climate effects combined



#######biological model(modeldata)########

mpg0 <- lmer(mean_offset ~ plant_group + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mpg0) # much less effect

mwood <- lmer(mean_offset ~ woodiness + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mwood) ### more effect

mlh0 <- lmer(mean_offset ~ leaf_habit + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlh0) ###

mls0 <- lmer(mean_offset ~ leaf_shape + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mls0) ## less effect

mwd0 <- lmer(mean_offset ~ wd + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mwd0) # looses effect

mRAP0 <- lmer(mean_offset ~ RAP + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mRAP0) # looses effect

mmyco0 <- lmer(mean_offset ~ myco + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmyco0) # looses effect

#biological combined

#######other (modeldata)######

mslt0 <- lmer(mean_offset ~ slt + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mslt0) #

mmesxyl0 <- lmer(mean_offset ~ xylem_measurement_method + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmesxyl0) # looses effect

mmessoil0 <- lmer(mean_offset ~ soil_measurement_method + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mmessoil0) ##

mlailv0 <- lmer(mean_offset ~ lai_lv + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlailv0) 

mlailvA0 <- lmer(mean_offset ~ laiLV_annual + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlailvA0) #?

mlaihv0 <- lmer(mean_offset ~ lai_hv + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlaihv0) #.

mlaihvA0 <- lmer(mean_offset ~ laiHV_annual + (1|study) , data=modeldata0,  weights = weightSWr, na.action = "na.omit",REML = TRUE)
summary(mlaihvA0) # looses effect

#####combined model(modeldata)########