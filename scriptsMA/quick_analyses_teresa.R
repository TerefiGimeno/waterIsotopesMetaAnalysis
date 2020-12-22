library(lme4)
library(tidyverse)
library(GGally)
library(MuMIn)
library(broom.mixed)
library(visreg)
# load this package to get p-values
library(lmerTest)
library(emmeans)

null <- lmer(mean_offset ~ (1|study) , data=modeldata, na.action = "na.omit",REML = TRUE)
# the overall offset is still significantly differnt from zero

myco <- lmer(mean_offset ~ myco + (1|study) , data=modeldata, na.action = "na.omit",REML = TRUE)
anova(myco)
# no significant differences in the offset among grous with/without mycos

wd <- myco <- lmer(mean_offset ~ wd + (1|study) , data=modeldata, na.action = "na.omit",REML = TRUE)
# no significant correlation with wood density

# analyses with monty data (only natural sites)
natural <- subset(modeldata, natural == 'natural')
temp <- lmer(mean_offset ~ temp_C + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(temp)
# singificant positive correlation between mean offset and temperature (the warmer the closer the offset is to zero)
sm1 <-lmer(mean_offset ~ smwl1 + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm1)
sm2 <-lmer(mean_offset ~ smwl2 + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm2)
sm3 <-lmer(mean_offset ~ smwl3 + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm3)
# significant negative correlation between the mean offset and soil moisture at 7, 28 and 100 cm depth
sm4 <-lmer(mean_offset ~ smwl4 + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm4)
# no significant correlation between mean offset and soil moisture at 200 cm
smA <-lmer(mean_offset ~ int_smwlA + (1|study), data=natural, na.action = "na.omit",REML = TRUE)
summary(smA)
smB <-lmer(mean_offset ~ int_smwlB + (1|study), data=natural, na.action = "na.omit",REML = TRUE)
summary(smB)
# significant negative correlation between the mean offset and soil content in the first 100 and 200 cm
smTemp <-lmer(mean_offset ~ temp_C * smwl1 + (1|study), data=natural, na.action = "na.omit",REML = TRUE)
summary(smTemp)
# the model with both temperature and smwl1 shows that the correlation is significant with temprature
# but not with smwl1, nor the interaction between them
# other varaibles
slt <- lmer(mean_offset ~ slt + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
# no significant effects of soil type
# there is something off with variables mper y mtpr
laiHv <- lmer(mean_offset ~ lai_hv + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
laiLv <- lmer(mean_offset ~ lai_lv + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
# no significant correlation with LAI

# analyses with annual trends
temp <- lmer(mean_offset ~ temp_annual + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(temp)
# no significant correlation with mean annual temperature
sm1 <- lmer(mean_offset ~ sm1_annual + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm1)
sm2 <- lmer(mean_offset ~ sm2_annual + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm2)
sm3 <- lmer(mean_offset ~ sm3_annual + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm3)
sm4 <- lmer(mean_offset ~ sm4_annual + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(sm4)
# significant negative correlation with sm2 and marginally significant with sm1
smA <- lmer(mean_offset ~ int_smwlA + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(smA)
smB <- lmer(mean_offset ~ int_smwlB + (1|study) , data=natural, na.action = "na.omit",REML = TRUE)
summary(smB)
# significant negative correatlin with both integrated soil moistures

group <- lmer(mean_offset ~ plant_group + (1|study) , data=modeldata, na.action = "na.omit",REML = TRUE)
summary(group)
group <- lmer(mean_offset ~ growth_form + (1|study) , data=modeldata, na.action = "na.omit",REML = TRUE)
anova(group)
group <- lmer(mean_offset ~ leaf_habit + (1|study) , data=modeldata, na.action = "na.omit",REML = TRUE)
