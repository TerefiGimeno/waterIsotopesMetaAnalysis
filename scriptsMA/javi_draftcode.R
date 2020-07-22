library(tidyverse)
metadata<- read.csv("metaprueba.csv", header= T,sep=";")
plantdata <- read.csv("plantaprueba.csv", header= T,sep=";")
sourcedata<- read.csv("fuenteprueba.csv", header=T, sep=";")

########## Obtainig the offset values
# Mean values of plant data

aggregate(cbind(plantdata$d18O_permil_plant+plantdata$d2H_permil_plant)~plantdata$id.campaign+plantdata$species,plantdata, mean)

# we create a dataframe for the analysis

PlantMeans<-aggregate(cbind(d18O_permil_plant,d2H_permil_plant)~id.campaign+species,plantdata, mean)

# we need the n value 
idcamp.sp<- plantdata %>%
  group_by(id.campaign,species) %>%
  summarise(id.campaign_species_count=n())

#Now we want to join plantmeans and idcamp.sp
joinedplant<- merge(idcamp.sp,PlantMeans, all.x = F, all.y=T)

#zSlope and intercept of the SWL

#We filter the data. We are only interested in soil data from label column 

soils<-sourcedata[sourcedata$label_class=="soil",]


lm(soils$deute~soils$Oxi, data=soils, subset = soils$id.campaign) #no me funciona

}
# First the values of slope

# Now the intercept

# Now we create a dataframe and joint it with PlantMeans

# We use the offset formula to