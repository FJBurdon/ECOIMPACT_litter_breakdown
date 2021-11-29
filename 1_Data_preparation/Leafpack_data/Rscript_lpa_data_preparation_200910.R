############################################################################################
##                                                                                        ##  
##                                LITTER BREAKDOWN                                        ##  
##                                                                                        ##
############################################################################################

##  Script: ECOIMPACT leafpack assay data preparation
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Prepare data for analysis by calculating breakdown coefficients
##  2. Correct or Remove erroneous values (i.e., human error, torn leafbags, etc.)
##  3. Calculate the invertebrate contribution to litter processing

##***************************************************************************************
##Load libraries you need
library(reshape)
library(reshape2)
library(plyr)
library(tidyverse)
library(readr)

##***************************************************************************************
## Load data - all litter types
leafpack<-read_csv("1_Input_data/DATA_lpa_all_200731.csv")
leafpack<- as.data.frame(leafpack)
str(leafpack)

## Load temperature data Marta calculated from ECOMIMPACT loggers
Temperature_data_Marta <- read_csv("1_Input_data/Temperature_data_Marta.csv")
Temperature_data_Marta <- as.data.frame(Temperature_data_Marta)

##***************************************************************************************
## Calculate temperature-related variables

## Check data
head(leafpack)

## Rationalise data
lpa_data <- leafpack[,c(1,5,4,6,7,19,8,14,13,10,11)]

## Calculate Final AFDM
lpa_data$Final_AFDM <- lpa_data$Dry-lpa_data$Ashed

## Merge Marta´s temperature data with leafpack data
lpa_data <- merge(lpa_data,Temperature_data_Marta[,c(1,3,6)],by=c("Site_code","Location"),all.x = T)

## Calculate mean temperature
lpa_data$T_mu <- lpa_data$DD/lpa_data$Days

## Rationalize and reorganise
lpa_data <- lpa_data[,c(3,1,2,4:8,14,13,9,12)]
head(lpa_data)

#****************************************************************************************
## Calculate mass loss (k)

lpa_data$k_ml   <- -1*log(lpa_data$Final_AFDM/lpa_data$Initial_AFDM)

##***************************************************************************************
## Calculate mass loss (k days)

lpa_data$k_day_ml   <- -1*log(lpa_data$Final_AFDM/lpa_data$Initial_AFDM)/lpa_data$Days 

##***************************************************************************************
## Calculate mass loss (k degree-days)

lpa_data$k_dd_ml   <- -1*log(lpa_data$Final_AFDM/lpa_data$Initial_AFDM)/lpa_data$DD 

##***************************************************************************************
## Separate Oak

leafpack_oak<-lpa_data[c(-which(lpa_data$Leaf == "alder")),] 
leafpack_oak_fine<-leafpack_oak[c(-which(leafpack_oak$Mesh == "coarse")),] 
leafpack_oak_coarse<-leafpack_oak[c(-which(leafpack_oak$Mesh == "fine")),]

##****************************************************************
## Separate Alder

leafpack_alder<-lpa_data[c(-which(lpa_data$Leaf == "oak")),] 
leafpack_alder_fine<-leafpack_alder[c(-which(leafpack_alder$Mesh == "coarse")),] 
leafpack_alder_coarse<-leafpack_alder[c(-which(leafpack_alder$Mesh == "fine")),]

#########################################################################################
##                                                                                     ##
##                                FINE ALDER 2013                                      ##
##                                                                                     ##
#########################################################################################

## Take 2013 data
lpa_alder_fine_2013<-leafpack_alder_fine[c(-which(leafpack_alder_fine$Year == "2014")),] 

## Assess basic trends in breakdown data
head(lpa_alder_fine_2013)
row.names(lpa_alder_fine_2013) <- NULL

## Screen for extreme values
qqnorm(lpa_alder_fine_2013$k_day_ml, pch = 1, frame = FALSE)
qqline(lpa_alder_fine_2013$k_day_ml, col = "steelblue", lwd = 2)

## Two values of concern - Field notes do not suggest a problem
lpa_alder_fine_2013[c(34,106),c(1:12)]
## 2013 COL US2 Unknown rep 4 Incorrect initial leaf mass or human error - too low
## 2013 KER DS Sandra rep 3 Hole in bag or human error - way too high

## See improvement in plot
qqnorm(lpa_alder_fine_2013$k_day_ml[-c(34,106)], pch = 1, frame = FALSE)
qqline(lpa_alder_fine_2013$k_day_ml[-c(34,106)], col = "steelblue", lwd = 2)

## Reorder locations for further analysis
lpa_alder_fine_2013$Sample_location <- NA
lpa_alder_fine_2013[which(lpa_alder_fine_2013$Location=="DS"),16] <- "3_D"
lpa_alder_fine_2013[which(lpa_alder_fine_2013$Location=="US1"),16] <- "2_U1"
lpa_alder_fine_2013[which(lpa_alder_fine_2013$Location=="US2"),16] <- "1_U2"

## Reorganise before writing output
lpa_alder_fine_2013 <- lpa_alder_fine_2013[,c(1:3,16,4:15)]
head(lpa_alder_fine_2013)

## Write data and outlier to file
write.csv(lpa_alder_fine_2013[c(34,106),],"lpa_alder_fine_2013_outliers.csv", row.names = F)
write.csv(lpa_alder_fine_2013[-c(34,106),],"lpa_alder_fine_2013_data.csv", row.names = F)

#########################################################################################
##                                                                                     ##
##                                FINE OAK 2013                                        ##
##                                                                                     ##
#########################################################################################

## Take 2013 data
lpa_oak_fine_2013<-leafpack_oak_fine

## Assess basic trends in breakdown data
head(lpa_oak_fine_2013)
row.names(lpa_oak_fine_2013) <- NULL

## Screen for extreme values
qqnorm(lpa_oak_fine_2013$k_day_ml, pch = 1, frame = FALSE)
qqline(lpa_oak_fine_2013$k_day_ml, col = "steelblue", lwd = 2)

## Plot by degree days to see data spread
plot(log(k_day_ml)~DD,data=lpa_oak_fine_2013)

## Four values of concern - Field notes do not suggest a problem
lpa_oak_fine_2013[c(53,90,123,87),c(1:12)]
## 2013 HER DS Sandra rep 5 Hole in bag or human error - way too high  
## 2013 HOR US1 Sandra rep 2 Incorrect initial leaf mass or human error - way too low
## 2013 MES DS Sandra rep 4 Incorrect initial leaf mass or human error - way too low 
## 2013 HOR DS Sandra rep 5 Incorrect initial leaf mass or human error - way too low

## See improvement in plots
qqnorm(lpa_oak_fine_2013$k_day_ml[-c(53,90,123,87)], pch = 1, frame = FALSE)
qqline(lpa_oak_fine_2013$k_day_ml[-c(53,90,123,87)], col = "steelblue", lwd = 2)

## Plot by degree days to see data spread
plot(log(k_day_ml)~DD,data=lpa_oak_fine_2013[-c(53,90,123,87),])

## Reorder locations for further analysis
lpa_oak_fine_2013$Sample_location <- NA
lpa_oak_fine_2013[which(lpa_oak_fine_2013$Location=="DS"),16] <- "3_D"
lpa_oak_fine_2013[which(lpa_oak_fine_2013$Location=="US1"),16] <- "2_U1"
lpa_oak_fine_2013[which(lpa_oak_fine_2013$Location=="US2"),16] <- "1_U2"

## Reorganise before writing output
lpa_oak_fine_2013 <- lpa_oak_fine_2013[,c(1:3,16,4:15)]
head(lpa_oak_fine_2013)

## Write data for further analysis
write.csv(lpa_oak_fine_2013[-c(53,90,123,87),],"lpa_oak_fine_2013_data.csv",row.names = F)
write.csv(lpa_oak_fine_2013[c(53,90,123,87),],"lpa_oak_fine_2013_outliers.csv",row.names = F)

#########################################################################################
##                                                                                     ##
##                                COARSE ALDER 2013                                    ##
##                                                                                     ##
#########################################################################################

## Take 2013 data for temperature plots
lpa_alder_coarse_2013<-leafpack_alder_coarse[which(leafpack_alder_coarse$Year=="2013"),]

## Assess basic trends in breakdown data
head(lpa_alder_coarse_2013)
row.names(lpa_alder_coarse_2013) <- NULL

## Screen for extreme values
qqnorm(lpa_alder_coarse_2013$k_day_ml, pch = 1, frame = FALSE)
qqline(lpa_alder_coarse_2013$k_day_ml, col = "steelblue", lwd = 2)

## Plot by degree days to see data spread
plot(log(k_day_ml)~DD,data=lpa_alder_coarse_2013)

## Three values of concern - Field notes suggest reps 1-3 NIE may have been partly exposed (out of water)
## Reps 4-6 seem ok
## 2013	NIE	US1	Marta	rep 3 - Too low
## 2013	NIE	US1	Marta	rep 1 - Too low
## 2013	NIE	US1	Marta	rep 2 - Too low
## Two values of concern - Field notes do not suggest a problem
## 2013 MES DS Marta rep 6 - Incorrect initial value or human error - too low
## 2013	ROT	US1	Sandra rep 5 - Hole in bag or human error - way too high

## See improvement in plots
qqnorm(lpa_alder_coarse_2013$k_day_ml[-c(125,143,148,147,182)], pch = 1, frame = FALSE)
qqline(lpa_alder_coarse_2013$k_day_ml[-c(125,143,148,147,182)], col = "steelblue", lwd = 2)

## Plot by degree days to see data spread
plot(log(k_day_ml)~DD,data=lpa_alder_coarse_2013[-c(125,143,148,147,182),])

## Reorder locations for further analysis
lpa_alder_coarse_2013$Sample_location <- NA
lpa_alder_coarse_2013[which(lpa_alder_coarse_2013$Location=="DS"),16] <- "3_D"
lpa_alder_coarse_2013[which(lpa_alder_coarse_2013$Location=="US1"),16] <- "2_U1"
lpa_alder_coarse_2013[which(lpa_alder_coarse_2013$Location=="US2"),16] <- "1_U2"

## Reorganise before writing output
lpa_alder_coarse_2013 <- lpa_alder_coarse_2013[,c(1:3,16,4:15)]
head(lpa_alder_coarse_2013)

## Write data for further analysis
write.csv(lpa_alder_coarse_2013[c(125,143,148,147,182),],"lpa_alder_coarse_2013_outliers.csv",row.names = T)
write.csv(lpa_alder_coarse_2013[-c(125,143,148,147,182),],"lpa_alder_coarse_2013_data.csv",row.names = T)

#########################################################################################
##                                                                                     ##
##                                COARSE OAK 2013                                      ##
##                                                                                     ##
#########################################################################################

## Take 2013 data for temperature plots
lpa_oak_coarse_2013<-leafpack_oak_coarse

## Assess basic trends in breakdown data
head(lpa_oak_coarse_2013)
row.names(lpa_oak_coarse_2013) <- NULL

## Screen for extreme values
qqnorm(lpa_oak_coarse_2013$k_day_ml, pch = 1, frame = FALSE)
qqline(lpa_oak_coarse_2013$k_day_ml, col = "steelblue", lwd = 2)

## Plot by degree days to see data spread
plot(log1p(k_day_ml)~DD,data=lpa_oak_coarse_2013)

## One value of concern
lpa_oak_coarse_2013[131,c(1:12)]
## 2013 MES US1 Pravin  rep 5 Incorrect initial leaf mass or human error - way too low

## See improvement in plot
qqnorm(lpa_oak_coarse_2013$k_day_ml[-131], pch = 1, frame = FALSE)
qqline(lpa_oak_coarse_2013$k_day_ml[-131], col = "steelblue", lwd = 2)

## Plot by degree days to see data spread
plot(log1p(k_day_ml)~DD,data=lpa_oak_coarse_2013[-131,])

## Reorder locations for further analysis
lpa_oak_coarse_2013$Sample_location <- NA
lpa_oak_coarse_2013[which(lpa_oak_coarse_2013$Location=="DS"),16] <- "3_D"
lpa_oak_coarse_2013[which(lpa_oak_coarse_2013$Location=="US1"),16] <- "2_U1"
lpa_oak_coarse_2013[which(lpa_oak_coarse_2013$Location=="US2"),16] <- "1_U2"

## Reorganise before writing output
lpa_oak_coarse_2013 <- lpa_oak_coarse_2013[,c(1:3,16,4:15)]
head(lpa_oak_coarse_2013)

## Write data for further analysis
write.csv(lpa_oak_coarse_2013[131,],"lpa_oak_coarse_2013_outliers.csv",row.names = F)
write.csv(lpa_oak_coarse_2013[-131,],"lpa_oak_coarse_2013_data.csv",row.names = F)

#########################################################################################
##                                                                                     ##
##                              INVERTEBRATE ALDER 2013                                ##
##                                                                                     ##
#########################################################################################

## Take 2013 data for temperature plots
lpa_alder_coarse_2013<-leafpack_alder_coarse[which(leafpack_alder_coarse$Year=="2013"),]

## Reorder locations for further analysis
lpa_alder_coarse_2013$Sample_location <- NA
lpa_alder_coarse_2013[which(lpa_alder_coarse_2013$Location=="DS"),16] <- "3_D"
lpa_alder_coarse_2013[which(lpa_alder_coarse_2013$Location=="US1"),16] <- "2_U1"
lpa_alder_coarse_2013[which(lpa_alder_coarse_2013$Location=="US2"),16] <- "1_U2"

## Create mean microbial breakdown values for each location (remove outlier values)
lpa_alder_fine_2013_mu <- data.frame(lpa_alder_fine_2013[-c(34,106),] %>%
                                  group_by(Year,Site_code,Location) %>% 
                                  summarise_at(vars("Final_AFDM"), mean, na.rm=T))
## Merge values
lpa_alder_invert_2013_mu <- merge(lpa_alder_coarse_2013,lpa_alder_fine_2013_mu[,-1],by=c("Site_code","Location"))

## Rename columns after merge
names(lpa_alder_invert_2013_mu)[12] <- "Final_AFDM"
names(lpa_alder_invert_2013_mu)[17] <- "Fine_mesh_AFDM"

## Calculate Invertebrate-mediated mass loss
lpa_alder_invert_2013_mu$k_invertebrate <- -1*log(lpa_alder_invert_2013_mu$Final_AFDM/lpa_alder_invert_2013_mu$Fine_mesh_AFDM)

##***************************************************************************************
## Calculate mass loss (k days)

lpa_alder_invert_2013_mu$k_invert_day_ml  <- -1*log(lpa_alder_invert_2013_mu$Final_AFDM/lpa_alder_invert_2013_mu$Fine_mesh_AFDM)/lpa_alder_invert_2013_mu$Days 

##***************************************************************************************
## Calculate mass loss (k degree-days)

lpa_alder_invert_2013_mu$k_invert_dd_ml  <- -1*log(lpa_alder_invert_2013_mu$Final_AFDM/lpa_alder_invert_2013_mu$Fine_mesh_AFDM)/lpa_alder_invert_2013_mu$DD 

## Check that it worked
head(lpa_alder_invert_2013_mu)

## Plot shows invertebrate mediated breakdown alway less than total breakdown
plot(k_invertebrate~k_ml,lpa_alder_invert_2013_mu)
abline(0,1,col="blue")

## Convert negative values to zero (not many)
lpa_alder_invert_2013_mu[c(18:20)] <- replace(lpa_alder_invert_2013_mu[c(18:20)] ,
                                              lpa_alder_invert_2013_mu[c(18:20)]  < 0, 0)

head(lpa_alder_invert_2013_mu)

## Bring together the microbial and invertebrate k rates
lpa_alder_invert_2013_trophic_levels <- merge(lpa_alder_fine_2013[,c(1:4,6,8:11,14:16)],
                                              lpa_alder_invert_2013_mu[,c(1:3,5,7,16,18:20)],
                                              by=c("Year","Site_code","Location","Sample_location","Leaf","Rep"),
                                              all.x = T)
## Reorganise 
head(lpa_alder_invert_2013_trophic_levels)

## Rename rows
names(lpa_alder_invert_2013_trophic_levels)[10] <- "k_microbe"
names(lpa_alder_invert_2013_trophic_levels)[11] <- "k_microbe_day"
names(lpa_alder_invert_2013_trophic_levels)[12] <- "k_microbe_dd"
names(lpa_alder_invert_2013_trophic_levels)[13] <- "k_invert"
names(lpa_alder_invert_2013_trophic_levels)[14] <- "k_invert_day"
names(lpa_alder_invert_2013_trophic_levels)[15] <- "k_invert_dd"

## Need to remove coarse outliers (since the same problems affect the data derived from coarse mesh)
## Samples flagged in the field
## 2013	NIE	US1	Marta	rep 3 - Too low
## 2013	NIE	US1	Marta	rep 1 - Too low
## 2013	NIE	US1	Marta	rep 2 - Too low
## Two values of concern - Field notes do not suggest a problem
## 2013 MES DS Marta rep 6 - Incorrect initial value or human error - too low
## 2013	ROT	US1	Sandra rep 5 - Hole in bag or human error - way too high

## Write output ## NIE US1 rep 1 missing?
write.csv(lpa_alder_invert_2013_trophic_levels[c(144,145,127,181),],"lpa_alder_invertebrate_2013_outliers.csv",row.names = F)
write.csv(lpa_alder_invert_2013_trophic_levels[-c(144,145,127,181),],"lpa_alder_invertebrate_2013_data.csv",row.names = F)

#########################################################################################
##                                                                                     ##
##                              INVERTEBRATE OAK   2013                                ##
##                                                                                     ##
#########################################################################################

## Take 2013 data for temperature plots
lpa_oak_coarse_2013<-leafpack_oak_coarse

## Reorder locations for further analysis
lpa_oak_coarse_2013$Sample_location <- NA
lpa_oak_coarse_2013[which(lpa_oak_coarse_2013$Location=="DS"),16] <- "3_D"
lpa_oak_coarse_2013[which(lpa_oak_coarse_2013$Location=="US1"),16] <- "2_U1"
lpa_oak_coarse_2013[which(lpa_oak_coarse_2013$Location=="US2"),16] <- "1_U2"

## Create mean microbial breakdown values for each location (remove outlier values)
lpa_oak_fine_2013_mu <- data.frame(lpa_oak_fine_2013[-c(53,90,123,87),] %>%
                                       group_by(Year,Site_code,Location) %>% 
                                       summarise_at(vars("Final_AFDM"), mean, na.rm=T))

## Merge data
lpa_oak_invert_2013_mu <- merge(lpa_oak_coarse_2013,lpa_oak_fine_2013_mu[,-1],by=c("Site_code","Location"))

## Rename columns after merge
names(lpa_oak_invert_2013_mu)[12] <- "Final_AFDM"
names(lpa_oak_invert_2013_mu)[17] <- "Fine_mesh_AFDM"

head(lpa_oak_invert_2013_mu)

##***************************************************************************************
## Calculate Invertebrate-mediated mass loss

lpa_oak_invert_2013_mu$k_invertebrate <- -1*log(lpa_oak_invert_2013_mu$Final_AFDM/lpa_oak_invert_2013_mu$Fine_mesh_AFDM)

##***************************************************************************************
## Calculate mass loss (k days)

lpa_oak_invert_2013_mu$k_invert_day_ml  <- -1*log(lpa_oak_invert_2013_mu$Final_AFDM/lpa_oak_invert_2013_mu$Fine_mesh_AFDM)/lpa_oak_invert_2013_mu$Days 

##***************************************************************************************
## Calculate mass loss (k degree-days)

lpa_oak_invert_2013_mu$k_invert_dd_ml  <- -1*log(lpa_oak_invert_2013_mu$Final_AFDM/lpa_oak_invert_2013_mu$Fine_mesh_AFDM)/lpa_oak_invert_2013_mu$DD 

## Check it worked
head(lpa_oak_invert_2013_mu)

## Plot and see outlier
plot(k_invertebrate~k_ml,lpa_oak_invert_2013_mu)
abline(0,1,col="blue")

## Convert negative values to zero (not many)
lpa_oak_invert_2013_mu[c(18:20)] <- replace(lpa_oak_invert_2013_mu[c(18:20)] ,
                                              lpa_oak_invert_2013_mu[c(18:20)]  < 0, 0)

head(lpa_oak_fine_2013)
## Rationalise and write output
## Bring together the microbial and invertebrate k rates
lpa_oak_invert_2013_trophic_levels <- merge(lpa_oak_fine_2013[,c(1:4,6,8:11,14:16)],
                                              lpa_oak_invert_2013_mu[,c(1:3,5,7,16,18:20)],
                                              by=c("Year","Site_code","Location","Sample_location","Leaf","Rep"),
                                              all.x = T)
head(lpa_oak_invert_2013_trophic_levels)
## Rename rows
names(lpa_oak_invert_2013_trophic_levels)[10] <- "k_microbe"
names(lpa_oak_invert_2013_trophic_levels)[11] <- "k_microbe_day"
names(lpa_oak_invert_2013_trophic_levels)[12] <- "k_microbe_dd"
names(lpa_oak_invert_2013_trophic_levels)[13] <- "k_invert"
names(lpa_oak_invert_2013_trophic_levels)[14] <- "k_invert_day"
names(lpa_oak_invert_2013_trophic_levels)[15] <- "k_invert_dd"

## Need to remove coarse outliers (since the same problems affect the data derived from coarse mesh)
## 2013 MES US1 Pravin rep 5 Incorrect initial leaf mass or human error - way too low

## Write output
write.csv(lpa_oak_invert_2013_trophic_levels[c(129),],"lpa_oak_invertebrate_2013_outliers.csv",row.names = F)
write.csv(lpa_oak_invert_2013_trophic_levels[-c(129),],"lpa_oak_invertebrate_2013_data.csv",row.names = F)

#########################################################################################
##                                                                                     ##
##                                COARSE ALDER All yrs                                 ##
##                                                                                     ##
#########################################################################################

## Take Allyrs data for temperature plots
lpa_alder_coarse_Allyrs<-leafpack_alder_coarse
row.names(lpa_alder_coarse_Allyrs) <- NULL

## Screen data for extreme values
qqnorm(lpa_alder_coarse_Allyrs$k_day_ml, pch = 1, frame = FALSE)
qqline(lpa_alder_coarse_Allyrs$k_day_ml, col = "steelblue", lwd = 2)

## Plot values against DD to test effects
plot(log1p(k_day_ml)~DD,data=lpa_alder_coarse_Allyrs)

## Three values of concern - Field notes suggest reps 1-3 NIE may have been partly exposed (out of water)
## Reps 4-6 seem ok
## 2013	NIE	US1	Marta	rep 3 - Too low
## 2013	NIE	US1	Marta	rep 1 - Too low
## 2013	NIE	US1	Marta	rep 2 - Too low
## Two values of concern - Field notes do not suggest a problem
## 2013 MES DS Marta rep 6 - Incorrect initial value or human error - too low
## 2013	ROT	US1	Sandra rep 5 - Hole in bag or human error - way too high
## 2nd hand account from Marta suggest that some reps US1 MUR may have been partly exposed 
## (out of water) - no record of reps, and breakdown values seem consistent with what we would expect
## Also some rep at US2 MAR were affected by sedimentation
## For these more broad site problems I think we have to accept that they are baked into the results
## Without a specific note about an exact rep I don´t want to drop an entire sampling location
## Two 2014 values with specific flags
## 2014 KNO DS Marta rep 1 Incorrect initial value or human error - too low (-ve)
## 2014 MAR US2 Marta  rep 1 clearly confounded, buried in fine sediment   
lpa_alder_coarse_Allyrs[c(159,190),c(1:12)]

## See improvement in plot
qqnorm(lpa_alder_coarse_Allyrs$k_day_ml[-c(159,190,197,231,235,236,288)], pch = 1, frame = FALSE)
qqline(lpa_alder_coarse_Allyrs$k_day_ml[-c(159,190,197,231,235,236,288)], col = "steelblue", lwd = 2)

## Plot values against DD to test effects
plot(log1p(k_day_ml)~DD,data=lpa_alder_coarse_Allyrs[-c(159,190,197,231,235,236,288),])

## Add new order for locations
lpa_alder_coarse_Allyrs$Sample_location <- NA
lpa_alder_coarse_Allyrs[which(lpa_alder_coarse_Allyrs$Location=="DS"),16] <- "3_D"
lpa_alder_coarse_Allyrs[which(lpa_alder_coarse_Allyrs$Location=="US1"),16] <- "2_U1"
lpa_alder_coarse_Allyrs[which(lpa_alder_coarse_Allyrs$Location=="US2"),16] <- "1_U2"

## Reorganise before writing output
lpa_alder_coarse_Allyrs <- lpa_alder_coarse_Allyrs[,c(1:3,16,4:15)]
head(lpa_alder_coarse_Allyrs)

## Write data for further analysis
write.csv(lpa_alder_coarse_Allyrs[c(159,190,197,231,235,236,288),],"lpa_alder_coarse_All_outliers.csv",row.names = F)
write.csv(lpa_alder_coarse_Allyrs[-c(159,190,197,231,235,236,288),],"lpa_alder_coarse_All_data.csv",row.names = F)

