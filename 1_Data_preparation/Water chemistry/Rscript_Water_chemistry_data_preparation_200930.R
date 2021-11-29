############################################################################################
##                                                                                        ##  
##                                   WATER CHEMISTRY                                      ##  
##                                                                                        ##
############################################################################################

##  Script: ECOIMPACT leafpack assay data preparation
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Demonstrate that there is no difference in water chemistry of the US sites

##***************************************************************************************
##Load libraries you need
library(reshape)
library(reshape2)
library(plyr)
library(tidyverse)
library(readr)
library(read_csv)

##***************************************************************************************
## Load data - all litter types
Water_chem <-read_csv("4_DATA_general_water_chemistry_160824.csv")
Water_chem <- as.data.frame(Water_chem)
str(Water_chem)

#############################################################################################
##                                                                                         ##
##                      TRANSFORM DATA AND RUN PERMANOVA                                   ##        
##                                                                                         ##    
#############################################################################################

Water_chem_lpa <- Water_chem[-c(which(Water_chem$Site=="AAD"),
                                     which(Water_chem$Site=="BIR"),
                                           which(Water_chem$Site=="UNT"),
                                                 which(Water_chem$Site=="ZUL")),]


Water_chem_lpa_log <- log(Water_chem_lpa[,-c(1:3)])

#############################################################################################
##                                                                                         ##
##                                       ADONIS                                            ##        
##                                                                                         ##    
#############################################################################################

library(vegan)
library(pairwiseAdonis)

## Create model from above data
dis_invert <- vegdist(Water_chem_lpa_log,"euc")

## Pairwise adonis: note strata only works for UBF vs FBF (site pairs)
P1 <- pairwise.adonis2(dis_invert ~ Location + Year, data = Water_chem_lpa[,c(1:3)], strata = 'Site')
P1

#############################################################################################
##                                                                                         ##
##                                     BETADISPER                                          ##        
##                                                                                         ##    
#############################################################################################

## Betadisper
B1 <- betadisper(dis_invert, Water_chem_lpa$Location, type = c("median"), bias.adjust = T)
## Test using permutation 
permutest(B1, permutations = 999, strata = Water_chem_lpa$Location)
permutest(B1, permutations = 999) ## Remove strata
(mod.HSD <- TukeyHSD(B1))
plot(mod.HSD)

##**********************************************************************
# Plot Betadisper results of Spiders
# Hellinger-transformed

png(filename="Water_chemistry_log_transformed_betadisper.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=16, 
    res=600)

par(mfcol = c(1,2), mar = c(4.6, 3, 4, 0.2), pty="s")

plot(B1,main="Water chemistry (Log-transform)")
boxplot(B1,xlab ="Location")

dev.off()

#############################################################################################
##                                                                                         ##
##                                       NMDS                                              ##        
##                                                                                         ##    
#############################################################################################

## Run NMDS on log-transformed distance matrix
mod1 <- metaMDS(dis_invert, dist="euc", trymax = 20, autotransform = T)

## Run envfit to test effect of location
envfit(mod1 ~ Location, data = Water_chem_lpa, p=999)

## Calculate distance between centroids
## Distance between U1 and U2
U1_U2 <- sqrt((0.7204-0.8932)^2+(-0.1458--0.0864)^2)
## 0.1827244

## Distance between D and U1
D_U1 <- sqrt((-1.6136-0.7204)^2+(0.2322--0.1458)^2)
## 2.364411

## Distance between D and U2
D_U2 <- sqrt((-1.6136-0.8932)^2+(0.2322--0.0864)^2)
## 2.526965

## Calculate difference as a %
D_U1/U1_U2*100
## 1294%

D_U2/U1_U2*100
## 1383%

#############################################################################################
##                                                                                         ##
##                                   EFFECT SIZE                                           ##        
##                                                                                         ##    
#############################################################################################

## Demonstrate how small the effect of change is
require(SingleCaseES)

U1_U2 <- Water_chem_lpa[-which(Water_chem_lpa$Location=="D"),]
D_U1 <- Water_chem_lpa[-which(Water_chem_lpa$Location=="U2"),]
D_U2 <- Water_chem_lpa[-which(Water_chem_lpa$Location=="U1"),]

colnames(Water_chem_lpa)

##********************************************************************************************
## COND
colnames(Water_chem_lpa)
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Cond, 
                              improvement = "increase")

U1_U2_Cond_Effsize <- data.frame("U1_U2","Cond",Unique_pair[3,])
names(U1_U2_Cond_Effsize)[1] <- "Contrast" 
names(U1_U2_Cond_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Cond, 
                              improvement = "increase")

D_U1_Cond_Effsize <- data.frame("D_U1","Cond",Unique_pair[3,])
names(D_U1_Cond_Effsize)[1] <- "Contrast" 
names(D_U1_Cond_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Cond, 
                              improvement = "increase")

D_U2_Cond_Effsize <- data.frame("D_U2","Cond",Unique_pair[3,])
names(D_U2_Cond_Effsize)[1] <- "Contrast" 
names(D_U2_Cond_Effsize)[2] <- "Response" 

Cond <- rbind(D_U1_Cond_Effsize,D_U2_Cond_Effsize,U1_U2_Cond_Effsize)

##********************************************************************************************
## pH
colnames(Water_chem_lpa[5])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = pH, 
                              improvement = "increase")

U1_U2_pH_Effsize <- data.frame("U1_U2","pH",Unique_pair[3,])
names(U1_U2_pH_Effsize)[1] <- "Contrast" 
names(U1_U2_pH_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = pH, 
                              improvement = "increase")

D_U1_pH_Effsize <- data.frame("D_U1","pH",Unique_pair[3,])
names(D_U1_pH_Effsize)[1] <- "Contrast" 
names(D_U1_pH_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = pH, 
                              improvement = "increase")

D_U2_pH_Effsize <- data.frame("D_U2","pH",Unique_pair[3,])
names(D_U2_pH_Effsize)[1] <- "Contrast" 
names(D_U2_pH_Effsize)[2] <- "Response" 

pH <- rbind(D_U1_pH_Effsize,D_U2_pH_Effsize,U1_U2_pH_Effsize)

##********************************************************************************************
## Alkal
colnames(Water_chem_lpa[6])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Alkal, 
                              improvement = "increase")

U1_U2_Alkal_Effsize <- data.frame("U1_U2","Alkal",Unique_pair[3,])
names(U1_U2_Alkal_Effsize)[1] <- "Contrast" 
names(U1_U2_Alkal_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Alkal, 
                              improvement = "increase")

D_U1_Alkal_Effsize <- data.frame("D_U1","Alkal",Unique_pair[3,])
names(D_U1_Alkal_Effsize)[1] <- "Contrast" 
names(D_U1_Alkal_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Alkal, 
                              improvement = "increase")

D_U2_Alkal_Effsize <- data.frame("D_U2","Alkal",Unique_pair[3,])
names(D_U2_Alkal_Effsize)[1] <- "Contrast" 
names(D_U2_Alkal_Effsize)[2] <- "Response" 

Alkal <- rbind(D_U1_Alkal_Effsize,D_U2_Alkal_Effsize,U1_U2_Alkal_Effsize)

##********************************************************************************************
## Hard
colnames(Water_chem_lpa[7])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Hard, 
                              improvement = "increase")

U1_U2_Hard_Effsize <- data.frame("U1_U2","Hard",Unique_pair[3,])
names(U1_U2_Hard_Effsize)[1] <- "Contrast" 
names(U1_U2_Hard_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Hard, 
                              improvement = "increase")

D_U1_Hard_Effsize <- data.frame("D_U1","Hard",Unique_pair[3,])
names(D_U1_Hard_Effsize)[1] <- "Contrast" 
names(D_U1_Hard_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Hard, 
                              improvement = "increase")

D_U2_Hard_Effsize <- data.frame("D_U2","Hard",Unique_pair[3,])
names(D_U2_Hard_Effsize)[1] <- "Contrast" 
names(D_U2_Hard_Effsize)[2] <- "Response" 

Hard <- rbind(D_U1_Hard_Effsize,D_U2_Hard_Effsize,U1_U2_Hard_Effsize)

##********************************************************************************************
## Na
colnames(Water_chem_lpa[8])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Na, 
                              improvement = "increase")

U1_U2_Sodium_Effsize <- data.frame("U1_U2","Sodium",Unique_pair[3,])
names(U1_U2_Sodium_Effsize)[1] <- "Contrast" 
names(U1_U2_Sodium_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Na, 
                              improvement = "increase")

D_U1_Sodium_Effsize <- data.frame("D_U1","Sodium",Unique_pair[3,])
names(D_U1_Sodium_Effsize)[1] <- "Contrast" 
names(D_U1_Sodium_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Na, 
                              improvement = "increase")

D_U2_Sodium_Effsize <- data.frame("D_U2","Sodium",Unique_pair[3,])
names(D_U2_Sodium_Effsize)[1] <- "Contrast" 
names(D_U2_Sodium_Effsize)[2] <- "Response" 

Sodium <- rbind(D_U1_Sodium_Effsize,D_U2_Sodium_Effsize,U1_U2_Sodium_Effsize)

##********************************************************************************************
## K
colnames(Water_chem_lpa[9])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = K, 
                              improvement = "increase")

U1_U2_Potassium_Effsize <- data.frame("U1_U2","Potassium",Unique_pair[3,])
names(U1_U2_Potassium_Effsize)[1] <- "Contrast" 
names(U1_U2_Potassium_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = K, 
                              improvement = "increase")

D_U1_Potassium_Effsize <- data.frame("D_U1","Potassium",Unique_pair[3,])
names(D_U1_Potassium_Effsize)[1] <- "Contrast" 
names(D_U1_Potassium_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = K, 
                              improvement = "increase")

D_U2_Potassium_Effsize <- data.frame("D_U2","Potassium",Unique_pair[3,])
names(D_U2_Potassium_Effsize)[1] <- "Contrast" 
names(D_U2_Potassium_Effsize)[2] <- "Response" 

Potassium <- rbind(D_U1_Potassium_Effsize,D_U2_Potassium_Effsize,U1_U2_Potassium_Effsize)

##********************************************************************************************
## Ca
colnames(Water_chem_lpa[10])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Ca, 
                              improvement = "increase")

U1_U2_Calcium_Effsize <- data.frame("U1_U2","Calcium",Unique_pair[3,])
names(U1_U2_Calcium_Effsize)[1] <- "Contrast" 
names(U1_U2_Calcium_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Ca, 
                              improvement = "increase")

D_U1_Calcium_Effsize <- data.frame("D_U1","Calcium",Unique_pair[3,])
names(D_U1_Calcium_Effsize)[1] <- "Contrast" 
names(D_U1_Calcium_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Ca, 
                              improvement = "increase")

D_U2_Calcium_Effsize <- data.frame("D_U2","Calcium",Unique_pair[3,])
names(D_U2_Calcium_Effsize)[1] <- "Contrast" 
names(D_U2_Calcium_Effsize)[2] <- "Response" 

Calcium <- rbind(D_U1_Calcium_Effsize,D_U2_Calcium_Effsize,U1_U2_Calcium_Effsize)

##********************************************************************************************
## Mg
colnames(Water_chem_lpa[11])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Mg, 
                              improvement = "increase")

U1_U2_Manganese_Effsize <- data.frame("U1_U2","Manganese",Unique_pair[3,])
names(U1_U2_Manganese_Effsize)[1] <- "Contrast" 
names(U1_U2_Manganese_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Mg, 
                              improvement = "increase")

D_U1_Manganese_Effsize <- data.frame("D_U1","Manganese",Unique_pair[3,])
names(D_U1_Manganese_Effsize)[1] <- "Contrast" 
names(D_U1_Manganese_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Mg, 
                              improvement = "increase")

D_U2_Manganese_Effsize <- data.frame("D_U2","Manganese",Unique_pair[3,])
names(D_U2_Manganese_Effsize)[1] <- "Contrast" 
names(D_U2_Manganese_Effsize)[2] <- "Response" 

Manganese <- rbind(D_U1_Manganese_Effsize,D_U2_Manganese_Effsize,U1_U2_Manganese_Effsize)

##********************************************************************************************
## NH4
colnames(Water_chem_lpa[12])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NH4, 
                              improvement = "increase")

U1_U2_NH4_Effsize <- data.frame("U1_U2","NH4",Unique_pair[3,])
names(U1_U2_NH4_Effsize)[1] <- "Contrast" 
names(U1_U2_NH4_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NH4, 
                              improvement = "increase")

D_U1_NH4_Effsize <- data.frame("D_U1","NH4",Unique_pair[3,])
names(D_U1_NH4_Effsize)[1] <- "Contrast" 
names(D_U1_NH4_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NH4, 
                              improvement = "increase")

D_U2_NH4_Effsize <- data.frame("D_U2","NH4",Unique_pair[3,])
names(D_U2_NH4_Effsize)[1] <- "Contrast" 
names(D_U2_NH4_Effsize)[2] <- "Response" 

NH4 <- rbind(D_U1_NH4_Effsize,D_U2_NH4_Effsize,U1_U2_NH4_Effsize)

##********************************************************************************************
## NO2
colnames(Water_chem_lpa[13])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NO2, 
                              improvement = "increase")

U1_U2_NO2_Effsize <- data.frame("U1_U2","NO2",Unique_pair[3,])
names(U1_U2_NO2_Effsize)[1] <- "Contrast" 
names(U1_U2_NO2_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NO2, 
                              improvement = "increase")

D_U1_NO2_Effsize <- data.frame("D_U1","NO2",Unique_pair[3,])
names(D_U1_NO2_Effsize)[1] <- "Contrast" 
names(D_U1_NO2_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NO2, 
                              improvement = "increase")

D_U2_NO2_Effsize <- data.frame("D_U2","NO2",Unique_pair[3,])
names(D_U2_NO2_Effsize)[1] <- "Contrast" 
names(D_U2_NO2_Effsize)[2] <- "Response" 

NO2 <- rbind(D_U1_NO2_Effsize,D_U2_NO2_Effsize,U1_U2_NO2_Effsize)

##********************************************************************************************
## NO3
colnames(Water_chem_lpa[14])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NO3, 
                              improvement = "increase")

U1_U2_NO3_Effsize <- data.frame("U1_U2","NO3",Unique_pair[3,])
names(U1_U2_NO3_Effsize)[1] <- "Contrast" 
names(U1_U2_NO3_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NO3, 
                              improvement = "increase")

D_U1_NO3_Effsize <- data.frame("D_U1","NO3",Unique_pair[3,])
names(D_U1_NO3_Effsize)[1] <- "Contrast" 
names(D_U1_NO3_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = NO3, 
                              improvement = "increase")

D_U2_NO3_Effsize <- data.frame("D_U2","NO3",Unique_pair[3,])
names(D_U2_NO3_Effsize)[1] <- "Contrast" 
names(D_U2_NO3_Effsize)[2] <- "Response" 

NO3 <- rbind(D_U1_NO3_Effsize,D_U2_NO3_Effsize,U1_U2_NO3_Effsize)

##********************************************************************************************
## TN
colnames(Water_chem_lpa[15])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TN, 
                              improvement = "increase")

U1_U2_TN_Effsize <- data.frame("U1_U2","TN",Unique_pair[3,])
names(U1_U2_TN_Effsize)[1] <- "Contrast" 
names(U1_U2_TN_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TN, 
                              improvement = "increase")

D_U1_TN_Effsize <- data.frame("D_U1","TN",Unique_pair[3,])
names(D_U1_TN_Effsize)[1] <- "Contrast" 
names(D_U1_TN_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TN, 
                              improvement = "increase")

D_U2_TN_Effsize <- data.frame("D_U2","TN",Unique_pair[3,])
names(D_U2_TN_Effsize)[1] <- "Contrast" 
names(D_U2_TN_Effsize)[2] <- "Response" 

TN <- rbind(D_U1_TN_Effsize,D_U2_TN_Effsize,U1_U2_TN_Effsize)

##********************************************************************************************
## SRP
colnames(Water_chem_lpa[16])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = SRP, 
                              improvement = "increase")

U1_U2_SRP_Effsize <- data.frame("U1_U2","SRP",Unique_pair[3,])
names(U1_U2_SRP_Effsize)[1] <- "Contrast" 
names(U1_U2_SRP_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = SRP, 
                              improvement = "increase")

D_U1_SRP_Effsize <- data.frame("D_U1","SRP",Unique_pair[3,])
names(D_U1_SRP_Effsize)[1] <- "Contrast" 
names(D_U1_SRP_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = SRP, 
                              improvement = "increase")

D_U2_SRP_Effsize <- data.frame("D_U2","SRP",Unique_pair[3,])
names(D_U2_SRP_Effsize)[1] <- "Contrast" 
names(D_U2_SRP_Effsize)[2] <- "Response" 

SRP <- rbind(D_U1_SRP_Effsize,D_U2_SRP_Effsize,U1_U2_SRP_Effsize)

##********************************************************************************************
## TP
colnames(Water_chem_lpa[17])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TP, 
                              improvement = "increase")

U1_U2_TP_Effsize <- data.frame("U1_U2","TP",Unique_pair[3,])
names(U1_U2_TP_Effsize)[1] <- "Contrast" 
names(U1_U2_TP_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TP, 
                              improvement = "increase")

D_U1_TP_Effsize <- data.frame("D_U1","TP",Unique_pair[3,])
names(D_U1_TP_Effsize)[1] <- "Contrast" 
names(D_U1_TP_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TP, 
                              improvement = "increase")

D_U2_TP_Effsize <- data.frame("D_U2","TP",Unique_pair[3,])
names(D_U2_TP_Effsize)[1] <- "Contrast" 
names(D_U2_TP_Effsize)[2] <- "Response" 

TP <- rbind(D_U1_TP_Effsize,D_U2_TP_Effsize,U1_U2_TP_Effsize)

##********************************************************************************************
## Chloride
colnames(Water_chem_lpa[18])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Ch, 
                              improvement = "increase")

U1_U2_Chloride_Effsize <- data.frame("U1_U2","Chloride",Unique_pair[3,])
names(U1_U2_Chloride_Effsize)[1] <- "Contrast" 
names(U1_U2_Chloride_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Ch, 
                              improvement = "increase")

D_U1_Chloride_Effsize <- data.frame("D_U1","Chloride",Unique_pair[3,])
names(D_U1_Chloride_Effsize)[1] <- "Contrast" 
names(D_U1_Chloride_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Ch, 
                              improvement = "increase")

D_U2_Chloride_Effsize <- data.frame("D_U2","Chloride",Unique_pair[3,])
names(D_U2_Chloride_Effsize)[1] <- "Contrast" 
names(D_U2_Chloride_Effsize)[2] <- "Response" 

Chloride <- rbind(D_U1_Chloride_Effsize,D_U2_Chloride_Effsize,U1_U2_Chloride_Effsize)

##********************************************************************************************
## Sulphur
colnames(Water_chem_lpa[19])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = S, 
                              improvement = "increase")

U1_U2_Sulphur_Effsize <- data.frame("U1_U2","Sulphur",Unique_pair[3,])
names(U1_U2_Sulphur_Effsize)[1] <- "Contrast" 
names(U1_U2_Sulphur_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = S, 
                              improvement = "increase")

D_U1_Sulphur_Effsize <- data.frame("D_U1","Sulphur",Unique_pair[3,])
names(D_U1_Sulphur_Effsize)[1] <- "Contrast" 
names(D_U1_Sulphur_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = S, 
                              improvement = "increase")

D_U2_Sulphur_Effsize <- data.frame("D_U2","Sulphur",Unique_pair[3,])
names(D_U2_Sulphur_Effsize)[1] <- "Contrast" 
names(D_U2_Sulphur_Effsize)[2] <- "Response" 

Sulphur <- rbind(D_U1_Sulphur_Effsize,D_U2_Sulphur_Effsize,U1_U2_Sulphur_Effsize)

##********************************************************************************************
## TOC
colnames(Water_chem_lpa[20])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TOC, 
                              improvement = "increase")

U1_U2_TOC_Effsize <- data.frame("U1_U2","TOC",Unique_pair[3,])
names(U1_U2_TOC_Effsize)[1] <- "Contrast" 
names(U1_U2_TOC_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TOC, 
                              improvement = "increase")

D_U1_TOC_Effsize <- data.frame("D_U1","TOC",Unique_pair[3,])
names(D_U1_TOC_Effsize)[1] <- "Contrast" 
names(D_U1_TOC_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = TOC, 
                              improvement = "increase")

D_U2_TOC_Effsize <- data.frame("D_U2","TOC",Unique_pair[3,])
names(D_U2_TOC_Effsize)[1] <- "Contrast" 
names(D_U2_TOC_Effsize)[2] <- "Response" 

TOC <- rbind(D_U1_TOC_Effsize,D_U2_TOC_Effsize,U1_U2_TOC_Effsize)

##********************************************************************************************
## DOC
colnames(Water_chem_lpa[21])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = DOC, 
                              improvement = "increase")

U1_U2_DOC_Effsize <- data.frame("U1_U2","DOC",Unique_pair[3,])
names(U1_U2_DOC_Effsize)[1] <- "Contrast" 
names(U1_U2_DOC_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = DOC, 
                              improvement = "increase")

D_U1_DOC_Effsize <- data.frame("D_U1","DOC",Unique_pair[3,])
names(D_U1_DOC_Effsize)[1] <- "Contrast" 
names(D_U1_DOC_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = DOC, 
                              improvement = "increase")

D_U2_DOC_Effsize <- data.frame("D_U2","DOC",Unique_pair[3,])
names(D_U2_DOC_Effsize)[1] <- "Contrast" 
names(D_U2_DOC_Effsize)[2] <- "Response" 

DOC <- rbind(D_U1_DOC_Effsize,D_U2_DOC_Effsize,U1_U2_DOC_Effsize)

##********************************************************************************************
## Silica
colnames(Water_chem_lpa[22])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Silica, 
                              improvement = "increase")

U1_U2_Silica_Effsize <- data.frame("U1_U2","Silica",Unique_pair[3,])
names(U1_U2_Silica_Effsize)[1] <- "Contrast" 
names(U1_U2_Silica_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Silica, 
                              improvement = "increase")

D_U1_Silica_Effsize <- data.frame("D_U1","Silica",Unique_pair[3,])
names(D_U1_Silica_Effsize)[1] <- "Contrast" 
names(D_U1_Silica_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = Silica, 
                              improvement = "increase")

D_U2_Silica_Effsize <- data.frame("D_U2","Silica",Unique_pair[3,])
names(D_U2_Silica_Effsize)[1] <- "Contrast" 
names(D_U2_Silica_Effsize)[2] <- "Response" 

Silica <- rbind(D_U1_Silica_Effsize,D_U2_Silica_Effsize,U1_U2_Silica_Effsize)

##********************************************************************************************
## SS
colnames(Water_chem_lpa[23])
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = U1_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = SS, 
                              improvement = "increase")

U1_U2_SS_Effsize <- data.frame("U1_U2","SS",Unique_pair[3,])
names(U1_U2_SS_Effsize)[1] <- "Contrast" 
names(U1_U2_SS_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat = D_U1,
                              #grouping = Site, 
                              condition = Location,
                              outcome = SS, 
                              improvement = "increase")

D_U1_SS_Effsize <- data.frame("D_U1","SS",Unique_pair[3,])
names(D_U1_SS_Effsize)[1] <- "Contrast" 
names(D_U1_SS_Effsize)[2] <- "Response" 
## Calculate effect size between site pairs
Unique_pair  <- batch_calc_ES(dat =D_U2,
                              #grouping = Site, 
                              condition = Location,
                              outcome = SS, 
                              improvement = "increase")

D_U2_SS_Effsize <- data.frame("D_U2","SS",Unique_pair[3,])
names(D_U2_SS_Effsize)[1] <- "Contrast" 
names(D_U2_SS_Effsize)[2] <- "Response" 

SS <- rbind(D_U1_SS_Effsize,D_U2_SS_Effsize,U1_U2_SS_Effsize)

##********************************************************************************************
## Combine for further analysis
colnames(Water_chem_lpa)

Water_chem_effects <- rbind(Cond,
                             pH,
                             Alkal,
                             Hard,
                             Sodium,
                             Potassium,
                             Calcium,
                             Manganese,
                             NH4,
                             NO2,
                             NO3,
                             TN,
                             SRP,
                             TP,
                             Chloride,
                             Sulphur,
                             TOC,
                             DOC,
                             Silica,
                             SS)

rownames(Water_chem_effects) <- NULL

head(Water_chem_effects)

## Spread from wide to long
Water_chem_effects_wide <- Water_chem_effects[,c(1,2,4)] %>% spread(Contrast,Est)

## Calculate mean effects - shows that there is virtually no change between U1 and U2
colMeans(Water_chem_effects_wide[,-1])
## > colMeans(Water_chem_effects_wide[,-1])
## D_U1             D_U2            U1_U2 
## 1.170456291      1.201212238     -0.008774653 
                             