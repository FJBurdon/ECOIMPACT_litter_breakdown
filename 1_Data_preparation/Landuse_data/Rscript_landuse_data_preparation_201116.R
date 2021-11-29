############################################################################################
##                                                                                        ##  
##                                LANDUSE (% COVER)                                       ##  
##                                                                                        ##
############################################################################################

##  Script: ECOIMPACT landuse data preparation
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Assess the correctness of the data
##  2. Remove unnecessary columns and write file for further analysis

##***************************************************************************************
##Load libraries you need
library(reshape)
library(reshape2)
library(plyr)
library(tidyverse)
library(readr)
require(read_csv)
library(readr)

##***************************************************************************************
## Load data - prepared by Christian Stamm to compare versions of the dataset for different
## ECOIMPACT publications (Burdon et al. 2006, Munz et al. 2017, Burdon et al. 2019,
##                            Burdon et al. 2020, Burdon et al. in prep)

## Tab EcoImpact Comparison Land use data		
## Christian Stamm	 16.11.2020	
## Rationale: Different EcoImpact publications provide different land use data; this table provides a compilation across the papers to get an overview		
## Burdon_16: Burdon et al. 2016: Environmental context and magnitude of disturbance influence trait-mediated community responses to wastewater in streams	 Ecology and Evolution	 6(12)	 3923-3939 (file: ece32165-sup-0001-appendixa-c.docx);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
## Munz_17: Munz et al.	 2017: Pesticides drive risk of micropollutants in wastewater-impacted streams during low flow conditions	 Water Research	110	 366-377 (file: Munz_2017_SI.docx);;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## Burdon_19: Burdon et al. 2019: Agriculture versus wastewater pollution as drivers of macroinvertebrate community structure in streams	 Science of the Total Environment	659	 1256-1265;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;	
## Burdon_20: Burdon et al.	 2020: Stream microbial communities and ecosystem functioning show complex responses to multiple stressors in wastewater	 Global Change Biology	 26(11)	 6363-6382;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
## Burdon_Leaflitter_draft		
## Site: EcoImpact study site		
## Year: Year of the field campaign		
## Catchment: Area of the hydrological catchment [km2]		
## Arable: Areal fraction of arable land in the catchment		
## Pasture: Areal fraction of pasture land in the catchment		
## Forest: Areal fraction of forest in the catchment		
## Urban: Areal fraction of urban areas in the catchment		
## Orchards: Areal fraction of orchards in the catchment		

Landuse_all <- read_delim("Tab_EcoIMpact_Landuse_Compare_v20201116.csv", 
                          ";", escape_double = FALSE, trim_ws = TRUE, 
                          skip = 17)
Landuse_all  <- as.data.frame(Landuse_all)

str(Landuse_all)

#############################################################################################
##                                                                                         ##
##                             CONFIRM KEY RELATIONSHIPS                                   ##        
##                                                                                         ##    
#############################################################################################

head(Landuse_all)

## Arable cropping
## Strong congruence except for Colombier
plot(Arable_4 ~ Arable_2, data = Landuse_all,
     xlab="Burdon et al. 2019",
     ylab="Burdon et al. in prep",
     main="Arable cropping")
abline(0,1,col="blue")

## Pasture
## Strong congruence
plot(Pasture_4 ~ Pasture_2, data = Landuse_all,
     xlab="Burdon et al. 2019",
     ylab="Burdon et al. in prep",
     main="Pasture")
abline(0,1,col="blue")

## Forest
## Strong congruence except for Colombier
plot(Forest_4 ~ Forest_2, data = Landuse_all,
     xlab="Burdon et al. 2019",
     ylab="Burdon et al. in prep",
     main="Forest")
abline(0,1,col="blue")

## Urban
## Mostly strong congruence - but tend to overestimate 
plot(Urban_4 ~ Urban_2, data = Landuse_all,
     xlab="Burdon et al. 2019",
     ylab="Burdon et al. in prep",
     main="Urban")
abline(0,1,col="blue")

## Orchards
## Weak congruence - but tend to underestimate 
plot(Orchards_4 ~ Orchards_2, data = Landuse_all,
     xlab="Burdon et al. 2019",
     ylab="Burdon et al. in prep",
     main="Urban")
abline(0,1,col="blue")

## Arable vs Orchards
plot(Orchards_2 ~ Arable_2, data = Landuse_all,
     xlab="Arable",
     ylab="Orchards",
     main="Burdon et al. 2019 - Arable vs Orchards")
abline(lm(Orchards_2 ~ Arable_2, data = Landuse_all),col="red")

## Catchment area
## Strong congruence - but tend to overestimate 
plot(Landuse_all$`Catchment area_3` ~ Landuse_all$`Catchment area_1`,
     xlab="Burdon et al. 2019",
     ylab="Burdon et al. in prep",
     main="Catchment area")
abline(0,1,col="blue")

#############################################################################################
##                                                                                         ##
##                    RATIONALISE DATA AND WRITE FOR FURTHER ANALYSIS                      ##        
##                                                                                         ##    
#############################################################################################

Landuse_temp <- Landuse_all[,c(1:3,16:21)]

head(Landuse_temp)

unique(Landuse_temp$Site)

## Add column and rearrange
Landuse_temp$SiteName <- NA
Landuse_temp <- Landuse_temp[,c(10,1:9)]
  
## Rename sites by three letter codes
Landuse_temp[which(Landuse_temp$Site=="Buttisholz"),1] <- "BUT"
Landuse_temp[which(Landuse_temp$Site=="Colombier"),1] <- "COL"
Landuse_temp[which(Landuse_temp$Site=="D\xfcrnten"),1] <- "DUR"
Landuse_temp[which(Landuse_temp$Site=="Herisau"),1] <- "HER"
Landuse_temp[which(Landuse_temp$Site=="Hochdorf"),1] <- "HOC"
Landuse_temp[which(Landuse_temp$Site=="Hornussen"),1] <- "HOR"
Landuse_temp[which(Landuse_temp$Site=="Kernenried"),1] <- "KER"
Landuse_temp[which(Landuse_temp$Site=="Messen"),1] <- "MES"
Landuse_temp[which(Landuse_temp$Site=="Niederdorf"),1] <- "NIE"
Landuse_temp[which(Landuse_temp$Site=="Romont"),1] <- "ROM"
Landuse_temp[which(Landuse_temp$Site=="Rothenturm"),1] <- "ROT"
Landuse_temp[which(Landuse_temp$Site=="S\xe9very"),1] <- "SEV"

Landuse_temp[which(Landuse_temp$Site=="Aadorf"),1] <- "AAD"
Landuse_temp[which(Landuse_temp$Site=="Birmensdorf"),1] <- "BIR"
Landuse_temp[which(Landuse_temp$Site=="Elgg"),1] <- "ELG"
Landuse_temp[which(Landuse_temp$Site=="Ellikon"),1] <- "ELL"
Landuse_temp[which(Landuse_temp$Site=="Knonau"),1] <- "KNO"
Landuse_temp[which(Landuse_temp$Site=="Marthalen"),1] <- "MAR"
Landuse_temp[which(Landuse_temp$Site=="Muri"),1] <- "MUR"
Landuse_temp[which(Landuse_temp$Site=="Reinach"),1] <- "REI"
Landuse_temp[which(Landuse_temp$Site=="Unterehrendingen"),1] <- "UNT"
Landuse_temp[which(Landuse_temp$Site=="Val-de-Ruz"),1] <- "VAL"
Landuse_temp[which(Landuse_temp$Site=="Villeret"),1] <- "VIL"
Landuse_temp[which(Landuse_temp$Site=="Zullwil"),1] <- "ZUL"

## Rename DUR and SEV since above did not work
Landuse_temp[3,1] <- "DUR"
Landuse_temp[12,1] <- "SEV"

## Add column and rearrange
Landuse_temp$Location <- NA
Landuse_temp <- Landuse_temp[,c(2,4,1,11,5,6,10,7:9)]

## Check reorganization
head(Landuse_temp)

## Rename columns
colnames(Landuse_temp) <- c("Site_ID","Year","Site","Location","CA","p_Crop","p_Orchard",
                             "p_Pasture","p_Forested","p_Urban")

## Create dummy variables for site data
Landuse_D   <- Landuse_temp
Landuse_U1  <- Landuse_temp
Landuse_U2  <- Landuse_temp
  
## Rename site location
Landuse_D$Location   <- "D"
Landuse_U1$Location   <- "U1"
Landuse_U2$Location   <- "U2"

## Combine data and reorganise
Landuse_final <- rbind(Landuse_D,Landuse_U1,Landuse_U2)
Landuse_final <- arrange(Landuse_final, Year, Site, Location)

## Tidyup workspace
rm(Landuse_D,Landuse_U1,Landuse_U2,Landuse_temp)

## Write output for use in SEM
write.csv(Landuse_final,"5_DATA_Landuse_201116.csv",row.names = F)
