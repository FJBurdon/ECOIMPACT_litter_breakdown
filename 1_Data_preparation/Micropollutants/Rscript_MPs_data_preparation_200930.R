############################################################################################
##                                                                                        ##  
##                                MICROPOLLUTANTS (MPs)                                   ##  
##                                                                                        ##
############################################################################################

##  Script: ECOIMPACT leafpack assay data preparation
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Assess the correctness of the data
##  2. Remove unnecessary columns and write file for further analysis
##  3. Compare data with that used in Munz et al. 2017

## NOTES: TUs summed for compound classes with different levels of characterisation 
## (see Meta-substances 1-5), average values for summed totals are used

##***************************************************************************************
##Load libraries you need
library(reshape)
library(reshape2)
library(plyr)
library(tidyverse)
library(readr)
require(read_csv)

##***************************************************************************************
## Load data - this was processed in Excel using pivot tabling
MPs<-read_csv("3_DATA_MPs_TU_metasubstance_VER02_161114.csv")
MPs<- as.data.frame(MPs)
str(MPs)

## Load raw data the above is based upon
MP_metadata<-read_csv("DATA_TU_metasubstance_201106.csv")
MP_metadata<- as.data.frame(MP_metadata)
str(MP_metadata)

## Load data used in Munz et al. 2017 (also includes above data)
MP_Ecoimpact <- read_csv("Munz et al. 2017/spear_median_sumTU_msPAF_workCS_201116.csv", 
                      skip = 4)

## Load data from Munz et al. 2017 - different version of the above
MP_Ecoimpact_sumTU_msPAF <- read_csv("Munz et al. 2017/spear_median_sumTU_msPAF.csv")
MP_Ecoimpact_sumTU_msPAF <- MP_Ecoimpact_sumTU_msPAF[,-1]

## Macroinvertebrate data
MP_Ecoimpact_invert <- read_delim("Munz et al. 2017/Tab_Delta_allMacro_v20181116.csv", 
                                           ";", escape_double = FALSE, trim_ws = TRUE, 
                                           skip = 19)

#############################################################################################
##                                                                                         ##
##                             CONFIRM KEY RELATIONSHIPS                                   ##        
##                                                                                         ##    
#############################################################################################

## Confirm that non-pesticide TUs are nested within all TUs, but make up signficant fraction
## See concentration at upper end of scale due to metals in 2014
plot(log(TU_avg_nonPest_MS1) ~ log(TUs), data = MPs)
abline(0,1,col="blue")

## But metals only measured in 2014
MPs_2013 <- MPs[which(MPs$Year=="2013"),]

## Without metals in 2013, non-pesticides fall off 1:1 line with total TUs
## and the correlation is weak
plot(log(TU_avg_nonPest_MS1) ~ log(TUs), data = MPs_2013[which(MPs_2013$Location=="D"),])
abline(0,1,col="blue")
summary(lm(log(TU_avg_nonPest_MS1) ~ log(TUs), data = MPs_2013[which(MPs_2013$Location=="D"),]))

## Without metals in 2013, pesticides are off the 1:1 line with total TUs
## but the correlation is strong
plot(log(TU_avg_Pest_MS1) ~ log(TUs), data = MPs_2013[which(MPs_2013$Location=="D"),])
abline(0,1,col="blue")
summary(lm(log(TU_avg_Pest_MS1) ~ log(TUs), data = MPs_2013[which(MPs_2013$Location=="D"),]))

## Check 2014 to see if the patterns change with metals included
## But metals only measured in 2014
MPs_2014 <- MPs[-which(MPs$Year=="2013"),]

## Without metals in 2013, non-pesticides are not far off 1:1 line with total TUs
## and the correlation is strong
plot(log(TU_avg_nonPest_MS1) ~ log(TUs), data = MPs_2014[which(MPs_2014$Location=="D"),])
abline(0,1,col="blue")
summary(lm(log(TU_avg_nonPest_MS1) ~ log(TUs), data = MPs_2014[which(MPs_2014$Location=="D"),]))

## Without metals in 2013, pesticides are way off the 1:1 line with total TUs
## and the correlation is weaker than in 2013
plot(log(TU_avg_Pest_MS1) ~ log(TUs), data = MPs_2014[which(MPs_2014$Location=="D"),])
abline(0,1,col="blue")
summary(lm(log(TU_avg_Pest_MS1) ~ log(TUs), data = MPs_2014[which(MPs_2014$Location=="D"),]))

## Looking at MS4 (Metals) and MS5 (Others)
## Shows that metals make up the majority of TUs when present
plot(log(TU_avg_Metals_MS4) ~ log(TUs), data = MPs,
     xlim=c(-6,-2))
abline(0,1,col="blue")

## See that metals are nested with others + includes industrial chemicals
plot(log(TU_avg_Other_MS5) ~ log(TU_avg_Metals_MS4), data = MPs,
     ylim=c(-5,-2))
abline(0,1,col="blue")

## See that metals are nested with others + includes industrial chemicals
boxplot(log(TU_avg_Insecticide_MS5) ~ Year, data = MPs[(MPs$Location=="U1"),])
summary(aov(log(TU_avg_Insecticide_MS5) ~ Year, data = MPs[(MPs$Location=="U1"),]))

#############################################################################################
##                                                                                         ##
##                    RATIONALISE DATA AND WRITE FOR FURTHER ANALYSIS                      ##        
##                                                                                         ##    
#############################################################################################

##***************************************************************************************
## Assess data strings
colnames(MPs)
##  "Year"                   "Site"                   "Location"               "TUs"                   
## "TU_avg_nonPest_MS1"     "TU_avg_Pest_MS1"        "TU_avg_nonPest_MS2"     "TU_avg_Pest_MS2"       
## "TU_avg_other_MS3"       "TU_avg_Pesticide_MS3"   "TU_avg_Pharma_MS3"      "TU_avg_Metals_MS4"     
## "TU_avg_Other_MS4"       "TU_avg_Pesticides_MS4"  "TU_avg_Pharma_MS4"      "TU_avg_Fungicide_MS5"  
## "TU_avg_Herbicide_MS5"   "TU_avg_Insecticide_MS5" "TU_avg_Other_MS5"       "TU_avg_Pharma_MS5"    

## Aim to keep:
##  "Year"                   "Site"                   "Location"               "TUs"    
## "TU_avg_Metals_MS4"
## "TU_avg_Fungicide_MS5"  
## "TU_avg_Herbicide_MS5"   "TU_avg_Insecticide_MS5" "TU_avg_Other_MS5"       "TU_avg_Pharma_MS5" 

write_MPs <- MPs[,c(1:4,12,16:20)]

## Shows that metasubstances 5 (Fung,Herb,Insect,Other,Pharma) sum to total TUs
plot(log(rowSums(write_MPs[,c(6:10)])) ~ log(TUs), data = write_MPs)
abline(0,1,col="blue")

write.csv(write_MPs,"3_DATA_MP_TUs_metasubstances.csv",row.names = F)

#############################################################################################
##                                                                                         ##
##                            RECREATE SUMMARY DATA USED ABOVE                             ##        
##                                                                                         ##    
#############################################################################################

## Since initial data was calculated using Excel pivot-tables it would be optimal to have
## summary data generated with R scripting for reproducibilty and verification. Also want to
## calculate summary concentrations for supplementary tables

## Aim: Sum totals (concentrations and TUs) for metasubstance classes: 
## 1) Total TUs,
## 2) M4 (Metals), and 
## 3) M5 (Pesticides)

## Need to correct duplication of Val-de-Ruz site names
MP_metadata[which(MP_metadata$SiteName=="Val-De-Ruz"),6] <- "Val-de-Ruz"

## Note that Thiacloprid should be included as an Insecticide
## Incorrectly assigned to Herbicides
MP_metadata[which(MP_metadata$ChemName=="Thiacloprid"),23] <- "Insecticide"

## Note that Simazin-2-hydroxy should be included as an Herbicide
## Incorrectly assigned to Pharmaceuticals
MP_metadata[which(MP_metadata$ChemName=="Simazin-2-hydroxy"),23] <- "Herbicide"
	
#############################################################################################
##                                                                                         ##
##                                      TOTAL TUs                                          ##        
##                                                                                         ##    
#############################################################################################

## Sum totals
Total_TUs <- MP_metadata %>%
              group_by(SiteName,LMR,Year) %>% 
              summarise_at(vars("TU_Final"),sum,na.rm=T)

## Calculate means
Total_TUs_mu <- Total_TUs %>%
  group_by(SiteName,LMR) %>% 
  summarise_at(vars("TU_Final"),mean,na.rm=T)

## Check header
head(Total_TUs_mu)
glimpse(Total_TUs_mu)

#############################################################################################
##                                                                                         ##
##                                      METAL TUs                                          ##        
##                                                                                         ##    
#############################################################################################

## Sum totals
Meta4_TUs <- MP_metadata %>%
  group_by(SiteName,LMR,Year,Meta_substance_4) %>% 
  summarise_at(vars("TU_Final"),sum,na.rm=T)

## Calculate means
Meta4_TUs_mu <- Meta4_TUs %>%
  group_by(SiteName,LMR,Meta_substance_4) %>% 
  summarise_at(vars("TU_Final"),mean,na.rm=T)

## Separate Metal TUs
Meta4_TUs_long <- Meta4_TUs_mu[which(Meta4_TUs_mu$Meta_substance_4=="Metals"),]

## Check header
head(Meta4_TUs_long)
glimpse(Meta4_TUs_long)

#############################################################################################
##                                                                                         ##
##                              PESTICIDE TUs (AND OTHERS)                                 ##        
##                                                                                         ##    
#############################################################################################

## Sum totals
Meta5_TUs <- MP_metadata %>%
  group_by(SiteName,LMR,Year,MS5_Pesticides) %>% 
  summarise_at(vars("TU_Final"),sum,na.rm=T)

## Calculate means
Meta5_TUs_mu <- Meta5_TUs %>%
  group_by(SiteName,LMR,MS5_Pesticides) %>% 
  summarise_at(vars("TU_Final"),mean,na.rm=T)

## Separate Pesticides etc TUs
Meta5_TUs_wide <- Meta5_TUs_mu %>% spread(MS5_Pesticides,TU_Final)

## Check header
head(Meta5_TUs_wide)
glimpse(Meta5_TUs_wide)

#############################################################################################
##                                                                                         ##
##                          MERGE DATA AND REMOVE EXTRA SITES                              ##        
##                                                                                         ##    
#############################################################################################

## Note: Could also remove effluent, since it is not used in the analyses

## Create dataframe for all data
TU_data <- Total_TUs_mu

## Merge total TUs and Metal TUs
TU_data <- merge(Total_TUs_mu, Meta4_TUs_long[,-3], by=c("SiteName","LMR"), all = T)

## Check header
head(TU_data)
glimpse(TU_data)

## Rename as you go
names(TU_data)[3]<-"Total"
names(TU_data)[4]<-"Metal"

## Merge TUs and Pesticide TUs
TU_data <- merge(TU_data, Meta5_TUs_wide, by=c("SiteName","LMR"), all = T)

## Check header
head(TU_data)
glimpse(TU_data)

## Rename as you go
names(TU_data)[3]<-"Total"
names(TU_data)[4]<-"Metal"

## Remove effluent data
TU_data <- TU_data[-which(TU_data$LMR=="effluent"),]

## Add dummy variables for U2 based on U1
TU_data_U2_dummy <- TU_data[-which(TU_data$LMR=="downstream"),]
TU_data_U2_dummy[which(TU_data_U2_dummy$LMR=="upstream1"),2] <- "U2"

## Rename D and U1
TU_data[which(TU_data$LMR=="downstream"),2] <- "D"
TU_data[which(TU_data$LMR=="upstream1"),2] <- "U1"

## Combine data
TU_data_all <- rbind(TU_data,TU_data_U2_dummy)

## First, drop extra 2014 sites (i.e, not used)
TU_data_all <- TU_data_all[-c(which(TU_data_all$SiteName=="Aadorf"),
                              which(TU_data_all$SiteName=="Birmensdorf"),
                              which(TU_data_all$SiteName=="Unterehrendingen"),
                              which(TU_data_all$SiteName=="Zulwil")),]

## Create column for year
TU_data_all$Year <- NA

## Give each site the year it was sampled
TU_data_all[which(TU_data_all$SiteName=="Buttisholz"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Colombier"),10] <- "2013"
#TU_data_all[which(TU_data_all$SiteName=="D<fc>rnten-Bubikon"),10] "2013"
TU_data_all[which(TU_data_all$SiteName=="Elgg"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Ellikon"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Herisau"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Hochdorf"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Hornussen"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Kernenried"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Knonau"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Marthalen"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Messen"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Muri"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Niederdorf"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Reinach"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Romont"),10] <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Rothenturm"),10] <- "2013"
#TU_data_all[which(TU_data_all$SiteName=="S<fc>very"),1]0 <- "2013"
TU_data_all[which(TU_data_all$SiteName=="Val-de-Ruz"),10] <- "2014"
TU_data_all[which(TU_data_all$SiteName=="Villeret"),10] <- "2014"

## Rename sites by three letter codes
TU_data_all[which(TU_data_all$SiteName=="Buttisholz"),1] <- "BUT"
TU_data_all[which(TU_data_all$SiteName=="Colombier"),1] <- "COL"
#TU_data_all[which(TU_data_all$SiteName=="D<fc>rnten-Bubikon"),1] <- "DUR"
TU_data_all[which(TU_data_all$SiteName=="Elgg"),1] <- "ELG"
TU_data_all[which(TU_data_all$SiteName=="Ellikon"),1] <- "ELL"
TU_data_all[which(TU_data_all$SiteName=="Herisau"),1] <- "HER"
TU_data_all[which(TU_data_all$SiteName=="Hochdorf"),1] <- "HOC"
TU_data_all[which(TU_data_all$SiteName=="Hornussen"),1] <- "HOR"
TU_data_all[which(TU_data_all$SiteName=="Kernenried"),1] <- "KER"
TU_data_all[which(TU_data_all$SiteName=="Knonau"),1] <- "KNO"
TU_data_all[which(TU_data_all$SiteName=="Marthalen"),1] <- "MAR"
TU_data_all[which(TU_data_all$SiteName=="Messen"),1] <- "MES"
TU_data_all[which(TU_data_all$SiteName=="Muri"),1] <- "MUR"
TU_data_all[which(TU_data_all$SiteName=="Niederdorf"),1] <- "NIE"
TU_data_all[which(TU_data_all$SiteName=="Reinach"),1] <- "REI"
TU_data_all[which(TU_data_all$SiteName=="Romont"),1] <- "ROM"
TU_data_all[which(TU_data_all$SiteName=="Rothenturm"),1] <- "ROT"
#TU_data_all[which(TU_data_all$SiteName=="S<fc>very"),1] <- "SEV"
TU_data_all[which(TU_data_all$SiteName=="Val-de-Ruz"),1] <- "VAL"
TU_data_all[which(TU_data_all$SiteName=="Villeret"),1] <- "VIL"

## Rename and reorganise
TU_data_all <- TU_data_all[,c(10,1:9)]
head(TU_data_all)
names(TU_data_all)[2] <- "Site"
names(TU_data_all)[3] <- "Location"

## Arrange files
TU_data_all <-arrange(TU_data_all, Year, Site, Location)

## Rename DUR and SEV since above did not work
TU_data_all[c(55:60),1] <- c(rep("2013",6))
TU_data_all[c(55:57),2] <- c(rep("DUR",3))
TU_data_all[c(58:60),2] <- c(rep("SEV",3))

## Arrange files
TU_data_all <-arrange(TU_data_all, Year, Site, Location)

## Compare with original data
MPs_pivot <-read_csv("3_DATA_MP_TUs_metasubstances.csv")

## First remove unwanted sites from original data
MPs_pivot <- MPs_pivot[-c(which(MPs_pivot$Site=="AAD"),
                              which(MPs_pivot$Site=="BIR"),
                              which(MPs_pivot$Site=="UNT"),
                              which(MPs_pivot$Site=="ZUL")),]

## Check Total TUs
plot(log(TU_data_all$Total)~log(MPs_pivot$TUs))
abline(0,1,col="blue")

## Check Metal TUs
plot(log(TU_data_all$Metal)~log(MPs_pivot$TU_avg_Metals_MS4))
abline(0,1,col="blue")

## Check Fungicide TUs
plot(log(TU_data_all$Fungicide)~log(MPs_pivot$TU_avg_Fungicide_MS5))
abline(0,1,col="blue")

## Check Herbicide TUs
plot(log(TU_data_all$Herbicide)~log(MPs_pivot$TU_avg_Herbicide_MS5))
abline(0,1,col="blue")

## Check Insecticide TUs
plot(log(TU_data_all$Insecticide)~log(MPs_pivot$TU_avg_Insecticide_MS5))
abline(0,1,col="blue")

## Very minor deviation with Thiacloprid included (see where present)
## Muri	    downstream	14_03
## Ellikon	upstream1	  14_07
## Ellikon	downstream  14_07
## Reinach	downstream  14_05
resid(lm(log(TU_data_all$Insecticide)~log(MPs_pivot$TU_avg_Insecticide_MS5)))

## Check Other TUs
plot(log(TU_data_all$Other)~log(MPs_pivot$TU_avg_Other_MS5)) 
abline(0,1,col="blue")

## Check Pharmaceutical TUs
plot(log(TU_data_all$Pharmaceutical)~log(MPs_pivot$TU_avg_Pharma_MS5)) 
abline(0,1,col="blue")

## Write as output validating the data generated from the Excel worksheet
write.csv(TU_data_all,"3_DATA_MP_TUs_metasubstances_final.csv",row.names = F)

#############################################################################################
##                                                                                         ##
##                               MICROPOLLUTANT CONCENTRATIONS                             ##        
##                                                                                         ##    
#############################################################################################

#############################################################################################
##                                                                                         ##
##                                     TOTAL CONCENTRATIONS                                ##        
##                                                                                         ##    
#############################################################################################

## Sum totals
Total_conc <- MP_metadata %>%
  group_by(SiteName,LMR,Year) %>% 
  summarise_at(vars("MeasuredValue"),sum,na.rm=T)

## Calculate means
Total_conc_mu <- Total_conc %>%
  group_by(SiteName,LMR) %>% 
  summarise_at(vars("MeasuredValue"),mean,na.rm=T)

## Calculate maximum
Total_conc_max <- Total_conc %>%
  group_by(SiteName,LMR) %>% 
  summarise_at(vars("MeasuredValue"),max,na.rm=T)

## Check header
head(Total_conc_mu)
glimpse(Total_conc_mu)

## Check header
head(Total_conc_max)
glimpse(Total_conc_max)

#############################################################################################
##                                                                                         ##
##                                  METAL CONCENTRATIONS                                   ##        
##                                                                                         ##    
#############################################################################################

## Sum totals
Meta4_conc <- MP_metadata %>%
  group_by(SiteName,LMR,Year,Meta_substance_4) %>% 
  summarise_at(vars("MeasuredValue"),sum,na.rm=T)

## Separate Metal TUs
Meta4_conc_long <- Meta4_conc[which(Meta4_conc$Meta_substance_4=="Metals"),]

## Calculate means
Meta4_conc_mu <- Meta4_conc %>%
  group_by(SiteName,LMR,Meta_substance_4) %>% 
  summarise_at(vars("MeasuredValue"),mean,na.rm=T)

## Calculate maximum
Meta4_conc_max <- Meta4_conc %>%
  group_by(SiteName,LMR,Meta_substance_4) %>% 
  summarise_at(vars("MeasuredValue"),max,na.rm=T)

## Separate Metal TUs
Meta4_conc_mu_long <- Meta4_conc_mu[which(Meta4_conc_mu$Meta_substance_4=="Metals"),]
Meta4_conc_max_long <- Meta4_conc_max[which(Meta4_conc_max$Meta_substance_4=="Metals"),]

## Check
glimpse(Meta4_conc_mu_long)
glimpse(Meta4_conc_max_long)

#############################################################################################
##                                                                                         ##
##                          PESTICIDE CONCENTRATIONS (AND OTHERS)                          ##        
##                                                                                         ##    
#############################################################################################

## Sum totals
Meta5_conc <- MP_metadata %>%
  group_by(SiteName,LMR,Year,MS5_Pesticides) %>% 
  summarise_at(vars("MeasuredValue"),sum,na.rm=T)

## Separate Other TUs
Meta5_conc_long <- Meta5_conc[which(Meta5_conc$MS5_Pesticides=="Other"),]

## Calculate means
Meta5_conc_mu <- Meta5_conc %>%
  group_by(SiteName,LMR,MS5_Pesticides) %>% 
  summarise_at(vars("MeasuredValue"),mean,na.rm=T)

## Calculate maximum
Meta5_conc_max <- Meta5_conc %>%
  group_by(SiteName,LMR,MS5_Pesticides) %>% 
  summarise_at(vars("MeasuredValue"),max,na.rm=T)

## Separate Pesticides etc TUs
Meta5_conc_mu_wide <- Meta5_conc_mu %>% spread(MS5_Pesticides,MeasuredValue)
Meta5_conc_max_wide <- Meta5_conc_max %>% spread(MS5_Pesticides,MeasuredValue)

## Check
glimpse(Meta5_conc_mu_wide)
glimpse(Meta5_conc_max_wide)

#############################################################################################
##                                                                                         ##
##                               OTHER CONCENTRATIONS                                      ##        
##                                                                                         ##    
#############################################################################################

## Need to correct "Others" above for "Metals" to distinguish the two in 2014 
Others_all <- merge(Meta5_conc_long, Meta4_conc_long, by=c("SiteName","LMR","Year"), all = T)

## Turn NAs to zeros
Others_all[is.na(Others_all)] <- 0

## Subtract Metal TUs from Others and calculate means and maximum values
Others_all$MeasuredValue <- Others_all[,5]-Others_all[,7]

## Calculate means
Other_conc_mu <- Others_all %>%
  group_by(SiteName,LMR) %>% 
  summarise_at(vars("MeasuredValue"),mean,na.rm=T)

## Calculate maximum
Other_conc_max <- Others_all %>%
  group_by(SiteName,LMR) %>% 
  summarise_at(vars("MeasuredValue"),max,na.rm=T)

## Check
glimpse(Other_conc_mu)
glimpse(Other_conc_max)

#############################################################################################
##                                                                                         ##
##                      MERGE DATA AND REMOVE EXTRA SITES: MEANS                           ##        
##                                                                                         ##    
#############################################################################################

## Note: Could also remove effluent, since it is not used in the analyses

## Create dataframe for all data
Conc_data_mu <- Total_conc_mu

## Merge total TUs and Metal TUs
Conc_data_mu <- merge(Total_conc_mu, Meta4_conc_mu_long[,-3], by=c("SiteName","LMR"), all = T)

## Check header
head(Conc_data_mu)
glimpse(Conc_data_mu)

## Rename as you go
names(Conc_data_mu)[3]<-"Total"
names(Conc_data_mu)[4]<-"Metal"

## Merge TUs and Pesticide TUs
Conc_data_mu <- merge(Conc_data_mu, Meta5_conc_mu_wide[,-6], by=c("SiteName","LMR"), all = T)

## Rename as you go
names(Conc_data_mu)[3]<-"Total"
names(Conc_data_mu)[4]<-"Metal"

## Check header
head(Conc_data_mu)
glimpse(Conc_data_mu)

## Merge TUs and Other TUs
Conc_data_mu <- merge(Conc_data_mu, Other_conc_mu, by=c("SiteName","LMR"), all = T)

## Rename as you go
names(Conc_data_mu)[3]<-"Total"
names(Conc_data_mu)[4]<-"Metal"
names(Conc_data_mu)[9]<-"Other"

## Check header
head(Conc_data_mu)
glimpse(Conc_data_mu)

## Remove effluent data
#Conc_data_mu <- Conc_data_mu[-which(Conc_data_mu$LMR=="effluent"),]

# Rename D and U1
Conc_data_mu[which(Conc_data_mu$LMR=="downstream"),2] <- "D"
Conc_data_mu[which(Conc_data_mu$LMR=="upstream1"),2] <- "U1"

## Create dataframe
Conc_data_mu_all <- Conc_data_mu

## First, drop extra 2014 sites (i.e, not used)
Conc_data_mu_all <- Conc_data_mu_all[-c(which(Conc_data_mu_all$SiteName=="Aadorf"),
                              which(Conc_data_mu_all$SiteName=="Birmensdorf"),
                              which(Conc_data_mu_all$SiteName=="Unterehrendingen"),
                              which(Conc_data_mu_all$SiteName=="Zulwil")),]

## Create column for year
Conc_data_mu_all$Year <- NA

## Give each site the year it was sampled
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Buttisholz"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Colombier"),10] <- "2013"
#Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="D<fc>rnten-Bubikon"),10] "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Elgg"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Ellikon"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Herisau"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Hochdorf"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Hornussen"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Kernenried"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Knonau"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Marthalen"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Messen"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Muri"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Niederdorf"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Reinach"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Romont"),10] <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Rothenturm"),10] <- "2013"
#Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="S<fc>very"),1]0 <- "2013"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Val-de-Ruz"),10] <- "2014"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Villeret"),10] <- "2014"

## Rename sites by three letter codes
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Buttisholz"),1] <- "BUT"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Colombier"),1] <- "COL"
#Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="D<fc>rnten-Bubikon"),1] <- "DUR"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Elgg"),1] <- "ELG"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Ellikon"),1] <- "ELL"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Herisau"),1] <- "HER"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Hochdorf"),1] <- "HOC"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Hornussen"),1] <- "HOR"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Kernenried"),1] <- "KER"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Knonau"),1] <- "KNO"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Marthalen"),1] <- "MAR"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Messen"),1] <- "MES"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Muri"),1] <- "MUR"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Niederdorf"),1] <- "NIE"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Reinach"),1] <- "REI"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Romont"),1] <- "ROM"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Rothenturm"),1] <- "ROT"
#Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="S<fc>very"),1] <- "SEV"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Val-de-Ruz"),1] <- "VAL"
Conc_data_mu_all[which(Conc_data_mu_all$SiteName=="Villeret"),1] <- "VIL"

## Rename and reorganise
Conc_data_mu_all <- Conc_data_mu_all[,c(10,1:9)]
head(Conc_data_mu_all)
names(Conc_data_mu_all)[2] <- "Site"
names(Conc_data_mu_all)[3] <- "Location"

## Arrange files
Conc_data_mu_all <- arrange(Conc_data_mu_all, Year, Site, Location)

## Rename DUR and SEV since above did not work
Conc_data_mu_all[c(55:60),1] <- c(rep("2013",6))
Conc_data_mu_all[c(55:57),2] <- c(rep("DUR",3))
Conc_data_mu_all[c(58:60),2] <- c(rep("SEV",3))

## Arrange files
Conc_data_mu_all  <-arrange(Conc_data_mu_all, Year, Site, Location)

## For different years
Conc_data_mu_year <- Conc_data_mu_all %>%
                          group_by(Location,Year) %>% 
                          summarise_at(vars("Total","Metal","Fungicide","Herbicide",
                                            "Insecticide","Pharmaceutical","Other"),mean,na.rm=T)

Conc_data_mu_year <- as.data.frame(t(Conc_data_mu_year))

write.csv(Conc_data_mu_year,"Concentration_means_years_2013_2014.csv",row.names = T)

## For both years together
Conc_data_mu_Allyrs <- Conc_data_mu_all %>%
  group_by(Location) %>% 
  summarise_at(vars("Total","Metal","Fungicide","Herbicide",
                    "Insecticide","Pharmaceutical","Other"),mean,na.rm=T)

Conc_data_mu_year <- as.data.frame(t(Conc_data_mu_Allyrs))

write.csv(Conc_data_mu_year,"Concentration_means_All_years.csv",row.names = T)

#############################################################################################
##                                                                                         ##
##                      MERGE DATA AND REMOVE EXTRA SITES: MAX                             ##        
##                                                                                         ##    
#############################################################################################

## Note: Could also remove effluent, since it is not used in the analyses

## Create dataframe for all data
Conc_data_max <- Total_conc_max

## Merge total TUs and Metal TUs
Conc_data_max <- merge(Total_conc_max, Meta4_conc_max_long[,-3], by=c("SiteName","LMR"), all = T)

## Check header
head(Conc_data_max)
glimpse(Conc_data_max)

## Rename as you go
names(Conc_data_max)[3]<-"Total"
names(Conc_data_max)[4]<-"Metal"

## Merge TUs and Pesticide TUs
Conc_data_max <- merge(Conc_data_max, Meta5_conc_max_wide[,-6], by=c("SiteName","LMR"), all = T)

## Rename as you go
names(Conc_data_max)[3]<-"Total"
names(Conc_data_max)[4]<-"Metal"

## Check header
head(Conc_data_max)
glimpse(Conc_data_max)

## Merge TUs and Other TUs
Conc_data_max <- merge(Conc_data_max, Other_conc_max, by=c("SiteName","LMR"), all = T)

## Rename as you go
names(Conc_data_max)[3]<-"Total"
names(Conc_data_max)[4]<-"Metal"
names(Conc_data_max)[9]<-"Other"

## Check header
head(Conc_data_max)
glimpse(Conc_data_max)

## Check data work (i.e., means should add up to Total)
Check <- rowSums(Conc_data_max[,c(4:9)],na.rm = T)
plot(Conc_data_max$Total~Check)

## Remove effluent data
#Conc_data_max <- Conc_data_max[-which(Conc_data_max$LMR=="effluent"),]

# Rename D and U1
Conc_data_max[which(Conc_data_max$LMR=="downstream"),2] <- "D"
Conc_data_max[which(Conc_data_max$LMR=="upstream1"),2] <- "U1"

## Create dataframe
Conc_data_max_all <- Conc_data_max

## First, drop extra 2014 sites (i.e, not used)
Conc_data_max_all <- Conc_data_max_all[-c(which(Conc_data_max_all$SiteName=="Aadorf"),
                                        which(Conc_data_max_all$SiteName=="Birmensdorf"),
                                        which(Conc_data_max_all$SiteName=="Unterehrendingen"),
                                        which(Conc_data_max_all$SiteName=="Zulwil")),]

## Create column for year
Conc_data_max_all$Year <- NA

## Give each site the year it was sampled
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Buttisholz"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Colombier"),10] <- "2013"
#Conc_data_max_all[which(Conc_data_max_all$SiteName=="D<fc>rnten-Bubikon"),10] "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Elgg"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Ellikon"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Herisau"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Hochdorf"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Hornussen"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Kernenried"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Knonau"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Marthalen"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Messen"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Muri"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Niederdorf"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Reinach"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Romont"),10] <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Rothenturm"),10] <- "2013"
#Conc_data_max_all[which(Conc_data_max_all$SiteName=="S<fc>very"),1]0 <- "2013"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Val-de-Ruz"),10] <- "2014"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Villeret"),10] <- "2014"

## Rename sites by three letter codes
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Buttisholz"),1] <- "BUT"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Colombier"),1] <- "COL"
#Conc_data_max_all[which(Conc_data_max_all$SiteName=="D<fc>rnten-Bubikon"),1] <- "DUR"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Elgg"),1] <- "ELG"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Ellikon"),1] <- "ELL"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Herisau"),1] <- "HER"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Hochdorf"),1] <- "HOC"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Hornussen"),1] <- "HOR"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Kernenried"),1] <- "KER"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Knonau"),1] <- "KNO"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Marthalen"),1] <- "MAR"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Messen"),1] <- "MES"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Muri"),1] <- "MUR"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Niederdorf"),1] <- "NIE"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Reinach"),1] <- "REI"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Romont"),1] <- "ROM"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Rothenturm"),1] <- "ROT"
#Conc_data_max_all[which(Conc_data_max_all$SiteName=="S<fc>very"),1] <- "SEV"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Val-de-Ruz"),1] <- "VAL"
Conc_data_max_all[which(Conc_data_max_all$SiteName=="Villeret"),1] <- "VIL"

## Rename and reorganise
Conc_data_max_all <- Conc_data_max_all[,c(10,1:9)]
head(Conc_data_max_all)
names(Conc_data_max_all)[2] <- "Site"
names(Conc_data_max_all)[3] <- "Location"

## Arrange files
Conc_data_max_all <- arrange(Conc_data_max_all, Year, Site, Location)

## Rename DUR and SEV since above did not work
Conc_data_max_all[c(55:60),1] <- c(rep("2013",6))
Conc_data_max_all[c(55:57),2] <- c(rep("DUR",3))
Conc_data_max_all[c(58:60),2] <- c(rep("SEV",3))

## Arrange files
Conc_data_max_all  <-arrange(Conc_data_max_all, Year, Site, Location)

## For different years
Conc_data_max_year <- Conc_data_max_all %>%
  group_by(Location,Year) %>% 
  summarise_at(vars("Total","Metal","Fungicide","Herbicide",
                    "Insecticide","Pharmaceutical","Other"),max,na.rm=T)

Conc_data_max_year <- as.data.frame(t(Conc_data_max_year))

write.csv(Conc_data_max_year,"Concentration_maximum_years_2013_2014.csv",row.names = T)

## For both years together
Conc_data_max_Allyrs <- Conc_data_max_all %>%
  group_by(Location) %>% 
  summarise_at(vars("Total","Metal","Fungicide","Herbicide",
                    "Insecticide","Pharmaceutical","Other"),max,na.rm=T)

Conc_data_max_year <- as.data.frame(t(Conc_data_max_Allyrs))

write.csv(Conc_data_max_year,"Concentration_maximum_All_years.csv",row.names = T)

#############################################################################################
##                                                                                         ##
##                          COMPARE DATA WITH MUNZ et al. 2017                             ##        
##                                                                                         ##    
#############################################################################################

#load "pheatmap" library
library(pheatmap)

head(MP_Ecoimpact[,c(3,4,6,7,11,39,41,53)])

#calculate the pearson correlation coefficient matrix
myMat.cor <-cor(log1p(MP_Ecoimpact[,c(3,4,6,7,11,39,41,53)]), method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)

## Plot Burdon version SUM TUs Insecticides vs SPEAR
plot(log(MP_Ecoimpact$SPEAR) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5),
     ylab="SPEAR",
     xlab="log sum TU insecticides average")
abline(lm(log(MP_Ecoimpact$SPEAR) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5)),col="red")

## Plot Burdon version SUM TUs Insecticides vs SPEAR (but see zero values?)
plot(log(MP_Ecoimpact$SPEAR[-c(34,46,48)]) ~ log(MP_Ecoimpact$median_sumTU_I_gm[-c(34,46,48)]),
     ylab="SPEAR",
     xlab="log sum TU insecticides median")
abline(lm(log(MP_Ecoimpact$SPEAR) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5)),col="red")

## Plot SUM TUs Insecticides (but see zero values?)
plot(log(MP_Ecoimpact$median_sumTU_I_gm[-c(34,46,48)]) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5[-c(34,46,48)]),
     ylab="log sum TU insecticides median",
     xlab="log sum TU insecticides average")
abline(0,1,col="blue")
abline(lm(log(MP_Ecoimpact$median_sumTU_I_gm[-c(34,46,48)]) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5[-c(34,46,48)])),col="red")

## Plot SUM TUs Insecticides median? (but see zero values?)
plot(log(MP_Ecoimpact$median_sumTU_I_med[-c(34,46,48)]) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5[-c(34,46,48)]),
     ylab="log sum TU insecticides median",
     xlab="log sum TU insecticides average")
abline(0,1,col="blue")
abline(lm(log(MP_Ecoimpact$median_sumTU_I_med[-c(34,46,48)]) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5[-c(34,46,48)])),col="red")

## Plot Burdon version SUM TUs Insecticides vs SPEAR
plot(log(SPEAR) ~ log(TU_avg_Insecticide_MS5),
     data=MP_Ecoimpact[which(MP_Ecoimpact$Location=="D"),],
     ylab="SPEAR",
     xlab="log sum TU insecticides average")
abline(lm(log(MP_Ecoimpact$SPEAR) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5)),col="red")

## Plot Burdon version SUM TUs Insecticides vs SPEAR
plot(log(SPEAR) ~ log(TU_avg_Insecticide_MS5),
     data=MP_Ecoimpact[which(MP_Ecoimpact$Location=="U1"),],
     ylab="SPEAR",
     xlab="log sum TU insecticides average")
abline(lm(log(MP_Ecoimpact$SPEAR) ~ log(MP_Ecoimpact$TU_avg_Insecticide_MS5)),col="red")

## Plot Burdon version SUM TUs Insecticides vs SPEAR
plot(log(SPEAR) ~ log(TU_avg_Insecticide_MS5),
     data=MP_Ecoimpact[which(MP_Ecoimpact$Location=="U1"),],
     ylab="SPEAR",
     xlab="log sum TU insecticides average")

## Combine input datasets to assess differences between years
MP_compare <- merge(MPs[,c(1:3,18)],MP_Ecoimpact[,c(3,11,37,38)],by=c("Site","Location"))
# Create new column filled with default colour
MP_compare$Colour="black"
# Set new column values to appropriate colours
MP_compare$Colour[MP_compare$Year=="2013"] <- "forestgreen"
MP_compare$Colour[MP_compare$Year=="2014"] <- "royalblue"

## Take U1 only
MP_compare_U1 <- MP_compare[which(MP_compare$Location=="U1"),]

## Plot Burdon version SUM TUs Insecticides vs SPEAR
## No consistent trends in TUs between years, although interesting bifucation in values across years 
plot(SPEAR ~ log(TU_avg_Insecticide_MS5),
     data= MP_compare_U1,
     col= MP_compare_U1$Colour,
     pch= 16, cex=2,
     ylab = "SPEAR",
     xlab = "log sum TU insecticides average")
abline(lm(MP_compare_U1 $SPEAR ~ log(MP_compare_U1$TU_avg_Insecticide_MS5)),col="red")

#############################################################################################
##                                                                                         ##
##                          COMPARE DATA WITH MUNZ et al. 2017 VER02                       ##
##                              &  BURDON et al. 2019 STOTEN                               ##
##                                                                                         ##    
#############################################################################################

## Repeats the analysis above
#load "pheatmap" library
library(pheatmap)

head(MP_Ecoimpact[,c(3,4,6,7,11,39,41,53)])

#calculate the pearson correlation coefficient matrix
myMat.cor <-cor(log1p(MP_Ecoimpact[,c(3,4,6,7,11,39,41,53)]), method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)

## Rationalise input data from Burdon et al. 2019 to include only TUs data 
MP_Ecoimpact_invert_TUs <- MP_Ecoimpact_invert[,c(1:4,26:39)]

## Take data generated from raw data (Metasubstances) - see updates for Thiacloprid
TU_data_all_mU2 <- TU_data_all[-which(TU_data_all$Location=="U2"),]

## Merge data sources for comparison
TU_data_Ecoimpact  <- merge(TU_data_all_mU2,MP_Ecoimpact_invert_TUs,by=c("Year","Site","Location"))

## Geometric mean sum Insecticides vs median sum Pesticides and Biocides
min(TU_data_Ecoimpact$sumTU_PB)
min(TU_data_Ecoimpact$Insecticide)
plot(log(Insecticide) ~ log(sumTU_PB), data=TU_data_Ecoimpact)
abline(lm(log(Insecticide) ~ log(sumTU_PB), data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(log(Insecticide) ~ log(sumTU_PB), data=TU_data_Ecoimpact))

## Geometric mean sum Insecticides vs "log_TU_PB_U1"??
min(TU_data_Ecoimpact$log_TU_PB_U1)
min(TU_data_Ecoimpact$Insecticide)
plot(log(Insecticide) ~ log_TU_PB_U1, data=TU_data_Ecoimpact)
abline(lm(log(Insecticide) ~ log_TU_PB_U1, data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(log(Insecticide) ~ log_TU_PB_U1, data=TU_data_Ecoimpact))

## Geometric mean sum Insecticides vs max sum Pesticides and Biocides
min(TU_data_Ecoimpact$maxTU_PB)
min(TU_data_Ecoimpact$Insecticide)
plot(log(Insecticide) ~ log(maxTU_PB), data=TU_data_Ecoimpact)
abline(lm(log(Insecticide) ~ log(maxTU_PB), data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(log(Insecticide) ~ log(maxTU_PB), data=TU_data_Ecoimpact))

## Geometric mean sum Insecticides vs "log_maxTU_PB"??
min(TU_data_Ecoimpact$log_maxTU_PB)
min(TU_data_Ecoimpact$Insecticide)
plot(log(Insecticide) ~ log_maxTU_PB, data=TU_data_Ecoimpact)
abline(lm(log(Insecticide) ~ log_maxTU_PB, data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(log(Insecticide) ~ log_maxTU_PB, data=TU_data_Ecoimpact))

## Now assess the data from Munz et al. 2017
head(MP_Ecoimpact_sumTU_msPAF)
#write.csv(MP_Ecoimpact_sumTU_msPAF,"MP_Ecoimpact_sumTU_msPAF.csv",row.names = F)

## Reload data with correct site info
MP_Ecoimpact_sumTU_msPAF <- read_csv("MP_Ecoimpact_sumTU_msPAF.csv")
glimpse(MP_Ecoimpact_sumTU_msPAF)

## Subset only Pesticides-Biocides and Insecticides
MP_Ecoimpact_sumTU_msPAF_PB_Insecticides <- MP_Ecoimpact_sumTU_msPAF[,c(1:4,8,12,15,20,23,28,31,36,39,44)] 

## Merge data sources for comparison
TU_data_Ecoimpact  <- merge(TU_data_all_mU2,MP_Ecoimpact_sumTU_msPAF_PB_Insecticides,by=c("Year","Site","Location"))

## Geometric mean sum Insecticides vs "median_msPAFraMix_PB" Pesticides and Biocides
min(TU_data_Ecoimpact$median_msPAFraMix_PB)
min(TU_data_Ecoimpact$Insecticide)
plot(log(Insecticide) ~ log(median_msPAFraMix_PB), data=TU_data_Ecoimpact)
abline(lm(log(Insecticide) ~ log(median_msPAFraMix_PB), data=TU_data_Ecoimpact),col="red")
summary(lm(log(Insecticide) ~ log(median_msPAFraMix_PB), data=TU_data_Ecoimpact))

## Geometric mean sum Insecticides vs "median_msPAFraMix_I" Insecticides
min(TU_data_Ecoimpact$median_msPAFraMix_I)
min(TU_data_Ecoimpact$Insecticide)
plot(sqrt(Insecticide) ~ sqrt(median_msPAFraMix_I), data=TU_data_Ecoimpact)
abline(lm(sqrt(Insecticide) ~ sqrt(median_msPAFraMix_I), data=TU_data_Ecoimpact),col="red")
summary(lm(sqrt(Insecticide) ~ sqrt(median_msPAFraMix_I), data=TU_data_Ecoimpact))

## Create plot for interpretation
png(filename="ECOIMPACT_Insecticides_TUs_D_magna_updated.png", 
    type="cairo",
    units="in", 
    width=8.75, 
    height=7, 
    pointsize=14, 
    res=600)

par(pty="s")

## Geometric mean sum Insecticides vs Geometric mean sum Insecticides?
min(TU_data_Ecoimpact$median_sumTU_I.DM_gm)
min(TU_data_Ecoimpact$Insecticide)
plot(sqrt(Insecticide) ~ sqrt(median_sumTU_I.DM_gm), 
     xlab="sqrt median_sumTU_I.DM_gm (Munz et al. 2017)",
     ylab="sqrt GM sum TUs Insecticides (Burdon et al. in prep)",
     main="Insecticides TUs Daphnia magna",
     data=TU_data_Ecoimpact)
abline(lm(sqrt(Insecticide) ~ sqrt(median_sumTU_I.DM_gm), data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(sqrt(Insecticide) ~ sqrt(median_sumTU_I.DM_gm), data=TU_data_Ecoimpact))

dev.off()

## See that my estimates of TUs never drops below 1:1 line, so could be due to differences 
## in Mode of Action definitions - but that appears not to be the case (I have checked)

## The difference could lie in how the data are used i.e. when does the median get applied to the GM?

## Note the additional values I added on the advice of Andreas Focks
## This is unlikely to be the source of the difference - none of the compounds are Insecticides
## Also Nicole likely used the updated version in her calculations

## See "Notes" tab in "DATA_TU_metasubstance_VER01_161111.xlsx"
## Chlorothalonil-4-hydroxy-Carbonsäureamid (TP611968)
## Irgarol-descyclopropyl
## Metazachlor-OXA
## Simazin-2-hydroxy
## Entered manually from Andreas data

## Additional notes for the above compounds
## Chlorothalonil-4-hydroxy-Carbonsäureamid: TP of Chlorothalonil; fungicide/biocide)
## Irgarol-descyclopropyl: TP of Irgarol and Terbutryn; Anti-fouling and Herbicide
## Metazachlor-OXA: TP of Metazachlor; herbicide
## Simazin-2-hydroxy: TP of Simazine and Simeton; Herbicides

## Create plot for interpretation
png(filename="ECOIMPACT_Insecticides_TUs_updated.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=7, 
    pointsize=14, 
    res=600)

par(mfrow=c(1,2)) 

## Geometric mean sum Insecticides vs Geometric mean sum Insecticides?
min(TU_data_Ecoimpact$median_sumTU_I_gm)
min(TU_data_Ecoimpact$Insecticide)
plot(sqrt(Insecticide) ~ sqrt(median_sumTU_I_gm), 
     xlab="sqrt median_sumTU_I_gm (Munz et al. 2017)",
     ylab="sqrt GM sum TUs Insecticides (Burdon et al. in prep)",
     main="Insecticides TUs (All - Munz et al. 2017)",
     data=TU_data_Ecoimpact)
abline(lm(sqrt(Insecticide) ~ sqrt(median_sumTU_I_gm), data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(sqrt(Insecticide) ~ sqrt(median_sumTU_I_gm), data=TU_data_Ecoimpact))

## Geometric mean sum Insecticides vs Geometric mean sum Insecticides?
min(TU_data_Ecoimpact$median_sumTU_I.DM_gm)
min(TU_data_Ecoimpact$Insecticide)
plot(sqrt(Insecticide) ~ sqrt(median_sumTU_I.DM_gm), 
     xlab="sqrt median_sumTU_I.DM_gm (Munz et al. 2017)",
     ylab="sqrt GM sum TUs Insecticides (Burdon et al. in prep)",
     main="Insecticides TUs Daphnia magna",
     data=TU_data_Ecoimpact)
abline(lm(sqrt(Insecticide) ~ sqrt(median_sumTU_I.DM_gm), data=TU_data_Ecoimpact),col="red")
abline(0,1,col="blue")
summary(lm(sqrt(Insecticide) ~ sqrt(median_sumTU_I.DM_gm), data=TU_data_Ecoimpact))

dev.off()

## Create plot for interpretation
png(filename="ECOIMPACT_Insecticides_TUs_SPEAR_updated.png", 
    type="cairo",
    units="in", 
    width=8.75, 
    height=7, 
    pointsize=14, 
    res=600)

#par(mfrow=c(1,2), mar = c(4,4,2,2) + 0.1) 

## Geometric mean sum Insecticides vs max sum Pesticides and Biocides
min(TU_data_Ecoimpact$SPEAR)
min(TU_data_Ecoimpact$Insecticide)
plot(log(Insecticide) ~ log(SPEAR), 
     xlab="log SPEAR (Munz et al. 2017)",
     ylab="log GM sum TUs Insecticides (Burdon et al. in prep)",
     main="Insecticides TUs Daphnia magna",
     data=TU_data_Ecoimpact)
abline(lm(log(Insecticide) ~ log(SPEAR), data=TU_data_Ecoimpact),col="red")
summary(lm(log(Insecticide) ~ log(SPEAR), data=TU_data_Ecoimpact))

dev.off()