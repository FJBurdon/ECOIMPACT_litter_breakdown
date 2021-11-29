############################################################################################
##                                                                                        ##  
##                                MICROPOLLUTANTS (MPs)                                   ##  
##                                                                                        ##
############################################################################################

##  Script: ECOIMPACT leafpack assay data preparation
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Combine three data files into one to be used in the SEM analyses
##  2. Remove unnecessary columns and write file for further analysis
##  3. Compare data for further quality assurance

## NOTE: There is a discrepancy in the "Tab_Delta_allMacro_v20181116" file - ELL and ELG mixed up?

## ADDITIONAL NOTES: From README file in folder 

##********************************************************
## ECOIMPACT Toxcicity data - MPs converted to Toxic Units
## README - Dr. Francis J. Burdon

## 1. A .csv file including max values for sum TUs for HM, PB, and Other TUs has been subsetted from the 
## data used in the second macroinvertebrate paper (Burdon et al. 2017). This is based on based on 
## tox data from all invertebrates, using the geometric mean for critical values. Values for U1 are 
## used as dummy variables for U2 - this is for data exploration purposes and not final analyses (use 
## only D and U1 locations):
  
##  - Output: MP_Ecoimpact_invert_TU_final.csv
##  -	Input: Tab_Delta_allMacro_v20181116.csv

## 2. A .csv file has similarly been generated from the data file containing median values for sum TUs 
## for HM, PB, and Other (also Insecticides, Herbicides, and Fungicides). The max values for the three 
## broad groups (HM, PB, and Other) are provided above. This data includes tox data from all invertebrates, 
## and Dapnia magna (DM) using the geometric mean and logarithmic mean. These data correspond to the data 
## provided by Nicole Munz. Values for U1 are used as dummy variables for U2 - this is for data exploration 
## purposes and not final analyses (use only D and U1 locations):
  
##  -	Output: MP_Ecoimpact_sumTU_msPAF_final.csv
##  -	Input: Spear_median_sumTU_msPAF.csv 

## The above input data file can be found the switchdrive (switchdrive\EcoImpact_Macro\Output Data\Overview_sumTU_msPAF\spear_median_sumTU_msPAF.csv).

## 3. A .csv file has similarly been generated from the raw data file used by Nicole Munz to generate the 
## above data using critical toxicity values for Dapnia magna. This data file includes mean, median, and 
## max values for sum TUs: total, HM, Insecticides, Herbicides, Fungicides, Pharmaceuticals, and Other. 
## This data includes tox data from Dapnia magna (DM) using the geometric mean. These data correspond to 
## the data provided by Nicole Munz. Values for U1 are used as dummy variables for U2 - this is for data 
## exploration purposes and not final analyses (use only D and U1 locations):
  
##  -	Output: 3_DATA_MP_TUs_metasubstances_final.csv
##  -	Input: 3_DATA_MPs_TU_metasubstance_VER02_161114.csv 

#The above output data (and output from 1 and 2) are generated in the R script "Rscript_MPs_data_preparation_201124.R" 

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

## Load data prepared by FJB for Burdon et al. in prep
MP_TUs_in_prep <- read_csv("3_DATA_MP_TUs_metasubstances_final.csv")

## Load data used in Munz et al. 2017 Water Research
MP_TUs_msPAF <- read_csv("MP_Ecoimpact_sumTU_msPAF_final.csv")

## Load data used in Burdon et al. 2019 STOTEN (based on Munz et al. 2017)
MP_TUs_Macroinvert <- read_csv("MP_Ecoimpact_invert_TU_final.csv")
MP_TUs_Macroinvert <- MP_TUs_Macroinvert [,-1]

#############################################################################################
##                                                                                         ##
##                            RATIONALISE AND COLLATE DATA                                 ##        
##                                                                                         ##    
#############################################################################################

glimpse(MP_TUs_in_prep)

## Subset data in "MP_TUs_msPAF" includes mean and median data for summed TUs. The TU values
## are based on "All" invertebrates (as oppposed to Daphnia magna in the above data). Values for 
## D. magna are included for the broad groups used in Burdon et al. 2019 STOTEN (HM, PB, Other, Org, and Total)
MP_TUs_msPAF_out <- MP_TUs_msPAF[,c(1:3,13:20,29:33)]

## Subset data in "MP_TUs_Macroinvert" focuses on max values for summed TUs. The TU values
## are based on "All" invertebrates
MP_TUs_Macroinvert_out <- MP_TUs_Macroinvert[,c(2:4,12:14)]

## First, merge selected data from "MP_TUs_msPAF"
MP_TUs_final <- merge(MP_TUs_in_prep,MP_TUs_msPAF_out,by=c("Year","Site","Location"))

## Second, merge selected data from "MP_TUs_Macroinvert"
MP_TUs_final <- merge(MP_TUs_final,MP_TUs_Macroinvert_out,by=c("Year","Site","Location"))

## Write output for use in SEM
write.csv(MP_TUs_final,"3_DATA_MP_TUs_final.csv",row.names = F)

#############################################################################################
##                                                                                         ##
##                             CONFIRM KEY RELATIONSHIPS                                   ##        
##                                                                                         ##    
#############################################################################################

## PB: Pesticides and Biocides
## Data in Burdon et al. 2019 relies on the geometric mean of critical toxicity values
## Then uses the median values (but recall reviewer demand to use max values)
plot(sqrt(MP_TUs_Macroinvert$sumTU_PB) ~ sqrt(MP_TUs_msPAF$median_sumTU_PB_gm))
abline(0,1,col="blue")

## Check deviant values
which.max(resid(lm(sqrt(MP_TUs_Macroinvert$sumTU_PB) ~ sqrt(MP_TUs_msPAF$median_sumTU_PB_gm))))
which.min(resid(lm(sqrt(MP_TUs_Macroinvert$sumTU_PB) ~ sqrt(MP_TUs_msPAF$median_sumTU_PB_gm))))

## Could these values be mixed up?
MP_TUs_Macroinvert[43,c(2:4)]
MP_TUs_msPAF[43,c(1:3)]

## Could these values be mixed up?
MP_TUs_Macroinvert[46,c(2:4)]
MP_TUs_msPAF[46,c(1:3)]

## This looks wrong - based off the raw data
MP_TUs_Macroinvert$sumTU_PB[46]-MP_TUs_Macroinvert$sumTU_PB[43]
## ELL D < ELG D -0.00010559
MP_TUs_Macroinvert$sumTU_PB[47]-MP_TUs_Macroinvert$sumTU_PB[44]
## ELL U1 < ELG U1 -9e-07

## This looks correct - based off the raw data - Munz et al. 2020
MP_TUs_msPAF$median_sumTU_PB_gm[46]-MP_TUs_msPAF$median_sumTU_PB_gm[43]
## ELL D > ELG D 0.00010559
MP_TUs_msPAF$median_sumTU_PB_gm[47]-MP_TUs_msPAF$median_sumTU_PB_gm[44]
## ELL U1 > ELG U1 9e-07

## PB: Pesticides and Biocides
## Data in Burdon et al. 2019 relies on the geometric mean of critical toxicity values
## Then uses the median values (but recall reviewer demand to use max values)
## Check max values to assess if those data are also mixed up in the Burdon et al. 2019 data
## Note entirely clear, but appear to be OK - no red flags in the plot
plot(sqrt(MP_TUs_Macroinvert$maxTU_PB) ~ sqrt(MP_TUs_msPAF$median_sumTU_PB_gm))
abline(0,1,col="blue")

MP_TUs_Macroinvert$maxTU_PB[43]-MP_TUs_msPAF$median_sumTU_PB_gm[43]
## Max > Med 0.000298141
MP_TUs_Macroinvert$maxTU_PB[44]-MP_TUs_msPAF$median_sumTU_PB_gm[44]
## Max > Med 0.000116385

MP_TUs_Macroinvert$maxTU_PB[46]-MP_TUs_msPAF$median_sumTU_PB_gm[46]
## Max > Med 0.000410339
MP_TUs_Macroinvert$maxTU_PB[47]-MP_TUs_msPAF$median_sumTU_PB_gm[47]
## Max > Med 0.000173255

## Other MPs
## Data in Burdon et al. 2019 relies on the geometric mean of critical toxicity values
## Then uses the median values (but recall reviewer demand to use max values)
plot(sqrt(MP_TUs_Macroinvert$sumTU_Oth) ~ sqrt(MP_TUs_msPAF$median_sumTU_Oth_gm))
abline(0,1,col="blue")

## Check deviant values
which.max(resid(lm(sqrt(MP_TUs_Macroinvert$sumTU_Oth) ~ sqrt(MP_TUs_msPAF$median_sumTU_Oth_gm))))
which.min(resid(lm(sqrt(MP_TUs_Macroinvert$sumTU_Oth) ~ sqrt(MP_TUs_msPAF$median_sumTU_Oth_gm))))

## Could these values be mixed up? Yes, see above
MP_TUs_Macroinvert[43,c(2:4)]
MP_TUs_msPAF[43,c(1:3)]

## Could these values be mixed up? Yes, see above
MP_TUs_Macroinvert[46,c(2:4)]
MP_TUs_msPAF[46,c(1:3)]


## Compare Insecticides
## Data in Burdon et al. in prep relies on the geometric mean of critical toxicity values
## Then uses the mean values
## First create data for comparison
MP_TUs_temp <- merge(MP_TUs_in_prep[,c(1:3,17)],MP_TUs_msPAF[,c(1:3,36)],by=c("Year","Site","Location"))

## This indicates the data used in Burdon et al. 2019 could be problematic if ELG and ELL are mixed up
## But helps confirm that the data generated for Burdon et al. in prep is valid
plot(sqrt(MP_TUs_temp$Insect_TU_med) ~ sqrt(MP_TUs_temp$median_sumTU_I.DM_gm))
abline(0,1,col="blue")
