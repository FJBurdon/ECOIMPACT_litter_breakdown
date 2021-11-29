#############################################################################################
##                                                                                         ##  
##                            STRUCTURAL EQUATION MODELLING                                ##   
##                                                                                         ##  
#############################################################################################

##  Script: ECOIMPACT leafpack assay data analysis
##  Author: Dr. Francis J. Burdon
##  VERSION: 3.0
##
##  Primary objective: Use structural equation modelling to explore the mechanisms driving 
##  differences in alder litter breakdown at sites. 
##
##  Version objective: using MEAN summed TU values for Daphnia magna
##  (see Burdon et al. 2019  STOTEN and Munz et al. 2017 for more information)
##
##  Specific objectives: 
##  1. Use data from all locations to assess potential predictors for further exploration in 
##     SEM (note that TUs from U1 are used as a dummy variable for U2)
##  2. Use data from all years at locations D and U1 to test the influence of key environmental 
##     predictors on invertebrate shredders and litter breakdown
##  3. Use data from 2014 (locations D and U1) to test the influence of key environmental predictors 
##     including Metal TUs on invertebrate shredders and litter breakdown
##  4. Create partial regression plots for main relationships elucidated in the SEM to help show the data
##     underpinning the main result (Figure 2, Main Text)   

##***************************************************************************************
## Load libraries you need

# Load piecewiseSEM from CRAN
require(piecewiseSEM) # Version 2.1.0
require(semEff) # Not available for piecewise SEM
require(packfor) # for forward-selection
require(qpcR) # to calculate AIC weights
require(lme4) # for forward-selection
require(nlme) # to calculate AIC weights

##Loads package used for Logit transformation (% data)
require(car)

##Loads vegan package for transformation
require(vegan)

##Loads reshape and plyr package for data manipulation
require(reshape)
require(reshape2)
require(plyr)

##Visualise colinearity
require(pheatmap)

## Forward Selection Blanchet et al. (2008)
require(packfor)

## Load data files
require(readr)

#############################################################################################
##                                                                                         ##  
##                                  LOAD DATA                                              ##
##                                                                                         ##
#############################################################################################

## Community - invertebrates
gammarid <- read.csv(here::here("1_Input_data","1_lpa_alder_coarse_allyrs_invert_abundances_data.csv")) 
invert_traits <- read.csv(here::here("1_Input_data","2_lpa_alder_coarse_allyrs_CWM_traits_data_all_leapacks.csv")) 

## Environmental
MPs <- read.csv(here::here("1_Input_data","3_DATA_MP_TUs_final.csv")) 
Waterchem <- read.csv(here::here("1_Input_data","4_DATA_general_water_chemistry_160824.csv")) 
Landuse <- read.csv(here::here("1_Input_data","5_DATA_Landuse_201116.csv")) 
Wastewater <- read.csv(here::here("1_Input_data","6_DATA_wastwater_dilution_factors_161010.csv")) 
Habitat <- read.csv(here::here("1_Input_data","7_DATA_habitat_160824.csv")) 

## Ecosystem
decomposition <- read.csv(here::here("1_Input_data","8_lpa_alder_coarse_All_data.csv")) 

#############################################################################################
##                                                                                         ##  
##                                  DECOMPOSITION                                          ##   
##                                                                                         ##  
#############################################################################################

## Sort data for coarse alder leafpacks
decomp_alder_coarse<-decomposition[which(decomposition$Leaf == "alder" & 
                                           decomposition$Mesh == "coarse"),] 

## Get mean values for further analysis
Mean_kdd_all_years<- ddply(decomp_alder_coarse, .(Site_code, Location, Year), summarise, 
                           mean = mean(k_dd_ml))

## Rename temperature corrected breakdown
names(Mean_kdd_all_years)[4]<-"Kdd_alder"

## Create new data object
Data_all <- Mean_kdd_all_years

## Sort data by Year, Site and Location
Data_all <- arrange(Data_all, Year, Site_code, Location)

## Tidy up workspace
remove(decomp_alder_coarse, Mean_kdd_all_years)

## Check data for normality (mean values)
shapiro.test(Data_all$Kdd_alder)
qqnorm(Data_all$Kdd_alder, pch = 1, frame = FALSE)
qqline(Data_all$Kdd_alder, col = "steelblue", lwd = 2) ## Need log-transforming

## Relabel the location
Data_all$Temp <- NA
Data_all[which(Data_all$Location=="DS"),5]   <- "D"
Data_all[which(Data_all$Location=="US1"),5]  <- "U1"
Data_all[which(Data_all$Location=="US2"),5]  <- "U2"

## Reorganise and relabel
Data_all <- Data_all[,c(1,5,3:4)]
names(Data_all)[2] <- "Location"

## Check output
head(Data_all)

#############################################################################################
##                                                                                         ##
##                                  INVERTEBRATES                                          ##  
##                                                                                         ##
#############################################################################################

head(gammarid[,c(1:3,6,16)])
head(invert_traits[,c(1:4,57,66)])

invertebrates <- merge(gammarid[,c(1:3,6,16)],invert_traits[,c(1:4,57,66)],
                       by=c("Site_code","Location","Rep","Year"))

head(invertebrates)

## Difference is because this data include zeros (no invertebrates found in leafpack)
plot(log1p(Food_3)~log1p(Gammaridae), data=invertebrates)
summary(lm(log1p(Food_3)~log1p(Gammaridae), data=invertebrates))
plot(log1p(Feeding_3)~log1p(Gammaridae), data=invertebrates)
summary(lm(log1p(Feeding_3)~log1p(Gammaridae), data=invertebrates))

require(tidyverse)
## Create mean values for each location
inverts_mu <- as.data.frame(invertebrates %>%
                group_by(Site_code,Location,Year) %>% 
                summarise_at(vars("Gammaridae","Food_3","Feeding_3"), 
                mean, na.rm=T))

## Merge data
Data_all <- merge(Data_all, inverts_mu, by=c("Site_code","Location","Year"))

##Test for normality
shapiro.test(Data_all$Gammaridae)
qqnorm(Data_all$Gammaridae, pch = 1, frame = FALSE)
qqline(Data_all$Gammaridae, col = "steelblue", lwd = 2) ## Need log+1 transforming

#############################################################################################
##                                                                                         ##  
##                                  ENVIRONMENTAL                                          ##
##                                                                                         ##  
#############################################################################################

## Arrange files (but depracated, because I will use the merge function)
MPs<-arrange(MPs, Year, Site, Location)
Waterchem<-arrange(Waterchem, Year, Site, Location)
Landuse<-arrange(Landuse, Year, Site, Location)
Wastewater<-arrange(Wastewater, Year, Site, Location)

## MICROPOLLUTANTS

## Calculate non-Insecticide TUs
MPs <- as.data.frame(MPs)

## Convert metals to zeros (NAs)
MPs[is.na(MPs)] <- 0

## Check column names for different groups of values
colnames(MPs)

## "Year"         "Site"          "Location"               
## "Total_TU_mu"  "Total_TU_med"  "Total_TU_max"           
## "Metal_TU_mu"  "Metal_TU_med"  "Metal_TU_max"           
## "Fung_TU_mu"   "Herb_TU_mu"    "Insect_TU_mu"   "Other_TU_mu"  "Pharma_TU_mu"          

## Exclude metals from non-insecticide TUs
MPs$MP_nonInsecticide <- NA
MPs$MP_nonInsecticide <- MPs[,4]-(MPs[,12]+MPs[,7])

## Change Data_all
names(Data_all)[1] <- "Site"

## Merge with data and rename columns (Drop TUs Other since correlated with TUs Metals)
Data_all <- merge(Data_all, MPs[,c(1:3,4,7,10:14,41)], by=c("Site","Location","Year"))
colnames(Data_all) <- c("Site","Location","Year","Kdd_alder","Gammaridae","Food_3","Feeding_3",
                "Total_TUs","Metal_TUs","Fungicide_TUs", 
                "Herbicide_TUs","Insecticide_TUs","Other_TUs",
                "Pharmaceutical_TUs","NonInsecticide_TUs")
## NUTRIENTS
                
## Calculate DIN
Waterchem <- as.data.frame(Waterchem)
Waterchem$DIN <- Waterchem[,c(12)]/1000+Waterchem[,c(13)]/1000+Waterchem[,c(14)]

## Assess correlation of DIN with NH4
plot(log(NH4) ~ log(DIN), data = Waterchem) ## Correlated
abline(lm(log(NH4) ~ log(DIN), data = Waterchem),col="red")
summary(lm(log(NH4) ~ log(DIN), data = Waterchem))
head(Waterchem)

## Merge with data and rename columns (if needed)
Data_all <- merge(Data_all, Waterchem[,c(1:3,16,24,12)], by=c("Site","Location","Year"))

## WASTEWATER

## Calculate % of WW in stream
Wastewater <- as.data.frame(Wastewater)
Wastewater$WW_fraction <- Wastewater[,5]/Wastewater[,4]*100

## Merge with data and rename columns (if needed)
Data_all <- merge(Data_all, Wastewater[,c(1:3,10)], by=c("Site","Location","Year"))

## LANDUSE

## Turn tibble in data frame
Landuse <- as.data.frame(Landuse)

## Check input data
head(Landuse[,c(2:4,6)])

## Merge with data and rename columns (if needed)
Data_all <- merge(Data_all, Landuse[,c(2:4,6)], by=c("Site","Location","Year"))

names(Data_all)[20] <- "Cropping"

## HABITAT (SEDIMENT)

## Turn tibble in data frame
Habitat <- as.data.frame(Habitat)

## Merge with data and rename columns (if needed)
Data_all <- merge(Data_all, Habitat[,c(1,3:6)], by=c("Site","Location","Year"))


#############################################################################################
##                                                                                         ##  
##                          INTIAL ANALYSIS USING ALL LOCATIONS                            ##
##                                                                                         ##  
#############################################################################################

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all

## Transformation
Data_final[,c(5:7,9)] <- log1p(Data_final[,c(5:7,9)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,8,10:18,21)] <- log(Data_final[,c(4,8,10:18,21)])# log transform Water chem + Kdd
Data_final[,c(19,22)] <- logit(Data_final[,c(19,22)]/100, adjust=0.025)# logit transform Cropping, %WW
Data_final[,c(20)] <- logit(Data_final[,c(20)], adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:22)] <- decostand(Data_final[,c(4:22)],"standardize")

#############################################################################################
##                                                                                         ##  
##                             MODEL PARAMETER SELECTION                                   ##
##                                                                                         ##  
#############################################################################################

#calculate the pearson correlation coefficient matrix
## Drop metals due to NAs
myMat.cor <- cor(Data_final[,-c(1:3)], method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 1 CWM FOOD                              ##
##                                                                                         ##  
#############################################################################################

## Assess hypothetical predictors to help guide model selection
## Fit global model for estimation of adjusted R2
## Note selection of TU data: insecticides, fungicides, and non-Insecticides (all hypothesis driver)

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Food_3 + Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result) ##0.1960471
vif.cca(rda_result)

global <- Data_final[,c(6,12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.1960471, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Food_3"     

## SHREDDERS

Shredders <- Data_final$Food_3

Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## 0.1343934
vif.cca(rda_result)

global <- Data_final[,c(12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1343934, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 2 CWM FEEDING                           ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Drop  + NonInsecticide_TUs
Gmm <- model.matrix( ~ Feeding_3 + Fungicide_TUs + Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result) ## 0.2204156
vif.cca(rda_result)

global <- Data_final[,c(7,12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2204156, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Feeding_3"     

## SHREDDERS

Shredders <- Data_final$Feeding_3

Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## 0.2155021
vif.cca(rda_result)

global <- Data_final[,c(12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.2155021, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##  
##                    FORWARD SELECTION: MODEL 3 GAMMARID ABUNDANCES                       ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Gammaridae + Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result) ## 0.2054795
vif.cca(rda_result)

global <- Data_final[,c(5,12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2054795, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Gammaridae"     

## SHREDDERS

Shredders <- Data_final$Gammaridae

Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN  + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## 0.3003787
vif.cca(rda_result)

global <- Data_final[,c(12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.3003787, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"
# "DIN"  

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all

## Transformation
Data_final[,c(5:7,9)] <- log1p(Data_final[,c(5:7,9)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,8,10:18,21)] <- log(Data_final[,c(4,8,10:18,21)])# log transform Water chem + Kdd
Data_final[,c(19,22)] <- logit(Data_final[,c(19,22)]/100, adjust=0.025)# logit transform Cropping, %WW
Data_final[,c(20)] <- logit(Data_final[,c(20)], adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
#Data_final[,c(4:22)] <- decostand(Data_final[,c(4:22)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1 CWM FOOD                                 ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Food_3 + SRP + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final ),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1 
                      ,NonInsecticide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% NonInsecticide_TUs
                      ,DIN %~~% NonInsecticide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M1)

# AICc
AIC(LPA_pSEM_M1,aicc = T)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_All_CEall_path_estimates.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M1)
#GOF_estimates <- rsquared(LPA_pSEM_M1)
#write.csv(GOF_estimates , "SEM_Food3_M1_nonInsect_All_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 162.253
## Without correlated error terms: BIC = 208.921

## Drop DIN ~ Cropping: BIC = 196.477
## Drop nonInsecticide TUs ~ Cropping: BIC = 159.74
## Drop Insecticide TUs ~ Cropping: BIC = 163.88
## Drop SRP ~ Cropping: BIC = 156.223

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Food3_M1_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1b CWM FOOD                                ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Try Fungicides instead of all non Insecticide TUs
# "Mode of action" hypothesis - fungicides expected to target microbes
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response   
  Leaf_response = lme(Kdd_alder ~ Food_3 + SRP + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final ),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Fungicide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final ),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1 
                      ,Fungicide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% Fungicide_TUs
                      ,DIN %~~% Fungicide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M1)

# AICc
AIC(LPA_pSEM_M1,aicc = T)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1b_Fung_All_CEall_path_estimates.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M1)
#GOF_estimates <- rsquared(LPA_pSEM_M1)
#write.csv(GOF_estimates , "SEM_Food3_M1b_Fung_All_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 164.315
## Without correlated error terms: BIC = 186.734

## Drop DIN ~ Cropping: BIC = 198.539
## Drop Fungicide TUs ~ Cropping: BIC = 163.57
## Drop Insecticide TUs ~ Cropping: BIC = 167.71
## Drop SRP ~ Cropping: BIC = 160.053

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Food3_M1b_Fung_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 2 CWM FEEDING                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response + SRP 
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + SRP + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Feeding_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                              random=~1|Year/Site, 
                              data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2 
                      ,NonInsecticide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% NonInsecticide_TUs
                      ,DIN %~~% NonInsecticide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M2)

# AICc
AIC(LPA_pSEM_M2,aicc = T)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Feeding3_M2_nonInsect_All_CEall_path_estimates.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M2)
#GOF_estimates <- rsquared(LPA_pSEM_M2)
#write.csv(GOF_estimates , "SEM_Feeding3_M2_nonInsect_All_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 161.536
## Without correlated error terms: BIC = 208.204

## Drop DIN ~ Cropping: BIC = 195.805
## Drop Fungicide TUs ~ Cropping: BIC = 159.025
## Drop Insecticide TUs ~ Cropping: BIC = 163.197
## Drop SRP ~ Cropping: BIC = 155.519

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Feeding3_M2_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                SEM: MODEL 2b CWM FEEDING                                ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Try Fungicide TUs instead of all non Insecticide TUs
# "Mode of action" hypothesis - fungicides expected to target microbes
# Note: model below is best-fitting identified by BIC scores

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + SRP + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Feeding_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Fungicide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2 
                      ,Fungicide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% Fungicide_TUs
                      ,DIN %~~% Fungicide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M2)

# AICc
AIC(LPA_pSEM_M2,aicc = T)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Feeding3_M2b_Fung_All_CEal_estimates.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M2)
#GOF_estimates <- rsquared(LPA_pSEM_M2)
#write.csv(GOF_estimates , "SEM_Feeding3_M2b_Fung_All_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 163.478
## Without correlated error terms: BIC = 185.897

## Drop DIN ~ Cropping: BIC = 197.747
## Drop Fungicide TUs ~ Cropping: BIC = 162.808
## Drop Insecticide TUs ~ Cropping: BIC = 166.98
## Drop SRP ~ Cropping: BIC = 159.302

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Feeding3_M2b_Fung_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                            SEM: MODEL 3 GAMMARID ABUNDANCES                             ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

head(Data_final)

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Gammaridae + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting Gammarid response
  Shredder_response = lme(Gammaridae ~ Insecticide_TUs + DIN + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticides_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                              random=~1|Year/Site, 
                              data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3 
                      ,NonInsecticide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% NonInsecticide_TUs
                      ,DIN %~~% NonInsecticide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M3)

# AICc
AIC(LPA_pSEM_M3,aicc = T)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Gammarid_M3_nonInsect_All_CEall_path_est.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M3)
#GOF_estimates <- rsquared(LPA_pSEM_M3)
#write.csv(GOF_estimates , "SEM_Gammarid_M3_nonInsect_All_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 166.482
## Without correlated error terms: BIC = 212.767

## Drop DIN ~ Cropping: BIC = 199.949
## Drop Fungicide TUs ~ Cropping: BIC = 163.647
## Drop Insecticide TUs ~ Cropping: BIC = 167.603
## Drop SRP ~ Cropping: BIC = 160.165

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Gammarid_M3_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                          SEM: MODEL 3b GAMMARID ABUNDANCES                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Try Fungicides instead of all non Insecticide TUs
# "Mode of action" hypothesis - fungicides expected to target microbes
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Gammaridae + SRP + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting Gammarid response
  Shredder_response = lme(Gammaridae ~ Insecticide_TUs + DIN + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticides_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                              random=~1|Year/Site, 
                              data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Fungicide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3 
                      ,Fungicide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% Fungicide_TUs
                      ,DIN %~~% Fungicide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M3)

# AICc
AIC(LPA_pSEM_M3,aicc = T)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Gammarid_M3b_Fung_All_CEall_path_est.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M3)
#GOF_estimates <- rsquared(LPA_pSEM_M3)
#write.csv(GOF_estimates , "SEM_Gammarid_M3b_Fung_All_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 165.008
## Without correlated error terms: BIC = 187.803

## Drop DIN ~ Cropping: BIC = 198.474
## Drop Fungicide TUs ~ Cropping: BIC = 164.943
## Drop Insecticide TUs ~ Cropping: BIC = 168.899
## Drop SRP ~ Cropping: BIC = 161.462

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Gammarid_M3b_Fung_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                  RANK INVERTEBRATE RESPONSES FOR VARIATION EXPLAINED                    ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Feeding_3 explains the most, but variation similar across invert predicto

## Feeding_3
Gmm <- model.matrix( ~ Feeding_3, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## $adj.r.squared
## [1] 0.2257114

## Food_3
Gmm <- model.matrix( ~ Food_3, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## $adj.r.squared
## [1] 0.2050562

## Gammaridae
Gmm <- model.matrix( ~ Gammaridae, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## $adj.r.squared
## [1] 0.2034623

#############################################################################################
##                                                                                         ##  
##                            FINAL ANALYSIS USED IN MANUSCRIPT                            ##
##                                                                                         ##  
#############################################################################################

#############################################################################################
##                                                                                         ##  
##                                 EXCLUDE LOCATION U2                                     ##
##                                                                                         ##  
#############################################################################################

## Here I will drop the Location U2 from the analysis to see how robust the findings are
## using the model structure selected above. This approach is more desirable because we only
## have MP concentrations (ergo MP TUs) for sampling locations U1 and D (U2 above uses dummy
## variables - i.e., the same as U1). Whilst the differences in water quality between U1 and 
## U2 are negligible, there is still variation in responses at sampling location U2 that could 
## be explained by actual data as opposed the dummy variables. For this reason the models below
## will be presented in the result section.

#############################################################################################
##                                                                                         ##  
##                                    TRANSFORMATION                                       ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all[-which(Data_all$Location=="U2"),]

## Transformation
Data_final[,c(5:7,9)] <- log1p(Data_final[,c(5:7,9)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,8,10:18,21)] <- log(Data_final[,c(4,8,10:18,21)])# log transform Water chem + Kdd
Data_final[,c(19,22)] <- logit(Data_final[,c(19,22)]/100, adjust=0.025)# logit transform Cropping, %WW
Data_final[,c(20)] <- logit(Data_final[,c(20)], adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
#Data_final[,c(4:22)] <- decostand(Data_final[,c(4:22)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                    CORRELATION PLOT                                     ##
##                                                                                         ##  
#############################################################################################

## Create plot for interpretation
png(filename="Figure_mU2_correlation_201221.png", 
    type="cairo",
    units="in", 
    width=8.75, 
    height=7.75, 
    pointsize=16, 
    res=600)

#calculate the pearson correlation coefficient matrix
## Drop metals due to NAs
myMat.cor <- cor(Data_final[,-c(1:3,9)], method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
#pheatmap(myMat.cor)

library(magrittr)

myMat.cor %>%
  set_rownames(c("Decomposition","Gammaridae","CWM Food","CWM Feeding","Total TUs","Fungicide TUs","Herbicide TUs","Insecticide TUs","Other TUs",
                 "Pharmaceutical TUs","Non-Insecticide TUs","SRP","DIN","NH4","% WW", "% Cropping","TSS","% Organic")) %>%
  set_colnames(c("Decomposition","Gammaridae","CWM Food","CWM Feeding","Total TUs","Fungicide TUs","Herbicide TUs","Insecticide TUs","Other TUs",
                 "Pharmaceutical TUs","Non-Insecticide TUs","SRP","DIN","NH4","% WW", "% Cropping","TSS","% Organic")) %>%
  pheatmap()

dev.off()

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 1 CWM FOOD                              ##
##                                                                                         ##  
#############################################################################################

## Assess hypothetical predictors to help guide model selection
## Fit global model for estimation of adjusted R2
## Note selection of TU data: insecticides, fungicides, and non-Insecticides (all hypothesis driver)

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Food_3 + Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result) ## 0.2376711
vif.cca(rda_result)

global <- Data_final[,c(6,12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2376711, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Food_3"     

## SHREDDERS

Shredders <- Data_final$Food_3

Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## 0.07906894
vif.cca(rda_result)

global <- Data_final[,c(12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.2679779, alpha = 0.3)
as.matrix(M1$variables) ## Selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "Insecticide_TUs"
## "SRP"            
## "DIN" 

## Use forward selection with appropriate criteria to assess which variables are essential
#M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.07906894, alpha = 0.05)
#as.matrix(M1$variables)
#rm(M1,global)
## No variables selected using adjusted R-sqr and alpha = 0.05


#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 2 CWM FEEDING                           ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Drop  + NonInsecticide_TUs
Gmm <- model.matrix( ~ Feeding_3 + Fungicide_TUs + Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result) ## 0.2291435
vif.cca(rda_result)

global <- Data_final[,c(7,12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2291435, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Feeding_3"     

## SHREDDERS

Shredders <- Data_final$Feeding_3

Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## 0.1208073
vif.cca(rda_result)

global <- Data_final[,c(12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.2786112, alpha = 0.3)
as.matrix(M1$variables) ## Selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "Insecticide_TUs"
## "SRP"            
## "DIN" 

## Use forward selection with appropriate criteria to assess which variables are essential
#M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1208073, alpha = 0.05)
#as.matrix(M1$variables)
#rm(M1,global)
## No variables selected using adjusted R-sqr and alpha = 0.05

#############################################################################################
##                                                                                         ##  
##                    FORWARD SELECTION: MODEL 3 GAMMARID ABUNDANCES                       ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Gammaridae + Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result) ## 0.2640745
vif.cca(rda_result)

global <- Data_final[,c(5,12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2640745, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Gammaridae"     

## SHREDDERS

Shredders <- Data_final$Gammaridae

Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + NonInsecticide_TUs + SRP + DIN  + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## 0.3003787
vif.cca(rda_result)

global <- Data_final[,c(12,10,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.4036372, alpha = 0.3)
as.matrix(M1$variables) ## Selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "Insecticide_TUs"
## "SRP"            
## "DIN" 
## "pORG"  

## Use forward selection with appropriate criteria to assess which variables are essential
#M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.2786112, alpha = 0.05)
#as.matrix(M1$variables)
#rm(M1,global)
## No variables selected using adjusted R-sqr and alpha = 0.05

#############################################################################################
##                                                                                         ##  
##                                    TRANSFORMATION                                       ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all[-which(Data_all$Location=="U2"),]

## Transformation
Data_final[,c(5:7,9)] <- log1p(Data_final[,c(5:7,9)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,8,10:18,21)] <- log(Data_final[,c(4,8,10:18,21)])# log transform Water chem + Kdd
Data_final[,c(19,22)] <- logit(Data_final[,c(19,22)]/100, adjust=0.025)# logit transform Cropping, %WW
Data_final[,c(20)] <- logit(Data_final[,c(20)], adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:22)] <- decostand(Data_final[,c(4:22)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1 CWM FOOD                                 ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Food_3 + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final ),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1 
                      ,NonInsecticide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% NonInsecticide_TUs
                      ,DIN %~~% NonInsecticide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M1)

# AICc
AIC(LPA_pSEM_M1,aicc = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
## Evaluate path significance using unstandardized coefficients (but standardized data)
Path_estimates <- coefs(LPA_pSEM_M1, standardize = "none", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_mU2_path_est_std.csv",row.names = T)

# Evaluate path significance using standardized coefficients
coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_mU2_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M1)
GOF_estimates <- rsquared(LPA_pSEM_M1)
write.csv(GOF_estimates , "SEM_Food3_M1_nonInsect_mU2_rsquared.csv",row.names = T)

# INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 145.719
## Without correlated error terms: BIC = 185.999

## Drop DIN ~ Cropping: BIC = 174.668
## Drop Fungicide TUs ~ Cropping: BIC = 143.816
## Drop Insecticide TUs ~ Cropping: BIC = 147.649
## Drop SRP ~ Cropping: BIC = 140.916

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
out <- akaike.weights(c(
135.092,
136.086,
138.839,
139.833,
140.916
))
write.csv(out$weights, "SEM_Food3_M1_nonInsect_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1b CWM FOOD                                ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Try Fungicides instead of all non Insecticide TUs
# "Mode of action" hypothesis - fungicides expected to target microbes
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M1b = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Food_3 + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final ),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Fungicide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final ),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M1b <- update(LPA_pSEM_M1b 
                       ,Fungicide_TUs %~~% Insecticide_TUs
                       ,SRP %~~% Fungicide_TUs
                       ,DIN %~~% Fungicide_TUs
                       ,DIN %~~% Insecticide_TUs
                       ,SRP %~~% Insecticide_TUs
                       ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1b, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M1b)

# AICc
AIC(LPA_pSEM_M1b,aicc = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
## Evaluate path significance using unstandardized coefficients (but standardized data)
#Path_estimates <- coefs(LPA_pSEM_M1b, standardize = "none", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1b_Fung_mU2_path_est_std.csv",row.names = T)

# Evaluate path significance using standardized coefficients
coefs(LPA_pSEM_M1b, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M1b, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Food3_M1b_Fung_mU2_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M1b)
GOF_estimates <- rsquared(LPA_pSEM_M1b)
write.csv(GOF_estimates , "SEM_Food3_M1b_Fung_mU2_rsquared.csv",row.names = T)

# INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 147.663
## Without correlated error terms: BIC = 168.23

## Drop DIN ~ Cropping: BIC = 176.612
## Drop Fungicide TUs ~ Cropping: BIC = 147.574
## Drop Insecticide TUs ~ Cropping: BIC = 151.408
## Drop SRP ~ Cropping: BIC = 144.675

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
out <- akaike.weights(c(
  138.483,
  139.845,
  140.808,
  142.102,
  143.464,
  144.675,
  145.364,
  145.426
))
write.csv(out$weights, "SEM_Food3_M1b_Fung_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 2 CWM FEEDING                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Feeding_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2 
                      ,NonInsecticide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% NonInsecticide_TUs
                      ,DIN %~~% NonInsecticide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M2)

# AICc
AIC(LPA_pSEM_M2,aicc = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
## Evaluate path significance using unstandardized coefficients (but standardized data)
Path_estimates <- coefs(LPA_pSEM_M2, standardize = "none", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Feed3_M2_nonInsect_mU2_path_est_std.csv",row.names = T)

# Evaluate path significance using standardized coefficients
coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Feed3_M2_nonInsect_mU2_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M2)
GOF_estimates <- rsquared(LPA_pSEM_M2)
write.csv(GOF_estimates , "SEM_Feed3_M2_nonInsect_mU2_rsquared.csv",row.names = T)

# INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 145.72
## Without correlated error terms: BIC = 185.999

## Drop DIN ~ Cropping: BIC = 174.766
## Drop Fungicide TUs ~ Cropping: BIC = 143.757
## Drop Insecticide TUs ~ Cropping: BIC = 147.675
## Drop SRP ~ Cropping: BIC = 140.858

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
out <- akaike.weights(c(
  134.447,
  136.337,
  138.141,
  140.032,
  140.858
))
write.csv(out$weights, "SEM_Feed3_M2_nonInsect_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                SEM: MODEL 2b CWM FEEDING                                ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Try Fungicide TUs instead of all non Insecticide TUs
# "Mode of action" hypothesis - fungicides expected to target microbes
# Note: model below is best-fitting identified by BIC scores

# Create component models and store in list
LPA_pSEM_M2b = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Feeding_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Fungicide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M2b <- update(LPA_pSEM_M2b 
                       ,Fungicide_TUs %~~% Insecticide_TUs
                       ,SRP %~~% Fungicide_TUs
                       ,DIN %~~% Fungicide_TUs
                       ,DIN %~~% Insecticide_TUs
                       ,SRP %~~% Insecticide_TUs
                       ,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M2b, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M2b)

# AICc
AIC(LPA_pSEM_M2b,aicc = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
## Evaluate path significance using unstandardized coefficients (but standardized data)
Path_estimates <- coefs(LPA_pSEM_M2b, standardize = "none", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Feed3_M2b_Fung_mU2_path_est_std.csv",row.names = T)

# Evaluate path significance using standardized coefficients
#coefs(LPA_pSEM_M2b, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M2b, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Feed3_M2b_Fung_mU2_path_est.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M2b)
GOF_estimates <- rsquared(LPA_pSEM_M2b)
write.csv(GOF_estimates , "SEM_Feed3_M2b_Fung_mU2_rsquared.csv",row.names = T)

# INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 147.743
## Without correlated error terms: BIC = 168.309

## Drop DIN ~ Cropping: BIC = 176.789
## Drop Fungicide TUs ~ Cropping: BIC = 147.666
## Drop Insecticide TUs ~ Cropping: BIC = 151.583
## Drop SRP ~ Cropping: BIC =

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
out <- akaike.weights(c(
  138.104,
  140.246,
  139.867,
  141.612,
  143.753,
  144.766,
  144.552,
  144.717
))
write.csv(out$weights, "SEM_Feed3_M2b_Fung_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                            SEM: MODEL 3 GAMMARID ABUNDANCES                             ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

head(Data_final)

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response  + WW_fraction
  Leaf_response = lme(Kdd_alder ~  Gammaridae + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting Gammarid response
  Shredder_response = lme(Gammaridae ~ Insecticide_TUs + SRP + DIN,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticides_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                              random=~1|Year/Site, 
                              data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Organic sediment response
  #Organics_response = lme(pORG ~ WW_fraction,
  #                              random=~1|Year/Site, 
  #                              data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3 
                      ,NonInsecticide_TUs %~~% Insecticide_TUs
                      ,SRP %~~% NonInsecticide_TUs
                      ,DIN %~~% NonInsecticide_TUs
                      ,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      ,DIN %~~% SRP
                      #,pORG %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M3)

# AICc
AIC(LPA_pSEM_M3,aicc = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
## Evaluate path significance using unstandardized coefficients (but standardized data)
Path_estimates <- coefs(LPA_pSEM_M3, standardize = "none", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Gammarid_M3_nonInsect_mU2_path_est_std.csv",row.names = T)

# Evaluate path significance using standardized coefficients
coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Gammarid_M3_nonInsect_mU2_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M3)
GOF_estimates <- rsquared(LPA_pSEM_M3)
write.csv(GOF_estimates , "SEM_Gammarid_M3_nonInsect_mU2_rsquared.csv",row.names = T)

# INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 176.132
## Without correlated error terms: BIC = 218.477

## Drop DIN ~ Cropping: BIC = 205.256
## Drop pORG ~ Cropping: BIC = 173.031
## Drop Fungicide TUs ~ Cropping: BIC = 171.159
## Drop Insecticide TUs ~ Cropping: BIC = 174.634
## Drop SRP ~ Cropping: BIC = 168.26

## Drop pORG ~ Gammaridae: BIC = 166.083

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
out <- akaike.weights(c(
  138.005,
  141.516,
  141.682,
  145.964,
  147.986,
  149.475
))
write.csv(out$weights, "SEM_Gammarid_M3_nonInsect_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                          SEM: MODEL 3b GAMMARID ABUNDANCES                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Try Fungicides instead of all non Insecticide TUs
# "Mode of action" hypothesis - fungicides expected to target microbes
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M3b = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Gammaridae + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting Gammarid response
  Shredder_response = lme(Gammaridae ~ Insecticide_TUs + SRP + DIN,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction ,
                     random=~1|Year/Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticides_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                              random=~1|Year/Site, 
                              data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Fungicide_TUs ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Organic sediment response
  #Organics_response = lme(pORG ~ WW_fraction,
  #                              random=~1|Year/Site, 
  #                              data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M3b <- update(LPA_pSEM_M3b 
                       ,Fungicide_TUs %~~% Insecticide_TUs
                       ,SRP %~~% Fungicide_TUs
                       ,DIN %~~% Fungicide_TUs
                       ,DIN %~~% Insecticide_TUs
                       ,SRP %~~% Insecticide_TUs
                       ,DIN %~~% SRP
                       #,pORG %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M3b, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M3b)

# AICc
AIC(LPA_pSEM_M3b,aicc = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
## Evaluate path significance using unstandardized coefficients (but standardized data)
Path_estimates <- coefs(LPA_pSEM_M3b, standardize = "none", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Gammarid_M3b_Fung_mU2_path_est_std.csv",row.names = T)

# Evaluate path significance using standardized coefficients
#coefs(LPA_pSEM_M3b, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M3b, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Gammarid_M3b_Fung_mU2_path_est.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M3b)
GOF_estimates <- rsquared(LPA_pSEM_M3b)
write.csv(GOF_estimates , "SEM_Gammarid_M3b_Fung_mU2_rsquared.csv",row.names = T)

# INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 177.955
## Without correlated error terms: BIC = 200.587

## Drop DIN ~ Cropping: BIC = 207.079
## Drop pORG ~ Cropping: BIC = 174.853
## Drop Fungicide TUs ~ Cropping: BIC = 174.726
## Drop Insecticide TUs ~ Cropping: BIC = 178.201
## Drop SRP ~ Cropping: BIC = 171.828

## Drop pORG ~ Gammaridae: BIC = 170.597

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
out <- akaike.weights(c(
  141.535,
  144.214,
  145.212,
  147.878,
  150.557,
  151.197,
  153.432,
  155.300
))
write.csv(out$weights, "SEM_Gammarid_M3b_Fung_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                  RANK INVERTEBRATE RESPONSES FOR VARIATION EXPLAINED                    ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Food_3 explains the most, but variation similar across invert predictors

## Food_3
Gmm <- model.matrix( ~ Food_3, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## $adj.r.squared
## [1] 0.3025651

## Gammaridae
Gmm <- model.matrix( ~ Gammaridae, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## $adj.r.squared
## [1] 0.297116

## Feeding_3
Gmm <- model.matrix( ~ Feeding_3, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## $adj.r.squared
## [1] 0.2903887

#############################################################################################
##                                                                                         ##  
##                           2014: EXPLICITLY CONSIDER METALS                              ##
##                                                                                         ##  
#############################################################################################

#############################################################################################
##                                                                                         ##  
##                              INCLUDE METALS EXPLICITLY                                  ##
##                              MODEL PARAMETER SELECTION                                  ##
##                                                                                         ##  
#############################################################################################

## Here I will take the 2014 data and explicitly test what influence Metals TU data have by: 
## 1) Using Metal TUs with other predictors to forward-select parameters using packfor
## 2) Using Metal TUs explicitly in the SEMs to see if they have any effects

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

## Select only 2014 sites
Data_final <- Data_all[-which(Data_all$Location=="U2"),]
Data_final <- Data_final[which(Data_final$Year=="2014"),]

## Transformation
Data_final[,c(5:7,9)] <- log1p(Data_final[,c(5:7,9)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,8,10:18,21)] <- log(Data_final[,c(4,8,10:18,21)])# log transform Water chem + Kdd
Data_final[,c(19,22)] <- logit(Data_final[,c(19,22)]/100, adjust=0.025)# logit transform Cropping, %WW
Data_final[,c(20)] <- logit(Data_final[,c(20)], adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:22)] <- decostand(Data_final[,c(4:22)],"standardize")

#############################################################################################
##                                                                                         ##  
##                             MODEL PARAMETER SELECTION                                   ##
##                                                                                         ##  
#############################################################################################

#calculate the pearson correlation coefficient matrix
## Drop metals due to NAs
myMat.cor <- cor(Data_final[,-c(1:3)], method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)

## Assess hypothetical predictors to help guide model selection
## Fit global model for estimation of adjusted R2

#############################################################################################
##                                                                                         ##
##                               CWM FOOD 3 MEAN ABUNDANCES                                ##
##                        INSECTICIDES + NON-INSECTICIDES + METALS                         ##  
##                                                                                         ##  
#############################################################################################

## See correlation with SRP among others
plot(SRP ~ Metal_TUs, Data_final)
abline(lm(SRP ~ Metal_TUs, Data_final),col="red")
summary(lm(SRP ~ Metal_TUs, Data_final)) ## highly correlated

##**********************
## Decompostion
Decomposition <- Data_final$Kdd_alder
## NonInsecticides (i.e., excluding Insecticide and Metals)
## Not include Metal_TUs due to VIF > 10 (i.e., 21.403029)
Gmm <- model.matrix( ~ Food_3 + Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.09224922

## Compared with Metal_TUs alone
Gmm <- model.matrix( ~ Metal_TUs, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.01139914

## Compared with Shredders
Gmm <- model.matrix( ~ Food_3, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.2730072

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(6,12,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.09224922, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"

##**********************************
Shredders <- Data_final$Food_3
## Shredders
## Not include Metal_TUs due to VIF > 10 (i.e., 20.305173)
Gmm <- model.matrix( ~ Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.3475552

## Metal_TUs
Gmm <- model.matrix( ~ Metal_TUs, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr -0.07141368

## Compared with Insecticide_TUs
Gmm <- model.matrix( ~ Insecticide_TUs, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.1808922

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(12,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.3475552, alpha = 0.05)
as.matrix(M1$variables)
#rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##
##                           INSECTICIDES + FUNGICIDES + METALS                            ##  
##                                                                                         ##  
#############################################################################################

## Decomposition
## Not include Metal_TUs due to VIF > 10 (i.e., 13.727341)
Gmm <- model.matrix( ~ Food_3 + Insecticide_TUs + Fungicide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.1298087

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(6,12,10,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.1298087, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"

##*****************
## Shredders
Shredders <- Data_final$Food_3
## Not include Metal_TUs due to VIF > 10 (i.e., 12.038968)
Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.3094362

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(12,10,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.3094362, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##
##                             CWM FEEDING 3 MEAN ABUNDANCES                               ##
##                        INSECTICIDES + NON-INSECTICIDES + METALS                         ##  
##                                                                                         ##  
#############################################################################################

##**********************
## Decompostion
Decomposition <- Data_final$Kdd_alder
## NonInsecticides (i.e., excluding Insecticide and Metals)
## Not include Metal_TUs due to VIF > 10 (i.e., 21.619592)
Gmm <- model.matrix( ~ Feeding_3 + Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.0853319

## Compared with Metal_TUs alone
Gmm <- model.matrix( ~ Metal_TUs, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.01139914

## Compared with Shredders
Gmm <- model.matrix( ~ Feeding_3, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.2664518

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(7,12,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.0853319, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"

##**********************************
Shredders <- Data_final$Feeding_3
## Shredders
## Not include Metal_TUs due to VIF > 10 (i.e., 20.305173)
Gmm <- model.matrix( ~ Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.3518244

## Metal_TUs
Gmm <- model.matrix( ~ Metal_TUs, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr -0.07140588

## Compared with Insecticide_TUs
Gmm <- model.matrix( ~ Insecticide_TUs, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.1794218

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(12,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.3518244, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##
##                           INSECTICIDES + FUNGICIDES + METALS                            ##  
##                                                                                         ##  
#############################################################################################

## Decomposition
Decomposition <- Data_final$Kdd_alder
## Not include Metal_TUs due to VIF > 10 (i.e., 13.671125)
Gmm <- model.matrix( ~ Feeding_3 + Insecticide_TUs + Fungicide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.1228746

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(6,12,10,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.1228746, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"

##*****************
## Shredders
Shredders <- Data_final$Feeding_3
## Not include Metal_TUs due to VIF > 10 (i.e., 12.038968)
Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.3149984

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(12,10,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.3149984, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"  

#############################################################################################
##                                                                                         ##
##                             GAMMARIDAE MEAN ABUNDANCES                                  ##
##                        INSECTICIDES + NON-INSECTICIDES + METALS                         ##  
##                                                                                         ##  
#############################################################################################

##**********************
## Decompostion
Decomposition <- Data_final$Kdd_alder
## NonInsecticides (i.e., excluding Insecticide and Metals)
## Not include Metal_TUs due to VIF > 10 (i.e., 21.758938)
Gmm <- model.matrix( ~ Gammaridae + Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.6360172

## Compared with Metal_TUs alone
Gmm <- model.matrix( ~ Metal_TUs, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.01139914

## Compared with Shredders
Gmm <- model.matrix( ~ Gammaridae, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.5573293

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(5,12,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.6360172, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Gammaridae"
# "DIN" 

##**********************************
Shredders <- Data_final$Gammaridae
## Shredders
## Not include Metal_TUs due to VIF > 10 (i.e., 20.305173)
Gmm <- model.matrix( ~ Insecticide_TUs + NonInsecticide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.4692595

## Metal_TUs
Gmm <- model.matrix( ~ Metal_TUs, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr -0.0687695

## Compared with Insecticide_TUs
Gmm <- model.matrix( ~ Insecticide_TUs, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
## adj.r.sqr 0.2788009

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(12,15,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh =  0.4692595, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##
##                           INSECTICIDES + FUNGICIDES + METALS                            ##  
##                                                                                         ##  
#############################################################################################

## Decomposition
Decomposition <- Data_final$Kdd_alder
## Not include Metal_TUs due to VIF > 10 (i.e., 13.777949 )
Gmm <- model.matrix( ~ Gammaridae + Insecticide_TUs + Fungicide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.6260059

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(5,12,10,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.6260059, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Gammaridae"

##*****************
## Shredders
Shredders <- Data_final$Gammaridae
## Not include Metal_TUs due to VIF > 10 (i.e., 12.038968)
Gmm <- model.matrix( ~ Insecticide_TUs + Fungicide_TUs + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)
## adj.r.sqr 0.3848767

## Metals strongly correlated with nutrients among others - so drop for forward selection
global <- Data_final[,c(12,10,16,17,18,21,22)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.3848767, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Insecticide_TUs"
# "SRP"   

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

## Select only 2014 sites
Data_final <- Data_all[-which(Data_all$Location=="U2"),]
Data_final <- Data_final[which(Data_final$Year=="2014"),]

## Transformation
Data_final[,c(5:7,9)] <- log1p(Data_final[,c(5:7,9)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,8,10:18,21)] <- log(Data_final[,c(4,8,10:18,21)])# log transform Water chem + Kdd
Data_final[,c(19,22)] <- logit(Data_final[,c(19,22)]/100, adjust=0.025)# logit transform Cropping, %WW
Data_final[,c(20)] <- logit(Data_final[,c(20)], adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
#Data_final[,c(4:22)] <- decostand(Data_final[,c(4:22)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 4 CWM FOOD                                 ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site as random effect)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M4 = psem(
  
  # Predicting LPA response + WW_fraction
  Leaf_response = lme(Kdd_alder ~ Food_3 + DIN,
                      random=~1|Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ Insecticide_TUs + SRP,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Site, 
                             data = Data_final),
  
  # Predicting non-Insecticide response
  Metals_response = lme(Metal_TUs ~ WW_fraction,
                   random=~1|Site, 
                  data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M4 <- update(LPA_pSEM_M4, 
                      Metal_TUs %~~% Insecticide_TUs,
                      SRP %~~% Metal_TUs,
                      DIN %~~% Metal_TUs,
                      DIN %~~% Insecticide_TUs,
                      SRP %~~% Insecticide_TUs,
                      DIN %~~% SRP
                      )

# Run goodness-of-fit tests
summary(LPA_pSEM_M4, standardize = "none", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M4)

# AICc
AIC(LPA_pSEM_M4,aicc = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 99.26
## Without correlated error terms: BIC = 134.832

## Drop DIN ~ Cropping: BIC = 104.407
## Drop Fungicide TUs ~ Cropping: BIC = 96.782
## Drop Insecticide TUs ~ Cropping: BIC = 97.67
## Drop SRP ~ Cropping: BIC = 98.534

out <- akaike.weights(c(
  93.260,
  95.370,
  95.676,
  95.963,
  96.000,
  96.782,
  97.587,
  97.784,
  98.011,
  98.299,
  100.221,
  100.628,
  104.512,
  105.173
))
write.csv(out$weights, "SEM_Food3_M4_Metals_2014_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                SEM: MODEL 5 CWM FEEDING 3                               ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M5 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + DIN,
                      random=~1|Site,
                      data = Data_final),
  
  # Predicting shredder response 
  Shredder_response = lme(Feeding_3 ~ Insecticide_TUs + SRP,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Site, 
                             data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Metal_TUs ~ WW_fraction,
                                random=~1|Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M5 <- update(LPA_pSEM_M5, 
                      Metal_TUs %~~% Insecticide_TUs,
                      SRP %~~% Metal_TUs,
                      DIN %~~% Metal_TUs,
                      DIN %~~% Insecticide_TUs,
                      SRP %~~% Insecticide_TUs,
                      DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M5, standardize = "none", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M5)

# AICc
AIC(LPA_pSEM_M5,aicc = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 99.264
## Without correlated error terms: BIC = 134.836

## Drop DIN ~ Cropping: BIC = 104.421
## Drop Fungicide TUs ~ Cropping: BIC = 96.764
## Drop Insecticide TUs ~ Cropping: BIC = 97.652
## Drop SRP ~ Cropping: BIC = 98.516

## Gammarids show evidence of negative indirect insecticide impact on decomposition
out <- akaike.weights(c(
  93.338,
  95.757,
  96.010,
  96.034,
  96.034,
  96.764,
  97.582,
  97.926,
  98.117,
  98.371,
  100.355,
  100.482,
  104.406,
  105.081
))
write.csv(out$weights, "SEM_Feed3_M5_Metals_2014_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                SEM: MODEL 6 GAMMARIDAE                                  ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M6 = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Gammaridae + DIN,
                      random=~1|Site,
                      data = Data_final),
  
  # Predicting shredder response 
  Shredder_response = lme(Gammaridae ~ Insecticide_TUs + SRP,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Site, 
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Site, 
                             data = Data_final),
  
  # Predicting non-Insecticide response
  nonInsecticide_response = lme(Metal_TUs ~ WW_fraction,
                                random=~1|Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M6 <- update(LPA_pSEM_M6, 
                      Metal_TUs %~~% Insecticide_TUs,
                      SRP %~~% Metal_TUs,
                      DIN %~~% Metal_TUs,
                      DIN %~~% Insecticide_TUs,
                      SRP %~~% Insecticide_TUs,
                      DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M6, standardize = "none", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M6)

# AICc
AIC(LPA_pSEM_M6,aicc = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 95.812
## Without correlated error terms: BIC = 131.384

## Drop DIN ~ Cropping: BIC = 102.596
## Drop Fungicide TUs ~ Cropping: BIC = 93.753
## Drop Insecticide TUs ~ Cropping: BIC = 95.2
## Drop SRP ~ Cropping: BIC = 92.841

out <- akaike.weights(c(
  90.401,
  92.483,
  92.696,
  92.841,
  93.441,
  93.486,
  94.484,
  94.697,
  95.523,
  99.091,
  99.863,
  101.010,
  101.167,
  101.367
))
write.csv(out$weights, "SEM_Gammarid_M6_Metals_2014_mU2_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                         SIMPLIFIED SEM: MODEL 1 CWM FOOD                                ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
## Could present if challenged about number of path lengths in above models
## Clearly shows same patterns with similar std. path estimates

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response  + WW_fraction
  Leaf_response = lme(Kdd_alder ~ Food_3 + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ Insecticide_TUs + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Insecticide response
  Insecticide_response = lme(Insecticide_TUs ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final ))
  
  # Predicting non-Insecticide response
  #nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction ,
  #                              random=~1|Year/Site, 
  #                              data = Data_final),
  
  # Predicting Nitrogen response
  #DIN_response = lme(DIN ~ WW_fraction + Cropping,
  #                   random=~1|Year/Site, 
  #                   data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1 
                      #,NonInsecticide_TUs %~~% Insecticide_TUs
                      #,SRP %~~% NonInsecticide_TUs
                      #,DIN %~~% NonInsecticide_TUs
                      #,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      #,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)

# Model adequacy
BIC(LPA_pSEM_M1)

# AICc
AIC(LPA_pSEM_M1,aicc = T)

# Evaluate path significance using unstandardized coefficients
coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_mU2_CEall_path_est.csv",row.names = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
#Path_estimates <- coefs(LPA_pSEM_M1, standardize = "none", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_mU2_CEall_path_est_stdized_errors.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M1)
#GOF_estimates <- rsquared(LPA_pSEM_M1)
#write.csv(GOF_estimates , "SEM_Food3_M1_nonInsect_mU2_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 134.667
## Without correlated error terms: BIC = 170.576

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path

#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Food3_M1_nonInsect_mU2_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                         SIMPLIFIED SEM: MODEL 1 CWM FOOD                                ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
## Could present if challenged about number of path lengths in above models
## Clearly shows same patterns with similar std. path estimates

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response  + WW_fraction
  Leaf_response = lmer(Kdd_alder ~ Food_3 + WW_fraction + (1|Site),
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lmer(Food_3 ~ Insecticide_TUs + SRP + (1|Site),
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lmer(SRP ~ WW_fraction + (1|Site),
                     data = Data_final),
  
  # Predicting Insecticide response
  Insecticide_response = lmer(Insecticide_TUs ~ WW_fraction + Cropping + (1|Site),
                             data = Data_final))

# Predicting non-Insecticide response
#nonInsecticide_response = lme(NonInsecticide_TUs ~ WW_fraction ,
#                              random=~1|Year/Site, 
#                              data = Data_final),

# Predicting Nitrogen response
#DIN_response = lme(DIN ~ WW_fraction + Cropping,
#                   random=~1|Year/Site, 
#                   data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1 
                      #,NonInsecticide_TUs %~~% Insecticide_TUs
                      #,SRP %~~% NonInsecticide_TUs
                      #,DIN %~~% NonInsecticide_TUs
                      #,DIN %~~% Insecticide_TUs
                      ,SRP %~~% Insecticide_TUs
                      #,DIN %~~% SRP
)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)

out <- suppressWarnings(semEff(LPA_pSEM_M1))

# Model adequacy
BIC(LPA_pSEM_M1)

# AICc
AIC(LPA_pSEM_M1,aicc = T)

# Evaluate path significance using unstandardized coefficients
coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_mU2_CEall_path_est.csv",row.names = T)

## Note: run final model with standardized data to get the standardized errors
## The scaled "Std.estimates" in the output above are the same as the Estimates using standardized data
#Path_estimates <- coefs(LPA_pSEM_M1, standardize = "none", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Food3_M1_nonInsect_mU2_CEall_path_est_stdized_errors.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M1)
#GOF_estimates <- rsquared(LPA_pSEM_M1)
#write.csv(GOF_estimates , "SEM_Food3_M1_nonInsect_mU2_CEall_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 134.667
## Without correlated error terms: BIC = 170.576

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path

#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Food3_M1_nonInsect_mU2_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                         SIMPLIFIED SEM: MODEL 1 CWM FOOD                                ## 
##          ATTEMPT TO ASSESS DIRECT AND INDIRECT EFFECTS WITH BOOTSTRAPPING               ##  
##                                                                                         ##  
#############################################################################################

# install.packages
#library(semEff)

## Investigate key correlations
#plot(Kdd_alder ~ WW_fraction,data = Data_final)
#cor.test(Data_final$Kdd_alder,Data_final$WW_fraction)

#plot(Food_3 ~ WW_fraction,data = Data_final)
#cor.test(Data_final$Food_3,Data_final$WW_fraction)

#plot(Food_3 ~ Cropping,data = Data_final)
#cor.test(Data_final$Food_3,Data_final$Cropping)

#plot(Kdd_alder ~ Cropping,data = Data_final)
#cor.test(Data_final$Kdd_alder, Data_final$Cropping)

## Assess total, direct, and indirect effects
## Need to fit piecewiseSEM with lmer4 for semEff package
## Note model fitted using lme4 package generates error ?isSingular (i.e., random effects are close to zero)
## Model respecified with only Site as the random effect
## However, semEff seems to struggle with assessing Insecticides as a mediator of WW effects, since
## % Cropping also contributes to the variation in Insecticide TUs
## Might be more effective to calculate the indirect effects by hand

# NOT RUN {
## Specification
#head(Data_final)
# }
# NOT RUN {
#Decomposition.SEM <- list(
#  "Decompostion" = lme4::lmer(Kdd_alder ~ Food_3 + WW_fraction + (1 | Site),
#                              data = Data_final),
#  "Shredder" = lme4::lmer(Food_3 ~ Insecticide_TUs + SRP +  (1 | Site),
#                          data = Data_final),
#  "SRP" = lme4::lmer(SRP ~ WW_fraction + (1 | Site),
#                        data = Data_final),
#  "Insecticide" = lme4::lmer(Insecticide_TUs ~ WW_fraction + Cropping + (1 | Site),
#                             data = Data_final)
#)
# }

## View model structure
#lapply(Decomposition.SEM, formula)

## Bootstrap model effects (10000 reps... can take a while)
#system.time(
#  Decomposition.SEM.Boot <- bootEff(Decomposition.SEM, ran.eff = "Site", seed = 53908)
# )

## Calculate SEM effects and CIs (use saved bootstrapped SEM)
#eff <- suppressWarnings(semEff(Decomposition.SEM.Boot))

## Summary of effects for response "Growth"
#eff$Summary$Decompostion
#tot <- totEff(eff)[["Decompostion"]]
#tot.b <- totEff(eff, type = "boot")[["Decompostion"]]

## Summary of effects for response "Growth"
#eff$Summary$Shredder
#tot <- totEff(eff)[["Shredder"]]
#tot.b <- totEff(eff, type = "boot")[["Shredder"]]


#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

## Create data frame for transformed and standardized data
Data_final <- Data_all

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

## Standardisation
Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

## Check data input
head(Data_final)


##***********************************************************************************************************
## Create plot with partial regressions for four main relationships shown in Fig.2 Main text
png(filename="ECOIMPACT_partial_regression_LME_insectidesTUs_dm_mu_211130.png", 
    type="cairo",
    units="in", 
    width=10, 
    height=8, 
    pointsize=20, 
    res=600)

##***********************************************************************************************************
## Detritivores | SRP ~ Insecticide TUs
plot(Food_3~SRP, data=Data_final)
abline(lm(Food_3~SRP, data=Data_final),col="red")
lm_Food3_SRP_bs<-summary(lm(Food_3~SRP, data=Data_final))
lm_Food3_SRP_residual<-lm_Food3_SRP_bs$residuals

## Create model with transformed data for plot
M1 <- lmer(lm_Food3_SRP_residual ~ Insecticide_TUs + (1|Site), 
           data=Data_final,
           lmerControl(optimizer = "Nelder_Mead"),
           REML=T)

## Assess parameter estimates
sjPlot::tab_model(M1)

# Create table with effects (parameters estimates)
effects_Insecticide_TUs <- effects::effect(term= "Insecticide_TUs", mod = M1)
summary(effects_Insecticide_TUs) #output of what the values are

# Save the effects values as a df:
x_Insecticide_TUs <- as.data.frame(effects_Insecticide_TUs)

## Check unique
unique(Data_final$Type_1)

## Head
head(Data_final)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Add column
Data_final$lm_Food3_SRP_residual <- lm_Food3_SRP_residual


## Create plot
Insecticide_TUs_plot <- ggplot() + 
  #1
  geom_line(data=x_Insecticide_TUs, aes(x=Insecticide_TUs, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_Insecticide_TUs, aes(x=Insecticide_TUs, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Data_final, aes(Insecticide_TUs, lm_Food3_SRP_residual, fill=Location), shape=21, size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  #scale_shape_manual(breaks=c("D", "U1", "U2"),
  #                   values=c(15, 17, 16),
  #                   labels = c("Unbuffered", "Buffered", "Forest"),
  #                   name = "Site type"
  #                   #, guide = FALSE
  #) +
  #5
  scale_fill_manual(breaks=c("D", "U1", "U2"),
                    values=c("firebrick","royalblue1", "royalblue4"),
                    labels = c("D", "U1", "U2"),
                    name = "Location"
                    #, guide = FALSE
  ) +
  #6
  labs(x="Insecticides (TUs)", y="Detritivores | Phosphorus") +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 13, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 13, face="bold", color="black"),
        axis.text.y=element_text(size = 13, color="black"),
        axis.text.x=element_text(size = 13, color="black"))

##***********************************************************************************************************
## Detritivores | Insecticide TUs ~ SRP
plot(Food_3~Insecticide_TUs, data=Data_final)
abline(lm(Food_3~Insecticide_TUs, data=Data_final),col="red")
lm_Food3_Insecticide_TUs_bs<-summary(lm(Food_3~Insecticide_TUs, data=Data_final))
lm_Food3_Insecticide_TUs_residual<-lm_Food3_Insecticide_TUs_bs$residuals

## Create model with transformed data for plot
M2 <- lmer(lm_Food3_Insecticide_TUs_residual ~ SRP + (1|Site), 
           data=Data_final,
           lmerControl(optimizer = "Nelder_Mead"),
           REML=T)

## Assess parameter estimates
sjPlot::tab_model(M2)

# Create table with effects (parameters estimates)
effects_SRP <- effects::effect(term= "SRP", mod = M2)
summary(effects_SRP) #output of what the values are

# Save the effects values as a df:
x_SRP <- as.data.frame(effects_SRP)

## Check unique
unique(Data_final$Type_1)

## Head
head(Data_final)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Add column
Data_final$lm_Food3_Insecticide_TUs_residual <- lm_Food3_Insecticide_TUs_residual

## Create plot
SRP_plot <- ggplot() + 
  #1
  geom_line(data=x_SRP, aes(x=SRP, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_SRP, aes(x=SRP, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Data_final, aes(SRP, lm_Food3_Insecticide_TUs_residual, fill=Location), shape=21, size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  #scale_shape_manual(breaks=c("D", "U1", "U2"),
  #                   values=c(15, 17, 16),
  #                   labels = c("Unbuffered", "Buffered", "Forest"),
  #                   name = "Site type"
  #                   #, guide = FALSE
  #) +
  #5
  scale_fill_manual(breaks=c("D", "U1", "U2"),
                    values=c("firebrick","royalblue1", "royalblue4"),
                    labels = c("D", "U1", "U2"),
                    name = "Location"
                    #, guide = FALSE
  ) +
  #6
  labs(x="Phosphorus", y="Detritivores | Insecticides (TUs)") +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 13, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 13, face="bold", color="black"),
        axis.text.y=element_text(size = 13, color="black"),
        axis.text.x=element_text(size = 13, color="black"))


head(Data_final)

##***********************************************************************************************************
## Decomposition | Detritivores ~ %Wastewater
plot(Kdd_alder~Food_3, data=Data_final)
abline(lm(Kdd_alder~Food_3, data=Data_final),col="red")
lm_kdd_Food3_bs<-summary(lm(Kdd_alder~Food_3, data=Data_final))
lm_kdd_Food3_residual<-lm_kdd_Food3_bs$residuals

## Create model with transformed data for plot
M3 <- lmer(lm_kdd_Food3_residual ~ WW_fraction + (1|Site), 
           data=Data_final,
           lmerControl(optimizer = "Nelder_Mead"),
           REML=T)

## Assess parameter estimates
sjPlot::tab_model(M3)

# Create table with effects (parameters estimates)
effects_WW_fraction <- effects::effect(term= "WW_fraction", mod = M3)
summary(effects_WW_fraction ) #output of what the values are

# Save the effects values as a df:
x_WW_fraction  <- as.data.frame(effects_WW_fraction)

## Check unique
unique(Data_final$Type_1)

## Head
head(Data_final)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Add column
Data_final$lm_kdd_Food3_residual <- lm_kdd_Food3_residual

## Create plot
WW_fraction_plot <- ggplot() + 
  #1
  geom_line(data=x_WW_fraction , aes(x=WW_fraction , y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_WW_fraction , aes(x=WW_fraction, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Data_final, aes(WW_fraction, lm_kdd_Food3_residual, fill=Location), shape=21, size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  #scale_shape_manual(breaks=c("D", "U1", "U2"),
  #                   values=c(15, 17, 16),
  #                   labels = c("Unbuffered", "Buffered", "Forest"),
  #                   name = "Site type"
  #                   #, guide = FALSE
  #) +
  #5
  scale_fill_manual(breaks=c("D", "U1", "U2"),
                    values=c("firebrick","royalblue1", "royalblue4"),
                    labels = c("D", "U1", "U2"),
                    name = "Location"
                    #, guide = FALSE
  ) +
  #6
  labs(x="% Wastewater", y="Decomposition | Detritivores") +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 13, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 13, face="bold", color="black"),
        axis.text.y=element_text(size = 13, color="black"),
        axis.text.x=element_text(size = 13, color="black"))


##***********************************************************************************************************
## Decomposition | %Wastewater ~ Detritivores
plot(Kdd_alder~WW_fraction, data=Data_final)
abline(lm(Kdd_alder~WW_fraction, data=Data_final),col="red")
lm_kdd_WWfract_bs<-summary(lm(Kdd_alder~WW_fraction, data=Data_final))
lm_kdd_WWfract_residual<-lm_kdd_WWfract_bs$residuals

##
head(Data_final)

## Create model with transformed data for plot
M4 <- lmer(lm_kdd_WWfract_residual ~ Food_3 + (1|Site), 
           data=Data_final,
           lmerControl(optimizer = "Nelder_Mead"),
           REML=T)

## Assess parameter estimates
sjPlot::tab_model(M4)

# Create table with effects (parameters estimates)
effects_Food_3 <- effects::effect(term= "Food_3", mod = M4)
summary(effects_Food_3) #output of what the values are

# Save the effects values as a df:
x_Food_3 <- as.data.frame(effects_Food_3)

## Check unique
unique(Data_final$Type_1)

## Head
head(Data_final)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Add column
Data_final$lm_kdd_WWfract_residual <- lm_kdd_WWfract_residual

## Create plot
Food_3_plot <- ggplot() + 
  #1
  geom_line(data=x_Food_3, aes(x=Food_3, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_Food_3, aes(x=Food_3, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Data_final, aes(Food_3,lm_kdd_WWfract_residual, fill=Location), shape=21, size=4, alpha=0.8) + 
  #4
  #geom_point(data=x_ripar, aes(x=PC1, y=fit), color="Grey20") +
  #scale_shape_manual(breaks=c("D", "U1", "U2"),
  #                   values=c(15, 17, 16),
  #                   labels = c("Unbuffered", "Buffered", "Forest"),
  #                   name = "Site type"
  #                   #, guide = FALSE
  #) +
  #5
  scale_fill_manual(breaks=c("D", "U1", "U2"),
                    values=c("firebrick","royalblue1", "royalblue4"),
                    labels = c("D", "U1", "U2"),
                    name = "Location"
                    #, guide = FALSE
  ) +
  #6
  labs(x="Detritivores", y="Decomposition | % Wastewater") +
  #
  theme(legend.position="left", legend.box = "vertical",
        #legend.position = c(0.2, 0.8),
        #panel.border = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.y=element_text(vjust=1.25, size = 13, face="bold", color="black"),
        axis.title.x=element_text(vjust=-0.65, size = 13, face="bold", color="black"),
        axis.text.y=element_text(size = 13, color="black"),
        axis.text.x=element_text(size = 13, color="black"))

#dis4_plot

## Arrange using grid plot
library(ggpubr)
ggarrange(SRP_plot, Insecticide_TUs_plot, WW_fraction_plot, Food_3_plot, common.legend = TRUE, legend="right", ncol = 2,
          nrow = 2)

library(grid)
library(gridExtra)
grid.text("(a)", x = unit(0.015, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(b)", x = unit(0.48, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(c)", x = unit(0.015, "npc"), y = unit(0.48, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(d)", x = unit(0.48, "npc"), y = unit(0.48, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))

dev.off()

