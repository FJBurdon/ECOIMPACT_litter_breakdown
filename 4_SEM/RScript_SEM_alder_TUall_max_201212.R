#############################################################################################
##                                                                                         ##  
##                            STRUCTURAL EQUATION MODELLING                                ##   
##                                                                                         ##  
#############################################################################################

##  Script: ECOIMPACT leafpack assay data analysis
##  Author: Dr. Francis J. Burdon
##  VERSION: 5.0
##
##  Primary objective: Use structural equation modelling to explore the mechanisms driving 
##  differences in alder litter breakdown at sites. 
##
##  Version objective: using MAXIMUM summed TU values for "all" invertebrates 
##  (see Burdon et al. 2019  STOTEN and Munz et al. 2017 for more information) - note that there 
##  are only three groups (i.e., "Heavy Metals","Pesticides and Biocides", and "Other" - no 
##  insecticides alone as with the other analysis versions 2, 3, and 5)
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

## Forward Selection: Blanchet et al. (2008)
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

## "maxTU_HM"  "maxTU_Oth"  "maxTU_PB" 

## Exclude metals from non-insecticide TUs
MPs$NonInsect_TU_mu <- NA
MPs$NonInsect_TU_mu <- MPs[,4]-(MPs[,7]+MPs[,12])

## Median
MPs$NonInsect_TU_med <- NA
MPs$NonInsect_TU_med <- MPs[,5]-(MPs[,8]+MPs[,17])

## Max (applied to non-Insecticide columns excluding metals)
MPs$NonInsect_TU_max <- NA
MPs$NonInsect_TU_max <- apply(MPs[,c(10,11,13,14)], 1, FUN=min)

## Now plot to assess the output - mean values
plot(log(NonInsect_TU_mu) ~ log(Total_TU_mu), data=MPs)
abline(0,1, col="blue")

## Median values
plot(log(NonInsect_TU_med) ~ log(Total_TU_med), data=MPs)
abline(0,1, col="blue")

## Maximum values
plot(log(NonInsect_TU_max) ~ log(Total_TU_max), data=MPs)
abline(0,1, col="blue")

## Change Data_all
names(Data_all)[1] <- "Site"

## Merge with data - here we will assess median values for all invertebrates
Data_all <- merge(Data_all, MPs[,c(1:3,38:40)], by=c("Site","Location","Year"))
colnames(Data_all) <- c("Site","Location","Year","Kdd_alder","Gammaridae","Food_3","Feeding_3",
                        "maxTU_HM","maxTU_Oth","maxTU_PB")
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
head(Landuse)

## Merge with data and rename columns (if needed)
Data_all <- merge(Data_all, Landuse[,c(2:4,6)], by=c("Site","Location","Year"))

names(Data_all)[15] <- "Cropping"

## HABITAT (SEDIMENT)

## Turn tibble in data frame
Habitat <- as.data.frame(Habitat)

## Merge with data and rename columns (if needed)
Data_all <- merge(Data_all, Habitat[,c(1,3:6)], by=c("Site","Location","Year"))

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all

colnames(Data_final)

min(Data_final$maxTU_HM)
min(Data_final$maxTU_Oth)
min(Data_final$maxTU_PB)

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

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

## Note PB and Organic MPs are virtually identical

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 1 CWM FOOD                              ##
##                                                                                         ##  
#############################################################################################

## Assess hypothetical predictors to help guide model selection
## Fit global model for estimation of adjusted R2
## Exclude metals since only measured in 2014

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

#head(Data_final)

Gmm <- model.matrix( ~ Food_3 + maxTU_Oth + maxTU_PB  + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(6,9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2241596, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
## "Food_3"  

## SHREDDERS

Shredders <- Data_final$Food_3

Gmm <- model.matrix( ~ maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final)
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Negative variation
vif.cca(rda_result)

global <- Data_final[,c(9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.06764464, alpha = 0.3)
as.matrix(M1$variables) ## Only selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "SRP"     
## "maxTU_PB"

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 2 CWM FEEDING                           ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Drop  + NonInsecticide_TUs
Gmm <- model.matrix( ~ Feeding_3 + maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

head(Data_final)
global <- Data_final[,c(7,9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2255398, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
## "Feeding_3"

## SHREDDERS

Shredders <- Data_final$Feeding_3

Gmm <- model.matrix( ~ maxTU_Oth + maxTU_PB  + SRP + DIN + NH4 + TSS + pORG, Data_final)
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Low variation but not negative
vif.cca(rda_result)

global <- Data_final[,c(9:13,16,17)]
head(global)

require(packfor)
## Blanchet et al. (2008)
## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1376445, alpha = 0.3)
as.matrix(M1$variables) ## Only selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "maxTU_PB"
## "SRP"     

#############################################################################################
##                                                                                         ##  
##                    FORWARD SELECTION: MODEL 3 GAMMARID ABUNDANCES                       ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Gammaridae + maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(5,9:13,16,17)]
head(global)

## Blanchet et al. (2008)
## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.1854535, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Gammaridae"     

## SHREDDERS

Shredders <- Data_final$Gammaridae

Gmm <- model.matrix( ~  maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final)
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Negative vartation
vif.cca(rda_result)

global <- Data_final[,c(9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1780834, alpha = 0.3)
as.matrix(M1$variables)## Only selected using unadjusted R-sqr and alpha = 0.3
## "SRP"     
## "maxTU_PB"
## "DIN"   

## Use same selection criteria as used for the litter breakdown rates
#M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.06744074, alpha = 0.05)
#as.matrix(M1$variables)
## No variables selected using adjusted R-sqr and alpha = 0.05

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

## Standardisation
#Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1 CWM FOOD                                 ##
##                                                                                         ##  
#############################################################################################

colnames(Data_final)

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
## "Food_3"  

## CWM trait: Food_3
## "SRP"     
## "maxTU_PB"

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Food_3 + SRP + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ maxTU_PB + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Pesticides response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                             random=~1|Year/Site, 
                             data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                                random=~1|Year/Site, 
                                data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1 
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth
                      ,DIN %~~% maxTU_Oth
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
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
## With correlated error terms: BIC = 165.397
## Without correlated error terms: BIC = 232.709

## Drop maxTU_Oth ~ Cropping: BIC = 161.932
## Drop maxTU_PB ~ Cropping: BIC = 170.388
## Drop SRP ~ Cropping: BIC = 158.984

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#
#))
#write.csv(out$weights, "SEM_Food3_M1_median_TUall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 2 CWM FEEDING                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
## "Feeding_3"

## Feeding_3
## "maxTU_PB"
## "SRP"     

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response  
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + SRP,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Feeding_3 ~  maxTU_PB,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Pesticides response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Year/Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Year/Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2 
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth 
                      ,DIN %~~% maxTU_Oth
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
                      ,DIN %~~% SRP
)

# Model adequacy
BIC(LPA_pSEM_M2)

# AICc
AIC(LPA_pSEM_M2,aicc = T)

# Run goodness-of-fit tests
summary(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Feeding3_M2_nonInsect_All_CEall_path_estimates.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M2)
#GOF_estimates <- rsquared(LPA_pSEM_M2)
#write.csv(GOF_estimates , "SEM_Feeding3_M2_nonInsect_All_CEall_rsquared.csv",row.names = T)

## Correlated error terms
## Without: BIC = 230.606
## With: BIC = 164.407

## Drop  maxTU_Oth <- Cropping: BIC = 160.939
## Drop  maxTU_PB <- Cropping: BIC = 169.441
## Drop  SRP <- Cropping: BIC = 157.987

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Feeding3_M2_nonInsect_All_CEall_BIC_weights.csv",row.names = T)


#############################################################################################
##                                                                                         ##  
##                            SEM: MODEL 3 GAMMARID ABUNDANCES                             ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
# "Gammaridae"    

## Gammarid
# "SRP"  
# "maxTU_PB"

head(Data_final)

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Gammaridae + SRP,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting shredder response 
  Shredder_response = lme(Gammaridae ~ maxTU_PB + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Pesticides response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Year/Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Year/Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth
                      ,DIN %~~% maxTU_Oth 
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
                      ,DIN %~~% SRP
)

# Model adequacy
BIC(LPA_pSEM_M3)

# AICc
AIC(LPA_pSEM_M3,aicc = T)

# Run goodness-of-fit tests
summary(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
#coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
#Path_estimates <- coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
#write.csv(Path_estimates, "SEM_Gammarid_M3_nonInsect_All_CEall_path_est.csv",row.names = T)

# Explore individual model fits
#rsquared(LPA_pSEM_M3)
#GOF_estimates <- rsquared(LPA_pSEM_M3)
#write.csv(GOF_estimates , "SEM_Gammarid_M3_nonInsect_All_CEall_rsquared.csv",row.names = T)

## Correlated error terms
## Without: BIC = 238.719
## With: BIC = 172.52

## Drop: maxTU_Oth <- Cropping: BIC = 169.183
## Drop: maxTU_PB <- Cropping: BIC = 177.405
## Drop: SRP <- Cropping: BIC = 166.332

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Gammarid_M3_nonInsect_All_CEall_BIC_weights.csv",row.names = T)


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
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all[-which(Data_all$Location=="U2"),]

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

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

## Note PB and Organic MPs are virtually identical

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 1 CWM FOOD                              ##
##                                                                                         ##  
#############################################################################################

## Assess hypothetical predictors to help guide model selection
## Fit global model for estimation of adjusted R2

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

#head(Data_final)

Gmm <- model.matrix( ~ Food_3 + maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(6,9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2741676, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Food_3"     

## SHREDDERS

Shredders <- Data_final$Food_3

Gmm <- model.matrix( ~ maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final)
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Negative variation
vif.cca(rda_result)

global <- Data_final[,c(9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1222631, alpha = 0.3)
as.matrix(M1$variables) ## Only selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "SRP"     
## "maxTU_PB"
## "DIN"    

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 2 CWM FEEDING                           ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Drop  + NonInsecticide_TUs
Gmm <- model.matrix( ~ Feeding_3 + maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(7,9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2486649, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Feeding_3"     

## SHREDDERS

Shredders <- Data_final$Feeding_3

Gmm <- model.matrix( ~ maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Negative variation
vif.cca(rda_result)

global <- Data_final[,c(9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1556532, alpha = 0.3)
as.matrix(M1$variables)  ## Only selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "maxTU_PB"
## "SRP"     
## "DIN" 

#############################################################################################
##                                                                                         ##  
##                    FORWARD SELECTION: MODEL 3 GAMMARID ABUNDANCES                       ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Gammaridae + maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(5,9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.2498132, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Gammaridae"     

## SHREDDERS

Shredders <- Data_final$Gammaridae

Gmm <- model.matrix( ~  maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(9:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because of negative variation for above
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.2330852, alpha = 0.3)
as.matrix(M1$variables)  ## Only selected using unadjusted R-sqr and alpha = 0.3
# "NH4"

## Note that selection criteria is relaxed because of negative variation for above
#M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.06532261, alpha = 0.05)
#as.matrix(M1$variables)
## No variable selected using adjusted R-sqr and alpha = 0.05

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all[-which(Data_all$Location=="U2"),]

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
#Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1 CWM FOOD                                 ##
##                                                                                         ##  
#############################################################################################

colnames(Data_final)

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
## Food_3

## Food_3
## "SRP"     
## "maxTU_PB"
## "DIN"    

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response + WW_fraction
  Leaf_response = lme(Kdd_alder ~ Food_3 + WW_fraction,
                      random=~1|Year/Site,
                      data = Data_final),
  
  # Predicting shredder response
  Shredder_response = lme(Food_3 ~ maxTU_PB + SRP,
                          random=~1|Year/Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Year/Site, 
                           data = Data_final ),
  
  # Predicting Metals response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Year/Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth 
                      ,DIN %~~% maxTU_Oth
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
                      ,DIN %~~% SRP
)

# Model adequacy
BIC(LPA_pSEM_M1)

# AICc
AIC(LPA_pSEM_M1,aicc = T)

# Run goodness-of-fit tests
summary(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M1, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Food3_TUall_max_M1b_path_estimates.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M1)
GOF_estimates <- rsquared(LPA_pSEM_M1)
write.csv(GOF_estimates , "SEM_Food3_TUall_max_M1b_rsquared.csv",row.names = T)

## INITIAL MODEL (GUIDED BY FORWARD-SELECTION)
## With correlated error terms: BIC = 145.847
## Without correlated error terms: BIC = 197.97

## Drop maxTU_HM ~ Cropping: BIC = 143.47
## Drop maxTU_PB ~ Cropping: BIC = 151.13
## Drop SRP ~ Cropping: BIC = 141.048

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#
#))
#write.csv(out$weights, "SEM_Food3_M1_median_TUall_BIC_weights.csv",row.names = T)


#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 2 CWM FEEDING                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
## Feeding_3

## CWM traits Feeding_3
## "maxTU_PB"
## "SRP"     
## "DIN" 

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response 
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  #Predicting Shredder response
  Shredder_response = lme(Feeding_3 ~ maxTU_PB + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Pesticides response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Year/Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction ,
                       random=~1|Year/Site, 
                       data = Data_final),
  
  # Predicting Nitrgen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth 
                      ,DIN %~~% maxTU_Oth 
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
                      ,DIN %~~% SRP
                      )

# Model adequacy
BIC(LPA_pSEM_M2)

# AICc
AIC(LPA_pSEM_M2,aicc = T)

# Run goodness-of-fit tests
summary(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M2, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Feeding3_TUall_max_M2b_path_estimates.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M2)
GOF_estimates <- rsquared(LPA_pSEM_M2)
write.csv(GOF_estimates , "SEM_Feeding3_TUall_max_M2b_rsquared.csv",row.names = T)

## Correlated error terms
## Without: BIC = 196.73
## With: BIC = 145.463

## Drop  maxTU_Oth <- Cropping: BIC = 143.121 
## Drop  maxTU_PB <- Cropping: BIC = 150.816
## Drop  SRP <- Cropping: BIC = 140.7

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Feeding3_M2_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                            SEM: MODEL 3 GAMMARID ABUNDANCES                             ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
## Feeding_3

## CWM traits Gammarid
## NH4 

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Gammaridae + WW_fraction,
                      random=~1|Year/Site, 
                      data = Data_final),
  
  # Predicting Gammarid response + NH4
  Shredder_response = lme(Gammaridae ~ maxTU_PB + DIN + SRP,
                          random=~1|Year/Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Year/Site, 
                     data = Data_final ),
  
  # Predicting Pestiicide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Year/Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Year/Site, 
                       data = Data_final),
  
  # Predicting NH4 response
  #NH4_response = lme(NH4 ~ WW_fraction ,
  #                     random=~1|Year/Site, 
  #                     data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Year/Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth
                      #,NH4 %~~% maxTU_Oth
                      ,DIN %~~% maxTU_Oth
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
                      #,NH4 %~~% maxTU_PB
                      ,DIN %~~% SRP
                      #,NH4 %~~% SRP
                      #,DIN %~~% NH4
)

# Model adequacy
BIC(LPA_pSEM_M3)

# AICc
AIC(LPA_pSEM_M3,aicc = T)

# Run goodness-of-fit tests
summary(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)

# Evaluate path significance using unstandardized coefficients
coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
Path_estimates <- coefs(LPA_pSEM_M3, standardize = "scale", intercepts = FALSE)
write.csv(Path_estimates, "SEM_Gammarid_TUall_max_M3b_path_est.csv",row.names = T)

# Explore individual model fits
rsquared(LPA_pSEM_M3)
GOF_estimates <- rsquared(LPA_pSEM_M3)
write.csv(GOF_estimates , "SEM_Gammarid_TUall_max_M3b_rsquared.csv",row.names = T)

## Correlated error terms
## Without: BIC = 249.282
## With: BIC = 179.46

## NH4 ~ Cropping = BIC: 176.086
## Drop maxTU_Oth ~ Cropping = BIC: 173.852
## Drop maxTU_PB ~ Cropping = BIC: 179.869
## Drop SRP ~ Cropping = BIC: 

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Gammarid_M3_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                  RANK INVERTEBRATE RESPONSES FOR VARIATION EXPLAINED                    ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Food_3 explains the most, but variation similar across invert predicto

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

Data_final <- Data_all[-which(Data_all$Location=="U2"),]
Data_final <- Data_final[-which(Data_final$Year=="2013"),]

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

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

## Note PB and Organic MPs are virtually identical

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 1 CWM FOOD                              ##
##                                                                                         ##  
#############################################################################################

## Assess hypothetical predictors to help guide model selection
## Fit global model for estimation of adjusted R2

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

#head(Data_final)

Gmm <- model.matrix( ~ Food_3 + maxTU_HM + maxTU_Oth + maxTU_PB + SRP + DIN +  NH4 + TSS  + pORG, Data_final)
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(6,8:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.1694445, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Food_3"

## SHREDDERS
## See massive correlations between MP predictors

Shredders <- Data_final$Food_3

Gmm <- model.matrix( ~ maxTU_HM + maxTU_Oth + maxTU_PB + SRP + DIN +  NH4 + TSS  + pORG, Data_final)
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result) ## See metals still > 10 VIF

global <- Data_final[,c(8:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.4244017, alpha = 0.3)
as.matrix(M1$variables) ## Only selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "maxTU_PB" 
## "SRP"      
## "maxTU_Oth"

#############################################################################################
##                                                                                         ##  
##                        FORWARD SELECTION: MODEL 2 CWM FEEDING                           ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Feeding_3 + maxTU_HM + maxTU_Oth + maxTU_PB + SRP + DIN +  NH4 + TSS  + pORG, Data_final)
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(7,8:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.1644378, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
# "Feeding_3"    

## SHREDDERS

Shredders <- Data_final$Feeding_3

Gmm <- model.matrix( ~  maxTU_HM + maxTU_Oth + maxTU_PB + SRP + DIN +  NH4 + TSS  + pORG, Data_final)
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Negative variation
vif.cca(rda_result)

global <-  Data_final[,c(8:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.423853, alpha = 0.3)
as.matrix(M1$variables) ## Only selected using unadjusted R-sqr and alpha = 0.3
rm(M1,global)
## "maxTU_PB" 
##  "SRP"      
## "maxTU_Oth"

#############################################################################################
##                                                                                         ##  
##                    FORWARD SELECTION: MODEL 3 GAMMARID ABUNDANCES                       ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

Gmm <- model.matrix( ~ Gammaridae +  maxTU_HM + maxTU_Oth + maxTU_PB + SRP + DIN +  NH4 + TSS  + pORG, Data_final) 
rda_result <- rda(Decomposition ~ Gmm)
RsquareAdj(rda_result)
vif.cca(rda_result)

global <- Data_final[,c(5,8:13,16,17)]
head(global)

## Blanchet et al. (2008)
## Use forward selection with appropriate criteria to assess which variables are essential
M1 <- forward.sel(Decomposition, global, nperm=999, adjR2thresh = 0.4562177, alpha = 0.05)
as.matrix(M1$variables)
rm(M1,global)
## "Gammaridae"

## SHREDDERS

Shredders <- Data_final$Gammaridae

Gmm <- model.matrix( ~  maxTU_HM + maxTU_Oth + maxTU_PB + SRP + DIN + NH4 + TSS + pORG, Data_final) 
rda_result <- rda(Shredders ~ Gmm)
RsquareAdj(rda_result) ## Negative vartation
vif.cca(rda_result)

global <- Data_final[,c(8:13,16,17)]
head(global)

## Use forward selection with appropriate criteria to assess which variables are essential
## Note that selection criteria is relaxed because otherwise no parameters are selected
M1 <- forward.sel(Shredders, global, nperm=999, adjR2thresh = 0.1766575, alpha = 0.3)
as.matrix(M1$variables) ## Only selected using unadjusted R-sqr and alpha = 0.3
## "maxTU_PB"
## "SRP"     

#############################################################################################
##                                                                                         ##  
##                          TRANSFORMATION AND STANDARDISATION                             ##
##                                                                                         ##  
#############################################################################################

Data_final <- Data_all[-which(Data_all$Location=="U2"),]
Data_final <- Data_final[-which(Data_final$Year=="2013"),]

## Transformation
Data_final[,c(5:8)] <- log1p(Data_final[,c(5:8)])# log+1 transform Gammarid abundances etc.
Data_final[,c(4,9:13,16)] <- log(Data_final[,c(4,9:13,16)])# log transform Water chem + Kdd
Data_final[,c(14,15,17)] <- logit(Data_final[,c(14,15,17)]/100, adjust=0.025)# logit transform Cropping, %WW

head(Data_final)

## Standardisation
Data_final[,c(4:17)] <- decostand(Data_final[,c(4:17)],"standardize")

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1 CWM FOOD                                 ##
##                                                                                         ##  
#############################################################################################

colnames(Data_final)

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Kdd
# "Food_3"  

# Food3
## "maxTU_PB" 
## "SRP"      
## "maxTU_Oth"

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Food_3 + DIN + SRP,
                      random=~1|Site,
                      data = Data_final),
  
  #Predicting shredder response 
  Shredder_response = lme(Food_3 ~ maxTU_PB + SRP,
                          random=~1|Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1
                      ,maxTU_Oth %~~% maxTU_PB
                      ,SRP %~~% maxTU_Oth
                      ,DIN %~~% maxTU_Oth 
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
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
## With correlated error terms: BIC = 101.236
## Without correlated error terms: BIC = 135.305

## Drop maxTU_Oth ~ Cropping: BIC = 101.052
## Drop maxTU_PB ~ Cropping: BIC = 101.235
## Drop SRP ~ Cropping: BIC = 101.382

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#
#))
#write.csv(out$weights, "SEM_Food3_M1_median_TUall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 1d CWM FOOD                                ##
##                                                                                         ##  
#############################################################################################

## TEST METALS EXPLICITLY

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M1 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Food_3 + DIN,
                      random=~1|Site,
                      data = Data_final),
  
  #Predicting shredder response 
  Shredder_response = lme(Food_3 ~ maxTU_PB + SRP,
                          random=~1|Site, 
                          data = Data_final ),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Site, 
                           data = Data_final ),
  
  # Predicting Metals response
  Metals_response = lme(maxTU_HM ~ WW_fraction,
                       random=~1|Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M1 <- update(LPA_pSEM_M1
                      ,maxTU_HM %~~% maxTU_PB
                      ,SRP %~~% maxTU_HM
                      ,DIN %~~% maxTU_HM 
                      ,DIN %~~% maxTU_PB
                      ,SRP %~~% maxTU_PB
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

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#
#))
#write.csv(out$weights, "SEM_Food3_M1_median_TUall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 2 CWM FEEDING                              ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Kdd
# "Feeding_3"    

# Feeding_3
## "maxTU_PB" 
## "SRP"      
## "maxTU_Oth"

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response  
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + DIN + SRP,
                      random=~1|Site, 
                      data = Data_final),
  
  #Predicting shredder response + Insecticide_TU_all_med
  Shredder_response = lme(Feeding_3 ~ maxTU_PB + SRP,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2 
                      ,maxTU_PB %~~% maxTU_Oth  
                      ,SRP %~~% maxTU_PB 
                      ,DIN %~~% maxTU_PB
                      ,DIN %~~% maxTU_Oth  
                      ,SRP %~~% maxTU_Oth  
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

## Correlated error terms
## Without: BIC = 136.205
## With: BIC = 101.283

## Drop  maxTU_Oth <- Cropping: BIC = 101.108
## Drop  maxTU_PB <- Cropping: BIC = 101.295
## Drop  SRP <- Cropping: BIC = 101.280

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Feeding3_M2_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                   SEM: MODEL 2d CWM FEEDING                             ##
##                                                                                         ##  
#############################################################################################

## TEST METALS EXPLICITLY

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

# Create component models and store in list
LPA_pSEM_M2 = psem(
  
  # Predicting LPA response  
  Leaf_response = lme(Kdd_alder ~ Feeding_3 + DIN + maxTU_HM,
                      random=~1|Site, 
                      data = Data_final),
  
  #Predicting Shredder response + Insecticide_TU_all_med
  Shredder_response = lme(Feeding_3 ~ maxTU_PB + SRP + maxTU_HM,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Site, 
                           data = Data_final ),
  
  # Predicting Metals response
  Metals_response = lme(maxTU_HM ~ WW_fraction,
                       random=~1|Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M2 <- update(LPA_pSEM_M2 
                      ,maxTU_PB %~~% maxTU_HM  
                      ,SRP %~~% maxTU_PB 
                      ,DIN %~~% maxTU_PB
                      ,DIN %~~% maxTU_HM 
                      ,SRP %~~% maxTU_HM  
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

## Correlated error terms
## Without: BIC = 143.155
## With: BIC = 102.913

## Drop  maxTU_HM <- Cropping: BIC = 100.75
## Drop  maxTU_PB <- Cropping: BIC = 100.938
## Drop  SRP <- Cropping: BIC = 100.922

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Feeding3_M2_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                            SEM: MODEL 3 GAMMARID ABUNDANCES                             ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

## Kdd
## "Gammaridae"          

## Gammaridae
## "maxTU_PB" 
## "SRP"      
## "maxTU_Oth"

head(Data_final)

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Gammaridae + WW_fraction,
                      random=~1|Site, 
                      data = Data_final),
  
  # Predicting Gammarid response maxTU_PB
  Shredder_response = lme(Gammaridae ~ SRP + maxTU_Oth,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction + Cropping,
                           random=~1|Site, 
                           data = Data_final ),
  
  # Predicting Other response
  Other_response = lme(maxTU_Oth ~ WW_fraction,
                       random=~1|Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3 
                      ,maxTU_Oth %~~% maxTU_PB  
                      ,SRP %~~% maxTU_PB  
                      ,DIN %~~% maxTU_PB  
                      ,DIN %~~% maxTU_Oth 
                      ,SRP %~~% maxTU_Oth 
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

## Correlated error terms
## Without: BIC = 134.129
## With: BIC = 99.206

## Drop maxTU_Oth ~ Cropping = BIC: 97.935
## Drop maxTU_PB ~ Cropping = BIC: 98.021
## Drop SRP ~ Cropping = BIC: 97.801

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Gammarid_M3_nonInsect_All_CEall_BIC_weights.csv",row.names = T)



#############################################################################################
##                                                                                         ##  
##                            SEM: MODEL 3d GAMMARID ABUNDANCES                            ##
##                                                                                         ##  
#############################################################################################

# Now fit piecewise model with random effects (Site nested in Year)
# Note: model below is best-fitting identified bv BIC scores

head(Data_final)

# Create component models and store in list
LPA_pSEM_M3 = psem(
  
  # Predicting LPA response
  Leaf_response = lme(Kdd_alder ~ Gammaridae + WW_fraction,
                      random=~1|Site, 
                      data = Data_final),
  
  # Predicting Gammarid response maxTU_PB
  Shredder_response = lme(Gammaridae ~ maxTU_PB + SRP,
                          random=~1|Site, 
                          data = Data_final),
  
  # Predicting Phosphorus response
  SRP_response = lme(SRP ~ WW_fraction,
                     random=~1|Site, 
                     data = Data_final ),
  
  # Predicting Pesticide response
  Pesticide_response = lme(maxTU_PB ~ WW_fraction,
                           random=~1|Site, 
                           data = Data_final ),
  
  # Predicting Metals response
  Metals_response = lme(maxTU_HM ~ WW_fraction,
                       random=~1|Site, 
                       data = Data_final),
  
  # Predicting Nitrogen response
  DIN_response = lme(DIN ~ WW_fraction + Cropping,
                     random=~1|Site, 
                     data = Data_final))

LPA_pSEM_M3 <- update(LPA_pSEM_M3 
                      ,maxTU_HM %~~% maxTU_PB  
                      ,SRP %~~% maxTU_PB  
                      ,DIN %~~% maxTU_PB  
                      ,DIN %~~% maxTU_HM 
                      ,SRP %~~% maxTU_HM 
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

## Correlated error terms
## Without: BIC = 146.028
## With: BIC = 105.574

## Drop maxTU_Oth ~ Cropping = BIC: 103.579
## Drop maxTU_PB ~ Cropping = BIC: 103.329
## Drop SRP ~ Cropping = BIC: 103.195

## See Excel table for model iterations and BIC scores used below for weightings
## Note: models were not considered if the test of directed separation indicated a significant
## unspecified path
#out <- akaike.weights(c(
#))
#write.csv(out$weights, "SEM_Gammarid_M3_nonInsect_All_CEall_BIC_weights.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                  RANK INVERTEBRATE RESPONSES FOR VARIATION EXPLAINED                    ##
##                                                                                         ##  
#############################################################################################

## DECOMPOSITION

Decomposition <- Data_final$Kdd_alder

## Food_3 explains the most, but variation similar across invert predicto

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
png(filename="ECOIMPACT_partial_regression_LME_211129.png", 
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
M1 <- lmer(lm_Food3_SRP_residual ~ maxTU_PB + (1|Site), 
            data=Data_final,
            lmerControl(optimizer = "Nelder_Mead"),
            REML=T)

## Assess parameter estimates
sjPlot::tab_model(M1)

# Create table with effects (parameters estimates)
effects_maxTU_PB <- effects::effect(term= "maxTU_PB", mod = M1)
summary(effects_maxTU_PB) #output of what the values are

# Save the effects values as a df:
x_maxTU_PB <- as.data.frame(effects_maxTU_PB)

## Check unique
unique(Data_final$Type_1)

## Head
head(Data_final)

## Set theme to classic
sjPlot::set_theme(base = theme_classic())

## Add column
Data_final$lm_Food3_SRP_residual <- lm_Food3_SRP_residual


## Create plot
maxTU_PB_plot <- ggplot() + 
  #1
  geom_line(data=x_maxTU_PB, aes(x=maxTU_PB, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_maxTU_PB, aes(x=maxTU_PB, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Data_final, aes(maxTU_PB, lm_Food3_SRP_residual, fill=Location), shape=21, size=4, alpha=0.8) + 
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
plot(Food_3~maxTU_PB, data=Data_final)
abline(lm(Food_3~maxTU_PB, data=Data_final),col="red")
lm_Food3_maxTU_PB_bs<-summary(lm(Food_3~maxTU_PB, data=Data_final))
lm_Food3_maxTU_PB_residual<-lm_Food3_maxTU_PB_bs$residuals

## Create model with transformed data for plot
M2 <- lmer(lm_Food3_maxTU_PB_residual ~ SRP + (1|Site), 
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
Data_final$lm_Food3_maxTU_PB_residual <- lm_Food3_maxTU_PB_residual

## Create plot
SRP_plot <- ggplot() + 
  #1
  geom_line(data=x_SRP, aes(x=SRP, y=fit), color="Grey10", linetype=1) +
  #2
  geom_ribbon(data= x_SRP, aes(x=SRP, ymin=lower, ymax=upper), alpha= 0.1, fill="Grey10") +
  #3
  geom_point(data=Data_final, aes(SRP, lm_Food3_maxTU_PB_residual, fill=Location), shape=21, size=4, alpha=0.8) + 
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
ggarrange(SRP_plot, maxTU_PB_plot, WW_fraction_plot, Food_3_plot, common.legend = TRUE, legend="right", ncol = 2,
           nrow = 2)

library(grid)
library(gridExtra)
grid.text("(a)", x = unit(0.02, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(b)", x = unit(0.48, "npc"), y = unit(0.98, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(c)", x = unit(0.02, "npc"), y = unit(0.48, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))
grid.text("(d)", x = unit(0.48, "npc"), y = unit(0.48, "npc"),gp=gpar(fontsize=16, fontface="bold", col="black"))

dev.off()



