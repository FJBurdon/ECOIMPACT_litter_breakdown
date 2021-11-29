#############################################################################################
##                                                                                         ##  
##                    IMPACTS OF WW: MIXED MODELS AND EFFECT SIZES                         ##  
##                                      ALDER                                              ##   
##                                                                                         ##  
#############################################################################################

##  Script: ECOIMPACT leafpack assay data analysis
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Calculate mean +/- SD for alder breakdown
##  2. Test changes in alder breakdown using linear mixed-models
##  3. Calculate effect sizes for alder breakdown

##***************************************************************************************
## Load libraries you need

##Loads package used for Logit transformation (% data)
require(car)
##Loads vegan package for transformation
library(vegan)
##Loads reshape and plyr package for data manipulation
library(reshape)
library(reshape2)
library(plyr)
## Load data files
require(readr)
## To calculate effect sizes
require(SingleCaseES)
## Calculate mixed models
require(lme4)
require(lmerTest)
require(blme)
require(emmeans)
## Load tidyverse
require(tidyverse)

##***************************************************************************************
## Load data

## All years coarse mesh decomposition
Alder_allyrs <-  read.csv(here::here("1_Alder_data_input","lpa_alder_coarse_All_data.csv")) 

## 2013 coarse mesh decomposition
Alder_2013 <- read.csv(here::here("1_Alder_data_input","lpa_alder_coarse_2013_data.csv")) 

## 2013 Microbe and Invertebrate decomposition
Microbial_breakdown_2013 <- read.csv(here::here("1_Alder_data_input","lpa_alder_fine_2013_data.csv")) 
Invert_breakdown_2013 <- read.csv(here::here("1_Alder_data_input","lpa_alder_invertebrate_2013_data.csv")) 

#############################################################################################
##                                                                                         ##  
##                      INVERTEBRATE EFFECT SIZES - 2013  ALDER                            ##   
##                                                                                         ##  
#############################################################################################

## Check input data
head(Invert_breakdown_2013)

## Test invertebrate-mediated breakdown k
k_invert_day_ml_M1 <- blmer(sqrt(k_invert_day) ~ Location + (1|Site_code), data = Invert_breakdown_2013, 
                            REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_invert_day_ml_M1)

## Check output
summary(k_invert_day_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_invert_day_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_invert_day_ml_M1), lower.tail = F)
k_invert_day_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_invert_day_ml_M1_coefficients,"K_day_alder_invert_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_invert_day_ml_M1_conf_intervals <- confint(k_invert_day_ml_M1, level = 0.95)
write.csv(k_invert_day_ml_M1_conf_intervals,"K_day_alder_invert_2013_LME_confint.csv",row.names = T)

## DF 2, 182.06
## Write ANOVA table
k_invert_day_ml_M1_anova <- anova(k_invert_day_ml_M1) 
k_invert_day_ml_M1_anova$p_value <- NA
k_invert_day_ml_M1_anova$p_value[1] <- pf(k_invert_day_ml_M1_anova[1,4], k_invert_day_ml_M1_anova[1,1], df.residual(k_invert_day_ml_M1), lower.tail = FALSE)
write.csv(k_invert_day_ml_M1_anova,"K_day_alder_invert_2013_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_invert_day_ml_M1)
n <- nrow(X)
Beta <- fixef(k_invert_day_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_invert_day_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0252
## Random 0.6591
## All 0.6843
write.csv(Model_variance,"K_day_alder_invert_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_invert_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_alder_invert_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_alder 2013
K_day_alder_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_day"), c(mean), na.rm=T))
K_day_alder_mu_invert_location <- data.frame(K_day_alder_mu_invert_site  %>%
                                               group_by(Location) %>% 
                                               summarise_at(vars("k_invert_day"), c(mean,sd), na.rm=T))
colnames(K_day_alder_mu_invert_location) <- c("Location","Mean","Stdev")
write.csv(K_day_alder_mu_invert_location,"K_day_alder_invert_2013_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate model for invertebrate K degree-day
k_invert_dd_ml_M1 <- blmer(sqrt(k_invert_dd) ~ Location + (1|Site_code), data = Invert_breakdown_2013, 
                           REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_invert_dd_ml_M1)

## Check output
summary(k_invert_dd_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_invert_dd_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_invert_dd_ml_M1), lower.tail = F)
k_invert_dd_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_invert_dd_ml_M1_coefficients,"K_degreeday_alder_invert_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_invert_dd_ml_M1_conf_intervals <- confint(k_invert_dd_ml_M1, level = 0.95)
write.csv(k_invert_dd_ml_M1_conf_intervals,"K_degreeday_alder_invert_2013_LME_confint.csv",row.names = T)

## DF 2, 182.08
## Write ANOVA table
k_invert_dd_ml_M1_anova <- anova(k_invert_dd_ml_M1) 
k_invert_dd_ml_M1_anova$p_value <- NA
k_invert_dd_ml_M1_anova$p_value[1] <- pf(k_invert_dd_ml_M1_anova[1,4], k_invert_dd_ml_M1_anova[1,1], df.residual(k_invert_dd_ml_M1), lower.tail = FALSE)
write.csv(k_invert_dd_ml_M1_anova,"K_degreeday_alder_invert_2013_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_invert_dd_ml_M1)
n <- nrow(X)
Beta <- fixef(k_invert_dd_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_invert_dd_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0363
## Random 0.5883
## All 0.6245
write.csv(Model_variance,"K_degreeday_alder_invert_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_invert_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_alder_invert_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_alder ALL
K_day_alder_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_dd"), c(mean), na.rm=T))
K_day_alder_mu_invert_location <- data.frame(K_day_alder_mu_invert_site  %>%
                                               group_by(Location) %>% 
                                               summarise_at(vars("k_invert_dd"), c(mean,sd), na.rm=T))
colnames(K_day_alder_mu_invert_location) <- c("Location","Mean","Stdev")
write.csv(K_day_alder_mu_invert_location,"K_degreeday_alder_invert_2013_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate effect sizes for invertebrate-mediated breakdown
Invert_2013 <- Invert_breakdown_2013
head(Invert_2013)
## Invertebrate k day 2013
Invertebrate_D_U1_kday <-  batch_calc_ES(dat = Invert_2013[-which(Invert_2013$Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_day, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_D_U2_kday  <-  batch_calc_ES(dat = Invert_2013[-which(Invert_2013$Location=="US1"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_day, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_U1_U2_kday  <-  batch_calc_ES(dat = Invert_2013[-which(Invert_2013$Location=="DS"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_invert_day, 
                                          improvement = "decrease",
                                          ES="SMD")

## Invertebrate k degree-day 2013
Invertebrate_D_U1_kdd  <-  batch_calc_ES(dat = Invert_2013[-which(Invert_2013$Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_dd, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_D_U2_kdd <-  batch_calc_ES(dat = Invert_2013[-which(Invert_2013$Location=="US1"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_dd, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_U1_U2_kdd <-  batch_calc_ES(dat = Invert_2013[-which(Invert_2013$Location=="D"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_invert_dd, 
                                          improvement = "decrease",
                                          ES="SMD")

## Bind output
Invert_kday <- rbind(Invertebrate_D_U1_kday,Invertebrate_D_U2_kday,Invertebrate_U1_U2_kday)
Invert_kdd  <- rbind(Invertebrate_D_U1_kdd,Invertebrate_D_U2_kdd,Invertebrate_U1_U2_kdd)
Invert_out  <- rbind(Invert_kday,Invert_kdd)
Invert_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)

## Calculate site means for alternative effect size
K_day_alder_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_day"), c(mean), na.rm=T))

K_dday_alder_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_dd"), c(mean), na.rm=T))

## Invertebrate k day 2013
Invertebrate_D_U1_kday <-  batch_calc_ES(dat = K_day_alder_mu_invert_site [-which(K_day_alder_mu_invert_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_day, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_D_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_invert_site [-which(K_day_alder_mu_invert_site $Location=="US1"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_invert_day, 
                                          improvement = "decrease",
                                          ES="SMD")

Invertebrate_U1_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_invert_site [-which(K_day_alder_mu_invert_site $Location=="DS"),],
                                           #grouping = Site, 
                                           condition = Location,
                                           outcome = k_invert_day, 
                                           improvement = "decrease",
                                           ES="SMD")

## Invertebrate k degree-day 2013
Invertebrate_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_alder_mu_invert_site [-which(K_dday_alder_mu_invert_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_dd, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_D_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_invert_site [-which(K_dday_alder_mu_invert_site $Location=="US1"),],
                                        #grouping = Site, 
                                        condition = Location,
                                        outcome = k_invert_dd, 
                                        improvement = "decrease",
                                        ES="SMD")

Invertebrate_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_invert_site [-which(K_dday_alder_mu_invert_site $Location=="DS"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_dd, 
                                         improvement = "decrease",
                                         ES="SMD")

## Bind output
Invert_kday <- rbind(Invertebrate_D_U1_kday,Invertebrate_D_U2_kday,Invertebrate_U1_U2_kday)
Invert_kdd  <- rbind(Invertebrate_D_U1_kdd,Invertebrate_D_U2_kdd,Invertebrate_U1_U2_kdd)
Invert_out  <- rbind(Invert_kday,Invert_kdd)
Invert_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)
Invert_out_mu  <- Invert_out

#############################################################################################
##                                                                                         ##  
##                        MICROBIAL EFFECT SIZES - 2013 ALDER                              ##   
##                                                                                         ##  
#############################################################################################

## Check input data
head(Microbial_breakdown_2013)

## Outliers
Alder_microbe_outliers_2013 <- Microbial_breakdown_2013[c(96,106,167),]
write.csv(Alder_microbe_outliers_2013,"Alder_microbe_2013_outliers.csv",row.names = F)

## Test microbially-mediated breakdown k
#k_microbe_day_ml_M1 <- blmer(sqrt(k_day_ml) ~ Location + (1|Site_code), data = Microbial_breakdown_2013[-c(96,106),], 
#                             REML=T, control = lmerControl(optimizer = "nloptwrap"))
k_microbe_day_ml_M1 <- blmer(sqrt(k_day_ml) ~ Location + (1|Site_code), data = Microbial_breakdown_2013[-c(96,106,167),], 
                             REML=T, control = lmerControl(optimizer = "nloptwrap"))
## Check residuals
plot(k_microbe_day_ml_M1)
Microbial_breakdown_2013[c(96,106,167),c(1:8,14)]## Human error - look to be larger than expected
## Year Site_code Location Sample_location Mesh  Leaf Person Rep      k_ml
## 96  2013       HOR      US1            2_U1 fine alder      S   4 0.5057105
## 106 2013       KER       DS             3_D fine alder      S   6 0.2154807
## 167 2013       ROM      US2            1_U2 fine alder      S   6 0.4855078

## Check output
summary(k_microbe_day_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_microbe_day_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_microbe_day_ml_M1), lower.tail = F)
k_microbe_day_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_microbe_day_ml_M1_coefficients,"K_day_alder_microbe_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_microbe_day_ml_M1_conf_intervals <- confint(k_microbe_day_ml_M1, level = 0.95)
write.csv(k_microbe_day_ml_M1_conf_intervals,"K_day_alder_microbe_2013_LME_confint.csv",row.names = T)

## DF 2, 186.06
## Write ANOVA table
k_microbe_day_ml_M1_anova <- anova(k_microbe_day_ml_M1) 
k_microbe_day_ml_M1_anova$p_value <- NA
k_microbe_day_ml_M1_anova$p_value[1] <- pf(k_microbe_day_ml_M1_anova[1,4], k_microbe_day_ml_M1_anova[1,1], df.residual(k_microbe_day_ml_M1), lower.tail = FALSE)
write.csv(k_microbe_day_ml_M1_anova,"K_day_alder_microbe_2013_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_microbe_day_ml_M1)
n <- nrow(X)
Beta <- fixef(k_microbe_day_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_microbe_day_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0190
## Random 0.5801
## All 0.5991
write.csv(Model_variance,"K_day_alder_microbe_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_microbe_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_alder_microbe_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_alder microbe
K_day_alder_mu_microbe_site <- data.frame(Microbial_breakdown_2013[-c(96,106,167),] %>%
                                            group_by(Site_code,Location) %>% 
                                            summarise_at(vars("k_day_ml"), c(mean), na.rm=T))
K_day_alder_mu_microbe_location <- data.frame(K_day_alder_mu_microbe_site  %>%
                                                group_by(Location) %>% 
                                                summarise_at(vars("k_day_ml"), c(mean,sd), na.rm=T))
colnames(K_day_alder_mu_microbe_location) <- c("Location","Mean","Stdev")
write.csv(K_day_alder_mu_microbe_location,"K_day_alder_microbe_2013_mean_std_location.csv",row.names = F)

##*****************************************************************************/
## Calculate model for microbial K degree day
k_microbe_dd_ml_M1 <- blmer(sqrt(k_dd_ml) ~ Location + (1|Site_code), data = Microbial_breakdown_2013[-c(96,106,167),], 
                            REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_microbe_dd_ml_M1)
which.max(residuals(k_microbe_dd_ml_M1))
Microbial_breakdown_2013[167,c(1:7,14)]
##    Year Site_code Location Sample_location Mesh  Leaf Person          k_ml
##167 2013       ROM      US2            1_U2 fine alder      S     0.4855078
## High value, probable human error

## Check output
summary(k_microbe_dd_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_microbe_dd_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_microbe_dd_ml_M1), lower.tail = F)
k_microbe_dd_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_microbe_dd_ml_M1_coefficients,"K_degreeday_alder_microbe_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_microbe_dd_ml_M1_conf_intervals <- confint(k_microbe_dd_ml_M1, level = 0.95)
write.csv(k_microbe_dd_ml_M1_conf_intervals,"K_degreeday_alder_microbe_2013_LME_confint.csv",row.names = T)

## DF 2, 186.41
## Write ANOVA table
k_microbe_dd_ml_M1_anova <- anova(k_microbe_dd_ml_M1) 
k_microbe_dd_ml_M1_anova$p_value <- NA
k_microbe_dd_ml_M1_anova$p_value[1] <- pf(k_microbe_dd_ml_M1_anova[1,4], k_microbe_dd_ml_M1_anova[1,1], df.residual(k_microbe_dd_ml_M1), lower.tail = FALSE)
write.csv(k_microbe_dd_ml_M1_anova,"K_degreeday_alder_microbe_2013_LME_ANOVA_200816.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_microbe_dd_ml_M1)
n <- nrow(X)
Beta <- fixef(k_microbe_dd_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_microbe_dd_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0001241
## Random 0.2596
## All 0.2597
write.csv(Model_variance,"K_degreeday_alder_microbe_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_microbe_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_alder_microbe_2013_LME_post_Location.csv",row.names = F)

## Calculate means for microbially-mediated K degree-day alder
K_day_alder_mu_microbe_site <- data.frame(Microbial_breakdown_2013[-c(96,106,167),] %>%
                                            group_by(Site_code,Location) %>% 
                                            summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))
K_day_alder_mu_microbe_location <- data.frame(K_day_alder_mu_microbe_site  %>%
                                                group_by(Location) %>% 
                                                summarise_at(vars("k_dd_ml"), c(mean,sd), na.rm=T))
colnames(K_day_alder_mu_microbe_location) <- c("Location","Mean","Stdev")
write.csv(K_day_alder_mu_microbe_location,"K_degreeday_alder_microbe_2013_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate effect sizes for microbially-mediated breakdown
Microbial_2013 <- Microbial_breakdown_2013[-c(96,106,167),]
## Microbe k day 2013
Microbe_D_U1_kday <-  batch_calc_ES(dat = Microbial_2013[-which(Microbial_2013$Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_day_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

Microbe_D_U2_kday  <-  batch_calc_ES(dat = Microbial_2013[-which(Microbial_2013$Location=="US1"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_day_ml, 
                                          improvement = "decrease",
                                          ES="SMD")

Microbe_U1_U2_kday  <-  batch_calc_ES(dat = Microbial_2013[-which(Microbial_2013$Location=="DS"),],
                                           #grouping = Site, 
                                           condition = Location,
                                           outcome = k_day_ml, 
                                           improvement = "decrease",
                                           ES="SMD")

## Microbe k degree-day 2013
Microbe_D_U1_kdd  <-  batch_calc_ES(dat = Microbial_2013[-which(Microbial_2013$Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_dd_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

Microbe_D_U2_kdd <-  batch_calc_ES(dat = Microbial_2013[-which(Microbial_2013$Location=="US1"),],
                                        #grouping = Site, 
                                        condition = Location,
                                        outcome = k_dd_ml, 
                                        improvement = "decrease",
                                        ES="SMD")

Microbe_U1_U2_kdd <-  batch_calc_ES(dat = Microbial_2013[-which(Microbial_2013$Location=="DS"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_dd_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

## Bind output
Microbe_kday <- rbind(Microbe_D_U1_kday,Microbe_D_U2_kday,Microbe_U1_U2_kday)
Microbe_kdd  <- rbind(Microbe_D_U1_kdd,Microbe_D_U2_kdd,Microbe_U1_U2_kdd)
Microbe_out  <- rbind(Microbe_kday,Microbe_kdd)
Microbe_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)

## Calculate site means for alternative effect size
K_day_alder_mu_microbe_site <- data.frame(Microbial_breakdown_2013[-c(96,106,167),] %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_day_ml"), c(mean), na.rm=T))

K_dday_alder_mu_microbe_site <- data.frame(Microbial_breakdown_2013 %>%
                                            group_by(Site_code,Location) %>% 
                                            summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))

## Microbe k day 2013
microbe_D_U1_kday <-  batch_calc_ES(dat = K_day_alder_mu_microbe_site [-which(K_day_alder_mu_microbe_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_day_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

microbe_D_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_microbe_site [-which(K_day_alder_mu_microbe_site $Location=="US1"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_day_ml, 
                                          improvement = "decrease",
                                          ES="SMD")

microbe_U1_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_microbe_site [-which(K_day_alder_mu_microbe_site $Location=="DS"),],
                                           #grouping = Site, 
                                           condition = Location,
                                           outcome = k_day_ml, 
                                           improvement = "decrease",
                                           ES="SMD")

## Microbe k degree-day 2013
microbe_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_alder_mu_microbe_site [-which(K_dday_alder_mu_microbe_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_dd_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

microbe_D_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_microbe_site [-which(K_dday_alder_mu_microbe_site $Location=="US1"),],
                                        #grouping = Site, 
                                        condition = Location,
                                        outcome = k_dd_ml, 
                                        improvement = "decrease",
                                        ES="SMD")

microbe_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_microbe_site [-which(K_dday_alder_mu_microbe_site $Location=="DS"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_dd_ml, 
                                         improvement = "decrease",
                                         ES="SMD")
## Bind output
microbe_kday <- rbind(microbe_D_U1_kday,microbe_D_U2_kday,microbe_U1_U2_kday)
microbe_kdd  <- rbind(microbe_D_U1_kdd,microbe_D_U2_kdd,microbe_U1_U2_kdd)
microbe_out  <- rbind(microbe_kday,microbe_kdd)
microbe_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)
microbe_out_mu  <- microbe_out

#############################################################################################
##                                                                                         ##  
##                            k ratio (invert/microbe)                                     ##   
##                                                                                         ##  
#############################################################################################

## Calculate K ratio
Invert_breakdown_2013$k_ratio <- Invert_breakdown_2013$k_invert/Invert_breakdown_2013$k_microbe

## Remove incomplete cases
Invert_breakdown_2013 <- Invert_breakdown_2013[complete.cases(Invert_breakdown_2013), ]

## Ensure row numbers are chronological
row.names(Invert_breakdown_2013) <- NULL

##*****************************************************************************/
## Calculate model for K ratio
#K_ratio_M1 <- lmer(sqrt(k_ratio) ~ Location + (1|Site_code), data = Invert_breakdown_2013, REML=T)
## remove microbial outliers where appropriate for consistency
K_ratio_M1 <- lmer(sqrt(k_ratio) ~ Location + (1|Site_code), data = Invert_breakdown_2013[-c(34,91,105,162),], REML=T)

## Check residuals
plot(K_ratio_M1) ## One exceptional value
which.max(residuals(K_ratio_M1))
Invert_breakdown_2013[34,c(1:6,10,13,16)]
## Could be due to mismatch between microbe and invertebrate estimates
## Removal does not alter result, but improves the fit
## Year Site_code Location Sample_location  Leaf Rep  k_microbe  k_invert  k_ratio
## 2013       COL      US2            1_U2 alder   4 0.03125254 0.2951408 9.443738
## Need to also exclude microbe outliers for consistency
Invert_breakdown_2013[c(91,105,162),c(1:6)]
## 91  2013       HOR      US1            2_U1 alder   4
## 105 2013       KER       DS             3_D alder   6
## 162 2013       ROM      US2            1_U2 alder   6
## Write outliers to file
Alder_kratio_2013_outliers <- Invert_breakdown_2013[-c(34,91,105,162),]
write.csv(Alder_kratio_2013_outliers,"Alder_kratio_2013_outliers.csv",row.names = T)

## Check output
summary(K_ratio_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(K_ratio_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(K_ratio_M1), lower.tail = F)
K_ratio_M1_coefficients <- cbind(ctable, pvals)
write.csv(K_ratio_M1_coefficients,"K_ratio_alder_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
K_ratio_M1_conf_intervals <- confint(K_ratio_M1, level = 0.95)
write.csv(K_ratio_M1_conf_intervals,"K_ratio_alder_2013_LME_confint.csv",row.names = T)

## DF 2, 178.09
## Write ANOVA table
K_ratio_M1_ANOVA <- anova(K_ratio_M1)
write.csv(K_ratio_M1_ANOVA ,"K_ratio_alder_2013_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(K_ratio_M1)
n <- nrow(X)
Beta <- fixef(K_ratio_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(K_ratio_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0389
## Random 0.5729
## All 0.6118
write.csv(Model_variance,"K_ratio_alder_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(K_ratio_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_ratio_alder_2013_post_Location.csv",row.names = F)

## Calculate mean values
K_ratio_data_mu_site <- data.frame(Invert_breakdown_2013[-c(34,91,105,162),] %>%
                                      group_by(Site_code,Location) %>% 
                                      summarise_at(vars("k_ratio"), c(mean), na.rm=T))
K_ratio_data_mu_location <- data.frame(K_ratio_data_mu_site %>%
                                         group_by(Location) %>% 
                                         summarise_at(vars("k_ratio"), c(mean,sd), na.rm=T))
colnames(K_ratio_data_mu_location) <- c("Location","Mean","Stdev")
write.csv(K_ratio_data_mu_location,"K_ratio_alder_2013_mean_std_location.csv",row.names = F)

##*************************************************************************
Ratio_data_2013 <- Invert_breakdown_2013[-c(34,91,105,162),]
## K_ratio - effect sizes
K_ratio_D_U1 <-  batch_calc_ES(dat = Ratio_data_2013 [-which(Ratio_data_2013 $Location=="US2"),],
                               #grouping = Site,
                               condition = Location,
                               outcome = k_ratio, 
                               improvement = "decrease",
                               ES="SMD")

K_ratio_D_U2 <-  batch_calc_ES(dat = Ratio_data_2013  [-which(Ratio_data_2013 $Location=="US1"),],
                               #grouping = Site,
                               condition = Location,
                               outcome = k_ratio, 
                               improvement = "decrease",
                               ES="SMD")

K_ratio_U1_U2 <-  batch_calc_ES(dat = Ratio_data_2013  [-which(Ratio_data_2013 $Location=="DS"),],
                                #grouping = Site,
                                condition = Location,
                                outcome = k_ratio, 
                                improvement = "decrease",
                                ES="SMD")
## Bind output
K_ratio_out <- rbind(K_ratio_D_U1,K_ratio_D_U2,K_ratio_U1_U2)
K_ratio_out$Contrast <- c("D_U1","D_U2","U1_U2")
K_ratio_out <- K_ratio_out[,c(6,1:5)]
write.csv(K_ratio_out,"K_ratio_alder_2013_effect_sizes.csv",row.names = F)

## K_ratio mean values
K_ratio_D_U1 <-  batch_calc_ES(dat = K_ratio_data_mu_site  [-which(K_ratio_data_mu_site $Location=="US2"),],
                               #grouping = Site,
                               condition = Location,
                               outcome = k_ratio, 
                               improvement = "decrease",
                               ES="SMD")

K_ratio_D_U2 <-  batch_calc_ES(dat = K_ratio_data_mu_site  [-which(K_ratio_data_mu_site $Location=="US1"),],
                               #grouping = Site,
                               condition = Location,
                               outcome = k_ratio, 
                               improvement = "decrease",
                               ES="SMD")

K_ratio_U1_U2 <-  batch_calc_ES(dat = K_ratio_data_mu_site  [-which(K_ratio_data_mu_site $Location=="DS"),],
                                #grouping = Site,
                                condition = Location,
                                outcome = k_ratio, 
                                improvement = "decrease",
                                ES="SMD")
## Bind output
K_ratio_out <- rbind(K_ratio_D_U1,K_ratio_D_U2,K_ratio_U1_U2)
K_ratio_out$Contrast <- c("D_U1","D_U2","U1_U2")
K_ratio_out <- K_ratio_out[,c(6,1:5)]
write.csv(K_ratio_out,"K_ratio_alder_2013_site_means_effect_sizes.csv",row.names = F)

#############################################################################################
##                                                                                         ##  
##                                    2013 ALDER                                           ##   
##                                                                                         ##  
#############################################################################################

## Select only 2013 sites
Alder_coarse_2013 <- Alder_allyrs[-which(Alder_allyrs$Year=="2014"),]

## Test 2013 litter breakdown k day
k_day_ml_M1 <- blmer(sqrt(k_day_ml) ~ Location + (1|Site_code), data = Alder_coarse_2013, 
                     REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_day_ml_M1)

## Check output
summary(k_day_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_day_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_day_ml_M1), lower.tail = F)
k_day_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_day_ml_M1_coefficients,"K_day_alder_Coarse_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_day_ml_M1_conf_intervals <- confint(k_day_ml_M1, level = 0.95)
write.csv(k_day_ml_M1_conf_intervals,"K_day_alder_Coarse_2013_LME_confint.csv",row.names = T)

## DF 2, 188.12
## Write ANOVA table
k_day_ml_M1_anova <- anova(k_day_ml_M1) 
k_day_ml_M1_anova$p_value <- NA
k_day_ml_M1_anova$p_value[1] <- pf(k_day_ml_M1_anova[1,4], k_day_ml_M1_anova[1,1], df.residual(k_day_ml_M1), lower.tail = FALSE)
write.csv(k_day_ml_M1_anova,"K_day_alder_Coarse_2013_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_day_ml_M1)
n <- nrow(X)
Beta <- fixef(k_day_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_day_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0207
## Random 0.6295
## All 0.6502
write.csv(Model_variance,"K_day_alder_Coarse_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_alder_Coarse_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_alder ALL
K_day_alder_mu_Coarse_2013_site <- data.frame(Alder_coarse_2013 %>%
                                        group_by(Site_code,Location) %>% 
                                        summarise_at(vars("k_day_ml"), c(mean), na.rm=T))
K_day_alder_mu_Coarse_2013_location <- data.frame(K_day_alder_mu_Coarse_2013_site  %>%
                                            group_by(Location) %>% 
                                            summarise_at(vars("k_day_ml"), c(mean,sd), na.rm=T))
colnames(K_day_alder_mu_Coarse_2013_location) <- c("Location","Mean","Stdev")
write.csv(K_day_alder_mu_Coarse_2013_location,"K_day_alder_Coarse_2013_mean_std_location.csv",row.names = F)

##*****************************************************************************
## Calculate model for coarse K degree-day
k_dd_ml_M1 <- blmer(sqrt(k_dd_ml) ~ Location + (1|Site_code), data = Alder_coarse_2013, 
                    REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_dd_ml_M1)

## Check output
summary(k_dd_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_dd_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_dd_ml_M1), lower.tail = F)
k_dd_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_dd_ml_M1_coefficients,"K_degreeday_alder_Coarse_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_dd_ml_M1_conf_intervals <- confint(k_dd_ml_M1, level = 0.95)
write.csv(k_dd_ml_M1_conf_intervals,"K_degreeday_alder_Coarse_2013_LME_confint.csv",row.names = T)

## DF 2, 188.19 
## Write ANOVA table
k_dd_ml_M1_anova <- anova(k_dd_ml_M1) 
k_dd_ml_M1_anova$p_value <- NA
k_dd_ml_M1_anova$p_value[1] <- pf(k_dd_ml_M1_anova[1,4], k_dd_ml_M1_anova[1,1], df.residual(k_dd_ml_M1), lower.tail = FALSE)
write.csv(k_dd_ml_M1_anova,"K_degreeday_alder_Coarse_2013_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_dd_ml_M1)
n <- nrow(X)
Beta <- fixef(k_dd_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_dd_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0425
## Random 0.4936
## All 0.5362
write.csv(Model_variance,"K_degreeday_alder_Coarse_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_alder_Coarse_2013_LME_post_Location.csv",row.names = F)

## Calculate means for k degree-day Alder 2013
K_dday_alder_mu_Coarse_2013_site <- data.frame(Alder_coarse_2013 %>%
                                         group_by(Site_code,Location) %>% 
                                         summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))
K_dday_alder_mu_Coarse_2013_location <- data.frame(K_dday_alder_mu_Coarse_2013_site  %>%
                                             group_by(Location) %>% 
                                             summarise_at(vars("k_dd_ml"), c(mean,sd), na.rm=T))
colnames(K_dday_alder_mu_Coarse_2013_location) <- c("Location","Mean","Stdev")
write.csv(K_dday_alder_mu_Coarse_2013_location,"K_degreeday_alder_Coarse_2013_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate effect sizes for litter breakdown
Alder_2013 <- Alder_coarse_2013
## Coarse k day
Coarse_D_U1_kday <-  batch_calc_ES(dat = Alder_2013 [-which(Alder_2013 $Location=="US2"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_D_U2_kday <-  batch_calc_ES(dat = Alder_2013 [-which(Alder_2013 $Location=="US1"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_U1_U2_kday <-  batch_calc_ES(dat = Alder_2013 [-which(Alder_2013 $Location=="DS"),],
                                    #grouping = c(Site_code,Year), 
                                    condition = Location,
                                    outcome = k_day_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

## Coarse k degree-day
Coarse_D_U1_kdd <-  batch_calc_ES(dat = Alder_2013 [-which(Alder_2013 $Location=="US2"),],
                                  #grouping = c(Site_code,Year), 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kdd <-  batch_calc_ES(dat = Alder_2013 [-which(Alder_2013 $Location=="US1"),],
                                  #grouping = c(Site_code,Year), 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_U1_U2_kdd <-  batch_calc_ES(dat = Alder_2013 [-which(Alder_2013 $Location=="DS"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_dd_ml, 
                                   improvement = "decrease",
                                   ES="SMD")
## Bind output
Coarse_kday <- rbind(Coarse_D_U1_kday,Coarse_D_U2_kday,Coarse_U1_U2_kday)
Coarse_kdd  <- rbind(Coarse_D_U1_kdd,Coarse_D_U2_kdd,Coarse_U1_U2_kdd)
Coarse_out  <- rbind(Coarse_kday,Coarse_kdd)
Coarse_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)
Coarse_2013 <- Coarse_out

## Calculate means
head(Alder_allyrs)

## Calculate site means for alternative effect size
K_day_alder_mu_ALL_site <- data.frame(Alder_coarse_2013 %>%
                                        group_by(Site_code,Location) %>% 
                                        summarise_at(vars("k_day_ml"), c(mean), na.rm=T))

K_dday_alder_mu_ALL_site <- data.frame(Alder_coarse_2013 %>%
                                         group_by(Site_code,Location) %>% 
                                         summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))

## Alder k day 2013
Coarse_D_U1_kday <-  batch_calc_ES(dat = K_day_alder_mu_ALL_site[-which(K_day_alder_mu_ALL_site$Location=="US2"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_day_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_ALL_site[-which(K_day_alder_mu_ALL_site$Location=="US1"),],
                                   #grouping = Site, 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_U1_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_ALL_site[-which(K_day_alder_mu_ALL_site$Location=="DS"),],
                                    #grouping = Site, 
                                    condition = Location,
                                    outcome = k_day_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

## Alder k degree-day 2013
Coarse_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_alder_mu_ALL_site [-which(K_dday_alder_mu_ALL_site $Location=="US2"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_ALL_site [-which(K_dday_alder_mu_ALL_site $Location=="US1"),],
                                 #grouping = Site, 
                                 condition = Location,
                                 outcome = k_dd_ml, 
                                 improvement = "decrease",
                                 ES="SMD")

Coarse_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_ALL_site [-which(K_dday_alder_mu_ALL_site $Location=="DS"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_kday <- rbind(Coarse_D_U1_kday,Coarse_D_U2_kday,Coarse_U1_U2_kday)
Coarse_kdd  <- rbind(Coarse_D_U1_kdd,Coarse_D_U2_kdd,Coarse_U1_U2_kdd)
Coarse_out  <- rbind(Coarse_kday,Coarse_kdd)
Coarse_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)
Coarse_out_mu  <- Coarse_out

#############################################################################################
##                                                                                         ##  
##                               ALL YEARS ALDER COARSE                                    ##   
##                                                                                         ##  
#############################################################################################

## Check input data
head(Alder_allyrs)

## Test litter breakdown k_day - fit with ML to improve estimation of random effects
k_day_ml_M1 <- blmer(log(k_day_ml) ~ Location + (1|Year/Site_code), data = Alder_allyrs, 
                     REML=F, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_day_ml_M1) ## One possible outlier
#which.min(residuals(k_day_ml_M1))
#Alder_allyrs[c(186),c(1:8,14)]
##     Year Site_code Location Sample_location   Mesh  Leaf Person Rep       k_ml
## 186 2014       MAR      US2            1_U2 coarse alder      M   6 0.09282493
## Probable sediment influence; no influence on result; leave in

## Check output
summary(k_day_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_day_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_day_ml_M1), lower.tail = F)
k_day_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_day_ml_M1_coefficients,"K_day_alder_Allyrs_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_day_ml_M1_conf_intervals <- confint(k_day_ml_M1, level = 0.95)
write.csv(k_day_ml_M1_conf_intervals,"K_day_alder_Allyrs_LME_confint.csv",row.names = T)

## DF: 2, 322.19
## Write ANOVA table
k_day_ml_M1_anova <- anova(k_day_ml_M1) 
k_day_ml_M1_anova$p_value <- NA
k_day_ml_M1_anova$p_value[1] <- pf(k_day_ml_M1_anova[1,4], k_day_ml_M1_anova[1,1], df.residual(k_day_ml_M1), lower.tail = FALSE)
write.csv(k_day_ml_M1_anova,"K_day_alder_Alyrs_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_day_ml_M1)
n <- nrow(X)
Beta <- fixef(k_day_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_day_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.01365
## Random 0.6280 
## All 0.6417
write.csv(Model_variance,"K_day_alder_Allyrs_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_alder_Allyrs_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_alder ALL
K_day_alder_mu_ALL_site <- data.frame(Alder_allyrs %>%
                                        group_by(Site_code,Location) %>% 
                                        summarise_at(vars("k_day_ml"), c(mean), na.rm=T))
K_day_alder_mu_ALL_location <- data.frame(K_day_alder_mu_ALL_site  %>%
                                            group_by(Location) %>% 
                                            summarise_at(vars("k_day_ml"), c(mean,sd), na.rm=T))
colnames(K_day_alder_mu_ALL_location) <- c("Location","Mean","Stdev")
write.csv(K_day_alder_mu_ALL_location,"K_day_alder_Allyrs_mean_std_location.csv",row.names = F)

##*****************************************************************************
## Calculate model for coarse K degree-day
k_dd_ml_M1 <- blmer(log(k_dd_ml) ~ Location + (1|Year/Site_code), data = Alder_allyrs, 
                    REML=F, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_dd_ml_M1)

## Check output
summary(k_dd_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_dd_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_dd_ml_M1), lower.tail = F)
k_dd_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_dd_ml_M1_coefficients,"K_degreeday_alder_Allyrs_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_dd_ml_M1_conf_intervals <- confint(k_dd_ml_M1, level = 0.95)
write.csv(k_dd_ml_M1_conf_intervals,"K_degreeday_alder_Allyrs_LME_confint.csv",row.names = T)

## DF: 2, 322.26
## Write ANOVA table
k_dd_ml_M1_anova <- anova(k_dd_ml_M1) 
k_dd_ml_M1_anova$p_value <- NA
k_dd_ml_M1_anova$p_value[1] <- pf(k_dd_ml_M1_anova[1,4], k_dd_ml_M1_anova[1,1], df.residual(k_dd_ml_M1), lower.tail = FALSE)
write.csv(k_dd_ml_M1_anova,"K_degreeday_alder_Allyrs_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(k_dd_ml_M1)
n <- nrow(X)
Beta <- fixef(k_dd_ml_M1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(k_dd_ml_M1)
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.0314
## Random 0.6092
## All 0.6406
write.csv(Model_variance,"K_degreeday_alder_Allyrs_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_alder_Allyrs_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_alder ALL
K_dday_alder_mu_ALL_site <- data.frame(Alder_allyrs %>%
                                         group_by(Site_code,Location) %>% 
                                         summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))
K_dday_alder_mu_ALL_location <- data.frame(K_dday_alder_mu_ALL_site  %>%
                                             group_by(Location) %>% 
                                             summarise_at(vars("k_dd_ml"), c(mean,sd), na.rm=T))
colnames(K_dday_alder_mu_ALL_location) <- c("Location","Mean","Stdev")
write.csv(K_dday_alder_mu_ALL_location,"K_degreeday_alder_Allyrs_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate effect sizes for litter breakdown
Alder_all <- Alder_allyrs
## Coarse k day
Coarse_D_U1_kday <-  batch_calc_ES(dat = Alder_all[-which(Alder_all$Location=="US2"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_D_U2_kday <-  batch_calc_ES(dat = Alder_all[-which(Alder_all$Location=="US1"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_U1_U2_kday <-  batch_calc_ES(dat = Alder_all[-which(Alder_all$Location=="DS"),],
                                    #grouping = c(Site_code,Year), 
                                    condition = Location,
                                    outcome = k_day_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

## Coarse k degree-day
Coarse_D_U1_kdd <-  batch_calc_ES(dat = Alder_all[-which(Alder_all$Location=="US2"),],
                                  #grouping = c(Site_code,Year), 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kdd <-  batch_calc_ES(dat = Alder_all[-which(Alder_all$Location=="US1"),],
                                  #grouping = c(Site_code,Year), 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_U1_U2_kdd <-  batch_calc_ES(dat = Alder_all[-which(Alder_all$Location=="DS"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_dd_ml, 
                                   improvement = "decrease",
                                   ES="SMD")
## Bind output
Coarse_kday <- rbind(Coarse_D_U1_kday,Coarse_D_U2_kday,Coarse_U1_U2_kday)
Coarse_kdd  <- rbind(Coarse_D_U1_kdd,Coarse_D_U2_kdd,Coarse_U1_U2_kdd)
Coarse_out  <- rbind(Coarse_kday,Coarse_kdd)
Coarse_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)
Coarse_All <- Coarse_out

## Calculate site means for alternative effect size
K_day_alder_mu_ALL_site <- data.frame(Alder_allyrs %>%
                                        group_by(Site_code,Location) %>% 
                                        summarise_at(vars("k_day_ml"), c(mean), na.rm=T))

K_dday_alder_mu_ALL_site <- data.frame(Alder_allyrs %>%
                                         group_by(Site_code,Location) %>% 
                                         summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))

## Alder k day 2013
Alder_D_U1_kday <-  batch_calc_ES(dat = K_day_alder_mu_ALL_site[-which(K_day_alder_mu_ALL_site$Location=="US2"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_day_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Alder_D_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_ALL_site[-which(K_day_alder_mu_ALL_site$Location=="US1"),],
                                   #grouping = Site, 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Alder_U1_U2_kday  <-  batch_calc_ES(dat = K_day_alder_mu_ALL_site[-which(K_day_alder_mu_ALL_site$Location=="DS"),],
                                    #grouping = Site, 
                                    condition = Location,
                                    outcome = k_day_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

## Alder k degree-day 2013
Alder_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_alder_mu_ALL_site [-which(K_dday_alder_mu_ALL_site $Location=="US2"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Alder_D_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_ALL_site [-which(K_dday_alder_mu_ALL_site $Location=="US1"),],
                                 #grouping = Site, 
                                 condition = Location,
                                 outcome = k_dd_ml, 
                                 improvement = "decrease",
                                 ES="SMD")

Alder_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_alder_mu_ALL_site [-which(K_dday_alder_mu_ALL_site $Location=="DS"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")
## Bind output
Alder_kday <- rbind(Alder_D_U1_kday,Alder_D_U2_kday,Alder_U1_U2_kday)
Alder_kdd  <- rbind(Alder_D_U1_kdd,Alder_D_U2_kdd,Alder_U1_U2_kdd)
Alder_out  <- rbind(Alder_kday,Alder_kdd)
Alder_out$Contrast <- rep(c("D_U1","D_U2","U1_U2"),2)
Alder_out_mu  <- Alder_out

#############################################################################################
##                                                                                         ##  
##        BRING TOGETHER INVERT, MICROBES, TOTAL (2013 and ALL YRS) AND EXPORT             ##   
##                                                                                         ##  
#############################################################################################

Effect_sizes <- rbind(Microbe_out,Invert_out,Coarse_2013,Coarse_All)
Effect_sizes$Response <- c(rep("K_day",3),rep("K_DD",3),rep("K_day",3),rep("K_DD",3),
                           rep("K_day",3),rep("K_DD",3),rep("K_day",3),rep("K_DD",3))
Effect_sizes$LPA <- c(rep("Microbe_2013",6),rep("Invert_2013",6),rep("Coarse_2013",6),rep("Coarse_All",6))
Effect_sizes <- Effect_sizes[,c(8,7,6,1:5)]
write.csv(Effect_sizes,"Alder_LPA_effect_sizes.csv",row.names = F)

#############################################################################################
##                                                                                         ##  
##        BRING TOGETHER INVERT, MICROBES, TOTAL (2013 and ALL YRS) AND EXPORT             ##   
##                                  MEAN VALUES                                            ##   
##                                                                                         ##  
#############################################################################################

Effect_sizes <- rbind(microbe_out_mu,Invert_out_mu,Coarse_out_mu,Alder_out_mu)
Effect_sizes$Response <- c(rep("K_day",3),rep("K_DD",3),rep("K_day",3),rep("K_DD",3),
                           rep("K_day",3),rep("K_DD",3),rep("K_day",3),rep("K_DD",3))
Effect_sizes$LPA <- c(rep("Microbe_2013",6),rep("Invert_2013",6),rep("Coarse_2013",6),rep("Coarse_All",6))
Effect_sizes <- Effect_sizes[,c(8,7,6,1:5)]
write.csv(Effect_sizes,"Alder_LPA_site_means_effect_sizes.csv",row.names = F)

#############################################################################################
##                                                                                         ##  
##                            Temperature (all years)                                      ##   
##                                                                                         ##  
#############################################################################################

## Calculate mean temperature per location
T_mu_ALL_site <- data.frame(Alder_allyrs %>%
                              group_by(Site_code,Location,Year) %>% 
                              summarise_at(vars("T_mu"), c(mean), na.rm=T))

## Test differences in temperature for Location and Year
T_mu_ALL_M1 <- blmer(log(T_mu) ~ Location + Year + (1|Site_code), data = T_mu_ALL_site, 
                     REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(T_mu_ALL_M1) ## One exceptional value
which.max(residuals(T_mu_ALL_M1))
T_mu_ALL_site[c(43:45),c(1:4)] ## Large difference at Reinach (i.e., nearly 5 degrees C warmer)

## Check output
summary(T_mu_ALL_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(T_mu_ALL_M1 ))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(T_mu_ALL_M1 ), lower.tail = F)
T_mu_ALL_M1_coefficients <- cbind(ctable, pvals)
write.csv(T_mu_ALL_M1_coefficients,"T_mu_Allyrs_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
T_mu_ALL_M1_conf_intervals <- confint(T_mu_ALL_M1 , level = 0.95)
write.csv(T_mu_ALL_M1_conf_intervals,"T_mu_Allyrs_LME_coefficients.csv",row.names = T)

## DF: 2, 38.002
## Write ANOVA table
T_mu_ALL_M1_anova <- anova(T_mu_ALL_M1) 
T_mu_ALL_M1_anova$p_value <- NA
T_mu_ALL_M1_anova$p_value[1] <- pf(T_mu_ALL_M1_anova[1,4], T_mu_ALL_M1_anova[1,1], df.residual(T_mu_ALL_M1 ), lower.tail = FALSE)
T_mu_ALL_M1_anova$p_value[2] <- pf(T_mu_ALL_M1_anova[2,4], T_mu_ALL_M1_anova[2,1], df.residual(T_mu_ALL_M1 ), lower.tail = FALSE)
write.csv(T_mu_ALL_M1_anova,"T_mu_Allyrs_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(T_mu_ALL_M1 )
n <- nrow(X)
Beta <- fixef(T_mu_ALL_M1 )
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(T_mu_ALL_M1 )
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.5992
## Random 0.3647
## All 0.9639
write.csv(Model_variance,"T_mu_Allyrs_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(T_mu_ALL_M1 , ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"T_mu_Allyrs_LME_post_Location.csv",row.names = F)

## Calculate means for Temperature All years
T_mu_ALL_site <- data.frame(Alder_allyrs %>%
                              group_by(Site_code,Location,Year) %>% 
                              summarise_at(vars("T_mu"), c(mean), na.rm=T))
T_mu_ALL_location <- data.frame(T_mu_ALL_site  %>%
                                  group_by(Location) %>% 
                                  summarise_at(vars("T_mu"), c(mean,sd), na.rm=T))
colnames(T_mu_ALL_location) <- c("Location","Mean","Stdev")
write.csv(T_mu_ALL_location,"T_mu_Allyrs_mean_std_location.csv",row.names = F)

## T_mu effect size
T_mu_ALL_D_U1 <-  batch_calc_ES(dat = T_mu_ALL_site[-which(T_mu_ALL_site$Location=="US2"),],
                               #grouping = Site,
                               condition = Location,
                               outcome = T_mu, 
                               improvement = "decrease",
                               ES="SMD")

T_mu_ALL_D_U2 <-  batch_calc_ES(dat = T_mu_ALL_site[-which(T_mu_ALL_site$Location=="US1"),],
                               #grouping = Site,
                               condition = Location,
                               outcome = T_mu, 
                               improvement = "decrease",
                               ES="SMD")

T_mu_ALL_U1_U2 <-  batch_calc_ES(dat = T_mu_ALL_site[-which(T_mu_ALL_site$Location=="DS"),],
                                #grouping = Site,
                                condition = Location,
                                outcome = T_mu, 
                                improvement = "decrease",
                                ES="SMD")

T_mu_ALL_out <- rbind(T_mu_ALL_D_U1,T_mu_ALL_D_U2,T_mu_ALL_U1_U2)
T_mu_ALL_out$Contrast <- c("D_U1","D_U2","U1_U2")
T_mu_ALL_out <- T_mu_ALL_out[,c(6,1:5)]
write.csv(T_mu_ALL_out,"T_mu_Allyrs_effect_sizes.csv",row.names = F)

#############################################################################################
##                                                                                         ##  
##                            Temperature (2013)                                           ##   
##                                                                                         ##  
#############################################################################################

## Use only data from 2013
T_mu_2013_site <- T_mu_ALL_site[which(T_mu_ALL_site$Year=="2013"),]

## Test differences in temperature for Location and Year
T_mu_2013_M1 <- blmer(log(T_mu) ~ Location + (1|Site_code), data = T_mu_2013_site, 
                     REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(T_mu_2013_M1)

## Check output
summary(T_mu_2013_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(T_mu_2013_M1 ))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(T_mu_2013_M1 ), lower.tail = F)
T_mu_2013_M1_coefficients <- cbind(ctable, pvals)
write.csv(T_mu_2013_M1_coefficients,"T_mu_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
T_mu_2013_M1_conf_intervals <- confint(T_mu_2013_M1 , level = 0.95)
write.csv(T_mu_2013_M1_conf_intervals,"T_mu_2013_LME_coefficients.csv",row.names = T)

## DF: 2, 22
## Write ANOVA table
T_mu_2013_M1_anova <- anova(T_mu_2013_M1) 
T_mu_2013_M1_anova$p_value <- NA
T_mu_2013_M1_anova$p_value[1] <- pf(T_mu_2013_M1_anova[1,4], T_mu_2013_M1_anova[1,1], df.residual(T_mu_2013_M1 ), lower.tail = FALSE)
write.csv(T_mu_2013_M1_anova,"T_mu_2013_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(T_mu_2013_M1 )
n <- nrow(X)
Beta <- fixef(T_mu_2013_M1 )
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(T_mu_2013_M1 )
Sl <- 
  sum(
    sapply(Sigma.list,
           function(Sigma)
           {
             Z <-X[,rownames(Sigma)]
             sum(diag(Z %*% Sigma %*% t(Z)))/n
           }))
Se <- attr(Sigma.list, "sc")^2
Sd <- 0
total.var <- Sf + Sl + Se + Sd

Marginal <- (Rsq.m <- Sf / total.var) 
Conditional <- (Rsq.c <- (Sf + Sl) / total.var) 
Random <- (Rsq.r<- Sl / total.var)

Model_variance <- rbind(cbind("Marginal",Marginal),
                        cbind("Conditional",Conditional),
                        cbind("Random",Random ))
## Fixed 0.02707
## Random 0.9453
## 2013 0.9724
write.csv(Model_variance,"T_mu_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(T_mu_2013_M1 , ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"T_mu_2013_LME_post_Location.csv",row.names = F)

## Calculate means for Temperature 2013
T_mu_2013_location <- data.frame(T_mu_2013_site  %>%
                                  group_by(Location) %>% 
                                  summarise_at(vars("T_mu"), c(mean,sd), na.rm=T))
colnames(T_mu_2013_location) <- c("Location","Mean","Stdev")
write.csv(T_mu_2013_location,"T_mu_2013_mean_std_location.csv",row.names = F)

## T_mu effect size
T_mu_2013_D_U1 <-  batch_calc_ES(dat = T_mu_2013_site[-which(T_mu_2013_site$Location=="US2"),],
                                #grouping = Site,
                                condition = Location,
                                outcome = T_mu, 
                                improvement = "decrease",
                                ES="SMD")

T_mu_2013_D_U2 <-  batch_calc_ES(dat = T_mu_2013_site[-which(T_mu_2013_site$Location=="US1"),],
                                #grouping = Site,
                                condition = Location,
                                outcome = T_mu, 
                                improvement = "decrease",
                                ES="SMD")

T_mu_2013_U1_U2 <-  batch_calc_ES(dat = T_mu_2013_site[-which(T_mu_2013_site$Location=="DS"),],
                                 #grouping = Site,
                                 condition = Location,
                                 outcome = T_mu, 
                                 improvement = "decrease",
                                 ES="SMD")

T_mu_2013_out <- rbind(T_mu_2013_D_U1,T_mu_2013_D_U2,T_mu_2013_U1_U2)
T_mu_2013_out$Contrast <- c("D_U1","D_U2","U1_U2")
T_mu_2013_out <- T_mu_2013_out[,c(6,1:5)]
write.csv(T_mu_2013_out,"T_mu_2013_effect_sizes.csv",row.names = F)
