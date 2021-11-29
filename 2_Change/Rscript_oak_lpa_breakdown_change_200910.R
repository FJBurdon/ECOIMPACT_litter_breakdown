#############################################################################################
##                                                                                         ##  
##                    IMPACTS OF WW: MIXED MODELS AND EFFECT SIZES                         ##  
##                                      OAK                                                ##   
##                                                                                         ##  
#############################################################################################

##  Script: ECOIMPACT leafpack assay data analysis
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Calculate mean +/- SD for oak breakdown
##  2. Test changes in oak breakdown using linear mixed-models
##  3. Calculate effect sizes for oak breakdown

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

## 2013 coarse mesh decomposition
Oak_2013 <- read.csv(here::here("2_Oak_data_input","lpa_oak_coarse_2013_data.csv")) 

## 2013 Microbe and Invertebrate decomposition
Microbial_breakdown_2013 <- read.csv(here::here("2_Oak_data_input","lpa_oak_fine_2013_data.csv")) 
Invert_breakdown_2013 <- read.csv(here::here("2_Oak_data_input","lpa_oak_invertebrate_2013_data.csv")) 

#############################################################################################
##                                                                                         ##  
##                      INVERTEBRATE EFFECT SIZES - 2013  OAK                              ##   
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
write.csv(k_invert_day_ml_M1_coefficients,"K_day_oak_invert_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_invert_day_ml_M1_conf_intervals <- confint(k_invert_day_ml_M1, level = 0.95)
write.csv(k_invert_day_ml_M1_conf_intervals,"K_day_oak_invert_2013_LME_confint.csv",row.names = T)

## DF 2, 183.13
## Write ANOVA table
k_invert_day_ml_M1_anova <- anova(k_invert_day_ml_M1) 
k_invert_day_ml_M1_anova$p_value <- NA
k_invert_day_ml_M1_anova$p_value[1] <- pf(k_invert_day_ml_M1_anova[1,4], k_invert_day_ml_M1_anova[1,1], df.residual(k_invert_day_ml_M1), lower.tail = FALSE)
write.csv(k_invert_day_ml_M1_anova,"K_day_oak_invert_2013_LME_ANOVA.csv",row.names = T)

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
## Fixed 0.0240
## Random 0.4090
## All 0.4329
write.csv(Model_variance,"K_day_oak_invert_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_invert_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_oak_invert_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_oak 2013
K_day_oak_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_day"), c(mean), na.rm=T))
K_day_oak_mu_invert_location <- data.frame(K_day_oak_mu_invert_site  %>%
                                               group_by(Location) %>% 
                                               summarise_at(vars("k_invert_day"), c(mean,sd), na.rm=T))
colnames(K_day_oak_mu_invert_location) <- c("Location","Mean","Stdev")
write.csv(K_day_oak_mu_invert_location,"K_day_oak_invert_2013_mean_std_location.csv",row.names = F)

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
write.csv(k_invert_dd_ml_M1_coefficients,"K_degreeday_oak_invert_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_invert_dd_ml_M1_conf_intervals <- confint(k_invert_dd_ml_M1, level = 0.95)
write.csv(k_invert_dd_ml_M1_conf_intervals,"K_degreeday_oak_invert_2013_LME_confint.csv",row.names = T)

## DF 2, 183.17
## Write ANOVA table
k_invert_dd_ml_M1_anova <- anova(k_invert_dd_ml_M1) 
k_invert_dd_ml_M1_anova$p_value <- NA
k_invert_dd_ml_M1_anova$p_value[1] <- pf(k_invert_dd_ml_M1_anova[1,4], k_invert_dd_ml_M1_anova[1,1], df.residual(k_invert_dd_ml_M1), lower.tail = FALSE)
write.csv(k_invert_dd_ml_M1_anova,"K_degreeday_oak_invert_2013_LME_ANOVA.csv",row.names = T)

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
## Fixed 0.0374
## Random 0.3306
## All 0.3680
write.csv(Model_variance,"K_degreeday_oak_invert_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_invert_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_oak_invert_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_oak ALL
K_day_oak_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_dd"), c(mean), na.rm=T))
K_day_oak_mu_invert_location <- data.frame(K_day_oak_mu_invert_site  %>%
                                               group_by(Location) %>% 
                                               summarise_at(vars("k_invert_dd"), c(mean,sd), na.rm=T))
colnames(K_day_oak_mu_invert_location) <- c("Location","Mean","Stdev")
write.csv(K_day_oak_mu_invert_location,"K_degreeday_oak_invert_2013_mean_std_location.csv",row.names = F)

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
K_day_oak_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_day"), c(mean), na.rm=T))

K_dday_oak_mu_invert_site <- data.frame(Invert_breakdown_2013 %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_invert_dd"), c(mean), na.rm=T))

## Invertebrate k day 2013
Invertebrate_D_U1_kday <-  batch_calc_ES(dat = K_day_oak_mu_invert_site [-which(K_day_oak_mu_invert_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_day, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_D_U2_kday  <-  batch_calc_ES(dat = K_day_oak_mu_invert_site [-which(K_day_oak_mu_invert_site $Location=="US1"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_invert_day, 
                                          improvement = "decrease",
                                          ES="SMD")

Invertebrate_U1_U2_kday  <-  batch_calc_ES(dat = K_day_oak_mu_invert_site [-which(K_day_oak_mu_invert_site $Location=="DS"),],
                                           #grouping = Site, 
                                           condition = Location,
                                           outcome = k_invert_day, 
                                           improvement = "decrease",
                                           ES="SMD")

## Invertebrate k degree-day 2013
Invertebrate_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_oak_mu_invert_site [-which(K_dday_oak_mu_invert_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_invert_dd, 
                                         improvement = "decrease",
                                         ES="SMD")

Invertebrate_D_U2_kdd <-  batch_calc_ES(dat = K_dday_oak_mu_invert_site [-which(K_dday_oak_mu_invert_site $Location=="US1"),],
                                        #grouping = Site, 
                                        condition = Location,
                                        outcome = k_invert_dd, 
                                        improvement = "decrease",
                                        ES="SMD")

Invertebrate_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_oak_mu_invert_site [-which(K_dday_oak_mu_invert_site $Location=="DS"),],
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
##                        MICROBIAL EFFECT SIZES - 2013 OAK                                ##   
##                                                                                         ##  
#############################################################################################

## Check input data
head(Microbial_breakdown_2013)

## Write outliers
Oak_microbe_2013_outliers <- Microbial_breakdown_2013[c(87,121,133,198),]
write.csv(Oak_microbe_2013_outliers,"Oak_microbe_2013_outliers.csv",row.names = F)

## Test microbially-mediated breakdown k
#k_microbe_day_ml_M1 <- blmer(sqrt(k_day_ml) ~ Location + (1|Site_code), data = Microbial_breakdown_2013, 
#                             REML=T, control = lmerControl(optimizer = "nloptwrap"))
k_microbe_day_ml_M1 <- blmer(sqrt(k_day_ml) ~ Location + (1|Site_code), data = Microbial_breakdown_2013[-c(87,121,133,198),], 
                              REML=T, control = lmerControl(optimizer = "nloptwrap"))
## Check residuals
plot(k_microbe_day_ml_M1)
which.max(residuals(k_microbe_day_ml_M1))
which.min(residuals(k_microbe_day_ml_M1))
#Microbial_breakdown_2013[c(87,133),c(1:8,14)]
##      Year Site_code Location Sample_location Mesh  Leaf Person Rep      k_ml
##  87  2013       HOR       DS             3_D fine  oak      P   4  0.4263141
##  133 2013       NIE       DS             3_D fine  oak      S   1  0.1321595
## Human error - #87 look to be much larger than expected, and #133 lower 
## Doesn´t change the result, but improved the model fit
## SEE BELOW FOR k degree-day (for consistent approach to data)
#Microbial_breakdown_2013[c(121,198),c(1:7,14)]
##    Year Site_code Location Sample_location Mesh  Leaf Person          k_ml
##121 2013       MES      US1            2_U1 fine  oak      S      0.2483309
##198 2013       SEV      US2            1_U2 fine  oak      P      0.2666130
## High values, probable human error; doesnt affect result but immproves fit

## Check output
summary(k_microbe_day_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_microbe_day_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_microbe_day_ml_M1), lower.tail = F)
k_microbe_day_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_microbe_day_ml_M1_coefficients,"K_day_oak_microbe_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_microbe_day_ml_M1_conf_intervals <- confint(k_microbe_day_ml_M1, level = 0.95)
write.csv(k_microbe_day_ml_M1_conf_intervals,"K_day_oak_microbe_2013_LME_confint.csv",row.names = T)

## DF 2, 181.1 
## Write ANOVA table
k_microbe_day_ml_M1_anova <- anova(k_microbe_day_ml_M1) 
k_microbe_day_ml_M1_anova$p_value <- NA
k_microbe_day_ml_M1_anova$p_value[1] <- pf(k_microbe_day_ml_M1_anova[1,4], k_microbe_day_ml_M1_anova[1,1], df.residual(k_microbe_day_ml_M1), lower.tail = FALSE)
write.csv(k_microbe_day_ml_M1_anova,"K_day_oak_microbe_2013_LME_ANOVA.csv",row.names = T)

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
## Fixed 0.0013
## Random 0.4901
## All 0.4914
write.csv(Model_variance,"K_day_oak_microbe_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_microbe_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_oak_microbe_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_oak microbe
K_day_oak_mu_microbe_site <- data.frame(Microbial_breakdown_2013[-c(87,121,133,198),] %>%
                                            group_by(Site_code,Location) %>% 
                                            summarise_at(vars("k_day_ml"), c(mean), na.rm=T))
K_day_oak_mu_microbe_location <- data.frame(K_day_oak_mu_microbe_site  %>%
                                                group_by(Location) %>% 
                                                summarise_at(vars("k_day_ml"), c(mean,sd), na.rm=T))
colnames(K_day_oak_mu_microbe_location) <- c("Location","Mean","Stdev")
write.csv(K_day_oak_mu_microbe_location,"K_day_oak_microbe_2013_mean_std_location.csv",row.names = F)

##*****************************************************************************/
## Calculate model for microbial K degree day
k_microbe_dd_ml_M1 <- blmer(sqrt(k_dd_ml) ~ Location + (1|Site_code), data = Microbial_breakdown_2013[-c(87,121,133,198),], 
                            REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_microbe_dd_ml_M1)
#which.max(residuals(k_microbe_dd_ml_M1))
#Microbial_breakdown_2013[c(121,198),c(1:7,14)]
##    Year Site_code Location Sample_location Mesh  Leaf Person          k_ml
##121 2013       MES      US1            2_U1 fine  oak      S      0.2483309
##198 2013       SEV      US2            1_U2 fine  oak      P      0.2666130
## High values, probable human error; doesnt affect result but immproves fit

## Check output
summary(k_microbe_dd_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_microbe_dd_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_microbe_dd_ml_M1), lower.tail = F)
k_microbe_dd_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_microbe_dd_ml_M1_coefficients,"K_degreeday_oak_microbe_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_microbe_dd_ml_M1_conf_intervals <- confint(k_microbe_dd_ml_M1, level = 0.95)
write.csv(k_microbe_dd_ml_M1_conf_intervals,"K_degreeday_oak_microbe_2013_LME_confint.csv",row.names = T)

## DF 2, 181.44
## Write ANOVA table
k_microbe_dd_ml_M1_anova <- anova(k_microbe_dd_ml_M1) 
k_microbe_dd_ml_M1_anova$p_value <- NA
k_microbe_dd_ml_M1_anova$p_value[1] <- pf(k_microbe_dd_ml_M1_anova[1,4], k_microbe_dd_ml_M1_anova[1,1], df.residual(k_microbe_dd_ml_M1), lower.tail = FALSE)
write.csv(k_microbe_dd_ml_M1_anova,"K_degreeday_oak_microbe_2013_LME_ANOVA_200816.csv",row.names = T)

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
## Fixed 0.0187
## Random 0.2681
## All 0.2868
write.csv(Model_variance,"K_degreeday_oak_microbe_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_microbe_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_oak_microbe_2013_LME_post_Location.csv",row.names = F)

## Calculate means for microbially-mediated K degree-day oak
K_day_oak_mu_microbe_site <- data.frame(Microbial_breakdown_2013[-c(87,121,133,198),] %>%
                                            group_by(Site_code,Location) %>% 
                                            summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))
K_day_oak_mu_microbe_location <- data.frame(K_day_oak_mu_microbe_site  %>%
                                                group_by(Location) %>% 
                                                summarise_at(vars("k_dd_ml"), c(mean,sd), na.rm=T))
colnames(K_day_oak_mu_microbe_location) <- c("Location","Mean","Stdev")
write.csv(K_day_oak_mu_microbe_location,"K_degreeday_oak_microbe_2013_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate effect sizes for microbially-mediated breakdown
Microbial_2013 <- Microbial_breakdown_2013[-c(87,121,133,198),]
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
K_day_oak_mu_microbe_site <- data.frame(Microbial_breakdown_2013[-c(87,121,133,198),] %>%
                                           group_by(Site_code,Location) %>% 
                                           summarise_at(vars("k_day_ml"), c(mean), na.rm=T))

K_dday_oak_mu_microbe_site <- data.frame(Microbial_breakdown_2013 %>%
                                            group_by(Site_code,Location) %>% 
                                            summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))

## Microbe k day 2013
microbe_D_U1_kday <-  batch_calc_ES(dat = K_day_oak_mu_microbe_site [-which(K_day_oak_mu_microbe_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_day_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

microbe_D_U2_kday  <-  batch_calc_ES(dat = K_day_oak_mu_microbe_site [-which(K_day_oak_mu_microbe_site $Location=="US1"),],
                                          #grouping = Site, 
                                          condition = Location,
                                          outcome = k_day_ml, 
                                          improvement = "decrease",
                                          ES="SMD")

microbe_U1_U2_kday  <-  batch_calc_ES(dat = K_day_oak_mu_microbe_site [-which(K_day_oak_mu_microbe_site $Location=="DS"),],
                                           #grouping = Site, 
                                           condition = Location,
                                           outcome = k_day_ml, 
                                           improvement = "decrease",
                                           ES="SMD")

## Microbe k degree-day 2013
microbe_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_oak_mu_microbe_site [-which(K_dday_oak_mu_microbe_site $Location=="US2"),],
                                         #grouping = Site, 
                                         condition = Location,
                                         outcome = k_dd_ml, 
                                         improvement = "decrease",
                                         ES="SMD")

microbe_D_U2_kdd <-  batch_calc_ES(dat = K_dday_oak_mu_microbe_site [-which(K_dday_oak_mu_microbe_site $Location=="US1"),],
                                        #grouping = Site, 
                                        condition = Location,
                                        outcome = k_dd_ml, 
                                        improvement = "decrease",
                                        ES="SMD")

microbe_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_oak_mu_microbe_site [-which(K_dday_oak_mu_microbe_site $Location=="DS"),],
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
K_ratio_M1 <- lmer(sqrt(k_ratio) ~ Location + (1|Site_code), data = Invert_breakdown_2013[-c(86,123,133,196),], REML=T)

## Check residuals
plot(K_ratio_M1)

## Write outliers to file
oak_kratio_2013_outliers <- Invert_breakdown_2013[c(86,123,133,196),]
write.csv(oak_kratio_2013_outliers,"Oak_kratio_2013_outliers.csv",row.names = T)

## Check output
summary(K_ratio_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(K_ratio_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(K_ratio_M1), lower.tail = F)
K_ratio_M1_coefficients <- cbind(ctable, pvals)
write.csv(K_ratio_M1_coefficients,"K_ratio_oak_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
K_ratio_M1_conf_intervals <- confint(K_ratio_M1, level = 0.95)
write.csv(K_ratio_M1_conf_intervals,"K_ratio_oak_2013_LME_confint.csv",row.names = T)

## DF 2, 178.09
## Write ANOVA table
K_ratio_M1_ANOVA <- anova(K_ratio_M1)
write.csv(K_ratio_M1_ANOVA ,"K_ratio_oak_2013_ANOVA.csv",row.names = T)

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
## Fixed 0.0307
## Random 0.2786
## All 0.3093
write.csv(Model_variance,"K_ratio_oak_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(K_ratio_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_ratio_oak_2013_post_Location.csv",row.names = F)

## Calculate mean values
K_ratio_data_mu_site <- data.frame(Invert_breakdown_2013[-c(34,91,105,162),] %>%
                                      group_by(Site_code,Location) %>% 
                                      summarise_at(vars("k_ratio"), c(mean), na.rm=T))
K_ratio_data_mu_location <- data.frame(K_ratio_data_mu_site %>%
                                         group_by(Location) %>% 
                                         summarise_at(vars("k_ratio"), c(mean,sd), na.rm=T))
colnames(K_ratio_data_mu_location) <- c("Location","Mean","Stdev")
write.csv(K_ratio_data_mu_location,"K_ratio_oak_2013_mean_std_location.csv",row.names = F)

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
write.csv(K_ratio_out,"K_ratio_oak_2013_effect_sizes.csv",row.names = F)

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
write.csv(K_ratio_out,"K_ratio_oak_2013_site_means_effect_sizes.csv",row.names = F)

#############################################################################################
##                                                                                         ##  
##                                  2013 OAK COARSE                                        ##   
##                                                                                         ##  
#############################################################################################

oak_coarse_2013 <- Oak_2013

## Test 2013 litter breakdown k day
k_day_ml_M1 <- blmer(sqrt(k_day_ml) ~ Location + (1|Site_code), data = oak_coarse_2013, 
                     REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_day_ml_M1)

## Check output
summary(k_day_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_day_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_day_ml_M1), lower.tail = F)
k_day_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_day_ml_M1_coefficients,"K_day_oak_Coarse_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_day_ml_M1_conf_intervals <- confint(k_day_ml_M1, level = 0.95)
write.csv(k_day_ml_M1_conf_intervals,"K_day_oak_Coarse_2013_LME_confint.csv",row.names = T)

## DF 2, 192.05
## Write ANOVA table
k_day_ml_M1_anova <- anova(k_day_ml_M1) 
k_day_ml_M1_anova$p_value <- NA
k_day_ml_M1_anova$p_value[1] <- pf(k_day_ml_M1_anova[1,4], k_day_ml_M1_anova[1,1], df.residual(k_day_ml_M1), lower.tail = FALSE)
write.csv(k_day_ml_M1_anova,"K_day_oak_Coarse_2013_LME_ANOVA.csv",row.names = T)

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
## Fixed  0.0161
## Random 0.5407
## All    0.5569
write.csv(Model_variance,"K_day_oak_Coarse_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_day_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_day_oak_Coarse_2013_LME_post_Location.csv",row.names = F)

## Calculate means for K_day_oak ALL
K_day_oak_mu_Coarse_2013_site <- data.frame(oak_coarse_2013 %>%
                                        group_by(Site_code,Location) %>% 
                                        summarise_at(vars("k_day_ml"), c(mean), na.rm=T))
K_day_oak_mu_Coarse_2013_location <- data.frame(K_day_oak_mu_Coarse_2013_site  %>%
                                            group_by(Location) %>% 
                                            summarise_at(vars("k_day_ml"), c(mean,sd), na.rm=T))
colnames(K_day_oak_mu_Coarse_2013_location) <- c("Location","Mean","Stdev")
write.csv(K_day_oak_mu_Coarse_2013_location,"K_day_oak_Coarse_2013_mean_std_location.csv",row.names = F)

##*****************************************************************************
## Calculate model for coarse K degree-day
k_dd_ml_M1 <- blmer(sqrt(k_dd_ml) ~ Location + (1|Site_code), data = oak_coarse_2013, 
                    REML=T, control = lmerControl(optimizer = "nloptwrap"))

## Check residuals
plot(k_dd_ml_M1)

## Check output
summary(k_dd_ml_M1)

## Write Parameter estimates and P-values
ctable <- coef(summary(k_dd_ml_M1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(k_dd_ml_M1), lower.tail = F)
k_dd_ml_M1_coefficients <- cbind(ctable, pvals)
write.csv(k_dd_ml_M1_coefficients,"K_degreeday_oak_Coarse_2013_LME_coefficients.csv",row.names = T)

## Write Confidence Intervals
k_dd_ml_M1_conf_intervals <- confint(k_dd_ml_M1, level = 0.95)
write.csv(k_dd_ml_M1_conf_intervals,"K_degreeday_oak_Coarse_2013_LME_confint.csv",row.names = T)

## DF 2, 192.2 
## Write ANOVA table
k_dd_ml_M1_anova <- anova(k_dd_ml_M1) 
k_dd_ml_M1_anova$p_value <- NA
k_dd_ml_M1_anova$p_value[1] <- pf(k_dd_ml_M1_anova[1,4], k_dd_ml_M1_anova[1,1], df.residual(k_dd_ml_M1), lower.tail = FALSE)
write.csv(k_dd_ml_M1_anova,"K_degreeday_oak_Coarse_2013_LME_ANOVA.csv",row.names = T)

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
## Fixed 0.0878
## Random 0.2349
## All 0.3227
write.csv(Model_variance,"K_degreeday_oak_Coarse_2013_LME_variance.csv",row.names = T)

## Write post-hoc test output
L.S <- pairs(lsmeans(k_dd_ml_M1, ~ Location))
Post_Location <- test(L.S, adjust = "tukey")
write.csv(Post_Location ,"K_degreeday_oak_Coarse_2013_LME_post_Location.csv",row.names = F)

## Calculate means for k degree-day oak 2013
K_dday_oak_mu_Coarse_2013_site <- data.frame(oak_coarse_2013 %>%
                                         group_by(Site_code,Location) %>% 
                                         summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))
K_dday_oak_mu_Coarse_2013_location <- data.frame(K_dday_oak_mu_Coarse_2013_site  %>%
                                             group_by(Location) %>% 
                                             summarise_at(vars("k_dd_ml"), c(mean,sd), na.rm=T))
colnames(K_dday_oak_mu_Coarse_2013_location) <- c("Location","Mean","Stdev")
write.csv(K_dday_oak_mu_Coarse_2013_location,"K_degreeday_oak_Coarse_2013_mean_std_location.csv",row.names = F)

##**************************************************************************
## Calculate effect sizes for litter breakdown
## Coarse k day
Coarse_D_U1_kday <-  batch_calc_ES(dat = Oak_2013 [-which(Oak_2013 $Location=="US2"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_D_U2_kday <-  batch_calc_ES(dat = Oak_2013 [-which(Oak_2013 $Location=="US1"),],
                                   #grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_U1_U2_kday <-  batch_calc_ES(dat = Oak_2013 [-which(Oak_2013 $Location=="DS"),],
                                    #grouping = c(Site_code,Year), 
                                    condition = Location,
                                    outcome = k_day_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

## Coarse k degree-day
Coarse_D_U1_kdd <-  batch_calc_ES(dat = Oak_2013 [-which(Oak_2013 $Location=="US2"),],
                                  #grouping = c(Site_code,Year), 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kdd <-  batch_calc_ES(dat = Oak_2013 [-which(Oak_2013 $Location=="US1"),],
                                  #grouping = c(Site_code,Year), 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_U1_U2_kdd <-  batch_calc_ES(dat = Oak_2013 [-which(Oak_2013 $Location=="DS"),],
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

## Calculate site means for alternative effect size
K_day_oak_mu_ALL_site <- data.frame(oak_coarse_2013 %>%
                                        group_by(Site_code,Location) %>% 
                                        summarise_at(vars("k_day_ml"), c(mean), na.rm=T))

K_dday_oak_mu_ALL_site <- data.frame(oak_coarse_2013 %>%
                                         group_by(Site_code,Location) %>% 
                                         summarise_at(vars("k_dd_ml"), c(mean), na.rm=T))

## oak k day 2013
Coarse_D_U1_kday <-  batch_calc_ES(dat = K_day_oak_mu_ALL_site[-which(K_day_oak_mu_ALL_site$Location=="US2"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_day_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kday  <-  batch_calc_ES(dat = K_day_oak_mu_ALL_site[-which(K_day_oak_mu_ALL_site$Location=="US1"),],
                                   #grouping = Site, 
                                   condition = Location,
                                   outcome = k_day_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_U1_U2_kday  <-  batch_calc_ES(dat = K_day_oak_mu_ALL_site[-which(K_day_oak_mu_ALL_site$Location=="DS"),],
                                    #grouping = Site, 
                                    condition = Location,
                                    outcome = k_day_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

## oak k degree-day 2013
Coarse_D_U1_kdd  <-  batch_calc_ES(dat = K_dday_oak_mu_ALL_site [-which(K_dday_oak_mu_ALL_site $Location=="US2"),],
                                  #grouping = Site, 
                                  condition = Location,
                                  outcome = k_dd_ml, 
                                  improvement = "decrease",
                                  ES="SMD")

Coarse_D_U2_kdd <-  batch_calc_ES(dat = K_dday_oak_mu_ALL_site [-which(K_dday_oak_mu_ALL_site $Location=="US1"),],
                                 #grouping = Site, 
                                 condition = Location,
                                 outcome = k_dd_ml, 
                                 improvement = "decrease",
                                 ES="SMD")

Coarse_U1_U2_kdd <-  batch_calc_ES(dat = K_dday_oak_mu_ALL_site [-which(K_dday_oak_mu_ALL_site $Location=="DS"),],
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
##        BRING TOGETHER INVERT, MICROBES, TOTAL (2013 and ALL YRS) AND EXPORT             ##   
##                                                                                         ##  
#############################################################################################

Effect_sizes <- rbind(Microbe_out,Invert_out,Coarse_2013)
Effect_sizes$Response <- c(rep("K_day",3),rep("K_DD",3),rep("K_day",3),rep("K_DD",3),
                           rep("K_day",3),rep("K_DD",3))
Effect_sizes$LPA <- c(rep("Microbe_2013",6),rep("Invert_2013",6),rep("Coarse_2013",6))
Effect_sizes <- Effect_sizes[,c(8,7,6,1:5)]
write.csv(Effect_sizes,"Oak_LPA_effect_sizes.csv",row.names = F)

#############################################################################################
##                                                                                         ##  
##        BRING TOGETHER INVERT, MICROBES, TOTAL (2013 and ALL YRS) AND EXPORT             ##   
##                                  MEAN VALUES                                            ##   
##                                                                                         ##  
#############################################################################################

Effect_sizes <- rbind(microbe_out_mu,Invert_out_mu,Coarse_out_mu)
Effect_sizes$Response <- c(rep("K_day",3),rep("K_DD",3),rep("K_day",3),rep("K_DD",3),
                           rep("K_day",3),rep("K_DD",3))
Effect_sizes$LPA <- c(rep("Microbe_2013",6),rep("Invert_2013",6),rep("Coarse_2013",6))
Effect_sizes <- Effect_sizes[,c(8,7,6,1:5)]
write.csv(Effect_sizes,"Oak_LPA_site_means_effect_sizes.csv",row.names = F)
