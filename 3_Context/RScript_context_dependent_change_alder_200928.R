#############################################################################################
##                                                                                         ##  
##                          EFFECT SIZES AND CONTEXT DEPENDENCY                            ##   
##                                                                                         ##  
#############################################################################################

##  Script: ECOIMPACT leafpack assay data analysis
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Calculate effect sizes for alder breakdown
##  2. Test context-dependent changes in alder breakdown using linear regression and LME
##      - specifcally the relationship between upstream trait abundances and litter breakdown
##  3. Plot results for manuscript

##***************************************************************************************
## Load libraries you need

##Loads package used for Logit transformation (% data)
require(car)

##Loads vegan package for transformation
library(vegan)

##Loads tudyverse, reshape and plyr package for data manipulation
library(reshape)
library(reshape2)
library(plyr)
require(tidyverse)

## Load data files
require(readr)

## To calculate effect sizes
require(SingleCaseES)

## To fit statistical models
require(blme)
require(mgcv)

## For results tables
require(sjPlot)

#############################################################################################
##                                                                                         ##
##                                  LOAD DATA                                              ##    
##                                                                                         ##
#############################################################################################

## All years decomposition
Alder_coarse_allyrs <-  read.csv(here::here("1_Data_input","lpa_alder_coarse_All_data.csv")) 

## Shredders
CWM_traits <- read.csv(here::here("1_Data_input","lpa_alder_coarse_allyrs_CWM_traits_data_all_leapacks.csv")) 
Taxa_abundances <- read.csv(here::here("1_Data_input","lpa_alder_coarse_allyrs_invert_abundances_data.csv")) 

#############################################################################################
##                                                                                         ##  
##                                  EFFECT SIZES - ALL YEARS                               ##   
##                                                                                         ##  
#############################################################################################

head(Alder_coarse_allyrs)

## Create copy of input data
Alder_allyrs  <- Alder_coarse_allyrs[!is.na(Alder_coarse_allyrs$k_dd_ml), ] 

## Ensure rows are in chronological order
rownames(Alder_allyrs) <- NULL

## Coarse k dd
Coarse_D_U1_site <-  batch_calc_ES(dat = Alder_allyrs [-which(Alder_allyrs$Location=="US2"),],
                                    grouping = c(Site_code,Year), 
                                    condition = Location,
                                    outcome = k_dd_ml, 
                                    improvement = "decrease",
                                    ES="SMD")

Coarse_D_U2_site <-  batch_calc_ES(dat = Alder_allyrs [-which(Alder_allyrs$Location=="US1"),],
                                   grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_dd_ml, 
                                   improvement = "decrease",
                                   ES="SMD")

Coarse_U1_U2_site <-  batch_calc_ES(dat = Alder_allyrs [-which(Alder_allyrs$Location=="DS"),],
                                   grouping = c(Site_code,Year), 
                                   condition = Location,
                                   outcome = k_dd_ml, 
                                   improvement = "decrease",
                                   ES="SMD")
## Create output
Coarse_D_U1_out <- cbind(rep("D_U1",20),Coarse_D_U1_site)
Coarse_D_U2_out <- cbind(rep("D_U2",20),Coarse_D_U2_site)
Coarse_U1_U2_out <- cbind(rep("U1_U2",20),Coarse_U1_U2_site)
colnames(Coarse_D_U1_out) <- c("Contrast","Site","Year","ES","Est","SE","CI_lower","CI_upper")
colnames(Coarse_D_U2_out) <- c("Contrast","Site","Year","ES","Est","SE","CI_lower","CI_upper")
colnames(Coarse_U1_U2_out) <- c("Contrast","Site","Year","ES","Est","SE","CI_lower","CI_upper")
Allyrs_effect_sizes <- rbind(Coarse_D_U1_out,Coarse_D_U2_out,Coarse_U1_U2_out)
write.csv(Allyrs_effect_sizes,"Data_Allyrs_effect_sizes_alder_coarse.csv",row.names = F)

median(Coarse_D_U1_site$Est)
median(Coarse_D_U2_site$Est)
median(Coarse_U1_U2_site$Est)

#############################################################################################
##                                                                                         ##  
##                            CHANGE VS DETRITIVORES (U1 and D)                            ##
##                                                                                         ##  
#############################################################################################

head(CWM_traits[,c(1:5)]) 

## Calculate mean values for CWM trait abundances - Foo and Feeding preference
CWM_traits_mu <- data.frame(CWM_traits %>%
                    group_by(Year,Site_code,Location,Mesh,Leaf) %>% 
                    summarise_at(vars("Food_3","Feeding_3"), mean, na.rm=T))

## Spread values by location
CWM_Food_3_mu_wide <- CWM_traits_mu[,c(1:3,6)]  %>% spread(Location, Food_3)
colnames(CWM_Food_3_mu_wide) <- c("Year","Site","Food_3_D","Food_3_U1","Food_3_U2")
Coarse_D_U1_out <- merge(Coarse_D_U1_out,CWM_Food_3_mu_wide,by=c("Site","Year"))

## Spread values by location
CWM_Feeding_3_mu_wide <- CWM_traits_mu[,c(1:3,7)]  %>% spread(Location, Feeding_3)
colnames(CWM_Feeding_3_mu_wide) <- c("Year","Site","Feeding_3_D","Feeding_3_U1","Feeding_3_U2")
Coarse_D_U1_out <- merge(Coarse_D_U1_out,CWM_Feeding_3_mu_wide,by=c("Site","Year"))

## Calculate mean values for gammrid abundances - dominant shredder taxon
Taxa_abundances_mu <- data.frame(Taxa_abundances %>%
                              group_by(Year,Site_code,Location,Mesh,Leaf) %>% 
                              summarise_at(vars("AFDM","Gammaridae"), mean, na.rm=T))

## Spread values by location
Gammaridae_mu_wide <- Taxa_abundances_mu[,c(1:3,7)]  %>% spread(Location, Gammaridae)
colnames(Gammaridae_mu_wide) <- c("Year","Site","Gammarid_D","Gammarid_U1","Gammarid_U2")
Coarse_D_U1_out <- merge(Coarse_D_U1_out,Gammaridae_mu_wide,by=c("Site","Year"))

##*************************************************************
## FOOD PREFERENCE - CWM trait abundances
## Initial plot
plot(Est ~ log1p(Food_3_U1), data=Coarse_D_U1_out)
abline(lm(Est ~ log1p(Food_3_U1), data=Coarse_D_U1_out),col="red")

## Initial analysis - linear regression with term for year
lm1 <- lm(Est ~ log1p(Food_3_U1) + Year, data=Coarse_D_U1_out)
tab_model(lm1)
BIC(lm1)
anova(lm1)

## Create summary object
Food_preference_U1_D_Allyrs <- summary(lm1)
## Write summary table
write.csv(Food_preference_U1_D_Allyrs$coefficients,"Food_preference_U1_D_Allyrs_Linear_summary.csv",row.names = T)
## Write GOF table
write.csv(data.frame(Food_preference_U1_D_Allyrs$r.squared,Food_preference_U1_D_Allyrs$adj.r.squared)
          ,"Food_preference_U1_D_Allyrs_Linear_Rsq.csv",row.names = T)
## Write ANOVA table
write.csv(anova(lm1),"Food_preference_U1_D_Allyrs_Linear_ANOVA.csv",row.names = T)

##*************************************************************
## FEEDING MODE PREFERENCE - CWM trait abundances
## Initial plot
plot(Est ~ log1p(Feeding_3_U1), data=Coarse_D_U1_out)
abline(lm(Est ~ log1p(Feeding_3_U1), data=Coarse_D_U1_out),col="red")

## Initial analysis - linear regression with term for year
lm2 <- lm(Est ~ log1p(Feeding_3_U1) + Year, data=Coarse_D_U1_out)
tab_model(lm2)
BIC(lm2)
anova(lm2)
## Create summary object
Feeding_preference_U1_D_Allyrs <- summary(lm2)
## Write summary table
write.csv(Feeding_preference_U1_D_Allyrs$coefficients,"Feeding_preference_U1_D_Allyrs_Linear_summary.csv",row.names = T)
## Write GOF table
write.csv(data.frame(Feeding_preference_U1_D_Allyrs$r.squared,Feeding_preference_U1_D_Allyrs$adj.r.squared)
          ,"Feeding_preference_U1_D_Allyrs_Linear_Rsq.csv",row.names = T)
## Write ANOVA table
write.csv(anova(lm2),"Feeding_preference_U1_D_Allyrs_Linear_ANOVA.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                            CHANGE VS DETRITIVORES (U2 and U1)                           ##  
##                                      NULL MODEL                                         ##
##                                                                                         ##  
#############################################################################################

## Spread values by location
CWM_Food_3_mu_wide <- CWM_traits_mu[,c(1:3,6)]  %>% spread(Location, Food_3)
colnames(CWM_Food_3_mu_wide) <- c("Year","Site","Food_3_D","Food_3_U1","Food_3_U2")
Coarse_U1_U2_out <- merge(Coarse_U1_U2_out,CWM_Food_3_mu_wide,by=c("Site","Year"))

## Spread values by location
CWM_Feeding_3_mu_wide <- CWM_traits_mu[,c(1:3,7)]  %>% spread(Location, Feeding_3)
colnames(CWM_Feeding_3_mu_wide) <- c("Year","Site","Feeding_3_D","Feeding_3_U1","Feeding_3_U2")
Coarse_U1_U2_out <- merge(Coarse_U1_U2_out,CWM_Feeding_3_mu_wide,by=c("Site","Year"))

## Spread values by location
Gammaridae_mu_wide <- Taxa_abundances_mu[,c(1:3,7)]  %>% spread(Location, Gammaridae)
colnames(Gammaridae_mu_wide) <- c("Year","Site","Gammarid_D","Gammarid_U1","Gammarid_U2")
Coarse_U1_U2_out <- merge(Coarse_U1_U2_out,Gammaridae_mu_wide,by=c("Site","Year"))

##*************************************************************
## FOOD PREFERENCE - CWM trait abundances
## Initial plot
plot(Est ~ log1p(Food_3_U2), data=Coarse_U1_U2_out)
abline(lm(Est ~ log1p(Food_3_U2), data=Coarse_U1_U2_out),col="red")

## Initial analysis - linear regression with term for year
lm3 <- lm(Est ~ log1p(Food_3_U2) + Year, data=Coarse_U1_U2_out)
tab_model(lm3)
BIC(lm3)
anova(lm3)

## Create summary object
Food_preference_U2_U1_Allyrs <- summary(lm3)
## Write summary table
write.csv(Food_preference_U2_U1_Allyrs$coefficients,"Food_preference_U2_U1_Allyrs_Linear_summary.csv",row.names = T)
## Write GOF table
write.csv(data.frame(Food_preference_U2_U1_Allyrs$r.squared,Food_preference_U2_U1_Allyrs$adj.r.squared)
          ,"Food_preference_U2_U1_Allyrs_Linear_Rsq.csv",row.names = T)
## Write ANOVA table
write.csv(anova(lm3),"Food_preference_U2_U1_Allyrs_Linear_ANOVA.csv",row.names = T)

##*************************************************************
## FEEDING MODE PREFERENCE - CWM trait abundances
## Initial plot
plot(Est ~ log1p(Feeding_3_U2), data=Coarse_U1_U2_out)
abline(lm(Est ~ log1p(Feeding_3_U2), data=Coarse_U1_U2_out),col="red")

## Initial analysis - linear regression with term for year
lm4 <- lm(Est ~ log1p(Feeding_3_U2) + Year, data=Coarse_U1_U2_out)
tab_model(lm4)
BIC(lm4)
anova(lm4)
## Create summary object
Feeding_preference_U2_U1_Allyrs <- summary(lm4)
## Write summary table
write.csv(Feeding_preference_U2_U1_Allyrs$coefficients,"Feeding_preference_U2_U1_Allyrs_Linear_summary.csv",row.names = T)
## Write GOF table
write.csv(data.frame(Feeding_preference_U2_U1_Allyrs$r.squared,Feeding_preference_U2_U1_Allyrs$adj.r.squared)
          ,"Feeding_preference_U2_U1_Allyrs_Linear_Rsq.csv",row.names = T)
## Write ANOVA table
write.csv(anova(lm4),"Feeding_preference_U2_U1_Allyrs_Linear_ANOVA.csv",row.names = T)

###*************************************************************************************
### Create figure showing context-dependent change - Food preferences
###*************************************************************************************

## Load packages to access images
library(png)
library(grid)
library(RCurl)

## Invertebrates
stonefly_url  <-"http://phylopic.org/assets/images/submissions/97804f49-f356-47d1-88e2-708855e4140d.original.png"
amphipod_image <-  readPNG("Amphipod.png")

## Get png
stonefly_logo <-  readPNG(getURLContent(stonefly_url))

## Turn into raster
st <- rasterGrob(stonefly_logo , interpolate=TRUE)# stonefly
am <- rasterGrob(amphipod_image, interpolate=TRUE)# amphipod

## Create publication quality figure
png(filename="Effect_size_food3_allyrs_comparison.png", 
    type="cairo",
    units="in", 
    width=11, 
    height=6, 
    pointsize=16, 
    res=600)

plot1 <-ggplot(Coarse_D_U1_out, aes(x=log1p(Food_3_U1))) +
  annotation_custom(st, xmin = 1.11, xmax = 1.65, ymin = 1.95, ymax = 6) +
  annotation_custom(am, xmin = 0.90, xmax = 1.25, ymin = 2.75, ymax = 6.3) +
  geom_hline(yintercept = 0, lty=3) +
  geom_hline(yintercept = 0, lty=3) +
  stat_smooth(aes(y = Est), method = "lm", formula = y ~ x,
              size = 1, colour = "black") +
  theme(axis.line=element_line(colour = "black", size = 0.5, 
                               linetype = "solid")) +
  geom_pointrange(aes(x=log1p(Food_3_U1),y=Est,ymin=CI_lower,ymax=CI_upper)) +
  geom_point(aes(x=log1p(Food_3_U1), y=Est), 
             fill = "white", 
             shape=21, size=5.5) +
  geom_point(aes(x=log1p(Food_3_U1), y=Est), 
             fill = "firebrick", 
             shape=21, size=5.5, alpha=0.7) +
  labs(x="",
    y="Effect of wastewater on decomposition") +
  expand_limits(y=c(-10,6), x=c(-0.01,1.51)) +
  theme_bw() +
  theme( 
    axis.title.y=element_text(size = 18, vjust=0.5),
    axis.title.x=element_text(size = 18, vjust=-0.5),
    axis.text.y=element_text(size = 14),
    axis.text.x=element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(), 
    text = element_text(size=11))

plot2 <-ggplot(Coarse_U1_U2_out, aes(x=log1p(Food_3_U2))) +
  annotation_custom(st, xmin = 1.11, xmax = 1.65, ymin = 1.95, ymax = 6) +
  annotation_custom(am, xmin = 0.90, xmax = 1.25, ymin = 2.75, ymax = 6.3) +
  geom_hline(yintercept = 0, lty=3) +
  stat_smooth(aes(y = Est), method = "lm", formula = y ~ x,
              size = 1, colour = "black") +
  theme(axis.line=element_line(colour = "black", size = 0.5, 
                               linetype = "solid")) +
  geom_pointrange(aes(x=log1p(Food_3_U2),y=Est,ymin=CI_lower,ymax=CI_upper)) +
  geom_point(aes(x=log1p(Food_3_U2), y=Est), 
             fill = "white", 
             shape=21, size=5.5) +
  geom_point(aes(x=log1p(Food_3_U2), y=Est), 
             fill = "royalblue", 
             shape=21, size=5.5, alpha=0.7) +
  labs(x="",
       y="Effect of location on decomposition") +
  expand_limits(y=c(-10,6), x=c(-0.01,1.51)) +
  theme_bw() +
  theme(
    axis.title.y=element_text(size = 18, vjust=0.5),
    axis.title.x=element_text(size = 18, vjust=-0.5),
    axis.text.y=element_text(size = 14),
    axis.text.x=element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(), 
    text = element_text(size=11))

grid.newpage()

require(cowplot)
cowplot::plot_grid(plot1, plot2, align = 'v')

grid.text(expression(paste(italic(F),""["1,17"]," = 4.55, ",italic(P)," < ",0.05,sep="")), 
          x = unit(0.1, "npc"), y = unit(0.25, "npc"), 
          gp=gpar(fontsize=16, col="black"), just = "left")
grid.text(expression(paste(italic(R^2),""[adj]," = 0.177",sep="")), 
          x = unit(0.1, "npc"), y = unit(0.2, "npc"), 
          gp=gpar(fontsize=16, col="black"), 
          just = "left")

grid.text(expression(paste(italic(F),""["1,17"]," = 0.035, ",italic(P)," = ",0.85,sep="")), 
          x = unit(0.6, "npc"), y = unit(0.25, "npc"), 
          gp=gpar(fontsize=16, col="black"), just = "left")
grid.text(expression(paste(italic(R^2),""[adj]," = 0.111",sep="")), 
          x = unit(0.6, "npc"), y = unit(0.2, "npc"), 
          gp=gpar(fontsize=16, col="black"), 
          just = "left")

grid.text(expression(paste(italic(log),"(",italic(x),"+1) CWM Detritivores")),
          x = unit(0.5, "npc"), y = unit(0.02, "npc"), 
          gp=gpar(fontsize=18, col="black"), just = "centre")

grid.text("a", x = unit(0.025, "npc"), y = unit(0.97, "npc"),gp=gpar(fontsize=20, fontface="bold", col="black"))
grid.text("b", x = unit(0.525, "npc"), y = unit(0.97, "npc"),gp=gpar(fontsize=20, fontface="bold", col="black"))

dev.off()

###*************************************************************************************
### Create figure showing context-dependent change - Feeding preferences
###*************************************************************************************

## Create publication quality figure
png(filename="Effect_size_feeding3_allyrs_comparison.png", 
    type="cairo",
    units="in", 
    width=11, 
    height=6, 
    pointsize=16, 
    res=600)

plot1 <-ggplot(Coarse_D_U1_out, aes(x=log1p(Feeding_3_U1))) +
  annotation_custom(st, xmin = 1.11, xmax = 1.65, ymin = 1.95, ymax = 6) +
  annotation_custom(am, xmin = 0.90, xmax = 1.25, ymin = 2.75, ymax = 6.3) +
  geom_hline(yintercept = 0, lty=3) +
  geom_hline(yintercept = 0, lty=3) +
  stat_smooth(aes(y = Est), method = "lm", formula = y ~ x,
              size = 1, colour = "black") +
  theme(axis.line=element_line(colour = "black", size = 0.5, 
                               linetype = "solid")) +
  geom_pointrange(aes(x=log1p(Feeding_3_U1),y=Est,ymin=CI_lower,ymax=CI_upper)) +
  geom_point(aes(x=log1p(Feeding_3_U1), y=Est), 
             fill = "white", 
             shape=21, size=5.5) +
  geom_point(aes(x=log1p(Feeding_3_U1), y=Est), 
             fill = "firebrick", 
             shape=21, size=5.5, alpha=0.7) +
  labs(x="",
       y="Effect of wastewater on decomposition") +
  expand_limits(y=c(-10,6), x=c(-0.01,1.51)) +
  theme_bw() +
  theme(
    axis.title.y=element_text(size = 18, vjust=0.5),
    axis.title.x=element_text(size = 18, vjust=-0.5),
    axis.text.y=element_text(size = 14),
    axis.text.x=element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(), 
    text = element_text(size=11))

plot2 <-ggplot(Coarse_U1_U2_out, aes(x=log1p(Feeding_3_U2))) +
  annotation_custom(st, xmin = 1.11, xmax = 1.65, ymin = 1.95, ymax = 6) +
  annotation_custom(am, xmin = 0.90, xmax = 1.25, ymin = 2.75, ymax = 6.3) +
  geom_hline(yintercept = 0, lty=3) +
  stat_smooth(aes(y = Est), method = "lm", formula = y ~ x,
              size = 1, colour = "black") +
  theme(axis.line=element_line(colour = "black", size = 0.5, 
                               linetype = "solid")) +
  geom_pointrange(aes(x=log1p(Feeding_3_U2),y=Est,ymin=CI_lower,ymax=CI_upper)) +
  geom_point(aes(x=log1p(Feeding_3_U2), y=Est), 
             fill = "white", 
             shape=21, size=5.5) +
  geom_point(aes(x=log1p(Feeding_3_U2), y=Est), 
             fill = "royalblue", 
             shape=21, size=5.5, alpha=0.7) +
  labs(x="",
       y="Effect of location on decomposition") +
  expand_limits(y=c(-10,6), x=c(-0.01,1.51)) +
  theme_bw() +
  theme(
    axis.title.y=element_text(size = 18, vjust=0.5),
    axis.title.x=element_text(size = 18, vjust=-0.5),
    axis.text.y=element_text(size = 14),
    axis.text.x=element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.background = element_blank(), 
    text = element_text(size=11))

grid.newpage()

require(cowplot)

cowplot::plot_grid(plot1, plot2, align = 'v')

grid.text(expression(paste(italic(F),""["1,17"]," = 3.51, ",italic(P)," = ",0.078,sep="")), 
          x = unit(0.1, "npc"), y = unit(0.25, "npc"), 
          gp=gpar(fontsize=16, col="black"), just = "left")
grid.text(expression(paste(italic(R^2),""[adj]," = 0.168",sep="")), 
          x = unit(0.1, "npc"), y = unit(0.2, "npc"), 
          gp=gpar(fontsize=16, col="black"), 
          just = "left")

grid.text(expression(paste(italic(F),""["1,17"]," = 0.079, ",italic(P)," = ",0.782,sep="")), 
          x = unit(0.6, "npc"), y = unit(0.25, "npc"), 
          gp=gpar(fontsize=16, col="black"), just = "left")
grid.text(expression(paste(italic(R^2),""[adj]," = 0.111",sep="")), 
          x = unit(0.6, "npc"), y = unit(0.2, "npc"), 
          gp=gpar(fontsize=16, col="black"), 
          just = "left")

grid.text(expression(paste(italic(log),"(",italic(x),"+1) CWM Detritivores")),
          x = unit(0.5, "npc"), y = unit(0.02, "npc"), 
          gp=gpar(fontsize=18, col="black"), just = "centre")

grid.text("a", x = unit(0.025, "npc"), y = unit(0.97, "npc"),gp=gpar(fontsize=20, fontface="bold", col="black"))
grid.text("b", x = unit(0.525, "npc"), y = unit(0.97, "npc"),gp=gpar(fontsize=20, fontface="bold", col="black"))

dev.off()

#############################################################################################
##                                                                                         ##  
##                            TEST IMPACT VS NULL MODEL                                    ##   
##                                MIXED MODEL APPROACH                                     ##         
##                                                                                         ##  
#############################################################################################

###******************************************************************************************
### FOOD_3
### Perform mixed model
head(Coarse_D_U1_out)
head(Coarse_U1_U2_out)

## Rationalise D_U1
D_U1_out <- Coarse_D_U1_out[,c(1:3,5,10)]
names(D_U1_out)[5]<-"Food_3"

## Rationalise U1_U2
U1_U2_out <- Coarse_U1_U2_out[,c(1:3,5,11)]
names(U1_U2_out)[5]<-"Food_3"

## Combine data
data_combined <- rbind(D_U1_out,U1_U2_out)

## Test slope (interaction) whilst accounting for year (Blocked) and Site (Random effect)
LME1 <- blmer(Est ~ Food_3:Contrast + Year + (1|Site), data = data_combined, REML=T)

## Check output
summary(LME1)

## Check ANOVA
car::Anova(LME1,type="III") ## signficant

## Write table for coefficients and p-values
ctable <- coef(summary(LME1))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(LME1), lower.tail = F)
LME1_coefficients <- cbind(ctable, pvals)
write.csv(LME1_coefficients,"Food_preference_U1_D_Allyrs_LME_coefficients.csv",row.names = T)

## Write confidence intervals for parameter estimates
LME1_conf_intervals <- confint(LME1, level = 0.95)
write.csv(LME1_conf_intervals,"Food_preference_U1_D_Allyrs_LME_confintervals.csv",row.names = T)

## Write ANOVA table
LME1_anova <- anova(LME1)
LME1_anova$p_value <- NA
LME1_anova$p_value[1] <- pf(LME1_anova[1,4], LME1_anova[1,1], df.residual(LME1), lower.tail = FALSE)
LME1_anova$p_value[2] <- pf(LME1_anova[2,4], LME1_anova[2,1], df.residual(LME1), lower.tail = FALSE)
write.csv(LME1_anova,"Food_preference_U1_D_Allyrs_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(LME1)
n <- nrow(X)
Beta <- fixef(LME1)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(LME1)
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
## Fixed 23.2%
## Random 17.5% 
## All 40.7%
write.csv(Model_variance,"Food_preference_U1_D_Allyrs_LME_variance.csv",row.names = T)

require(phia)
## This will test the slope for WW impact vs the Null model
LME1_post_output <- testInteractions(LME1, pairwise="Contrast", slope="Food_3", adjust="BH")
write.csv(LME1_post_output ,"Food_preference_U1_D_Allyrs_LME_slopes_comparison.csv",row.names = T)

###******************************************************************************************
### FEEDING_3
### Perform mixed model
head(Coarse_D_U1_out)
head(Coarse_U1_U2_out)

## Rationalise D_U1
D_U1_out <- Coarse_D_U1_out[,c(1:3,5,13)]
names(D_U1_out)[5]<-"Feeding_3"

## Rationalise U1_U2
U1_U2_out <- Coarse_U1_U2_out[,c(1:3,5,14)]
names(U1_U2_out)[5]<-"Feeding_3"

## Combine data
data_combined <- rbind(D_U1_out,U1_U2_out)

## Test slope (interaction) whilst accounting for year (Blocked) and Site (Random effect)
LME2 <- blmer(Est ~ Feeding_3:Contrast + Year + (1|Site), data = data_combined, REML=T)

## Check output
summary(LME2)

## Check ANOVA
car::Anova(LME2,type="III")

## Write table for coefficients and p-values
ctable <- coef(summary(LME2))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(LME2), lower.tail = F)
LME2_coefficients <- cbind(ctable, pvals)
write.csv(LME2_coefficients,"Feeding_preference_U1_D_Allyrs_LME_coefficients.csv",row.names = T)

## Write confidence intervals for parameter estimates
LME2_conf_intervals <- confint(LME2, level = 0.95)
write.csv(LME2_conf_intervals,"Feeding_preference_U1_D_Allyrs_LME_confintervals.csv",row.names = T)

## Write ANOVA table
LME2_anova <- anova(LME2)
LME2_anova$p_value <- NA
LME2_anova$p_value[1] <- pf(LME2_anova[1,4], LME2_anova[1,1], df.residual(LME2), lower.tail = FALSE)
LME2_anova$p_value[2] <- pf(LME2_anova[2,4], LME2_anova[2,1], df.residual(LME2), lower.tail = FALSE)
write.csv(LME2_anova,"Feeding_context_depend_LME2_LME_ANOVA_200814.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(LME2)
n <- nrow(X)
Beta <- fixef(LME2)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(LME2)
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
## Fixed 21.9%
## Random 17.9% 
## All 39.9%
write.csv(Model_variance,"Feeding_preference_U1_D_Allyrs_LME_variance.csv",row.names = T)

## This will test the slope for WW impact vs the Null model
LME2_post_output <- testInteractions(LME2, pairwise="Contrast", slope="Feeding_3", adjust="BH")
write.csv(LME2_post_output ,"Feeding_preference_U1_D_Allyrs_LME_slopes_comparison.csv",row.names = T)

#############################################################################################
##                                                                                         ##  
##                                GAMMARID ABUNDANCES                                      ##
##                              EST IMPACT VS NULL MODEL                                   ##   
##                      LINEAR REGRESSION and MIXED MODEL APPROACH                         ##         
##                                                                                         ##  
#############################################################################################

##*************************************************************
## GAMMARID ABUNDANCES
## Initial plot
plot(Est ~ log1p(Gammarid_U1), data=Coarse_D_U1_out)
abline(lm(Est ~ log1p(Gammarid_U1), data=Coarse_D_U1_out),col="red")

## Initial analysis - linear regression with term for year
lm5 <- lm(Est ~ log1p(Gammarid_U1) + Year, data=Coarse_D_U1_out)
tab_model(lm5)
BIC(lm5)
anova(lm5)

## Create summary object
Gammarid_abundance_U1_D_Allyrs <- summary(lm5)
## Write summary table
write.csv(Gammarid_abundance_U1_D_Allyrs$coefficients,"Gammarid_abundance_U1_D_Allyrs_Linear_summary.csv",row.names = T)
## Write GOF table
write.csv(data.frame(Gammarid_abundance_U1_D_Allyrs$r.squared,Gammarid_abundance_U1_D_Allyrs$adj.r.squared)
          ,"Gammarid_abundance_U1_D_Allyrs_Linear_Rsq.csv",row.names = T)
## Write ANOVA table
write.csv(anova(lm5),"Gammarid_abundance_U1_D_Allyrs_Linear_ANOVA.csv",row.names = T)

##*************************************************************
## GAMMARID ABUNDANCES
## Initial plot
plot(Est ~ log1p(Gammarid_U2), data=Coarse_U1_U2_out)
abline(lm(Est ~ log1p(Gammarid_U2), data=Coarse_U1_U2_out),col="red")

## Initial analysis - linear regression with term for year
lm6 <- lm(Est ~ log1p(Gammarid_U2) + Year, data=Coarse_U1_U2_out)
tab_model(lm6)
BIC(lm6)
anova(lm6)

## Create summary object
Gammarid_abundance_U2_U1_Allyrs <- summary(lm6)
## Write summary table
write.csv(Gammarid_abundance_U2_U1_Allyrs$coefficients,"Gammarid_abundance_U2_U1_Allyrs_Linear_summary.csv",row.names = T)
## Write GOF table
write.csv(data.frame(Gammarid_abundance_U2_U1_Allyrs$r.squared,Gammarid_abundance_U2_U1_Allyrs$adj.r.squared)
          ,"Gammarid_abundance_U2_U1_Allyrs_Linear_Rsq.csv",row.names = T)
## Write ANOVA table
write.csv(anova(lm6),"Gammarid_abundance_U2_U1_Allyrs_Linear_ANOVA.csv",row.names = T)

##*************************************************************
## GAMMARID ABUNDANCES
head(Coarse_D_U1_out)
head(Coarse_U1_U2_out)

## Rationalise D_U1
D_U1_out <- Coarse_D_U1_out[,c(1:3,5,16)]
names(D_U1_out)[5]<-"Gammarid"

## Rationalise U1_U2
U1_U2_out <- Coarse_U1_U2_out[,c(1:3,5,17)]
names(U1_U2_out)[5]<-"Gammarid"

## Combine data
data_combined <- rbind(D_U1_out,U1_U2_out)

## Test slope (interaction) whilst accounting for year (Blocked) and Site (Random effect)
LME3 <- blmer(Est ~ Gammarid:Contrast + Year + (1|Site), data = data_combined, REML=T)

## Check output
summary(LME3)

## Check ANOVA
car::Anova(LME3,type="III")

## Write table for coefficients and p-values
ctable <- coef(summary(LME3))
pvals <- 2 * pt(abs(ctable[, "t value"]), df.residual(LME3), lower.tail = F)
LME3_coefficients <- cbind(ctable, pvals)
write.csv(LME3_coefficients,"Gammarid_abundance__U1_D_Allyrs_LME_coefficients.csv",row.names = T)

## Write confidence intervals for parameter estimates
LME3_conf_intervals <- confint(LME3, level = 0.95)
write.csv(LME3_conf_intervals,"Gammarid_abundance__U1_D_Allyrs_LME_confintervals.csv",row.names = T)

## Write ANOVA table
LME3_anova <- anova(LME3)
LME3_anova$p_value <- NA
LME3_anova$p_value[1] <- pf(LME3_anova[1,4], LME3_anova[1,1], df.residual(LME3), lower.tail = FALSE)
LME3_anova$p_value[2] <- pf(LME3_anova[2,4], LME3_anova[2,1], df.residual(LME3), lower.tail = FALSE)
write.csv(LME3_anova,"Gammarid_abundance__U1_D_Allyrs_LME_ANOVA.csv",row.names = T)

## Estimate marginal and conditional R-squares
X <- model.matrix(LME3)
n <- nrow(X)
Beta <- fixef(LME3)
Sf <- var(X %*% Beta)
Sigma.list <- VarCorr(LME3)
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
## Fixed 7.0%
## Random 14.2% 
## All 21.2%
write.csv(Model_variance,"Gammarid_abundance__U1_D_Allyrs_LME_variance.csv",row.names = T)

## This will test the slope for WW impact vs the Null model
LME3_post_output <- testInteractions(LME3, pairwise="Contrast", slope="Gammarid", adjust="BH")
write.csv(LME3_post_output ,"Gammarid_abundance__U1_D_Allyrs_LME_slopes_comparison.csv",row.names = T)



