############################################################################################
##                                                                                        ##  
##                                LEAFPACK INVERTEBRATES                                  ##  
##                                                                                        ##
############################################################################################

##  Script: ECOIMPACT leafpack assay invertebrate data preparation
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Prepare invertebrate data for analysis by calculating CWM trait abundances
##  2. Explore underlying patterns in data supporting analyses used in the Main Text

##***************************************************************************************
##Load libraries you need
library(reshape)
library(reshape2)
library(plyr)
library(tidyverse)
library(readr)

#Load packages needed for trait analyses
require(magic)
require(geometry)
require(ade4)
library(FD)

##***************************************************************************************
## Load data - site information
sites <- read_csv("1_Input_data/DATA_sites_locations_names_200802.csv")

## Load data - all data for leafpacks
leafpack <- read_csv("1_Input_data/DATA_Taxa_LPA_updated_200802.csv")
lpa_taxa <- merge(leafpack,sites[,c(7,2)],by="Site")
lpa_taxa <- as.data.frame(lpa_taxa[,c(2,28,3:27)])

############################################################################################
##                                                                                        ##  
##                          LEAFPACK INVERTEBRATES TRAITS                                 ##  
##                                                                                        ##
############################################################################################

LPA_taxa_traits_raw_200803 <- read_csv("1_Input_data/LPA_taxa_traits_raw_200803.csv")

colnames(LPA_taxa_traits_raw_200803)

require(tidyverse)
## Create mean values for each taxa (Family) from the Tachet trait database
## This is important as it undepins the trait-based values used in subsequent analyses
Traits_mu <- data.frame(LPA_taxa_traits_raw_200803 %>%
                          group_by(Taxa_code) %>% 
                          summarise_at(vars("Size_1","Size_2","Size_3","Size_4","Size_5","Size_6","Size_7",
                                            "Lifecycle_1","Lifecycle_2","Fecundity_1","Fecundity_2","Fecundity_3","Lifestages_1","Lifestages_2",
                                            "Lifestages_3","Lifestages_4","Reproduction_1","Reproduction_2","Reproduction_3",
                                            "Reproduction_4","Reproduction_5","Reproduction_6","Reproduction_7","Reproduction_8",
                                            "Dispersal_1","Dispersal_2","Dispersal_3","Dispersal_4","Resistance_1",
                                            "Resistance_2","Resistance_3","Resistance_4","Resistance_5","Respiration_1",
                                            "Respiration_2","Respiration_3","Respiration_4","Respiration_5","Movement_1",
                                            "Movement_2","Movement_3","Movement_4","Movement_5","Movement_6",
                                            "Movement_7","Movement_8","Food_1","Food_2","Food_3",
                                            "Food_4","Food_5","Food_6","Food_7","Food_8",
                                            "Food_9","Feeding_1","Feeding_2","Feeding_3","Feeding_4",
                                            "Feeding_5","Feeding_6","Feeding_7","Feeding_8","Transversal_1",
                                            "Transversal_2","Transversal_3","Transversal_4","Transversal_5","Transversal_6",
                                            "Transversal_7","Longitudinal_1","Longitudinal_2","Longitudinal_3","Longitudinal_4",
                                            "Longitudinal_5","Longitudinal_6","Longitudinal_7","Longitudinal_8","Elevation_1",
                                            "Elevation_2","Elevation_3","Substrate_1","Substrate_2","Substrate_3",
                                            "Substrate_4","Substrate_5","Substrate_6","Substrate_7","Substrate_8",
                                            "Substrate_9","Current_1","Current_2","Current_3","Current_4",
                                            "Trophic_1","Trophic_2","Trophic_3","Salinity_1","Salinity_2",
                                            "Temperature_1","Temperature_2","Temperature_3","Saprobity_1","Saprobity_2",
                                            "Saprobity_3","Saprobity_4","Saprobity_5","pH_1","pH_2",
                                            "pH_3","pH_4","pH_5","pH_6"),mean,na.rm=T))

write.csv(Traits_mu,"lpa_taxa_traits_mean_data.csv",row.names = F)

############################################################################################
##                                                                                        ##  
##                            LEAFPACK INVERTEBRATES - ALDER                              ##  
##                                                                                        ##
############################################################################################

## Rationalise data - need Coarse mesh Alder leaves
lpa_alder_coarse_new <- lpa_taxa[which(lpa_taxa$Mesh=="Coarse" & lpa_taxa$Leaf=="Alder"),]

## Rationalise data - remove NAs (where no AFDM data available)
lpa_alder_coarse_new <- lpa_alder_coarse_new[!is.na(lpa_alder_coarse_new$AFDM),]

## Combine taxa synonymns (Limnephillidae + Trichoptera)
lpa_alder_coarse_new$Limnephillid <- lpa_alder_coarse_new$Limnephillidae+lpa_alder_coarse_new$Trichoptera

## Remove taxa synonyms
lpa_alder_coarse_new <- lpa_alder_coarse_new[,-c(10,17)]

## Check output
head(lpa_alder_coarse_new)

## Re-label headers to ensure taxa families are correct
colnames(lpa_alder_coarse_new) <- c("Year","Site_code","Location","Mesh","Leaf","Rep","ID","AFDM",
                                    "Gammaridae","Baetidae","Nemouridae","Caloptyerigidae",
                                    "Chironomidae","Simuliidae","Asellidae","Rhyacophilidae","Elmidae","Dytiscidae",   
                                    "Tipulidae","Planorbidae","Lymnaeidae", "Oligochaete","Hirudinea","Acari","Other",
                                    "Limnephilidae") 
## Subset and sort columns by Family
taxa <- lpa_alder_coarse_new[,c(9:26)] 
taxa <- taxa[,order(colnames(taxa))]

## Recombine with site information and write output
invert_abundance_out <- data.frame(lpa_alder_coarse_new[,-c(9:26)],taxa)
write.csv(invert_abundance_out,"lpa_alder_coarse_allyrs_invert_abundances_data.csv",row.names = F)

##***************************************************************************
## Running "FD" package in R
## See: Laliberte, E., Legendre, P., and B. Shipley. (2014). FD: measuring functional
## diversity from multiple traits, and other tools for functional ecology. R package
## version 1.0-12.
##***************************************************************************

## For FRic, FEve, and FDiv see:
## Villeger, S., et al. (2008). New multidimensional functional diversity indices for a 
## multifaceted framework in functional ecology. Ecology 89(8): 2290-2301.

## For FDis see:
## Laliberte, E. and P. Legendre (2010). A distance-based framework for measuring functional 
## diversity from multiple traits." Ecology 91(1): 299-305.

head(invert_abundance_out)

## Remove unidentified taxa "Other"
invert_abund_out <- invert_abundance_out[,-22]

## Remove leafpacks with no observations
invert_abund_out <- invert_abund_out[rowSums(invert_abund_out[-(1:8)]) !=0, ]

## Remove abundances row header for FD function
abundances <- invert_abund_out[,-c(1:8)]

## Remove taxa with no observations
abundances <- abundances[,which(colSums(abundances) !=0)]## Remove missing taxa

## Remove taxa with no observations from trait matrix
Traits <- Traits_mu[-c(which(Traits_mu$Taxa_code=="Acari"),
                       which(Traits_mu$Taxa_code=="Dytiscidae"),
                       which(Traits_mu$Taxa_code=="Elmidae"),
                       which(Traits_mu$Taxa_code=="Simuliidae")),]

## Remove abundances row header for FD function
rownames(Traits) <- Traits[,1]
traits<-Traits[,-1]

## Calculate FD for all
FD_all <- dbFD(traits, abundances, w.abun = TRUE, stand.x = TRUE, ord = "metric", 
               asym.bin = NULL, corr = "sqrt", 
               calc.FRic = TRUE, m = "max", stand.FRic = T, scale.RaoQ = FALSE, 
               calc.FGR = TRUE, clust.type = "ward", km.inf.gr = 2, km.sup.gr = nrow(x) - 1, 
               km.iter = 100, km.crit = "calinski", calc.CWM = TRUE, CWM.type = "dom", 
               calc.FDiv = TRUE, dist.bin = 2, print.pco = T, messages = TRUE)

#g - separate by group
#8 - number of functional groups

## Extract FD estimates
FD_metrics_all_invert   <- data.frame(FD_all$nbsp,
                                      FD_all$sing.sp,
                                      FD_all$qual.FRic,
                                      FD_all$FRic,
                                      FD_all$FEve,
                                      FD_all$FDiv,
                                      FD_all$FDis,
                                      FD_all$RaoQ,
                                      FD_all$FGR)
## Relabel FD estimates
col_headings <- c('nbsp','sing.sp',"qual.FRic",'FRic','FEve','FDiv','FDis','RaoQ','FGR')
names(FD_metrics_all_invert) <- col_headings

## Write FD estimates to file
write_data <- data.frame(invert_abund_out[,c(1:8)],FD_metrics_all_invert)
write.csv(write_data,file = "lpa_alder_coarse_allyrs_FD_metrics_raw_data.csv",row.names = F)

## Create mean values for each location
FD_mu <- data.frame(write_data[,-c(4:7)] %>%
                      group_by(Year,Site_code,Location) %>% 
                      summarise_at(vars("AFDM","nbsp","sing.sp","qual.FRic",
                                        "FRic","FEve","FDiv","FDis","RaoQ","FGR"),mean,na.rm=T))
## Convert NA to zero
FD_mu[is.na(FD_mu)] <- 0
write.csv(FD_mu,"lpa_alder_coarse_allyrs_FD_metrics_mean_data.csv",row.names = F)

#load "pheatmap" library
library(pheatmap)

#calculate the pearson correlation coefficient matrix
myMat.cor <- cor(FD_mu[,-c(1:3,7)], method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)
dev.off()

## Note potential for postive functional diversity effects
plot(sqrt(AFDM)~log1p(FEve),data=FD_mu)
abline(lm(sqrt(AFDM)~log1p(FEve),data=FD_mu),col="red")
summary(lm(sqrt(AFDM)~log1p(FEve),data=FD_mu))

## Note potential for postive functional diversity effects FDis looks promising
plot(sqrt(AFDM)~FDis,data=FD_mu)
abline(lm(sqrt(AFDM)~FDis,data=FD_mu),col="red")
summary(lm(sqrt(AFDM)~FDis,data=FD_mu))

##******************************************************************************************
## COMMUNITY-WEIGHTED MEANS
## See: Lavorel, S. and E. Garnier (2002). Predicting changes in community composition and 
## ecosystem functioning from plant traits: revisiting the Holy Grail. 
## Functional Ecology 16(5): 545-556.

## Write CMW (community weighted means) estimates to file
write_data <- data.frame(invert_abund_out[,c(1:8)],FD_all$CWM)
write.csv(write_data,file = "lpa_alder_coarse_allyrs_CWM_traits_data.csv",row.names = F)
rm(write_data)

## Explore patterns underpinning CMW trait predictors
CWM_traits <- data.frame(invert_abund_out[,c(1:8)],FD_all$CWM)

## Plot feeding traits and other key traits and Gammarid abundances 
png(filename="lpa_alder_coarse_allyrs_CWM_Gammarid_regressions.png", 
    type="cairo",
    units="in", 
    width=11, 
    height=10, 
    pointsize=16, 
    res=600)

par(mfrow=c(2,2), mar = c(4,4,2,2) + 0.1) 

plot(log(CWM_traits$Food_3)~log1p(invert_abund_out$Gammaridae),
     main = "Food: plant detritus > 1mm",
     ylab = "log Trait abundance",
     xlab = "")
abline(lm(log(CWM_traits$Food_3)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log(CWM_traits$Food_3)~log1p(invert_abund_out$Gammaridae)))
text(3.4,-0.73,bquote(italic(italic(R^2)[adj])*~"= 0.2721"))

plot(log(CWM_traits$Feeding_3)~log1p(invert_abund_out$Gammaridae),
     main = "Feeding habits: shredder",
     ylab = "log Trait abundance",
     xlab = "")
abline(lm(log(CWM_traits$Feeding_3)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log(CWM_traits$Feeding_3)~log1p(invert_abund_out$Gammaridae)))
text(3.4,-0.1,bquote(italic(italic(R^2)[adj])*~"= 0.2848"))

plot(log1p(CWM_traits$Size_5)~log1p(invert_abund_out$Gammaridae),
     main = "Maximal size: >2-4 cm",
     ylab = "log Trait abundance",
     xlab = "log+1 Gammarid abundance")
abline(lm(log1p(CWM_traits$Size_5)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log1p(CWM_traits$Size_5)~log1p(invert_abund_out$Gammaridae)))
text(3.4,0.1,bquote(italic(italic(R^2)[adj])*~"= 0.3045"))

plot(log1p(CWM_traits$Lifestages_4)~log1p(invert_abund_out$Gammaridae),
     main = "Aquatic stages: adult",
     ylab = "log Trait abundance",
     xlab = "log+1 Gammarid abundance")
abline(lm(log1p(CWM_traits$Lifestages_4)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log1p(CWM_traits$Lifestages_4)~log1p(invert_abund_out$Gammaridae)))
text(3.4,0.1,bquote(italic(italic(R^2)[adj])*~"= 0.3838"))

dev.off()

## Assess potential for differences between site locations
## Seems equivocal

CWM_traits[,c(2,3,57)] %>%
  group_by(Location) %>%
  summarize(mean=mean(Food_3),
            sd=sd(Food_3),
            median=median(Food_3))

CWM_traits[,c(2,3,66)] %>%
  group_by(Location) %>%
  summarize(mean=mean(Feeding_3),
            sd=sd(Feeding_3),
            median=median(Feeding_3))

##***************************************************************************************
## Explore invertebrate data to see how dominant gammarids are and who the other players
## are
##***************************************************************************************

## Rationalise data - need Coarse mesh Alder leaves
lpa_alder_coarse_new <- lpa_taxa_new[which(lpa_taxa$Mesh=="Coarse" & lpa_taxa$Leaf=="Alder"),]

## Rationalise data - remove NAs (where no AFDM data available)
lpa_alder_coarse_new <- lpa_alder_coarse_new[!is.na(lpa_alder_coarse_new$AFDM),]

## Combine taxa synonymns (Limnephillidae + Trichoptera)
lpa_alder_coarse_new$Limnephillid <- lpa_alder_coarse_new$Limnephillidae+lpa_alder_coarse_new$Trichoptera

## Remove taxa synonyms
lpa_alder_coarse_new <- lpa_alder_coarse_new[,-c(10,17)]

## Check output
head(lpa_alder_coarse_new)

## Re-label headers to ensure taxa families are correct
colnames(lpa_alder_coarse_new) <- c("Year","Site_code","Location","Mesh","Leaf","Rep","ID","AFDM",
                                    "Gammaridae","Baetidae","Nemouridae","Caloptyerigidae",
                                    "Chironomidae","Simuliidae","Asellidae","Rhyacophilidae","Elmidae","Dytiscidae",   
                                    "Tipulidae","Planorbidae","Lymnaeidae", "Oligochaete","Hirudinea","Acari","Other",
                                    "Limnephilidae") 
## Subset and sort columns by Family
taxa <- lpa_alder_coarse_new[,c(9:26)] 
taxa <- taxa[,order(colnames(taxa))]

## Recombine with site information and write output
invert_abundance_out <- data.frame(lpa_alder_coarse_new[,-c(9:26)],taxa)
row.names(invert_abundance_out) <- NULL 

## Exclude erroneous values
## 2013	NIE	US1	Marta	rep 3 - Too low
## 2013	NIE	US1	Marta	rep 1 - Too low
## 2013	NIE	US1	Marta	rep 2 - Too low
## 2013 MES DS Marta rep 6 - Incorrect initial value or human error - too low
## 2013	ROT	US1	Sandra rep 5 - Hole in bag or human error - way too high
## 2014 KNO DS Marta rep 1 Incorrect initial value or human error - too low (-ve)
## 2014 MAR US2 Marta  rep 1 clearly confounded, buried in fine sediment 

invert_abundance_out <- invert_abundance_out[-c(168,174,197,231:233,288),]

taxa <- invert_abundance_out[,c(9:26)]

##*****************************************************************************
## ABUNDANCES and proportions - to see the main players   
## Sum column abundances
Total_inverts <-  colSums(taxa, na.rm = FALSE, dims = 1)

head(Total_inverts)

## Calculate total proportions for dominant Families
Chironomid_abundance    <- Total_inverts[5]/sum(Total_inverts)*100
Gammarid_abundance    <- Total_inverts[8]/sum(Total_inverts)*100
Limnephilid_abundance <- Total_inverts[10]/sum(Total_inverts)*100
Nemourid_abundance    <- Total_inverts[12]/sum(Total_inverts)*100

##*****************************************************************************
## DENSITIES and relative abundance
## Convert abundances to densities
invert_density <- taxa/invert_abundance_out$AFDM

## Sum column densities
Total_invert_densities <-  colSums (invert_density, na.rm = FALSE, dims = 1)

## Calculate total proportions for dominant Families
Chironomid_abundance    <- Total_invert_densities[5]/sum(Total_invert_densities)*100
Gammarid_abundance    <- Total_invert_densities[8]/sum(Total_invert_densities)*100
Limnephilid_abundance <- Total_invert_densities[10]/sum(Total_invert_densities)*100
Nemourid_abundance    <- Total_invert_densities[12]/sum(Total_invert_densities)*100

## Write file for invertebrate densities
## Recombine with site information and write output
invert_density_out <- data.frame(invert_abundance_out[,-c(9:26)],invert_density)

write.csv(invert_density_out,"lpa_alder_coarse_allyrs_densities_data.csv",row.names = F)

############################################################################################
##                                                                                        ##  
##                     ALDER LEAFPACK INVERTEBRATES COMMUNITIES                           ##  
##                                                                                        ##
############################################################################################

## Load vegan for multivariate analysis
require(vegan)

## Hellinger transform density data and run NMDS
dinvert <- decostand(invert_density, "hel")##Abundance data Hellinger-transformed
dinvert  <- dinvert [,colSums(invert_density)!=0]## Remove missing taxa
mod1 <- metaMDS(dinvert, dist="euc", trymax = 20, autotransform =F)

## Extracts scores for plotting below
nmsc <- scores(mod1, display = "sites")
species <- scores(mod1, display = "species")

#write.csv(data.frame(lpa_alder_coarse[,-c(7:23)],nmsc),"NMDS_sites_hell_200802.csv",row.names = F)
#write.csv(lpa_alder_coarse,"Coarse_alder_invertebrates_200802.csv",row.names = F)

## Compare gammarid densities and NMDS2
gammarid_density <- data.frame(invert_density[,8],nmsc)
names(gammarid_density)[1] <- "Gammarid_density"
head(gammarid_density)

## Compare gammarid abundances and NMDS2
gammarid_abund <- data.frame(taxa[,8],nmsc)
names(gammarid_abund)[1] <- "Gammarid_abundance"
head(gammarid_abund)

## Compare gammarid abundances/densities and NMDS2
png(filename="lpa_alder_coarse_allyrs_NMDS2_gammarid_regressions.png", 
    type="cairo",
    units="in", 
    width=9, 
    height=6, 
    pointsize=13, 
    res=600)

par(mfrow=c(1,2), mar = c(4,4,2,2) + 0.1, pty="s")

plot(NMDS2 ~ log1p(Gammarid_abundance), gammarid_abund, 
     main = "Gammaridae abundance",
     xlab = "log+1 Gammarid abundance")
abline(lm(NMDS2 ~ log1p(Gammarid_abundance), gammarid_abund),col="red")
summary(lm(NMDS2 ~ log1p(Gammarid_abundance), gammarid_abund))
text(3.15,0.7,bquote(italic(italic(R^2)[adj])*~"= 0.5078 "))

plot(NMDS2 ~ log1p(Gammarid_density), gammarid_density, 
     main = "Gammaridae density",
     xlab = "log+1 Gammarid density")
abline(lm(NMDS2 ~ log1p(Gammarid_density), gammarid_density),col="red")
summary(lm(NMDS2 ~ log1p(Gammarid_density), gammarid_density))
text(4.4,0.7,bquote(italic(italic(R^2)[adj])*~"= 0.4102 "))

dev.off()

## Plot contour plot to better undestand the influence of dominant taxa
png(filename="lpa_alder_coarse_allyrs_invert_density_NMDS_contour_plot.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=12, 
    res=600)

par(mfrow=c(2,2), mar = c(4,4,2,2) + 0.1, pty="s")

ordisurf(mod1 ~ Chironomidae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Chironomidae (58%)")
ordisurf(mod1 ~ Gammaridae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Gammaridae (21%)")
ordisurf(mod1 ~ Limnephilidae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Limnephilidae (8%)")
ordisurf(mod1 ~ Nemouridae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Nemouridae (6%)")

dev.off()


## Test abundances of gammarids
Gammarid_envfit <- envfit(mod1~taxa$Gammaridae, strata=invert_abundance_out$Site_code, perm=999)
Gammarid_vector <- Gammarid_envfit$vectors[1]

## Test abundances of chironomids
Chironomid_envfit <- envfit(mod1~taxa$Chironomidae, strata=invert_abundance_out$Site_code, perm=999)
Chironomid_vector <- Chironomid_envfit$vectors[1]

## Test abundances of caddis
envfit(mod1~taxa$Limnephilidae, strata=invert_abundance_out$Site_code, perm=999)

## Test abundances of stoneflies
envfit(mod1~taxa$Nemouridae, strata=invert_abundance_out$Site_code, perm=999)

## Test litter breakdown
Mass_remaining <- invert_abundance_out$AFDM/3.215
Mass_remaining <- asin(sqrt(Mass_remaining))
Mass_remaining_envfit <- envfit(mod1~Mass_remaining, strata=invert_abundance_out$Site_code, perm=999)
Mass_remaining_vector <- Mass_remaining_envfit$vectors[1]


##############################
##    NMDS plot

png(filename="lpa_alder_coarse_allyrs_NMDS_invert_densities_plot.png", 
    type="cairo",
    units="in", 
    width=7, 
    height=7, 
    pointsize=16, 
    res=600)

head(sites)
sites  <- invert_density_out[,colSums(invert_density_out[c(9:26)])!=0]## Remove missing taxa

max(data.scores$NMDS1)
min(data.scores$NMDS1)

max(data.scores$NMDS2)
min(data.scores$NMDS2)

#dinvert <- decostand(invert_density, "hel")##Abundance data Hellinger-transformed
#dinvert  <- dinvert [,colSums(invert_density)!=0]## Remove missing taxa
#mod1 <- metaMDS(dinvert, dist="euc", trymax = 20, autotransform =F)

data.scores <- as.data.frame(scores(mod1 ))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$number <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$site <- sites[,1]  
data.scores$location <- sites[,2]  
#data.scores$rep <- sites[,6]  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(mod1 , "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

#bact_species_scores <- species.scores

require(plyr)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ site, data = data.scores, FUN = mean)
segs <- merge(data.scores, setNames(cent, c('site','oNMDS1','oNMDS2')),
              by = 'site', sort = FALSE)

find_hull <- function(data.scores) data.scores[chull(data.scores$NMDS1, data.scores$NMDS2), ]
hulls <- ddply(data.scores, "location", find_hull)

p1 <-  ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=location)) + 
  #theme(panel.background = element_rect(fill='white', colour='black')) + 
  #theme(panel.grid.major = none, panel.grid.minor = none) + 
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2), 
               linetype = 2, size = 0.5, colour="gray40", show.legend = NA) + # spiders for sites
  #stat_ellipse(data=data.scores,aes(x=NMDS1,y=NMDS2, colour=location),type = "t",level = 0.5,
  #             show.legend = FALSE) + #, linetype = 2
  geom_polygon(data = hulls, alpha = 0.4) +
  #geom_point(shape=21, size=3, alpha = 0.6) + 
  scale_fill_manual(values=c("grey20","grey60","grey80")) +
  scale_color_manual(values=c("grey20","grey60","grey80"),name = NULL) +
  geom_segment(aes(x = 0, y = 0, xend = Gammarid_vector$arrows[1], yend = Gammarid_vector$arrows[2]), col="grey10") + 
  geom_segment(aes(x = 0, y = 0, xend = Chironomid_vector$arrows[1], yend = Chironomid_vector$arrows[2]), col="grey10") + 
  geom_segment(aes(x = 0, y = 0, xend = Mass_remaining_vector$arrows[1], yend =  Mass_remaining_vector$arrows[2]), col="grey10") + 
  labs(fill = "Sampling \n location") +
  scale_x_continuous(breaks=c(-1.0,-0.5,0,0.5,1.0)) +
  scale_y_continuous(breaks=c(-1.0,-0.5,0,0.5,1.0)) +
  expand_limits(x=c(-1.1,1.1),y=c(-1.1,1.1)) +
  #coord_equal(ratio=1) +
  #coord_fixed() +
  #annotate("text", x = bact_up$NMDS1, y = bact_up$NMDS2, label = "italic(Lacihabitans)",parse=T) +
  #annotate("text", x = bact_down$NMDS1, y = bact_down$NMDS2, label = "italic(Flavobacterium)",parse=T) +
  theme_bw() + 
  theme(legend.position = c(0.1,0.2),
        legend.title = element_text(color = "black", size = 11),
        legend.text = element_text(color = "black", size = 8),
        #legend.direction = "horizontal",
        #axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        #axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=16, vjust = c(0.5), face="bold"), # remove x-axis labels
        axis.title.y = element_text(size=16, vjust = c(0.5), face="bold"), # remove y-axis labels
        axis.text.x = element_text(size=14), # remove x-axis labels
        axis.text.y = element_text(size=14), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

print(p1)

require(grid)
grid.text("Stress = 0.104", x = unit(0.85, "npc"), y = unit(0.15, "npc"),gp=gpar(fontsize=12, col="black"))
grid.text("Gammaridae", x = unit(0.48, "npc"), y = unit(0.14, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text("Chironomidae", x = unit(0.9, "npc"), y = unit(0.475, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text("% litter mass remaining", x = unit(0.52, "npc"), y = unit(0.94, "npc"), gp=gpar(fontsize=12, col="black"))

dev.off()

##************************************************************************************************
## PERMANOVA - test for differences in community composition

head(invert_abundance_out)
Year <- invert_abundance_out$Year
Site_code <- invert_abundance_out$Site_code
Location <- invert_abundance_out$Location

require(pairwiseAdonis)
permanova_2 <- pairwise.adonis2(vegdist(dinvert, "euc") ~ Location + Year, data=invert_abundance_out, strata = 'Site_code')
permanova_2 ## Signficant difference between U1 and U2 could be confounded by dispersion

##************************************************************************************************
## BETAPART - test for differences in community dispersion

## Test effect of Location
Mod_betapart_location <- betadisper(vegdist(dinvert, "euc"), Location, type = c("median"))
Mod_betapart_location ## Signficant difference between U1 and U2 could be confounded by dispersion
anova(Mod_betapart_location)
## S3 method for class 'betadisper'
TukeyHSD(Mod_betapart_location, which = "group", ordered = FALSE, conf.level = 0.95)

## Test effect of Year
Mod_betapart_year <- betadisper(vegdist(dinvert, "euc"), Year, type = c("median"))
Mod_betapart_year ## Signficant difference between U1 and U2 could be confounded by dispersion
anova(Mod_betapart_year)
## S3 method for class 'betadisper'
TukeyHSD(Mod_betapart_year, which = "group", ordered = FALSE, conf.level = 0.95)

png(filename="lpa_alder_coarse_allyrs_invert_density_betadisper.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=12, 
    res=600)

par(mfrow=c(2,2), mar = c(4,4,2,2) + 0.1, pty="s")

## Plot the results from above 
plot(Mod_betapart_location, main = "Location")
boxplot(Mod_betapart_location, ylab = "Distance to centroid",xlab = "Location")

## Plot the results from above 
plot(Mod_betapart_year, main = "Year")
boxplot(Mod_betapart_year, ylab = "Distance to centroid",xlab = "Year")

dev.off()

## Comment: Probably need to use mean values, and missing values could be problematic
## Also consider using abundances instead of densities

############################################################################################
##                                                                                        ##  
##                            LEAFPACK INVERTEBRATES - OAK                                ##  
##                                                                                        ##
############################################################################################

## Rationalise data - need Coarse mesh Oak leaves
lpa_oak_coarse_new <- lpa_taxa[which(lpa_taxa$Mesh=="Coarse" & lpa_taxa$Leaf=="Oak"),]

## Rationalise data - remove NAs (where no AFDM data available)
lpa_oak_coarse_new <- lpa_oak_coarse_new[!is.na(lpa_oak_coarse_new$AFDM),]

## Combine taxa synonymns (Limnephillidae + Trichoptera)
lpa_oak_coarse_new$Limnephillid <- lpa_oak_coarse_new$Limnephillidae+lpa_oak_coarse_new$Trichoptera

## Remove taxa synonyms
lpa_oak_coarse_new <- lpa_oak_coarse_new[,-c(10,17)]

## Check output
head(lpa_oak_coarse_new)

## Re-label headers to ensure taxa families are correct
colnames(lpa_oak_coarse_new) <- c("Year","Site_code","Location","Mesh","Leaf","Rep","ID","AFDM",
                                    "Gammaridae","Baetidae","Nemouridae","Caloptyerigidae",
                                    "Chironomidae","Simuliidae","Asellidae","Rhyacophilidae","Elmidae","Dytiscidae",   
                                    "Tipulidae","Planorbidae","Lymnaeidae", "Oligochaete","Hirudinea","Acari","Other",
                                    "Limnephilidae") 
## Subset and sort columns by Family
taxa <- lpa_oak_coarse_new[,c(9:26)] 
taxa <- taxa[,order(colnames(taxa))]

## Recombine with site information and write output
invert_abundance_out <- data.frame(lpa_oak_coarse_new[,-c(9:26)],taxa)
write.csv(invert_abundance_out,"lpa_oak_coarse_2013_invert_abundances_data.csv",row.names = F)

##***************************************************************************
## Running "FD" package in R
## See: Laliberte, E., Legendre, P., and B. Shipley. (2014). FD: measuring functional
## diversity from multiple traits, and other tools for functional ecology. R package
## version 1.0-12.
##***************************************************************************

## For FRic, FEve, and FDiv see:
## Villeger, S., et al. (2008). New multidimensional functional diversity indices for a 
## multifaceted framework in functional ecology. Ecology 89(8): 2290-2301.

## For FDis see:
## Laliberte, E. and P. Legendre (2010). A distance-based framework for measuring functional 
## diversity from multiple traits." Ecology 91(1): 299-305.

head(invert_abundance_out)

## Remove unidentified taxa "Other"
invert_abund_out <- invert_abundance_out[,-22]

## Remove leafpacks with no observations
invert_abund_out <- invert_abund_out[rowSums(invert_abund_out[-(1:8)]) !=0, ]

## Remove abundances row header for FD function
abundances <- invert_abund_out[,-c(1:8)]

## Remove taxa with no observations
abundances <- abundances[,which(colSums(abundances) !=0)]## Remove missing taxa

head(invert_abund_out)
head(abundances)

## Remove taxa with no observations from trait matrix
Traits <- Traits_mu[-c(which(Traits_mu$Taxa_code=="Acari"),
                       which(Traits_mu$Taxa_code=="Tipulidae")),]

## Remove abundances row header for FD function
rownames(Traits) <- Traits[,1]
traits<-Traits[,-1]

## Calculate FD for all
FD_all <- dbFD(traits, abundances, w.abun = TRUE, stand.x = TRUE, ord = "metric", 
               asym.bin = NULL, corr = "sqrt", 
               calc.FRic = TRUE, m = "max", stand.FRic = T, scale.RaoQ = FALSE, 
               calc.FGR = TRUE, clust.type = "ward", km.inf.gr = 2, km.sup.gr = nrow(x) - 1, 
               km.iter = 100, km.crit = "calinski", calc.CWM = TRUE, CWM.type = "dom", 
               calc.FDiv = TRUE, dist.bin = 2, print.pco = T, messages = TRUE)

#g - separate by group
#9 - number of functional groups

## Extract FD estimates
FD_metrics_all_invert   <- data.frame(FD_all$nbsp,
                                      FD_all$sing.sp,
                                      FD_all$qual.FRic,
                                      FD_all$FRic,
                                      FD_all$FEve,
                                      FD_all$FDiv,
                                      FD_all$FDis,
                                      FD_all$RaoQ,
                                      FD_all$FGR)
## Relabel FD estimates
col_headings <- c('nbsp','sing.sp',"qual.FRic",'FRic','FEve','FDiv','FDis','RaoQ','FGR')
names(FD_metrics_all_invert) <- col_headings

## Write FD estimates to file
write_data <- data.frame(invert_abund_out[,c(1:8)],FD_metrics_all_invert)
write.csv(write_data,file = "lpa_oak_coarse_2013_FD_metrics_raw_data.csv",row.names = F)

## Create mean values for each location
FD_mu <- data.frame(write_data[,-c(4:7)] %>%
                      group_by(Year,Site_code,Location) %>% 
                      summarise_at(vars("AFDM","nbsp","sing.sp","qual.FRic",
                                        "FRic","FEve","FDiv","FDis","RaoQ","FGR"),mean,na.rm=T))
## Convert NA to zero
FD_mu[is.na(FD_mu)] <- 0
write.csv(FD_mu,"lpa_oak_coarse_2013_FD_metrics_mean_data.csv",row.names = F)

#load "pheatmap" library
library(pheatmap)

#calculate the pearson correlation coefficient matrix
myMat.cor <- cor(FD_mu[,-c(1:3,7)], method=c("pearson"))

#plot a heatmap using the calculated correlation matrix
pheatmap(myMat.cor)
dev.off()

## Note potential for postive functional diversity effects
plot(sqrt(AFDM)~log1p(FEve),data=FD_mu)
abline(lm(sqrt(AFDM)~log1p(FEve),data=FD_mu),col="red")
summary(lm(sqrt(AFDM)~log1p(FEve),data=FD_mu))

## Note potential for postive functional diversity effects FDiv looks promising
plot(sqrt(AFDM)~FDiv,data=FD_mu)
abline(lm(sqrt(AFDM)~FDiv,data=FD_mu),col="red")
summary(lm(sqrt(AFDM)~FDiv,data=FD_mu))

##******************************************************************************************
## COMMUNITY-WEIGHTED MEANS
## See: Lavorel, S. and E. Garnier (2002). Predicting changes in community composition and 
## ecosystem functioning from plant traits: revisiting the Holy Grail. 
## Functional Ecology 16(5): 545-556.

## Write CMW (community weighted means) estimates to file
write_data <- data.frame(invert_abund_out[,c(1:8)],FD_all$CWM)
write.csv(write_data,file = "lpa_oak_coarse_2013_CWM_traits_data.csv",row.names = F)
rm(write_data)

## Explore patterns underpinning CMW trait predictors
CWM_traits <- data.frame(invert_abund_out[,c(1:8)],FD_all$CWM)

## Plot feeding traits and other key traits and Gammarid abundances 
png(filename="lpa_oak_coarse_2013_CWM_Gammarid_regressions.png", 
    type="cairo",
    units="in", 
    width=11, 
    height=10, 
    pointsize=16, 
    res=600)

par(mfrow=c(2,2), mar = c(4,4,2,2) + 0.1) 

plot(log(CWM_traits$Food_3)~log1p(invert_abund_out$Gammaridae),
     main = "Food: plant detritus > 1mm",
     ylab = "log Trait abundance",
     xlab = "")
abline(lm(log(CWM_traits$Food_3)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log(CWM_traits$Food_3)~log1p(invert_abund_out$Gammaridae)))
text(4,-1.5,bquote(italic(italic(R^2)[adj])*~"= 0.5327"))

plot(log(CWM_traits$Feeding_3)~log1p(invert_abund_out$Gammaridae),
     main = "Feeding habits: shredder",
     ylab = "log Trait abundance",
     xlab = "")
abline(lm(log(CWM_traits$Feeding_3)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log(CWM_traits$Feeding_3)~log1p(invert_abund_out$Gammaridae)))
text(4,-1.2,bquote(italic(italic(R^2)[adj])*~"= 0.4994"))

plot(log1p(CWM_traits$Size_5)~log1p(invert_abund_out$Gammaridae),
     main = "Maximal size: >2-4 cm",
     ylab = "log Trait abundance",
     xlab = "log+1 Gammarid abundance")
abline(lm(log1p(CWM_traits$Size_5)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log1p(CWM_traits$Size_5)~log1p(invert_abund_out$Gammaridae)))
text(4,0.2,bquote(italic(italic(R^2)[adj])*~"= 0.6775"))

plot(log1p(CWM_traits$Lifestages_4)~log1p(invert_abund_out$Gammaridae),
     main = "Aquatic stages: adult",
     ylab = "log Trait abundance",
     xlab = "log+1 Gammarid abundance")
abline(lm(log1p(CWM_traits$Lifestages_4)~log1p(invert_abund_out$Gammaridae)),col="red")
summary(lm(log1p(CWM_traits$Lifestages_4)~log1p(invert_abund_out$Gammaridae)))
text(4,0.2,bquote(italic(italic(R^2)[adj])*~"= 0.6624"))

dev.off()

## Assess potential for differences between site locations
## Stronger with Oak, but that could also be the Year (2013)

CWM_traits[,c(2,3,57)] %>%
  group_by(Location) %>%
  summarize(mean=mean(Food_3),
            sd=sd(Food_3),
            median=median(Food_3))

CWM_traits[,c(2,3,66)] %>%
  group_by(Location) %>%
  summarize(mean=mean(Feeding_3),
            sd=sd(Feeding_3),
            median=median(Feeding_3))

##***************************************************************************************
## Explore invertebrate data to see how dominant gammarids are and who the other players
## are
##***************************************************************************************

## Rationalise data - need Coarse mesh Oak leaves
lpa_oak_coarse_new <- lpa_taxa[which(lpa_taxa$Mesh=="Coarse" & lpa_taxa$Leaf=="Oak"),]

## Rationalise data - remove NAs (where no AFDM data available)
lpa_oak_coarse_new <- lpa_oak_coarse_new[!is.na(lpa_oak_coarse_new$AFDM),]

## Combine taxa synonymns (Limnephillidae + Trichoptera)
lpa_oak_coarse_new$Limnephillid <- lpa_oak_coarse_new$Limnephillidae+lpa_oak_coarse_new$Trichoptera

## Remove taxa synonyms
lpa_oak_coarse_new <- lpa_oak_coarse_new[,-c(10,17)]

## Check output
head(lpa_oak_coarse_new)

## Re-label headers to ensure taxa families are correct
colnames(lpa_oak_coarse_new) <- c("Year","Site_code","Location","Mesh","Leaf","Rep","ID","AFDM",
                                    "Gammaridae","Baetidae","Nemouridae","Caloptyerigidae",
                                    "Chironomidae","Simuliidae","Asellidae","Rhyacophilidae","Elmidae","Dytiscidae",   
                                    "Tipulidae","Planorbidae","Lymnaeidae", "Oligochaete","Hirudinea","Acari","Other",
                                    "Limnephilidae") 
## Subset and sort columns by Family
taxa <- lpa_oak_coarse_new[,c(9:26)] 
taxa <- taxa[,order(colnames(taxa))]

## Recombine with site information and write output
invert_abundance_out <- data.frame(lpa_oak_coarse_new[,-c(9:26)],taxa)
row.names(invert_abundance_out) <- NULL 

## 2013 MES US1 Pravin rep 5 Incorrect initial leaf mass or human error - way too low

invert_abundance_out <- invert_abundance_out[-c(131),]

taxa <- invert_abundance_out[-184,c(9:26)]

##*****************************************************************************
## ABUNDANCES and proportions - to see the main players   
## Sum column abundances
Total_inverts <-  colSums(taxa, na.rm = FALSE, dims = 1)

head(Total_inverts)

## Calculate total proportions for dominant Families
Chironomid_abundance    <- Total_inverts[5]/sum(Total_inverts)*100
Gammarid_abundance    <- Total_inverts[8]/sum(Total_inverts)*100
Limnephilid_abundance <- Total_inverts[10]/sum(Total_inverts)*100
Nemourid_abundance    <- Total_inverts[12]/sum(Total_inverts)*100

##*****************************************************************************
## DENSITIES and relative abundance
## Convert abundances to densities
invert_density <- taxa/invert_abundance_out$AFDM[-184]

## Sum column densities
Total_invert_densities <-  colSums (invert_density, na.rm = FALSE, dims = 1)

## Calculate total proportions for dominant Families
Chironomid_density    <- Total_invert_densities[5]/sum(Total_invert_densities)*100
Gammarid_density    <- Total_invert_densities[8]/sum(Total_invert_densities)*100
Limnephilid_density <- Total_invert_densities[10]/sum(Total_invert_densities)*100
Nemourid_density    <- Total_invert_densities[12]/sum(Total_invert_densities)*100

## Write file for invertebrate densities
## Recombine with site information and write output
invert_density_out <- data.frame(invert_abundance_out[-184,-c(9:26)],invert_density)

write.csv(invert_density_out,"lpa_oak_coarse_2013_densities_data.csv",row.names = F)

############################################################################################
##                                                                                        ##  
##                     OAK LEAFPACK INVERTEBRATES COMMUNITIES                             ##  
##                                                                                        ##
############################################################################################

## Load vegan for multivariate analysis
require(vegan)

## Generate data
invert_density <- taxa/invert_abundance_out$AFDM[-184]

## remove taxa missing representation
invert_density <- invert_density[,colSums(invert_density)!=0]

## Remove leafpacks with no observations
invert_density <- invert_density[rowSums(invert_density) !=0, ]

## Replace NAs with zeros (if any)
invert_density[is.na(invert_density)] <- 0

## Hellinger transform density data and run NMDS
dinvert <- decostand(invert_density, "hell")
mod1    <- metaMDS(dinvert , dist="euc", trymax = 20, autotransform = F)

## Extracts scores for plotting below
nmsc <- scores(mod1, display = "sites")
species <- scores(mod1, display = "species")

#write.csv(data.frame(lpa_oak_coarse[,-c(7:23)],nmsc),"NMDS_sites_hell_200802.csv",row.names = F)
#write.csv(lpa_oak_coarse,"Coarse_oak_invertebrates_200802.csv",row.names = F)

## Compare gammarid densities and NMDS2
gammarid_density <- data.frame(invert_density[,7],nmsc)
names(gammarid_density)[1] <- "Gammarid_density"
head(gammarid_density)

## remove taxa missing representation
taxa <- taxa[,colSums(taxa)!=0]

## Remove leafpacks with no observations
taxa <- taxa[rowSums(taxa) !=0, ]

## Replace NAs with zeros (if any)
taxa[is.na(taxa)] <- 0

## Compare gammarid abundances and NMDS2
gammarid_abund <- data.frame(taxa[,7],nmsc)
names(gammarid_abund)[1] <- "Gammarid_abundance"
head(gammarid_abund)

## Compare gammarid abundances/densities and NMDS1
png(filename="lpa_oak_coarse_2013_NMDS1_gammarid_regressions.png", 
    type="cairo",
    units="in", 
    width=9, 
    height=6, 
    pointsize=13, 
    res=600)

par(mfrow=c(1,2), mar = c(4,4,2,2) + 0.1, pty="s")

plot(NMDS1 ~ log1p(Gammarid_abundance), gammarid_abund, 
     main = "Gammaridae abundance",
     xlab = "log+1 Gammarid abundance")
abline(lm(NMDS1 ~ log1p(Gammarid_abundance), gammarid_abund),col="red")
summary(lm(NMDS1 ~ log1p(Gammarid_abundance), gammarid_abund))
text(4,0.4,bquote(italic(italic(R^2)[adj])*~"= 0.5273"))

plot(NMDS1 ~ log1p(Gammarid_density), gammarid_density, 
     main = "Gammaridae density",
     xlab = "log+1 Gammarid density")
abline(lm(NMDS1 ~ log1p(Gammarid_density), gammarid_density),col="red")
summary(lm(NMDS1 ~ log1p(Gammarid_density), gammarid_density))
text(3.5,0.4,bquote(italic(italic(R^2)[adj])*~"= 0.4993"))

dev.off()

## Plot contour plot to better undestand the influence of dominant taxa
png(filename="lpa_oak_coarse_2013_invert_density_NMDS_contour_plot.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=12, 
    res=600)

par(mfrow=c(2,2), mar = c(4,4,2,2) + 0.1, pty="s")

ordisurf(mod1 ~ Chironomidae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Chironomidae (49%)")
ordisurf(mod1 ~ Gammaridae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Gammaridae (37%)")
ordisurf(mod1 ~ Limnephilidae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Limnephilidae (3%)")
ordisurf(mod1 ~ Nemouridae, dinvert, bubble =4, method = "REML", bs = "ts",
         main="Nemouridae (4%)")

dev.off()

## Create file for sites
sites <- invert_density_out[,c(1:26)]

## remove taxa missing representation
sites<- sites[,colSums(sites[,c(9:26)])!=0]

## Remove leafpacks with no observations
sites <- sites[rowSums(sites[,c(9:22)]) !=0, ]

## Replace NAs with zeros (if any)
sites[is.na(sites)] <- 0

## Test abundances of gammarids
Gammarid_envfit <- envfit(mod1~taxa$Gammaridae, strata=sites$Site_code, perm=999)
Gammarid_vector <- Gammarid_envfit$vectors[1]

## Test abundances of chironomids
Chironomid_envfit <- envfit(mod1~taxa$Chironomidae, strata=sites$Site_code, perm=999)
Chironomid_vector <- Chironomid_envfit$vectors[1]

## Test abundances of caddis
Limnephilid_envfit <- envfit(mod1~taxa$Limnephilidae, strata=sites$Site_code, perm=999)
Limnephilid_vector <- Limnephilid_envfit$vectors[1]

## Test abundances of stoneflies
Nemourid_envfit <- envfit(mod1~taxa$Nemouridae, strata=sites$Site_code, perm=999)
Nemourid_vector <- Nemourid_envfit$vectors[1]

## Generate data
taxa <- invert_abundance_out[-184,c(9:26)]
invert_density <- taxa/invert_abundance_out$AFDM[-184]

## Create data frame
invert_density_out <- data.frame(invert_abundance_out[-184,-c(9:26)],invert_density)
head(invert_density_out)

data <- invert_density_out[,c(8:26)]

## remove taxa missing representation
data <- data[,colSums(data)!=0]

## Remove leafpacks with no observations
data <- data[rowSums(data[,-1]) !=0, ]

## Replace NAs with zeros (if any)
data[is.na(data)] <- 0

## Test litter breakdown
Mass_remaining <- data$AFDM/3.520
Mass_remaining <- asin(sqrt(Mass_remaining))
Mass_remaining_envfit <- envfit(mod1~Mass_remaining, strata=sites$Site_code, perm=999)
Mass_remaining_vector <- Mass_remaining_envfit$vectors[1]

##############################
##    NMDS plot

png(filename="lpa_oak_coarse_allyrs_NMDS_invert_densities_plot.png", 
    type="cairo",
    units="in", 
    width=7, 
    height=7, 
    pointsize=16, 
    res=600)

sites <- invert_density_out[,c(1:26)]

## remove taxa missing representation
sites<- sites[,colSums(sites[,c(9:26)])!=0]

## Remove leafpacks with no observations
sites <- sites[rowSums(sites[,c(9:22)]) !=0, ]

## Replace NAs with zeros (if any)
sites[is.na(sites)] <- 0

max(data.scores$NMDS1)
min(data.scores$NMDS1)

max(data.scores$NMDS2)
min(data.scores$NMDS2)

#dinvert <- decostand(invert_density, "hel")##Abundance data Hellinger-transformed
#dinvert  <- dinvert [,colSums(invert_density)!=0]## Remove missing taxa
#mod1 <- metaMDS(dinvert, dist="euc", trymax = 20, autotransform =F)

data.scores <- as.data.frame(scores(mod1 ))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$number <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
data.scores$site <- sites[,1]  
data.scores$location <- sites[,2]  
#data.scores$rep <- sites[,6]  #  add the grp variable created earlier
head(data.scores)  #look at the data

species.scores <- as.data.frame(scores(mod1 , "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
head(species.scores)  #look at the data

#bact_species_scores <- species.scores

require(plyr)
cent <- aggregate(cbind(NMDS1, NMDS2) ~ site, data = data.scores, FUN = mean)
segs <- merge(data.scores, setNames(cent, c('site','oNMDS1','oNMDS2')),
              by = 'site', sort = FALSE)

find_hull <- function(data.scores) data.scores[chull(data.scores$NMDS1, data.scores$NMDS2), ]
hulls <- ddply(data.scores, "location", find_hull)

p1 <-  ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=location)) + 
  #theme(panel.background = element_rect(fill='white', colour='black')) + 
  #theme(panel.grid.major = none, panel.grid.minor = none) + 
  geom_segment(data = segs, mapping = aes(xend = oNMDS1, yend = oNMDS2), 
               linetype = 2, size = 0.5, colour="gray40", show.legend = NA) + # spiders for sites
  #stat_ellipse(data=data.scores,aes(x=NMDS1,y=NMDS2, colour=location),type = "t",level = 0.5,
  #             show.legend = FALSE) + #, linetype = 2
  geom_polygon(data = hulls, alpha = 0.4) +
  #geom_point(shape=21, size=3, alpha = 0.6) + 
  scale_fill_manual(values=c("grey20","grey60","grey80")) +
  scale_color_manual(values=c("grey20","grey60","grey80"),name = NULL) +
  geom_segment(aes(x = 0, y = 0, xend = Gammarid_vector$arrows[1], yend = Gammarid_vector$arrows[2]), col="grey10") + 
  geom_segment(aes(x = 0, y = 0, xend = Chironomid_vector$arrows[1], yend = Chironomid_vector$arrows[2]), col="grey10") + 
  geom_segment(aes(x = 0, y = 0, xend = Nemourid_vector$arrows[1], yend = Nemourid_vector$arrows[2]), col="grey10") + 
  geom_segment(aes(x = 0, y = 0, xend = Limnephilid_vector$arrows[1], yend = Limnephilid_vector$arrows[2]), col="grey10") + 
  geom_segment(aes(x = 0, y = 0, xend = Mass_remaining_vector$arrows[1], yend =  Mass_remaining_vector$arrows[2]), col="grey10") + 
  labs(fill = "Sampling \n location") +
  scale_x_continuous(breaks=c(-0.8,-0.4,0,0.4,0.8)) +
  scale_y_continuous(breaks=c(-0.8,-0.4,0,0.4,0.8)) +
  expand_limits(x=c(-0.9,1.1),y=c(-0.9,1.1)) +
  #coord_equal(ratio=1) +
  #coord_fixed() +
  #annotate("text", x = bact_up$NMDS1, y = bact_up$NMDS2, label = "italic(Lacihabitans)",parse=T) +
  #annotate("text", x = bact_down$NMDS1, y = bact_down$NMDS2, label = "italic(Flavobacterium)",parse=T) +
  theme_bw() + 
  theme(legend.position = c(0.1,0.85),
        legend.title = element_text(color = "black", size = 11),
        legend.text = element_text(color = "black", size = 8),
        #legend.direction = "horizontal",
        #axis.text.x = element_blank(),  # remove x-axis text
        #axis.text.y = element_blank(), # remove y-axis text
        #axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=16, vjust = c(0.5), face="bold"), # remove x-axis labels
        axis.title.y = element_text(size=16, vjust = c(0.5), face="bold"), # remove y-axis labels
        axis.text.x = element_text(size=14), # remove x-axis labels
        axis.text.y = element_text(size=14), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

print(p1)

require(grid)
grid.text("Stress = 0.123", x = unit(0.85, "npc"), y = unit(0.16, "npc"),gp=gpar(fontsize=12, col="black"))
grid.text("Gammaridae", x = unit(0.19, "npc"), y = unit(0.21, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text("Chironomidae", x = unit(0.84, "npc"), y = unit(0.26, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text("Limnephilidae", x = unit(0.6, "npc"), y = unit(0.92, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text("Nemouridae", x = unit(0.45, "npc"), y = unit(0.94, "npc"), gp=gpar(fontsize=12, col="black"))
grid.text("% litter mass\n remaining", x = unit(0.87, "npc"), y = unit(0.72, "npc"), gp=gpar(fontsize=12, col="black"))

dev.off()

Gammarid_vector
Chironomid_vector
Limnephilid_vector
Nemourid_vector
Mass_remaining_vector

##************************************************************************************************
## PERMANOVA - test for differences in community composition

head(sites)
Year <- sites$Year
Site_code <- sites$Site_code
Location <- sites$Location

require(pairwiseAdonis)
permanova_2 <- pairwise.adonis2(vegdist(dinvert, "euc") ~ Location, data=sites, strata = 'Site_code')
permanova_2 ## Signficant difference between U1 and U2 could be confounded by dispersion

##************************************************************************************************
## BETAPART - test for differences in community dispersion

## Test effect of Location
Mod_betapart_location <- betadisper(vegdist(dinvert, "euc"), Location, type = c("median"))
Mod_betapart_location ## Signficant difference between U1 and U2 could be confounded by dispersion
anova(Mod_betapart_location)
## S3 method for class 'betadisper'
TukeyHSD(Mod_betapart_location, which = "group", ordered = FALSE, conf.level = 0.95)

png(filename="lpa_oak_coarse_2013_invert_density_betadisper.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=12, 
    res=600)

par(mfrow=c(1,2), mar = c(4,4,2,2) + 0.1, pty="s")

## Plot the results from above 
plot(Mod_betapart_location, main = "Location")
boxplot(Mod_betapart_location, ylab = "Distance to centroid",xlab = "Location")

dev.off()
