#############################################################################################
##                                                                                         ##  
##                               SANKEY FIGURE FOR MPS                                     ##   
##                                                                                         ##  
#############################################################################################

##  Script: ECOIMPACT leafpack assay data analysis
##  Author: Dr. Francis J. Burdon
##  Objectives:
##  1. Use Sankey diagram to describe different MPs classes due to Agriculture (US) and
##     amd WW

##***************************************************************************************
## Load libraries you need

##Loads package used for Logit transformation (% data)
require(car)

##Loads vegan package for transformation
require(vegan)

##Loads reshape and plyr package for data manipulation
require(reshape)
require(reshape2)
require(plyr)

## Load packages for Sankey 
require(networkD3)
require(dplyr)
require(tidyverse)
require(data.table)
require(janitor)

## Load package to convert pdf to png
require(pdftools)

## Load data files
require(readr)

#############################################################################################
##                                                                                         ##  
##                                  LOAD DATA                                              ##
##                                                                                         ##
#############################################################################################

## Environmental
MPs <- read.csv(here::here("1_Input_data","3_DATA_MPs_TU_metasubstance_VER02_161114.csv")) 

#############################################################################################
##                                                                                         ##  
##                                  ENVIRONMENTAL                                          ##
##                                                                                         ##  
#############################################################################################

## Arrange files (but depracated, because I will use the merge function)
MPs<-arrange(MPs, Year, Site, Location)

## MICROPOLLUTANTS

## Calculate non-Insecticide TUs
MPs <- as.data.frame(MPs)
MPs$MP_nonInsecticide <- NA
MPs$MP_nonInsecticide <- MPs[,4]-MPs[,18]

head(MPs[,c(1:6,16:18,20,21)])
## Change Data_all

## Merge with data and rename columns
Data_all <-MPs[,c(1:6,16:18,20,21)]
colnames(Data_all) <- c("Year","Site","Location",
                        "Total_TUs","NonPest_TUs","Pesticide_TUs","Fungicide_TUs", 
                        "Herbicide_TUs","Insecticide_TUs",
                        "Pharmaceutical_TUs","NonInsecticide_TUs")

## Visual plots to check data makes sense
plot(log(NonPest_TUs)~log(Total_TUs),Data_all)
abline(0,1,col="blue")
plot(log(Pesticide_TUs)~log(Total_TUs),Data_all)
abline(0,1,col="blue")
plot(log(Fungicide_TUs)~log(Pesticide_TUs),Data_all)
abline(0,1,col="blue")
plot(log(Herbicide_TUs)~log(Pesticide_TUs),Data_all)
abline(0,1,col="blue")
plot(log(Insecticide_TUs)~log(Pesticide_TUs),Data_all)
abline(0,1,col="blue")## See step function - must be WW
plot(log(Pharmaceutical_TUs)~log(Total_TUs),Data_all)
abline(0,1,col="blue")
plot(log(NonInsecticide_TUs)~log(Total_TUs),Data_all)
abline(0,1,col="blue")

#############################################################################################
##                                                                                         ##  
##                       CREATE MEAN VALUES FOR SANKEY DIAGRAM                             ##
##                                                                                         ##  
#############################################################################################

## Calculate mean values
MPs_mu <- data.frame(Data_all %>%
                      group_by(Location) %>% 
                      summarise_at(vars("Fungicide_TUs","Herbicide_TUs","Insecticide_TUs",
                                        "Pharmaceutical_TUs"),mean,na.rm=T))
## Transpose data
MPs_mu <- data.frame(t(as.matrix(MPs_mu )))
## Convert row names to column
MPs_mu <- setDT(MPs_mu, keep.rownames = TRUE)[]
## Convert first row to column header
MPs_mu <- MPs_mu %>% row_to_names(row_number = 1)
## Drop U2 to focus on U1
MPs_mu <- MPs_mu[,-4]
## Rorder rows
row.names(MPs_mu) <- NULL
## Relabel columns
colnames(MPs_mu) <- c("Location","D","Agriculture")
## Convert to data frame
MPs_mu <- as.data.frame(MPs_mu)
## Convert data to numeric form
MPs_mu[,2] <- as.numeric(as.character(MPs_mu$D))
MPs_mu[,3] <- as.numeric(as.character(MPs_mu$Agriculture))
## Create estimated influence of Agriculture and WW
Agriculture <- MPs_mu[,3]
Wastewater <- MPs_mu[,2]-MPs_mu[,3]
## Create vectors               
Ag <- data.frame(MPs_mu[,1],Agriculture)
WW <- data.frame(MPs_mu[,1],Wastewater)
## Relable 
colnames(Ag) <- c("guide","value")
colnames(WW) <- c("guide","value")
## Relable 
data <- rbind(Ag,WW)

#############################################################################################
##                                                                                         ##  
##                                  SANKEY CHART                                           ##
##                                                                                         ##
#############################################################################################

## Create values from the estimates calculated above
values <- data$value
values <- log1p(values) ## Log scale since Insecticides dominate TUs

# A connection data frame is a list of flows with intensity for each flow
links <- data.frame(
  source=c("Agriculture","Agriculture","Agriculture","Agriculture","Wastewater","Wastewater","Wastewater","Wastewater"), 
  target=c("Fungicides","Herbicides","Insecticides","Pharmaceuticals"), 
  value=values
)

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, units = "TUs", fontSize = 20,
                   fontFamily = "Arial", nodeWidth = 30)
p
## See viewer - export to Chrome and save as .pdf

# save the widget
# library(htmlwidgets)
# saveWidget(p, file=paste0( getwd(), "sankeyBasic1.html"))

##*************************************************
## Convert .pdf to .png for further use
pdf_convert("Sankey_TUs.pdf",
                        format = "png",
                        pages = 1,
                        filenames = NULL,
                        dpi = 1000,
                        antialias = TRUE,
                        opw = "",
                        upw = "",
                        verbose = TRUE)
