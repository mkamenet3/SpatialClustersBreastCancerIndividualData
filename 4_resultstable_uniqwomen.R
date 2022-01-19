#M.Kamenetsky
#Updated 2021-5-22
#Table of unique RR for unique women controls
#############################
#load required packages
#############################
rm(list=ls())

library(spdep)
library(ggplot2)
library(rgdal)
library(rgeos)
library(cowplot)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(sf)
library(raster)
library(gt)

table.out <- NULL

'%ni%' <- Negate('%in%')
make_circles <- function(centers, radius, nPoints = 100){
    # centers: the data frame of centers with ID
    # radius: radius measured in kilometer
    #
    meanLat <- mean(centers$latcentroids)
    # length per longitude changes with lattitude, so need correction
    radiusLon <- radius /111 / cos(meanLat/57.3)
    radiusLat <- radius / 111
    circleDF <- data.frame(ID = rep(centers$centerID, each = nPoints))
    angle <- seq(0,2*pi,length.out = nPoints)
    
    circleDF$lon <- unlist(lapply(centers$loncentroids, function(x) x + radiusLon * cos(angle)))
    circleDF$lat <- unlist(lapply(centers$latcentroids, function(x) x + radiusLat * sin(angle)))
    return(circleDF)
}


############################################################################################
#load results
############################################################################################
miresults <- read.csv("../../Results/Analysis/Wisconsin/tables/dat_pooledmiresults_20210517.csv")
clustersdf_wisc_bic <- read.csv("../../Results/Analysis/Wisconsin/tables/clustersdf_wisc_bic.csv")%>%
    mutate(centerID = oldclusterid) 


################################################################################################################################
################################################################################################################################
################################################################################################################################
#PART 3: Make Table
################################################################################################################################
################################################################################################################################
################################################################################################################################


###############################################
#Fully-Adjusted
###############################################

geocodes <- clustersdf_wisc_bic %>%
    mutate(clusterid = oldclusterid) %>%
    dplyr::select(clusterid,loncentroids, latcentroids, r) %>%
    mutate(City = if_else(clusterid==164, "Lake Superior",
                          if_else(clusterid==12378, "Weston",
                                  if_else(clusterid==283, "Platteville",
                                          if_else(clusterid==91588, "Greenfield",
                                                  if_else(clusterid==64, "Alden",
                                                          if_else(clusterid==4375, "Hudson", 
                                                                  if_else(clusterid==100, "Platteville",
                                                                          if_else(clusterid==850, "Smelser",
                                                                                  if_else(clusterid==72, "Lake Superior",
                                                                                          if_else(clusterid==106, "Lake Superior",
                                                                                                  if_else(clusterid==1450, "Superior",
                                                                                                          if_else(clusterid==1541, "Superior",
                                                                                                                  if_else(clusterid==1625, "Superior",
                                                                                                                          if_else(clusterid==68729, "West Allis",
                                                                                                                                  if_else(clusterid==72637, "West Allis","NA")))))))))))))))) %>%
    mutate(sigclusters = ifelse(clusterid %in% c(12378, 91588, 64, 4375),1,0),
           clusterid2= paste0(City, "_", clusterid),
           clusterid_orig = paste0("cluster", clusterid))

# #####################################################
#Because these are fitted values, a woman can belong to multiple clusters
sumcasesctrl <- miresults %>%
    dplyr::select(cluster100:cluster4375,miest_age,miSEest_age,
                  miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                  miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,bc) %>%
    pivot_longer(cluster100:cluster4375, names_to="clusterid", values_to = "count") %>%
    filter(count==1)%>%
    group_by(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)  %>% mutate(uniqRRID = cur_group_id()) %>% ungroup()%>%
    group_by(uniqRRID) %>%
    summarize(sumcases = sum(bc),
              sumcontrols = sum(ifelse(bc==1,0,1)))
    
    
fullyadjusted_uniq <- miresults %>%

    dplyr::select(cluster100:cluster4375,miest_age,miSEest_age,
                  miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                  miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,bc) %>%
    pivot_longer(cluster100:cluster4375, names_to="clusterid", values_to = "count") %>%
    filter(count==1)%>%
    group_by(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)  %>% mutate(uniqRRID = cur_group_id()) %>% ungroup()%>%
    
    dplyr::filter(bc==0) %>%
    dplyr::select(-bc) %>%

    
    distinct(miest_age, miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks, .keep_all = TRUE) %>%
    mutate(modfull = miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
           LBfull = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)-1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks),
           UBfull = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)+1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    mutate(modage = miest_age,
           LBage = exp(log(miest_age)-1.96*miSEest_age),
           UBage = exp(log(miest_age)+1.96*miSEest_age))

fullyadjusted_uniq_casesctrl <- merge(fullyadjusted_uniq, sumcasesctrl, by="uniqRRID")    

fullyadjusted_uniqout <- merge(fullyadjusted_uniq_casesctrl, geocodes, by.x="clusterid", by.y="clusterid_orig") %>%
    mutate(countynum = if_else(City=="Superior",1,
                               if_else(City=="Lake Superior",2,
                                       if_else(City %in% c("Alden", "Hudson"),3,
                                               if_else(City %in% c("Platteville","Smelser"),4,
                                                       if_else(City=="Weston",5,
                                                               if_else(City %in% c("West Allis","Greenfield"),6, 0))))))) %>%
    mutate(countyID = if_else(countynum==1,"Douglas County",
                              if_else(countynum==2,"Ashland & Bayfield Counties",
                                      if_else(countynum==3,"Polk & St. Croix Counties",
                                              if_else(countynum==4,"Grant & Lafayette Counties",
                                                      if_else(countynum==5,"Marathon County",
                                                              if_else(countynum==6,"Waukesha & Milwaukee Counties", NA_character_))))))) %>%
    mutate_at(c("modage", "LBage", "UBage", "modfull", "LBfull", "UBfull"), function(x) round(x,2)) %>%
    arrange(countynum, modage) %>%
    mutate(outputage = paste0("(", LBage, ", ",UBage, ")"),
           outputfull = paste0("(", LBfull, ", ",UBfull, ")")) %>%
    relocate(City, modage, outputage, modfull, outputfull) %>%
    dplyr::select(countyID,City, modage, outputage, modfull, outputfull, sumcases, sumcontrols)

###########################################
#MAKE TABLE
###########################################
fullyadjusted_uniqout %>%
    gt(groupname_col = "countyID") %>%
    cols_label(City= md("**City**"),
               modage = md("**OR**"),
               outputage = md("**(95% CI)**"),
               modfull= md("**OR**"),
               outputfull = md("**(95% CI)**"),
               sumcases = md("**Cases**"),
               sumcontrols = md("**Controls**")) %>%
        tab_spanner(label=md("**Age-Adjusted**"), columns=vars(modage, outputage)) %>%
        tab_spanner(label=md("**Fully-Adjusted**"), columns=vars(modfull, outputfull)) %>%
        tab_header(title="Fitted Results from Participant-Level Spatial Cluster Model, Wisconsin, 1997-2000") %>%
    gtsave("../../manuscript/tables/modelresultstable_fitted.tex")





