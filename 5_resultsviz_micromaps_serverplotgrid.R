#M.Kamenetsky
#Updated 2020-06-09
#Analyse cluster results and create logistic model + viz based on BIC cluster selection
#Update 2020-6-27 plot representative cluster ppl instead using micromaps
#This script run on the server because the maps take too much memory
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
library(xtable)
library(geosphere)
library(forcats)
library(mapproj)
library(maptools)
library(grid)
library(gridExtra)
library(ggsn)
library(ggspatial)

table.out <- NULL
############################################################################################
#Helper Functions
############################################################################################
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
#Load results
############################################################################################
miresults <- read.csv("dat_pooledmiresults_20210517.csv")
clustersdf_wisc_bic <- read.csv("data/clustersdf_wisc_bic.csv")%>%
  mutate(centerID = oldclusterid) 


######################################################################################
#######################################################################################
#######################################################################################
#PART 1: Identify ALL Clusters for Plotting
#######################################################################################
#######################################################################################
#######################################################################################

################################################################################################################################
################################################################################################################################
################################################################################################################################
#PART 3: PLOTS
################################################################################################################################
################################################################################################################################
################################################################################################################################
######################################
#PlotRaw Case and Control Points
######################################
sigclusters <- c(12378, 91588, 64, 4375)
#Load maps
cases <- miresults %>% filter(bc==1) 
controls <- miresults %>% filter(bc==0)
undercounties <- c("Barron", "Bayfield","Buffalo","Burnett",
                   "Dunn","Eau Claire", "Pepin","Pierce","Polk", "St. Croix")
wi <- readOGR(dsn=".", layer="wi", verbose = TRUE)


wi_map <- fortify(wi)
wi_undercount <- subset(wi, wi$NAME %in% undercounties)
undercounties_map <- fortify(wi_undercount)

p1 <- ggplot() + geom_map(data=wi_map, map=wi_map, 
                          aes(x=long, y=lat, map_id=id),
                          color="grey", fill="white", size=0.25) +
    coord_map() +
    theme_map() 

p2 <- p1 + geom_map(data=undercounties_map, map=undercounties_map,
                    aes(x=long, y=lat, map_id=id),
                    color="black",fill="white", size=1)



###############################################
#Fully-Adjusted
###############################################

set.seed(2)
representativeclusterppl_ctrls <- miresults %>%
    filter(bc==0) %>%
    group_by(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    filter(row_number()==sample(1:n(),1)) %>%
    mutate(expxbetahat = miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
           LB = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)-1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks),
           UB = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)+1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    dplyr::select(expxbetahat, LB, UB, lat, lon) %>%
    distinct(expxbetahat, LB, UB, .keep_all = TRUE) 

ctrls <- miresults %>%
    filter(bc==0) %>%
    group_by(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    mutate(expxbetahat = miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
           LB = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)-1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks),
           UB = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)+1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    dplyr::select(expxbetahat, LB, UB, lat, lon) 

representativeclusterppl_ctrls_sf <- st_as_sf(as.data.frame(representativeclusterppl_ctrls), 
                                              coords=c("lon", "lat"), 
                                              crs="+proj=longlat + ellps=WGS84")

ctrls_sf <- st_as_sf(as.data.frame(ctrls), 
                     coords=c("lon", "lat"), 
                     crs="+proj=longlat + ellps=WGS84")

out <- ctrls
out <- out %>% arrange(-expxbetahat)


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
  


#####################################################
########################################
#Over all cluster OR
fullyadjusted <- miresults %>%
    dplyr::filter(bc==0) %>%
    mutate(inacluster= ifelse(inacluster==1,0,1)) %>%
    pivot_longer(cluster100:inacluster, names_to="clusterid", values_to = "count") %>%
    dplyr::filter(count==1) %>%
    dplyr::mutate(clusterid=recode(clusterid, inacluster="No Cluster", .default=levels(clusterid)))
fullyadjusted<- merge(fullyadjusted, geocodes, by.x="clusterid", by.y="clusterid_orig", all.x=TRUE)


#merge in geocodes with cities

fullyadjusted_all<- fullyadjusted %>% 
    mutate(expxbetahat = miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
           LB = exp(log(miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)-1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks),
           UB = exp(log(miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)+1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    mutate(womanid = paste0(City," ",X)) %>%
    arrange(expxbetahat)


fullyadjusted_uniq<- fullyadjusted %>%
    mutate(expxbetahat = miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks) %>%
    distinct(expxbetahat,  .keep_all = TRUE) %>%
    mutate(LB = exp(log(expxbetahat) -1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks),
           UB = exp(log(expxbetahat) +1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)) %>%
    mutate(UB = ifelse(UB >3,3,UB)) %>%
    arrange(-expxbetahat) %>%
    mutate(womanid = paste0(City, " ", X))




################################################################
#By Groups of 6
################################################################
makemicropanels <- function(fullyadjusted_all,fullyadjusted_unique,clustersdf_wisc_bic, pcounties,wi_undercount, wi, colorsii, map){
  #browser()
    ##################################
    #Plot panel counties
    ##################################
    wi_counties <- subset(wi, wi$NAME %in% pcounties)
    wi_map <- fortify(wi)
    counties_map <- fortify(wi_counties)
    
    clusters_bic_counties <- subset(clustersdf_wisc_bic, nb %in% pcounties)
    clusters_bic_c <- lapply(1:nrow(clusters_bic_counties ), function(x) make_circles(clusters_bic_counties[x,c(13,14,3)], radius=clusters_bic_counties [x,6]))
    
    ##################################
    #Plot (if any) under-reporting counties
    ##################################
    p1 <- ggplot() + 
        geom_map(data=counties_map, map=counties_map,
                 aes(x=long, y=lat, map_id=id),
                 color="grey",fill="white", size=1) +
        #coord_equal() +
        theme_bw()
        pcounties_map <- ggplot() + 
            geom_map(data=counties_map, map=counties_map,
                     aes(x=long, y=lat, map_id=id),
                     color="grey",fill="white", size=1) +
            coord_equal() 
            theme_bw() 

    ###################################################################
    #Create Counties Map (counties identified on WI map by color block)
    wi_zoom<- subset(wi, wi$NAME %in% pcounties)
    zoomcounties_map <- fortify(wi_zoom)
   # browser()
    if(pcounties=="Douglas"){
      countymap <- ggplot() + 
        geom_map(data=wi_map, map=wi_map,
                 aes(x=long, y=lat, map_id=id),
                 color="grey",fill="white", size=1) +
        geom_map(data=zoomcounties_map, map=zoomcounties_map,
                 aes(x=long, y=lat, map_id=id),fill="grey60", size=1, color="grey80")+
        geom_map(data=undercounties_map, map=undercounties_map,
                 aes(x=long, y=lat, map_id=id),
                 color="black",fill=NA, size=1)  +
        coord_fixed(1.5)+
        theme_map() +
        ggsn::scalebar(wi_map, dist=50, location="topright", st.size=5, transform=TRUE, dist_unit="km",
                       st.dist = 0.1)
    } else {
      countymap <- ggplot() + 
        geom_map(data=wi_map, map=wi_map,
                 aes(x=long, y=lat, map_id=id),
                 color="grey",fill="white", size=1) +
        geom_map(data=zoomcounties_map, map=zoomcounties_map,
                 aes(x=long, y=lat, map_id=id),fill="grey60", size=1, color="grey80")+
        geom_map(data=undercounties_map, map=undercounties_map,
                 aes(x=long, y=lat, map_id=id),
                 color="black",fill=NA, size=1)  +
        coord_fixed(1.5)+
        theme_map()
      
    }
 
    ###################################################################
    ##################################
    #Extract women in selected polygons
    ##################################
    fullyadjusted_all <- merge(fullyadjusted_all, subset(fullyadjusted_uniq, select=c(womanid, clusterid2)),
                               by=c("womanid"), all.x=TRUE)
    #need to fill in by expxbetahat
    fullyadjusted_all<- fullyadjusted_all %>%
        group_by(expxbetahat) %>%
        arrange(expxbetahat) %>%
        fill(clusterid2.y, .direction = "updown")
    ##ALL
    coords <- cbind(fullyadjusted_all$lon, fullyadjusted_all$lat)
    #this next step creates the spatialpointsdataframe
    fullyadjusted_all_spdf <- SpatialPointsDataFrame(coords = coords,
                                                     data=fullyadjusted_all,
                                                     proj4string=crs(wi_counties))
    pts_intersect <- over(fullyadjusted_all_spdf, as(wi_counties, "SpatialPolygons"))
    fullyadjusted_all_counties <- fullyadjusted_all[which(!is.na(pts_intersect)),]
    fullyadjusted_all_counties$clusterid2.y <- as.factor(fullyadjusted_all_counties$clusterid2.y)

    fullyadjusted_all_counties <- fullyadjusted_all_counties %>% filter(expxbetahat!=1)%>% group_by(expxbetahat)  %>% mutate(uniqRRID = cur_group_id()) %>%ungroup()
    uniqRR <- table(fullyadjusted_all_counties$uniqRRID)
    colmap <- unlist(sapply(1:length(uniqRR), function(i) rep(colorsii[i], each=uniqRR[i])))
    namesmap <- unlist(lapply(1:length(uniqRR), function(i) rep(attributes(uniqRR)$dimnames[[1]][i], uniqRR[i])))
    colnam_df <- cbind.data.frame(color = colmap, names=namesmap) %>% mutate(color = as.factor(color),
                                                                             names = as.factor(names)) %>%
      distinct()
    
    fullyadjusted_all_counties <-merge(fullyadjusted_all_counties, colnam_df, by.x="uniqRRID", by.y="names")
    
    fullyadjusted_uniq_counties <- fullyadjusted_all_counties %>%
          distinct(expxbetahat,  .keep_all = TRUE) %>%
          arrange(-expxbetahat) %>%
          mutate(clusterid2.y = ifelse(is.na(clusterid2.y), "No Cluster", clusterid2.y),
                 womanid = ifelse(clusterid=="No Cluster", "No Cluster", womanid),
                 UBw = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)+(1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)),
                 LBw = exp(log(miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks)-(1.96*miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks))) %>%
        filter(womanid!="No Cluster")

    ####################################################################
    print(paste0("Length unique colors: ", length(unique(fullyadjusted_all_counties$clusterid2.y))))
    ##################################
    #Create plots
    ##################################
    fullyadjusted_all_counties_null <- fullyadjusted_all_counties %>% filter(expxbetahat==1) %>% droplevels()
    
    fullyadjusted_all_counties_or <- fullyadjusted_all_counties %>% filter(expxbetahat!=1) %>% droplevels() %>%
        arrange(-expxbetahat) %>%
      filter(color!="grey")
    print(paste0("Uniq womanid: ",sum(table(unique(fullyadjusted_uniq_counties$womanid))) ))
    
    if(map==TRUE){
        pout <- pcounties_map +
          lapply(1:length( clusters_bic_c), function(i) geom_polygon(data =  clusters_bic_c[[i]], aes(lon, lat, group = ID),
                                                                     color="grey50",
                                                                     fill=NA)) +
            geom_point(data=fullyadjusted_all_counties_or,
                       aes(x=lon+rnorm(n=1, mean=0,sd=0.001), y=lat+rnorm(n=1, mean=0,sd=0.001), color=color), alpha=0.4,size=1) +
            scale_colour_identity() + 
            theme_map() +
            theme(plot.margin = unit(c(0,0,0,0), "cm"),
                  legend.position = "none") +
            coord_fixed(ratio=1.5)
    } else if(map==FALSE) {
      if(any(fullyadjusted_uniq_counties$UBw > 3) | any(fullyadjusted_uniq_counties$LBw < 0.3)){
        print("UBmax or LBmax")
        #set max to 3
        ix <- which(fullyadjusted_uniq_counties$UBw>3)
        fullyadjusted_uniq_counties$UBw[ix] <- 3
        fullyadjusted_uniq_counties$upperout <- NA
        fullyadjusted_uniq_counties$upperout[ix] <- 3
        #set min to 1/3 (0.33)
        ix <- which(fullyadjusted_uniq_counties$LBw<0.3)
        fullyadjusted_uniq_counties$LBw[ix] <- 0.3
        fullyadjusted_uniq_counties$lowerout <- NA
        fullyadjusted_uniq_counties$lowerout[ix] <- 0.3
        #if estimate lower than 0.3
        ix <- which(fullyadjusted_uniq_counties$expxbetahat<0.3)
        fullyadjusted_uniq_counties$expxbetahat[ix] <- 0.3
        fullyadjusted_uniq_counties$estout <- NA
        fullyadjusted_uniq_counties$estout[ix] <- 0.3
        
        pout <- ggplot(data=fullyadjusted_uniq_counties) +
          geom_hline(yintercept=1, col="grey10", size=1.5) +
          geom_point(aes(x=fct_inorder(womanid), y=expxbetahat, color = color, size=3)) +  
          geom_point(aes(x=fct_inorder(womanid), y=estout, color = color, size=3)) +  
          geom_errorbar(aes(x=fct_inorder(womanid),ymin=LBw, ymax=UBw), colour="black", width=0) +
          geom_segment(aes(y = expxbetahat,x=fct_inorder(womanid), xend=fct_inorder(womanid), yend=upperout), 
                       arrow=arrow(length = unit(0.25, "cm"), type="closed"))+
          geom_segment(aes(y = expxbetahat,x=fct_inorder(womanid), xend=fct_inorder(womanid), yend=lowerout),
                       arrow=arrow(length = unit(0.25, "cm"), type="closed"))+

          scale_colour_identity() + 
          coord_fixed(ratio=0.5)+
          coord_flip() +
          

          theme_classic() +
          ylab("OR") +
          xlab("") +
          labs(fill=" ") +
          ggtitle("")+
          guides(size=FALSE) +
          theme(plot.title = element_text(hjust = 0.5, size=12),
                plot.caption=element_text(hjust=0.5, size=10),
                legend.position = "none",
                axis.text.y=element_blank(),
                axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))  +
          scale_y_continuous(breaks=c(0.3,0.5,  1,  2, 3),trans = scales::log_trans(),
                             limits=c(0.3,3)) 
       
        
      } else{
        pout <-ggplot(data=fullyadjusted_uniq_counties) +
          geom_hline(yintercept=1, col="grey10", size=1.5) +
          geom_point(aes(x=fct_inorder(womanid), y=expxbetahat, color = color, size=3)) +  
          geom_errorbar(aes(x=fct_inorder(womanid),ymin=LBw, ymax=UBw), colour="black", width=0) +
          scale_colour_identity() + 
          coord_fixed(ratio=0.5)+
          coord_flip() +
          theme_classic() +
          ylab("OR") +
          xlab("") +
          labs(fill=" ") +
          ggtitle("")+
          theme(plot.title = element_text(hjust = 0.5, size=12),
                plot.caption=element_text(hjust=0.5, size=10),
                legend.position = "none",
                axis.text.y=element_blank(),
                axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))  +
          scale_y_continuous(breaks=c(0.3,0.5,  1,  2, 3),trans = scales::log_trans(),
                             limits=c(0.3,3))  
      }
                
    } else {
      pout <-   countymap
    }
    return(pout)
}


# 
counties1 <- c("Douglas")
counties2 <- c("Ashland","Bayfield")
counties3 <- c("Polk","St. Croix")
counties4 <- c("Grant", "Lafayette")
counties5 <- c("Marathon")
counties6 <- c("Waukesha", "Milwaukee")

brewer.pal(n=8, name="Dark2")

colorsi <- brewer.pal(n=8, name="Dark2")
colorsi[7] <- "#666666"
colorsi[8] <- "#A6761D"


panel1_coef <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                               counties1, wi_undercount, wi, colorsi, map=FALSE)
panel1_map <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties1, wi_undercount, wi, colorsi, map=TRUE)
panel1_cty <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties1, wi_undercount, wi, colorsi, map="cty")
#10 unq cols
panel2_coef <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq,clustersdf_wisc_bic,
                               counties2, wi_undercount, wi, colorsi, map=FALSE)
panel2_map <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties2, wi_undercount, wi, colorsi, map=TRUE)
panel2_cty <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties2, wi_undercount, wi, colorsi, map="cty")
#6colors
panel3_coef <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                               counties3, wi_undercount, wi, colorsi, map=FALSE)
panel3_map <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties3, wi_undercount, wi, colorsi, map=TRUE)
panel3_cty <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties3, wi_undercount, wi, colorsi, map="cty")
#8 colors
panel4_coef <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                               counties4, wi_undercount, wi, colorsi, map=FALSE)
panel4_map <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties4, wi_undercount, wi, colorsi, map=TRUE)
panel4_cty <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties4, wi_undercount, wi, colorsi, map="cty")
#2

#panel5: 4 unique colors
panel5_coef <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq,clustersdf_wisc_bic,
                               counties5, wi_undercount, wi, colorsi, map=FALSE)
panel5_map <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq,clustersdf_wisc_bic,
                              counties5, wi_undercount, wi, colorsi, map=TRUE)
panel5_cty <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties5, wi_undercount, wi, colorsi, map="cty")


#panel6: 3 unique colors
#colorsii <- c(colorsi[1:2], "grey")
panel6_coef <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                               counties6, wi_undercount, wi, colorsi, map=FALSE)
panel6_map <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties6, wi_undercount, wi, colorsi, map=TRUE)
panel6_cty <- makemicropanels(fullyadjusted_all, fullyadjusted_uniq, clustersdf_wisc_bic,
                              counties6, wi_undercount, wi, colorsi, map="cty")


plt <- plot_grid(panel1_coef, panel1_map, panel1_cty,
                 panel2_cty, panel2_map, panel2_coef,
          panel3_coef, panel3_map, panel3_cty,
          panel4_cty, panel4_map, panel4_coef,
          panel5_coef, panel5_map, panel5_cty,
          panel6_cty, panel6_map, panel6_coef,
          ncol=6, 
          labels=c('A) Douglas', '','',
                   'B) Bayfield & Ashland', '','',
                   'C) Polk & St. Croix', '','',
                   'D) Grant & Lafayette', '','',
                   'E) Marathon', '','',
                   'F) Waukesha & Milwaukee', '',''))
plt


ggsave(file="pointplotuniqurelativerisksbycluster_6panel_low2_20210521.png",
       path = "../../manuscript/figures", dpi=72, width=14, height=10, units="in")

ggsave(file="pointplotuniqurelativerisksbycluster_6panel_low2_20210521.tiff",
       path = "../../manuscript/figures", dpi=300, width=14, height=10, units="in", compression = "lzw")

ggsave(file="pointplotuniqurelativerisksbycluster_6panel_low2_20210521_dpi300.png",
       path = "../../manuscript/figures", dpi=300, width=14, height=10, units="in")

ggsave(file="pointplotuniqurelativerisksbycluster_6panel_low2_20210521.eps",
       path = "../../manuscript/figures", dpi=72, width=14, height=10, units="in")

ggsave(file="pointplotuniqurelativerisksbycluster_6panel_low2_20210521.pdf",
       path = "../../manuscript/figures", dpi=72, width=14, height=10, units="in")


##############################################################################################################################################################################################################################################################################
rm(list=ls())