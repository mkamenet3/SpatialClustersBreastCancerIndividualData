rm(list=ls())
########################################################
#M.Kamenetsky
#2020-02-01
#The goal of this script is to explore how to
#create neighbors for each county and block
#by each county (and its neighbors) and perform
#the same SER analysis.

#this script is executed from the main BreastCancerWI directory
#Update 2020-02-08: 
#   - rerun analysis because need new and dat created by
#   each county in order to merge with original bc
#   dataset at the person level (identify which people are in
#   which grid cell, but grid cell id is only created by
#   each neighborhood)
#Update 2020-02-14:
#   - Create a map that summarizes teh RR clusters from the model
#       - For each control location, calculate a cluster effect (elevated RR could come from cluster overlap, even if
#          CI's for the individual cluster span 1)
#       - Need to pick some sort of background location
#       - there should be a smaller set of CI's for people characterized by the cluster
#       - find all unique risks for points in clusters; locations for control persons & plot; summarize over them
#       - extract coefs and vcov matrix of logistic model; set other coefs to zero (that are not the cluster)


########################################################
#############################
#load required packages
#############################
library(spdep)
library(tidyverse)
library(clusso)
library(rgdal)
library(rgeos)
library(cowplot)
library(knitr)
library(RColorBrewer)
library(gridExtra)
library(ggforce)
library(sf)
library(raster)
library(Matrix)


table.out <- NULL

'%ni%' <- Negate('%in%')

#############################
#Helper functions
#############################
gridplots <- function(lassoresult.ic, summary_xkm, grid_lim, title, gridmax, rmax, ic, centercountyname){
    maxrr <- 2
    minrr <- 0.5
    estrr.qic <- lassoresult.ic
    cluster_ix <- redblue(log(maxrr *  pmax(minrr, pmin(estrr.qic, maxrr)))/log(maxrr^2))
    rr.qic <- cbind.data.frame(rr=lassoresult.ic, 
                               id=summary_xkm$id, 
                               colors = cluster_ix)
    names(cluster_ix) <- rr.qic$id
    
    
    
    merged.qic <- merge(summary_xkm, rr.qic, by="id")
    row.names(merged.qic) <- merged.qic$id
    merged.qic_spdf <- SpatialPolygonsDataFrame(as(grid_lim, "SpatialPolygons"), data=merged.qic, match.ID = TRUE)
    #here, merge to 
    datdf.qic <- st_as_sf(merged.qic_spdf)
    #RRmap
    ggplot()+ geom_sf(data=datdf.qic, aes(fill=id)) +
        scale_fill_manual(values=cluster_ix)+
        theme_bw()+
        theme(legend.position = "none",
              plot.title = element_text(hjust=0.5)) +
        ggtitle(title)
    ggsave(paste0(centercountyname,"grid",gridmax,"_rad", rmax,"_ic",ic, "_space.PNG"),
           path="../../Results/Analysis/Wisconsin/figures")
    #ident map
    ggplot()+ geom_sf(data=subset(datdf.qic, rr==1), fill='transparent') +
        geom_sf(data=subset(datdf.qic, rr<1), color="blue") +
        geom_sf(data=subset(datdf.qic, rr>1), color="red") +
        theme_bw()+
        theme(legend.position = "none",
              plot.title = element_text(hjust=0.5)) +
        ggtitle(title)
    #ggtitle("Quasi-Binomial, Space: (Q)BIC") 
    ggsave(paste0(centercountyname,"grid",gridmax,"_rad", rmax,"_ic", ic,"_space_ident.PNG"),
           path="../../Results/Analysis/Wisconsin/figures")
}

#############################
#load data
#############################
bc <- read.csv(".") #load in WWHS data (omitted here due to privacy)
#clean bc data a bit
bc$phase<- recode_factor(bc$phase, `regvar` = "phase 0")
bc$phase <- relevel(bc$phase, ref="phase 0")
bc <- bc %>%
    filter(!is.na(county))


#######################################################################################
#######################################################################################
#######################################################################################
#PART 1: clusso Each County Neighborhood in Wisconsin
#######################################################################################
#######################################################################################
#######################################################################################
#load Wisconsin map
wi <- readOGR(dsn="../../data", layer="wi", verbose = TRUE)
row.names(wi) <- row.names(wi@data)
head(wi@data)
wi$id <- 1:nrow(wi)
wi_map <- fortify(wi)
p1 <- ggplot() + geom_map(data=wi_map, map=wi_map, 
                          aes(x=long, y=lat, map_id=id),
                          color="black", fill="white", size=0.25) +
    coord_map() +
    theme_bw()
p1
#extract centroids from map
#convert to simple features sf
widatsf <- st_as_sf(wi)
coords <- as.data.frame(st_coordinates(st_centroid(widatsf)))
head(coords)
p1 + geom_point(data=coords,aes(x=X, y=Y), color="red", shape=3)


pborders <- ggplot() + geom_map(data=wi_map, map=wi_map, 
                                aes(x=long, y=lat, map_id=id),
                                color="black", fill="grey90", size=1, alpha=0.5) +
    coord_map()

#############################
#create queen's case neighborhood
#############################
wi_nb <- poly2nb(wi, row.names = wi$NAME)
str(wi_nb)
#create row-standardized weights
listw_wi <- nb2listw(wi_nb, style="W")
plot(wi);plot(listw_wi, coords, col="blue", add=TRUE)

#For example, Grant is neighbors with Crawford (4), Lafayette,(47) and some others
#View(cbind(wi@data, 1:72))

#######################################################################################
#For each county, we want to perform cluster analysis
##Therefore, first extract a single county, and then I can loop over that
#######################################################################################
countiesgridexpand <- NULL
countynames <- as.vector(wi@data$NAME)
clusterslistqaic <- vector("list", 72)
clusterslistbic <- vector("list", 72)



# Start the clock!
ptm <- proc.time()
for (i in 1:nrow(wi@data)) {
    block1 <- wi[c(i,wi_nb[[i]]),]
    print(centercountyname <- as.vector(wi@data[i,"NAME"]))
    #subset to only counties of interest
    counties1 <- as.vector(block1@data$NAME)
    
    bc_focal <- bc %>%
        filter(county %in% centercountyname ) %>%
        droplevels()
    
    
    bc_quad <- bc %>%
        filter(county %in% counties1) %>%
        droplevels()

    #create nbd spdf
    xyfocal <- cbind(bc_focal$x, bc_focal$y)
    bc_focal_spdf <- SpatialPointsDataFrame(coords = xyfocal,
                                      data=bc_focal,
                                      proj4string = CRS("+init=epsg:3070"))
    xy <- cbind(bc_quad$x, bc_quad$y)
    bc_spdf <- SpatialPointsDataFrame(coords = xy,
                                            data=bc_quad,
                                            proj4string = CRS("+init=epsg:3070"))
    #create focal county spdf 
    a <- bbox(bc_spdf)
    x <- seq(from = a[1,1]-1000, to = a[1,2]+1000, by = 1000)
    y <- seq(from = a[2,1]-1000 , to = a[2,2]+1000, by = 1000)
    xy <- expand.grid(x = x, y = y)
    
    grid.pts<-SpatialPointsDataFrame(coords= xy, data=xy, 
                                     proj4string = CRS("+init=epsg:3070"))
    class(grid.pts)
    gridded(grid.pts) <- TRUE
    grid <- as(grid.pts, "SpatialPolygons")
    summary(grid)
    gridspdf <- SpatialPolygonsDataFrame(grid, data=data.frame(id=row.names(grid),
                                                               row.names=row.names(grid)))
    names.grd<-sapply(gridspdf@polygons, function(x) slot(x,"ID"))
    text(coordinates(gridspdf), labels=sapply(slot(gridspdf, "polygons"), function(m) slot(m,
                                                                                           "ID")), cex=0.3)
    o = over(bc_spdf,gridspdf) #overlaying grid on top of nbd
    oo = over(bc_focal_spdf,gridspdf) 
    head(o) #these are grid IDs for the whole nbd
    head(oo) #these are grid IDs for the focal county only
    ix <- which(o$id%in%oo$id)
    new = cbind(bc_spdf@data, o) #ok this is a the sub-data set of all grid cells in nbd
    new$focalcell <-0
    new$focalcell[ix] <-1
    head(new)
    #add if statement, if o is not empty, next to expand bbox xy
    lengthoverlap <- length(which(is.na(new$id))) #want this to be zero
    print(lengthoverlap)
    if (lengthoverlap!=0){
        message(paste0("length overlap not equal to zero: ", centercountyname))
        countiesgridexpand <- c(countiesgridexpand, centercountyname)
        next
    }
    
    #find centroid of each grid cell from the focal county 
    new$controls <- ifelse(new$bc==1,0,1)
    truecentroids <- gCentroid(grid, byid = TRUE)
    cent_df <- cbind.data.frame(as.data.frame(truecentroids@coords),
                              id = row.names(truecentroids@coords))
    #subset cent_df to only centers in focal county
    summary_xkm<- new%>%
        group_by(id, focalcell) %>%
        summarise(ncases = sum(bc),
                  ncontrols = sum(controls)) %>%
        mutate(ntrials = ncases + ncontrols,
               phase = "phase1")%>%
        ungroup()
    
    dat <- merge(summary_xkm, cent_df , by.x='id', by.y='id')
    dat$id <- as.factor(dat$id)
    ############################
    #perform clusso
    ############################
    #limit to only focal cells as cluster ids
    datsub <- dat %>% filter(focalcell==1) %>% droplevels()
    #set params for cluster detection via clusso
    rMax <- 10
    Time <- 1
    xx <- dat$x/1000
    yy <- dat$y/1000
    subx <- datsub$x/1000
    suby <- datsub$y/1000
    focalcellsix <- which(dat$focalcell==1)

    resclusso <- clusso(df=dat,
                        expected = ntrials,
                        observed = ncases,
                        timeperiod = phase,
                        id = id,
                        covars = FALSE,
                        x= xx,
                        y= yy,
                        rMax = rMax,
                        utm=TRUE,
                        analysis = "space",
                        model="binomial",
                        maxclust = 15, collapsetime = TRUE, focalcells = focalcellsix,subxy = cbind(subx, suby))
    
    save(resclusso, dat, new, file = paste0("../../Results/Analysis/Wisconsin/RDataFiles/",centercountyname, ".RData"))
    (clussoprettyout <- clussopretty(resclusso, analysis="space", model="binomial",clusteridentify=FALSE))
    clussoprettyout$centercounty <- centercountyname
    table.out <- rbind(table.out, clussoprettyout)
    print(table.out)
    write.csv(table.out, file=paste0("../../Results/Analysis/Wisconsin/tables/counties/",centercountyname,".csv"),
              row.names = FALSE)
    #######################################################################################
    #IDENTIFY PEOPLE IN CLUSTERS and BUILD X matrix based on BIC RESULTS
    #######################################################################################
    ###############################
    #BIC
    ###############################
    #selected PC:
    identPC.bic <- resclusso$lassoresult.p.s$coefs.lasso.all[, resclusso$lassoresult.p.s$selections$select.qbic] #select PC
    gridincluster.bic <- resclusso$sparseMAT %*% matrix(identPC.bic,ncol=1)
    gridincluster_bin.bic <- ifelse(gridincluster.bic!=0,1,0)
    #these are which grid cells are inside the identified cluster area
    #The same order is preserved
    ix_gridincluster.bic <- which(gridincluster_bin.bic!=0)
    dat_ix_gridincluster.bic <- dat[ix_gridincluster.bic,]
    dat_ix_gridincluster.bic <- droplevels(dat_ix_gridincluster.bic) #these are the grid cells that are in a cluster
    ####################3
    #Ok, now need to idnetify which PC people belong to.
    #for now, just do it this way -> create wide dataframe
    cident_ix.bic <- which(identPC.bic!=0)
    clist.bic<- lapply(1:length(cident_ix.bic), function(x) rep(0,length(identPC.bic)))
    clistdf.bic <- vector("list", length(cident_ix.bic))
    print(cident_ix.bic)
    if (length(cident_ix.bic)==0){
        print(paste0("no BIC-based clusters found in ", centercountyname))
        clusterslistbic[[i]] <-  NULL
    }
    else{ for (j in 1:length(cident_ix.bic)){
        clist.bic[[j]][cident_ix.bic[[j]]] <-1
    }
        for (k in 1:length(cident_ix.bic)){
            c_gridincluster.bic <- resclusso$sparseMAT %*% matrix(clist.bic[[k]],ncol=1)
            c_ix_gridincluster.bic <- which(c_gridincluster.bic@x!=0)
            c_dat_ix_gridincluster.bic <- dat[c_ix_gridincluster.bic,]
            c_dat_ix_gridincluster.bic <- droplevels(c_dat_ix_gridincluster.bic)
            
            c_dat_ix_gridincluster.bic$clusterid <- paste0("cluster", cident_ix.bic[[k]])
            clistdf.bic[[k]] <- c_dat_ix_gridincluster.bic
        }
        clusters_bic <- dplyr::bind_rows(clistdf.bic, .id="column_label")
        clusters_bic$clusteryes <- 1
        clusters_bic$countyname <- countynames[[i]]
        #bind list of dataframes
        clusters_bic$clusterid <- as.factor(clusters_bic$clusterid)
        #go long to wide
        clusters_bic_w <- clusters_bic %>% 
            dplyr::select(id, phase,  clusterid, clusteryes) %>%
            pivot_wider(names_from = clusterid,
                        values_from = clusteryes,
                        values_fill = list(clusteryes=0))
        ccs_clusters_bic_w <- merge(new,clusters_bic_w, by.x="id", by.y="id") 
        #add to list that will then be collapsed
        clusterslistbic[[i]] <-  ccs_clusters_bic_w
    }
    #####################
    #Create GRIDPLOTS
    #####################
    gridplots
    gridplots(resclusso$lassoresult.p.s$E.qbic, summary_xkm, grid_lim = grid_lim,
              title="Binomial, Space: BIC", gridmax = 1, rmax=10, ic="bic", centercountyname)
}
# Stop the clock
proc.time() - ptm
print(countiesgridexpand)
save(clusterslistbic, file="../../Results/Analysis/Wisconsin/RDataFiles/clusterslistbic.RData")

###############################
#Complete dataframe
###############################
#####################
#BIC
#####################
#collapse list into df
clusters_bicdf <- dplyr::bind_rows(clusterslistbic) %>%
    mutate_at(vars(cluster100:cluster4375), ~replace_na(.,0)) 
#replace any NA with 0's
write.csv(clusters_bicdf,file="../../Results/Analysis/Wisconsin/tables/clusters_bicdf.csv")

#########################################################################################################
#POST-PROCESSING
#Since I look at the neighborhood level, there will be duplicate idnum but unique id
#The dataset should be unique at the idnum level
#Therefore, by idnum, carryforward all indicators of 1 in cluster variables
#then drop id (which is the grid id) and drop duplicates.
#This will then tell me whenever a person was in a cluster id, regardless of if it was the 
#center county or just a neighboring county.
#########################################################################################################
# #####################
# #BIC
# #####################
# 
clusters_bicdf <- read.csv("../../Results/Analysis/Wisconsin/tables/clusters_bicdf.csv")

namestofill <- names(clusters_bicdf)[33:47]
clusters_bicdf_uniqidnum <- clusters_bicdf %>%
    group_by(idnum) %>%
    fill(namestofill) %>%
    fill(namestofill, .direction = "up") %>%
    dplyr::select(-id) %>%
    distinct() %>%
    mutate_at(vars(cluster100:cluster4375), ~replace_na(.,0)) %>%
write.csv(clusters_bicdf_uniqidnum,file="../../Results/Analysis/Wisconsin/tables/clusters_bicdf_uniqidnum.csv")


#########################################################################################################
#Combine counties into main csv called allcounties.csv
#########################################################################################################
allcounties0 <- read.csv("../../Results/Analysis/Wisconsin/tables/counties/Sauk.csv") 
#Sauk was the last county in the loop and contains all the other county results
write.csv(allcounties0, file="../../Results/Analysis/Wisconsin/tables/allcounties.csv")


