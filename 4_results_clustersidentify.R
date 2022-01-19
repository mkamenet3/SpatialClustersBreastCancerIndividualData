#M.Kamenetsky
#2021-5-21
#Identify the clusters from the models
#Not used in final manuscript
######################################################################################
#######################################################################################
#######################################################################################
#PART 1: Identify ALL Clusters for Plotting
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################
miresults <- read.csv("../../Results/Analysis/Wisconsin/tables/dat_pooledmiresults_20210517.csv")

###################################################################################################
#BIC RESULTS
###################################################################################################
clusterids <- names(miresults)[17:31]
clusterids <- as.numeric(unlist(regmatches(clusterids, gregexpr("[[:digit:]]+", clusterids))))


countieswithclusters <- read.csv("../../Results/Analysis/Wisconsin/tables/allcounties.csv") %>%
    dplyr::filter(model=="Binomial" & numclust.BIC !=0) %>%
    dplyr::select(centercounty)

#load - to loop this
clusteridents <- NULL
potentialclusteridents <- NULL
for (i in 1:nrow(countieswithclusters)){
    load(paste0("../../Results/Analysis/Wisconsin/RDataFiles/", countieswithclusters$centercounty[i],".RData"))
    resclussoSparseMAT <- resclusso$sparseMAT
    resclussodatsub <- dat %>% filter(focalcell==1) %>% droplevels()
    subx <- resclussodatsub$x/1000
    suby <- resclussodatsub$y/1000
    
    ix_x <- which(resclusso$clustersdf$x %in% subx)
    ix_y <- which(resclusso$clustersdf$y %in% suby)
    ix <- intersect(ix_x,ix_y)
    clusterdfsub <- resclusso$clustersdf[ix,]
    
    #identified clusteres
    selected <- resclusso$lassoresult.p.s$selections$select.qbic
    clusterix <- which(resclusso$lassoresult.p.s$coefs.lasso.all[,selected]!=0)
    cluster <- cbind(clusterdfsub[clusterix,], nb=as.character(countieswithclusters$centercounty[i]),
                     oldclusterid = row.names(resclusso$clustersdf[clusterix,]))
    clusteridents <- rbind(clusteridents, cluster)
    #all potential clusters

}
#IDENTIFIED CLUSTERS
#Convert WTM coordinates to lat/lon
#create spatialpoints dataframe
clusteridents$x1000 <- clusteridents$x*1000
clusteridents$y1000 <- clusteridents$y*1000
clusteridents.points <- SpatialPointsDataFrame(cbind(clusteridents$x1000, clusteridents$y1000),
                                               clusteridents, proj4string = CRS("+init=epsg:3070"))
clusteridents.latlon <- spTransform(clusteridents.points, CRS("+proj=longlat"))

#######################################################################################
#Create circular clusters to overlay on top of maps (see Plotting below)
#######################################################################################
#IDENTIFIED CLUSTERS
clusteridents$loncentroids <- coordinates(clusteridents.latlon)[,1]
clusteridents$latcentroids <- coordinates(clusteridents.latlon)[,2]
clusteridents$clusterid <- row.names(clusteridents)
write.csv(clusteridents, file="../../Results/Analysis/Wisconsin/tables/clustersdf_wisc_bic.csv")

