########################################################
#Additional Figures
#M.Kamenetsky
#update 2021-12-28
#added county names in Figure 3 that has clusters in them
#2021-02-01
########################################################
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
library(cowplot)
library(ggsn)
library(GISTools)
library(ggspatial)

table.out <- NULL
#############################
#helper funcs
#############################
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


coefestclusters <- function(coefdf, clustersdf_wisc_bic, idx){
    clustername <- clustersdf_wisc_bic$centerID[idx]
    mycircle <- make_circles(clustersdf_wisc_bic[idx,c(13,14,3)], radius=clustersdf_wisc_bic[idx,6])
    mycircle$clusterID <- paste0("cluster",mycircle$ID)
    mycircle$expest <- exp(coefdf$coefest[which(coefdf$clusterid==mycircle$clusterID[1])])
    return(mycircle)
}


#############################
#load data
#############################
bc <- read.csv(".")  #load in WWHS data (omitted here due to privacy)
#clean bc data a bit
bc$phase<- recode_factor(bc$phase, `regvar` = "phase 0")
bc$phase <- relevel(bc$phase, ref="phase 0")
bc <- bc %>%
    filter(!is.na(county))
########################################################
########################################################
#PLOT 1:
#Gridded milwaukee with circles from radius 0 through 10km
########################################################
########################################################
wi <- readOGR(dsn="../../data", layer="wi", verbose = TRUE)

#change crs to utm
wi_utm <- spTransform(x = wi, CRS("+init=epsg:3070"))
coords_utm <- coordinates(wi_utm)
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
undercounties <- c("Barron", "Bayfield","Buffalo","Burnett",
                   "Dunn","Eau Claire", "Pepin","Pierce","Polk", "St. Croix")
wi_undercount <- subset(wi, wi$NAME %in% undercounties)
undercounties_map <- fortify(wi_undercount)
p1 <- ggplot() + geom_map(data=wi_map, map=wi_map, 
                          aes(x=long, y=lat, map_id=id),
                          color="grey", fill="white", size=0.25) +
    coord_map() +
    theme_map() +
    ggsn::scalebar(data=wi_map, dist = 50, dist_unit = "km",
                          transform = TRUE, model = "WGS84", location = "topright", st.size=3.5)

#create layer with county names that have clusters in them
clustercounties <- c("Douglas", "Ashland", "Bayfield", "Polk", "St. Croix", "Grant",
                     "Lafayette", "Marathon", "Waukesha", "Milwaukee")
wi_clustercounties <- subset(wi, wi$NAME %in% clustercounties)
wi_clustercounties_sf<- st_as_sf(wi_clustercounties)




p2 <- p1 + geom_map(data=undercounties_map, map=undercounties_map,
                    aes(x=long, y=lat, map_id=id),
                    color="black",fill="white", size=1)

#############################################################################
mkenbd <- subset(wi, wi$NAME %in% c("Milwaukee", "Ozaukee", "Washington", "Waukesha", "Racine"))
mkenbd_map <- fortify(mkenbd)
mke <- subset(wi, wi$NAME %in% c("Milwaukee"))
mke_map <- fortify(mke)



wimkemap <- p2 + geom_map(data=mkenbd_map, map=mkenbd_map,
              aes(x=long, y=lat, map_id=id),
              fill="grey40", color="grey80") +
    geom_map(data=mke_map, map=mke_map,
             aes(x=long, y=lat, map_id=id),
             fill="black", color="grey80") +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering,
                           height= unit(1,"cm"),
                           width= unit(1,"cm"))  +
    ggtitle("Milwaukee County Neighborhood") +
    theme(plot.title = element_text(hjust = 0.5))
    


#######################################################################################
#For each county, we want to perform the cluster analysis
##Therefore, first extract a single county, and then I can loop over that
#######################################################################################
countiesgridexpand <- NULL
countynames <- as.vector(wi@data$NAME)
clusterslistqaic <- vector("list", 72)
clusterslistbic <- vector("list", 72)

wi_nb <- poly2nb(wi_utm, row.names = wi$NAME)
block1 <- wi[c(49, wi_nb[[49]]),]


print(centercountyname <- as.vector(wi@data[49,"NAME"]))
#subset to only counties of interest
counties1 <- as.vector(block1@data$NAME)
bc_quad <- bc %>%
    filter(county %in% counties1) %>%
    droplevels()

bc_focal <- bc %>%
    filter(county %in% centercountyname ) %>%
    droplevels()

#create nbd spdf
xyfocal <- cbind(bc_focal$x, bc_focal$y)
bc_focal_spdf <- SpatialPointsDataFrame(coords = xyfocal,
                                        data=bc_focal,
                                        proj4string = CRS("+init=epsg:3070"))

#join bc_sf and grid_50km
xy <- cbind(bc_quad$x, bc_quad$y)


#create spatial points dataframe
bc_spdf <- SpatialPointsDataFrame(coords = xy,
                                  data=bc_quad,
                                  
                                  proj4string = CRS("+init=epsg:3070"))
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
o = over(bc_spdf,gridspdf)
oo = over(bc_focal_spdf,gridspdf) 
head(o)
new = cbind(bc_spdf@data, o)
head(new)
ix <- which(o$id%in%oo$id)
new$focalcell <-0
new$focalcell[ix] <-1
#add if statement, if o is not empty, next to expand bbox xy
lengthoverlap <- length(which(is.na(new$id))) #want this to be zero
print(lengthoverlap)
if (lengthoverlap!=0){
    message(paste0("length overlap not equal to zero: ", centercountyname))
    countiesgridexpand <- c(countiesgridexpand, centercountyname)
    next
}

#find centroid of each grid cell
new$controls <- ifelse(new$bc==1,0,1)
truecentroids <- gCentroid(grid, byid = TRUE)
cent_df <- cbind.data.frame(as.data.frame(truecentroids@coords),
                            id = row.names(truecentroids@coords))


summary_xkm<- new%>%
    group_by(id, focalcell) %>%
    summarise(ncases = sum(bc),
              ncontrols = sum(controls)) %>%
    mutate(ntrials = ncases + ncontrols,
           phase = "phase1")%>%
    ungroup()



dat <- merge(summary_xkm, cent_df , by.x='id', by.y='id')
grid_df <-data.frame(id=character(), stringsAsFactors=FALSE )
for (m in grid@polygons ) { grid_df <- rbind(grid_df, data.frame(id=m@ID, stringsAsFactors=FALSE))  }
row.names(grid_df) <- grid_df$id
grid_spdf <- SpatialPolygonsDataFrame(grid, grid_df)
grid_lim <- grid_spdf[grid_spdf$id %in% dat$id,]



########################################
#convert spdf to sf
grid_lim_sf <- st_as_sf(grid_lim)

grid_summary <- merge(grid_lim_sf, summary_xkm, by.x="id", by.y="id")
grid_summary <- st_transform(grid_summary, st_crs(widatsf))
#find centroids of all cells
cellcentroids <- st_coordinates(st_centroid(grid_summary))
grid_summary$x <- cellcentroids[,1]
grid_summary$y <- cellcentroids[,2]
ixg3655 <- which(grid_summary$id=="g4335")
grid_summary <- grid_summary %>%
    mutate(centerID = id,
           latcentroids = y,
           loncentroids = x)

#make circles around grid cell "g3655"
rad1 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=1)
rad2 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=2)
rad3 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=3)
rad4 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=4)
rad5 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=5)
rad6 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=6)
rad7 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=7)
rad8 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=8)
rad9 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=9)
rad10 <- make_circles(grid_summary[ixg3655,c(10,11,12)], radius=10)

mkenbdsubset <- subset(widatsf, NAME %in%  c("Milwaukee", "Ozaukee", "Washington", "Waukesha", "Racine"))
mkenbd_bbox <- st_bbox(mkenbdsubset)
xrange <- mkenbd_bbox$xmax - mkenbd_bbox$xmin
yrange <- mkenbd_bbox$ymax - mkenbd_bbox$ymin

mkenbd_bbox[4] <- mkenbd_bbox[4] + (0.1 * yrange) # ymax - top
mkenbd_bbox[3] <- mkenbd_bbox[3] + (0.1 * xrange) # xmax - top
attr(mkenbd_bbox, "class") <- "bbox"
attr(st_geometry(mkenbdsubset), "bbox") <- mkenbd_bbox


pmke_grid <- ggplot() + 
        geom_sf(data=mkenbdsubset , 
                fill=NA, color="black", size=1)+
    geom_point(data=subset(grid_summary, focalcell==1), aes(x=x, y=y), col="black", size=0.5) +
    geom_sf(data=grid_summary,fill=NA) +
    geom_point(data=subset(grid_summary, id== "g4335"), color="black", aes(x=x,y=y)) +
    geom_polygon(data = rad1, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad2, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad3, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad4, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad5, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad6, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad7, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad8, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad9, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    geom_polygon(data = rad10, aes(lon, lat, group = ID), color = "black", alpha = 0, size=0.75) +
    theme_map()+
    ggsn::scalebar(data=mkenbdsubset, dist = 5, dist_unit = "km",
              transform = TRUE, model = "WGS84", location = "topright",st.size = 3) +
    annotation_north_arrow(location = "br", which_north = "true", 
                               pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
                               style = north_arrow_fancy_orienteering,
                               height= unit(1,"cm"),
                               width= unit(1,"cm")) +
    ggtitle("Potential Clusters for Single Grid Cell") +
    theme(plot.title = element_text(hjust = 0.5))

# ##########################################################
#Combine MKE NBD plot and grid plot
plot_grid(wimkemap, pmke_grid, labels = c('A)', 'B)'), ncol=2)    
ggsave(filename = "mkenbdgridcells.eps",path="../../manuscript/figures/", height = 5, width = 12, units="in", dpi=300)
ggsave(filename = "mkenbdgridcells.pdf",path="../../manuscript/figures/", height = 6, width = 12.5, units="in", dpi=300)

# ##########################################################
###########################################################
#from external: E:\BreastCancerWI\Results\Analysis\Wisconsin\tables
#also exists in C:\Users\Maria\Google Drive\RESEARCH\BreastCancerWI\Results\Analysis\Wisconsin\tables
miresults <- read.csv("../../Results/Analysis/Wisconsin/tables/dat_pooledmiresults_20210517.csv")
clustersdf_wisc_bic <- read.csv("../../Results/Analysis/Wisconsin/tables/clustersdf_wisc_bic.csv") %>%
    mutate(centerID = oldclusterid)
#############################################################################################
load("../../Results/Analysis/Wisconsin/RDataFiles/mi_coefse_age.RData")
load("../../Results/Analysis/Wisconsin/RDataFiles/mi_coefse_agefhnumbirthagefirstbagemenoraceeducbmidrinks.RData")

#CLUSTERS IDENTIFIED by BIC
p3 <- p1 + geom_map(data=undercounties_map, map=undercounties_map,
                    aes(x=long, y=lat, map_id=id),
                    color="black",fill="grey95") 

#micoefest_age
coefest_age_clusters <- data.frame(coefest=coefest_age[-c(1,17)])
coefest_age_clusters$clusterid <- attributes(coefest_age[-c(1,17)])$names
clusters_bic_age <- lapply(1:nrow(clustersdf_wisc_bic),
                           function(x) coefestclusters(coefest_age_clusters, clustersdf_wisc_bic, x))

#mi_coefse_agefhnumbirthagefirstbagemenoraceeducbmidrinks
coefest_full_clusters <- data.frame(coefest=coefest_agefhnumbirthagefirstbagemenoraceeducbmidrinks[-c(1:14,30:34)])
coefest_full_clusters$clusterid <- attributes(coefest_agefhnumbirthagefirstbagemenoraceeducbmidrinks[-c(1:14,30:34)])$names
clusters_bic_full <- lapply(1:nrow(clustersdf_wisc_bic),
                           function(x) coefestclusters(coefest_full_clusters, clustersdf_wisc_bic, x))

#fully adjusted
pbic_full<- p3 +
    lapply(1:length(clusters_bic_full), function(i) geom_polygon(data = clusters_bic_full[[i]], aes(lon, lat, group = ID,
                                                                                               fill=expest,
                                                                                               color=expest), 
                                                            alpha = 0.5)) +
    ggtitle("Fully-Adjusted Model Results") + 
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering,
                           height= unit(1,"cm"),
                           width= unit(1,"cm")) +
    labs(fill="OR")+
    theme(plot.title= element_text(hjust=0.5),
          legend.position = "bottom") +
    guides(color=FALSE) +
    scale_fill_gradient2(midpoint=0, low="darkblue", mid="white", high="red3",
                         breaks=c(0.15,0.25, 0.5,0.75, 1, 1.5, 2), limits=c(0.1, 2.1),trans = scales::log_trans()) +
    scale_color_gradient2(midpoint=0, low="darkblue", mid="white", high="red3", limits=c(0.1, 2.1),trans = scales::log_trans()) +
    geom_sf_text(data=subset(wi_clustercounties_sf, 
                             NAME %ni% c("Bayfield", "Milwaukee", "Waukesha")), 
                 aes(label=NAME)) +
    geom_sf_text(data=subset(wi_clustercounties_sf, 
                             NAME %in% c("Bayfield")), 
                 aes(label=NAME), nudge_y = 0.25)+
    geom_sf_text(data=subset(wi_clustercounties_sf, 
                             NAME %in% c("Milwaukee")), 
                 aes(label=NAME), nudge_y = 0.1, nudge_x = 0.2) +
    geom_sf_text(data=subset(wi_clustercounties_sf, 
                             NAME %in% c("Waukesha")), 
                 aes(label=NAME), nudge_x = -0.1, nudge_y = -0.1)

legend <- get_legend(pbic_full +theme(legend.position = "bottom", legend.key.width= unit(1, "cm"),
                                      legend.key.height = unit(0.25, "cm"),
                                legend.justification="center" ,
                                legend.text = element_text(angle = 45, vjust = 0.5, hjust=1)))
plot_row <- plot_grid(pbic_full + theme(legend.position="none"),
                      ncol=1,
                      align = 'h',
                      hjust=-1)

plot_grid(plot_row, legend, ncol=1,rel_heights = c(1,0.1))
ggsave("clustersidentified_fullyadjusted.pdf",
       path="../../manuscript/figures/", dpi=300,
       width=7, height=8, units="in")

