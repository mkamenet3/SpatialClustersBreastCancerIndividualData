#M.Kamenetsky
#EDA

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
library(xtable)
library(geosphere)
library(ggsn)
library(ggspatial)

table.out <- NULL

'%ni%' <- Negate('%in%')

############################################################################################
#Load in csvs to process
############################################################################################
bc <- read.csv(".")  #load in WWHS data (omitted here due to privacy)
#clean bc data a bit
bc$phase<- recode_factor(bc$phase, `regvar` = "phase 0")
bc$phase <- relevel(bc$phase, ref="phase 0")
bc <- bc %>%
    filter(!is.na(county))

###########
################################################################################################################################
#PART 3: PLOTS
################################################################################################################################
################################################################################################################################
################################################################################################################################
######################################
#Plot Raw Case and Control Points
######################################
#Load maps
bc$agegroups <- cut(bc$age, breaks = 5)
cases <- bc %>% filter(bc==1)  %>% drop_na(c("lat","lon"))
controls <- bc %>% filter(bc==0) %>% drop_na(c("lat","lon"))
undercounties <- c("Barron", "Bayfield","Buffalo","Burnett",
                   "Dunn","Eau Claire", "Pepin","Pierce","Polk", "St. Croix")
wi <- readOGR(dsn="../../data", layer="wi", verbose = TRUE)

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
controlpoints <- p2 + geom_point(data=controls, aes(x=lon, y=lat),col="darkslateblue", size=1, alpha=0.1) +
    ggtitle("Controls") +
    theme(plot.title = element_text(hjust = 0.5, size=22)) 

casepoints <- p2 + geom_point(data=cases, aes(x=lon, y=lat),col="firebrick", size=1, alpha=0.1) +
    ggtitle("Cases") +
    theme(plot.title = element_text(hjust = 0.5, size=22)) 

plotout <- grid.arrange(casepoints, controlpoints, ncol = 2)
ggsave(filename = "caseandcontrollocationsacrossphase.pdf", plot=plotout,path="../../manuscript/figures/", height = 5, width = 6, units="in", dpi=300)
ggsave(filename = "caseandcontrollocationsacrossphase.eps", plot=plotout,path="../../manuscript/figures/", height = 5, width = 6, units="in", dpi=300)


######################################
#Update maps for cases and controls
######################################
#convert wi to wisf
wisf <- st_transform(st_as_sf(wi),"+init=epsg:3070")
wisf_centroids <- st_coordinates(st_transform(st_centroid(wisf),"+init=epsg:4326"))
wisf$x <- wisf_centroids[,1]
wisf$y <- wisf_centroids[,2]
#sum cases and controls, resp
countycases <- cases %>%
    group_by(county) %>%
    summarize(casesTot =  sum(bc))
countycontrols <- controls %>%
    group_by(county) %>%
    summarize(controlsTot =  sum(ifelse(bc==0,1,0)))
#merge into wisf
wisf <- merge(wisf, countycases, by.x="NAME", by.y="county")
wisf <- merge(wisf, countycontrols, by.x="NAME", by.y="county")


#make 
pcases <- p2 +
    geom_point(data=wisf, aes(x=x, y=y, size=casesTot), shape=16) +
    labs(size="", title="Cases") +
    theme( legend.position = 'bottom',
           legend.spacing.x = unit(0.75, 'cm'),
          plot.title = element_text(hjust = 0.5)) +
    ggsn::scalebar(wi_map, dist=50, location="topright", st.size=2, transform=TRUE, dist_unit="km") +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) 
pcontrols <- p2 +
    geom_point(data=wisf, aes(x=x, y=y, size=controlsTot), shape=15) +
    labs(size="", title="Controls") +
    theme( legend.position = 'bottom',
           legend.spacing.x = unit(1, 'cm'),
           plot.title = element_text(hjust = 0.5)) +
    scalebar(wi_map, dist=50, location="topright", st.size=2, transform=TRUE, dist_unit="km") +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) 
plot_grid(pcases, pcontrols, labels = c('A)', 'B)'))    
ggsave(filename = "casesandcontrols_circlesquare.eps",path="../../manuscript/figures/", height = 5, width = 12, units="in", dpi=300)
ggsave(filename = "casesandcontrols_circlesquare.pdf",path="../../manuscript/figures/", height = 5, width = 12, units="in", dpi=300)



