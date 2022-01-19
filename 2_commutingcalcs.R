#M.Kamenetsky
#2021-10-02
#Purpose: This code calculates the mean and median distance travelled in Wisconsin
#Data Downloaded from Federal Highway Administration U.S. Department of Transporation. 2017 National Household Travel Survey. Technical report, 2017.
#https://nhts.ornl.gov/

#load packages
library(tidyverse)
library(survey)

#load dataset
dat <- read.csv("../../data/csv/perpub.csv")
wts <- read.csv("../../data/ReplicatesCSV/perwgt.csv")

wts2 <- wts%>%
    dplyr::select(HOUSEID,PERSONID, WTPERFIN)
dat2 <- dat %>% 
    dplyr::select(HOUSEID, PERSONID,GCDWORK,HBHUR, HHSTATE) %>%
    filter(GCDWORK !=-9 & HBHUR!="-9" & HHSTATE=="WI") 

main <-merge(dat2, wts2, by=c("HOUSEID", "PERSONID")) %>%
    mutate(ID = paste0(HOUSEID,PERSONID))

nhts <- svydesign(id = ~ID,
                  strata = NULL,
                  weights = ~WTPERFIN,
                  nest = FALSE,
                  data = main)
summary(nhts)

svymean(~GCDWORK, nhts)
svyquantile(~GCDWORK, nhts, quantile=c(0.5), ci=TRUE)

main%>%
    summarize(mean= mean(GCDWORK),
              median=median(GCDWORK)) %>%
    arrange(median)
