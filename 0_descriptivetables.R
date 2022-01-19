rm(list=ls())

########################################################
#2020-5-01
#M.Kamenetsky
#create descriptive tables (both imputed and not imputed)
#create overall tables and tables separated out by CCS status

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
library(gt)
library(tableone)
library(geosphere)

table.out <- NULL

'%ni%' <- Negate('%in%')

############################################################################################
#Load in csvs to process
############################################################################################
bc <- read.csv(".") #load in WWHS data (omitted here due to privacy)

#what are the response rates?
test <- bc %>%
    group_by(bc) %>%
    mutate_at(c(15:27), .funs = funs(response =  case_when( !is.na(.)  ~ 1,
                                                                             is.na(.)  ~ 0 ))) %>%
    mutate(response = rowSums(across(c(27:39)))) %>%
    mutate(response = ifelse(response!=0,1,0)) %>%
    summarize(response_pct = mean(response))


#clean bc data a bit
bc$phase<- recode_factor(bc$phase, `regvar` = "phase 0")
bc$phase <- relevel(bc$phase, ref="phase 0")
table(bc$bc)
length(which(is.na(bc$county)))

bc <- bc %>%
    filter(!is.na(county))
table(bc$bc)
round(mean(bc$age),2)
round(sd(bc$age),2)
round(mean(bc$famhist, na.rm=TRUE),2)
round(mean(bc$number_births[bc$gave_birth==1],na.rm=TRUE),2)
round(sd(bc$number_births[bc$gave_birth==1],na.rm=TRUE),2)
round(mean(bc$age_first_birth[bc$gave_birth==1],na.rm=TRUE),2)
round(sd(bc$age_first_birth[bc$gave_birth==1],na.rm=TRUE),2)
round(length(which(bc$race=="white"))/sum(unlist(table(bc$race))),2)
round(length(which(bc$educ=="grade 12"))/sum(unlist(table(bc$educ))),2)


#by BC status
tapply(bc$age, bc$bc, function(x) round(mean(x),1))
tapply(bc$age, bc$bc, function(x) round(sd(x),1))
tapply(bc$famhist, bc$bc, function(x) round(mean(x, na.rm=TRUE),1))
tapply(bc$number_births, bc$bc, function(x) round(mean(x[bc$gave_birth==1],na.rm=TRUE),1))
tapply(bc$number_births, bc$bc, function(x) round(sd(x[bc$gave_birth==1],na.rm=TRUE),1))
tapply(bc$age_first_birth, bc$bc, function(x) round(mean(x[bc$gave_birth==1],na.rm=TRUE),1))
tapply(bc$age_first_birth, bc$bc, function(x) round(sd(x[bc$gave_birth==1],na.rm=TRUE),1))

round(length(which(bc$race=="white" & bc$bc==1))/length(which(bc$bc==1)),2)
round(length(which(bc$race=="white" & bc$bc==0))/length(which(bc$bc==0)),2)
round(length(which(bc$educ=="grade 12" & bc$bc==1))/length(which(bc$bc==1)),2)
round(length(which(bc$educ=="grade 12" & bc$bc==0))/length(which(bc$bc==0)),2)


############################################################################################
#Table 1: with imputed 
##########################################################################
bc_clean <- bc %>%
    dplyr::select(bc, stage, age, race, educ, bmi, weight, height, drinker, drinks_per_week,
                  famhist, menopause, age_menopause, hrt,
                  gave_birth, number_births, age_first_birth)
bc_clean$gave_birth <- as.factor(bc_clean$gave_birth)


bc_clean$stage <- factor(bc_clean$stage, labels=c("Control",
                                                  "Distant",
                                                  "Local",
                                                  "Regional"))
bc_clean$race <- factor(bc_clean$race, labels=c("Black",
                                                "Hispanic",
                                                "Other",
                                                "White"))
bc_clean$educ <- factor(bc_clean$educ, labels=c("College 1-3",
                                                "College 4",
                                                "Doctor",
                                                "Grade 1-7",
                                                "Grade 12",
                                                "Grade 8",
                                                "Grade 9-11",
                                                "Masters",
                                                "None"))
bc_clean$educ <- ordered(bc_clean$educ, 
                         levels = c("None",
                                    "Grade 1-7",
                                    "Grade 8",
                                    "Grade 9-11",
                                    "Grade 12",
                                    "College 1-3",
                                    "College 4",
                                    "Masters",
                                    "Doctor"))

bc_clean$hrt <- factor(bc_clean$hrt, labels=c("Current",
                                              "Former",
                                              "Never"))
bc_clean$gave_birth <- factor(bc_clean$gave_birth, labels=c("No",
                                                            "Yes"))
bc_clean$bc <- factor(bc_clean$bc, labels=c("Controls", "Cases"))

bc_clean$drinker <- factor(bc_clean$drinker,
                                   labels=c("No", "Yes"))

#clean up labels
labelled::var_label(bc_clean$stage) <- "Stage"
labelled::var_label(bc_clean$age) <- "Age (years)"
labelled::var_label(bc_clean$race) <- "Race"
labelled::var_label(bc_clean$educ) <- "Education"
labelled::var_label(bc_clean$bmi) <- "BMI (kg/m2)"
labelled::var_label(bc_clean$weight) <-"Weight (kg)"
labelled::var_label(bc_clean$height) <- "Height (m)"
labelled::var_label(bc_clean$age_menopause) <- "Age at Menopause (years)"
labelled::var_label(bc_clean$hrt) <- "HRT"
labelled::var_label(bc_clean$gave_birth) <- "Gave Birth"
labelled::var_label(bc_clean$number_births) <- "Parity"
labelled::var_label(bc_clean$age_first_birth) <- "Age at First Birth (years)"
labelled::var_label(bc_clean$drinker) <- "Drinker"
labelled::var_label(bc_clean$drinks_per_week) <- "Drinks per Week"
labelled::var_label(bc_clean$bc) <- "Case/Control Status"
# 

############################################################################################
#Everything
############################################################################################

contvars <- c("age", "bmi", "height", "weight", "age_menopause",
              "number_births", "age_first_birth", "drinks_per_week", "stage",
              "race", "educ", "hrt", "gave_birth", "drinker")
catVars <- c("stage","race", "educ", "hrt", "gave_birth", "drinker")

tab1 <- CreateTableOne(vars=contvars, data=bc_clean, factorVars=catVars,
                       strata="bc", test=FALSE)
out1 <- print(tab1, showAllLevels=TRUE, varLabels=TRUE, missing=TRUE, noSpaces=TRUE, printToggle=FALSE)

kable(out1, format="latex")


############################################################################################
#Table 1 w/o imputations
############################################################################################
#Create small tables to be smooshed together based on problem variables
#Vars: 
#subt1
#- age_menopause*menopause
#- hrt*menopause

#subt2
#- number_births*gave_birth
#- age_first_birth*gave_birth

#subt3
#- drinks_per_week*drinker


nonconditional <-  c("age", "bmi", "height", "weight", 
                     "stage",
                     "race", "educ")
catVars <- c("stage","race", "educ")
taba <- CreateTableOne(vars=nonconditional, data=bc_clean, factorVars=catVars,
                       strata="bc", test=FALSE)
out_nonconditional <- print(taba, showAllLevels=TRUE, varLabels=TRUE, missing=TRUE, noSpaces=TRUE, printToggle=FALSE)
ix <-c(seq(10,50,2),53,54)
lines_uncond <- kable(out_nonconditional, format="latex",
                      caption = "Descriptive statistics for Wisconsin Women's Breast Cancer Study.") %>%
                capture.output %>%  .[-ix] 

#conditional on menopause
bc_clean_meno <- bc_clean %>% dplyr::filter(menopause==1)
meno <-  c("age_menopause","hrt")
catVars <- c("hrt")
tabb <- CreateTableOne(vars=meno, data=bc_clean_meno, factorVars=catVars,
                       strata="bc", test=FALSE)
out_meno <- print(tabb, showAllLevels=TRUE, varLabels=TRUE, missing=TRUE, noSpaces=TRUE, printToggle=FALSE)
ix <-c(1,2,3,seq(7,13,2),16)
lines_meno <- kable(out_meno, format="latex") %>%capture.output %>%  .[-ix]

#condiional on gave_birth
bc_clean_birth <- bc_clean %>% dplyr::filter(gave_birth=="Yes")
birth <-  c("number_births","age_first_birth")
tabc <- CreateTableOne(vars=birth, data=bc_clean_birth, 
                       strata="bc", test=FALSE)
out_birth<- print(tabc, showAllLevels=TRUE, varLabels=TRUE, missing=TRUE, noSpaces=TRUE, printToggle=FALSE)
ix <-c(1,2,3,seq(7,9,2),12)
lines_birth<- kable(out_birth, format="latex") %>%capture.output%>%  .[-ix]


#conditional on drinker
bc_clean_drink <- bc_clean %>% dplyr::filter(drinker=="Yes")
drink <-  c("drinks_per_week")
tabd <- CreateTableOne(vars=drink, data=bc_clean_drink, 
                       strata="bc", test=FALSE)
out_drink<- print(tabd, showAllLevels=TRUE, varLabels=TRUE, missing=TRUE, noSpaces=TRUE, printToggle=FALSE)
ix <- c(1:6,10)
lines_drink<- kable(out_drink, format="latex", caption="", title="test") %>%capture.output %>%  .[-ix]

#combine tables
combo <- c(lines_uncond,lines_meno,lines_birth, lines_drink) 
#clean combo
headers_ix <- c(7, 32,40, 46)
header <- "  &  & Controls & Cases & Missing\\\\" 
combo[headers_ix] <- header
#alter specific headers now
h_ix <-7
hmeno_ix <- 32
hbirth_ix <-40
hdrink_ix <- 46

h_all <- "   \\textbf{All Women} &  & \\textbf{Controls} & \\textbf{Cases} & \\textbf{Missing}\\\\" 
hmeno <- "   \\textbf{Post-Menopausal Women} &  & \\textbf{Controls} & \\textbf{Cases} & \\textbf{Missing}\\\\" 
hbirth <- "   \\textbf{Gave Birth} &  & \\textbf{Controls} & \\textbf{Cases} & \\textbf{Missing}\\\\" 
hdrink <- "   \\textbf{Self-Reported Drinkers} &  & \\textbf{Controls} & \\textbf{Cases} & \\textbf{Missing}\\\\" 

combo[hmeno_ix] <-hmeno
combo[hbirth_ix] <- hbirth
combo[hdrink_ix] <- hdrink
combo[h_ix] <- h_all
#######################################
#Output
#######################################
combo %>%
    cat(file="../../manuscript/tables/table1.tex",sep = "\n")



