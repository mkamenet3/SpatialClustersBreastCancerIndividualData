#M.Kamenetsky
#Updated 2020-03-31
#Impute Missing values for key covariates
#This script is run from the root BreastCancerWI directory
#############################
#load required packages
#############################
library(Hmisc)
library(tidyverse)
table.out <- NULL
set.seed(04082020)
'%ni%' <- Negate('%in%')

############################################################################################
#Helper Functions
############################################################################################
cleanimpute <- function(bc, bc_mi, imputeid){
    bc_Imp<- bc
    imputed<- impute.transcan(bc_mi,imputation=imputeid,data=bc,list.out=TRUE,pr=FALSE,check=FALSE)
    bc_Imp[names(imputed)] <- imputed
    
    bc_Imp <- bc_Imp %>%   
        mutate(drinks_per_week=ifelse(drinker==1,drinks_per_week,2.8),
               number_births=ifelse(gave_birth==1,number_births,3),
               age_first_birth=ifelse(gave_birth==1,age_first_birth,24),
               age_menopause=ifelse(menopause==1,age_menopause,48),
               hrt=factor(ifelse(menopause==1,hrt,1),labels=c("never","former","current")))
    return(bc_Imp)
}

checkimpute <- function(impi, groupvar, sumvar){
    impi %>%
        group_by(!!! groupvar) %>%
        summarise(mean = mean(!!! sumvar),
                  sd = sd(!!!sumvar),
                  median = median(!!!sumvar),
                  q25 = quantile(!!!sumvar, 0.25),
                  q75 = quantile(!!!sumvar, 0.75))
}

checkimpbc <- function(impi){
    print(checkimpute(impi, quos(drinker), quos(drinks_per_week)))
    print(checkimpute(impi, quos(menopause), quos(age_menopause)))
    print(with(impi, table(menopause, hrt)))
    print(checkimpute(impi, quos(gave_birth), quos(number_births)))
    print(checkimpute(impi, quos(gave_birth), quos(age_first_birth)))
}
############################################################################################
#Load in csvs to process
############################################################################################

bc <- read.csv(".")#load in WWHS data (omitted here due to privacy)
#clean bc data a bit
bc$phase<- recode_factor(bc$phase, `regvar` = "phase 0")
bc$phase <- relevel(bc$phase, ref="phase 0")

#exclude if county missing:
bc <- bc %>%
    dplyr::filter(!is.na(county))

#figure out % missing for key vars
bc %>%
    group_by(drinker) %>%
    summarise(napercent_drinksperweek = mean(is.na(drinks_per_week))) 
bc %>%
    group_by(menopause) %>%
    summarise(napercent_agemeno = mean(is.na(age_menopause)),
              napercent_hrt = mean(is.na(hrt))) 
bc %>%
    group_by(gave_birth) %>%
    summarise(napercent_numbirths = mean(is.na(number_births)),
              napercent_agefirstbirth = mean(is.na(age_first_birth))) 

############################################################################################
#Imputation
#aregImpute = multiple imputation using additive regression, bootstrapping, and predictive mean matching
#binary variables are restricted to be linear
#default is predictive mean matching
#14 imputations because 14% highest missing variable
############################################################################################
#Impute
bc_mi <- aregImpute(~x+y+phase+bc+stage+refyr+age+race+educ+bmi+drinker+drinks_per_week+famhist+menopause+age_menopause+hrt+
                                       gave_birth+number_births+age_first_birth, 
                                   n.impute=14, tlinear=FALSE, data=bc)
#add imputations to dataframe
bc_Imp1 <- cleanimpute(bc, bc_mi, 1)
bc_Imp2 <- cleanimpute(bc, bc_mi, 2)
bc_Imp3 <- cleanimpute(bc, bc_mi, 3)
bc_Imp4 <- cleanimpute(bc, bc_mi, 4)
bc_Imp5 <- cleanimpute(bc, bc_mi, 5)
bc_Imp6 <- cleanimpute(bc, bc_mi, 6)
bc_Imp7 <- cleanimpute(bc, bc_mi, 7)
bc_Imp8 <- cleanimpute(bc, bc_mi, 8)
bc_Imp9 <- cleanimpute(bc, bc_mi, 9)
bc_Imp10 <- cleanimpute(bc, bc_mi, 10)
bc_Imp11 <- cleanimpute(bc, bc_mi, 11)
bc_Imp12 <- cleanimpute(bc, bc_mi, 12)
bc_Imp13 <- cleanimpute(bc, bc_mi, 13)
bc_Imp14 <- cleanimpute(bc, bc_mi, 14)

############################################################################################
#Check Missing
##1: check % missing
##2: Check that everyone who has drinks per week is not missing
###These variables are: 
#           - drinker-drinks_per_week
#           - menopause-age_menopause
#           - hrt - menopause
#           - gave_birth - number_births
#           - gave_birth - age_first_birth
############################################################################################
#bc_Imp1
checkimpbc(bc_Imp1)
checkimpbc(bc_Imp2)
checkimpbc(bc_Imp3)
checkimpbc(bc_Imp4)
checkimpbc(bc_Imp5)
checkimpbc(bc_Imp6)
checkimpbc(bc_Imp7)
checkimpbc(bc_Imp8)
checkimpbc(bc_Imp9)
checkimpbc(bc_Imp10)
checkimpbc(bc_Imp11)
checkimpbc(bc_Imp12)
checkimpbc(bc_Imp13)
checkimpbc(bc_Imp14)

############################################################################################
#save output
############################################################################################
save(bc_Imp1,file="data/imputations/bc_Imp1.RData")
save(bc_Imp2,file="data/imputations/bc_Imp2.RData")
save(bc_Imp3,file="data/imputations/bc_Imp3.RData")
save(bc_Imp4,file="data/imputations/bc_Imp4.RData")
save(bc_Imp5,file="data/imputations/bc_Imp5.RData")
save(bc_Imp6,file="data/imputations/bc_Imp6.RData")
save(bc_Imp7,file="data/imputations/bc_Imp7.RData")
save(bc_Imp8,file="data/imputations/bc_Imp8.RData")
save(bc_Imp9,file="data/imputations/bc_Imp9.RData")
save(bc_Imp10,file="data/imputations/bc_Imp10.RData")
save(bc_Imp11,file="data/imputations/bc_Imp11.RData")
save(bc_Imp12,file="data/imputations/bc_Imp12.RData")
save(bc_Imp13,file="data/imputations/bc_Imp13.RData")
save(bc_Imp14,file="data/imputations/bc_Imp14.RData")






































