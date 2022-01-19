#M.Kamenetsky
#2021-09-27
#Performs individal-level analysis with identified clusters
#This script is run on the server

rm(list=ls())
#############################
#load required packages
#############################
library(dplyr)
library(tidyr)
table.out <- NULL


#############################
#helper functions
#############################
'%ni%' <- Negate('%in%')

prepdata <- function(bcimpi, clusters_bicdf_uniqidnum){
  #clean up duplicate IDs (due to focal cell)
  clusters_bicdf_uniqidnum0 <- clusters_bicdf_uniqidnum %>%
    dplyr::select(-c(focalcell, X.1, X)) %>%
    group_by(idnum) %>%
    mutate(counter = n()) 
  
  
  dedups <- clusters_bicdf_uniqidnum0 %>%
    filter(counter!=1) %>%
    na_if(0) %>% 
    tidyr::fill(cluster100:cluster4375, .direction="downup") %>% 
    mutate_at(vars(cluster100:cluster4375), ~replace_na(.,0)) %>%
    distinct() %>%
    ungroup()
  
  singles <- clusters_bicdf_uniqidnum0 %>%
    filter(counter==1) 
  
  #combine back
  clusters_bicdf_uniqidnum_clean <- rbind.data.frame(singles, dedups)
  
  clusters_bicdf <-merge(bcimpi, clusters_bicdf_uniqidnum_clean[,c(1,30:44)], by.x="idnum",by.y="idnum", all.x=TRUE)
  #BIC
  clusters_bicdf_selectedvars <- clusters_bicdf %>%
    dplyr::select(lat, lon, bc, age, gave_birth,famhist, number_births, age_first_birth,menopause, age_menopause, race, educ, 
                  bmi, drinks_per_week,
                  hrt,
                  cluster100:cluster4375) 
  uniqcolnames <- attributes(unique(as.matrix(clusters_bicdf_selectedvars), MARGIN = 2))$dimnames[[2]]
  ix <- which(names(clusters_bicdf_selectedvars) %in% uniqcolnames)
  clusters_bicdf_selectedvars <- clusters_bicdf_selectedvars[, ix]
  
  clusters_bicdf_selectedvars_allppl <- clusters_bicdf_selectedvars %>%
    dplyr::select(lat, lon, bc, age, gave_birth,famhist, number_births, race, educ, bmi, drinks_per_week,
                  hrt,
                  age_first_birth,menopause, age_menopause,cluster100:cluster4375) %>%
    mutate_at(vars(cluster100:cluster4375), ~replace_na(.,0))
  clusters_bicdf_selectedvars_allppl$inacluster <-  if_else(rowSums(dplyr::select(clusters_bicdf_selectedvars_allppl,
                                                                                  starts_with("cluster"))) > 0, 1, 0)
  clusters_bicdf_selectedvars_allppl_complete <- clusters_bicdf_selectedvars_allppl
  return(clusters_bicdf_selectedvars_allppl_complete)
  
}

logitmodels.bic <- function(bcimpi, clusters_bicdf_uniqidnum){
    #######################################
    #prep data
    #######################################
    clusters_bicdf_selectedvars_allppl_complete <- prepdata(bcimpi, clusters_bicdf_uniqidnum)
    #######################################
    #Age-adjusted model
    #######################################
    print("AGE-ADJUSTED MODEL")
    clusters_bicdf_selectedvars_allppl_complete$age_cent <- (clusters_bicdf_selectedvars_allppl_complete$age-mean(clusters_bicdf_selectedvars_allppl_complete$age, na.rm=TRUE))/sd(clusters_bicdf_selectedvars_allppl_complete$age, na.rm=TRUE)
    clusters_bicdf_selectedvars_allppl_complete$age_menopause_cent <- (clusters_bicdf_selectedvars_allppl_complete$age_menopause-mean(clusters_bicdf_selectedvars_allppl_complete$age_menopause, na.rm=TRUE))/sd(clusters_bicdf_selectedvars_allppl_complete$age_menopause, na.rm=TRUE)
    clusters_bicdf_selectedvars_allppl_complete$age_first_birth_cent <- (clusters_bicdf_selectedvars_allppl_complete$age_first_birth-mean(clusters_bicdf_selectedvars_allppl_complete$age_first_birth, na.rm=TRUE))/sd(clusters_bicdf_selectedvars_allppl_complete$age_first_birth, na.rm=TRUE)
    clusters_bicdf_selectedvars_allppl_complete$bmi_cent <- (clusters_bicdf_selectedvars_allppl_complete$bmi - mean(clusters_bicdf_selectedvars_allppl_complete$bmi,na.rm=TRUE))/sd(clusters_bicdf_selectedvars_allppl_complete$bmi, na.rm=TRUE)
    
    m_age <- glm(bc ~ . -lat -lon -famhist -race -gave_birth  -number_births -hrt
                            - age - age_menopause_cent -age_first_birth_cent -bmi -bmi_cent -drinks_per_week -age_first_birth -educ -menopause -age_menopause -inacluster , 
                 data=clusters_bicdf_selectedvars_allppl_complete,
                 family="binomial")
    print(summary(m_age))
    #all ppl
    clusterMAT <- model.matrix(m_age)
    coefsvec <- coef(m_age)
    print(coefsvec)
    coefsvec[c(1,17)] <- 0
    print(coefsvec)
    vec <- (as.vector(clusterMAT%*% coefsvec))
    
    clusters_bicdf_selectedvars_allppl_complete$georisk_age <-  vec
    clusters_bicdf_selectedvars_allppl_complete$SEgeorisk_age <-  sqrt(diag(clusterMAT%*%vcov(m_age)%*%t(clusterMAT)))
  
    #######################################
    #Age, FamHist, Parity, AgeFirstBirth,AgeMenopause,Race, Education-adjusted model
    #######################################
     clusters_bicdf_selectedvars_allppl_complete$race <- relevel( clusters_bicdf_selectedvars_allppl_complete$race, ref="white")
     clusters_bicdf_selectedvars_allppl_complete$educ <- relevel( clusters_bicdf_selectedvars_allppl_complete$educ, ref="grade 12")
    print("AGE, FH, Parity, Age First Birth,agemeno, Race, Education-ADJUSTED MODEL")
    m_agefhnumbirthagefirstbagemenoraceeduc <- glm(bc ~ . + gave_birth:number_births +
                                                       gave_birth:age_first_birth_cent + menopause:age_menopause_cent +
                                                       race + educ - gave_birth - number_births - age_first_birth -menopause -age_menopause
                                                   - age_menopause_cent -age_first_birth_cent
                                                       -age  -bmi -bmi_cent -drinks_per_week -hrt,
                                                   data=subset(clusters_bicdf_selectedvars_allppl_complete,
                                                               select=-c(lat, lon,
                                                                         georisk_age,SEgeorisk_age,
                                                                         inacluster)),
                                                   family="binomial")

    print(summary(m_agefhnumbirthagefirstbagemenoraceeduc))
    #all ppl
    clusterMAT <- model.matrix(m_agefhnumbirthagefirstbagemenoraceeduc)
    coefsvec <- coef(m_agefhnumbirthagefirstbagemenoraceeduc)
    print(coefsvec)
    coefsvec[c(1:13, 29:32)] <-0
    print(coefsvec)
    vec <- (as.vector(clusterMAT%*% coefsvec))
    
    clusters_bicdf_selectedvars_allppl_complete$georisk_agefhnumbirthagefirstbagemenoraceeduc <-  vec
    clusters_bicdf_selectedvars_allppl_complete$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc <-  sqrt(diag(clusterMAT%*%vcov(m_agefhnumbirthagefirstbagemenoraceeduc)%*%t(clusterMAT)))
    
    #######################################
    #Age, FamHist, Parity, AgeFirstBirth,AgeMenopause,Race, Education, BMI, Drinks_per_week-adjusted model
    #######################################
    print("AGE, FH, Parity, Age First Birth,agemeno, Race, Education, BMI, Drinks_per_week-ADJUSTED MODEL")
    m_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- glm(bc ~ . + gave_birth:number_births +
                                                       gave_birth:age_first_birth_cent + menopause:age_menopause_cent +
                                                       race + bmi_cent + drinks_per_week + educ  
                                                      -age  - gave_birth - bmi - number_births - age_first_birth -menopause -age_menopause -hrt
                                                      -age_menopause_cent -age_first_birth_cent,
                                                   data=subset(clusters_bicdf_selectedvars_allppl_complete,
                                                               select=-c(lat, lon,
                                                                         georisk_age,SEgeorisk_age,
                                                                         georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                                         SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc, 
                                                                         inacluster)),
                                                   family="binomial")
    
    print(summary(m_agefhnumbirthagefirstbagemenoraceeducbmidrinks))
    #all ppl
    clusterMAT <- model.matrix(m_agefhnumbirthagefirstbagemenoraceeducbmidrinks)
    coefsvec <- coef(m_agefhnumbirthagefirstbagemenoraceeducbmidrinks)
    print(coefsvec)
    coefsvec[c(1:14, 30:34)] <-0
    print(coefsvec)
    vec <- (as.vector(clusterMAT%*% coefsvec))
    
    clusters_bicdf_selectedvars_allppl_complete$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks <-  vec
    clusters_bicdf_selectedvars_allppl_complete$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks <-  sqrt(diag(clusterMAT%*%vcov(m_agefhnumbirthagefirstbagemenoraceeducbmidrinks)%*%t(clusterMAT)))
    
    
    #######################################
    #Age, FamHist, Parity, AgeFirstBirth,AgeMenopause,Race, Education, BMI, Drinks_per_week, HRT
    #-adjusted model
    #######################################
    print("AGE, FH, Parity, Age First Birth,agemeno, Race, Education, BMI, Drinks_per_week, HRT-ADJUSTED MODEL")

    m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- glm(bc ~ . + gave_birth:number_births +
                                                              gave_birth:age_first_birth_cent + menopause:age_menopause_cent +
                                                              race + bmi_cent + drinks_per_week + educ + hrt 
                                                            -age - gave_birth - number_births - bmi - age_first_birth -menopause -age_menopause
                                                            -age_menopause_cent -age_first_birth_cent,
                                                            data=subset(clusters_bicdf_selectedvars_allppl_complete,
                                                                        select=-c(lat, lon,
                                                                                   georisk_age,SEgeorisk_age,
                                                                                  georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                                                  SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc, 
                                                                                  georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                                                  SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                                                  inacluster)),
                                                            family="binomial")
    
    print(summary(m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt))
    #all ppl
    clusterMAT <- model.matrix(m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt)
    coefsvec <- coef(m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt)
    print(coefsvec)
    coefsvec[c(1:16, 32:36)] <-0
    print(coefsvec)
    vec <- (as.vector(clusterMAT%*% coefsvec))
    clusters_bicdf_selectedvars_allppl_complete$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <-  vec
    clusters_bicdf_selectedvars_allppl_complete$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <-  sqrt(diag(clusterMAT%*%vcov(m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt)%*%t(clusterMAT)))
    
    #return everything

    return(list(clusters_bicdf_selectedvars_allppl_complete = clusters_bicdf_selectedvars_allppl_complete,
                m_age = m_age,
                m_agefhnumbirthagefirstbagemenoraceeduc = m_agefhnumbirthagefirstbagemenoraceeduc,
                m_agefhnumbirthagefirstbagemenoraceeducbmidrinks = m_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt = m_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt))
}

############################################################################################
#Load in csvs to process
############################################################################################
#####################
#BIC
#####################
#load imputation 1-5
load("data/imputations/bc_Imp1.RData")
load("data/imputations/bc_Imp2.RData")
load("data/imputations/bc_Imp3.RData")
load("data/imputations/bc_Imp4.RData")
load("data/imputations/bc_Imp5.RData")
load("data/imputations/bc_Imp6.RData")
load("data/imputations/bc_Imp7.RData")
load("data/imputations/bc_Imp8.RData")
load("data/imputations/bc_Imp9.RData")
load("data/imputations/bc_Imp10.RData")
load("data/imputations/bc_Imp11.RData")
load("data/imputations/bc_Imp12.RData")
load("data/imputations/bc_Imp13.RData")
load("data/imputations/bc_Imp14.RData")

#load clusters df
clusters_bicdf_uniqidnum<-  read.csv("data/clusters_bicdf_uniqidnum.csv")

resbic.imp1 <- logitmodels.bic(bc_Imp1, clusters_bicdf_uniqidnum)
resbic.imp2 <- logitmodels.bic(bc_Imp2, clusters_bicdf_uniqidnum)
resbic.imp3 <- logitmodels.bic(bc_Imp3, clusters_bicdf_uniqidnum)
resbic.imp4 <- logitmodels.bic(bc_Imp4, clusters_bicdf_uniqidnum)
resbic.imp5 <- logitmodels.bic(bc_Imp5, clusters_bicdf_uniqidnum)
resbic.imp6 <- logitmodels.bic(bc_Imp6, clusters_bicdf_uniqidnum)
resbic.imp7 <- logitmodels.bic(bc_Imp7, clusters_bicdf_uniqidnum)
resbic.imp8 <- logitmodels.bic(bc_Imp8, clusters_bicdf_uniqidnum)
resbic.imp9 <- logitmodels.bic(bc_Imp9, clusters_bicdf_uniqidnum)
resbic.imp10 <- logitmodels.bic(bc_Imp10, clusters_bicdf_uniqidnum)
resbic.imp11 <- logitmodels.bic(bc_Imp11, clusters_bicdf_uniqidnum)
resbic.imp12 <- logitmodels.bic(bc_Imp12, clusters_bicdf_uniqidnum)
resbic.imp13 <- logitmodels.bic(bc_Imp13, clusters_bicdf_uniqidnum)
resbic.imp14 <- logitmodels.bic(bc_Imp14, clusters_bicdf_uniqidnum)

############################################################
############################################################
############################################################
#COMBINE IMPUTATION ESTIMATES
############################################################
############################################################
############################################################
####################################
#Age-adjusted
####################################
############
#Coef estimates & SEs
############
coefests_age <- cbind(coef(resbic.imp1[[2]]),
                              coef(resbic.imp2[[2]]),
                              coef(resbic.imp3[[2]]),
                              coef(resbic.imp4[[2]]),
                              coef(resbic.imp5[[2]]),
                              coef(resbic.imp6[[2]]),
                              coef(resbic.imp7[[2]]),
                              coef(resbic.imp8[[2]]),
                              coef(resbic.imp9[[2]]),
                              coef(resbic.imp10[[2]]),
                              coef(resbic.imp11[[2]]),
                              coef(resbic.imp12[[2]]),
                              coef(resbic.imp13[[2]]),
                              coef(resbic.imp14[[2]]))
coefest_age <- apply(coefests_age,1,mean)

sescoef_age <- cbind(coef(summary(resbic.imp1[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp2[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp3[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp4[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp5[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp6[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp7[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp8[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp9[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp10[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp11[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp12[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp13[[2]]))[, "Std. Error"],
                            coef(summary(resbic.imp14[[2]]))[, "Std. Error"])
secoef_age <- sqrt(apply(coefests_age,1,var)+(1+1/14)*apply(sescoef_age,1,function(x) mean(x^2)))
#save results to RData file
save(coefest_age,secoef_age , file="mi_coefse_age.RData")

############
#per person
############
ests.georisk_age <- cbind(resbic.imp1[[1]]$georisk_age,
                                             resbic.imp2[[1]]$georisk_age,
                                             resbic.imp3[[1]]$georisk_age,
                                             resbic.imp4[[1]]$georisk_age,
                                             resbic.imp5[[1]]$georisk_age,
                                             resbic.imp6[[1]]$georisk_age,
                                             resbic.imp7[[1]]$georisk_age,
                                             resbic.imp8[[1]]$georisk_age,
                                             resbic.imp9[[1]]$georisk_age,
                                             resbic.imp10[[1]]$georisk_age,
                                             resbic.imp11[[1]]$georisk_age,
                                             resbic.imp12[[1]]$georisk_age,
                                             resbic.imp13[[1]]$georisk_age,
                                             resbic.imp14[[1]]$georisk_age)
est.georisk_age <- apply(ests.georisk_age,1,mean)

se.ests <- cbind(resbic.imp1[[1]]$SEgeorisk_age,
                 resbic.imp2[[1]]$SEgeorisk_age,
                 resbic.imp3[[1]]$SEgeorisk_age,
                 resbic.imp4[[1]]$SEgeorisk_age,
                 resbic.imp5[[1]]$SEgeorisk_age,
                 resbic.imp6[[1]]$SEgeorisk_age,
                 resbic.imp7[[1]]$SEgeorisk_age,
                 resbic.imp8[[1]]$SEgeorisk_age,
                 resbic.imp9[[1]]$SEgeorisk_age,
                 resbic.imp10[[1]]$SEgeorisk_age,
                 resbic.imp11[[1]]$SEgeorisk_age,
                 resbic.imp12[[1]]$SEgeorisk_age,
                 resbic.imp13[[1]]$SEgeorisk_age,
                 resbic.imp14[[1]]$SEgeorisk_age)
se.georisk_age <- sqrt(apply(ests.georisk_age,1,var)+(1+1/14)*apply(se.ests,1,function(x) mean(x^2)))


####################################
#Age, family history, parity, age first birth, age menopause,race, education-adjusted
####################################
############
#Coef estimates & SEs
############
coefests_agefhnumbirthagefirstbagemenoraceeduc <- cbind(coef(resbic.imp1[[3]]),
                                                 coef(resbic.imp2[[3]]),
                                                 coef(resbic.imp3[[3]]),
                                                 coef(resbic.imp4[[3]]),
                                                 coef(resbic.imp5[[3]]),
                                                 coef(resbic.imp6[[3]]),
                                                 coef(resbic.imp7[[3]]),
                                                 coef(resbic.imp8[[3]]),
                                                 coef(resbic.imp9[[3]]),
                                                 coef(resbic.imp10[[3]]),
                                                 coef(resbic.imp11[[3]]),
                                                 coef(resbic.imp12[[3]]),
                                                 coef(resbic.imp13[[3]]),
                                                 coef(resbic.imp14[[3]]))
coefest_agefhnumbirthagefirstbagemenoraceeduc <- apply(coefests_agefhnumbirthagefirstbagemenoraceeduc,1,mean)

sescoef_agefhnumbirthagefirstbagemenoraceeduc <- cbind(coef(summary(resbic.imp1[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp2[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp3[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp4[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp5[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp6[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp7[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp8[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp9[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp10[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp11[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp12[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp13[[3]]))[, "Std. Error"],
                                                coef(summary(resbic.imp14[[3]]))[, "Std. Error"])
secoef_agefhnumbirthagefirstbagemenoraceeduc <- sqrt(apply(coefests_agefhnumbirthagefirstbagemenoraceeduc,1,var)+(1+1/14)*apply(sescoef_agefhnumbirthagefirstbagemenoraceeduc,1,function(x) mean(x^2)))
#save results to RData file
save(coefest_agefhnumbirthagefirstbagemenoraceeduc,secoef_agefhnumbirthagefirstbagemenoraceeduc , file="mi_coefse_agefhnumbirthagefirstbagemenoraceeduc.RData")

############
#per person
############
ests.georisk_agefhnumbirthagefirstbagemenoraceeduc <- cbind(resbic.imp1[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp2[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp3[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp4[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp5[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp6[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp7[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp8[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp9[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp10[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp11[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp12[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp13[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc,
                                                     resbic.imp14[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeduc)
est.georisk_agefhnumbirthagefirstbagemenoraceeduc <- apply(ests.georisk_agefhnumbirthagefirstbagemenoraceeduc,1,mean)

se.ests <- cbind(resbic.imp1[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp2[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp3[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp4[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp5[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp6[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp7[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp8[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp9[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp10[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp11[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp12[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp13[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc,
                 resbic.imp14[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeduc)
se.georisk_agefhnumbirthagefirstbagemenoraceeduc <- sqrt(apply(ests.georisk_agefhnumbirthagefirstbagemenoraceeduc,1,var)+(1+1/14)*apply(se.ests,1,function(x) mean(x^2)))


####################################
#Age, family history, parity, age first birth, age menopause,race, education, bmi, drinks_per_week-adjusted
####################################
############
#Coef estimates & SEs
############
coefests_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- cbind(coef(resbic.imp1[[4]]),
                                                        coef(resbic.imp2[[4]]),
                                                        coef(resbic.imp3[[4]]),
                                                        coef(resbic.imp4[[4]]),
                                                        coef(resbic.imp5[[4]]),
                                                        coef(resbic.imp6[[4]]),
                                                        coef(resbic.imp7[[4]]),
                                                        coef(resbic.imp8[[4]]),
                                                        coef(resbic.imp9[[4]]),
                                                        coef(resbic.imp10[[4]]),
                                                        coef(resbic.imp11[[4]]),
                                                        coef(resbic.imp12[[4]]),
                                                        coef(resbic.imp13[[4]]),
                                                        coef(resbic.imp14[[4]]))
coefest_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- apply(coefests_agefhnumbirthagefirstbagemenoraceeducbmidrinks,1,mean)

sescoef_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- cbind(coef(summary(resbic.imp1[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp2[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp3[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp4[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp5[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp6[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp7[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp8[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp9[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp10[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp11[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp12[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp13[[4]]))[, "Std. Error"],
                                                       coef(summary(resbic.imp14[[4]]))[, "Std. Error"])
secoef_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- sqrt(apply(coefests_agefhnumbirthagefirstbagemenoraceeducbmidrinks,1,var)+(1+1/14)*apply(sescoef_agefhnumbirthagefirstbagemenoraceeducbmidrinks,1,function(x) mean(x^2)))
#save results to RData file
save(coefest_agefhnumbirthagefirstbagemenoraceeducbmidrinks,secoef_agefhnumbirthagefirstbagemenoraceeducbmidrinks , file="mi_coefse_agefhnumbirthagefirstbagemenoraceeducbmidrinks.RData")

############
#per person
############
ests.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- cbind(resbic.imp1[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp2[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp3[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp4[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp5[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp6[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp7[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp8[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp9[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp10[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp11[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp12[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp13[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                                                            resbic.imp14[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks)
est.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- apply(ests.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,1,mean)

se.ests <- cbind(resbic.imp1[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp2[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp3[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp4[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp5[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp6[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp7[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp8[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp9[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp10[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp11[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp12[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp13[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,
                 resbic.imp14[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks)
se.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- sqrt(apply(ests.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks,1,var)+(1+1/14)*apply(se.ests,1,function(x) mean(x^2)))

####################################
#Age, family history, parity, age first birth, age menopause,race, education, bmi, drinks_per_week, HRT-adjusted
####################################
############
#Coef estimates & SEs
############
coefests_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- cbind(coef(resbic.imp1[[5]]),
                                                                 coef(resbic.imp2[[5]]),
                                                                 coef(resbic.imp3[[5]]),
                                                                 coef(resbic.imp4[[5]]),
                                                                 coef(resbic.imp5[[5]]),
                                                                 coef(resbic.imp6[[5]]),
                                                                 coef(resbic.imp7[[5]]),
                                                                 coef(resbic.imp8[[5]]),
                                                                 coef(resbic.imp9[[5]]),
                                                                 coef(resbic.imp10[[5]]),
                                                                 coef(resbic.imp11[[5]]),
                                                                 coef(resbic.imp12[[5]]),
                                                                 coef(resbic.imp13[[5]]),
                                                                 coef(resbic.imp14[[5]]))
coefest_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- apply(coefests_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,1,mean)

sescoef_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- cbind(coef(summary(resbic.imp1[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp2[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp3[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp4[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp5[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp6[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp7[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp8[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp9[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp10[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp11[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp12[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp13[[5]]))[, "Std. Error"],
                                                                coef(summary(resbic.imp14[[5]]))[, "Std. Error"])
secoef_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- sqrt(apply(coefests_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,1,var)+(1+1/14)*apply(sescoef_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,1,function(x) mean(x^2)))
#save results to RData file
save(coefest_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,secoef_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt , file="mi_coefse_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt.RData")

############
#per person
############
ests.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- cbind(resbic.imp1[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp2[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp3[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp4[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp5[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp6[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp7[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp8[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp9[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp10[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp11[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp12[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp13[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                                                                     resbic.imp14[[1]]$georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt)
est.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- apply(ests.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,1,mean)

se.ests <- cbind(resbic.imp1[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp2[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp3[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp4[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp5[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp6[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp7[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp8[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp9[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp10[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp11[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp12[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp13[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,
                 resbic.imp14[[1]]$SEgeorisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt)
se.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- sqrt(apply(ests.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt,1,var)+(1+1/14)*apply(se.ests,1,function(x) mean(x^2)))



####################################
#add mi averaged estimates and standard errors to dataframe
####################################
dat <- resbic.imp1[[1]]
dat$inacluster <-  if_else(rowSums(dplyr::select(dat,starts_with("cluster"))) > 0, 1, 0)

#add in age
dat$miest_age<- exp(est.georisk_age)
dat$miSEest_age<- se.georisk_age
#add in age, family history, parity, age first birth, agemenopause, race, education
dat$miest_agefhnumbirthagefirstbagemenoraceeduc  <- exp(est.georisk_agefhnumbirthagefirstbagemenoraceeduc)
dat$miSEest_agefhnumbirthagefirstbagemenoraceeduc <- se.georisk_agefhnumbirthagefirstbagemenoraceeduc

#add in age, family history, parity, age first birth, agemenopause, race, education, bmi, drinksperweek
dat$miest_agefhnumbirthagefirstbagemenoraceeducbmidrinks  <- exp(est.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks)
dat$miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinks <- se.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinks

#add in age, family history, parity, age first birth, agemenopause, race, education, bmi, drinksperweek, hrt
dat$miest_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt  <- exp(est.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt)
dat$miSEest_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt <- se.georisk_agefhnumbirthagefirstbagemenoraceeducbmidrinkshrt


#save
write.csv(dat, file="dat_pooledmiresults_20210517.csv")




