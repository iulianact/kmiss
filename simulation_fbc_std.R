#----------------------------------------------------------------------------------
# full, ignore, compl + standardisation 
#----------------------------------------------------------------------------------
rm(list=ls())

library(mice)
library(stdReg)
library(sandwich)
library(survival)
library(boot)

source("functions.R")
load("general_params.RData")
load("fit_params.RData")

#CHOOSE SCENARIO
choose_scenario <- 1 #1 or 3
choose_offset <- 0 #0 for censtype c
choose_censtype <- "a" #"a" or "c"
choose_misspec <- FALSE #TRUE or FALSE

#load chosen scenario parameters
scen <- paste0("simdat",choose_scenario,choose_censtype,"_o",choose_offset)
load(paste0(scen,"_params.RData"))

#load data
load(paste0(scen,"_simdat.RData"))

#models to fit
if(choose_misspec){form_std <- form_std_mis} else {form_std <- form_std_correct}
if(choose_misspec){form_ps <- form_ps_mis} else {form_ps <- form_ps_correct}

#----------------------------------------------------------------------------------
set.seed(1)

full <- vector(mode="list",length=nrep)
for(i in 1:nrep){
        dat0 <- simdat[[i]]
        mfull <- coxph(form_std, dat=dat0, method="breslow")
        fit.std_full <- stdCoxph(fit=mfull, data=dat0, X="A",t=setTime)
        ACEfull_std_est <- Sdiff(fit.std_full$est)
        ACEfull_std_var <- sapply(1:length(setTime),FUN=function(x){t(c(-1,1)) %*% fit.std_full$vcov[[x]] %*% c(-1,1)})
        
        print(paste("full",i))
        
        full[[i]] <- list(est=ACEfull_std_est,
                          var=ACEfull_std_var,
                          estall=fit.std_full$est,
                          varall=data.frame(`A0`=sapply(fit.std_full$vcov,`[[`,1),`A1`=sapply(fit.std_full$vcov,`[[`,4)))
}

bias <- vector("list",nrep)
for(i in 1:nrep){
        dat0 <- simdat[[i]]
        mbias <- coxph(formBias_std, dat=dat0, method="breslow")
        fit.std_bias <- stdCoxph(fit=mbias, data=dat0, X="A",t=setTime)
        ACEbias_std_est <- Sdiff(fit.std_bias$est)
        ACEbias_std_var <- sapply(1:length(setTime),FUN=function(x){t(c(-1,1)) %*% fit.std_bias$vcov[[x]] %*% c(-1,1)})
        
        bias[[i]] <- list(est=ACEbias_std_est,
                          var=ACEbias_std_var,
                          estall=fit.std_bias$est,
                          varall=data.frame(`A0`=sapply(fit.std_bias$vcov,`[[`,1),`A1`=sapply(fit.std_bias$vcov,`[[`,4)))
}

compl <- vector("list",nrep)
for(i in 1:nrep){
        datc <- simdat[[i]][simdat[[i]]$D>d0,]  
        mcompl <- coxph(form_std, dat=datc, method="breslow")
        fit.std_compl <- stdCoxph(fit=mcompl, data=datc, X="A",t=setTime)
        ACEcompl_std_est <- Sdiff(fit.std_compl$est)
        ACEcompl_std_var <- sapply(1:length(setTime),FUN=function(x){t(c(-1,1)) %*% fit.std_compl$vcov[[x]] %*% c(-1,1)})
        
        compl[[i]] <- list(est=ACEcompl_std_est,
                           var=ACEcompl_std_var,
                           estall=fit.std_compl$est,
                           varall=data.frame(`A0`=sapply(fit.std_compl$vcov,`[[`,1),`A1`=sapply(fit.std_compl$vcov,`[[`,4)))
}
