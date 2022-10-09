#----------------------------------------------------------------------------------
# mice1, mice2, mice3 and mice4 imputation + standardisation
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
load(paste0(scen,"_impdat_mice1.RData"))
load(paste0(scen,"_impdat_mice2.RData"))
load(paste0(scen,"_impdat_mice3.RData"))
load(paste0(scen,"_impdat_mice4.RData"))

#models to fit
if(choose_misspec){form_std <- form_std_mis} else {form_std <- form_std_correct}
if(choose_misspec){form_ps <- form_ps_mis} else {form_ps <- form_ps_correct}

#----------------------------------------------------------------------------------
set.seed(21); impdat_mice_current <- impdat_mice1; time_mice_current <- time_mice1
# set.seed(22); impdat_mice_current <- impdat_mice2; time_mice_current <- time_mice2
# set.seed(23); impdat_mice_current <- impdat_mice3; time_mice_current <- time_mice3
# set.seed(24); impdat_mice_current <- impdat_mice4; time_mice_current <- time_mice4

mice_std_M <- vector("list",nrep)
for(i in 1:nrep){
  impdat <- impdat_mice_current[[i]]
  timemice <- time_mice_current[[i]]
  
  miceACE2_est <- vector(mode = "list", length = numimp)
  miceACE2_estall <- vector(mode = "list", length = numimp)
  miceACE2_var <- vector(mode = "list", length = numimp)
  miceACE2_varall <- vector(mode = "list", length = numimp)
  for(j in 1:numimp){
    tmpdat <- cbind(time=timemice,complete(impdat,j))
    tmpfit <- coxph(form_std, data=tmpdat, method = "breslow")
    fit.std_tmp <- stdCoxph(fit=tmpfit, data=tmpdat, X="A",t=setTime)
    miceACE2_est[[j]] <- Sdiff(fit.std_tmp$est)
    miceACE2_estall[[j]] <- fit.std_tmp$est
    miceACE2_var[[j]] <- sapply(1:length(setTime),FUN=function(x){t(c(-1,1)) %*% fit.std_tmp$vcov[[x]] %*% c(-1,1)})
    miceACE2_varall[[j]] <- data.frame(`A0`=sapply(fit.std_tmp$vcov,`[[`,1),`A1`=sapply(fit.std_tmp$vcov,`[[`,4))
  }
  
  ACEmice_std_est <- Reduce("+", miceACE2_est)/numimp
  ACEmice_std_estall <- Reduce("+", miceACE2_estall)/numimp
  
  between <- numeric(length(setTime))
  betweenall <- data.frame(matrix(0,nrow=length(setTime),ncol=2))
  colnames(betweenall) <- c("A=0","A=1")
  rownames(betweenall) <- paste0("t.",setTime)
  for(j in 1:numimp){
    between <- between + (miceACE2_est[[j]]-ACEmice_std_est)^2
    betweenall <- betweenall + (miceACE2_estall[[j]]-ACEmice_std_estall)^2
  }
  W <- Reduce("+", miceACE2_var)/numimp 
  B <- between/(numimp-1)
  ACEmice_std_var <- W + (1+1/numimp)*B 
  
  Wall <- Reduce("+", miceACE2_varall)/numimp  
  Ball <- betweenall/(numimp-1) 
  ACEmice_std_varall <- Wall + (1+1/numimp)*Ball 
  
  mice_std_M[[i]] <- list(est=ACEmice_std_est,
                           estall=ACEmice_std_estall,
                           var=ACEmice_std_var,
                           varall=ACEmice_std_varall,
                           fdeg=(numimp-1)*(1+numimp*W/((numimp+1)*B))^2,
                           fdegall=(numimp-1)*(1+numimp*Wall/((numimp+1)*Ball))^2,
                           B=B,W=W,Ball=Ball,Wall=Wall)
}
