#----------------------------------------------------------------------------------
#  mice3 + other
#----------------------------------------------------------------------------------
rm(list=ls())

library(mice)
library(stdReg)
library(sandwich)
library(survival)
library(boot)
library(eventglm)

source("functions.R")
load("general_params.RData")
load("fit_params.RData")

#CHOOSE SCENARIO
choose_scenario <- 3 #1 or 3
choose_offset <- 0 #0 for censtype c
choose_censtype <- "c" #"a" or "c"
choose_misspec <- FALSE #TRUE or FALSE

choose_stratyear <- 5 #1 or 5
stratyear <- choose_stratyear

#load chosen scenario parameters
scen <- paste0("simdat",choose_scenario,choose_censtype,"_o",choose_offset)
load(paste0(scen,"_params.RData"))

#load data
load(paste0(scen,"_simdat.RData"))
load(paste0(scen,"_impdat_mice3.RData"))

#models to fit
if(choose_misspec){form_std <- form_std_mis} else {form_std <- form_std_correct}
if(choose_misspec){form_ps <- form_ps_mis} else {form_ps <- form_ps_correct}

#----------------------------------------------------------------------------------
#PSCox
set.seed(53)
mice3_PSCox <- vector("list",nrep)
for(i in 1:nrep){
  impdat <- impdat_mice3[[i]]
  timemice <- time_mice3[[i]]
  
  miceACE2_est <- vector(mode = "list", length = numimp)
  miceACE2_estall <- vector(mode = "list", length = numimp)
  for(j in 1:numimp){
    datmice <- cbind(time=timemice,complete(impdat,j))
    
    tmpdat <- datmice
    
    mfull_ipw <- glm(form_ps,family="binomial",data=tmpdat)
    predfull_ipw <- mfull_ipw$fitted.values
    wfull <- 1/(tmpdat$A*predfull_ipw+(1-tmpdat$A)*(1-predfull_ipw))
    
    tmpfit <- coxph(Surv(time,event) ~ A, data=tmpdat, method = "breslow", weights = wfull)
    tmpfit$weights <- rep(1,length(tmpfit$y))
    fit.std_tmp <- stdCoxph(fit=tmpfit, data=tmpdat, X="A",t=setTime)
    
    miceACE2_est[[j]] <- Sdiff(fit.std_tmp$est)
    miceACE2_estall[[j]] <- fit.std_tmp$est
  }
  
  ACEmice_PS_est <- Reduce("+", miceACE2_est)/numimp
  ACEmice_PS_estall <- Reduce("+", miceACE2_estall)/numimp
  
  mice3_PSCox[[i]] <- list(est=ACEmice_PS_est,
                           estall=ACEmice_PS_estall)
}

#----------------------------------------------------------------------------------
#PSCox-entry
set.seed(53)
mice3_PSCox_withD <- vector("list",nrep)
for(i in 1:nrep){
  impdat <- impdat_mice3[[i]]
  timemice <- time_mice3[[i]]
  
  miceACE2_est <- vector(mode = "list", length = numimp)
  miceACE2_estall <- vector(mode = "list", length = numimp)
  for(j in 1:numimp){
    datmice <- cbind(time=timemice,complete(impdat,j))
    
    tmpdat <- datmice
    
    mfull_ipw <- glm(form_ps,family="binomial",data=tmpdat)
    predfull_ipw <- mfull_ipw$fitted.values
    wfull <- 1/(tmpdat$A*predfull_ipw+(1-tmpdat$A)*(1-predfull_ipw))
    
    tmpfit <- coxph(Surv(time,event) ~ A + D, data=tmpdat, method = "breslow", weights = wfull)
    tmpfit$weights <- rep(1,length(tmpfit$y))
    fit.std_tmp <- stdCoxph(fit=tmpfit, data=tmpdat, X="A",t=setTime)
    
    miceACE2_est[[j]] <- Sdiff(fit.std_tmp$est)
    miceACE2_estall[[j]] <- fit.std_tmp$est
  }
  
  ACEmice_PS_est <- Reduce("+", miceACE2_est)/numimp
  ACEmice_PS_estall <- Reduce("+", miceACE2_estall)/numimp
  
  mice3_PSCox_withD[[i]] <- list(est=ACEmice_PS_est,
                           estall=ACEmice_PS_estall)
}

#----------------------------------------------------------------------------------
#PSCox-full
set.seed(53)
mice3_PSCox_full <- vector("list",nrep)
for(i in 1:nrep){
  impdat <- impdat_mice3[[i]]
  timemice <- time_mice3[[i]]
  
  miceACE2_est <- vector(mode = "list", length = numimp)
  miceACE2_estall <- vector(mode = "list", length = numimp)
  for(j in 1:numimp){
    datmice <- cbind(time=timemice,complete(impdat,j))
    
    tmpdat <- datmice
    
    mfull_ipw <- glm(form_ps,family="binomial",data=tmpdat)
    predfull_ipw <- mfull_ipw$fitted.values
    wfull <- 1/(tmpdat$A*predfull_ipw+(1-tmpdat$A)*(1-predfull_ipw))
    
    tmpfit <- coxph(form_std, data=tmpdat, method = "breslow", weights = wfull)
    tmpfit$weights <- rep(1,length(tmpfit$y))
    fit.std_tmp <- stdCoxph(fit=tmpfit, data=tmpdat, X="A",t=setTime)
    
    miceACE2_est[[j]] <- Sdiff(fit.std_tmp$est)
    miceACE2_estall[[j]] <- fit.std_tmp$est
  }
  
  ACEmice_PS_est <- Reduce("+", miceACE2_est)/numimp
  ACEmice_PS_estall <- Reduce("+", miceACE2_estall)/numimp
  
  mice3_PSCox_full[[i]] <- list(est=ACEmice_PS_est,
                           estall=ACEmice_PS_estall)
}


#----------------------------------------------------------------------------------
#PSstrat5 or PSstrat1
if(stratyear==5){setTime <- c(1,5)}

set.seed(413)
mice3_PSstrat_boot_easy <- vector("list",nrep)
for(i in 1:nrep){
  impdat <- impdat_mice3[[i]]
  timemice <- time_mice3[[i]]
  
  miceACE2_est <- vector(mode = "list", length = numimp)
  for(j in 1:numimp){
    datmice <- cbind(time=timemice,complete(impdat,j))
    tmpdat <- datmice
    
    if(stratyear==1){
      Dorig <- tmpdat$D
      Dstrat <- factor(floor(Dorig))
    }
    if(stratyear==5){
      Dorig <- tmpdat$D
      Dstrat <- numeric(length(Dorig))
      Dstrat[Dorig < 8] <- "1-7"
      Dstrat[Dorig >= 8 & Dorig < 13] <- "8-12"
      Dstrat[Dorig >= 13 & Dorig < 18] <- "13-17"
      Dstrat[Dorig >= 18 & Dorig < 23] <- "18-22"
      Dstrat[Dorig >= 23 & Dorig <= 28] <- "23-27"
      Dstrat <- factor(Dstrat,levels=c("1-7","8-12","13-17","18-22","23-27"))
    }
    
    miceACE_boot <- data.frame(matrix(NA,nrow=length(levels(Dstrat)),ncol=length(setTime)))
    for(d in 1:length(levels(Dstrat))){
      idd <- which(Dstrat==levels(Dstrat)[d])
      datidd <- tmpdat[idd,]
      
      mfull_ipw <- glm(form_ps,family="binomial",data=datidd)
      predfull_ipw <- mfull_ipw$fitted.values
      wfull <- 1/(datidd$A*predfull_ipw+(1-datidd$A)*(1-predfull_ipw))
      
      survfull_ipw_times <- survfit(Surv(time,event) ~ A, data=datidd, weights = wfull)
      ss_times <- summary(survfull_ipw_times, times=setTime)
      df_times <- data.frame(time=ss_times$time,A=ss_times$strata,surv=ss_times$surv)
      df_less1 <- df_times$surv[df_times$A=="A=1"]
      df_less0 <- df_times$surv[df_times$A=="A=0"]
      
      miceACE_boot[d,] <- df_less1 - df_less0
    }
    
    miceACE2_est[[j]] <- colMeans(miceACE_boot) 
  }
  
  ACEmice_PS_est <- Reduce("+", miceACE2_est)/numimp
  
  mice3_PSstrat_boot_easy[[i]] <- list(est=ACEmice_PS_est)
}

#----------------------------------------------------------------------------------
#PO
set.seed(60)
mice3_PO <- vector("list",nrep)
for(i in 1:nrep){
  impdat <- impdat_mice3[[i]]
  timemice <- time_mice3[[i]]
  
  miceACE2_est <- vector(mode = "list", length = numimp)
  warnall <- data.frame(matrix(NA,nrow=numimp,ncol=length(setTime)))
  for(j in 1:numimp){
    datmice <- cbind(time=timemice,complete(impdat,j))
    
    dA0 <- datmice; dA0$A <- 0
    dA1 <- datmice; dA1$A <- 1
    
    miceACE_boot <- numeric(length(setTime))
    for(k in 1:length(setTime)){
      tmpfit <- cumincglm(form_std,time=setTime[k],link="cloglog",data=datmice)
      warnall[j,k] <- 1-tmpfit$converged
      miceACE_boot[k] <- -(mean(predict(tmpfit,newdata = dA1,type="response")) - mean(predict(tmpfit,newdata = dA0,type="response"))) #since cloglog=log(-log(1-x))
    }
    
    miceACE2_est[[j]] <- miceACE_boot
  }
  
  ACEmice_PO_est <- Reduce("+", miceACE2_est)/numimp
  
  mice3_PO[[i]] <- list(est=ACEmice_PO_est,
                        warnall=warnall)
}

