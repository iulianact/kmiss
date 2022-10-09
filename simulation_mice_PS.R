#----------------------------------------------------------------------------------
# mice1, mice2, mice3 and mice4 + PS 
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

numcores <- 2 

#----------------------------------------------------------------------------------
set.seed(331); impdat_mice_current <- impdat_mice1; time_mice_current <- time_mice1
# set.seed(332); impdat_mice_current <- impdat_mice2; time_mice_current <- time_mice2
# set.seed(333); impdat_mice_current <- impdat_mice3; time_mice_current <- time_mice3
# set.seed(334); impdat_mice_current <- impdat_mice4; time_mice_current <- time_mice4

mice_PS_boot_easy <- vector("list",nrep)
for(i in 1:nrep){
        impdat <- impdat_mice_current[[i]]
        timemice <- time_mice_current[[i]]
        
        miceACE2_est <- vector(mode = "list", length = numimp)
        miceACE2_estall <- vector(mode = "list", length = numimp)
        miceACE2_var <- vector(mode = "list", length = numimp)
        miceACE2_varall <- vector(mode = "list", length = numimp)
        for(j in 1:numimp){
                datmice <- cbind(time=timemice,complete(impdat,j))
                
                bootcall <- boot(data=datmice,
                                 statistic=easy_boot_PS,
                                 R=nboot,
                                 form_ps=form_ps,
                                 setTime=setTime,
                                 parallel = "multicore", ncpus=numcores)
                
                estdf <- data.frame(matrix(bootcall$t0,nrow=length(setTime),ncol=3))
                colnames(estdf) <- c("A=0","A=1","Diff")
                rownames(estdf) <- paste0("t.",setTime)
                
                vardf <- data.frame(matrix(apply(bootcall$t,2,var),nrow=length(setTime),ncol=3))
                colnames(vardf) <- c("A=0","A=1","Diff")
                rownames(vardf) <- paste0("t.",setTime)
                
                miceACE2_est[[j]] <- estdf[,3] 
                miceACE2_estall[[j]] <- estdf[,1:2] 
                miceACE2_var[[j]] <- vardf[,3]
                miceACE2_varall[[j]] <- vardf[,1:2]
        }
        
        ACEmice_PS_est <- Reduce("+", miceACE2_est)/numimp
        ACEmice_PS_estall <- Reduce("+", miceACE2_estall)/numimp
        
        between <- numeric(length(setTime))
        betweenall <- data.frame(matrix(0,nrow=length(setTime),ncol=2))
        colnames(betweenall) <- c("A=0","A=1")
        rownames(betweenall) <- paste0("t.",setTime)
        for(j in 1:numimp){
                between <- between + (miceACE2_est[[j]]-ACEmice_PS_est)^2
                betweenall <- betweenall + (miceACE2_estall[[j]]-ACEmice_PS_estall)^2
        }
        W <- Reduce("+", miceACE2_var)/numimp  
        B <- between/(numimp-1) 
        ACEmice_PS_var <- W + (1+1/numimp)*B 
        
        Wall <- Reduce("+", miceACE2_varall)/numimp  
        Ball <- betweenall/(numimp-1) 
        ACEmice_PS_varall <- Wall + (1+1/numimp)*Ball 
        
        mice_PS_boot_easy[[i]] <- list(est=ACEmice_PS_est,
                                        estall=ACEmice_PS_estall,
                                        var=ACEmice_PS_var,
                                        varall=ACEmice_PS_varall,
                                        fdeg=(numimp-1)*(1+numimp*W/((numimp+1)*B))^2,
                                        fdegall=(numimp-1)*(1+numimp*Wall/((numimp+1)*Ball))^2,
                                        B=B,W=W,Ball=Ball,Wall=Wall)
}
