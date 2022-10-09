#----------------------------------------------------------------------------------
# full, ignore, complete + PS 
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

numcores <- 2 

#----------------------------------------------------------------------------------
set.seed(1)

fullPS <- vector(mode="list",length=nrep)
for(i in 1:nrep){
        dat0 <- simdat[[i]]
        
        bootcall <- boot(data=dat0,
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
        
        fullPS[[i]] <- list(est=estdf[,3],
                            estall=estdf[,1:2],
                            var=vardf[,3],
                            varall=vardf[,1:2])
}

biasPS <- vector("list",nrep)
for(i in 1:nrep){
        dat0 <- simdat[[i]]
        
        bootcall <- boot(data=dat0,
                          statistic=easy_boot_PS,
                          R=nboot,
                          form_ps=formBias_ps,
                          setTime=setTime,
                          parallel = "multicore", ncpus=numcores)
        
        estdf <- data.frame(matrix(bootcall$t0,nrow=length(setTime),ncol=3))
        colnames(estdf) <- c("A=0","A=1","Diff")
        rownames(estdf) <- paste0("t.",setTime)
        
        vardf <- data.frame(matrix(apply(bootcall$t,2,var),nrow=length(setTime),ncol=3))
        colnames(vardf) <- c("A=0","A=1","Diff")
        rownames(vardf) <- paste0("t.",setTime)
        
        biasPS[[i]] <- list(est=estdf[,3],
                            estall=estdf[,1:2],
                            var=vardf[,3],
                            varall=vardf[,1:2])
}

complPS <- vector("list",nrep)
for(i in 1:nrep){
        datc <- simdat[[i]][simdat[[i]]$D>d0,]  
        
        bootcall <- boot(data=datc,
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
        
        complPS[[i]] <- list(est=estdf[,3],
                             estall=estdf[,1:2],
                             var=vardf[,3],
                             varall=vardf[,1:2])
}
