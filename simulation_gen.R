#----------------------------------------------------------------------------------
# Generate full and imputed data
#----------------------------------------------------------------------------------
rm(list=ls())

library(mice)
library(stdReg)
library(sandwich)
library(survival)
library(boot)

source("functions.R")
load("general_params.RData")

#CHOOSE SCENARIO
choose_scenario <- 1 #1 or 3
choose_offset <- 0 #0 for censtype c
choose_censtype <- "a" #"a" or "c"

#load chosen scenario parameters
scen <- paste0("simdat",choose_scenario,choose_censtype,"_o",choose_offset)
load(paste0(scen,"_params.RData"))

#----------------------------------------------------------------------------------
# Data
#----------------------------------------------------------------------------------
set.seed(42)

simdat <- lapply(1:nrep, FUN=function(x){generate_surv(N=N,d0=d0,
                                                       maxD=maxD,minD=minD,
                                                       cH=cH,cC=cC,cA=cA,
                                                       cT=cT,lambdaT=lambdaT,nuT=nuT,
                                                       cCens=cCens,lambdaC=lambdaC,admcens=admcens,offset=offset)})

rATE <- getreal(mrep=1000,M=10^6,
                maxD=maxD,minD=minD,
                cH=cH,cC=cC,cA=cA,
                cT=cT,lambdaT=lambdaT,nuT=nuT,
                setTime=setTime)

params <- list(N=N,d0=d0,
               maxD=maxD,minD=minD,
               cH=cH,cC=cC,cA=cA,
               cT=cT,lambdaT=lambdaT,nuT=nuT,
               cCens=cCens,lambdaC=lambdaC,admcens=admcens,offset=offset)

# save(simdat,params,rATE,file=paste0(scen,"_simdat.RData"))

#----------------------------------------------------------------------------------
# Impute and save mice1
#----------------------------------------------------------------------------------
set.seed(91)

impdat_mice1 <- list("vector",nrep)
time_mice1 <- list("vector",nrep)
for(i in 1:nrep){
        datm <- simdat[[i]]
        datm$C[datm$D<=d0] <- NA
        datmice <- datm
        datmice$C <- factor(datmice$C)
        datmice$logtime <- log(datmice$time)
        
        timemice <- datmice$time
        datmice$time <- NULL
        
        impdat_mice1[[i]] <- mice(datmice, m=numimp, maxit = 5, method = 'logreg', print=FALSE)
        time_mice1[[i]] <- timemice
}

# save(impdat_mice1,time_mice1,file=paste0(scen,"_impdat_mice1.RData"))

#----------------------------------------------------------------------------------
# Impute and save mice2
#----------------------------------------------------------------------------------
set.seed(92)

impdat_mice2 <- list("vector",nrep)
time_mice2 <- list("vector",nrep)
for(i in 1:nrep){
        datm <- simdat[[i]]
        datm$C[datm$D<=d0] <- NA
        datmice <- datm
        datmice$C <- factor(datmice$C)
        
        sfit <- survfit(Surv(datmice$time,datmice$event)~1,se.fit=FALSE)
        sumsurv <- summary(sfit,times=datmice$time)
        datmice$HT <- sumsurv$cumhaz[match(datmice$time,sumsurv$time)]
        
        timemice <- datmice$time
        datmice$time <- NULL
        
        impdat_mice2[[i]] <- mice(datmice, m=numimp, maxit = 5, method = 'logreg', print=FALSE)
        time_mice2[[i]] <- timemice
}

# save(impdat_mice2,time_mice2,file=paste0(scen,"_impdat_mice2.RData"))

#----------------------------------------------------------------------------------
# Impute and save mice3 (A, D, H, event indicator, H_0(T), H_0(T) x {A,D,H})
#----------------------------------------------------------------------------------
set.seed(93)

impdat_mice3 <- list("vector",nrep)
time_mice3 <- list("vector",nrep)
for(i in 1:nrep){
        datm <- simdat[[i]]
        datm$C[datm$D<=d0] <- NA
        datmice <- datm
        datmice$C <- factor(datmice$C)
        
        sfit <- survfit(Surv(datmice$time,datmice$event)~1,se.fit=FALSE)
        sumsurv <- summary(sfit,times=datmice$time)
        datmice$HT <- sumsurv$cumhaz[match(datmice$time,sumsurv$time)]
        
        datmice$HTD <- datmice$D * datmice$HT
        datmice$HTH <- datmice$H * datmice$HT
        datmice$HTA <- datmice$A * datmice$HT
        timemice <- datmice$time
        datmice$time <- NULL

        impdat_mice3[[i]] <- mice(datmice, m=numimp, maxit = 5, method = 'logreg', print=FALSE)
        time_mice3[[i]] <- timemice
}

# save(impdat_mice3,time_mice3,file=paste0(scen,"_impdat_mice3.RData"))

#----------------------------------------------------------------------------------
# Impute and save mice4 (A, D, H, event indicator, logT, logT x {A,D,H})
#----------------------------------------------------------------------------------
set.seed(94)

impdat_mice4 <- list("vector",nrep)
time_mice4 <- list("vector",nrep)
for(i in 1:nrep){
        datm <- simdat[[i]]
        datm$C[datm$D<=d0] <- NA
        datmice <- datm
        datmice$C <- factor(datmice$C)
        datmice$logtime <- log(datmice$time)
        
        datmice$logTD <- datmice$D * datmice$logtime
        datmice$logTH <- datmice$H * datmice$logtime
        datmice$logTA <- datmice$A * datmice$logtime
        timemice <- datmice$time
        datmice$time <- NULL
        
        impdat_mice4[[i]] <- mice(datmice, m=numimp, maxit = 5, method = 'logreg',print=FALSE)
        time_mice4[[i]] <- timemice
}

# save(impdat_mice4,time_mice4,file=paste0(scen,"_impdat_mice4.RData"))
