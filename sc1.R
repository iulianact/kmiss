#----------------------------------------------------------------------------------
# Scenario I
#----------------------------------------------------------------------------------
rm(list=ls())
set.seed(42)

admcens <- TRUE
offset <- 0 

cCens <- c(cCens_D = 0) 
lambdaC <- 1/15 

censtype <- ifelse(admcens,"a",ifelse(cCens==0,"c","b"))
offset <- ifelse(censtype=="a",offset,0)
scen <- paste0("simdat1",censtype)

scen <- paste0(scen,"_o",offset)

#----------------------------------------------------------------------------------
N <- 2*10^4

maxD <- 27; minD <- 1

cH <- c(offset=15,alpha=3,beta=1,scale=60) 

cC <- c(cC_0=-6.5, cC_H=0.1, cC_D=0, cC_D2=0) 

cA <- c(cA_0=3,  
        cA_H=-0.12, 
        cA_D=0.01,  
        cA_C=-0.5, 
        cA_DC=0.05, 
        cA_DH=0.002) 

cT <- c(cT_H=0.005, 
        cT_D=0, 
        cT_C=2,
        
        cT_DC=0,
        cT_DH=0, 
        cT_CH=0,  
        
        cT_A=-1.2,
        cT_AD=0, 
        cT_AH=0, 
        cT_AC=0, 
        cT_ACD=0)  

lambdaT <- 1/300 
nuT <- 2.5

d0 <- 7

#----------------------------------------------------------------------------------
Dcont <- runif(N,min=minD,max=maxD+1)
D <- Dcont 
H <- cH["offset"] + rbeta(N, shape1=cH["alpha"],shape2=cH["beta"]) * cH["scale"]
C <- rbinom(N,1,plogis(cC["cC_0"] + cC["cC_H"]*H + cC["cC_D"]*D + cC["cC_D2"]*D^2))
A <- rbinom(N,1,plogis(cA["cA_0"] + cA["cA_H"]*H + cA["cA_D"]*D + cA["cA_C"]*C + cA["cA_DC"]* D*C +cA["cA_DH"]*D*H))

linT <- cT["cT_H"]*H + cT["cT_D"]*D + cT["cT_C"]*C + cT["cT_A"]*A + cT["cT_DC"]*D*C + cT["cT_DH"]*D*H + cT["cT_CH"]*C*H + cT["cT_AD"]*A*D + cT["cT_AH"]*A*H + cT["cT_AC"]*A*C + cT["cT_ACD"]*A*C*D
timeT <- (-log(runif(N,min=0,max=1))/(lambdaT*exp(linT)))^(1/nuT)

timeC <- maxD+1+offset-Dcont
time <- pmin(timeT,timeC)
event <- ifelse(timeT<timeC,1,0)

linC <- cCens["cCens_D"]*D
timeCR <- (-log(runif(N,min=0,max=1))/(lambdaC*exp(linC)))
timeR <- pmin(timeT,timeCR)
eventR <- ifelse(timeT<timeCR,1,0)

if(admcens){timedat <- time} else {timedat <- timeR}
if(admcens){eventdat <- event} else {eventdat <- eventR}
dat0 <- data.frame(time=timedat,event=eventdat,A=A,C=C,D=D,H=H)

# save(N,maxD,minD,cH,cC,cA,cT,
#      lambdaT,nuT,
#      cCens,lambdaC,
#      admcens,offset,d0,
#      file=paste0(scen,"_params.RData"))
