generate_surv <- function(N,
                          d0,
                          maxD,minD,
                          cH,cC,cA,
                          cT,lambdaT,nuT,
                          cCens,lambdaC,admcens,offset){
  Dcont <- runif(N,min=minD,max=maxD+1)
  D <- Dcont 
  H <- cH["offset"] + rbeta(N, shape1=cH["alpha"],shape2=cH["beta"]) * cH["scale"]
  C <- rbinom(N,1,plogis(cC["cC_0"] + cC["cC_H"]*H + cC["cC_D"]*D + cC["cC_D2"]*D^2))
  A <- rbinom(N,1,plogis(cA["cA_0"] + cA["cA_H"]*H + cA["cA_D"]*D + cA["cA_C"]*C + cA["cA_DC"]* D*C +cA["cA_DH"]*D*H))
  
  linT <- cT["cT_H"]*H + cT["cT_D"]*D + cT["cT_C"]*C + cT["cT_A"]*A + cT["cT_DC"]*D*C + cT["cT_DH"]*D*H + cT["cT_CH"]*C*H + cT["cT_AD"]*A*D + cT["cT_AH"]*A*H + cT["cT_AC"]*A*C + cT["cT_ACD"]*A*C*D
  timeT <- (-log(runif(N,min=0,max=1))/(lambdaT*exp(linT)))^(1/nuT)
  
  if(admcens){
    timeC <- maxD+1+offset-Dcont
  } else {
    linC <- cCens["cCens_D"]*D
    timeC <- (-log(runif(N,min=0,max=1))/(lambdaC*exp(linC)))
  }
  
  time <- pmin(timeT,timeC)
  event <- ifelse(timeT<timeC,1,0)
  
  dat0 <- data.frame(time=time,event=event,A=A,C=C,D=D,H=H)
  
  return(dat0)
}

#----------------------------------------------------------------------------------
getreal <- function(mrep,M,
                    maxD,minD,
                    cH,cC,cA,
                    cT,lambdaT,nuT,
                    setTime){
  
  rS <- data.frame(matrix(NA,nrow=mrep,ncol=length(setTime)))
  rS1 <- data.frame(matrix(NA,nrow=mrep,ncol=length(setTime)))
  rS0 <- data.frame(matrix(NA,nrow=mrep,ncol=length(setTime)))
  
  for(k in 1:mrep){
    Dcont <- runif(M,min=minD,max=maxD+1)
    D <- Dcont 
    H <- cH["offset"] + rbeta(M, shape1=cH["alpha"],shape2=cH["beta"]) * cH["scale"]
    C <- rbinom(M,1,plogis(cC["cC_0"] + cC["cC_H"]*H + cC["cC_D"]*D + cC["cC_D2"]*D^2))
    A <- rbinom(M,1,plogis(cA["cA_0"] + cA["cA_H"]*H + cA["cA_D"]*D + cA["cA_C"]*C + cA["cA_DC"]* D*C +cA["cA_DH"]*D*H))
    
    linT1 <- cT["cT_H"]*H + cT["cT_D"]*D + cT["cT_C"]*C + cT["cT_A"]*1 + cT["cT_DC"]*D*C + cT["cT_DH"]*D*H + cT["cT_CH"]*C*H + cT["cT_AD"]*1*D + cT["cT_AH"]*1*H + cT["cT_AC"]*1*C + cT["cT_ACD"]*1*C*D
    timeT1 <- (-log(runif(M,min=0,max=1))/(lambdaT*exp(linT1)))^(1/nuT)
    linT0 <- cT["cT_H"]*H + cT["cT_D"]*D + cT["cT_C"]*C + cT["cT_A"]*0 + cT["cT_DC"]*D*C + cT["cT_DH"]*D*H + cT["cT_CH"]*C*H + cT["cT_AD"]*0*D + cT["cT_AH"]*0*H + cT["cT_AC"]*0*C + cT["cT_ACD"]*0*C*D
    timeT0 <- (-log(runif(M,min=0,max=1))/(lambdaT*exp(linT0)))^(1/nuT)

    S1 <- data.frame(matrix(NA,nrow=M,ncol=length(setTime)))
    S0 <- data.frame(matrix(NA,nrow=M,ncol=length(setTime)))
    for(i in 1:length(setTime)){
      S1[,i] <- timeT1 > setTime[i]
      S0[,i] <- timeT0 > setTime[i]
    }
    rS[k,] <- colMeans(S1) - colMeans(S0)
    rS1[k,] <- colMeans(S1)
    rS0[k,] <- colMeans(S0)

  }
  rATE <- colMeans(rS)
  
  pS1 <- colMeans(rS1)
  pS0 <- colMeans(rS0)
  
  list(rATE=rATE,
       pS1=pS1,
       pS0=pS0)
}

#----------------------------------------------------------------------------------
Sdiff <- function(est){
  return(est[,2]-est[,1]) 
} 

#----------------------------------------------------------------------------------
easy_boot_PS <- function(data,indices,form_ps,setTime){
  tmpdat <- data[indices,]
  
  mfull_ipw <- glm(form_ps,family="binomial",data=tmpdat)
  predfull_ipw <- mfull_ipw$fitted.values
  wfull <- 1/(tmpdat$A*predfull_ipw+(1-tmpdat$A)*(1-predfull_ipw))
  
  survfull_ipw_times <- survfit(Surv(time,event) ~ A, data=tmpdat, weights = wfull)
  ss_times <- summary(survfull_ipw_times, times=setTime)
  df_times <- data.frame(time=ss_times$time,A=ss_times$strata,surv=ss_times$surv)
  df_less1 <- df_times$surv[df_times$A=="A=1"]
  df_less0 <- df_times$surv[df_times$A=="A=0"]
  
  res <- c(df_less0, df_less1, df_less1 - df_less0)
  names(res) <- paste0(rep(paste0("t.",setTime,"."),3),rep(c("A=0","A=1","Diff"),each=length(setTime)))
  return(res)
}
