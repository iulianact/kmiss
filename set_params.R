#----------------------------------------------------------------------------------
# Set parameters
#----------------------------------------------------------------------------------
nrep <- 500 
numimp <- 10
nboot <- 200
setTime <- c(1,5,10) 
stratyear <- 5

form_std_correct <- formula(Surv(time,event) ~ A + H + D + C + D:C + A:D) 
form_ps_correct <- formula(A ~ H + D + C + D:C + D:H)

form_std_mis <- formula(Surv(time,event) ~ A + H + D + C)
form_ps_mis <- formula(A ~ H + D + C)

formBias_std <- formula(Surv(time,event) ~ A + H + D + A:D)
formBias_ps <- formula(A ~ H + D + D:H)

# save(nrep,numimp,nboot,
#      setTime,stratyear,
#      file="general_params.RData")
#
# save(form_std_correct,form_ps_correct,
#      form_std_mis,form_ps_mis,
#      formBias_std,formBias_ps,
#      file="fit_params.RData")
