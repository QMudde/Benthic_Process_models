## ----setup, include=FALSE, warning=FALSE, message=FALSE-------------------------------------------
options(width = 100)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##>",
  fig.width=8,
  fig.height=8
)

library(deSolve)
library(BFFMod)

## -------------------------------------------------------------------------------------------------
ls("package:BFFMod")

## -------------------------------------------------------------------------------------------------
args(mussel_run_debind)

## -------------------------------------------------------------------------------------------------
args(mussel_run_debcohort)

## -------------------------------------------------------------------------------------------------
args(mussel_get_debparms)

## -------------------------------------------------------------------------------------------------
args(mussel_get_debvars)

## -------------------------------------------------------------------------------------------------
args(mussel_get_debstates)

## -------------------------------------------------------------------------------------------------
mod1 <- mussel_run_debind()
plot(mod1, which = c("STRUCT","RESERVE","REPROD","Length"))

## -------------------------------------------------------------------------------------------------
plot(mod1, which = c("f_Phyto","f_Temp","f_Detritus","Length"))
plot(mod1, which = c("RESERVE","Gonads","Gametes","Length"))

## -------------------------------------------------------------------------------------------------
HtHf <- mussel_run_debind(times = 0:(3*365), f_Phyto = 10, f_Temp = 20) 
HtLf <- mussel_run_debind(times = 0:(3*365), f_Phyto = 5, f_Temp = 20)
LtHf <- mussel_run_debind(times = 0:(3*365), f_Phyto = 10, f_Temp = 5) 
LtLf <- mussel_run_debind(times = 0:(3*365), f_Phyto = 5, f_Temp = 5)
plot(HtHf,HtLf,LtHf,LtLf, 
     which = c("STRUCT","RESERVE","REPROD","Length"), 
     lwd = 2)
legend("bottomright", col = c(1:4), legend = c("HtHf", "HtLf","LtHf","LtLf"), lwd = 2, cex=0.5)

## ----fig.width=8,fig.height=4---------------------------------------------------------------------
HtHf[,"Spawn_rate"] <- cumsum(HtHf[,"Spawn_rate"])
HtLf[,"Spawn_rate"] <- cumsum(HtLf[,"Spawn_rate"])
LtHf[,"Spawn_rate"] <- cumsum(LtHf[,"Spawn_rate"])
LtLf[,"Spawn_rate"] <- cumsum(LtLf[,"Spawn_rate"])

plot(HtHf,HtLf,LtHf,LtLf, main = c("Length", "Total reproductive output"),
     which = c("Length","Spawn_rate"), 
     lwd = 2)
legend("bottomright", col = c(1:4), legend = c("HtHf", "HtLf","LtHf","LtLf"), lwd = 2, cex=0.5)
  

## -------------------------------------------------------------------------------------------------
f_dynaTEMP = function(t) {(10 * sin((2*pi*(t-120)/365)) + 17)}
f_dynaPHYTO = function(t) {(5 * sin((2*pi*(t-120)/365)) + 6)}
f_dynaDT = function(t) {(20 * sin((2*pi*(t-200)/365)) + 30)}
DYNAmod <- mussel_run_debind(times = 0:(2*365), 
                             f_Phyto = f_dynaPHYTO,
                             f_Temp = f_dynaTEMP,
                             f_Detritus =f_dynaDT)
plot(DYNAmod, which = c("f_Phyto","f_Temp","f_Detritus","Length"))
plot(DYNAmod, which = c("STRUCT","RESERVE","REPROD","Length"))

## -------------------------------------------------------------------------------------------------
DYNAmod_HoxHsim <- mussel_run_debind(times = 0:(2*365), 
                             f_Phyto = f_dynaPHYTO,
                             f_Temp = f_dynaTEMP,
                             f_Detritus = f_dynaDT,
                             f_O2 = 400, f_Sim = 10)
DYNAmod_LoxHsim <- mussel_run_debind(times = 0:(2*365), 
                             f_Phyto = f_dynaPHYTO,
                             f_Temp = f_dynaTEMP,
                             f_Detritus = f_dynaDT,
                             f_O2 = 80, f_Sim = 10)
DYNAmod_HoxLsim <- mussel_run_debind(times = 0:(2*365), 
                             f_Phyto = f_dynaPHYTO,
                             f_Temp = f_dynaTEMP,
                             f_Detritus =f_dynaDT,
                             f_O2 = 400, f_Sim = 0.5)
DYNAmod_LoxLsim <- mussel_run_debind(times = 0:(2*365), 
                             f_Phyto = f_dynaPHYTO,
                             f_Temp = f_dynaTEMP,
                             f_Detritus =f_dynaDT,
                             f_O2 = 80, f_Sim = 0.5)

plot(DYNAmod_HoxHsim,DYNAmod_LoxHsim,DYNAmod_HoxLsim,DYNAmod_LoxLsim,
     which = c("STRUCT","RESERVE","REPROD","Length"),
     lwd = 2)
legend("bottomright", col = c(1:4), legend = c("HoxHsim", "LoxHsim","HoxLsim","LoxLsim"),
       lwd = 2, cex = 0.5)


## -------------------------------------------------------------------------------------------------
plot(HtHf,HtLf,LtHf,LtLf, 
     which = c("Clearance_rate","Assimilation","Respiration","Length"), 
     lwd = 2)
legend("bottomright", col = c(1:4), legend = c("HtHf", "HtLf","LtHf","LtLf"), lwd = 2, cex=0.7)

plot(DYNAmod_HoxHsim,DYNAmod_LoxHsim,DYNAmod_HoxLsim,DYNAmod_LoxLsim,
     which = c("Clearance_rate","Assimilation","Respiration","Length"),
     lwd = 2)
legend("bottomright", col = c(1:4), legend = c("HoxHsim", "LoxHsim","HoxLsim","LoxLsim"),
       lwd = 2, cex = 0.7)

## -------------------------------------------------------------------------------------------------
pop <- mussel_run_debcohort(times = 0:(5*365))
plot(pop, which = c("Total_POP_m2","Total_WW_m2","LARV_DENS","LARV_BIOM"))

## -------------------------------------------------------------------------------------------------
pop_HtHf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 10, f_Temp = 20) 
pop_HtLf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 5, f_Temp = 20)
pop_LtHf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 10, f_Temp = 5) 
pop_LtLf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 5, f_Temp = 5)

plot(pop_HtHf,pop_HtLf,pop_LtHf,pop_LtLf,
     which = c("Total_POP_m2","Total_WW_m2","LARV_DENS","LARV_BIOM"))

## -------------------------------------------------------------------------------------------------
pop_HtHf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 10, f_Temp = 20,dynamic_environment = FALSE) 
pop_HtLf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 5, f_Temp = 20,dynamic_environment = FALSE)
pop_LtHf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 10, f_Temp = 5,dynamic_environment = FALSE) 
pop_LtLf <- mussel_run_debcohort(times = 0:(5*365), f_Phyto = 5, f_Temp = 5,dynamic_environment = FALSE)

plot(pop_HtHf,pop_HtLf,pop_LtHf,pop_LtLf,
     which = c("Total_POP_m2","Total_WW_m2","LARV_DENS","LARV_BIOM"))

## -------------------------------------------------------------------------------------------------
pop_dyna <- mussel_run_debcohort(times = 100:(5*365),parms = c(water_renew = 1),
                                f_Phyto = f_dynaPHYTO,
                                f_Temp = f_dynaTEMP,
                                f_Detritus = f_dynaDT) 
pop_dyna2 <- mussel_run_debcohort(times = 100:(5*365),parms = c(water_renew = 2),
                                f_Phyto = f_dynaPHYTO,
                                f_Temp = f_dynaTEMP,
                                f_Detritus = f_dynaDT)
pop_dyna3 <- mussel_run_debcohort(times = 100:(5*365),parms = c(water_renew = 3),
                                f_Phyto = f_dynaPHYTO,
                                f_Temp = f_dynaTEMP,
                                f_Detritus = f_dynaDT)


plot(pop_dyna,pop_dyna2,pop_dyna3,
     which = c("Total_POP_m2","Total_WW_m2","LARV_DENS","LARV_BIOM"))

