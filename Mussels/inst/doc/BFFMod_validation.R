## ----setup, warning=FALSE, include = FALSE--------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "##>",
  fig.width=6,
  fig.height=6
)
library(deSolve)
library(FME)
library(BFFMod)

## -------------------------------------------------------------------------------------------------
# default parameter values
fullpars <- mussel_get_debparms()
pars <- fullpars$default
names(pars)<-fullpars$names
# Datasets:
# physiological rate snapshots
getwd()
load("../data/M_gallo_Feeding_data.rda") 
# growth trajectories
load("../data/Sinop_data.rda")
load("../data/Marm_data.rda")
load("../data/Crimea_data.rda")

## -------------------------------------------------------------------------------------------------
## Sinop data 
sinop_o2_f <- make_ff(day = Sinop_data$day, values = Sinop_data$Oxygen_bamhbi, recalc = 'none')
Sinop_phyc <- make_ff(day = Sinop_data$day, values = Sinop_data$ChlorophyllA.mug_L, recalc = 'chl_mugl')
Sinop_DtC <- make_ff(day = Sinop_data$day, values = Sinop_data$POM.mg_L, recalc = 'pom_g')
Sinop_DtC[,2] <- Sinop_DtC[,2] - Sinop_phyc[,2] # substract phytoplankton from poc to get detritus
Sinop_SIM <- make_ff(day = Sinop_data$day, values = Sinop_data$PIM.mg_L, recalc = 'none')
Sinop_temp <- make_ff(day = Sinop_data$day, values = Sinop_data$Temperature.C, recalc = 'none')
Sinop_times <- seq(from = min(Sinop_data$day), to = max(Sinop_data$day), by = 1 )
# estimate initial conditions
Sinop_INI <- get_INI_indiv(Length = 4.03,WW = 1.52)
Sinop_obs <- data.frame(time = Sinop_data$day,
                        Length = Sinop_data$length_mean,
                        WW_ind = Sinop_data$WW_mean)

## Marmara data
Marm_o2_f <- make_ff(day = Marm_data$day, values = Marm_data$Oxygen_est, recalc = 'none')
Marm_DtC <- make_ff(day = Marm_data$day, values = Marm_data$POM.mg_L, recalc = 'pom_g')
Marm_phyc <- make_ff(day = Marm_data$day, values = Marm_data$ChlorA.mug_L, recalc = 'chl_mugl')
Marm_SIM <- make_ff(day = Marm_data$day, values = Marm_data$PIM.mg_L, recalc = 'none')
Marm_temp <- make_ff(day = Marm_data$day, values = Marm_data$Temperature.C, recalc = 'none')
Marm_DtC[,2] <- Marm_DtC[,2] - Marm_phyc[,2] # substract phytoplankton from poc to get detritus
Marm_times <- seq(from = min(Marm_data$day), to = max(Marm_data$day), by = 1 )
# estimate initial conditions: 
Marm_INI = get_INI_indiv(Length = 0.9819)
# Marmara weights are in total weight, recalculate to wet weight estimate
gTW_gWW <- mussel_get_auxiliary()[["gTW_gWW"]]
Marm_obs <- data.frame(time = Marm_data$day,
                       Length = Marm_data$length.mm/10,
                       WW_ind = Marm_data$tot.ww.g/gTW_gWW)

## Crimea data (2 growth curves)
Crim_o2_f <- make_ff(day = Crimea_data$times, values = Crimea_data$ox_mod, recalc = 'none')
Crim_DtC <- make_ff(day = Crimea_data$times, values = Crimea_data$dt_mod, recalc = 'none')
Crim_phyc <- make_ff(day = Crimea_data$times, values = Crimea_data$phyc_mod , recalc = 'none')
Crim_temp <- make_ff(day = Crimea_data$times, values = Crimea_data$temperature, recalc = 'none')
POM_g_to_mmolC <- mussel_get_auxiliary()[["POM_g_to_mmolC"]]
Crim_SIM <- cbind(Crimea_data$times,(Crim_DtC[,2]+Crim_phyc[,2])/POM_g_to_mmolC) # equal weights SIM & POM 
Crim_times <- seq(from = min(Crimea_data$times), to = max(Crimea_data$times), by = 1 )
# initial conditions
Crim_INI_35 = get_INI_indiv(Length = 3.5, Rel_E_dens  = 0.5) # half relative energy density because we start after winter
Crim_INI_60 = get_INI_indiv(Length = 6, Rel_E_dens  = 0.5) # half relative energy density because we start after winter
# sets obs dat
Crim_35_obs <- data.frame(time = Crimea_data$times,
                          Length = Crimea_data$ini_35_size)
Crim_60_obs <- data.frame(time = Crimea_data$times,
                          Length = Crimea_data$ini_60_size)

## -------------------------------------------------------------------------------------------------
FEEDPs <- c("max_clear", "max_phyto_filt", "max_DT_filt", "max_spm_filt",
            "rho_Phyto","rho_DT","rho_SIM","max.ing.Phy","max.ing.DT",
            "max_SIM_ing", "aeff")

DEBPs <- c("kappa", "maint","cat_s","eg", "mat_pub")
ENVPs <- c("ks_O2","crit_O2","T_a","T_l","T_h","T_ah","T_al")
senspars <- c(FEEDPs,DEBPs)
# sensitivity of 10 %
SENS_RANGE <- data.frame(min = pars[senspars]*0.9, max = pars[senspars]*1.1)

## ----fig.width=6, fig.height=6--------------------------------------------------------------------

sensrun_sinop <-function(pars, outtimes)
{
  P <- mussel_run_debind(parms = pars, yini = Sinop_INI,
                         times = Sinop_times,
                         f_Phyto = Sinop_phyc ,
                         f_Temp = Sinop_temp,
                         f_Detritus = Sinop_DtC,
                         f_O2 = sinop_o2_f,
                         f_Sim = Sinop_SIM)
  return(as.data.frame(P))
}
SENS_SINOP <- sensRange(func = sensrun_sinop, parms = pars, dist = "unif",
                        parRange =SENS_RANGE, sensvar=c("Length","WW_ind","STRUCT","Gonads","Gametes","RESERVE"),
                        num=50)
Gsens_SINOP<-summary(SENS_SINOP)
plot(Gsens_SINOP,which = c("WW_ind","STRUCT","RESERVE","Gametes"))

# SINOP LENGTH PLOT

plot(Gsens_SINOP,which = c("Length"),
     ylab = c("cm"),
     xlab = c("days"),
     main = c("Shell Length : SINOP"))
points(x = Sinop_data$day , y = Sinop_data$length_mean, pch = 19, cex = 1)
arrows(Sinop_data$day, Sinop_data$length_upper, 
       Sinop_data$day, Sinop_data$length_lower, 
       length=0.05, angle=90, code=3)

# SINOP WW PLOT

plot(Gsens_SINOP,which = c("WW_ind"),
     ylab = c("gram"),
     xlab = c("days"),
     main = c("Wet Weight : SINOP"))

points(x = Sinop_data$day , y = Sinop_data$WW_mean, pch = 19, cex = 1)
arrows(Sinop_data$day, Sinop_data$WW_upper, 
       Sinop_data$day, Sinop_data$WW_lower, 
       length=0.05, angle=90, code=3)

## ----fig.width=6, fig.height=6--------------------------------------------------------------------
sensrun_marmara <-function(pars, outtimes)
{
  P <- mussel_run_debind(parms = pars, yini = Marm_INI,
                         times = Marm_times,
                         f_Phyto = Marm_phyc ,
                         f_Temp = Marm_temp,
                         f_Detritus = Marm_DtC,
                         f_Sim = Marm_SIM,
                         f_O2 = Marm_o2_f)
  return(as.data.frame(P))
}
SENS_MARM <- sensRange(func = sensrun_marmara, parms = pars, dist = "unif",
                       parRange =SENS_RANGE, sensvar=c("STRUCT","RESERVE","REPROD","Length","WW_ind"),
                       num=50)
Gsens_MARM<-summary(SENS_MARM)
plot(Gsens_MARM,which = c("STRUCT","RESERVE","REPROD"))
# Marmara Length
plot(Gsens_MARM,which = c("Length"), 
     ylab = c("cm shell"),
     xlab = c("days"))
points(x = Marm_data$day , y = Marm_data$length.mm/10, pch = 19, cex = 1)
arrows(Marm_data$day, Marm_data$length.mm.min/10, 
       Marm_data$day, Marm_data$length.mm.max/10, 
       length=0.05, angle=90, code=3)
# plot marmara weight
plot(Gsens_MARM,which = c("WW_ind"), 
     ylab = c("gram"),
     xlab = c("days"))
points(x = Marm_data$day , y = Marm_data$ww_g_est, pch = 19, cex = 1)
arrows(Marm_data$day, Marm_data$ww_g_est - Marm_data$ww_g_est_sd, 
       Marm_data$day, Marm_data$ww_g_est + Marm_data$ww_g_est_sd, 
       length=0.05, angle=90, code=3)

## ----fig.width=6, fig.height=6, warning = FALSE---------------------------------------------------
sensrun_Crim35 <-function(pars, outtimes)
{
  P <- mussel_run_debind(parms = pars, yini = Crim_INI_35,
                         times = Crim_times,
                         f_Phyto = Crim_phyc,
                         f_Temp = Crim_temp,
                         f_Detritus = Crim_DtC,
                         f_O2 = Crim_o2_f,
                         f_Sim = Crim_SIM)
  return(as.data.frame(P))
}

SENS_CRIM35 <- sensRange(func = sensrun_Crim35, parms = pars, dist = "unif",
                         parRange =SENS_RANGE, sensvar=c("STRUCT","RESERVE","REPROD","Length","WW_ind"),
                         num=50)
Gsens_CRIM35<-summary(SENS_CRIM35)
plot(Gsens_CRIM35,which = c("Length"), 
     ylab = "cm shell", xlab = "days", main = "Length, 35 mm start")
points(x = Crimea_data$times , y = Crimea_data$ini_35_size, pch = 19, cex = 1)
arrows(Crimea_data$times, Crimea_data$ini_35_size - Crimea_data$ini_35_sd, 
       Crimea_data$times, Crimea_data$ini_35_size + Crimea_data$ini_35_sd, 
       length=0.05, angle=90, code=3)


# dataset 2 mussel growing from 60mm
sensrun_Crim60 <-function(pars, outtimes)
{
  P <- mussel_run_debind(parms = pars, yini = Crim_INI_60,
                         times = Crim_times,
                         f_Phyto = Crim_phyc ,
                         f_Temp = Crim_temp,
                         f_Detritus = Crim_DtC,
                         f_O2 = Crim_o2_f,
                         f_Sim = Crim_SIM)
  return(as.data.frame(P))
}
# run sensitivity Crim 35 
SENS_CRIM60 <- sensRange(func = sensrun_Crim60, parms = pars, dist = "unif",
                         parRange =SENS_RANGE, sensvar=c("STRUCT","RESERVE","REPROD","Length","WW_ind"),
                         num=50)
Gsens_CRIM60<-summary(SENS_CRIM60)
plot(Gsens_CRIM60,which = c("Length"), 
     ylab = "cm shell", xlab = "days", main = "Length, 60 mm start")
points(x = Crimea_data$times , y = Crimea_data$ini_60_size, pch = 19, cex = 1)
arrows(Crimea_data$times, Crimea_data$ini_60_size - Crimea_data$ini_60_sd, 
       Crimea_data$times, Crimea_data$ini_60_size + Crimea_data$ini_60_sd, 
       length=0.05, angle=90, code=3)

## ----fig.width=6, fig.height=6--------------------------------------------------------------------
RATES.FIT <- data.frame(index = NA, Clear = NA, OrgFilt = NA,Reject = NA, OrgIng = NA, Assim = NA)
RATES.low <- data.frame(index = NA, Clear = NA, OrgFilt = NA,Reject = NA, OrgIng = NA, Assim = NA)
RATES.high <- data.frame(index = NA, Clear = NA, OrgFilt = NA,Reject = NA, OrgIng = NA, Assim = NA)
for(i in 1:nrow(Feeding_data)){
  
  
  
  F_INI = get_INI_indiv(Length = Feeding_data$Length[i])
  
  fphy = cbind(1:2, rep(Feeding_data$PhyC[i],2))
  ftemp = cbind(1:2, rep(Feeding_data$Temp[i],2))
  fsim = cbind(1:2, rep(Feeding_data$SIM[i],2))
  fdet =   cbind(1:2, rep(Feeding_data$DT[i],2))
  
  OUT= mussel_run_debind(times = 1:2,
                         yini = F_INI,
                         f_Phyto = fphy,
                         f_Temp = ftemp,
                         f_Sim = fsim,
                         f_Detritus = fdet,
                         parms = pars)
  
  
  minvals= SENS_RANGE$min
  names(minvals)<- rownames(SENS_RANGE)
  maxvals= SENS_RANGE$max
  names(maxvals)<- rownames(SENS_RANGE)
  
  fitP_low <- pars
  fitP_low[names(minvals)] <- minvals
  fitP_high <- pars
  fitP_high[names(maxvals)] <- maxvals
  
  OUT_h= mussel_run_debind(times = 1:2,
                           yini = F_INI,
                           f_Phyto = fphy,
                           f_Temp = ftemp,
                           f_Sim = fsim,
                           f_Detritus = fdet,
                           parms = fitP_high)
  
  
  
  OUT_l= mussel_run_debind(times = 1:2,
                           yini = F_INI,
                           f_Phyto = fphy,
                           f_Temp = ftemp,
                           f_Sim = fsim,
                           f_Detritus = fdet,
                           parms = fitP_low)
  
  
  
  
  # need specific unit transformations to compare with dataset
  
  outdf = as.data.frame(OUT)
  outdf_h = as.data.frame(OUT_h)
  outdf_l = as.data.frame(OUT_l)
  
  
  RATES.FIT[i,]<- c(i,
                    outdf$Clearance_rate[1], #[clearance rate m3/d -> l/h]
                    outdf$Filtration_C[1],
                    outdf$Pseudofaeces_C[1]/(0.001+outdf$Filtration_C[1])*100,# % of filtred C rejected
                    outdf$Ingestion[1],
                    outdf$Assimilation[1]/pars["aeff"])
  
  RATES.high[i,]<- c(i,
                     outdf_h$Clearance_rate[1], #[clearance rate m3/d -> l/h]
                     outdf_h$Filtration_C[1],
                     outdf_h$Pseudofaeces_C[1]/(0.001+outdf_h$Filtration_C[1])*100,# % of filtred C rejected
                     outdf_h$Ingestion[1],
                     outdf_h$Assimilation[1]/fitP_high["aeff"])
  
  RATES.low[i,]<- c(i,
                    outdf_l$Clearance_rate[1], #[clearance rate m3/d -> l/h]
                    outdf_l$Filtration_C[1],
                    outdf_l$Pseudofaeces_C[1]/(0.001+outdf_l$Filtration_C[1])*100,# % of filtred C rejected
                    outdf_l$Ingestion[1],
                    outdf_l$Assimilation[1]/fitP_low["aeff"])
}

RATES.FIT <- na.omit(RATES.FIT)
RATES.high <- na.omit(RATES.high)
RATES.low <- na.omit(RATES.low)

## ----fig.width=6, fig.height=6--------------------------------------------------------------------
# plot for clearance rate
plot(x = Feeding_data$index, y = Feeding_data$Clear, ylim = c(0,0.1), cex = 1.5, lwd = 2,pch = 19,
     xlab = "index", ylab = "Clear rate [m3/d]",main = "Clearance rate")
arrows(Feeding_data$index, Feeding_data$Clear + Feeding_data$Clear_sd*2, 
       Feeding_data$index, Feeding_data$Clear - Feeding_data$Clear_sd*2, 
       length=0.05, angle=90, code=3)

points(x = RATES.FIT$index, y = RATES.FIT$Clear, col = "blue", pch = 4,cex = 1.5, lwd = 2)
arrows(RATES.FIT$index, RATES.low$Clear, 
       RATES.FIT$index, RATES.high$Clear, 
       length=0.05, angle=90, code=3, col = "blue")



legend("topright", pch = c(19,4), col = c("black","blue"), cex = 0.8,
       legend = c("Data","Modeled values"))

## ----fig.width=6, fig.height=6--------------------------------------------------------------------
# plot for filtration rate
plot(x = Feeding_data$index, y = Feeding_data$OrgFilt, ylim = c(0,4), cex = 1.5, lwd = 2,pch = 19,
     xlab = "index", ylab = "C filtration rate [mmolC/d]",main = "Organic Filtration")
arrows(Feeding_data$index, Feeding_data$OrgFilt + Feeding_data$OrgFilt_sd*2, 
       Feeding_data$index, Feeding_data$OrgFilt - Feeding_data$OrgFilt_sd*2, 
       length=0.05, angle=90, code=3)

points(x = RATES.FIT$index, y = RATES.FIT$OrgFilt, col = "blue", pch = 4,cex = 1.5, lwd = 2)
arrows(RATES.FIT$index, RATES.low$OrgFilt, 
       RATES.FIT$index, RATES.high$OrgFilt, 
       length=0.05, angle=90, code=3, col = "blue")
legend("topright", pch = c(19,4), col = c("black","blue"), cex = 0.8,
       legend = c("Data","Modeled values"))

## REJECTION 
#  plot(x = Feeding_data$index, y = Feeding_data$Reject, ylim = c(0,50), cex = 1.5, lwd = 2,pch = 19,
#       xlab = "index", ylab = "Rejection [%]",main = "Rejection")
#  points(x = Feeding_data$index, y = Feeding_data$Reject + Feeding_data$Reject_sd)
#  points(x = Feeding_data$index, y = Feeding_data$Reject - Feeding_data$Reject_sd)

#  points(x = default.run$index, y = default.run$Reject, col = "red", pch = 2,cex = 1.5, lwd = 2)
#  points(x = FIT.run$index, y = FIT.run$Reject, col = "blue", pch = 3,cex = 1.5, lwd = 2)
#  legend("topright", pch = c(19,2,3), col = c("black","red","blue"), cex = 0.8,
#         legend = c("DATA","MOD_Guess","MOD_FIT"))


## ----fig.width=6, fig.height=6--------------------------------------------------------------------

#
plot(x = Feeding_data$index, y = Feeding_data$OrgIng, ylim = c(0,3), cex = 1.5, lwd = 2,pch = 19,
     xlab = "index", ylab = "C Ingestion rate [mmolC/d]",main = "Organic Ingeston")
arrows(Feeding_data$index, Feeding_data$OrgIng + Feeding_data$OrgIng_sd*2, 
       Feeding_data$index, Feeding_data$OrgIng - Feeding_data$OrgIng_sd*2, 
       length=0.05, angle=90, code=3)


points(x = RATES.FIT$index, y = RATES.FIT$OrgIng, col = "blue", pch = 4,cex = 1.5, lwd = 2)
arrows(RATES.FIT$index, RATES.low$OrgIng, 
       RATES.FIT$index, RATES.high$OrgIng, 
       length=0.05, angle=90, code=3, col = "blue")
legend("topright", pch = c(19,4), col = c("black","blue"), cex = 0.8,
      legend = c("Data","Modeled values"))

#
plot(x = Feeding_data$index, y = Feeding_data$Assim, ylim = c(0,3), cex = 1.5, lwd = 2,pch = 19,
     xlab = "index", ylab = "Assimilation rate [mmolC/d]",main = "Assimilation")
arrows(Feeding_data$index, Feeding_data$Assim + Feeding_data$Assim_sd*2, 
       Feeding_data$index, Feeding_data$Assim - Feeding_data$Assim_sd*2, 
       length=0.05, angle=90, code=3)

# ! /aeff because feeding respiration is not counted into assimilation in model
points(x = RATES.FIT$index, y = (RATES.FIT$Assim/pars["aeff"]), col = "blue", pch = 4,cex = 1.5, lwd = 2)
arrows(RATES.FIT$index, RATES.low$Assim/pars["aeff"], 
       RATES.FIT$index, RATES.high$Assim/pars["aeff"], 
       length=0.05, angle=90, code=3, col = "blue")
legend("topright", pch = c(19,4), col = c("black","blue"), cex = 0.8,
       legend = c("Data","Modeled values"))

