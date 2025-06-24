## ----setup, include=FALSE, warning=FALSE, message=FALSE-------------------------------------------
options(width = 100)
require(dtPosidonia)
require(marelac)

## -------------------------------------------------------------------------------------------------
ls("package:dtPosidonia")

## -------------------------------------------------------------------------------------------------
args(Posidonia_run)

## -------------------------------------------------------------------------------------------------
args(Posidonia_get_parms)

## -------------------------------------------------------------------------------------------------
args(Posidonia_get_vars)

## ----include=FALSE--------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -------------------------------------------------------------------------------------------------
Parms  <- Posidonia_get_parms(as.vector = TRUE)
States <- Posidonia_get_states(as.vector = TRUE)

## -------------------------------------------------------------------------------------------------
extract_data <- function(data,  # data.frame with 3 columns (variable, day, value)
                         name){
  
  data.extract   <- subset(data, 
                           subset = variable == name, 
                           select = c(day, value))
  
  colnames(data.extract)[2] <- name
  
  return(data.extract)
}

## -------------------------------------------------------------------------------------------------
data.oxygen   <- extract_data(obs_calvi, "OXYGEN")
data.leaves   <- extract_data(obs_calvi, "LEAVES")
data.roots    <- extract_data(obs_calvi, "ROOTS")  
data.rhizomes <- extract_data(obs_calvi, "RHIZOMES")

## -------------------------------------------------------------------------------------------------
data.temp     <- extract_data(forcs_calvi, "f_Temperature")
data.light    <- extract_data(forcs_calvi, "f_Light")
data.wind     <- extract_data(forcs_calvi, "f_Wind")

funtemp       <- approxfun(x = data.temp$day, 
                           y = data.temp$f_Temperature, rule = 2)
funlight      <- approxfun(x = data.light$day,       
                           y = data.light$f_Light,      rule = 2)
funwind       <- approxfun(x = data.wind$day,        
                           y = data.wind$f_Wind,        rule = 2)

Day   <- trunc(data.wind$day)
Wind  <- aggregate(x   = data.wind$f_Wind, 
                   by  = list(day = Day), 
                   FUN = mean)
fwind <- approxfun(x    = Wind$day, 
                   y    = Wind$x, 
                   rule = 2)

# Estimated with marelac functions
Time <- seq(from = 0., to = 365, by = 1)

ftemp     <- cbind(Time, funtemp (Time))
flight    <- cbind(Time, funlight(Time))
fwind     <- cbind(Time, funwind (Time))

funpiston  <- approxfun(x = Time, 
                        y = gas_transfer(t   = ftemp[,2], 
                                         u10 = fwind[,2], 
                                         species = "O2",
                                         method = "Nightingale")*86400)  # m/d
funSatox   <- approxfun(x = Time, 
                        y = gas_O2sat(S = Parms["salinity"], 
                                      t = ftemp[,2]) / molweight("O2")) #in mol/m3
fpiston    <- cbind(Time, funpiston(Time))
fSatox     <- cbind(Time, funSatox(Time))


## -------------------------------------------------------------------------------------------------
times  <- seq(from = 0, to = 365, by = 1)

 P    <- Posidonia_run(times = times, yini = States, f_Temperature = ftemp, 
                       f_Light = flight, f_Wind = fwind, f_Piston = fpiston,
                       f_SatO2 = fSatox)

 
Parms2 <- list(maxPhotosynthesis= 0.05)   # /D 

 P2    <- Posidonia_run(parms = Parms2, 
                       times = times, yini = States, f_Temperature = ftemp, 
                       f_Light = flight, f_Wind = fwind, f_Piston = fpiston,
                       f_SatO2 = fSatox)

Parms3 <- list(maxPhotosynthesis= 0.005)   # /D 
 P3    <- Posidonia_run(parms = Parms3, 
                       times = times, yini = States, f_Temperature = ftemp, 
                       f_Light = flight, f_Wind = fwind, f_Piston = fpiston,
                       f_SatO2 = fSatox)


## ----fig.width = 8, fig.height=10-----------------------------------------------------------------
plot(P, P2, P3, mfrow = c(4, 3), which = 1:12)
plot(P, P2, P3, mfrow = c(4, 3), which = 13:24)


## ----fig.width = 8, fig.height=10-----------------------------------------------------------------
par  (mfrow = c(2,2))
plot(P, P2, P3, obs = obs_calvi, 
     which = c("LEAVES", "ROOTS", "RHIZOMES", "OXYGEN"),
     xlab = "Days", ylab = c("molC/m2", "molC/m2", "molC/m2", "mol O2/m2"))

## -------------------------------------------------------------------------------------------------
Q  <- as.data.frame(P)

## CARBON budget 
F.DICtoLEAVES               <- mean(Q$GrossPrimaryProduction)-mean(Q$LeafTotalRespiration)
F.LEAVEStoDETRITUS          <- mean(Q$LeafBreaking)
F.LEAVEStoDEADLEAVES        <- mean(Q$LeafMortality)
F.ASSIMILATEStoLEAVES       <- mean(Q$LeafAssimilation)
F.LEAVEStoASSIMILATES       <- mean(Q$GrossPrimaryProduction)+mean(Q$LeafNetReallocation)
F.DEADLEAVEStoDETRITUS      <- mean(Q$DeadLeafAbscission)+mean(Q$DeadLeafBreaking)
F.DEADLEAVEStoDIC           <- mean(Q$DeadLeafDecay)
F.RHIZOMEtoDIC              <- mean(Q$RhizomeTotalRespiration)
F.ASSIMILATEStoRHIZOMES     <- mean(Q$RhizomeAssimilation)
F.RHIZOMEStoASSIMILATES     <- mean(Q$RhizomeMobilisation)
F.RHIZOMEStoDEADBELOWGROUND <- mean(Q$RhizomeMortality)
F.ROOTtoDIC                 <- mean(Q$RootTotalRespiration)
F.ASSIMILATEStoROOTS        <- mean(Q$RootAssimilation)
F.ROOTtoDEADBELOWGROUND     <- mean(Q$RootMortality)
F.DEADBELOWGROUNDtoDIC      <- mean(Q$BelowgroundDecay)
F.DETRITUStoDIC             <- mean(Q$DetritusDecay)
F.DETRITUStoEXPORT          <- mean(Q$DetritusExport)

##Others
AnnualGPP                    <- mean(Q$GrossPrimaryProduction)
AnnualNPP                    <- mean(Q$NetO2production)
AnnualTotalRespiration       <- mean(Q$TotalRespiration)


## ----echo = FALSE, fig.width=10, fig.height=10----------------------------------------------------
require(diagram)
#Carbon Flows
names        <-c("DIC", "ASSIMILATES", "LEAVES", "DEADLEAVES", 
                 "RHIZOMES", "ROOTS", "DEADBELOWGROUND", "DETRITUS")
CarbonFlows  <- matrix(nrow=8, ncol=8, byrow=TRUE, data=c(
+ # DIC            ASSIMILATES          LEAVES           DEADLEAVES     RHIZOMES                 ROOTS       DEADBELOWGROUND     DETRITUS
+    0,               0,           F.DICtoLEAVES,           0,             0,                      0,              0,                 0,          #DIC
+    0,               0,       F.ASSIMILATEStoLEAVES,       0,  F.ASSIMILATEStoRHIZOMES, F.ASSIMILATEStoROOTS,     0,                 0,          #ASSIMILATES
+    0,   F.LEAVEStoASSIMILATES,           0,       F.LEAVEStoDEADLEAVES,  0,                      0,              0,      F.LEAVEStoDETRITUS,    #LEAVES
+ F.DEADLEAVEStoDIC,  0,                   0,               0,             0,                      0,              0,     F.DEADLEAVEStoDETRITUS, #DEADLEAVES
+ F.RHIZOMEtoDIC, F.RHIZOMEStoASSIMILATES, 0,               0,             0,                      0,    F.RHIZOMEStoDEADBELOWGROUND, 0,          #RHIZOMES
+ F.ROOTtoDIC,        0,                   0,               0,             0,                      0,    F.ROOTtoDEADBELOWGROUND,     0,          #ROOTS
+ F.DEADBELOWGROUNDtoDIC, 0,               0,               0,             0,                      0,              0,                 0,          #DEADBELOWGROUND
+ F.DETRITUStoDIC,    0,                   0,               0,             0,                      0,              0,                 0           #DETRITUS
))

dimnames(CarbonFlows) <- list(names,names)
par(mar=c(1,1,1,1), mfrow=c(1,1))
plotweb(CarbonFlows, main="Posidonia Carbon Flow", sub="molC/m2/day", val=TRUE, 
        budget=TRUE, lab.size = 1.2, maxflow = 0.03, minflow = 0.001)

