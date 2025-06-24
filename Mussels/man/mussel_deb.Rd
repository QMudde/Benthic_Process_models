\name{Mussel dynamics for one indidual based on dynamic energy budget}

\alias{mussel_run_debind}

\title{
  Physiological model of the dynamics of an average mussel individual based on Dynamic Energy Budget.
}
\description{
  \code{mussel_run_debind} is a dynamic energy budget (\code{DEB}) model that describes the dynamics of mussel structure, reserve and reproductive carbon in one \code{ind}ividual, under environmental forcings.  The physiology of the mussel is described with three state variables:

\itemize{
  \item  RESERVE,                  Reserve carbon per individual, [mmolC/ind]
  \item   STRUCT                Structural carbon per individual, [mmolC/ind]
  \item   REPROD  Reproductive carbon (gonads + gametes) per ind, [mmolC/ind]
}  
The feeding mechanism is adapted from \code{Saraiva et al., 2011} and is specifically developed to 
model benthic filter feeders in a dynamic environments with multiple particle types. Environmental 
forcings that the model works with dynamically are: 

\itemize{
  \item Temperature                       [dgC]
  \item Phytoplankton                     [mmolC/m3]
  \item Detritus                          [mmolC/m3]
  \item Suspended inorganic matter (SIM)  [g/m3]
  \item Oxygen                            [mmolO2/m3]
}
}

\usage{

mussel_run_debind(parms = list(), times = 0:365, yini = NULL, 
                  f_Phyto = 10, f_Temp = 20, f_Detritus = 10, f_Sim = 1, f_O2 = 360)
           
}
\arguments{

  \item{f_Phyto, f_Temp, f_Detritus, f_Sim, f_O2}{ Forcing functions that represent external environmental conditions. The default is to keep the forcing functions constant, but it is also allowed to pass a data set consisting of (time, value). f_Phyto represents phytoplankton concentrations [mmolC/m3], f_Temp is the temperature [dgC], f_Detritus is detrital organic carbon, [mmolC/m3], f_Sim is suspended inorganic matter, [g/m3], f_O2 is oxygen concentration, [mmolO2/m3].
  }  
  \item{yini}{ Initial condition. If NULL, then the initial condition is estimated based on the parameters that define the weight of settling larvae (weight_settle) and the fraction that goes to reserves (pReserve_settle) as follows: c(RESERVE = weight_settle*pReserve_settle, STRUCT = weight_settle*(1-pReserve_settle), REPROD = 0.).
  }
  \item{times}{ simulation time specification, in days. Usually a vector. If one single value, then the model is is called once, and the derivatives and output variables are returned.
  }
  \item{parms}{ vector or list with parameter values that overrule the default parameter values. See examples for parameters and their default values.
  }
}
\value{
The return values consists of an object of class * dtDyn* and of class *deSolve* as returned by the \link[deSolve]{ode} function of package deSolve. The first column consists of the output \code{times}. Then several output variables are returned. 

The object has several useful attributes, such as the units of all output variables, and the parameter settings.
}

\author{
  Karline Soetaert and Quinten Mudde
}
\details{ 
The parameters for the mussel deb models, and their default values and units can be found via the function \code{mussel_parms}. 
}
\examples{

# ---------------------------------------------
# ---------------------------------------------
# Simple DEB dynamics
# ---------------------------------------------
# ---------------------------------------------

yini  <- c(0.1, 0.1, 0.1)
times <- 0:1095

out  <- mussel_run_debind(times = times, 
                          yini  = yini)

out2 <- mussel_run_debind(times   = times, 
                          yini    = yini, 
                          f_Phyto = 0)

out3 <- mussel_run_debind(times  = times, 
                          yini   = yini, 
                          f_Temp = 10)

# incrase eg to make assimilaiton of structure very costly
out4  <- mussel_run_debind(times = times, 
                           yini  = yini, 
                           parms = list(eg = 0.25))
                    
plot(out, out2, out3, out4, 
     lty = 1, lwd = 2)

# ---------------------------------------------
# fluctuating forcing functions
# ---------------------------------------------

temp.fluct = function(t) (   6 * sin((2*pi*(t-50)/365)) + 20) 
phy.fluct  = function(t) (12.5 * sin((2*pi*(t   )/365)) + 12.5)

out5 <- mussel_run_debind(times  = times, 
                          yini    = yini, 
                          f_Temp  = temp.fluct, 
                          f_Phyto = phy.fluct)

plot(out5, 
     which = c("RESERVE", "STRUCT", "REPROD", 
     "Gametes", "Length", "Spawn_rate", "C_ind", 
     "f_Phyto", "f_Temp"))

# ---------------------------------------------
# Sensitivity analysis
# ---------------------------------------------

times   <- seq(from = 0, to = 5*365, by = 10)
Phy.seq <- seq(from = 0, to = 100, by = 1)
parms   <- list(spawn_eff = 0) # no spawning

out <- list()
SV <- NULL
for (i in 1: length(Phy.seq)){
  out[[i]] <- mussel_run_debind(times = times, yini = yini, 
                                parms = parms,
                                f_Phyto = Phy.seq[i])
                            
  SV <- rbind(SV, out[[i]][nrow(out[[i]]), 2:8])                       
}

head(SV)
par(mfrow = c(2,2))

plot(Phy.seq, SV[,"RESERVE"]/ SV[,"STRUCT"], 
     xlab = "Phytoplankton C", main = "R/S", 
     ylab = "molC/molC")
     
plot(Phy.seq, SV[,"Length"], 
     xlab = "Phytoplankton C", main = "Length", 
     ylab = "cm")

plot(out[[1]], out[[50]], out[[100]], 
     which = c("REPROD", "Length"), 
     xlab = "day", lty=1, lwd = 2, mfrow = NULL)
     
legend("bottomright", legend = Phy.seq[c(1, 50, 100)], 
       col=1:3, lty = 1)     
       

}


\keyword{ utilities }
