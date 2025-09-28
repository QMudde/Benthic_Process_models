\name{Macroalgae model}

\alias{Macroalgae_run}

\title{
  Model of the dynamics of seagrasses in the environment, based on Phyllophora.
}
\description{
  \code{Macroalgae_run} runs a mechanistic model that describes the dynamics of seagrass leaves, roots and rhizomes, dead leaves, dead belowground tissue, under environmental forcings.  
  
  The seagrasses are described with the following state variables:

\itemize{
  \item  STRUCT_C: structural C of macroalgae, mmolC/m2
  \item  RESERVE_C: reserve C of macroalgae, mmolC/m2
  \item  RESERVE_N: reserve N of macroalgae, mmolN/m2
  \item  RESERVE_P: reserve P of macroalgae, mmolP/m2
  \item  NO3 nitrate concentration in water, mmolN/m3
  \item  NH4: ammonium concentration in water, mmolN/m3
  \item  PO4: phosphate concentration in water, mmolP/m3
  \item  O2: oxygen concentration in water, mmolO2/m3
  \item  C_CLOSING: all C that goes in and out of the Box model, mmolC/m2
  \item  N_CLOSING: all N that goes in and out of the Box model, mmolN/m2
  \item  P_CLOSING: all P that goes in and out of the BOX model, mmolP/m2
  }  
}

\usage{

Macroalgae_run(parms = list(), times = 0:365, yini = NULL, 
           f_Temperature = 20, f_Current = 1, f_Light = 100,  
           f_PHYTO = 10, f_ExtNO3 = 20, f_ExtNH4 = 1, 
           f_ExtPO4 = 1, f_ExtO2 = 200)    
           
}
\arguments{

  \item{f_Temperature, f_Light, f_Current, f_PHYTO, 
  f_ExtNO3, f_ExtNH4, f_ExtPO4, f_ExtO2}{ Forcing 
  functions that represent external environmental conditions. 
  The default is to keep the forcing functions constant, but it is also allowed to 
  pass a data set consisting of (time, value). 
  f_Temperature is the temperature of the water [dgC], 
  f_Current is the current velocity [m/s],
  f_Light is the light intensity [W/m2], 
  f_PHYTO is the phytoplankton concentration in the water, [mmolC/m3], 
  f_ExtNO3, f_ExtNH4 are the external nitrate and ammoniume concentration, [mmol N/m3] 
  f_ExtPO4 is the external phosphate concentration, [mmol P/m3] 
  f_ExtO2 is the external oxygen concentration, [mol O2/m3].
  }  
  \item{yini}{ Initial condition. If NULL, then a default value is used (see \link{Macroalgae_get_states} to extract those).
  Else a vector or list with the values of one of ("STRUCT_C", "RESERVE_C", "RESERVE_N", 
  "RESERVE_P", "NO3", "NH4", "PO4", "O2", "C_CLOSING", "N_CLOSING", "P_CLOSING").
  If the names in \code{yini} are not provided, it is assumed that all state variables are passed, in the correct order as given above.
  }
  \item{times}{ simulation time specification, in days. 
  If one value, the model is called once, and the derivatives and output variables are returned.
  }
  \item{parms}{ vector or list with parameter values that overrule the default 
  parameter values (see \link{Macroalgae_get_parms} to extract those). 
  }
}
\value{
The return values consists of an object of class *dtMacroalgae* and of class *deSolve* as returned by the \link[deSolve]{ode} function of package deSolve. The first column consists of the output \code{times}. Then several output variables are returned. 

The object has several useful attributes, such as the units of all output variables, and the parameter settings.
}

\author{
  Karline Soetaert, Quinten Mudde, Joran de Gang, Marilaure Gregoire
}
\details{ 
The parameters for the model, and their default values and units can be found via the function \code{Macroalgae_get_parms}. 

The names and units of the state variables and their default initial conditions can be found via function \link{Macroalgae_get_states}.
}
\examples{

# ---------------------------------------------
# ---------------------------------------------
# default parameters, states and forcings
# ---------------------------------------------
# ---------------------------------------------

Macroalgae_get_parms()
Macroalgae_get_states()
Macroalgae_get_forcings()

# ---------------------------------------------
# ---------------------------------------------
# default model run
# ---------------------------------------------
# ---------------------------------------------

out1 <- Macroalgae_run()       # default run
plot(out1, which = 1:6)

parms2 <- c(pMax = 0.05)       # double the default
out2 <- Macroalgae_run(parms = parms2)

out3 <- Macroalgae_run(f_Temperature = 10)   # half the default
plot(out1, out2, out3, which = 1:6)

out4 <- Macroalgae_run(f_Light = 1000)
plot(out1, out2, out3, out4, which = 1:4)

# ---------------------------------------------
# fluctuating forcing functions
# ---------------------------------------------

t <- 0:730
temp.fluct   =  20  + 6 * sin((2*pi*(t-50)/365)) 
light.fluct  =  170 - 150*sin((2*pi*(t   )/365)) 

out5 <- Macroalgae_run(times = t,
                      f_Temperature  = cbind(t, temp.fluct), 
                      f_Light        = cbind(t, light.fluct))

plot(out5, which = c("Temperature", "Light", 
                     "STRUCT_C", "RESERVE_C", "RESERVE_N", "RN_SCr" ))


# ---------------------------------------------
# Sensitivity analysis
# ---------------------------------------------


}
\references{
  ... to be added?
}

\keyword{ utilities }
