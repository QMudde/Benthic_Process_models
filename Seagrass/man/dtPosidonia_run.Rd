\name{Posidonia model}

\alias{Posidonia_run}

\title{
  Model of the dynamics of seagrasses in the environment, based on Posidonia oceanica.
}
\description{
  \code{Posidonia_run} runs a mechanistic model that describes the dynamics of seagrass leaves, roots and rhizomes, dead leaves, dead belowground tissue, under environmental forcings.  
  
  The seagrasses are described with the following state variables:

\itemize{
  \item  LEAVES            C biomass of leaves  molC/m2
  \item  ROOTS             C biomass of roots molC/m2
  \item  RHIZOMES          C biomass of rhizomes molC/m2
  \item  NRESERVES         Nutrient reserves in rhizomes  molN/m2
  \item  DEADLEAVES        dead aboveground C biomass molC/m2
  \item  DEADBELOWGROUND   dead belowground C biomass molC/m2
  \item  BURIED            total C buried  molC/m2
  \item  DETRITUS          detritus biomass molC/ind
  \item  OXYGEN            pelagic oxygen concentration molO2/m3
  \item  DIN               pelagic dissolved inorganic nitrogen concentration  molN/m3
  \item  DINSED            benthic dissolved inorganic nitrogen concentration  molN/m2
  \item  CLOSING           variable to close the C budget (mainly DIC) molC/m2
  }  
}

\usage{

Posidonia_run(parms = list(), times = 0:365, yini = NULL, 
           f_Temperature = 20, f_Light = 100, f_Daylength = 12,
           f_Wind = 1, f_Piston = 1, f_SatO2 = 0.36)
           
}
\arguments{

  \item{f_Temperature, f_Light, f_Daylength, f_Wind, f_Piston, f_SatO2}{ Forcing 
  functions that represent external environmental conditions. 
  The default is to keep the forcing functions constant, but it is also allowed to 
  pass a data set consisting of (time, value). 
  f_Temperature is the temperature of the water [dgC], 
  f_Light is the light intensity [W/m2], f_Wind is the wind speed, [m/s], 
  f_Daylength is the daylength, in [hours].
  f_Piston is piston velocity (air-sea exchange), [m/d], 
  f_SatO2 is the saturated oxygen concentration, [molO2/m3].
  }  
  \item{yini}{ Initial condition. If NULL, then a default value is used (see \link{Posidonia_get_states} to extract those).
  Else a vector or list with the values of one of ("LEAVES", "ROOTS" , "RHIZOMES",
  "NRESERVES", "DEADLEAVES", "DEADBELOWGROUND", "BURIED", "DETRITUS", "OXYGEN", "DIN", "DINSED", "CLOSING").
  If the names in \code{yini} are not provided, it is assumed that all state variables are passed, in the correct order as given above.
  }
  \item{times}{ simulation time specification, in days. 
  If one value, the model is called once, and the derivatives and output variables are returned.
  }
  \item{parms}{ vector or list with parameter values that overrule the default 
  parameter values (see \link{Posidonia_get_parms} to extract those). 
  }
}
\value{
The return values consists of an object of class *dtPosidonia* and of class *deSolve* as returned by the \link[deSolve]{ode} function of package deSolve. The first column consists of the output \code{times}. Then several output variables are returned. 

The object has several useful attributes, such as the units of all output variables, and the parameter settings.
}

\author{
  Karline Soetaert, Quinten Mudde, Joran de Gang, Marilaure Gregoire
}
\details{ 
The parameters for the model, and their default values and units can be found via the function \code{Posidonia_get_parms}. 

The names and units of the state variables and their default initial conditions can be found via function \link{Posidonia_get_states}.
}
\examples{

# ---------------------------------------------
# ---------------------------------------------
# default parameters, states and forcings
# ---------------------------------------------
# ---------------------------------------------

Posidonia_get_parms()
Posidonia_get_states()
Posidonia_get_forcings()

# ---------------------------------------------
# ---------------------------------------------
# default model run
# ---------------------------------------------
# ---------------------------------------------

out1 <- Posidonia_run()
plot(out1, which = 1:6)

parms2 <- c(maxPhotosynthesis = 0.02)       # double the default
out2 <- Posidonia_run(parms = parms2)

out3 <- Posidonia_run(f_Temperature = 10)   # half the default
plot(out1, out2, out3, which = 1:6)

out4 <- Posidonia_run(f_Light = 1000)
plot(out1, out2, out3, out4, which = 1:6)

# ---------------------------------------------
# fluctuating forcing functions
# ---------------------------------------------

t <- 0:730
temp.fluct   =  20  + 6 * sin((2*pi*(t-50)/365)) 
light.fluct  =  170 - 150*sin((2*pi*(t   )/365)) 

# formula for daylength; latitude = 50
lat = 50
daylength = 12*(1+2/pi*asin(tan(lat*pi/180) * 
                            tan(23.45*pi/180*sin(2*pi*(284+t)/365))))

out5 <- Posidonia_run(times = t,
                      f_Temperature  = cbind(t, temp.fluct), 
                      f_Light        = cbind(t, light.fluct),
                      f_Daylength    = cbind(t, daylength))

plot(out5, which = 1:6)


# ---------------------------------------------
# Sensitivity analysis
# ---------------------------------------------


}
\references{
  ... to be added?
}

\keyword{ utilities }
