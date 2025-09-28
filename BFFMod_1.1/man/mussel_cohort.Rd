\name{Mussel cohort DEB model}

\alias{mussel_run_debcohort}

\title{
  Model that describes the physiology of a mussel bed, expressed as several cohorts per m2
}
\description{
  \code{mussel_run_debcohort} combines the physiological deb formulation with population dynamics in a number of cohorts. It also includes the interaction with the environment. New cohort birth is determined
  by the collective spawning and larvae dynamics. Cohort mortality is determined by the state of the
  average DEB individual in each cohort. If dynamic environments are enabled, competition for 
  resources drives the population.
  
  The following state variables are described (where N is 1: the maximum number of cohorts)

\itemize{
  \item  RESERVE(N),                  Reserve carbon per individual, [mmolC/ind]
  \item   STRUCT(N)                Structural carbon per individual, [mmolC/ind]
  \item   REPROD(N)  Reproductive carbon (gonads + gametes) per ind, [mmolC/ind]
  \item POPULATION(N)                            Population density,    [ind/m2]
  \item  LARV_DENS                       Density of pelagic larvae, [number/m3]
  \item  LARV_BIOM                       Biomass of pelagic larvae,  [mmolC/m3]
  \item      PHYTO                     Phytoplankton concentration,  [mmolC/m3]
  \item   DETRITUS                          Detritus concentration,  [mmolC/m3]
  \item        SIM  Suspended inorganic matter (spm) concentration,      [g/m3]
  \item        O2                              Oxygen concentration, [mmolO2/m3]
}  
  
}

\usage{
mussel_run_debcohort(parms = list(), times = 0:365, yini = NULL, 
           f_Phyto = 10, f_Temp = 20, f_Detritus = 10, f_Sim = 1, f_O2 = 360,
           max_cohort = 12, dynamic_environment = TRUE)
}
\arguments{

  \item{f_Phyto, f_Temp, f_Detritus, f_Sim, f_O2}{ Forcing functions that represent external environmental conditions. The default is to keep the forcing functions constant, but it is also allowed to pass a data set consisting of (time, value). f_Phyto represents phytoplankton concentrations [mmolC/m3], f_Temp is the temperature [dgC], f_Detritus is detrital organic carbon, [mmolC/m3], f_Sim is suspended inorganic matter, [g/m3], f_O2 is oxygen concentration, [mmolO2/m3].
  }  
  \item{yini}{ Initial condition. If NULL, then the initial condition is set =0 for all mussel tissues (for all cohorts), the density of pelagic larvae, LARV_DENS, is set to 1e6 [ind/m3], and the total biomass of pelagic larvae, LARV_BIOM, is set equal to 1e6*parms["weight_egg"]; the state variables PHYTO, DETRITUS, Sim and O2 are set equal to the value of the corresponding forcing function at the initial times.
  }
  \item{times}{ simulation time specification, in days. If one single value, then the model is is called once, and the derivatives and output variables are returned.
  }
  \item{parms}{ vector or list with parameter values that overrule the default parameter values. See examples for parameters and their default values.
  }
  \item{max_cohort}{ The maximum number of cohorts that are modeled.
  }  
  \item{dynamic_environment}{ If TRUE, then the environmental conditions will be dynamically updated by the mussel dynamics.
  }  
}
\value{
The return values consists of an object of class *dtDyn* and of class *deSolve* as returned by the \link[deSolve]{ode} function of package deSolve. The first column consists of the output \code{times}. Then several output variables are returned. 

The object has several useful attributes, such as the units of all output variables, and the parameter settings. The units are used to label the y-axis from plots.
}

\author{
  Karline Soetaert and Quinten Mudde
}
\details{ 
The parameters for the mussel deb models, and their default values and units can be found via the function \code{mussel_get_debparms()}. 
}
\examples{

# ---------------------------------------------
# ---------------------------------------------
# Simple cohort DEB dynamics
# ---------------------------------------------
# ---------------------------------------------

times <- 0:(12*365)

cout  <- mussel_run_debcohort(times = times)
cout2 <- mussel_run_debcohort(times = times, f_Phyto = 0)
cout3 <- mussel_run_debcohort(times = times, f_Temp = 10)

# water renew, as proxy for phytoplankton growth. increase to increase population. 
cout4  <- mussel_run_debcohort(times = times, 
                               parms = list(water_renew = 10))
                    
plot(cout, cout2, cout3, cout4, 
     which = c(seq(1, 49, by = 4 ), 50:54),
     lty = 1, lwd = 2)

# ---------------------------------------------
# fluctuating forcing functions
# ---------------------------------------------

temp.fluct = function(t) (   6 * sin((2*pi*(t-50)/365)) + 20) 
phy.fluct  = function(t) (12.5 * sin((2*pi*(t   )/365)) + 12.5)

cout5 <- mussel_run_debcohort(
                      times   = times, 
                      f_Temp  = temp.fluct, 
                      f_Phyto = phy.fluct)

plot(cout5, 
     which = c("RESERVE_1", "RESERVE_6", "RESERVE_12",
               "STRUCT_1", "STRUCT_6", "STRUCT_12",
               "REPROD_1", "REPROD_6", "REPROD_12", 
               "Gametes_1", "Gametes_6", "Gametes_12", 
               "Length_1",  "Length_6", "Length_12", 
               "f_Phyto", "f_Temp"))


}


\keyword{ utilities }
