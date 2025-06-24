\name{observed data and forcinngs}
\alias{obs_calvi}
\alias{forcs_calvi}

\title{
  Observed data and forcing functions from the Bay of Calvi (Corsica, France) to run the Posidonia model.
}
  
\description{

  \code{forcs_calvi} contains forcing function data, obtained from EMWCF and from U Liege. This comprises:
  \itemize{
  \item{f_Temperature}, water temperature in dgC (U Liege), 2007-2016
  \item{f_Light}, Surface direct short-wave (solar) radiation in W/m2 (EMWCF), 2007-2016.
  \item{f_Wind}, wind speed in m/s (EMWCF), 2007-2016.
  }
  
  \code{obs_calvi} contains observed data, obtained from ??? for the year 2007. This comprises:
  \itemize{
  \item{LEAVES, ROOTS, RHIZOMES}, obsrved leaf, root and rhizome biomass in molC/m2.
  \item{OXYGEN}, water-column oxygen concentration, molO2/m3.
  }  
  
  The data.frames contain 3 columns:
  \itemize{
  \item{variables},  the name of the forcing or data variable
  \item{day},  the day nr (day 1 = 01-01-2007)
  \item{value},  the value
  }  
}

\usage{
  data(obs_calvi)
  data(forcs_calvi)
}


\author{
  Karline Soetaert and Quinten Mudde
}

\examples{

## ==================
## Observed data
## ==================

par(mfrow = c(2,2))
oxygen <- subset(obs_calvi, 
                 subset = variable == "OXYGEN",
                 select = c(day, value))
plot(oxygen, ylab = "mol O2/m3", main = "OXYGEN", type = "l")                 


leaves <- subset(obs_calvi, 
                 subset = variable == "LEAVES",
                 select = c(day, value))
plot(leaves, ylab = "mol C/m2", main = "LEAVES", type = "l")                 

roots <- subset(obs_calvi, 
                 subset = variable == "ROOTS",
                 select = c(day, value))
plot(roots, ylab = "mol C/m2", main = "ROOTS", type = "l")                 

rhizomes <- subset(obs_calvi, 
                 subset = variable == "RHIZOMES",
                 select = c(day, value))
plot(rhizomes, ylab = "mol C/m2", main = "RHIZOMES", type = "l")                 



## ==================
## Forcing data
## ==================

par(mfrow = c(2,2))
temp <- subset(forcs_calvi, 
                 subset = variable == "f_Temperature",
                 select = c(day, value))
plot(temp, ylab = "dgC", main = "f_Temperature", type = "l")                 


light <- subset(forcs_calvi, 
                 subset = variable == "f_Light",
                 select = c(day, value))
plot(light, ylab = "W/m2", main = "f_Light", type = "l")                 

wind <- subset(forcs_calvi, 
                 subset = variable == "f_Wind",
                 select = c(day, value))
plot(wind, ylab = "dgC", main = "f_Wind", type = "l")                 


}
\references{
  ... to be added?
}

\keyword{ data }
