\name{make_ff}
\alias{make_ff}
\alias{BFFMod make_ff}
\title{
  Utility functions in BFFMod package for working with fied data or analysing output : \code{make_ff}
}
\description{
Additional functions of the BFFMod package are described: 
\itemize{
  \item \code{get_INI_indiv}: Provides initial DEB conditions of individual based on length & biomass
  \item \code{get_INI_coh}: Creates vector of initial conditions for \code{mussel_run_debcohort} based on 
  condions of the first cohort
  \item \code{make_ff}: Creates a forcing function based on data, sets the correct format and units
  \item \code{mussel_pop_struct}: Restructures the \code{mussel_run_debcohort} output to the desired     population structure
}}

\usage{
        make_ff(day, values, recalc)
}

\arguments{
      \item{day}{vector of days from time series data}
      \item{values}{vector of values from time series data}
      \item{recalc}{recalculation of units, default is no recaculation. 
      supported units are:
      phytoplankton in grams: "phyc_g", 
      particulate organic matter in grams: "pom_g", 
      chlorophyll a in micrograms: "chl_mugl"}
}

\examples{
  ## function recalculates certian units to units used in the model (mmolC)
      make_ff(day = c(1,2,3,4,5), 
              values = c(2,3,4,5,3), 
              recalc = "chl_mugl")
              
  ## input days and values do not need to be sequential
  ## by default no recalculation are done
  
      make_ff(day = c(10,30,100,1000), 
              values = c(2,3,4,5,3))        

}

\keyword{ utilities }
