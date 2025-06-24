\name{get_INI_indiv}
\alias{get_INI_indiv}
\alias{BFFMod get_INI_indiv}
\title{
  Utility functions in BFFMod package for working with fied data or analysing output : \code{get_INI_indiv}
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
 get_INI_indiv(Length,
               AFDW       = NA,
               DW         = NA,
               WW         = NA,
               Rel_E_dens = 1,
               matured = TRUE, 
               parms = list())
 }

\arguments{
      \item{Length}{Shell length of mussel at start of model run}
      \item{AFDW}{inital ash free dry weight in grams, if known}
      \item{DW}{inital dry weight in grams, if known}
      \item{WW}{inital shell free wet weight in grams, if known}
      \item{Rel_E_dens}{inital reserves density relative to the maximum (Reserves/Structure). This can differ between summer and winter}
      \item{matured}{set to FALSE if mussel is too young to be at full maturity}
      \item{parms}{change parameters if applied to distinct species or population with known DEB parameters}
}

\examples{
## simple inputs to construct vector of DEB compartments
  get_INI_indiv(Length= 5, WW = 1.5)
  
## Depending on available data and season, we can manipulate. 
  get_INI_indiv(Length= 3, AFDW = 0.2,
  Rel_E_dens = 0.5,
  matured = FALSE) 
         
}

\keyword{ utilities }
