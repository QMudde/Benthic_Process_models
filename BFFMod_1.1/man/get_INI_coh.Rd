\name{get_INI_coh}
\alias{get_INI_coh}
\alias{BFFMod get_INI_coh}
\title{
  Utility functions in BFFMod package for working with fied data or analysing output : \code{get_INI_coh}
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
 get_INI_coh(Indiv_ini, Dens_ini,
             max_coh = 12,
             PhyC_ini=10, DT_ini=10, Sim_ini=1, O2_ini=360)
 }

\arguments{
      \item{Indiv_ini}{vector of c(RESERVE = X, STRUCT = Y, REPROD = Z).
      usually direct output from \code{get_INI_indiv}}
      \item{Dens_ini}{inital density of first cohort}
      \item{max_coh}{maximum number of cohorts in required model run.
      this must be an integer between 2 and 20}
      \item{PhyC_ini}{Phytoplankton Concentration at first timestep}
      \item{DT_ini}{Organic detritus Concentration at first timestep}
      \item{Sim_ini}{Suspended inorganic matter concentration at first timestep}
      \item{O2_ini}{Oxygen concentration at first timestep}
}

\examples{
## simple inputs to construct vector for population with cohorts. 
## the first cohort will be filled according to the input
  get_INI_coh(Indiv_ini = c(RESERVE = 1, STRUCT = 1, REPROD = 0),
  Dens_ini = 100)
  
## First timestep of dynanimc forcings can be provided to avoid rapid 
## initial changes in potentially short model runs 
## set max cohort value to allign with future model runs
  get_INI_coh(Indiv_ini = c(RESERVE = 1, STRUCT = 1, REPROD = 0),
  max_coh = 8,
  Dens_ini = 100, PhyC_ini = 2,DT_ini = 30,O2_ini=100)

## it is easy to provide output from \code{get_INI_indiv}
  get_INI_coh(Indiv_ini =   get_INI_indiv(Length= 3, WW = 1.5),
  Dens_ini = 100)
}

\keyword{ utilities }
