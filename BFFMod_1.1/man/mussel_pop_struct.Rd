\name{mussel_pop_struct}
\alias{mussel_pop_struct}
\alias{BFFMod mussel_pop_struct}
\title{
  Utility functions in BFFMod package for working with fied data or analysing output : \code{mussel_pop_struct}
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
 mussel_pop_struct(model, 
                   age_classes = c(0,1,2,3), 
                   length_classes=c(0,1,2.5,5))
}

\arguments{
      \item{model}{model output from \code{mussel_run_debcohort}}
      \item{age_classes}{vector of age class cutoffs (must start with 0)}
      \item{length_classes}{vector of length class cutoffs (must start with 0)}
}


\examples{
## get cohort model output of eg 6 years
out <- mussel_run_debcohort(times = 0:(6*365))

## the population structure is given by cohort number, but not age or length.
## we can quickly transform to get population structure. 

out_struct <- mussel_pop_struct(out)

## custom structure can be required.
## just 3 age classes (between 0-1, 1-2, 2-older) 
## and 2 size classes (between 0-1.5 cm, 1.5-longer)
## keep in mind classes must *start* at 0

out_struct <- mussel_pop_struct(model = out,
                                age_classes = c(0,1,2),
                                length_classes = c(0,1.5))
         
}

\keyword{ utilities }
