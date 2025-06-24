\name{parameters and variables}
\alias{Posidonia_get_parms}
\alias{Posidonia_get_states}
\alias{Posidonia_get_vars}
\alias{Posidonia_get_forcings}

\title{
  Functions to retrieve parameters, variables and state variables, and forcing functions for the Posidonia dynamic models
}
  
\description{

  \code{Posidonia_get_parms} retrieve the parameter names, description, units, default values and values from model solutions.

  \code{Posidonia_get_forcings} retrieve the forcing function names, description, and units.

  \code{Posidonia_get_vars, and Posidonia_get_states} retrieve the names, description, units and mean values for the output variables and the state variables respectively.
}

\usage{
  Posidonia_get_parms(out, as.vector = FALSE, which = NULL)
  
  Posidonia_get_forcings(out, as.vector = FALSE, which = NULL)
  
  Posidonia_get_vars (out, as.vector = FALSE, which = NULL)
  
  Posidonia_get_states(out, as.vector = FALSE, which = NULL)
}

\arguments{

  \item{out }{an output object returned by \link{Posidonia_run}. If absent, \code{Posidonia_get_parms} will return the default (parameter) values. For \code{Posidonia_get_vars} and \code{Posidonia_get_states}, an absent argument \code{out} will return just the names, description and units.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter, forcing or variable vector, else a data.frame that also contains the units and description.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}

\value{

\code{Posidonia_get_vars} and \code{Posidonia_get_states} return the output or state variables of the solution as a vector or data.frame. If a dynamic run is passed (argument \code{out}), the output, averaged over the run, is also added.
 
\code{Posidonia_get_parms} return the parameter settings of a run
}

\author{
  Karline Soetaert and Quinten Mudde
}

\examples{

# ---------------------------------------------------
# runs to work with
# ---------------------------------------------------
times <- seq(0 : (365*3))

# ----------------------------------
# 3 cohorts, two runs

cout  <- Posidonia_run(times      = times)     # default parameters
                         
cout2 <- Posidonia_run(times      = times,
                       parms      = list(maxPhotosynthesis = 0.05))

plot(cout, cout2, which = 1:9)

# ---------------------------------------------------
# Definition, units and default parameter settings
# ---------------------------------------------------

Posidonia_get_parms()

# ---------------------------------------------------
# Parameter settings of a run
# ---------------------------------------------------

Posidonia_get_parms(cout, which = "maxPhotosynthesis")
Posidonia_get_parms(cout2, which = "maxPhotosynthesis")

# ---------------------------------------------------
# Definition, units of output variables
# ---------------------------------------------------

Posidonia_get_vars()

Posidonia_get_vars(cout)   # more variables

# ---------------------------------------------------
# Definition, units of state variables
# ---------------------------------------------------

Posidonia_get_states()
Posidonia_get_states(cout)   

# Note: the mean is averaged over all cohorts. 
# To get the mean for each cohort:

summary(cout)[4,]

}
\references{
  ... to be added?
}

\keyword{ utilities }
