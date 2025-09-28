\name{parameters and variables}
\alias{Macroalgae_get_parms}
\alias{Macroalgae_get_states}
\alias{Macroalgae_get_vars}
\alias{Macroalgae_get_forcings}

\title{
  Functions to retrieve parameters, variables and state variables, and forcing functions for the Macroalgae dynamic models
}
  
\description{

  \code{Macroalgae_get_parms} retrieve the parameter names, description, units, default values and values from model solutions.

  \code{Macroalgae_get_forcings} retrieve the forcing function names, description, and units.

  \code{Macroalgae_get_vars, and Macroalgae_get_states} retrieve the names, description, units and mean values for the output variables and the state variables respectively.
}

\usage{
  Macroalgae_get_parms(out, as.vector = FALSE, which = NULL)
  
  Macroalgae_get_forcings(out, as.vector = FALSE, which = NULL)
  
  Macroalgae_get_vars (out, as.vector = FALSE, which = NULL)
  
  Macroalgae_get_states(out, as.vector = FALSE, which = NULL)
}

\arguments{

  \item{out }{an output object returned by \link{Macroalgae_run}. If absent, \code{Macroalgae_get_parms} will return the default (parameter) values. For \code{Macroalgae_get_vars} and \code{Macroalgae_get_states}, an absent argument \code{out} will return just the names, description and units.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter, forcing or variable vector, else a data.frame that also contains the units and description.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}

\value{

\code{Macroalgae_get_vars} and \code{Macroalgae_get_states} return the output or state variables of the solution as a vector or data.frame. If a dynamic run is passed (argument \code{out}), the output, averaged over the run, is also added.
 
\code{Macroalgae_get_parms} return the parameter settings of a run
}

\author{
  Karline Soetaert and Quinten Mudde
}

\examples{

# ---------------------------------------------------
# runs to work with
# ---------------------------------------------------
times <- seq(0 : (365*3))

cout  <- Macroalgae_run(times      = times)     # default parameters
                         
cout2 <- Macroalgae_run(times      = times,
                       parms       = list(pMax = 0.05))

plot(cout, cout2, which = 1:9)

# ---------------------------------------------------
# Definition, units and default parameter settings
# ---------------------------------------------------

Macroalgae_get_parms()

# ---------------------------------------------------
# Parameter settings of a run
# ---------------------------------------------------

Macroalgae_get_parms(cout,  which = "Photosynthesis")
Macroalgae_get_parms(cout2, which = "Photosynthesis")

# ---------------------------------------------------
# Definition, units of output variables
# ---------------------------------------------------

Macroalgae_get_vars()

Macroalgae_get_vars(cout)   # more variables

# ---------------------------------------------------
# Definition, units of state variables
# ---------------------------------------------------

Macroalgae_get_states()
Macroalgae_get_states(cout)   

# To get the mean, we can also do this:

summary(cout)[4,]

}
\references{
  ... to be added?
}

\keyword{ utilities }
