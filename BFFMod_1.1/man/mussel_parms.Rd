\name{parameters and variables}
\alias{mussel_get_debparms}
\alias{mussel_get_debstates}
\alias{mussel_get_debvars}
\alias{mussel_get_debforcs}
\alias{mussel_get_auxiliary}

\title{
  Functions to retrieve parameters, forcing functions, variables and state variables for the mussel dynamic models based on deb}
  
\description{

  \code{mussel_get_debparms} retrieve the parameter names, description, units, default values and values from model solutions obtained with \link{mussel_run_debind} or \link{mussel_run_debcohort}.

  \code{mussel_get_debvars}, \code{mussel_get_debstates} and \code{mussel_get_debforcs} retrieve the names, description, units and mean values for the output variables, the state variables and the forcing functions respectively.
}

\usage{
  mussel_get_debparms (out, as.vector = FALSE, which = NULL)
  mussel_get_debvars  (out, as.vector = FALSE, which = NULL)
  mussel_get_debstates(out, as.vector = FALSE, which = NULL)
  mussel_get_debforcs (out, as.vector = FALSE, which = NULL)
  mussel_get_auxiliary(out, as.vector = FALSE, which = NULL)
}

\arguments{

  \item{out }{an output object returned by \link{mussel_run_debind} or \link{mussel_run_debcohort}. If absent, \code{mussel_get_debparms} will return the default (parameter) values. For \code{mussel_get_debvars} and \code{mussel_get_debstates}, an absent argument \code{out} will return just the names, description and units.
  }
  \item{as.vector }{if \code{TRUE} will return the parameter or variable vector, else a data.frame that also contains the units and description.
  }
  \item{which }{if not \code{NULL}, a vector with names of the variables/parameters to return.
  }
}

\value{

\code{mussel_get_debvars} and \code{mussel_get_debstates} return the output or state variables of the solution as a vector or data.frame. If a dynamic run is passed (argument \code{out}), the output, averaged over the run, is also added.
 
\code{mussel_get_debparms} return the parameter settings of a run
}

\author{
  Karline Soetaert and Quinten Mudde
}

\examples{

# ---------------------------------------------------
# runs to work with
# ---------------------------------------------------
times <- seq(0 : (365*2))

# ----------------------------------
# 3 cohorts, two runs

cout  <- mussel_run_debcohort(
                      max_cohort = 4,     # default parameters
                      times      = times)
cout2 <- mussel_run_debcohort(
                      max_cohort = 4,
                      times      = times,
                      parms      = list(aeff = 0.5))

plot(cout, cout2, which = 1:9)

# ----------------------------------
# single individual model, two runs

dout  <- mussel_run_debind(times = times)

dout2 <- mussel_run_debind(times = times,
                           parms = list(aeff = 0.5))

plot(dout, dout2, which = 1:9)

# ---------------------------------------------------
# Definition, units and default parameter settings
# ---------------------------------------------------

mussel_get_debparms()

# ---------------------------------------------------
# Parameter settings of a run
# ---------------------------------------------------

mussel_get_debparms(cout, which = "aeff")
mussel_get_debparms(cout2, which = "aeff")

# ---------------------------------------------------
# Definition, units of output variables
# ---------------------------------------------------

mussel_get_debvars()

mussel_get_debvars(dout)
mussel_get_debvars(cout)   # more variables

# ---------------------------------------------------
# Definition, units of state variables
# ---------------------------------------------------

mussel_get_debstates()
mussel_get_debstates(dout)
mussel_get_debstates(cout)   # more state variables

# Note: the mean is averaged over all cohorts. 
# To get the mean for each cohort:

summary(cout)[4,]

# ---------------------------------------------------
# Definition, units of forcing functions
# ---------------------------------------------------

mussel_get_debforcs()

mussel_get_debforcs(dout)
mussel_get_debforcs(cout)   


}


\keyword{ utilities }
