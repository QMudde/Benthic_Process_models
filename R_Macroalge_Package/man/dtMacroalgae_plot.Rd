\name{plotting functions}
\alias{plot.dtMacroalgae}

\title{
  Utility functions for plotting dynamic model output.
}
  
\description{
Plots dynamic output with labels on the y-axis. 
}

\usage{
\method{plot}{dtMacroalgae}(x, ..., 
                           select = NULL, which = select, ylab = NULL, 
                           lty = 1, las = 1) 
}

\arguments{

  \item{x }{an output object returned by \link{Macroalgae_run}. 
  }

  \item{which, select}{ The name(s) of the output variables to be plotted versus time.
  }
  \item{ylab}{ The labels of the y-axis. The default uses the units of the variables.
  }
  \item{lty, las}{ The type of the lines used for plotting (default = solid line) and the label orientation.
  }
  \item{...}{ Any other dynamic model output object.
  }
}  

\note{
 The function adds the units as the labels for the y-axis, and then calls the plotting function as defined for \code{deSolve} objects.
}

\value{

 returns nothing
}

\author{
  Karline Soetaert
}


\examples{

# ---------------------------------------------------
# runs to work with
# ---------------------------------------------------
times <- seq(0 : (365*3))

# ----------------------------------
# two runs

cout  <- Macroalgae_run ()
cout2 <- Macroalgae_run (parms = list(pMax = 0.05))

# ----------------------------------
# plot one variable for one run

plot(cout, 
     which = "Total_C",
     xlab  = "days")

# ----------------------------------
# plot several variables for several runs at once

plot(cout, cout2, 
     which = 1:9, 
     las   = 0)   # default for plotting


# ----------------------------------
# both models in same plot

plot(cout, cout2, 
     which = c("STRUCT_C", "RESERVE_C"),
#     obs = obs_calvi,
     mfrow = c(2,2))

plot(cout, cout2, 
     which = c("O2"),
     mfrow = NULL)       # to avoid overruling mfrow


}
\references{
  ... to be added?
}

\keyword{ utilities }
