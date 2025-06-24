\name{plotting functions}
\alias{plot.dtDyn}
\alias{matplot.mussel}

\title{
  Utility functions for plotting dynamic model output.
}
  
\description{
Plots dynamic output with labels on the y-axis. 
}

\usage{
\method{plot}{dtDyn}(x, ..., select = NULL, which = select, ylab = NULL, 
                     lty = 1, las = 1) 

matplot.mussel(x, which = NULL, 
               ylab = NULL, lty = 1, las = 1, 
               legend = list(x = "topleft", cex=0.5), ...)
}

\arguments{

  \item{x }{an output object returned by \link{mussel_run_debind} or \link{mussel_run_debcohort}. 
  }

  \item{which, select}{ The name(s) of the output variables to be plotted versus time.
  }
  \item{ylab}{ The labels of the y-axis. The default uses the units of the variables.
  }
  \item{lty, las}{ The type of the lines used for plotting (default = solid line) and the label orientation.
  }
  \item{legend}{ Specifications of the legend to be written for matplot.mussel.
  If \code{NULL}, no legend will be written.
  }
  
  \item{...}{ Any other dynamic model output object.
  }
}  

\note{
 The \code{plot} function adds the units as the labels for the y-axis, and then calls the plotting function as defined for \code{deSolve} objects.
 
 The \code{matplot.mussel} function allows to plot several variables in one plot.
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
# 3 cohorts, two runs

cout  <- mussel_run_debcohort(max_cohort = 4,     # default parameters
                              times      = times)
cout2 <- mussel_run_debcohort(max_cohort = 4,
                              times      = times,
                              parms      = list(aeff = 0.9))

# ----------------------------------
# plot one variable for one run

plot(cout, 
     which = "Total_C_m2",
     xlab  = "days")

# ----------------------------------
# plot several variables for several runs at once

plot(cout, cout2, 
     which = 1:9, 
     las   = 0)   # default for plotting

# ----------------------------------
# single individual model, two runs

dout  <- mussel_run_debind(times = times)

dout2 <- mussel_run_debind(times = times,
                           parms = list(aeff = 0.9))

plot(dout, dout2, which = 1:9) 

# ----------------------------------
# both models in same plot

plot(dout, dout2, 
     which = c("RESERVE", "STRUCT"),
     mfrow = c(2,2))

plot(cout, cout2, 
     which = c("RESERVE_1", "STRUCT_1"),
     mfrow = NULL)       # to avoid overruling mfrow

# ----------------------------------
# using matplot.mussel

matplot.mussel(dout, 
               which = 1:3, 
               lwd = 2, 
               legend = list(cex = 1, lwd  = 2))

matplot.mussel(cout, 
               which = "RESERVE", 
               lwd = 2)
               

}


\keyword{ utilities }
