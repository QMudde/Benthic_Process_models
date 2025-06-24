\name{BFFMod_data}
\docType{data}
\alias{BFFMod_data}
\alias{Crimea_data}
\alias{M_gallo_Feeding_data}
\alias{Marm_data}
\alias{Sinop_data}

\title{
Growth and physiology data of Mytilus galloprovincialis at Sinop Bay
}
\description{
The dataset \code{Sinop_data} contains shell length and wet weight growth data of Mytilus galloprovincialis, together with the associated environmental values measuered in the field. This dataset is used for the calibration of the feeding parameters in \code{mussel_get_debparms}

  Meta information on the \code{Sinop_data}:
  \itemize{
      \item The data originate from the paper by \code{Karayucel et al., 2010}
      \item measurments were done betweeh may 2005 and may 2006 on a raft culture, with ropes hanging 13 meters deep near Sinop Bay, Northern Turkey.
      \item Included model data originates from \code{copernicus marine service}
      \item Coordinates, associated with data are -> Lat: 42.002 ;  Lon : 35.154
}
  \code{Sinop_data} contains: 
  \itemize{
    \item \code{time, date, (julian) day} the days of the measurements 
    \item \code{length mean, upper, lower} The mean shell length of the population and the mean plus/minus the standard error. 
    \item \code{total weight (TW) mean, upper, lower} The mean total wet weight (including shell) of the population and the mean plus/minus the standard error.
    \item \code{wet weight (WW) mean, upper, lower} The mean wet weight (exlcuding shell) of the population and the mean plus/minus the standard error.
    \item \code{Temperature, Salinity} Measured physical conditions of seawater at depth of the raft. Temperatue expressed in dgC, salinity in psu. 
    \item \code{Chlorophyll A} Measured cholrophyll A concentration of seawater at time of sampling, expressed in microgram per liter.
    \item \code{POM} Measured total particulate organic matter concentration of seawater at time of sampling, expressed in milligram per liter.
    \item \code{PIM} Measured total particulate inorganic matter concentration of seawater at time of sampling, expressed in milligram per liter.
    \item \code{dens} Populaiton of mussels (nr of individuals/meter of rope)
    \item \code{Oxygen_bamhbi} Modeled dissolved oxygen values, expressed in mmolO2/m3, extracted from a NEMO-BAMHBI coupled hydrodynamic model. Values are associated with the the provided coordinates, depth and date.
}




The dataset \code{Marm_data} contains shell length and total weight growth data of Mytilus galloprovincialis, together with the associated environmental values measuered in the field. This dataset is used for the calibration of the feeding parameters in \code{mussel_get_debparms}

  Meta information on the \code{Marm_data}:
  \itemize{
      \item The data originate from the paper by \code{Vural et al., 2015}
      \item measurments were done between august 2008 and march 2009 hanging net bag populations. Cardak Lagoon, Dardanelles, Marmara Sea avg depth of 5 m.
      \item Oxygen values are a rough estimeate based on measurements by \code{Deniz & Tas 2009} in the Marmara Sea.
      \item Provided values are an average of 3 growth trajectories at different depths (0.5, 5 and 10m)
}
  \code{Marm_data} contains: 
  \itemize{
    \item \code{time, date, (julian) day} the days of the measurements 
    \item \code{length mean, SE, min, max} The mean shell length of the population and the standard error, expressed in mm. Also, the minimum and maximum observed values are proveded for each date. 
    \item \code{total weight mean, SE, min, max} The mean total wet weight (including shell) of the population, their standard error, and the minimum and maximum value expressed in grams.
    \item \code{wet weight estimate (WW) mean, SE} The mean and standard error of wet weight (exlcuding shell), estimed from the TW with a conversion factor.
    \item \code{Temperature, Salinity} Measured physical conditions of seawater at depth of the raft. Temperatue expressed in dgC, salinity in promille. 
    \item \code{Chlorophyll A} Measured cholrophyll A concentration of seawater at time of sampling, expressed in microgram per liter.
    \item \code{POM} Measured total particulate organic matter concentration of seawater at time of sampling, expressed in milligram per liter.
    \item \code{PIM} Measured total particulate inorganic matter concentration of seawater at time of sampling, expressed in milligram per liter.
    \item \code{tot.pop} Populaiton of mussels (total nr of individuals in all 3 netbags)
    \item \code{Oxygen_est} Rough estimation of dissolved oxygen concentration mmolO2/m3, calculated by a sine funciton based on measurements by \code{Deniz & Tas 2009}.
    marm_o2_f<- function(t){(62.5 * sin((2*pi*(t)/365)) +250)}
}

The dataset \code{Sinop_data} contains shell length and wet weight growth data of Mytilus galloprovincialis, together with the associated environmental values measuered in the field. This dataset is used for the calibration of the feeding parameters in \code{mussel_get_debparms}

  Meta information on the \code{Sinop_data}:
  \itemize{
      \item The data originate from the paper by \code{Karayucel et al., 2010}
      \item measurments were done betweeh may 2005 and may 2006 on a raft culture, with ropes hanging 13 meters deep near Sinop Bay, Northern Turkey.
      \item Included model data originates from \code{copernicus marine service}
      \item Coordinates, associated with data are -> Lat: 42.002 ;  Lon : 35.154
}
  \code{Sinop_data} contains: 
  \itemize{
    \item \code{time, date, (julian) day} the days of the measurements 
    \item \code{length mean, upper, lower} The mean shell length of the population and the mean plus/minus the standard error. 
    \item \code{total weight (TW) mean, upper, lower} The mean total wet weight (including shell) of the population and the mean plus/minus the standard error.
    \item \code{wet weight (WW) mean, upper, lower} The mean wet weight (exlcuding shell) of the population and the mean plus/minus the standard error.
    \item \code{Temperature, Salinity} Measured physical conditions of seawater at depth of the raft. Temperatue expressed in dgC, salinity in psu. 
    \item \code{Chlorophyll A} Measured cholrophyll A concentration of seawater at time of sampling, expressed in microgram per liter.
    \item \code{POM} Measured total particulate organic matter concentration of seawater at time of sampling, expressed in milligram per liter.
    \item \code{PIM} Measured total particulate inorganic matter concentration of seawater at time of sampling, expressed in milligram per liter.
    \item \code{dens} Populaiton of mussels (nr of individuals/meter of rope)
    \item \code{Oxygen_bamhbi} Modeled dissolved oxygen values, expressed in mmolO2/m3, extracted from a NEMO-BAMHBI coupled hydrodynamic model. Values are associated with the the provided coordinates, depth and date.
}


\usage{
BFFMod_data(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
