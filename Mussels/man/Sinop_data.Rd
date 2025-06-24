\name{Sinop_data}
\docType{data}
\alias{Sinop_data}
\alias{BFFMod_data}

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
}

\usage{
load("./data/Sinop_data.rda")
}

\references{
Karayücel, S., Çelik, M. Y., Karayücel, İ., & Erik, G. (2010).
Growth and production of raft cultivated Mediterranean mussel (Mytilus galloprovincialis Lamarck, 1819) in Sinop, Black sea.
Turkish Journal of Fisheries and Aquatic Sciences, 10(1).
}

\examples{
## data structure
ncol(Sinop_data)## number of cols
nrow(Sinop_data)## number of rows
## plot data on shell length
plot(x=Sinop_data$day, y = Sinop_data$length_mean)
points(x=Sinop_data$day, y = Sinop_data$length_upper)
points(x=Sinop_data$day, y =Sinop_data$length_lower)
## plot data on wet weights
plot(x=Sinop_data$day, y = Sinop_data$WW_mean)
points(x=Sinop_data$day, y = Sinop_data$WW_upper)
points(x=Sinop_data$day, y =Sinop_data$WW_lower)
}

