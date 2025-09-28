\name{Marm_data}
\docType{data}
\alias{Marm_data}
\alias{BFFMod_data}

\title{
Growth and weight data of Mytilus galloprovincialis at Cardak Lagoon, Dardanelles
}
\description{
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
}

\usage{
load("./data/Marm_data.rda")
}

\references{
Vural, P., Yildiz, H., & Acarli, S. (2015).
Growth and survival performances of Mediterranean mussel (Mytilus galloprovincialis, Lamarck, 1819) on different depths in Cardak lagoon, Dardanelles.
Marine Science and Technology Bulletin, 4(1), 7-12.

Deniz, N., & Taş, S. (2009). 
Seasonal variations in the phytoplankton community in the north-eastern Sea of Marmara and a species list.
Journal of the Marine Biological Association of the United Kingdom, 89(2), 269-276.
}

\examples{
## data structure
ncol(Marm_data)## number of cols
nrow(Marm_data)## number of rows
## plot data on shell length
plot(x= Marm_data$day, y = Marm_data$length.mm)
points(x= Marm_data$day, y = Marm_data$length.mm + 2*Marm_data$length.mm.sd)
points(x= Marm_data$day, y = Marm_data$length.mm - 2*Marm_data$length.mm.sd)
## plot data on wet weights
plot(x= Marm_data$day, y = Marm_data$tot.ww.g)
points(x= Marm_data$day, y = Marm_data$tot.ww.g + 2*Marm_data$tot.ww.g.sd)
points(x= Marm_data$day, y = Marm_data$tot.ww.g - 2*Marm_data$tot.ww.g.sd)
}

