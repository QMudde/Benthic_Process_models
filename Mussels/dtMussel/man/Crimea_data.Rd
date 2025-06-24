\name{Crimea_data}
\docType{data}
\alias{Crimea_data}
\alias{BFFMod_data}

\title{
Shell length data of two Mytilus galloprovincialis populations at Sevastopol Bay
}
\description{
The dataset \code{Crimea_data} contains shell length Mytilus galloprovincialis, growing in Sevastopol bay, crimea. Two populations, based on their initial mean shell length (35mm and 60mm), are distinguished and their mean size trajectories are proveded as data.
This dataset is used for the validation of the individual growth model \code{mussel_run_debind()} with the parameters shown in \code{mussel_get_debparms}

Meta information on the \code{Crimea_data}:
  \itemize{
      \item The data originate from the paper by \code{Chelyadina & Popov, 2020}
      \item Measurments were done between december 2018 and march 2020. Populations were cultures in hanging oyster cages at 3m depth near a mussel-oyster culture farm at the opening of Sevastopol Bay, Crimea.
      \item Original paper only provided seawater temperature as environmental variable, additional variables are model data originating from \code{copernicus marine service}
      \item Coordinates, associated with data are -> Lat: 44.620389 ;  Lon : 33.503778
}
  \code{Crimea_data} contains: 
  \itemize{
    \item \code{times as days since jan01 2018} the days of the measurements 
    \item \code{ini_35_size/SE, ini_60_size/SE} The mean and standard error of shell length of the two followed populations expressed in cm. One with intial size of 35mm and the second with initial size of 60mm 
    \item \code{Temperature} Measured temperatuere in dgC at given time 
    \item \code{phyc_mod, dt_mod, ox_mod} Modeled Phytoplankton [mmolC/m3], Detritus [mmolC/m3] and Oxygen [mmolO2/m3] concentrations, extracted from \code{Copernicus marine service} at the given times and coordinates. 
}
}

\usage{
load("./data/Crimea_data.rda")
}

\references{
Chelyadina, N., & Popov, M. (2020). 
Individual variability of growth of the cultivated Mytilus galloprovincialis Lam. in various habitat conditions of the Black sea coastal zone. 
Aquaculture Research, 51(12), 4846-4851.
}

\examples{
## data structure
ncol(Crimea_data)## number of cols
nrow(Crimea_data)## number of rows
## plot data of 35mm population
plot(x=Crimea_data$times, y = Crimea_data$ini_35_size)
points(x=Crimea_data$times, y = Crimea_data$ini_35_size+2*Crimea_data$ini_35_sd)
points(x=Crimea_data$times, y = Crimea_data$ini_35_size-2*Crimea_data$ini_35_sd)
## plot data of 60 mm population
plot(x=Crimea_data$times, y = Crimea_data$ini_60_size)
points(x=Crimea_data$times, y = Crimea_data$ini_60_size+2*Crimea_data$ini_60_sd)
points(x=Crimea_data$times, y = Crimea_data$ini_60_size-2*Crimea_data$ini_60_sd)
}

