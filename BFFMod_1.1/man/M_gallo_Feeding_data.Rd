\name{M_gallo_Feeding_data}
\docType{data}
\alias{M_gallo_Feeding_data}
\alias{BFFMod_data}

\title{
Phydiological feeding rate data of Mytilus galloprovincialis clearance and filtration, Alfacs Bay, Spain
}
\description{
The dataset \code{M_gallo_Feeding_data} contains data on physiological feeding rates of Mytilus galloprovincialis. The data were measured on a spanish population in Alfacs Bay, and have been used to calibrate order of magnitude levels for clearance and filtration rates in \code{mussel_get_debparms}

  Meta information on the \code{M_gallo_Feeding_data}:
  \itemize{
      \item The data originate from the paper by \code{Galimany et al., 2015}
      \item Mesurements are from experimental setup fed with bay water. which were carried out in November 2006 (4) and February (4), April (3) and July 2007(3).
      \item The original data contains POM, which is roughly separated into Detritus and Phytoplankton. The expected ratio between Detritus and Phytoplankton in winter (Nov&Feb), spring (April) and summer (July) were based on values published by \code{Delgado 1987}
}

  \code{M_gallo_Feeding_data} contains: 
  \itemize{
    \item \code{index} number/index of the experiment. The enries are split into 4 groups, whith all the same average mussel size, and same average temperature but different particle concentrations. November 2006 (1-4), February 2007 (5-8), April 2007 (9-11), July 2007 (10-12) 
    \item \code{Temp} Average seawater temperature in dgC per experiment group.
    \item \code{SIM} Suspended inorganic matter particle concentration in gram per liter.
    \item \code{POM, PhyC, DT} Particular organic matter (POM) expressed in gram/liter. POM values were manually seperated, witch ratios based on \code{Delgado 1987}, into phytoplankton [mmolC/m3] and detritus [mmolC/m3]. the PhyC and DT values are not from the original data. 
    \item \code{Length, DW} Average seawater Length and dry weight in cm and grams per experiment group.
    \item \code{Clear, Clear_sd} Clearance rate means and their standard error, expressed in m3/day/ind.
    \item \code{Filt, Filt_sd} Filtration rate means and their standard error, expressed in g/day/ind, this value contains both organic and inorganic particles.
    \item \code{OrgFilt, OrgFilt_sd} Organic filtration rate means and their standard error, expressed in mmolC/day/ind, this value contains both organic and inorganic particles.
    \item \code{Reject, Reject_sd} Mean and SE of percentage of total particles rejected trough pseudofeaces. 
    \item \code{OrgIng, OrgIng_sd} Mean and SE of ingested organic particles into gut system in mmolC/day/ind 
    \item \code{OrgIng, OrgIng_sd} Mean and SE of ingested organic particles into gut system in mmolC/day/ind
    \item \code{Assim, Assim_sd} Mean and SE of carbon assimilation rate. 
    \item \code{Clear_w, OrgFilt_w,Reject_w,OrgIng_w,Assim_w} Weight factor of physiological rates that is used for fitting purposes. This is defined by the mean value of the timeseries divided by the SE of each entry.  
}
}

\usage{
load("./data/M_gallo_Feeding_data.rda")
}

\references{
Galimany, E., Ramón, M., & Ibarrola, I. (2011).
Feeding behavior of the mussel Mytilus galloprovincialis (L.) in a Mediterranean estuary: a field study.
Aquaculture, 314(1-4), 236-243..

Delgado, M. (1987).
Fitoplancton de las bahías del delta del Ebro.
}

\examples{
## data structure
ncol(M_gallo_Feeding_data)## number of cols
nrow(M_gallo_Feeding_data)## number of rows

## physiological rates are different trough the year
plot(x= M_gallo_Feeding_data$day, y = M_gallo_Feeding_data$Clear)
points(x= Marm_data$index,
y = M_gallo_Feeding_data$Clear + 2*M_gallo_Feeding_data$Clear_sd)
points(x= Marm_data$index,
y = M_gallo_Feeding_data$Clear.mm - 2*M_gallo_Feeding_data$Clear_sd)

## plot data on wet weights
plot(x= M_gallo_Feeding_data$day, y = M_gallo_Feeding_data$OrgFilt)
points(x= Marm_data$index, 
y = M_gallo_Feeding_data$Clear - 2*M_gallo_Feeding_data$OrgFilt_sd)
points(x= Marm_data$index,
y = M_gallo_Feeding_data$Clear - 2*M_gallo_Feeding_data$OrgFilt_sd)
}

