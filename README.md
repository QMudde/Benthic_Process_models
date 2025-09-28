# Benthic_Process_models
Benthic Process Model code for NECCTON (WP6)

Repository contians code for mechanistic process models of benthic lifeforms. 
Process models are developed for Mussel beds (Mytilus galloprovincialis), Seaweeds (Phyllophora) and Seagrass (Posidonia). Model dynamics are based on detailed, species specific ecophysiologial processes and rates.
Models are formulated as a complex systems of differential equations, they are written in Fortran and can be called in R. Calibration is based on rates provided by literature and the associated calibration datasets. 

In the scope of the NECCTON project, the models are used to predict growth and reproduction potential of shellfish in the Black sea, and predict growth of Phyllophora in competition with Cystoseira seaweed assemblages in the Black sea.
For the Black sea, the models can be one-way coupled with BAMHBI-NEMO, meaning that the models can make use of direct input from BAMHBI-NEMO, but cannot provide feedback. 
Black Sea reanalysis data, obtained from Copernicus Marine, can be used as forcings.

# Instructions & Usage: 
Each finished model is accompanied by a .Rmd vignette to instruct the user about the model dynamics and its usage.


## References 
published data from literatre are used to calibrate the models
### Mussels 
- Karayucel et al., 2010
- Chelyadina & Popov, 2020
- Vural et al., 2015
- Galimany et al., 2013
### Seagrass
- Champenois & Borger, 2021
### Macroalgae 
- Afanasyef et al., 2017
- Blindova & Trishna, 1990
- Pankeeva & Mironova, 2019

#### Dependencies 
The R code structure makes use of several pacakges: 
- library(deSolve) # solving differential equations in Fortran & R
- library(FME) # providing fitting and cost fuctions for calibration of parameters
