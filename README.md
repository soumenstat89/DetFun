[README.md](https://github.com/soumenstat89/DetFun/files/7032571/README.md)
# DetFun

####################################################################################################
#----#----#----#----#----  SINGLE-SEASON SCR: VARIOUS DETECTION FUNCTIONS  ----#----#----#----#----#
####################################################################################################
# LOADING R LIBRARY PACKAGES

require(rgeos)                     # Import the geospatial analysis libraries
require(rgdal)                      # Import the spatial input/ouput libraries
require(raster)                     # Import the raster spatial package
require(coda)                       # Import the MCMC diagnostic tools
require(nimble)                    # Import the NIMBLE subroutines
require(nimbleSCR)             #  Importing custom NIMBLE functions for SCR analysis 
require(abind)                     # Import the library for manipulating multidimensional arrays


# THIS SCRIPT IS FOR RUNNING THE MODEL BY SOURCING THE R SCRIPT "DetFun.R".
source("DetFun.R")

#THE FOLLOWING R FUNCTIONS ARE LOADED BY SOURCING "DetFun.R"

# SimulateDetection_detfun( ): Function to simulate individual detections within an SCR framework, under the influence of a specified detection function.
# UTMToGrid( ): Function to convert UTM coordinates to GRID coordinates for both habitat and detectors
# MakeAugmentation( ): increases the dimensions of an object along the individual and/or the year dimension. It returns a Vector or Matrix object with the expanded object.
# MakeInitXY( ): Function to set initial values of XY coordinates of ACS
# dbernPP( ): Probability Density of a Bernoulli Point Process 
# dbin_LESSCachedAllSparse_detfun_hn( ): Function to create a NIMBLE custom distribution for faster SCR model runs with respect to half-normal detection function.
# dbin_LESSCachedAllSparse_detfun_hnp( ): Function to create a NIMBLE custom distribution for faster SCR model runs with respect to half-normal plateau detection function.
# dbin_LESSCachedAllSparse_detfun_ex( ): Function to create a NIMBLE custom distribution for faster SCR model runs with respect to exponential detection function.
# HRA_HNP_nimble_array( ): Function to compute home range radius and area for a specified circular quantile with respect to half-normal plateau detection function.
# HRA_EX_nimble_array( ): Function to compute home range radius and area for a specified circular quantile with respect to exponential detection function.
# HRA_detfun( ): Function to compute home range radius and area for a specified circular quantile with respect to any given detection function.
# simulateFromPosterior( ): Function to draw posterior replicates of the SCR data for a given MCMC sample of parameters.
# pvaluefn( ): Function to compute Bayesian p-value for a given discrepancy metric.
# runDetFun( ): Master function for simulation of SCR data set, run MCMC, compute home range radius and area, Bayesian p-value for the simulated SCR data set with given detection function for simulation and model fitting

myVars <- list(
  HABITAT =    list(extent     = c(20), # Extent of habitat in each side
                    resolution = 1,  # Habitat cell resolution
                    buffer     = 4.5),    # Buffer width 
  DETECTORS =  list(resolution = 1),  # Detector resolution     
  POPULATION = list(N = 200,               # true population size
                    M = 400),              # size of the augmented population
  # From below, only use the list of parameter values that corresponds to the detection function to be used for simulation
  #-- Half normal
  HN = list(p0  = 0.3,
            sigma = 1.5), 
  #-- Half normal plateau
  HNP = list(p0  = 0.25,
             sigma = 1,
             wd = 1.5), 
  #-- Donut
  DN = list(p0 = 0.25, 
            sigma = 1, 
            sigma.b = 1.5,
            wd = 1.5), 
  # Asymmetric logistic 
  AL = list(p0 = 0.3, 
            sigma = 2,
            slope = 5, 
            slopecon = 1), 
  # Exponential 
  EX = list(p0 = 0.3, 
            sigma = 1.5),
  #Bimodal 
  BI = list(p0 = 0.25, 
            sigma = 0.5, 
            p0.b = 0.15, 
            sigma.b = 1, 
            wd = 2) 
)

#============= POSSIBLE CHOICES FOR DETECTION FUNCTIONS ==========
# detfun.names.sim <- c("HN", "EX", "HNP", "AL", "DN", "BI") # FOR SIMULATIONS
# detfun.names.fit <- c("HN", "HNP", "EX") # FOR MODEL FITTING


#============= RUN THE MCMC AND OBTAIN ESTIMATES ==========
runDetFun(
  detfunsim = "AL"             # Detection function used for simulation of SCR data
  , detfunfit = "EX"              # Detection function used for SCR model fitting
  , burnin = 1                      # MCMC burnin period
  , niter =   100                   # No. of MCMC iterations
  , nchains =  2                   # No. of MCMC chains
  , pr = 0.95                        # The quantile probability in order to compute the home range radius 
  , WD = "ABC/DEF/GHI"   # Working Directory to save the figures and model output.
  , myVars = myVars          # Details of Habitat design and detector array details
  , plot.save.check = T)      # Indicator to save the plots in working directory WD 



