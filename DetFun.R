# THIS SCRIPT CONTAINS ALL THE FUNCTIONS AND MODEL CODES THAT ARE USED DURING THE ANALYSIS.

library(nimble)                    # Import the NIMBLE subroutines
library(rgeos)                     # Import the geospatial analysis libraries
library(rgdal)                     # Import the spatial input/ouput libraries
library(raster)                    # Import the raster spatial package
library(coda)                      # Import the MCMC diagnostic tools
library(nimble)                    # Import the NIMBLE subroutines
library(nimbleSCR)
library(abind)                     # Import the library for manipulating multidimensional arrays

###############################################################
##################### SimulateDetection_detfun ################
###############################################################
#' @title Function to simulate individual detections within an SCR framework, under the influence of 
#' different detection functions.
#'
#' @description
#' \code{SimulateDetection_detfun} returns a list object with \code{} 
#' 
#' @param params \code{Numeric} vector denoting the parameter values of specified detection function
#' @param AC.sp A \code{SpatialPointsDataFrame} object with individual activity center locations.
#' @param detector.sp  A \code{SpatialPointsDataFrame} object with detectors' locations.
#' @param detfun.type A \code{Character} value denoting the type of detection function in the SCR model.
#'  Avalaible choices are: "HN", "EX", "HNP", "AL", "DN", "BI".
#' @param seed A seed number for results to be reproducible (default is NULL)
#' @param plot.check A \code{logical} for whether (\code{TRUE}) or not (\code{FALSE}) plots are to be generated during simulations.

#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' det.sim<-SimulateDetection_detfun(params= c(0.6,1200),AC.sp=AC.sp,detector.sp=detector.sp,
#' detfun.type = "HN", seed = NULL, plot.check = TRUE)

SimulateDetection_detfun <- function( params = NULL
                                         , AC.sp
                                         , detector.sp
                                         , detfun.type = "HN" 
                                         , seed = NULL
                                         , plot.check = T
)
{
  
  ##------------------------------------------------------------------------------------------------------------------- 
  ## ==== CLEAN & SET-UP THE INPUT DATA ====
  p0 <- params[1]
  sigma <- params[2]
  projection(AC.sp) <- projection(detector.sp)
  
  myP0 <- p0                                  ## For plotting purpose
  mySIGMA <- sigma                                              ## For plotting purpose
  n.trials <- matrix(1, length(AC.sp), length(detector.sp), byrow = TRUE) 
  
  #--- Create separate detectors and sub-detectors SpatialPointsDataFrame when needed
  #--- Subset detector.sp to main detectors only
  sub.detector.sp <- detector.sp
  temp <- unique(detector.sp@data[ ,c("main.cell.id","main.cell.x","main.cell.y")])
  temp <- temp[order(temp$main.cell.id), ]
  detector.sp <- SpatialPointsDataFrame(temp[ ,c("main.cell.x","main.cell.y")], data = data.frame(temp), proj4string = CRS(projection(detector.sp)))
  
  #--- Occasionally rownames pose a problem; fix now:
  dimnames(detector.sp@data)[[1]] <- dimnames(detector.sp@coords)[[1]] <- 1:length(detector.sp)
  dimnames(AC.sp@data)[[1]] <- dimnames(AC.sp@coords)[[1]] <- 1:length(AC.sp)
  
  sigma <- matrix(sigma, length(AC.sp), length(detector.sp), byrow = FALSE)
  individual.sigma <- sigma[,1]
  p0 <- matrix(p0, length(AC.sp), length(detector.sp))
  individual.p0 <- p0
  
  D <- gDistance(detector.sp, AC.sp, byid=TRUE)
  
  ## HALF NORMAL 
  if(detfun.type == "HN"){ P <- p0 * exp(- D*D / (2.0 * sigma * sigma)) }
  
  ## HALF NORMAL PLATEAU
  if(detfun.type == "HNP"){
    wd <- params[3]
    dens <- as.matrix(D) #rep(NA, length(D))
    dens[D<=wd] <- 1
    dens[D>wd] <- exp(-(D[D>wd]-wd)*(D[D>wd]-wd)/(2*sigma[D>wd]*sigma[D>wd]))
    P <- p0 * dens
  }
  
  # DONUT
  if(detfun.type == "DN"){
    sigma.b <- params[3]
    wd <- params[4]
    dens <- as.matrix(D) # rep(NA, length(D))
    dens[D<=wd] <- exp(-(D[D<=wd]-wd)*(D[D<=wd]-wd)/(2.0*sigma.b*sigma.b))
    dens[D>wd] <- exp(-(D[D>wd]-wd)*(D[D>wd]-wd)/(2.0*sigma[D>wd]*sigma[D>wd])) # sigma is a matrix with same dim as D
    P <- p0 * dens
  }
  
  ## Asymmetrical logistic 
  if(detfun.type == "AL"){
    slope <- params[3]
    slopecon <- params[4]
    Anuslope <- (2*abs(slope)*slopecon)/(1+slopecon)
    fx <- 1/ (1+(D/sigma)^Anuslope)
    den <- 1+fx*((D/sigma)^(slope))+(1-fx)*((D/sigma)^(slope*slopecon))
    P <- p0/den
  }
  
  ## EXPONENTIAL 
  if(detfun.type == "EX"){P <- p0 * exp(-D/sigma)}
  
  ## BIMODAL - DETECTION KERNEL WITH A TWO PEAKS 
  if(detfun.type == "BI"){
    p0.b <- params[3]
    sigma.b <- params[4]
    wd <- params[5]
    densa <- p0 *  exp(-D*D/(2.0*sigma*sigma)) 
    densb <-  p0.b * exp(-(D-wd)*(D-wd)/(2.0*sigma.b*sigma.b))
    P <- densa + densb
  }
  
  
  ##-------------------------------------------------------------------------------------------------------------------
  ## ==== BERNOULLI or BINOMIAL DETECTION PROCESS ====
  #-- Set the seed for replication
  y.seed <- P
  if(!is.null(seed)){set.seed(seed)}
  y.seed[] <- sample(1000000, length(n.trials), replace = TRUE)
  #-- Sample number of individual detections per main detector
  temp <- abind(n.trials, P, y.seed, along = 3)
  y <- apply(temp, c(1,2), function(x){
    set.seed(x[3])
    rbinom(1, x[1], x[2])})
  # }#if              
  
  ##------------------------------------------------------------------------------------------------------------------- 
  ## ==== OUTPUT LIST ====
  #--- Fix the names
  if(!is.vector(y)){dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])}  #---DO NOT REMOVE! :)
  dimnames(AC.sp@coords) <- list(c(1:length(AC.sp)), c("x","y"))
  
  #--- Slim to individuals detected
  y.all <- y
  if(length(dim(y)) > 1){
    detected <- apply(y, 1, max) > 0
    if(!any(detected)){detected <- rep(TRUE, dim(y)[1])}
    y <- y[detected, ]
  }
  if(length(dim(y)) == 0){
    detected <- which(y < dim(P)[2])
    if(!any(detected)){detected <- rep(TRUE, length(y))}
    y <- y[detected]
  }
  
  #--- Output list
  out <- list( y = y        
               , y.all = y.all
               , P = P
               , D = D
               , p0 = p0
               , sigma = sigma
               , n.trials = n.trials[1, ])
  
  ##------------------------------------------------------------------------------------------------------------------- 
  ## ==== PLOT DETECTIONS ====
  if(plot.check){
    
    plot(AC.sp)
    plot(detector.sp, col = "gray40", pch = 19, cex = 0.6, add = TRUE)
    plot(AC.sp, add = TRUE)
    
    plot(AC.sp[detected, ], col = "red", pch = 19, add = TRUE)
    
    lapply(which(detected), function(x){
      this.row <- y.all[x, ]
      this.det <- detector.sp[this.row>0, ]
      if(length(this.det)>0){
        segments( coordinates(this.det)[,1], coordinates(this.det)[,2]
                  , coordinates(AC.sp[x,])[,1], coordinates(AC.sp[x,])[,2]
                  , col = "red")
      }})
  }
  
  
  # }#if
  ##------------------------------------------------------------------------------------------------------------------- 
  return(out)
}


###############################################################
##################### UTMToGrid ###############################
###############################################################
#' @title Function to convert UTM coordinates to GRID coordinates for both habitat and detectors
#'
#' @description
#' \code{UTMToGrid} returns a list object with \code{} 
#' 
#' @param data.sp  A \code{SpatialPointsDataFrame} object with detectors' locations.
#' @param grid.sp A \code{SpatialPointsDataFrame} object with habitat grid cell locations.

#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' scaled<-UTMToGrid(grid.sp=habitat.sp,detector.sp=detector.sp)

UTMToGrid <- function(data.sp = NULL,
                      grid.sp = NULL){
  # PREPARE THE DATA
  grid.xy <- as.array(coordinates(grid.sp))
  dimnames(grid.xy) <- list(1:length(grid.sp), c("x","y"))
  data.xy <- as.array(coordinates(data.sp))
  dimnames(data.xy) <- list(1:length(data.sp), c("x","y"))
  # CALCULATE THE RESOLUTION
  resolution <- min(diff(unique(sort(grid.xy[ ,"x"]))))#---assumes square grid cells and utm projection (units in meters or km; not latlong!)
  ## obtain x and y min
  start0.y <- max(grid.xy[ ,"y"]) + resolution/2 #---because we are moving from top to bottom
  start0.x <- min(grid.xy[ ,"x"]) - resolution/2 #---because we are moving from left to right
  ##---- Re-projecting the grid cell centers
  grid.scaled.xy <- grid.xy
  grid.scaled.xy[ ,"y"] <- (start0.y - grid.xy[ ,"y"])/resolution
  grid.scaled.xy[ ,"x"] <- (grid.xy[ ,"x"] - start0.x)/resolution
  ##---- REPROJECTING THE DATA
  data.scaled.xy <- data.xy
  data.scaled.xy[ ,"y"] <- (start0.y - data.xy[ ,"y"])/resolution
  data.scaled.xy[ ,"x"] <- (data.xy[ ,"x"] - start0.x)/resolution 
  out <- list(grid.scaled.xy = grid.scaled.xy,
              grid.xy = grid.xy,
              data.scaled.xy = data.scaled.xy,
              data.xy = data.xy)
  return(out)
}


###############################################################
##################### MakeAugmentation ########################
###############################################################
#' @title MakeAugmentation function
#'
#' @description \code{MakeAugmentation} increases the dimensions of an object along
#'  the individual and/or the year dimension. It returns a \code{Vector} or \code{Matrix} 
#'  object with the expanded object.
#'
#' @param y A \code{Vector} or \code{Matrix} object containing the individual detection histories.
#' @param aug.factor A \code{numeric} object defining the augmentation factor to be used.
#' @param replace.value A \code{numeric} object defining the value to be repeated for augmented individuals.
#' 
#' @return A \code{Vector} or \code{Matrix} object containing the augmented y.

MakeAugmentation <- function( y,
                              aug.factor= NULL,
                              aug.years = NULL,
                              replace.value = NA){
  ## Vector Data augmentation
  if(is.vector(y)){
    if(is.null(names(y))){ 
      names(y) <- 1:length(y) 
    }
    if(!is.null(aug.factor)){
      y.aug <- c(y, rep(replace.value, round(length(y) * aug.factor)))
      names(y.aug) <- c(names(y), rep("Augmented", round(length(y) * aug.factor)))
      y <- y.aug
    }
  }
  
  ## Matrix augmentation
  if(is.matrix(y)){
    if(is.null(dimnames(y))){
      dimnames(y) <- list(1:dim(y)[1], 1:dim(y)[2])
    }
    if(!is.null(aug.factor)){
      n.tot <- round(dim(y)[1]*(1 + aug.factor))
      y.aug <- matrix(replace.value, n.tot, dim(y)[2])
      y.aug[1:dim(y)[1], ] <- y
      dimnames(y.aug) <- list(c( dimnames(y)[[1]], rep("Augmented", n.tot - dim(y)[1])),
                              dimnames(y)[[2]])
      y <- y.aug
    }
    
  }  
  
  return (y)
}


###############################################################
##################### MakeInitXY ##############################
###############################################################

#' @title Function to set initial values of XY coordinates of ACS
#'
#' @description
#' \code{InitXY} returns a matrix object with with the x coordinates ([,1]) and y coordinates  ([,2]). it returns the location of a detection for an detected individual and a random location for augmented individuals. 
#' 
#' @param y \code{matrix}  with individual detections. row=Ids, col= detectors. 
#' @param habitat.mx \code{matrix} of corrdinates from the buffered habitat grid
#' @param detector.xy \code{matrix} with coordinates of detectors scaled
#' @param IDCells.mx \code{matrix} the ID of the cells of the habitat
#' @param grid.xy \code{matrix} with the scaled coordinates of the habitat 
#' @examples
#' # Generate simulated detection histories within an SCR framework:
#' MakeInitXY(   y = data$y
#'              , habitat.mx = data$habitat.mx
#'              , detector.xy = data$detector.xy
#'              , grid.xy = grid.xy)

MakeInitXY <- function(    y = y
                           , habitat.mx = habitat.mx
                           , detector.xy = detector.xy
                           , IDCells.mx = IDCells.mx
                           , grid.xy = grid.xy){
  # , xy.bounds = NULL ){
  
  # ---- STEP 1: GET THE OBJECT READY TO STORE THE DATA ----- 
  # Make it general to work with the time dimention 
  n.years <- 1
  # n.years <- ifelse(length(dim(y))>=3, dim(y)[3], 1)
  n.individuals <- dim(y)[1]
  n.detectors <- dim(detector.xy)[1]
  
  # if(length(dim(y)) == 2){y <- array(y, c(n.individuals, n.detectors, n.years))}
  # if(length(dim(detector.xy)) == 2){detector.xy <- array(detector.xy, c(n.detectors, 2, n.years))}
  # n.year <- dim(y)[3]
  
  
  # STORE THE DIMENSION OF THE WINDOW 
  # IF NO FRAGMENTATION IS USED 
  dim.mx <- array(0, c(n.individuals, 2, 2))#, n.year))
  # for(t in 1:n.year){
  dim.mx[,1,1] <- 1 # min x
  dim.mx[,2,1] <- 1 # min y
  dim.mx[,1,2] <- dim(habitat.mx)[2]# max x
  dim.mx[,2,2] <- dim(habitat.mx)[1]# max y 
  
  # if(!is.null(xy.bounds)){
  #     dim.mx <- xy.bounds
  #     dim.mx[,1,1] <-  floor(dim.mx[,1,1]) + 1 # Add one to match the habitat 
  #     dim.mx[,2,1] <-  floor(dim.mx[,2,1]) + 1 # Add one to match the habitat
  # }
  # }
  
  # IF  FRAGMENTATION IS USED 
  
  
  
  # empty <- array(numeric(),c(n.individuals, 2, n.year))
  empty <- array(numeric(),c(n.individuals, 2))
  ids <- lapply(c(1:n.individuals), function(x) x)
  
  # ---- STEP 2: FIND SXY FOR ALL INDIVIDUALS  ----- 
  
  # for(t in 1:n.year){
  listt <- apply(y, 1, function(x) which(x>0, arr.ind = T))# obtain detectors with detections for each ids 
  
  if(length(listt)==0){listt <- list()
  for(i in 1:n.individuals){
    listt[[i]] <- integer()
  }
  }
  
  # empty[,,t] <- do.call(rbind, lapply(ids, function(i) {
  empty <- do.call(rbind, lapply(ids, function(i) {
    
    #print(i)
    x <-  listt[[i]]
    # if only one detections 
    if(length(x)==1){ # If only one detection, use that detection as a starting value
      detector.xy[x,]
    }else{
      
      # if several detections    
      if(length(x)>1){# if more than 1 detection, use average value
        mean.det <- colMeans(detector.xy[x,])
        # if doesnt end up in habitat, use a random coordinate within habitat 
        if(habitat.mx[floor(mean.det)[2]+1,floor(mean.det)[1]+1]==1){
          mean.det 
        }else{
          min.x <- dim.mx[i,1,1]
          max.x <- dim.mx[i,1,2]
          min.y <- dim.mx[i,2,1]
          max.y <- dim.mx[i,2,2]
          
          window.mx.ind <- habitat.mx[c(min.y:max.y), c(min.x:max.x)]
          window.ID.ind <- IDCells.mx[c(min.y:max.y), c(min.x:max.x)]
          
          id.cell <- window.ID.ind[window.mx.ind==1]
          sxy.grid <- matrix(grid.xy[id.cell,], nrow=length(id.cell),ncol=2)
          
          matrix(sxy.grid[floor(runif(1, 1, length(id.cell))),], ncol=2, nrow=1)                     }
      }else{
        
        # if no detections (augmented IDS)
        if(length(x)==0){
          
          min.x <- dim.mx[i,1,1]
          max.x <- dim.mx[i,1,2]
          min.y <- dim.mx[i,2,1]
          max.y <- dim.mx[i,2,2]
          
          window.mx.ind <- habitat.mx[c(min.y:max.y), c(min.x:max.x)]
          window.ID.ind <- IDCells.mx[c(min.y:max.y), c(min.x:max.x)]
          
          id.cell <- window.ID.ind[window.mx.ind==1]
          sxy.grid <- matrix(grid.xy[id.cell,], nrow=length(id.cell),ncol=2)
          
          matrix(sxy.grid[floor(runif(1, 1, length(id.cell))),], ncol=2, nrow=1)
        }
      }}      
  })) 
  
  # }
  return(empty)
  
}


###############################################################
##################### Modelcodes ##############################
###############################################################

#########################################################################
##################### FIT modelCode: Half-normal ########################
#########################################################################
modelCode_hn_LESSCachedAllSparse <- nimbleCode({
  
  #----- SPATIAL PROCESS 
  sumHabInt <- sum(mu[1:numHabWindows])
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernPP(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      intensityWeights = mu[1:numHabWindows],
      sumIntensity = sumHabInt,
      numPoints = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max]
      
    )
  }#i
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  
  ##---- DETECTION PROCESS
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 40) # 50)
  # alpha <- -1 / (2 * sigma * sigma)
  
  for (i in 1:n.individuals){
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    ## USE z[i] TO TURN OFF CALCULATIONS WHEN NOT NECESSARY 
    y[i, 1:nMaxDetectors] ~ dbin_LESSCachedAllSparse_detfun_hn(pZero = p0,
                                                                  sxy = sxy[i,1:2],
                                                                  sigma = sigma,
                                                                  nbDetections = nbDetections[i],
                                                                  yDets = yDets[i,1:nMaxDetectors],
                                                                  detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                  trials = trials[1:n.detectors],
                                                                  detectorIndex = detectorIndex[1:n.cells,1:maxNBDets], 
                                                                  nDetectorsLESS = nDetectorsLESS[1:n.cells],
                                                                  ResizeFactor = ResizeFactor,
                                                                  maxNBDets = maxNBDets,
                                                                  habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                  maxDist=MaxDist,
                                                                  indicator = z[i]) #,
  
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})
modelCode_hn_LESS_Cached <- nimbleCode({
  
  #----- SPATIAL PROCESS 
  sumHabInt <- sum(mu[1:numHabWindows])
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernPP(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      intensityWeights = mu[1:numHabWindows],
      sumIntensity = sumHabInt,
      numPoints = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max]
      
    )
    
  }#i
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  ##---- DETECTION PROCESS
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 40) # 50)
  
  alpha <- -1 / (2 * sigma * sigma)
  
  for (i in 1:n.individuals){
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    ## USE z[i] TO TURN OFF CALCULATIONS WHEN NOT NECESSARY 
    
    # # Calculate distance between i-th AC location and each detector
    d2[i,1:n.detectors] <- pow(sxy[i,1]-detector.xy[1:n.detectors,1],2) +
      pow(sxy[i,2]-detector.xy[1:n.detectors,2],2)
    # # Calculate detection prob. of i-th individual at each detector
    p[i,1:n.detectors] <- p0*exp(alpha * d2[i,1:n.detectors])
    # # Calculate prob. density of CR data for i-th individual at each detector
    
    
    y[i, 1:n.detectors] ~ dbin_LESS_Cached_detfun_hn(pZero = p0,
                                                                 sxy = sxy[i,1:2],
                                                                 sigma = sigma,
                                                                 detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                 trials = trials[1:n.detectors],
                                                                 detectorIndex = detectorIndex[1:n.cells,1:maxNBDets],
                                                                 nDetectorsLESS = nDetectorsLESS[1:n.cells],
                                                                 ResizeFactor = ResizeFactor,
                                                                 habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                 indicator = z[i]) #,
  }#i
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
 
})

#########################################################################
##################### FIT modelCode: Half-normal plateau ################
#########################################################################
modelCode_hnplateau_LESSCachedAllSparse <- nimbleCode({
  
  #----- SPATIAL PROCESS 
  sumHabInt <- sum(mu[1:numHabWindows])
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernPP(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      intensityWeights = mu[1:numHabWindows],
      sumIntensity = sumHabInt,
      numPoints = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max]
      
    )
   
  }#i
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  
  ##---- DETECTION PROCESS
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 40) # 50)
  wd ~ dunif(0, 20)
  # alpha <- -1 / (2 * sigma * sigma)
  
  for (i in 1:n.individuals){
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    ## USE z[i] TO TURN OFF CALCULATIONS WHEN NOT NECESSARY 
    y[i, 1:nMaxDetectors] ~ dbin_LESSCachedAllSparse_detfun_hnplateau(pZero = p0,
                                                                         sxy = sxy[i,1:2],
                                                                         sigma = sigma,
                                                                         nbDetections = nbDetections[i],
                                                                         yDets = yDets[i,1:nMaxDetectors],
                                                                         detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                         trials = trials[1:n.detectors],
                                                                         detectorIndex = detectorIndex[1:n.cells,1:maxNBDets], 
                                                                         nDetectorsLESS = nDetectorsLESS[1:n.cells],
                                                                         ResizeFactor = ResizeFactor,
                                                                         maxNBDets = maxNBDets,
                                                                         habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                         maxDist=MaxDist,
                                                                         indicator = z[i],
                                                                         wd = wd) 
    
    
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])

  
})
modelCode_hnplateau_LESS_Cached <- nimbleCode({
  
  #----- SPATIAL PROCESS 
  sumHabInt <- sum(mu[1:numHabWindows])
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernPP(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      intensityWeights = mu[1:numHabWindows],
      sumIntensity = sumHabInt,
      numPoints = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max]
      
    )
  
  }#i
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  ##---- DETECTION PROCESS
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  # wd <- 2
  wd ~ dunif(0, 20)
  
  alpha <- -1 / (2 * sigma * sigma)
  
  for (i in 1:n.individuals){
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    ## USE z[i] TO TURN OFF CALCULATIONS WHEN NOT NECESSARY 
    
     # # Calculate distance between i-th AC location and each detector
    D[i,1:n.detectors] <- pow(pow(sxy[i,1]-detector.xy[1:n.detectors,1],2) +
                                pow(sxy[i,2]-detector.xy[1:n.detectors,2],2), 0.5)
    # 
    # #   # Calculate detection prob. of i-th individual at each detector
    # # c1[i,1:n.detectors] <- D[i,1:n.detectors] <= wd
    for (j in 1:n.detectors){
      c1[i,j] <- D[i,j] <= wd
      # dens[i,j] <- c1[i,j] + (1-c1[i,j])*exp(alpha *(D[i,j]-wd)*(D[i,j]-wd))
    }
    dens[i,1:n.detectors] <- c1[i,1:n.detectors] +
      (1-c1[i,1:n.detectors])*exp(alpha *(D[i,1:n.detectors]-wd)*(D[i,1:n.detectors]-wd))
    # 
    p[i,1:n.detectors] <- p0*dens[i,1:n.detectors]
    # # Calculate prob. density of CR data for i-th individual at each detector
     
    y[i, 1:n.detectors] ~ dbin_LESS_Cached_detfun_hnplateau(pZero = p0,
                                                                        sxy = sxy[i,1:2],
                                                                        sigma = sigma,
                                                                        detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                        trials = trials[1:n.detectors],
                                                                        detectorIndex = detectorIndex[1:n.cells,1:maxNBDets],
                                                                        nDetectorsLESS = nDetectorsLESS[1:n.cells],
                                                                        ResizeFactor = ResizeFactor,
                                                                        habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                        indicator = z[i],
                                                                        wd = wd)
  }#i
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})


#########################################################################
##################### FIT modelCode: Exponential ########################
#########################################################################
modelCode_ex_LESSCachedAllSparse <- nimbleCode({
  
  #----- SPATIAL PROCESS 
  sumHabInt <- sum(mu[1:numHabWindows])
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernPP(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      intensityWeights = mu[1:numHabWindows],
      sumIntensity = sumHabInt,
      numPoints = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max]
      
    )
    
  }#i
  
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  
  ##---- DETECTION PROCESS
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  # wd ~ dunif(0, 20)
  # alpha <- -1 / (2 * sigma * sigma)
  
  for (i in 1:n.individuals){
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    ## USE z[i] TO TURN OFF CALCULATIONS WHEN NOT NECESSARY 
    y[i, 1:nMaxDetectors] ~ dbin_LESSCachedAllSparse_detfun_ex(pZero = p0,
                                                                  sxy = sxy[i,1:2],
                                                                  sigma = sigma,
                                                                  nbDetections = nbDetections[i],
                                                                  yDets = yDets[i,1:nMaxDetectors],
                                                                  detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                  trials = trials[1:n.detectors],
                                                                  detectorIndex = detectorIndex[1:n.cells,1:maxNBDets], 
                                                                  nDetectorsLESS = nDetectorsLESS[1:n.cells],
                                                                  ResizeFactor = ResizeFactor,
                                                                  maxNBDets = maxNBDets,
                                                                  habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                  maxDist=MaxDist,
                                                                  indicator = z[i]) #,
     
  }#i
  
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})
modelCode_ex_LESS_Cached <- nimbleCode({
  
  #----- SPATIAL PROCESS 
  sumHabInt <- sum(mu[1:numHabWindows])
  for(i in 1:n.individuals){
    sxy[i, 1:2] ~ dbernPP(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      intensityWeights = mu[1:numHabWindows],
      sumIntensity = sumHabInt,
      numPoints = 1,
      habitatGrid = habitatGrid[1:y.max,1:x.max]
      
    )
 
  }#i
  ##---- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  for (i in 1:n.individuals){
    z[i] ~ dbern(psi)
  }#i
  ##---- DETECTION PROCESS
  p0 ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  
  alpha <- -1 / sigma# (2 * sigma * sigma)
  
  for (i in 1:n.individuals){
    ## LESS APPROACH AND USE OF SPARSE MATRICES
    ## USE z[i] TO TURN OFF CALCULATIONS WHEN NOT NECESSARY 
    
     # # Calculate distance between i-th AC location and each detector
    D[i,1:n.detectors] <- pow(pow(sxy[i,1]-detector.xy[1:n.detectors,1],2) +
                                pow(sxy[i,2]-detector.xy[1:n.detectors,2],2), 0.5)
    
    # # Calculate detection prob. of i-th individual at each detector
    p[i,1:n.detectors] <- p0*exp(alpha * D[i,1:n.detectors])
    # # Calculate prob. density of CR data for i-th individual at each detector
    
    y[i, 1:n.detectors] ~ dbin_LESS_Cached_detfun_ex(pZero = p0,
                                                                 sxy = sxy[i,1:2],
                                                                 sigma = sigma,
                                                                 detector.xy =  detector.xy[1:n.detectors,1:2],
                                                                 trials = trials[1:n.detectors],
                                                                 detectorIndex = detectorIndex[1:n.cells,1:maxNBDets],
                                                                 nDetectorsLESS = nDetectorsLESS[1:n.cells],
                                                                 ResizeFactor = ResizeFactor,
                                                                 habitatID = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                                 indicator = z[i]) #,
  }#i
  ##---- DERIVED QUANTITIES
  N <- sum(z[1:n.individuals])
  
})


###############################################################
##################### dbernPP - Bernoulli Point Process #######
###############################################################

# Density function
#' @title Probability Density of a Bernoulli Point Process 
#' 
#' @description Probability density function for a Bernoulli point process (a binomial point process 
#' with one point only) observed within a set of observation windows.
#' 
#' @param x A vector of coordinates of one point generated from the point process.
#' @param lowerCoords A matrix containing the lower coordinates of the observation windows. 
#' Each row represents an observation window.
#' @param upperCoords A matrix containing the upper coordinates of the observation windows. 
#' Each row represents an observation window.
#' @param intensityWeights A vector of length equal to the number of observation windows.
#' Sets the baseline intensity for each observation window.
#' @param sumIntensity A value calculated by integrating the intensity function over the whole spatial domain. 
#' A discrete summation when the domain is divided into grid cells. This value is obtained by summing up
#' the product of \code{intensityWeights} and \code{windowSizes}; however providing this as an argument 
#' helps to reduce some calculations when the function is used repeatedly in a loop in some cases.
#' @param numPoints The number of points generated by the point process. 
#' If this value is non-zero then one point will be generated. 
#' This is useful for NIMBLE models where \code{dbinomPP} is nested in a loop and 
#' where there may be zero points in some cases.
#' @param habitatGrid  A \code{matrix} object with ID of each habitat grid cells
#' @param numWindows The number of observation windows. If this value is less than zero
#' then the number of windows is taken from \code{dim(lowerCoords)[1]}. The value
#' is then used to truncate \code{lowerCoords} and \code{upperCoords} (so that extra rows
#' beyond \code{numWindows} are ignored). This is useful when using this function in
#' NIMBLE models where \code{dbernPP} is nested in a loop and where may be zero or one
#' observation window in some cases.
#' @param parsCheck If \code{TRUE} then check the validity of the input parameter(s).
#' @param log If \code{TRUE} then return the log density.
#' 
#' @return A scalar containing the (log) probability density.
#' 
#' @author 
#' @export
#' 
dbernPP <- nimbleFunction(
  run = function(
    x                = double(1),
    lowerCoords      = double(2),
    upperCoords      = double(2),
    intensityWeights = double(1),
    sumIntensity     = double(0),
    numPoints        = integer(0),
    habitatGrid      = double(2),
    log              = integer(0, default = 0)
  ){
    ## Likelihood is one if there are zero points
    if(numPoints == 0) {
      if(log){
        return(0.0) 
      } 
      else {
        return(1.0)
      } 
    }
    nrowMax <- dim(habitatGrid)[1]
    ncolMax <- dim(habitatGrid)[2]
    # MAKE SURE THE POINT FALLS WITHIN THE HABITAT 
    if(x[2] < 0 |  x[2] > ncolMax | x[1] < 0 | x[1] > nrowMax){
      if(log) {
        return(-Inf)
      } else {
        return(0.0)
      }
      
    }else{
      ## Find which window the point x falls within
      windowInd <- habitatGrid[trunc(x[2])+1, trunc(x[1])+1]
      pointIntensity <- intensityWeights[windowInd]
    }
    
    ## Log probability density 
    outProb <- log(pointIntensity) - log(sumIntensity)
    if(log == 0) {
      outProb <- exp(outProb)
    }
    return(outProb)
    returnType(double(0))
  }
)
###############################################################
# simulation function 
#' @title Generate a Sample from a Bernoulli Point Process
#' 
#' @description Function to generate a sample (i.e. one point) from a Bernoulli
#' point process for a series of observation windows.
#' 
#' @param n The number of samples to generate (this should always be set to 1)
#' @param lowerCoords A matrix containing the lower coordinates of the observation
#' windows. Each row represents an observation window.
#' @param upperCoords A matrix containing the upper coordinates of the observation
#' windows. Each row represents an observation window.
#' @param intensityWeights A vector of length equal to the number of observation windows.
#' Sets the baseline intensity for each observation window.
#' @param windowSizes A vector of all the observation window sizes.
#' @param sumIntensity A value calculated by integrating the intensity function over the whole spatial domain. 
#' A discrete summation when the domain is divided into grid cells. This value is obtained by summing up
#' the product of \code{intensityWeights} and \code{windowSizes}; however providing this as an argument 
#' helps to reduce some calculations when the function is used repeatedly in a loop in some cases. 
#' @param numPoints The number of points generated by the point process. 
#' If this value is non-zero then one point will be generated. 
#' This is useful for NIMBLE models where \code{dbinomPP} is nested in a loop and 
#' where there may be zero points in some cases.
#' @param numWindows The number of observation windows. If this value is less than zero
#' then the number of windows is taken from \code{dim(lowerCoords)[1]}. The value
#' is then used to truncate \code{lowerCoords} and \code{upperCoords} (so that extra rows
#' beyond \code{numWindows} are ignored). This is useful when using this function in
#' NIMBLE models where \code{dbernPP} is nested in a loop and where may be zero or one
#' observation window in some cases.
#' @param parsCheck If \code{TRUE} then check the validity of the input parameter(s).
#' 
#' @return A vector of coordinates for a point generated from the Bernoulli point process
#' 
#' @author Wei Zhang
#' @export
#' 
rbernPP <- nimbleFunction(
  run = function(
    n                = integer(0),
    lowerCoords      = double(2),
    upperCoords      = double(2),
    intensityWeights = double(1),
    sumIntensity     = double(0),
    numPoints        = integer(0),
    habitatGrid      = double(2)
    
  ){
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("rbernPP only allows n = 1; using n = 1")
    }  
    ## Ensure that there is at least one window with non-zero size         
    # if(sum(windowSizes) == 0){ 
    #   stop("at least one of the observation windows should have a non-zero size")
    # }
    ## If no points to be generated, return a null vector
    if(numPoints == 0) {
      outCoordinates <- numeric(length = 0)
      return(outCoordinates)
    } 
    ## Intensities of all windows
    areaIntensities <- intensityWeights #* windowSizes
    ## Observation window index
    obsWindowInd <- rcat(1, areaIntensities)
    ## Generate coordinates within the relevant observation window
    targetLowerCoords <- lowerCoords[obsWindowInd,]
    targetUpperCoords <- upperCoords[obsWindowInd,]
    numDims <- dim(lowerCoords)[2]
    outCoordinates <- targetLowerCoords[1:numDims] + runif(numDims, 0.0, 1.0) * (targetUpperCoords[1:numDims] - targetLowerCoords[1:numDims])
    
    return(outCoordinates)
    returnType(double(1)) 
  }
)


registerDistributions( list(
###############################################################
  
  ## Define the dbernPP distribution 
  dbernPP = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dbernPP(lowerCoords, upperCoords, intensityWeights, sumIntensity, numPoints, habitatGrid)",
    # Set the input and output types and dimension structure
    types = c("value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)", "intensityWeights = double(1)", 
              "sumIntensity = double(0)", "numPoints = double(0)","habitatGrid = double(2)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
)
)

###############################################################
##################### dbin_LESSCachedAllSparse_detfun_hn ######
###############################################################

#' @title Function to create a NIMBLE custom distribution for faster SCR model runs with respect to half-normal detection function.
#'
#' @description
#' \code{dbin_LESSCachedAllSparse_detfun_hn} returns the likelihood of a given individual & spatial binary detection history y with respect to half-normal detection function
#' 
#' @param x \code{Vector} containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param nbDetections A \code{integer} giving the number of detector locations where the individual is detected.
#' @param yDets A \code{Vector} giving the detector IDs where detection took place
#' @param detector.xy A \code{Matrix} giving the scaled detector locations
#' @param trials A \code{Vector} giving the number of detection sessions for each detector
#' @param detectorIndex A \code{Matrix}  with the detectors ID (along columns) assigned to each cells (along rows) of the habitat under the LESS approach.
#' @param nDetectorsLESS A \code{Vector}  of dimensions n.cells  with the number of detectors assigned to each cells of the habitat under the LESS approach.
#' @param ResizeFactor A \code{Numeric} giving scaling factor for the habitat grid cells
#' @param maxNBDets A \code{Numeric} giving the maximum number of possible detectors within a radius of any habitat grid cell
#' @param habitatID A \code{Matrix} giving the IDs of each habitat grid cell 
#' @param maxDist A \code{Numeric} with the maximum distance to the AC where detections are considered possible. This applies a local evaluation of the state space (Milleret et al. 2018 Ecology and Evolution)
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)

#### Density function
dbin_LESSCachedAllSparse_detfun_hn <- nimbleFunction(run = function( x = double(1)
                                                                        , pZero = double(0)
                                                                        , sxy = double(1)
                                                                        , sigma = double(0, default = 1.0)
                                                                        , nbDetections = double(0)
                                                                        , yDets = double(1)
                                                                        , detector.xy = double(2)
                                                                        , trials = double(1)
                                                                        , detectorIndex = double(2)
                                                                        , nDetectorsLESS = double(1)
                                                                        , ResizeFactor = double(0, default = 1)
                                                                        , maxNBDets = double(0)
                                                                        , habitatID = double(2)                      
                                                                        , maxDist = double(0, default = 0.0)
                                                                        , indicator = double(0, default = 1.0)
                                                                        , log = integer(0, default = 0)
){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ## CHECK INPUT DIMENSIONS
  nDetectors <- length(trials)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(nbDetections == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  ## GET DETECTOR INDEX FROM THE HABITATID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  index <- detectorIndex[sxyID,1:nDetectorsLESS[sxyID]]
  ## GET NECESSARY INFO
  n.detectors <- length(index)
  maxDist_squared <- maxDist*maxDist
  
  ## RECREATE Y
  y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
  if(nbDetections > 0){
    for(j in 1:nbDetections){
      y[yDets[j]] <- x[j]
      ## check if a detection is out of the "detection window"
      if(sum(yDets[j]==index)==0){
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
  }
  
  
  ## CALCULATE THE LIKELIHOOD
  alpha <- -1.0 / (2.0 * sigma * sigma)
  # if(detfun != 3.0) alpha <- -1.0 / (2.0 * sigma * sigma)
  
  
  logProb <- 0.0 
  count <- 1 
  index1 <- c(index, 0) # so the count is not an issue
  
  for(j in 1:nDetectors){
    if(index1[count] == j){ # IF J IS EQUAL TO THE RELEVANT DETECTOR 
      D <- pow(pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2), 0.5)
      dens <- exp(alpha * D*D)
      
      p <- pZero * dens 
      logProb <- logProb + dbinom(y[j], prob = p, size = trials[j], log = TRUE)
      count <- count + 1
    }#else{logProb <- logProb + dbinom(y[j], prob = 0, size = trials[j], log = TRUE)}
  }
  
  if(log)return(logProb)
  else return(exp(logProb))
  
})

#### Registration 
registerDistributions(list(
  dbin_LESSCachedAllSparse_detfun_hn = list(
    BUGSdist = "dbin_LESSCachedAllSparse_detfun_hn(pZero, sxy, sigma, nbDetections, yDets, 
                detector.xy, trials, detectorIndex, nDetectorsLESS, ResizeFactor, maxNBDets, habitatID, maxDist, indicator)",
    types = c( "value = double(1)", "pZero = double(0)","sxy = double(1)", "sigma = double(0)",
               "nbDetections = double(0)", "yDets = double(1)", "detector.xy = double(2)","trials = double(1)","detectorIndex = double(2)" ,
               "nDetectorsLESS = double(1)", "ResizeFactor = double(0)",
               "maxNBDets = double(0)","habitatID= double(2)",
               "maxDist = double(0)" ,"indicator = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE)))



###############################################################
##################### dbin_LESSCachedAllSparse_detfun_hnp #####
###############################################################

#' @title Function to create a NIMBLE custom distribution for faster SCR model runs with respect to half-normal plateau detection function.
#'
#' @description
#' \code{dbin_LESSCachedAllSparse_detfun_hnplateau} returns the likelihood of a given individual & spatial binary detection history y with respect to half-normal plateau detection function
#' 
#' @param x \code{Vector} containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal plateau detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal plateau detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param nbDetections A \code{integer} giving the number of detector locations where the individual is detected.
#' @param yDets A \code{Vector} giving the detector IDs where detection took place
#' @param detector.xy A \code{Matrix} giving the scaled detector locations
#' @param trials A \code{Vector} giving the number of detection sessions for each detector
#' @param detectorIndex A \code{Matrix}  with the detectors ID (along columns) assigned to each cells (along rows) of the habitat under the LESS approach.
#' @param nDetectorsLESS A \code{Vector}  of dimensions n.cells  with the number of detectors assigned to each cells of the habitat under the LESS approach.
#' @param ResizeFactor A \code{Numeric} giving scaling factor for the habitat grid cells
#' @param maxNBDets A \code{Numeric} giving the maximum number of possible detectors within a radius of any habitat grid cell
#' @param habitatID A \code{Matrix} giving the IDs of each habitat grid cell 
#' @param maxDist A \code{Numeric} with the maximum distance to the AC where detections are considered possible. This applies a local evaluation of the state space (Milleret et al. 2018 Ecology and Evolution)
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param wd A \code{Numeric} giving the length of the plateaufor the half-normal plateau detection function
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)

#### Density function
dbin_LESSCachedAllSparse_detfun_hnplateau <- nimbleFunction(run = function( x = double(1)
                                                                               , pZero = double(0)
                                                                               , sxy = double(1)
                                                                               , sigma = double(0, default = 1.0)
                                                                               , nbDetections = double(0)
                                                                               , yDets = double(1)
                                                                               , detector.xy = double(2)
                                                                               , trials = double(1)
                                                                               , detectorIndex = double(2)
                                                                               , nDetectorsLESS = double(1)
                                                                               , ResizeFactor = double(0, default = 1)
                                                                               , maxNBDets = double(0)
                                                                               , habitatID = double(2)                      
                                                                               , maxDist = double(0, default = 0.0)
                                                                               , indicator = double(0, default = 1.0)
                                                                               , wd = double(0, default = 3.0)
                                                                               , log = integer(0, default = 0)
){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ## CHECK INPUT DIMENSIONS
  nDetectors <- length(trials)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(nbDetections == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  ## GET DETECTOR INDEX FROM THE HABITATID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  index <- detectorIndex[sxyID,1:nDetectorsLESS[sxyID]]
  ## GET NECESSARY INFO
  n.detectors <- length(index)
  maxDist_squared <- maxDist*maxDist
  
  ## RECREATE Y
  y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
  if(nbDetections > 0){
    for(j in 1:nbDetections){
      y[yDets[j]] <- x[j]
      ## check if a detection is out of the "detection window"
      if(sum(yDets[j]==index)==0){
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
  }
  
  
  ## CALCULATE THE LIKELIHOOD
  alpha <- -1.0 / (2.0 * sigma * sigma)
  # if(detfun != 3.0) alpha <- -1.0 / (2.0 * sigma * sigma)
  
  
  logProb <- 0.0 
  count <- 1 
  index1 <- c(index, 0) # so the count is not an issue
  
  for(j in 1:nDetectors){
    if(index1[count] == j){ # IF J IS EQUAL TO THE RELEVANT DETECTOR 
      D <- pow(pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2), 0.5)
      if(D <= wd) dens <- 1.0
      if(D >  wd) dens <- exp(alpha *(D-wd)*(D-wd))#/(2.0*sigma*sigma))
      p <- pZero * dens 
      logProb <- logProb + dbinom(y[j], prob = p, size = trials[j], log = TRUE)
      count <- count + 1
    }#else{logProb <- logProb + dbinom(y[j], prob = 0, size = trials[j], log = TRUE)}
  }
  
  if(log)return(logProb)
  else return(exp(logProb))
  
})

#### Registration 
registerDistributions(list(
  dbin_LESSCachedAllSparse_detfun_hnplateau = list(
    BUGSdist = "dbin_LESSCachedAllSparse_detfun_hnplateau(pZero, sxy, sigma, nbDetections, yDets, 
                detector.xy, trials, detectorIndex, nDetectorsLESS, ResizeFactor, maxNBDets, habitatID, maxDist, indicator,
                wd)",
    types = c( "value = double(1)", "pZero = double(0)","sxy = double(1)", "sigma = double(0)",
               "nbDetections = double(0)", "yDets = double(1)", "detector.xy = double(2)","trials = double(1)","detectorIndex = double(2)" ,
               "nDetectorsLESS = double(1)", "ResizeFactor = double(0)",
               "maxNBDets = double(0)","habitatID= double(2)",
               "maxDist = double(0)" ,"indicator = double(0)",
               "wd = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE)))



###############################################################
##################### dbin_LESSCachedAllSparse_detfun_ex ######
###############################################################

#' @title Function to create a NIMBLE custom distribution for faster SCR model runs with respect to exponential detection function.
#'
#' @description
#' \code{dbin_LESSCachedAllSparse_detfun_ex} returns the likelihood of a given individual & spatial binary detection history y with respect to exponential detection function
#' 
#' @param x \code{Vector} containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the exponential detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the exponential detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param nbDetections A \code{integer} giving the number of detector locations where the individual is detected.
#' @param yDets A \code{Vector} giving the detector IDs where detection took place
#' @param detector.xy A \code{Matrix} giving the scaled detector locations
#' @param trials A \code{Vector} giving the number of detection sessions for each detector
#' @param detectorIndex A \code{Matrix}  with the detectors ID (along columns) assigned to each cells (along rows) of the habitat under the LESS approach.
#' @param nDetectorsLESS A \code{Vector}  of dimensions n.cells  with the number of detectors assigned to each cells of the habitat under the LESS approach.
#' @param ResizeFactor A \code{Numeric} giving scaling factor for the habitat grid cells
#' @param maxNBDets A \code{Numeric} giving the maximum number of possible detectors within a radius of any habitat grid cell
#' @param habitatID A \code{Matrix} giving the IDs of each habitat grid cell 
#' @param maxDist A \code{Numeric} with the maximum distance to the AC where detections are considered possible. This applies a local evaluation of the state space (Milleret et al. 2018 Ecology and Evolution)
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)

#### Density function
dbin_LESSCachedAllSparse_detfun_ex <- nimbleFunction(run = function( x = double(1)
                                                                     , pZero = double(0)
                                                                     , sxy = double(1)
                                                                     , sigma = double(0, default = 1.0)
                                                                     , nbDetections = double(0)
                                                                     , yDets = double(1)
                                                                     , detector.xy = double(2)
                                                                     , trials = double(1)
                                                                     , detectorIndex = double(2)
                                                                     , nDetectorsLESS = double(1)
                                                                     , ResizeFactor = double(0, default = 1)
                                                                     , maxNBDets = double(0)
                                                                     , habitatID = double(2)                      
                                                                     , maxDist = double(0, default = 0.0)
                                                                     , indicator = double(0, default = 1.0)
                                                                     , log = integer(0, default = 0)
){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ## CHECK INPUT DIMENSIONS
  nDetectors <- length(trials)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(nbDetections == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  ## GET DETECTOR INDEX FROM THE HABITATID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  index <- detectorIndex[sxyID,1:nDetectorsLESS[sxyID]]
  ## GET NECESSARY INFO
  n.detectors <- length(index)
  maxDist_squared <- maxDist*maxDist
  
  ## RECREATE Y
  y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
  if(nbDetections > 0){
    for(j in 1:nbDetections){
      y[yDets[j]] <- x[j]
      ## check if a detection is out of the "detection window"
      if(sum(yDets[j]==index)==0){
        if(log == 0) return(0.0)
        else return(-Inf)
      }
    }
  }
  
  
  ## CALCULATE THE LIKELIHOOD
  alpha <- -1.0 / sigma # (2.0 * sigma * sigma)
  # if(detfun != 3.0) alpha <- -1.0 / (2.0 * sigma * sigma)
  
  
  logProb <- 0.0 
  count <- 1 
  index1 <- c(index, 0) # so the count is not an issue
  
  for(j in 1:nDetectors){
    if(index1[count] == j){ # IF J IS EQUAL TO THE RELEVANT DETECTOR 
      D <- pow(pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2), 0.5)
      dens <- exp(alpha * D)
      p <- pZero * dens 
      # # p <- pZero * exp(alpha * d2)
      logProb <- logProb + dbinom(y[j], prob = p, size = trials[j], log = TRUE)
      count <- count + 1
    }#else{logProb <- logProb + dbinom(y[j], prob = 0, size = trials[j], log = TRUE)}
  }
  
  if(log)return(logProb)
  else return(exp(logProb))
  
})

#### Registration 
registerDistributions(list(
  dbin_LESSCachedAllSparse_detfun_ex = list(
    BUGSdist = "dbin_LESSCachedAllSparse_detfun_ex(pZero, sxy, sigma, nbDetections, yDets, 
                detector.xy, trials, detectorIndex, nDetectorsLESS, ResizeFactor, maxNBDets, habitatID, maxDist, indicator)",
    types = c( "value = double(1)", "pZero = double(0)","sxy = double(1)", "sigma = double(0)",
               "nbDetections = double(0)", "yDets = double(1)", "detector.xy = double(2)","trials = double(1)","detectorIndex = double(2)" ,
               "nDetectorsLESS = double(1)", "ResizeFactor = double(0)",
               "maxNBDets = double(0)","habitatID= double(2)",
               "maxDist = double(0)" ,"indicator = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE)))

###############################################################
##################### dbin_LESS_Cached_detfun_hn ##############
###############################################################
#' @title Function to create a NIMBLE custom distribution for faster SCR model runs with respect to half-normal detection function.
#'
#' @description
#' \code{dbin_LESS_Cached_detfun_hn} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] with respect to half-normal detection function
#' 
#' @param x \code{Vector} of length n.detectors containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param detector.xy A \code{Matrix}  of dimensions n.detectors x 2 with the locations of each detector.
#' @param trials A \code{Vector} giving the number of detection sessions for each detector.
#' @param detectorIndex A \code{Matrix}  with the detectors ID (along columns) assigned to each cells (along rows) of the habitat under the LESS approach.
#' @param nDetectorsLESS A \code{Vector}  of dimensions n.cells  with the number of detectors assigned to each cells of the habitat under the LESS approach.
#' @param ResizeFactor A \code{Numeric} giving scaling factor for the habitat grid cells
#' @param habitatID A \code{Matrix} giving the IDs of each habitat grid cell 
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)

#### Density function 
dbin_LESS_Cached_detfun_hn <- nimbleFunction(run = function( x = double(1),
                                                                pZero = double(0),
                                                                sxy = double(1),
                                                                sigma = double(0),
                                                                detector.xy = double(2),
                                                                trials = double(1),
                                                                detectorIndex = double(2),
                                                                nDetectorsLESS = double(1),
                                                                ResizeFactor = double(0, default = 1),
                                                                habitatID = double(2),
                                                                indicator = double(0, default = 1.0),
                                                                log = integer(0, default = 0)){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ## CHECK INPUT DIMENSIONS
  nDetectors <- length(x)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(sum(x[1:nDetectors]) == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  ## GET NECESSARY INFO
  # n.detectors <- length(x) #NOT NEEDED
  alpha <- -1.0 / (2.0 * sigma * sigma)
  
  ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  
  ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
  index <- detectorIndex[sxyID, 1:nDetectorsLESS[sxyID]]
  
  # ## RECREATE Y ==> MOVE TO NIMBLE MODEL CODE
  # y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
  # if(nbDetections > 0){
  #   for(j in 1:nbDetections){
  #     y[yDets[j]] <- x[j]
  #     ## check if a detection is out of the "detection window"
  #     if(sum(yDets[j]==index)==0){
  #       if(log == 0) return(0.0)
  #       else return(-Inf)
  #     }
  #   }
  # }
  
  ## CALCULATE THE LIKELIHOOD OF THE DETECTION VECTOR
  logProb <- 0.0 
  count <- 1 
  index1 <- c(index, 0) # so the count is not an issue
  
  for(j in 1:nDetectors){
    if(x[j] > 0){
      ## Check if the detection is out of the "detection window"
      if(sum(j == index) == 0){
        if(log == 0) return(0.0) else return(-Inf)
      }
    }
    ## If j is equal to the relevant detector
    if(index1[count] == j){
      d2 <- pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2)
      p <- pZero * exp(alpha * d2)
      logProb <- logProb + dbinom(x[j], prob = p, size = trials[j], log = TRUE)
      count <- count + 1
    }
  }#j
  
  if(log)return(logProb)
  else return(exp(logProb))
})

#### Sampling function 
rbin_LESS_Cached_detfun_hn <- nimbleFunction(run = function( n = integer(0),
                                                                pZero = double(0),
                                                                sxy = double(1),
                                                                sigma = double(0),
                                                                detector.xy = double(2),
                                                                trials = double(1),
                                                                detectorIndex = double(2),
                                                                nDetectorsLESS = double(1),
                                                                ResizeFactor = double(0, default = 1),
                                                                habitatID = double(2),
                                                                indicator = double(0, default = 1.0)){
  # RETURN TYPE DECLARATION
  returnType(double(1))
  
  ## GET NECESSARY INFO
  alpha <- -1.0 / (2.0 * sigma * sigma)
  n.detectors <- dim(detector.xy)[1]
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){return(rep(0.0, n.detectors))}
  
  ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  
  ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
  index <- detectorIndex[sxyID, 1:nDetectorsLESS[sxyID]]
  
  ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS 
  detectOut <- rep(0, n.detectors)
  
  ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
  for(j in 1:length(index)){
    # for(j in index){
    d2 <- pow(detector.xy[index[j],1] - sxy[1], 2) + pow(detector.xy[index[j],2] - sxy[2], 2)
    p <- pZero * exp(alpha * d2)
    # Draw the observation at detector j from a binomial distribution with probability p
    detectOut[index[j]] <- rbinom(1, trials[index[j]], p)    
  }#j
  
  ## OUTPUT
  return(detectOut)
})

#### Registration 
registerDistributions(list(
  dbin_LESS_Cached_detfun_hn = list(
    BUGSdist = "dbin_LESS_Cached_detfun_hn(pZero, sxy, sigma,
    detector.xy, trials, detectorIndex, nDetectorsLESS,
    ResizeFactor, habitatID, indicator)",
    types = c( "value = double(1)", "pZero = double(0)","sxy = double(1)", "sigma = double(0)",
               "detector.xy = double(2)", "trials = double(1)", "detectorIndex = double(2)", "nDetectorsLESS = double(1)",
               "ResizeFactor = double(0)", "habitatID= double(2)", "indicator = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
    )))




###############################################################
##################### dbin_LESS_Cached_detfun_hnp #############
###############################################################
#' @title Function to create a NIMBLE custom distribution for faster SCR model runs with respect to half-normal plateau detection function.
#'
#' @description
#' \code{dbin_LESS_Cached_detfun_hnplateau} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] with respect to half-normal plateau detection function
#' 
#' @param x \code{Vector} of length n.detectors containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal plateau detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal plateau detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param detector.xy A \code{Matrix}  of dimensions n.detectors x 2 with the locations of each detector.
#' @param trials A \code{Vector} giving the number of detection sessions for each detector.
#' @param detectorIndex A \code{Matrix}  with the detectors ID (along columns) assigned to each cells (along rows) of the habitat under the LESS approach.
#' @param nDetectorsLESS A \code{Vector}  of dimensions n.cells  with the number of detectors assigned to each cells of the habitat under the LESS approach.
#' @param ResizeFactor A \code{Numeric} giving scaling factor for the habitat grid cells
#' @param habitatID A \code{Matrix} giving the IDs of each habitat grid cell 
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param wd A \code{Numeric} giving the length of the plateaufor the half-normal plateau detection function
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)

#### Density function
dbin_LESS_Cached_detfun_hnplateau <- nimbleFunction(run = function( x = double(1),
                                                                       pZero = double(0),
                                                                       sxy = double(1),
                                                                       sigma = double(0),
                                                                       detector.xy = double(2),
                                                                       trials = double(1),
                                                                       detectorIndex = double(2),
                                                                       nDetectorsLESS = double(1),
                                                                       ResizeFactor = double(0, default = 1),
                                                                       habitatID = double(2),
                                                                       indicator = double(0, default = 1.0),
                                                                       wd = double(0, default = 3.0),
                                                                       log = integer(0, default = 0)){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ## CHECK INPUT DIMENSIONS
  nDetectors <- length(x)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(sum(x[1:nDetectors]) == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  ## GET NECESSARY INFO
  # n.detectors <- length(x) #NOT NEEDED
  alpha <- -1.0 / (2.0 * sigma * sigma)
  
  ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  
  ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
  index <- detectorIndex[sxyID, 1:nDetectorsLESS[sxyID]]
  
  # ## RECREATE Y ==> MOVE TO NIMBLE MODEL CODE
  # y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
  # if(nbDetections > 0){
  #   for(j in 1:nbDetections){
  #     y[yDets[j]] <- x[j]
  #     ## check if a detection is out of the "detection window"
  #     if(sum(yDets[j]==index)==0){
  #       if(log == 0) return(0.0)
  #       else return(-Inf)
  #     }
  #   }
  # }
  
  ## CALCULATE THE LIKELIHOOD OF THE DETECTION VECTOR
  logProb <- 0.0 
  count <- 1 
  index1 <- c(index, 0) # so the count is not an issue
  
  for(j in 1:nDetectors){
    if(x[j] > 0){
      ## Check if the detection is out of the "detection window"
      if(sum(j == index) == 0){
        if(log == 0) return(0.0) else return(-Inf)
      }
    }
    ## If j is equal to the relevant detector
    if(index1[count] == j){
      # d2 <- pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2)
      # p <- pZero * exp(alpha * d2)
      D <- pow(pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2), 0.5)
      if(D <= wd) dens <- 1.0
      if(D >  wd) dens <- exp(alpha *(D-wd)*(D-wd))#/(2.0*sigma*sigma))
      p <- pZero * dens 
      
      logProb <- logProb + dbinom(x[j], prob = p, size = trials[j], log = TRUE)
      count <- count + 1
    }
  }#j
  
  if(log)return(logProb)
  else return(exp(logProb))
})

#### Sampling function 
rbin_LESS_Cached_detfun_hnplateau <- nimbleFunction(run = function( n = integer(0),
                                                                       pZero = double(0),
                                                                       sxy = double(1),
                                                                       sigma = double(0),
                                                                       detector.xy = double(2),
                                                                       trials = double(1),
                                                                       detectorIndex = double(2),
                                                                       nDetectorsLESS = double(1),
                                                                       ResizeFactor = double(0, default = 1),
                                                                       habitatID = double(2),
                                                                       indicator = double(0, default = 1.0),
                                                                       wd = double(0, default = 3.0)){
  # RETURN TYPE DECLARATION
  returnType(double(1))
  
  ## GET NECESSARY INFO
  alpha <- -1.0 / (2.0 * sigma * sigma)
  n.detectors <- dim(detector.xy)[1]
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){return(rep(0.0, n.detectors))}
  
  ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  
  ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
  index <- detectorIndex[sxyID, 1:nDetectorsLESS[sxyID]]
  
  ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS 
  detectOut <- rep(0, n.detectors)
  
  ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
  for(j in 1:length(index)){
    # for(j in index){
    # d2 <- pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2)
    # p <- pZero * exp(alpha * d2)
    D <- pow(pow(detector.xy[index[j],1] - sxy[1], 2) + pow(detector.xy[index[j],2] - sxy[2], 2), 0.5)
    if(D <= wd) dens <- 1.0
    if(D >  wd) dens <- exp(alpha *(D-wd)*(D-wd))#/(2.0*sigma*sigma))
    p <- pZero * dens 
    # Draw the observation at detector j from a binomial distribution with probability p
    detectOut[index[j]] <- rbinom(1, trials[index[j]], p)    
  }#j
  
  ## OUTPUT
  return(detectOut)
})

#### Registration 
registerDistributions(list(
  dbin_LESS_Cached_detfun_hnplateau = list(
    BUGSdist = "dbin_LESS_Cached_detfun_hnplateau(pZero, sxy, sigma,
    detector.xy, trials, detectorIndex, nDetectorsLESS,
    ResizeFactor, habitatID, indicator, wd)",
    types = c( "value = double(1)", "pZero = double(0)","sxy = double(1)", "sigma = double(0)",
               "detector.xy = double(2)", "trials = double(1)", "detectorIndex = double(2)", "nDetectorsLESS = double(1)",
               "ResizeFactor = double(0)", "habitatID= double(2)", "indicator = double(0)",
               "wd = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )))

###############################################################
##################### dbin_LESS_Cached_detfun_ex ##############
###############################################################
#' @title Function to create a NIMBLE custom distribution for faster SCR model runs with respect to exponential detection function.
#'
#' @description
#' \code{dbin_LESS_Cached_detfun_ex} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] with respect to exponential detection function
#' 
#' @param x \code{Vector} of length n.detectors containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the exponential detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the exponential detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param detector.xy A \code{Matrix}  of dimensions n.detectors x 2 with the locations of each detector.
#' @param trials A \code{Vector} giving the number of detection sessions for each detector.
#' @param detectorIndex A \code{Matrix}  with the detectors ID (along columns) assigned to each cells (along rows) of the habitat under the LESS approach.
#' @param nDetectorsLESS A \code{Vector}  of dimensions n.cells  with the number of detectors assigned to each cells of the habitat under the LESS approach.
#' @param ResizeFactor A \code{Numeric} giving scaling factor for the habitat grid cells
#' @param habitatID A \code{Matrix} giving the IDs of each habitat grid cell 
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)

#### Density function 
dbin_LESS_Cached_detfun_ex <- nimbleFunction(run = function( x = double(1),
                                                             pZero = double(0),
                                                             sxy = double(1),
                                                             sigma = double(0),
                                                             detector.xy = double(2),
                                                             trials = double(1),
                                                             detectorIndex = double(2),
                                                             nDetectorsLESS = double(1),
                                                             ResizeFactor = double(0, default = 1),
                                                             habitatID = double(2),
                                                             indicator = double(0, default = 1.0),
                                                             log = integer(0, default = 0)){
  # RETURN TYPE DECLARATION
  returnType(double(0))
  
  ## CHECK INPUT DIMENSIONS
  nDetectors <- length(x)
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){
    if(sum(x[1:nDetectors]) == 0){
      if(log == 0) return(1.0)
      else return(0.0)
    }else{
      if(log == 0) return(0.0)
      else return(-Inf)
    }
  }
  
  ## GET NECESSARY INFO
  # n.detectors <- length(x) #NOT NEEDED
  alpha <- -1.0 / sigma # (2.0 * sigma * sigma)
  
  ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  
  ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
  index <- detectorIndex[sxyID, 1:nDetectorsLESS[sxyID]]
  
  # ## RECREATE Y ==> MOVE TO NIMBLE MODEL CODE
  # y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
  # if(nbDetections > 0){
  #   for(j in 1:nbDetections){
  #     y[yDets[j]] <- x[j]
  #     ## check if a detection is out of the "detection window"
  #     if(sum(yDets[j]==index)==0){
  #       if(log == 0) return(0.0)
  #       else return(-Inf)
  #     }
  #   }
  # }
  
  ## CALCULATE THE LIKELIHOOD OF THE DETECTION VECTOR
  logProb <- 0.0 
  count <- 1 
  index1 <- c(index, 0) # so the count is not an issue
  
  for(j in 1:nDetectors){
    if(x[j] > 0){
      ## Check if the detection is out of the "detection window"
      if(sum(j == index) == 0){
        if(log == 0) return(0.0) else return(-Inf)
      }
    }
    ## If j is equal to the relevant detector
    if(index1[count] == j){
      D <- pow(pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2), 0.5)
      p <- pZero * exp(alpha * D)
      logProb <- logProb + dbinom(x[j], prob = p, size = trials[j], log = TRUE)
      count <- count + 1
    }
  }#j
  
  if(log)return(logProb)
  else return(exp(logProb))
})

#### Sampling function 
rbin_LESS_Cached_detfun_ex <- nimbleFunction(run = function( n = integer(0),
                                                             pZero = double(0),
                                                             sxy = double(1),
                                                             sigma = double(0),
                                                             detector.xy = double(2),
                                                             trials = double(1),
                                                             detectorIndex = double(2),
                                                             nDetectorsLESS = double(1),
                                                             ResizeFactor = double(0, default = 1),
                                                             habitatID = double(2),
                                                             indicator = double(0, default = 1.0)){
  # RETURN TYPE DECLARATION
  returnType(double(1))
  
  ## GET NECESSARY INFO
  alpha <- -1.0 / sigma # (2.0 * sigma * sigma)
  n.detectors <- dim(detector.xy)[1]
  
  ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
  if(indicator == 0){return(rep(0.0, n.detectors))}
  
  ## RETRIEVE THE ID OF THE HABITAT WINDOW THE CURRENT sxy FALLS IN FROM THE HABITAT_ID MATRIX
  sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
  
  ## RETRIEVE THE IDs OF THE RELEVANT DETECTORS
  index <- detectorIndex[sxyID, 1:nDetectorsLESS[sxyID]]
  
  ## INITIALIZE THE OUTPUT VECTOR OF DETECTIONS 
  detectOut <- rep(0, n.detectors)
  
  ## SAMPLE THE DETECTION HISTORY (FOR RELEVANT DETECTORS ONLY)
  for(j in 1:length(index)){
    # for(j in index){
    D <- pow(pow(detector.xy[index[j],1] - sxy[1], 2) + pow(detector.xy[index[j],2] - sxy[2], 2), 0.5)
    p <- pZero * exp(alpha * D)
    # Draw the observation at detector j from a binomial distribution with probability p
    detectOut[index[j]] <- rbinom(1, trials[index[j]], p)    
  }#j
  
  ## OUTPUT
  return(detectOut)
})

#### Registration 
registerDistributions(list(
  dbin_LESS_Cached_detfun_ex = list(
    BUGSdist = "dbin_LESS_Cached_detfun_ex(pZero, sxy, sigma,
    detector.xy, trials, detectorIndex, nDetectorsLESS,
    ResizeFactor, habitatID, indicator)",
    types = c( "value = double(1)", "pZero = double(0)","sxy = double(1)", "sigma = double(0)",
               "detector.xy = double(2)", "trials = double(1)", "detectorIndex = double(2)", "nDetectorsLESS = double(1)",
               "ResizeFactor = double(0)", "habitatID= double(2)", "indicator = double(0)"),
    pqAvail = FALSE,
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )))

###############################################################
##################### HOME RANGE AREA - HN Plateau ############
###############################################################
#' @title Function to compute home range radius and area for a specified circular quantile with respect to half-normal plateau detection function.
#'
#' @description
#' \code{HRA_HNP_nimble_array} returns the estimates of home range radius and area for a given array of parameters with respect to half-normal plateau detection function using bisection algorithm
#' 
#' @param xlim \code{Vector} of length 2 giving the range along x-axis (used in setup part of function)
#' @param ylim \code{Vector} of length 2 giving the range along y-axis (used in setup part of function)
#' @param ng \code{Numeric} variable denoting the number of breaks along an axis (used in setup part of function)
#' @param tol A \code{Numeric} variable denoting the allowed tolerance in the radius estimate (used in setup part of function)
#' @param niter A \code{Numeric}  giving the maximum number of iterations in bisection algorithm (used in setup part of function)
#' @param p0Vec A \code{Vector}  of values for 'p0' variable. (used in run part of function)
#' @param sigmaVec A \code{Vector}  of values for 'sigma' variable. (used in run part of function)
#' @param wdVec  A \code{Vector}  of values for 'wd' variable. (used in run part of function)
#' @param pr A \code{Numeric}  variable denoting the quantile probability in order to compute the radius (used in run part of function)
#' @param d A \code{Numeric} giving an inital value of the radius (used in run part of function)

#Bisection algorithm
HRA_HNP_nimble_array <- nimbleFunction(
  setup = function(xlim, ylim, ng, tol, niter){
    s <- c(sum(xlim)/2, sum(ylim)/2) # center of the rectangular grid
    x1 <- rep(seq(xlim[1], xlim[2], length.out = ng), ng) # vector (ng^2x 1)
    x2.temp <- seq(ylim[1], ylim[2], length.out = ng)
    x2 <- nimNumeric(length = ng^2, value = 0, init = TRUE)
    for(i in 1:ng){
      x2[((i-1)*ng+1):(i*ng)] <- rep(x2.temp[i], ng)
    }
    
    n <- length(x1) # ng^2
    D <- sqrt((s[1] - x1[1:n])^2 + (s[2] - x2[1:n])^2 ) # vector (n x 1)
    
  },
  run = function( p0Vec = double(1), 
                  sigmaVec = double(1),
                  wdVec = double(1),
                  pr = double(0),
                  d = double(0)){
 
    mcmc.niter <- length(p0Vec)
    out <- matrix(NA, mcmc.niter, 2)
    
    for(mcmc.iter in 1:mcmc.niter ){
      if (mcmc.iter%%500 == 0)  cat('..... MCMC iter #', mcmc.iter, '\n')
      
      p0 <- p0Vec[mcmc.iter]
      sigma <- sigmaVec[mcmc.iter]
      wd <- wdVec[mcmc.iter]
      
      p <- nimNumeric(length = n, value = p0, init = TRUE)
      for(i in 1:n){
        if(D[i] > wd) {p[i] <- p0*exp(-(D[i] - wd)*(D[i] - wd)/(2*sigma*sigma))}
      }
      psi <- p[1:n]/sum(p[1:n]) # vector (n x 1)
      
      d.lower <- min(D)
      d.upper <- max(D)
      
      temp <- 0
      for(i in 1:n){
        if(D[i] <= d) temp <- temp + psi[i]
      }
      iter <- 1 #0
      
      while(iter <= niter &
            abs(d.upper-d.lower) > d*tol){
        
        if(temp >= pr){
          d <- (d.lower+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp < pr) d.lower <- d
        }
        if(temp < pr){
          d <- (d.upper+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp >= pr) d.upper <- d
        }
        iter <- iter + 1
        
      }#while
      if(temp < pr) {d <- d*(1+tol) }
      radius <- d
      area <- pi * radius^2
      out[mcmc.iter, 1:2] <- c(radius, area)
      
    }#mcmc.iter
    
    returnType(double(2))
    
    return(out) 
  })



###############################################################
##################### HOME RANGE AREA - Exponential ###########
###############################################################
#' @title Function to compute home range radius and area for a specified circular quantile with respect to exponential detection function.
#'
#' @description
#' \code{HRA_EX_nimble_array} returns the estimates of home range radius and area for a given array of parameters with respect to exponential detection function using bisection algorithm
#' 
#' @param xlim \code{Vector} of length 2 giving the range along x-axis (used in setup part of function)
#' @param ylim \code{Vector} of length 2 giving the range along y-axis (used in setup part of function)
#' @param ng \code{Numeric} variable denoting the number of breaks along an axis (used in setup part of function)
#' @param tol A \code{Numeric} variable denoting the allowed tolerance in the radius estimate (used in setup part of function)
#' @param niter A \code{Numeric}  giving the maximum number of iterations in bisection algorithm (used in setup part of function)
#' @param p0Vec A \code{Vector}  of values for 'p0' variable. (used in run part of function)
#' @param sigmaVec A \code{Vector}  of values for 'sigma' variable. (used in run part of function)
#' @param pr A \code{Numeric}  variable denoting the quantile probability in order to compute the radius (used in run part of function)
#' @param d A \code{Numeric} giving an inital value of the radius (used in run part of function)

#Bisection algorithm
HRA_EX_nimble_array <- nimbleFunction(
  setup = function(xlim, ylim, ng, tol, niter){
    s <- c(sum(xlim)/2, sum(ylim)/2) # center of the rectangular grid
    x1 <- rep(seq(xlim[1], xlim[2], length.out = ng), ng) # vector (ng^2x 1)
    x2.temp <- seq(ylim[1], ylim[2], length.out = ng)
    x2 <- nimNumeric(length = ng^2, value = 0, init = TRUE)
    for(i in 1:ng){
      x2[((i-1)*ng+1):(i*ng)] <- rep(x2.temp[i], ng)
    }
    
    n <- length(x1) # ng^2
    D <- sqrt((s[1] - x1[1:n])^2 + (s[2] - x2[1:n])^2 ) # vector (n x 1)
    
  },
  run = function( p0Vec = double(1), 
                  sigmaVec = double(1),
                  pr = double(0),
                  d = double(0)){
    
    mcmc.niter <- length(p0Vec)
    out <- matrix(NA, mcmc.niter, 2)
    
    for(mcmc.iter in 1:mcmc.niter ){
      if (mcmc.iter%%500 == 0)  cat('..... MCMC iter #', mcmc.iter, '\n')
      
      p0 <- p0Vec[mcmc.iter]
      sigma <- sigmaVec[mcmc.iter]
      
      p <- p0*exp(-D[1:n]/sigma)
      psi <- p[1:n]/sum(p[1:n]) # vector (n x 1)
      
      d.lower <- min(D)
      d.upper <- max(D)
      
      temp <- 0
      for(i in 1:n){
        if(D[i] <= d) temp <- temp + psi[i]
      }
      iter <- 1 #0
      
      while(iter <= niter &
            abs(d.upper-d.lower) > d*tol){
        
        if(temp >= pr){
          d <- (d.lower+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp < pr) d.lower <- d
        }
        if(temp < pr){
          d <- (d.upper+d)/2
          temp <- 0
          for(i in 1:n){
            if(D[i] <= d) temp <- temp + psi[i]
          }
          if(temp >= pr) d.upper <- d
        }
        iter <- iter + 1
        
      }#while
      if(temp < pr) {d <- d*(1+tol) }
      radius <- d
      area <- pi * radius^2
      out[mcmc.iter, 1:2] <- c(radius, area)
      
    }#mcmc.iter
    
    returnType(double(2))
    
    return(out) 
  })


#############################################################################
##################### HOME RANGE AREA - any detection function ##############
#############################################################################
#' @title Function to compute home range radius and area for a specified circular quantile with respect to a given detection function.
#'
#' @description
#' \code{HRA_detfun} returns the estimates of home range radius and area for a given array of parameters with respect to a given detection function using bisection algorithm
#' 
#' @param detfun \code{Character} value denoting the type of detection function in the SCR model.
#'  Avalaible choices are: "HN", "EX", "HNP", "AL", "DN", "BI".
#' @param param A \code{Vector} of values for detection function parameters. Names of the different values should be according the corresponding parameter.
#' @param pr A \code{Numeric}  variable denoting the quantile probability in order to compute the radius
#' @param xlim \code{Vector} of length 2 giving the range along x-axis
#' @param ylim \code{Vector} of length 2 giving the range along y-axis
#' @param ng \code{Numeric} variable denoting the number of breaks along an axis
#' @param tol A \code{Numeric} variable denoting the allowed tolerance in the radius estimate
#' @param niter A \code{Numeric}  giving the maximum number of iterations in bisection algorithm

# Home range area for any detection function
# Bisection algorithm
HRA_detfun <- function (detfun, param, pr = 0.95, xlim, ylim, 
                        ng = 1200, tol = 1E-5, niter = 2000) 
{
  # HALF-NORMAL
  if(detfun == 'HN'){
    sigma <- param['sigma']
    q<-qchisq(pr,2)#--FOLLOWING ROYLE BOOK PAGE 136
    radius<-sigma*sqrt(q)
    area<-pi*radius^2
    out <- c(radius, area)
    names(out) <- c("HRradius", "HRarea")
    return(out)
  }
  
  if(detfun != 'HN'){
  # HALF-NORMAL PLATEAU
  if(detfun == 'HNP'){
    func <- function(param, D){ 
    sigma <- param['sigma']; p0 <- param['p0']; wd <- param['wd']
    dens <- D # rep(NA, length(D))
    dens[D<=wd] <- 1
    dens[D>wd] <- exp(-(D[D>wd]-wd)*(D[D>wd]-wd)/(2.0*sigma*sigma)) 
    return(p0*dens)
  }
  }
  # DONUT
  if(detfun == 'DN'){
    func <- function(param, D){ 
    sigma <- param['sigma']; sigma.b <- param['sigma.b'] 
    p0 <- param['p0']; wd <- param['wd']
    dens <- D # rep(NA, length(D))
    dens[D<=wd] <- exp(-(D[D<=wd]-wd)*(D[D<=wd]-wd)/(2.0*sigma.b*sigma.b))
    dens[D>wd] <- exp(-(D[D>wd]-wd)*(D[D>wd]-wd)/(2.0*sigma*sigma)) 
    return(p0*dens)
  }
  }
  # Asymmetrical logistic 
  if(detfun == 'AL'){
    func <- function(param, D){ 
      sigma <- param['sigma']; p0 <- param['p0']
    slope <- param['slope']; slopecon <- param['slopecon']
    Anuslope <- (2*abs(slope)*slopecon)/(1+slopecon)
    fx <- 1/ (1+(D/sigma)^Anuslope)
    den <- 1+fx*((D/sigma)^(slope))+(1-fx)*((D/sigma)^(slope*slopecon))
    detp <- p0/den
    return(detp)
  }
  }
  # LAPLACE OR EXPONENTIAL
  if(detfun == 'EX'){
    func <- function(param, D){ 
      sigma <- param['sigma']; p0 <- param['p0']
    dens <- exp(-D/sigma)
    return(p0*dens)
  }
  }
  # Bimodal 
  if(detfun == 'BI'){
    func <- function(param, D){ 
    sigma <- param['sigma']  # Scale parameter of the first peak (could be different between the two peaks but less params with one)
    sigma.b <- param['sigma.b']  # Baseline detection probability of the second peak
    p0 <- param['p0']     # Baseline detection probability of the first peak (could be different between the two peaks but less params with one)
    p0.b <- param['p0.b']     # Baseline detection probability of the second peak
    wd <- param['wd']       # Distance between the home-range center and the "doughnut" center
    densa <- p0 *  exp(-D*D/(2.0*sigma*sigma)) 
    densb <-  p0.b * exp(-(D-wd)*(D-wd)/(2.0*sigma.b*sigma.b))
    dens <- densa + densb
    return(dens)
  }
  }
  
  s <- c(sum(xlim)/2, sum(ylim)/2) # center of the rectangular grid
  x1 <- rep(seq(xlim[1], xlim[2], length.out = ng), ng) # vector (ng^2x 1)
  x2 <- sort(rep(seq(ylim[1], ylim[2], length.out = ng), ng)) # vector (ng^2x 1)
  
 
  # delta <- min(diff(x1[1:ng])) # minimum of the lag-1 differences
  # x1 <- rep(seq(xlim[1] - delta/2, xlim[2] + delta/2, delta), ng+1) # vector (ng(ng+1)x 1)
  # delta <- min(diff(sort(unique(x2)))) # minimum of the lag-1 differences
  # x2 <- sort(rep(seq(ylim[1] - delta/2, ylim[2] + delta/2, delta), ng+1)) # vector (ng(ng+1)x 1)
  
  D <- sqrt((s[1] - x1)^2 + (s[2] - x2)^2) # vector (length(x1) x 1)
  p <- func(param, D) # vector (length(x1) x 1)
  
  # Normed detection probability density
  psi <- p/sum(p) # vector (length(x1) x 1)
  
    obj <- function(d){ sum(psi[D <= d]) - pr}
    d.lower <- min(D)
    d.upper <- max(D)
    d <- (d.lower+d.upper)/2 
    del <- obj(d)
    iter <- 0
    repeat {
      iter <- iter + 1
      # print(iter)
      if(del >= 0){
        d <- (d.lower+d)/2
        del <- obj(d)
        if(del < 0) d.lower <- d
      }
      if(del < 0){
        d <- (d.upper+d)/2
        del <- obj(d)
        if(del >= 0) d.upper <- d
      }
      if(iter == niter | abs(d.upper-d.lower) <= d*tol){break}
    }
    if(obj(d) < 0) {d <- d*(1+tol) }
    radius <- d
    area <- pi * radius^2
    out <- c(radius, area)
    names(out) <- c("HRradius", "HRarea")
    return(out) 
  } # if not HN
  
}#HRA_detfun




#############################################################################
##################### SIMULATING POSTERIOR REPLICATES OF THE DATA ###########
#############################################################################
#' @title Function to draw posterior replicates of the data for a given MCMC sample of parameters.
#'
#' @description
#' \code{simulateFromPosterior} returns posterior replicates of the data for a given MCMC sample of parameters.
#' 
#' @param model A \code{model} object (used in setup part of function)
#' @param nodes A \code{Vector} of nodes (used in setup part of function)
#' @param posteriorSamples A \code{matrix} of dimension nIter x nNodes with MCMC samples for each 'nodes' variable along its columns. (used in run part of function)

simulateFromPosterior <- nimbleFunction(
  setup = function(model, nodes) { # 'setup' function  is written is R language
    mv <- modelValues(model)
    allVars <- model$getVarNames()
    downNodes <- model$getDependencies(nodes = nodes, self = FALSE, includeData = TRUE)
  },
  run = function(posteriorSamples = double(2)) { # 'run' function is written is nimble language
    nIter <- dim(posteriorSamples)[1]
    nNodes <- dim(posteriorSamples)[2]
    resize(mv, nIter)
    for(ii in 1:nIter){
      if (ii%%500 == 0)  cat('..... drawing sample #', ii, '\n')
      values(model, nodes) <- posteriorSamples[ii,1:nNodes]
      model$simulate(nodes = downNodes, includeData = TRUE) # Specified 'nodes' (by MCMC samples) should not get updated 
      model$calculate(nodes = downNodes) # Log-likelihood p(yrep | theta_iter) of the simulated data 
      copy(from = model, nodes = allVars, to = mv, rowTo = ii, logProb = TRUE)
    }
  }) 


##############################################################
##################### COMPUTING BAYESIAN P-VALUES ############
##############################################################
#' @title Function to compute Bayesian p-value for a given discrepancy metric.
#'
#' @description
#' \code{pvaluefn} returns the estimates of Bayesian p-value for a given discrepancy metric and a given MCMC sample of parameters and SCR data replicates
#' 
#' @param repData A \code{List} of niter number of matrices, each of which is a replicate of single season SCR data (dimension M x ntrap) based on MCMC sample of the model parameters and latent variables.
#' @param obsData A \code{Matrix} of observed single season SCR data (dimension M x ntrap)
#' @param detprob A \code{List} of niter number of matrices, each of which is detection probability matrix (dimension M x ntrap)  based on MCMC sample of the model parameters and latent variables.
#' @param trials A \code{Vector} giving the number of detection sessions for each detector
#' @param metric value denoting the type of discrepancy metric to compute Bayesian p-value.
#'  Avalaible choices are: 'ALL', FT', 'PearsonsChiSq', 'FT_detector', 'FT_individual'. 
#'  If 'ALL', then all Bayesian p-value is computed with all the four metrics: FT', 'PearsonsChiSq', 'FT_detector', 'FT_individual'.

pvaluefn <- function(repData = NULL,# list (niter) # components: matrix (M x ntrap) 
                     obsData = NULL,# M x ntrap
                     detprob = NULL, # list (niter) # components: matrix (M x ntrap) 
                     trials = NULL, # no. of occasions or trial
                     metric = 'FT')  # c('FT', 'PearsonsChiSq', 'FT_detector', 'FT_individual')
  
{
  
  statnames <- c('FT', 'PearsonsChiSq', 'FT_detector', 'FT_individual') 
                 
  if(sum(metric %in% 'ALL') == 0) metric <- intersect(statnames, metric)
  eps <- 1e-10
  Tstat.sim <- Tstat.obs <- NULL
  pvalue.out <- c()
  
  #------------------
  # FT
  #------------------
  if(sum(metric %in% c('FT', 'ALL')) > 0){
    # disfn <- function(y1, pr){ sum((sqrt(y1) - sqrt(pr*ntrial))^2) }
    
    if(sum(trials)==length(trials)){
      disfn <- function(y1, pr){ 
        A1 <- sqrt(y1)
        C1 <- sum((A1-sqrt(pr))^2)
        return(C1)}
    }else{
      disfn <- function(y1, pr){ 
        A1 <- sqrt(y1)
        B1 <- sqrt(t(apply(pr, 1, function(x,y){x*y}, y = trials)))
        C1 <- sum((A1-B1)^2)
        return(C1)}
    }
    
    Tstat.sim.FT <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.FT <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT)
    if(sum(metric %in% c('FT')) > 0){
      pvalue.FT <- mean(Tstat.sim.FT>Tstat.obs.FT)  
      pvalue.out <- c(pvalue.out, pvalue.FT)
    }
  } #FT
  
 
  #------------------
  # Pearson's ChiSq
  #------------------ 
  if(sum(metric %in% c('PearsonsChiSq', 'ALL')) > 0){
    if(sum(trials)==length(trials)){
      disfn <- function(y1, pr){ 
        A1 <- y1
        #B1 <- t(apply(pr, 1, function(x,y){x*y}, y = trials))
        C1 <- sum(((A1-pr)^2)/(pr+eps))
        return(C1)}
    }else{
      disfn <- function(y1, pr){ 
        A1 <- y1
        B1 <- t(apply(pr, 1, function(x,y){x*y}, y = trials))
        C1 <- sum(((A1-B1)^2)/(B1+eps))
        return(C1)}
    }
    Tstat.sim.omni <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.omni <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.omni)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.omni)
    if(sum(metric %in% c('PearsonsChiSq')) > 0){
      pvalue.omni <- mean(Tstat.sim.omni>Tstat.obs.omni)  
      pvalue.out <- c(pvalue.out, pvalue.omni)
    }
  } #Omnibus
  
  #------------------
  # FT_detector
  #------------------ 
  if(sum(metric %in% c('FT_detector', 'ALL')) > 0){
    # disfn <- function(y1, pr){ sum((sqrt(colSums(y1)) - sqrt(colSums(pr)*ntrial))^2) }
    disfn <- function(y1, pr){ sum((sqrt(colSums(y1)) - sqrt(colSums(pr)*trials))^2) }
    Tstat.sim.FT.det <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.FT.det <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.det)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.det)
    if(sum(metric %in% c('FT_detector')) > 0){
      pvalue.FT.det <- mean(Tstat.sim.FT.det>Tstat.obs.FT.det)  
      pvalue.out <- c(pvalue.out, pvalue.FT.det)
    }
  } #FT_detector
  
  #------------------
  # FT_individual
  #------------------ 
  if(sum(metric %in% c('FT_individual', 'ALL')) > 0){
    if(sum(trials)==length(trials)){
      disfn <- function(y1, pr){ 
        A1 <- sqrt(rowSums(y1))
        C1 <- sum((A1-sqrt(rowSums(pr)))^2)
        return(C1)}
    }else{
      disfn <- function(y1, pr){ 
        A1 <- sqrt(rowSums(y1))
        B1 <- sqrt(rowSums(t(apply(pr, 1, function(x,y){x*y}, y = trials))))
        C1 <- sum((A1-B1)^2)
        return(C1)}
    }
    
    Tstat.sim.FT.ind <- mapply(disfn, y1 = repData, pr = detprob, SIMPLIFY = T) # niter x 1
    Tstat.obs.FT.ind <- unlist(lapply(detprob, disfn, y1 = obsData)) # niter x 1
    Tstat.sim <- cbind(Tstat.sim, Tstat.sim.FT.ind)
    Tstat.obs <- cbind(Tstat.obs, Tstat.obs.FT.ind)
    if(sum(metric %in% c('FT_individual')) > 0){
      pvalue.FT.ind <- mean(Tstat.sim.FT.ind>Tstat.obs.FT.ind)  
      pvalue.out <- c(pvalue.out, pvalue.FT.ind)
    }
  } #FT_individual
  
  #=========================================================
  condition <- sum(metric %in% 'ALL') == 0 && sum(metric %in% statnames) < length(statnames)
  if(condition){ names(pvalue.out) <- metric }
  if(metric == 'ALL' || sum(metric %in% statnames) == length(statnames)){
    pvalue.out <- mapply(function(tsim,tobs){mean(tsim>tobs)},
                         tsim = split2list(Tstat.sim, 2),
                         tobs = split2list(Tstat.obs, 2), SIMPLIFY = T)
    names(pvalue.out) <- statnames
  }
  # out <- list(Tstat.sim = Tstat.sim, Tstat.obs = Tstat.obs, pvalue = pvalue.out)
  out <- list(pvalue = pvalue.out)
}# pvaluefn function

split2list <- function(x, margin){ # x is an array # margin is a scalar
  if(length(dim(x))>2)  xx <- lapply(split(x, slice.index(x, margin)), function(a){array(a, dim(x)[-margin])})
  if(length(dim(x))==2 && margin == 1)  xx <- lapply(snow::splitRows(x, dim(x)[[1]]), function(a){c(a)}) 
  if(length(dim(x))==2 && margin == 2)  xx <- lapply(snow::splitCols(x, dim(x)[[2]]), function(a){c(a)})
  if(is.vector(x)) xx<- split(x, cut(x, length(x), labels = FALSE)) 
  names(xx) <- NULL
  return(xx)
}

################################################################################################### 
##################### RUNNING THE MCMC, computing HR radius and Bayesian p-values #################
###################################################################################################

#' @title Function for simulation of SCR data set, run MCMC, compute home range radius and area, Bayesian p-value for the simulated SCR data set with given detection function for simulation and model fitting
#'
#' @description
#' \code{runDetFun} returns the MCMC samples of monitored variables (myNimbleOutput), summary statistics of the parameters of ineterest (summary.output),
#' samples of home range radius and area (HR.mcmc.sample) for a given quantile probability and MCMC samples of model parameters, estimates of home range radius and area corresponding to the
#' detection function used for SCR data simulation  (HR.true), the quantile probability to calculate home range radius and area (pr), Bayesian p-values (pval),
#' MCMC run time (Runtime), numer of MCMC iterations (niter), burnin period (burnin), no. of parallel MCMC chains (nchains).
#' 
#' @param detfunsim \code{Vector} containing detection/non detections 
#' @param detfunfit \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param burnin \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param nchains A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param pr A \code{Numeric}  variable denoting the quantile probability in order to compute the radius (used in run part of function)
#' @param WD A \code{Character} giving the working directory where the output will be saved.
#' @param myVars A \code{List} giving various parameters for the simulation: HABITAT (extent, resolution, buffer), DETECTORS (resolution),
#' POPULATION (N - population size, M - an upper bound of the population size) and detection function parameters.
#' @param plot.save.check A \code{LOGICAL} variable to indicate whether the plots need to be saved in the working directory (WD) or not.

runDetFun <- function(
  detfunsim = "HN"
  , detfunfit = "EX"
  , burnin = 1
  , niter =   100 #10000 # 20000 # 10000
  , nchains =  2
  , pr = 0.95
  , WD = NA
  , myVars = NULL
  , plot.save.check = FALSE){ ##--DO ALL
  
  if(is.na(WD)){WD <- getwd()}
  if(is.null(myVars)){
    myVars <- list(
      HABITAT =    list(extent     = c(20),
                        resolution = 1,
                        buffer     = 4.5),     
      DETECTORS =  list(resolution = 1),       
      POPULATION = list(N = 200,               # true population size
                        M = 400),              # size of the augmented population
      HN = list(p0  = 0.3,
                sigma = 1.5), 
      HNP = list(p0  = 0.25,
                 sigma = 1,
                 wd = 1.5), 
      DN = list(p0 = 0.25, 
                sigma = 1, 
                sigma.b = 1.5,
                wd = 1.5), 
      AL = list(p0 = 0.3, 
                sigma = 2,
                slope = 5, 
                slopecon = 1), 
      EX = list(p0 = 0.3, 
                sigma = 1.5),
      BI = list(p0 = 0.25, 
                sigma = 0.5, 
                p0.b = 0.15, 
                sigma.b = 1, 
                wd = 2) 
    )
  }
  
  require(rgeos)                     # Import the geospatial analysis libraries
  require(rgdal)                     # Import the spatial input/ouput libraries
  require(raster)                    # Import the raster spatial package
  require(coda)                      # Import the MCMC diagnostic tools
  require(nimble)                    # Import the NIMBLE subroutines
  require(nimbleSCR)
  require(abind)                     # Import the library for manipulating multidimensional arrays
  
  modelName = paste("DetFun_sim", detfunsim,"_fit",detfunfit, sep = "")## [SD] 
  dir.create(file.path(WD, modelName))
  
  if(plot.save.check == TRUE){
    graphics.off()
    path <- file.path(WD,modelName, paste("PLOTS_sim", detfunsim,"_fit",detfunfit, ".pdf", sep = ""))
    pdf(file=path, width = 10, height = 10)
  }
  ## ----------------------------------------------------------------------------------------------
  ## ------ 1. SET-UP HABITAT AND DETECTORS -----
  ## ----------------------------------------------------------------------------------------------
  
  ### ==== 1.1. GENERATE HABITAT ====
  buffer <- myVars$HABITAT$buffer
  grid.size <- c(myVars$HABITAT$extent)
  
  coords <- matrix(c(buffer            , buffer ,
                     grid.size + buffer, buffer ,
                     grid.size + buffer, grid.size + buffer,
                     buffer            , grid.size + buffer,
                     buffer            , buffer 
  ), ncol = 2, byrow = TRUE)
  
  P1 <- Polygon(coords)
  myStudyArea <- SpatialPolygons(list(Polygons(list(P1), ID = "a")),
                                 proj4string = CRS("+proj=utm +zone=33 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"))

  poly <- myStudyArea
  resolution <- myVars$HABITAT$resolution
  
  ## ----- Store the original study area polygon -----
  polygon.orig <- poly
  
  ## ----- Create the buffer zone around the study area -----
  buff <- rgeos::gBuffer(poly, width = buffer)
  poly <- raster::union(poly, buff)
  
  ## ----- Create Habitat raster (study area + buffer)-----   
  r <- raster(extent(poly))                                         ## Rasterise polygon
  res(r) <- resolution  ## Change resolution
  
  polygon.r <- rasterize(as(poly, "SpatialLines"), r, field=1)      ## rasterise as lines and then polygons so all cells touched by the polygon are covered 
  polygon.r <- rasterize(poly, polygon.r, update=T, field=1)        ## Add field=1 so its 1 for study area and NA for the rest.
  
  ## ----- Create Habitat matrix (habitat : 1 and non-habitat: 0) -----
  habitat.r <- polygon.r
  habitat.r[is.na(habitat.r)] <- 0                                     ## Give 0 values to raster cells outside the study area
  habitat.mx <- as.matrix(habitat.r)                                   ## Convert to matrix
  
  # ## ----- Give unique IDs to cells ----- 
  IDCells.r <- polygon.r
  IDCells.r[] <- 1:length(IDCells.r)							         ## Cell ID starts from the top left corner
  IDCells.mx <- as.matrix(IDCells.r)							               ## Convert to matrix
  
  ## ----- Obtain xy coordinates of cells -----   
  habitat.xy <- xyFromCell(polygon.r, 1:ncell(polygon.r))
  dimnames(habitat.xy) <- list(1:length(habitat.xy[,"x"]), c("x","y"))
  habitat.sp <- SpatialPointsDataFrame(data.frame(habitat.xy[,c("x","y")]), data=data.frame(habitat.xy), proj4string=CRS(projection(poly)))
  habitat.index <- which(!is.na(over(as(habitat.sp, "SpatialPoints"), as(poly,"SpatialPolygons"))))
  habitat.clip.sp <- habitat.sp[!is.na(over(as(habitat.sp, "SpatialPoints"), as(poly, "SpatialPolygons"))), ]
  
  ## -- Obtain lower and upper cell coordinates
  lower.hab.sp <- data.frame(coordinates(habitat.sp) - resolution/2)
  upper.hab.sp <- data.frame(coordinates(habitat.sp) + resolution/2)
  coordinates(upper.hab.sp) <- upper.hab.sp
  coordinates(lower.hab.sp) <- lower.hab.sp
  proj4string(lower.hab.sp) <- CRS(projection(poly))
  proj4string(upper.hab.sp) <- CRS(projection(poly))
  
  myHabitat <- list(habitat.sp = habitat.sp, #
                    habitat.clip.sp = habitat.clip.sp,
                    habitat.xy = habitat.xy,
                    IDCells.mx = IDCells.mx,#
                    habitat.r = habitat.r, #
                    habitat.mx = habitat.mx,
                    resolution = resolution,
                    buffered.habitat.poly = poly,#
                    habitat.poly = polygon.orig,#
                    habitat.index = habitat.index,
                    upper.hab.sp = upper.hab.sp, ##
                    lower.hab.sp = lower.hab.sp ##
  )
  
  par(mfrow=c(1,1))
  
  plot(habitat.r, main = 'Buffered habitat', legend = FALSE) ## The habitat and polygons 
  plot(poly, add=TRUE, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 0.5))						
  plot(polygon.orig, add=TRUE, col = rgb(red = 102/255,green = 102/255,blue = 102/255,alpha = 1))
 
  lowerHabCoords <- coordinates(myHabitat$lower.hab.sp)
  upperHabCoords <- coordinates(myHabitat$upper.hab.sp)
  nHabCells <- dim(lowerHabCoords)[1]
  habQualityDims <- dim(myHabitat$habitat.r)[2:1]
  
  ### ====  GENERATE DETECTORS ====
  
  data <- myStudyArea
  resolution.det <- myVars$DETECTORS$resolution
  center <- TRUE
  
  detector.r1 <- raster(extent(data), resolution = resolution.det, crs = proj4string(data))
  maindetector.r <- rasterize(data, detector.r1)
  
  ### ==== CREATE POLYGONS FROM RASTER  ====
  ## Main detectors 
  temp.r <- raster(maindetector.r)
  temp.r[] <- 1:length(temp.r)
  maindetector.poly <- rasterToPolygons(temp.r, dissolve = TRUE)
  
  ### ==== OBTAIN SPATIALPOINTS FROM DETECTORS ====
  ## Main detectors 
  main.detector.xy <- xyFromCell(maindetector.r, 1:ncell(maindetector.r))
  main.detector.sp <- SpatialPointsDataFrame( data.frame(main.detector.xy[,c("x","y")]),
                                              data=data.frame(main.detector.xy),
                                              proj4string=CRS(projection(data)))
  names(main.detector.sp@data) <- c("main.cell.x","main.cell.y")
  main.detector.sp@data$main.cell.id <- 1:length(main.detector.sp)
  
  ### ==== MAKE A PLOTTING FUNCTION ====
  
  plot(data, col=grey(0.3), main = 'DETECTORS')
  plot(maindetector.poly, add=TRUE, lwd=3)
  plot(main.detector.sp, pch=19, cex=1, col=as.numeric(main.detector.sp@data$main.cell.id), add=TRUE)
  
  myDetectors <- list( detector.sp = main.detector.sp, 
                       grid.poly = maindetector.poly)
  
  ### ==== SCALE DETECTORS & HABITAT & UPPER/LOWER COORDINATES FOR HABITAT WINDOWS ====
  
  # SCALE DETECTORS
  scaled <- UTMToGrid(data.sp = myDetectors$detector.sp,
                      grid.sp = myHabitat$habitat.sp)
  ScaledLowerCoords <- UTMToGrid(data.sp = SpatialPoints(lowerHabCoords),
                                 grid.sp = myHabitat$habitat.sp)$data.scaled.xy
  ScaledUpperCoords <- UTMToGrid(data.sp = SpatialPoints(upperHabCoords),
                                 grid.sp = myHabitat$habitat.sp)$data.scaled.xy
  
  ScaledUpperCoords[,2] <- ScaledUpperCoords[,2]+1
  ScaledLowerCoords[,2] <- ScaledLowerCoords[,2]-1
  
  
  ### ====  CREATE CACHED OBJECTS FOR LESS RESTRICTION ====
  # Creating a set of detector IDs within maxDist radius for each habitat cells
  DetectorIndexLESS <- nimbleSCR::getLocalTraps(
    habitatMask=myHabitat$habitat.mx,
    trapCoords=scaled$data.scaled.xy,
    dmax=15/res(myHabitat$habitat.r)[1],
    resizeFactor = 1,
    plot.check = TRUE
  )
  
  ### ==== SIMULATE INDIVIDUAL AC LOCATIONS ====
  ## Sample random (uniform) activity center locations
  ## sample ACs from the the same habitat extent than the one used in the model (based on habitat.r; i.e. not the rounded edges of the bufferered habitat polygon).
  simulated.ACS <- spsample( x = raster::aggregate(rasterToPolygons(myHabitat$habitat.r,function(x)x>0)),
                             n = myVars$POPULATION$N,
                             type="random")
  simulated.ACS$id <- 1:length(simulated.ACS)
  
  plot(myHabitat$buffered.habitat.poly, main = 'Simulated activity centres')
  plot(raster::aggregate(rasterToPolygons(myHabitat$habitat.r,function(x)x>0)), add=TRUE)
  plot(simulated.ACS, add=TRUE)
  
  ### ==== SIMULATE DETECTION ARRAY : y.ar ====
  
  params <- unlist(myVars[[detfunsim]])
  
  detectionSimulationOut <- SimulateDetection_detfun( 
    params = params,
    AC.sp = simulated.ACS,
    detector.sp =  myDetectors$detector.sp,
    detfun.type = detfunsim, 
    plot.check = T
  ) #,#)#,
  title(main = 'Simulated SCR data')
  plot(myStudyArea, add=TRUE)
  plot(raster::aggregate(rasterToPolygons(myHabitat$habitat.r,function(x)x>0)), add=TRUE)
  
  ### ====  SUBSET TO INDIVIDUALS DETECTED AT LEAST ONE YEAR/DATA SET ==== 
  y <- detectionSimulationOut$y
  
  ## CHECK THE NUMBER OF INIDVIDUALS DETECTED
  detected <- apply( detectionSimulationOut$y.all,1, function(x) sum(x)>0)
  
  
  ### ====  AUGMENT DATA ====
  
  ##--augmentation factor that ensures that the total number of individuals 
  # (alive + available) is always the same, regardless of the simulation
  this.AF <- myVars$POPULATION$M/sum(detected)-1
  
  #---check it: 
  sum(detected) * (1+ this.AF) ==  myVars$POPULATION$M
  
  y <- MakeAugmentation(y = y, aug.factor = this.AF, replace.value = 0)
  z <- MakeAugmentation( y = rep(1, nrow(detectionSimulationOut$y)),
                         aug.factor = this.AF, replace.value = NA)
  
  
  ### ==== TRANSFORM Y TO SPARSE MATRICES =====
  
  SparseY <- nimbleSCR::getSparseY(y)
  
  ### ==== SET THE INPUT FOR NIMBLE ====
  
  ### ==== Define the constants to be used in the model ====
  # Set the list of model constants
  nimConstants <- list(
    x.max = habQualityDims[1],
    y.max =  habQualityDims[2],
    n.detectors = dim(y)[2], 
    n.individuals = dim(y)[1],
    numHabWindows = nHabCells,
    nMaxDetectors = SparseY$maxDetNums, 
    y.maxDet = dim(DetectorIndexLESS$habitatGrid)[1], 
    x.maxDet = dim(DetectorIndexLESS$habitatGrid)[2], 
    ResizeFactor = DetectorIndexLESS$resizeFactor, 
    maxNBDets = DetectorIndexLESS$numLocalTrapsMax,
    n.cells = dim(DetectorIndexLESS$localTrapsIndices)[1],
    MaxDist = params['sigma'] *  5 
  )
  
  ### ==== Define the data to be used in the model ====
  # Set the list of data
  nimData <- list(
    lowerHabCoords = ScaledLowerCoords,
    upperHabCoords = ScaledUpperCoords,
    z = z,
    y = SparseY$y[,,1],
    detector.xy = scaled$data.scaled.xy,
    trials = rep(1, nrow(scaled$data.scaled.xy)),  
    mu = rep(1, nHabCells),
    yDets = SparseY$detIndices[,,1], 
    nbDetections = SparseY$detNums[,1], 
    detectorIndex = DetectorIndexLESS$localTrapsIndices, 
    nDetectorsLESS = DetectorIndexLESS$numLocalTraps, 
    habitatIDDet = DetectorIndexLESS$habitatGrid, 
    habitatGrid = myHabitat$IDCells.mx
  )
  
  ### ==== Define the initial values to be used in the model  ====
  
  # Initialise AC locations
  sxy <- MakeInitXY( y = y,
                     habitat.mx = myHabitat$habitat.mx,
                     detector.xy = scaled$data.scaled.xy,
                     IDCells.mx = myHabitat$IDCells.mx,
                     grid.xy = scaled$grid.scaled.xy)
  
  # Initialise z values
  z.init <- ifelse(!is.na(z), NA, 1)
  z.init[!is.na(z.init)] <- rbinom(sum(!is.na(z.init)), size = 1, prob = 0.5)
  
  nimInits <- list(p0 = runif(1,0,1),
                   sigma = runif(1,0,20), #runif(1,0,10),
                   z = z.init, #nimInits$z,           
                   psi = runif(1,0,1),
                   sxy = sxy #[ , ,1],
  )
  ### ==== Nimble model and the parameters to be monitored ====
  
  if(detfunfit == "HN"){
    nimParams <- c("N", "sxy", "z", "p0", "sigma")#Add NIn if you want to get the number of ids excluding buffer area
    modelCode <- modelCode_hn_LESSCachedAllSparse
  }
  
  if(detfunfit == "HNP"){
    nimInits$wd <- runif(1,0,10)
    nimParams <- c("N", "sxy", "z", "p0", "sigma", "wd")#Add NIn if you want to get the number of ids excluding buffer area
    modelCode <- modelCode_hnplateau_LESSCachedAllSparse
 }
  
  if(detfunfit == "EX"){
    nimParams <- c("N", "sxy", "z", "p0", "sigma")#Add NIn if you want to get the number of ids excluding buffer area
    modelCode <- modelCode_ex_LESSCachedAllSparse
  }
  
  
  ### ==== Save the NIMBLE input data (and simulation data) ====
  file.name <- paste( "NimbleInFile_sim", detfunsim,"_fit",detfunfit, ".RData", sep = "")
  this.path <- file.path(WD, modelName, file.name)
  
  save( nimData,
        nimConstants,
        nimParams,
        modelCode,
        nimInits,
        detected,
        simulated.ACS,
        myVars,
        params,
        detectionSimulationOut,
        detfunsim,
        detfunfit,
        file = this.path
  )
  
  # BUILDING AND COMPILING THE NIMBLE MODEL 
  model <- nimbleModel( code = modelCode,
                        constants = nimConstants,
                        data = nimData,
                        inits = nimInits,
                        check = F,       
                        calculate = F)
  print(model$calculate())
  
  cmodel <- compileNimble(model)
  
  print(cmodel$calculate())  
  
  #################################
  MCMCconf <- configureMCMC(model,
                            control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                            useConjugacy = FALSE)
  MCMCconf$addMonitors(nimParams)
  
  #################################
  # MCMCconf$printSamplers()
  
  if(detfunfit == "HNP"){
    # retrieve samplerConf list
    samplerConfList <- MCMCconf$getSamplers()
    samplerConfList[1:5]
    # change the 'adaptive' elements of the control list of the third sampler
    
    # sigma
    control <- samplerConfList[[3]]$control
    control$log <- F
    control$reflective <- T
    control$adaptive <- T
    control$scale <- 1
    samplerConfList[[3]]$setControl(control)
    
    # wd
    control <- samplerConfList[[4]]$control
    control$log <- F
    control$reflective <- T
    control$adaptive <- T
    control$scale <- 1
    samplerConfList[[4]]$setControl(control)
    
    # use this modified list of samplerConf objects in the MCMC configuration
    MCMCconf$setSamplers(samplerConfList)
    
    ##-- RECYCLING GIBBS SAMPLER
    # retrieve the current ordering of sampler execution
    ordering <- MCMCconf$getSamplerExecutionOrder()
    ordering
    len <- length(ordering)
    MCMCconf$setSamplerExecutionOrder(c(ordering[1:2], 
                                        rep(ordering[3], 10), # sigma
                                        rep(ordering[4], 10), # wd
                                        ordering[5:len]))
  }
  #################################
  
  MCMC <- buildMCMC(MCMCconf)
  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
  
  
  ### ====  RUN THE MCMC ====
  
  Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC, 
                                                    nburnin = burnin, 
                                                    niter = niter, 
                                                    nchains = nchains, #6, 
                                                    samplesAsCodaMCMC = TRUE))
  
  print(Runtime)
  
  ### ====  PLOT THE OUTPUT ====
  
  if(detfunfit == "HNP"){paramnames <- c('N', 'p0', 'sigma', 'wd')}
  if(detfunfit != "HNP"){paramnames <- c('N', 'p0', 'sigma')}
  sim.values <- c(myVars$POPULATION$N, params)
  

  par(mfrow=c(1,2))
  for(i in 1:length(paramnames))
  {
    traceplot(myNimbleOutput[ ,paramnames[i]])
    abline(h = sim.values[i], col = "black", lwd = 3, lty = 2)
    plot(density(unlist(myNimbleOutput[ ,paramnames[i]])), main = paramnames[i])
    abline(v = sim.values[i], col = "red", lwd = 2)
  }

  myNimbleOutput <- as.mcmc.list(myNimbleOutput) 
  if(length(myNimbleOutput) > 1){posteriorSamples <- do.call(rbind, myNimbleOutput)} 
  if(length(myNimbleOutput) == 1){posteriorSamples <- as.matrix(myNimbleOutput)} 
  
  summary.output <- do.call(rbind, lapply(paramnames, function(x){
    this.posterior <- posteriorSamples[, x]
    out <- list(MEAN = mean(this.posterior),
                CV = sd(this.posterior)/mean(this.posterior),
                LCI = quantile(this.posterior, 0.025),
                UCI = quantile(this.posterior, 0.975))
    return(unlist(out))
  }))
  dimnames(summary.output)[[1]] <- paramnames
  
  
  ### ==== Home range area calculation ====
  
  if(detfunfit == "HN"){
    sigchain <- posteriorSamples[,'sigma'] # vector (n.iter x 1)
    q<-qchisq(pr,2)#--FOLLOWING ROYLE BOOK PAGE 136
    radius<-sigchain*sqrt(q)
    area<-pi*radius^2
    HRest <- cbind(radius, area)
    dimnames(HRest)[[2]] <- c("HRradius", "HRarea")
  }
  
  if(detfunfit == "HNP"){
    HRAnim.arr <- HRA_HNP_nimble_array(xlim = c(0, habQualityDims[1]), ylim = c(0, habQualityDims[2]),
                                       ng = 1200, tol = 1E-5, niter = 2000)
    cHRAnim.arr <- compileNimble(HRAnim.arr, resetFunctions = TRUE)
    HRest <- cHRAnim.arr$run(p0Vec = posteriorSamples[,'p0'], sigmaVec = posteriorSamples[,'sigma'],
                             wdVec = posteriorSamples[,'wd'], pr = pr, d = 6)
    dimnames(HRest)[[2]] <- c("HRradius", "HRarea")
  }
  
  if(detfunfit == "EX"){
    HRAnim.arr <- HRA_EX_nimble_array(xlim = c(0, habQualityDims[1]), ylim = c(0, habQualityDims[2]),
                                      ng = 1200, tol = 1E-5, niter = 2000)
    cHRAnim.arr <- compileNimble(HRAnim.arr, resetFunctions = TRUE)
    HRest <- cHRAnim.arr$run(p0Vec = posteriorSamples[,'p0'], sigmaVec = posteriorSamples[,'sigma'],
                             pr = pr, d = 6)
    dimnames(HRest)[[2]] <- c("HRradius", "HRarea")
    
  }
  
  HR.true <- HRA_detfun(detfun = detfunsim, param = params, pr = pr,
                        xlim = c(0, habQualityDims[1]), ylim = c(0, habQualityDims[2]),
                        ng = 1200, tol = 1E-5, niter = 2000) 
  
  
  paramnames.hr <- c("HRradius", "HRarea")
  
  HR.mcmc.sample <- list()
  par(mfrow=c(1,2))
  for(i in 1:length(paramnames.hr)){
    this.posterior <- HRest[,paramnames.hr[i]]
    this.posterior.array <- matrix(this.posterior,
                                   nrow=length(this.posterior)/nchains,
                                   ncol=nchains,
                                   byrow = FALSE)
    HR.mcmc.sample[[i]] <- as.mcmc.list(lapply(1:nchains,function(x) mcmc(this.posterior.array[,x])))
    names(HR.mcmc.sample[[i]]) <- paste("chain", 1:nchains, sep = "")
    traceplot(HR.mcmc.sample[[i]])
    abline(h = HR.true[paramnames.hr[i]], col = "black", lwd = 3, lty = 2)
    plot(density(this.posterior), main = paramnames.hr[i])
    abline(v = HR.true[paramnames.hr[i]], col = "red", lwd = 2)
  }
  names(HR.mcmc.sample) <- paramnames.hr
  
  if(plot.save.check == TRUE){ graphics.off()}
  
  
  
  ### ==== Bayesian p-vaues ====
  
  ## ADD TO NIMDATA
  nimData$y <- y
  
  ## REMOVE SPARSE OBJECT TO AVOID WARNING MESSAGES ABOUT UNUSED DATA
  nimData <- nimData[!names(nimData) %in% c("yDets", "nbDetections")]
  nimConstants <- nimConstants[!names(nimConstants) %in% c("nMaxDetectors")]
  
  # nimData <- nimData[!names(nimData) %in% c("y", "yDets", "nbDetections", "trials",
  # "detectorIndex", "nDetectorsLESS", "habitatIDDet")]
  
  ## CREATE NIMBLE MODEL OBJECT
  if(detfunfit == "HN"){modelCode <- modelCode <- modelCode_hn_LESS_Cached}
  if(detfunfit == "HNP"){modelCode <- modelCode <- modelCode_hnplateau_LESS_Cached}
  if(detfunfit == "EX"){modelCode <- modelCode <- modelCode_ex_LESS_Cached}
  
  ### BUILDING AND COMPILING THE NIMBLE MODEL 
  
  model <- nimbleModel( code = modelCode, # UPDATED FROM INPUT FILES
                        constants = nimConstants,
                        data = nimData, # UPDATED FROM INPUT FILES
                        inits = nimInits,
                        check = F,       
                        calculate = F)
  
  ## COMPILE NIMBLE MODEL
  cmodel <- compileNimble(model)
  ## SET-UP THE FUNCTION 
  simFromPost <- simulateFromPosterior(model = model, nodes = dimnames(posteriorSamples)[[2]]) # Just using the model and the tracked nodes
  
  ## COMPILE THE FUNCITON (building c++ equivalent version of the code)
  csimFromPost <- compileNimble(simFromPost, project = cmodel, resetFunctions = TRUE)
  
  ## RUN THE COMPILED VERSION (much faster) (Note that, here we are using the already obtained MCMC samples of the tracked model parameters)
  csimFromPost$run(posteriorSamples = posteriorSamples) # Using the MCMC samples for all the chains
  
  
  ## STORING POSTERIOR REPLICATES OF THE DATA AND OTHER CHAINS
  
  ## STORE REPLICATED DATASETS
  zchain <- csimFromPost$mv[["z"]]  # list (niter) # components : vectors (M x 1)
  detprobchain<- csimFromPost$mv[["p"]] # list (niter) # components: matrix (M x ntrap)
  
  trials <- nimData$trials
  
  replicatedData <- csimFromPost$mv[["y"]] # list (niter) # components: matrix (M x ntrap) 
  
  n.iter <- dim(posteriorSamples)[1] #length(p0chain)
  names(replicatedData) = names(zchain) = names(detprobchain) = paste('iter', 1:n.iter, sep='')
  
  ### GETTING P-VALUES 
  pvalue_metric <- 'ALL'
  all_pvalue_metrics <- c('FT', 'PearsonsChiSq', 'FT_detector', 'FT_individual')
  if(sum(pvalue_metric %in% 'ALL') == 0){ 
    pvalue_metric <- intersect(all_pvalue_metrics, pvalue_metric)
  }
  pval <- pvaluefn(repData = replicatedData,
                   obsData = nimData$y,
                   detprob = detprobchain,
                   trials = trials,
                   metric = pvalue_metric)$pvalue
  
  ################################################################################################### 
  
  cat('No. of detected individuals = ',sum(detected), '\n', sep = '')
  
  cat('Posterior summary output:', '\n', sep = '')
  print(summary.output)
  
  cat("Simulated detection function: ", detfunsim, fill = TRUE)
  cat("Radius to achieve ", 100*pr, "% of home range area: ", HR.true["HRradius"], fill = TRUE, sep = "") # [SD]
  cat("Home range area: ", HR.true["HRarea"], fill = TRUE)
  
  cat('Bayesian p-values with different metrics:', '\n', sep = '')
  print(pval)
  
  ################################################################################################### 
  
  file.name <- paste( "NimbleOutFile_sim", detfunsim,"_fit",detfunfit, ".RData", sep = "")
  this.path <- file.path(WD, modelName, file.name)
  
  save(myNimbleOutput, summary.output,
       HR.mcmc.sample, HR.true, pr, pval,
       Runtime, niter, burnin, nchains,
       detfunsim, detfunfit,
       file=this.path)
  
  
}#runDetFun
################################################################################################### 
###################################################################################################
###################################################################################################
