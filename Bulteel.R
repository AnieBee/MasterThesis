###########################################################################################################
#
# This code is an adaptation from Matlab code into R. The original code was written by Kirsten Bulteel.
# The original code is supplementary material of the following article:                
# 
# Bulteel, K., Tuerlinckx, F., Brose, A., & Ceulemans, E. (2016). Clustering vector autoregressive
#    models: capturing qualitative differences in within-person dynamics. Frontiers In Psychology,
#    7, 1540. doi:10.3389/fpsyg.2016.01540
#
###########################################################################################################


library(cluster) #agnes() agglomerative hierarchical clustering
library(MASS) #ginv()


### Subfunctions ### ------------------------

## Caclculate fit ##
CalculateFit <- function(X, Y, Pmr_vec, nObs_per_pers){
  
  nIndepVar <- ncol(X)
  nDepVar <- ncol(Y)
  nClusters <- max(Pmr_vec)
  nObs <- sum(nObs_per_pers)
  
  RegrCoeff <- array(0, dim=c(nIndepVar+1, nDepVar*nClusters))
  Y_pred <- array(0, dim=c(nObs, nDepVar))
  
  # Indicate cluster of each observation
  index <- array(0, dim=c(1, nObs))
  index[c(1,(cumsum(nObs_per_pers[1:(nrow(nObs_per_pers) - 1),]) + 1))] = 1
  Part_obs <- Pmr_vec[cumsum(index), ]
  
  # Calculate slope of each cluster and estimate Y
  for(cl in 1:nClusters){
    ind_obs = which(Part_obs == cl)
    X_cte = cbind(rep(1, length(ind_obs)), X[ind_obs,])
    RegrCoeff_cl = ginv(t(X_cte)%*%X_cte)%*%t(X_cte)%*%Y[ind_obs,]
    Y_pred[ind_obs,] = X_cte%*%RegrCoeff_cl
    RegrCoeff[,(nDepVar*(cl-1)+1):(cl*nDepVar)] = RegrCoeff_cl
  }
  NewFit = sum((Y-Y_pred)^2)
  
  CalculateFitList = list("NewFit" = NewFit,
                          "RegrCoeff" = RegrCoeff,
                          "Y_pred" = Y_pred)
  return(CalculateFitList)
}

## RANDPARTITION ##
# Generate random partition matrix for the one-mode clustering problem
randpartition <- function(nPersons, nClusters){

  IM = diag(nClusters)
  #Possible rows of the partition matrix (the number of rows equals the number of clusters)
  # Preallocate the partition matrix
  Pmr = array(0, dim=c(nPersons, nClusters))
  #As long as there are empty clusters, search process continues
  while (sum(colSums(Pmr) == 0) > 0){
    Pmr_vec = t(ceiling(nClusters*t(runif(nPersons))))
    Pmr = IM[Pmr_vec[,1],]
  }
  randpartitionList = list("Pmr" = Pmr,
                           "Pmr_vec" =  Pmr_vec)
  return(randpartitionList)
}

## CreateInitialSolution ##
#This function creates an initial partition matrix 
#and corresponding fit, based on a random assignment to clusters
CreateInitialSolution <- function(X, Y, nStim_per_pers, nClusters){
  
  nPersons = length(nStim_per_pers)
  randpartitionList = randpartition(nPersons, nClusters)
  Pmr = randpartitionList$Pmr
  Part_vec = randpartitionList$Pmr_vec
  #Distributes the persons randomly over the clusters
  CalculateFitList = CalculateFit(X, Y, Part_vec, nStim_per_pers)
  Fit = CalculateFitList$NewFit
  RegrCoeff = CalculateFitList$RegrCoeff
  Y_pred = CalculateFitList$Y_pred 
  
  CreateInitialSolutionList = list("Pmr" = Pmr,
                                   "NewFit" = Fit,
                                   "RegrCoeff" = RegrCoeff,
                                   "Y_pred" = Y_pred)
  return(CreateInitialSolutionList)
}


## upd_Pmr ###-----------------------------
#Update partition matrix
upd_Pmr <- function(X, Y, nObs_per_pers, Pmr){
  nPersons = nrow(nObs_per_pers)
  nClusters = ncol(Pmr)
  IM = diag(nClusters)
  
  #Possible rows of the partition matrix (the number of rows
  # equals the number of clusters)
  for(pers in 1:nPersons){
    fit_cl =  as.matrix(rep(0, nClusters))
    for(cl in 1:nClusters){
      # Put person in every cluster (1 per 1) and recalculate the fit
      Pmr[pers,] = IM[cl,]
      Part_vec_trial = as.matrix(rowSums(Pmr%*%diag(1:nClusters)))
      CalculateFitList <- CalculateFit(X, Y, Part_vec_trial, nObs_per_pers)
      NewFit_trial = CalculateFitList$NewFit
      fit_cl[cl, 1] = NewFit_trial
      
    }
    # Put person in cluster that gives best fit
    best = which(fit_cl == min(fit_cl))
    Pmr[pers, ] = IM[best[1], ] 
  }
  Pmr_vec = as.matrix(rowSums(Pmr%*%diag(1:nClusters)), ncol=1)
  
  CalculateFitList = CalculateFit(X, Y, Pmr_vec, nObs_per_pers)
  Fit = CalculateFitList$NewFit
  RegrCoeff = CalculateFitList$RegrCoeff
  Y_pred = CalculateFitList$Y_pred 
  
  upd_PmrList = list("NewFit" = Fit,
                     "RegrCoeff" = RegrCoeff,
                     "Pmr" = Pmr,
                     "Y_pred" = Y_pred)
  return(upd_PmrList)
}


# Main function --------------------------------

# INPUT #
# X <- is 2 dim to allow for unequal # of obs per person
#       rows= timepoints*#person; cols = variables;
# Y <- outcome variable: same as X but from time = 2 to time = T+1


# nObs_per_pers =   each row contains # of observations for participant of that number
# (length(nObs_per_pers)= #participants)
# needs to be vector of dimension: npers x 1
# example: nObs_per_pers = matrix(vector, length(vector), 1)

# nObs <- sum(nObs_per_pers)
# nClusters <- # of clusters
# nRandomStarts <- number of random starts
# SmartStart: numerical logical:
# 1 = use of hierarchical clustering as smart start
# 0 = no use of a smart start


# OUTPUT #
#Pmr: nPersons * nClusters partition matrix 
#Pmrvex: pmr in vector format
#RegrCoeff: matrices for each cluster containing cluster-specific intercepts and slopes
#Fit: SS of prediction errors
#varDr: total SS
#Y_pred: predicted values for outcome variables

FitClVAR1model <- function(X, Y, nObs_per_pers, nPersons, nClusters, nObs, nRandomStarts, SmartStart){
  nIndepVar <- ncol(X)
  nDepVar <- ncol(Y)
  
  # Preallocate the output matrices
  allFit <- matrix(0, (nRandomStarts + SmartStart), 1)
  Pmr <- array(0, dim= c(nPersons, nClusters, (nRandomStarts + SmartStart)))
  RegrCoeff <- array(0, dim=c(nIndepVar + 1, nClusters*nDepVar, (nRandomStarts + SmartStart))) 
  GOF <- array(0, dim=c((nRandomStarts + SmartStart),1))
  Y_pred <- array(0, dim=c(nObs,nDepVar, (nRandomStarts + SmartStart)))
  AttractionRate <- 0 
  
  # Calculate the total SS ----
  varDr = sum((Y-mean(mean(Y)))^2) 
  
  # If the number of clusters equals one, an ordinary least squares instead
  # of an alternating least squares (ALS) algorithm is sufficient.
  if(nClusters == 1){
    Pmr <- matrix(1, nPersons, 1)
    Pmrvec <- Pmr
    CalculateFitList  <- CalculateFit(X, Y, Pmrvec, nObs_per_pers)
    Fit = CalculateFitList$NewFit
    RegrCoeff = CalculateFitList$RegrCoeff
    Y_pred = CalculateFitList$Y_pred

    GOF = 1 - Fit/(varDr)
  }else{ #Perform the ALS algorithm
    #Random starts
    for(run in 1:nRandomStarts){
      # cat ("\n ", nClusters,"  ", run," ")  
      CreateInitialSolutionList <- CreateInitialSolution(X, Y, nObs_per_pers, nClusters)
      Pmr_run = CreateInitialSolutionList$Pmr
      NewFit = CreateInitialSolutionList$NewFit 
      diff = 1;
      while(diff > 0.000001){ #Keep updating as long as CurFit really improves
        CurFitTotal = NewFit
        upd_PmrList <- upd_Pmr(X, Y, nObs_per_pers, Pmr_run)
        NewFit = upd_PmrList$NewFit
        RegrCoeff_run = upd_PmrList$RegrCoeff
        Pmr_run = upd_PmrList$Pmr
        Y_pred_run = upd_PmrList$Y_pred
        diff = CurFitTotal-NewFit
      }
      GOF_run= 1 - NewFit/(varDr)
      allFit[run] = NewFit
      Pmr[, , run] = Pmr_run
      RegrCoeff[, , run] = RegrCoeff_run
      GOF[run] = GOF_run
      Y_pred[ , , run] = Y_pred_run
    }
    
    if(SmartStart == 1){ # Using a rational start
      
      # Indicate to which person each observation belongs
      index <- matrix(0, 1, nObs)
      Index <- nObs_per_pers[1:(length(nObs_per_pers)-1)]
      index[c(1, (apply(matrix(Index, length(Index), 1), 2, cumsum) +1))] = 1
      Part_obs <- matrix(apply(index, 1 , cumsum), 1, nObs)
      
      # Compute the parameters of a VAR(1) model per person
      RegrModel_per_pers = matrix(0, nPersons, nIndepVar*nDepVar)
      for(pers in 1:nPersons){
        # cat(pers, nPersons, "\n")
        ind_obs = which(Part_obs == pers)
        X_cte = cbind(rep(1, length(ind_obs)), X[ind_obs,]) 
        RegrCoeff_pers = ginv(t(X_cte)%*%X_cte)%*%t(X_cte)%*%Y[ind_obs, ]
        RegrCoeff_pers = RegrCoeff_pers[-1,]
        RegrModel_per_pers[pers, ] = as.vector(RegrCoeff_pers)
      }
      
      #Hierarchical clustering on the individual VAR(1) weights
      Storage = agnes(RegrModel_per_pers, metric = "euclidean", method = "ward") 
      Part_HierarchicalClustering = cutree(Storage, nClusters)
      
      # Run the ALS algorithm starting from the solution of the
      # hierarchical clustering procedure
      CalculateFitList = CalculateFit(X, Y, matrix(Part_HierarchicalClustering, ncol= 1), nObs_per_pers)
      Fit_initial = CalculateFitList$NewFit
      IM = diag(nClusters)
      Pmr_run = IM[Part_HierarchicalClustering, ]
      CurFitTotal = Fit_initial
      diff = 1
      
      # Keep updating as long as CurFit really improves
      while (diff > 0.000001){
        upd_PmrList = upd_Pmr(X, Y, nObs_per_pers, Pmr_run)
        NewFit = upd_PmrList$NewFit
        RegrCoeff_run = upd_PmrList$RegrCoeff
        Pmr_run = upd_PmrList$Pmr
        Y_pred_run = upd_PmrList$Y_pred
        diff = CurFitTotal - NewFit
        CurFitTotal = NewFit
      }

      GOF_run = 1 - NewFit/(varDr)
      allFit[length(allFit)] = NewFit
      Pmr[ , , (nRandomStarts + SmartStart)] = Pmr_run 
      RegrCoeff[ , , length(allFit)] = RegrCoeff_run
      GOF[(nRandomStarts + SmartStart)] = GOF_run
      Y_pred[ , , (nRandomStarts + SmartStart)] =  Y_pred_run
      
    }
    
    # Select the model with the best fit
    Fit = min(allFit)
    Index_RunWithBestFit = which.min(allFit) # here it is determined which run has lowest fit
    AttractionRate = (length(which(allFit == Fit))/ length(allFit))
    Pmr = Pmr[ , , Index_RunWithBestFit]
    RegrCoeff = RegrCoeff[ , ,Index_RunWithBestFit]
    GOF = GOF[Index_RunWithBestFit]
    Y_pred = Y_pred[ , ,Index_RunWithBestFit]
    Pmrvec = as.matrix(rowSums(Pmr%*%diag(1:nClusters)), ncol=1)
    
  }
  FitClVAR1modelList <- list("Pmr" = Pmr, "Pmrvec" = Pmrvec, "RegrCoeff" = RegrCoeff, "Fit" = Fit,
                             "GOF" = GOF, "varDr"= varDr, "Y_pred"= Y_pred, "AttractionRate" = AttractionRate)
  return(FitClVAR1modelList)
}



