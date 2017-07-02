###########################################################################################################
#
# Cite this code:                
# 
# Ernst, A. F. (2017). Probabilistic Time Series Clustering by Vector Autoregressive Metric - R code.
#      https://github.com/AnieBee/MasterThesis. GitHub.
#
#
########################################################################################################

library(mclust) # Mclust() 
library(MASS) #ginv()

### fit Gussian mixture model to data

ProbClustVAR1 <- function(X, Y, nObs_per_pers, nPersons, nClusters, nObs, nRandomStarts, SmartStart){
  nIndepVar <- ncol(X)
  nDepVar <- ncol(Y)  
  # stopifnot( SmartStart == 0 | SmartStart == 1 )
  
  # Allocation of output 
  allFit <- matrix(NA, (nRandomStarts + SmartStart), 1)
  Pmrvec <- array(NA, dim= c(nPersons, (nRandomStarts + SmartStart)))
  RegrCoeff <- array(NA, dim=c(nIndepVar+1, nClusters*nDepVar, (nRandomStarts + SmartStart))) 
  AttractionRate <- 0 
  
  # Indicate to which person each observation belongs
  index <- matrix(0, 1, nObs)
  Index <- nObs_per_pers[1:(length(nObs_per_pers)-1)]
  index[c(1, (apply(matrix(Index, length(Index), 1), 2, cumsum) +1))] = 1
  Part_obs <- matrix(apply(index, 1 , cumsum), 1, nObs)
  RegrModel_per_pers = matrix(0, nPersons, ((nIndepVar+1)*nDepVar))
  
  # Calculate VAR(1) model for each person
  for(pers in 1:nPersons){
    ind_obs = which(Part_obs == pers)
    X_cte = cbind(rep(1, length(ind_obs)), X[ind_obs, ]) 
    RegrCoeff_pers = ginv(t(X_cte)%*%X_cte)%*%t(X_cte)%*%Y[ind_obs, ]
    RegrModel_per_pers[pers, ] = as.vector(RegrCoeff_pers)
  }
  
  # Use Mclust() to cluster using random partitioning as start
  for(run in 1:nRandomStarts){
    mod <- Mclust(RegrModel_per_pers, G = nClusters, 
                   initialization = list(hcPairs = randomPairs(RegrModel_per_pers)))
    
    ProbRegrCoeff <- matrix(summary(mod, parameters = T)$mean , nrow = (numberVaribales + 1))
    
    allFit[run] = summary(mod, parameters = T)$bic
    Pmrvec[, run] = mod$classification
    RegrCoeff[, , run] =  ProbRegrCoeff
    
  }
  # Use Mclust() to cluster using model-based hierarchical agglomerative clustering results to initialize EM
  if(SmartStart){
    mod <- Mclust(RegrModel_per_pers, G = nClusters)
    
    # Make ProbRegrCoeff have dimensions nVar+1 by (nVAr*nCluster)
    ProbRegrCoeff <- matrix(summary(mod, parameters = T)$mean , nrow = (numberVaribales + 1))
    
    allFit[(nRandomStarts + SmartStart)] = mod$bic
    Pmrvec[ , (nRandomStarts + SmartStart)] = mod$classification 
    RegrCoeff[ , , (nRandomStarts + SmartStart)] = ProbRegrCoeff
  }

  Fit = max(na.omit(allFit))
  Index_RunWithBestFit = which.max(na.omit(allFit)) # which run has best fit
  AttractionRate = (length(which(allFit == Fit))/ length(allFit))
  
  ProbClustVAR1List <- list("Pmrvec" =   Pmrvec[ , Index_RunWithBestFit], 
                            "RegrCoeff" =  RegrCoeff[ , ,Index_RunWithBestFit],
                            "AttractionRate" = AttractionRate)
  return(ProbClustVAR1List)
}

