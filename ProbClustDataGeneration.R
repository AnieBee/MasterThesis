###########################################################################################################
#
# Cite this code:                
# 
# Ernst, A. F. (2017). Probabilistic Time Series Clustering by Vector Autoregressive Metric - R code.
#      https://github.com/AnieBee/MasterThesis. GitHub.
#
#
########################################################################################################

library(portes) #varima.sim()
library(MASS) #mvrnorm()
set.seed(705)


nPhikqGenerated <- 0
nPhikqUsed <- 0 

Phi_ks <- array(1, dim = c(numberVaribales, numberVaribales, 
                           length(nPersons), # index:i
                           length(nClusters), # index:j
                           length(ClusterSize), # index: k
                           length(Observations), #index: l
                           length(Distance), # index: m
                           max(nClusters), #index n
                           Replications #index r
                    ))



dims = (Observations[1] + 1)*nPersons[1]
ILcounter <- matrix(1:(length(nPersons)*length(Observations)), 
                    nrow= length(nPersons), ncol= length(Observations), byrow= T)
# ILcounter[i, l] gives list index for DataSets[[index]] to return i = nPersons, l = Observations 
# Allows array indexing properties on the Dataset despite it being a list of arrays of unequal dimensions

for(i in 1:length(nPersons)){
  for(l in 1:length(Observations)){
    if(! (i == 1 & l == 1)){
      dims2 = (Observations[l] + 1)*nPersons[i]
      dims = c(dims, dims2)
    }
  }
}

DataSets <- list(1:length(dims))

for(i in 1:length(dims)){
  DataSets[[i]] <- array(i, dim = c(dims[i], 
                                    numberVaribales + 1 , # numberVariables + 1 col for true cluster membership
                                    # length(nPersons), # index:i
                                    length(nClusters), # index:j
                                    length(ClusterSize), # index: k
                                    # length(Observations), #index: l
                                    length(Distance), # index: m
                                    Replications #index r
  ))
}


#Main Function ----
for(i in 1:length(nPersons)){
  for(j in 1:length(nClusters)){
    for(k in 1:length(ClusterSize)){ #Determine # of persons in cluster = nPerClust
      nPerClust <- matrix((1-ClusterSize[k])/(nClusters[j]-as.numeric(ClusterSize[k]!=0)), nrow=nClusters[j], ncol=1)
      #nPerClust is vector of length nClusters[j]
      #contains the remaining percentages distributed equally
      if(ClusterSize[k]){ # if first cluster has different percentage, replace equal precentage with it
        nPerClust[1] = ClusterSize[k]
      }
      nPerClust = nPerClust*nPersons[i]
      if(nPerClust[1] == 7.5){
        nPerClust[sample(1:4, 2, replace = FALSE)] = 8
        nPerClust = floor(nPerClust)
      }
      cat("Pers", nPersons[i], "nClust", nClusters[j], "percentage", ClusterSize[k], nPerClust,"\n")
      
      for(l in 1:length(Observations)){ # Observations
        for(m in 1:length(Distance)){ # Distance between clusters
          for(r in 1:Replications){
            
            # SunnySideUp contains data for all individuals of a cluster after one another.
            # To have a non-ordered dataset where individuals of the same cluster are
            # not presented after one another, persons are scrambled unorderly. The result is 
            # ScrambledData which is the final data set for the condition that is saved in DataSets.
            SunnySideUpData <- array(NA,
                                     dim= c((Observations[l] + 1),  numberVaribales + 1, nPersons[i]
                                            # numberVaribales  +1 so you can save which cluster an indiviual came from
                                            # Observations[l] + 1 because first Y and last X is not used
                                     ))
            
            ScrambledData <- array(NA, dim= c((Observations[l] + 1)*nPersons[i],
                                              numberVaribales + 1 # +1 so you can save which cluster an individual came from
                                    )) 
            
            PeopleAlreadyInThisDataSet <- 0
            for(n in 1:c(nClusters[j])){# For every cluster determine Phi_k--------
                cat(i, j , k, l, m, r, n , nClusters[j], Distance[m], Observations[l], "\n")
              
                Phi_k <- matrix(NA, ncol = numberVaribales, nrow = numberVaribales)
                diag(Phi_k) <- runif(numberVaribales, min= .7, max= .9)
              
                if(Distance[m] == 2){# similiar cluster condition
                    Phi_k[upper.tri(Phi_k)] <- sample(c(runif(numbPhiOffDiag, min= .3, max= .5), 
                                                    runif(numbPhiOffDiag, min= 0, max= .2)), 
                                                  numbPhiOffDiag)
                    Phi_k[lower.tri(Phi_k)] <- sample(c(runif(numbPhiOffDiag, min= .3, max= .5), 
                                                    runif(numbPhiOffDiag, min= 0, max= .2)), 
                                                  numbPhiOffDiag)
                    # rescaling to ensure stationarity
                    Phi_k = Phi_k*(.99/max(abs(eigen(Phi_k)$values)))
                }else{ # highley similiar or disimilar  
                    Phi_k[upper.tri(Phi_k)] <- runif(numbPhiOffDiag, min= .3, max= .5)
                    Phi_k[lower.tri(Phi_k)] <- runif(numbPhiOffDiag, min= .3, max= .5)
                    # rescaling to ensure stationarity
                    Phi_k = Phi_k*(.99/max(abs(eigen(Phi_k)$values)))
                
                    if(Distance[m] == 3){ # highley disimilar: half of cross regressive coeffs are assigned negative signs
                        indexNegNumbs <- c(sample(x = c(which(upper.tri(Phi_k)), which(lower.tri(Phi_k))), size = numbPhiOffDiag, replace = FALSE))
                        Phi_k[indexNegNumbs] <- Phi_k[indexNegNumbs]*(-1)
                  
                      }
                } #------
                Phi_ks[ , , i, j, k, l, m, n, r] <- Phi_k # save Phi_k
                Sigma <- diag(numberVaribales)
                #Sigma is the variance/covariance of the white noise series: set to identity
               
                
                for(q in 1:nPerClust[n]){# n gives the number of current evaluated cluster
                    cat(q, nPerClust[n], "\n") # q starts at 1 for new cluster, q is used as counter for people in the current cluster
                    Phi_kq <-  array(1, dim = c(numberVaribales, numberVaribales, 1)) 
                    while( sum(abs(Phi_kq) >= 1) || max(Mod(eigen(Phi_kq[ , ,1], symmetric = FALSE, only.values = TRUE)$value)) >= 1 ){  # create variation within clusters, sample participants individual Phis from MultiNormal of Phi_k 
                      # Don't leave loop until all elements of Phi_kq have an absolute value smaller than 1 
                      # and until the largest Modulus of eigenvalues of Phi_kq is smaller than 1
                      Phi_kq <- array(mvrnorm(n = 1, mu = as.vector(Phi_k), 
                                              Sigma = (diag(numberVaribales*numberVaribales)*VariancePhi_kqs),
                                              empirical = FALSE),
                                      dim = c(numberVaribales, numberVaribales, 1))
                      nPhikqGenerated = nPhikqGenerated + 1
                    }
                    nPhikqUsed =  nPhikqUsed + 1
                    
                    # the assertions below are always met because of the while loop above
                    # stopifnot(max(Mod(eigen(Phi_kq[ , ,1], symmetric = FALSE, only.values = TRUE)$value )) < 1 )
                    # stopifnot( sum(abs(Phi_kq)>= 1) == 0)
                    
                    Data <- varima.sim(list(ar = Phi_kq, ma = NULL), n = (Observations[l] +1), k = numberVaribales, 
                                       sigma = Sigma)
                    # Data is array[rows hold numberVariables, col hold Observations]
                    Data <- cbind(Data, rep(n,  (Observations[l] +1)))
                    SunnySideUpData[ , , PeopleAlreadyInThisDataSet + q] <- Data
                }
                PeopleAlreadyInThisDataSet <- PeopleAlreadyInThisDataSet + q
            }
            # Scramble Sunny Side up and save in DataSets so people appear unordered in the final dataset
            IndexDataScramble <- c(sample(1:nPersons[i], nPersons[i], replace = FALSE))
            ScrambledData <- apply(SunnySideUpData[ , , IndexDataScramble], 2, cbind)
            DataSets[[ILcounter[i, l]]][ , , j, k, m, r] <- ScrambledData 
            cat(i, j, k, l, m, n , r, "\n")
          }
        }
      }
    }
  }
}

save(DataSets, file= "data/ProbClustDataGenerationDataSets.RData")
save(Phi_ks, file= "data/ProbClustDataGenerationPhi_ks.RData")
PhikqProportion <- list("nPhikqGenerated" = nPhikqGenerated, "nPhikqUsed" = nPhikqUsed)
save(PhikqProportion, file= "data/PhikqProportion.RData")
