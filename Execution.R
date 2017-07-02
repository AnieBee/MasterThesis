###########################################################################################################
#
# Cite this code:                
# 
# Ernst, A. F. (2017). Probabilistic Time Series Clustering by Vector Autoregressive Metric - R code.
#      https://github.com/AnieBee/MasterThesis. GitHub.
#
#
########################################################################################################

library(mclust) # adjustedRandIndex() 
library(e1071) # matchClasses() 

nRandomStarts <-100 # of random starts: 1 hierchical and 100 random starts
SmartStart <- 1 #numerical logical 1 = use of hierarchical clustering as smart start
numberVaribales <- 4
BoolProbClustDataUsed <- 1 # Boolean that runs from 0 to 1 in outer loop, 0 = BulteelDataGenerationData is analyzed; 1 = ProbClustDataGeneration Data is analyzed

### Vary ### 
nPersons <- c(30, 60, 120) # index: i 
nClusters <- c(2, 4) # index: j 
ClusterSize <- c(0, .1, .6) #index: k
Observations <- c(50, 100, 500)  #index: l. 1 Observation gets added in varima.sim because in X you delete last Observation
#and in Y you need to delete first Observation
Distance <- c(1, 2, 3) # index: m
# 1 = small diff with all positive cross-regressive coefs
# 2 =  larg diffs between clusters with all positive cross-regressive coefs
# 3 = largere diffs between clusters with pos and neg cross-regressive coefs
Replications <- 10 #index: r
numbPhiOffDiag <- ((numberVaribales^2)-numberVaribales)/2 
VariancePhi_kqs <- .025


### Execution and Analysis -----------
source("Bulteel.R") # Analyses Data using the Bulteel non-probabilistic clustering model
source("ProbClust.R") # Analyses Data using the probabilistic GMM model
for(a in 0:BoolProbClustDataUsed){ 
  
  if(a){#ProbClust Data is analyzed
    source("ProbClustDataGeneration.R")
  }else{#Bulteel Data is analyzed
    source("BulteelDataGeneration.R")
  }
  ### Sotrage of output ###
  NonProb_Prob <- 2 # Used as index. First dimension: Non-probalistic (Bulteel) clustering Model was used, 2nd dimension: Prob Clust Model employed
  # Adjusted Rand Indeces
  ARI <- array(NA, dim =c(NonProb_Prob, length(nPersons), length(Observations),  length(ClusterSize), length(nClusters), length(Distance), Replications))
  # ARI[ , i, l, k, j, m, r]
  # In this case this is Euclidean distance
  MahalanobisDistance <- array(NA, dim = c(NonProb_Prob, length(nPersons), length(Observations),  
                                           length(ClusterSize), length(nClusters), length(Distance), Replications))
  # MahalanobisDistance[ , i, l, k, j, m, r]
  AttractionRates <- array(NA, dim = c(NonProb_Prob, length(nPersons), length(Observations),  
                                       length(ClusterSize), length(nClusters), length(Distance), Replications))
  # Final estimates of VAR(1) slopes 
  SlopeEstimates <- list(length(nClusters))  # list 1 = 2 clusters , list 2 = 4 clusters #index: j
  for(runner in 1:length(nClusters)){
    SlopeEstimates[[runner]] <- array(NA, dim = c(numberVaribales, (nClusters[runner]*numberVaribales), 
                                                  NonProb_Prob,  # First dimension: Non-probalistic (Bulteel) clustering Model was used, 2nd dimension: Prob Clust Model employed
                                                  length(nPersons), # index:i
                                                  length(ClusterSize), # index: k
                                                  length(Observations), #index: l
                                                  length(Distance), # index: m
                                                  Replications #index r
                                                  #SlopeEstimates[[j]][  , , NonProb_Prob, i, k, l, m, r]
    ))
  }
  # Final estimates of cluster membership
  ClusterMembershipEstimates <-list(length(nPersons)) # list 1 = 30 people , list 2 = 60 people....#index: i
  for(runner in 1:length(nPersons)){
    ClusterMembershipEstimates[[runner]] <- array(NA, dim = c(nPersons[runner], 
                                                              NonProb_Prob,  # First dimension: Non-probalistic (Bulteel) clustering Model was used, 2nd dimension: Prob Clust Model employed
                                                              length(nClusters), # index:j
                                                              length(ClusterSize), # index: k
                                                              length(Observations), #index: l
                                                              length(Distance), # index: m
                                                              Replications #index r
                                                              #ClusterMembershipEstimates[[i]][  , NonProb_Prob, j, k, l, m, r]
    ))
  }
  
  for(i in 1:length(nPersons)){
    for(l in 1:length(Observations)){
      
      nObs_per_pers = matrix(c(rep(Observations[l], nPersons[i])), ncol= 1) #index[i, 1]
      nObs <- sum(nObs_per_pers)
      
      for(k in 1:length(ClusterSize)){
        for(j in 1:length(nClusters)){
          cat("\n", "Pers", nPersons[i], "Obs", Observations[l], "nClust", nClusters[j], "percentage", ClusterSize[k], "\n")
          for(m in 1:length(Distance)){ # Distance between clusters
            for(r in 1:Replications){
              # delete last obs of every participant in X and first obs in Y
              X <- DataSets[[ILcounter[i, l]]][-c(seq((Observations[l] + 1), ((Observations[l] + 1)*nPersons[i]), by = (Observations[l] + 1))) , 1:numberVaribales, j, k, m, r]
              Y <- DataSets[[ILcounter[i, l]]][ -c(seq(1, ((Observations[l] + 1)*nPersons[i]), by = (Observations[l] + 1))), 1:numberVaribales, j, k, m, r]
              TruePartition <- DataSets[[ILcounter[i, l]]][c(seq((Observations[l] + 1), ((Observations[l] + 1)*nPersons[i]), by = (Observations[l] + 1))) , numberVaribales + 1, j, k, m, r] 
              
              
              ### Non-Prob Clustering Model Employed for Analysis ### --------------------------
              FitClVAR1modelList <- FitClVAR1model(X, Y, nObs_per_pers, nPersons[i], nClusters[j], nObs, nRandomStarts, SmartStart)
              # determine the permutation of labels that yields the highest proportion of agreement between labels
              SwapTruePartition <-  c(as.vector(matchClasses(table(TruePartition, FitClVAR1modelList$Pmrvec), method = "exact")))
              
              ### Save Estimated Model ###
              ClusterMembershipEstimates[[i]][  , 1, j, k, l, m, r] <- FitClVAR1modelList$Pmrvec
              SlopeEstimates[[j]][  , , 1, i, k, l, m, r] <- FitClVAR1modelList$RegrCoeff[-1, ]
              
              ### ARI ###
              ARI[ 1, i, l, k, j, m, r] <- adjustedRandIndex(TruePartition, FitClVAR1modelList$Pmrvec)
              
              ### Mahalanobis distance ###
              # distance from every clusters estimated parameters to the "true GMM" parameters 
              MahalanobisSum <- 0 
              for(n in 1:nClusters[j]){
                Increment <- dist(rbind(matrix(SlopeEstimates[[j]][  , , 1, i, k, l, m, r], 
                                               nrow = (numberVaribales*(numberVaribales)))[, SwapTruePartition[n]], 
                                        as.vector(t(Phi_ks[ , , i, j, k, l, m, n, r]))), method = "euclidean")
                MahalanobisSum = MahalanobisSum + Increment
              }
              MahalanobisDistance[ 1, i, l, k, j, m, r] <- MahalanobisSum
              
              ### Attraction Rate ### 
              AttractionRates[ 1, i, l, k, j, m, r] <- FitClVAR1modelList$AttractionRate
              
              ### ProbClust Model Employed for Analysis ###------------------------------------
              ProbClustMod <- ProbClustVAR1(X, Y, nObs_per_pers, nPersons[i], nClusters[j], nObs, nRandomStarts, SmartStart)
              # determine the permutation of labels that yields the highest proportion of agreement between labels
              SwapTruePartition <-  c(as.vector(matchClasses(table(TruePartition, ProbClustMod$Pmrvec), method = "exact")))
              
              
              ### Save Estimated Model ###
              ClusterMembershipEstimates[[i]][  , 2, j, k, l, m, r] <- ProbClustMod$Pmrvec
              SlopeEstimates[[j]][  , , 2, i, k, l, m, r] <- ProbClustMod$RegrCoeff[-1, ]
              
              #ARI
              ARI[ 2, i, l, k, j, m, r] <- adjustedRandIndex(TruePartition, ProbClustMod$Pmrvec)
              
              ### Mahalanobis distance ###
              # distance from every clusters estimated parameters to the true parameters 
              MahalanobisSum <- 0 
              for(n in 1:nClusters[j]){
                Increment <- dist(rbind(matrix(SlopeEstimates[[j]][  , , 2, i, k, l, m, r], 
                                                 nrow = (numberVaribales*(numberVaribales)))[, SwapTruePartition[n]], 
                                          as.vector(t(Phi_ks[ , , i, j, k, l, m, n, r]))), method = "euclidean")
                
                MahalanobisSum = MahalanobisSum + Increment
                MahalanobisSum = MahalanobisSum + Increment
              }
              MahalanobisDistance[ 2, i, l, k, j, m, r] <- MahalanobisSum
              
              ### Attraction Rate ### 
              AttractionRates[ 2, i, l, k, j, m, r] <- ProbClustMod$AttractionRate
            }
          }
        }
      }
    }
  }
  # Divide the Euclidean Distance by the number of clusters to get mean values
  MahalanobisDistance[ , , , , 1, , ] <- MahalanobisDistance[ , , , , 1, , ]*.5
  MahalanobisDistance[ , , , , 2, , ] <- MahalanobisDistance[ , , , , 2, , ]*.25
  
  ### Save Results ###
  if(a){# ProbClust Data was analyzed
    ProbClustDataResults <- list("ARI" = ARI, "AttractionRates" = AttractionRates, 
                                 "MahalanobisDistance" = MahalanobisDistance,
                                 "ClusterMembershipEstimates" = ClusterMembershipEstimates, 
                                 "SlopeEstimates" = SlopeEstimates)
    save(ProbClustDataResults, file= "data/ProbClustDataResults.RData")
  }else{# Bulteel Data was analyzed
    BulteelDataResults <- list("ARI" = ARI, "AttractionRates" = AttractionRates, 
                               "MahalanobisDistance" = MahalanobisDistance,
                               "ClusterMembershipEstimates" = ClusterMembershipEstimates, 
                               "SlopeEstimates" = SlopeEstimates)
    save(BulteelDataResults, file= "data/BulteelDataResults.RData")
  }
}
