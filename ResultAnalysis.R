library(reshape)
library(heplots)
library(Hmisc)
options(scipen = 1)
options(digits = 3)
nRandomStarts <- 100 # of random starts: 1 hierchical and 100 random starts
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
# MahalanobisCovariance <- set to I*VariancePhi_kqs when ProbClust Data is analyzed (Mahalanobis distance)
# set to I when Bulteel Data is analyzed, (Mahalanobis distance reduceds to Euclidean distance)
NonProb_Prob <- 2 # First dimension: Non-probalistic (Bulteel) clustering Model was used, 2nd dimension: Prob Clust Model employed

#ARI[ NonProbProb, i, l, k, j, m, r]
### First Steps
length(which(is.na(ProbClustDataResults$ARI)))
length(which(is.na(BulteelDataResults$ARI)))

## Data Cleaning
BoolProbClustDataUsed <- # Assign 0 or 1 

if(BoolProbClustDataUsed){#ProbClust Data is used
  load("data/ProbClustDataResults.RData")
  ARI <- ProbClustDataResults$ARI
  MahalanobisMean <- ProbClustDataResults$MahalanobisDistance
  AttractionRates <- ProbClustDataResults$AttractionRates
}else{#Bulteel Data is used
  load("data/BulteelDataResults.RData")
  ARI <- BulteelDataResults$ARI
  MahalanobisMean <- BulteelDataResults$MahalanobisDistance
  AttractionRates <- BulteelDataResults$AttractionRates
}
# Divide the Overall Mahalanobis Distances by the number of clusters to get the Mean value

## ANOVA ARI ##-----------------
dimnames(ARI)[[1]] <- c("NonProbModel", "ProbModel")
dimnames(ARI)[[2]] <- c("30Pers", "60Pers", "120Pers") #nPersons
dimnames(ARI)[[3]] <- c("50Obs", "100Obs", "500Obs") #Observations
dimnames(ARI)[[4]] <- c("EqualProportion", "MinorityCluster", "MajorityCluster") #Clustersize
dimnames(ARI)[[5]] <- c("2CL", "4CL") #nClusters
dimnames(ARI)[[6]] <- c("SmallDistance", "MediumDistance", "LargeDistance") #Distance

MeltARI <- melt.array(ARI)
names(MeltARI) <- c("ClustModel", "nPers", "nObs", "ClusterSize", "nClusters", "Distance", "Rep", "value" )
str(MeltARI)

ANOVA <- aov(value ~ ClustModel + nPers + nObs + ClusterSize + nClusters + Distance +
               ClustModel*nPers + ClustModel*nObs + ClustModel*ClusterSize + ClustModel*nClusters + ClustModel*Distance +
               nPers*nObs + nPers*ClusterSize + nPers*nClusters + nPers*Distance +
               nObs*ClusterSize + nObs*nClusters + nObs*Distance +
               ClusterSize*nClusters + ClusterSize*Distance +
               nClusters*Distance
             , data = MeltARI)

summary(ANOVA)
anova(ANOVA)
Anova(ANOVA)
model.tables(ANOVA, type = "effects")
model.tables(ANOVA, type = "means")
etasq(ANOVA)
which(etasq(ANOVA)>.1)


## ANOVA ARI[1, ..] ##-----------------
ARI <- BulteelDataResults$ARI[2, , , , , , ]
dimnames(ARI)[[1]] <- c("30Pers", "60Pers", "120Pers") #nPersons
dimnames(ARI)[[2]] <- c("50Obs", "100Obs", "500Obs") #Observations
dimnames(ARI)[[3]] <- c("EqualProportion", "MinorityCluster", "MajorityCluster") #Clustersize
dimnames(ARI)[[4]] <- c("2CL", "4CL") #nClusters
dimnames(ARI)[[5]] <- c("SmallDistance", "MediumDistance", "LargeDistance") #Distance

MeltARI <- melt.array(ARI)
names(MeltARI) <- c("nPers", "nObs", "ClusterSize", "nClusters", "Distance", "Rep", "value" )
str(MeltARI)

ANOVA <- aov(value ~ nPers + nObs + ClusterSize + nClusters + Distance +
               nPers*nObs + nPers*ClusterSize + nPers*nClusters + nPers*Distance +
               nObs*ClusterSize + nObs*nClusters + nObs*Distance +
               ClusterSize*nClusters + ClusterSize*Distance +
               nClusters*Distance
             , data = MeltARI)

summary(ANOVA)
anova(ANOVA)
Anova(ANOVA)
model.tables(ANOVA, type = "effects")
model.tables(ANOVA, type = "means")
etasq(ANOVA)
which(etasq(ANOVA)>.1)





## Mahalanobis ARI ##-----------------
dimnames(MahalanobisMean)[[1]] <- c("NonProbModel", "ProbModel")
dimnames(MahalanobisMean)[[2]] <- c("30Pers", "60Pers", "120Pers") #nPersons
dimnames(MahalanobisMean)[[3]] <- c("50Obs", "100Obs", "500Obs") #Observations
dimnames(MahalanobisMean)[[4]] <- c("EqualProportion", "MinorityCluster", "MajorityCluster") #Clustersize
dimnames(MahalanobisMean)[[5]] <- c("2CL", "4CL") #nClusters
dimnames(MahalanobisMean)[[6]] <- c("SmallDistance", "MediumDistance", "LargeDistance")#Distance

MeltMahalanobis <- melt.array(MahalanobisMean)
names(MeltMahalanobis) <- c("ClustModel", "nPers", "nObs", "ClusterSize", "nClusters", "Distance", "Rep", "value" )
str(MeltMahalanobis)

ANOVA <- aov(value ~ ClustModel + nPers + nObs + ClusterSize + nClusters + Distance +
               ClustModel*nPers + ClustModel*nObs + ClustModel*ClusterSize + ClustModel*nClusters + ClustModel*Distance +
               nPers*nObs + nPers*ClusterSize + nPers*nClusters + nPers*Distance +
               nObs*ClusterSize + nObs*nClusters + nObs*Distance +
               ClusterSize*nClusters + ClusterSize*Distance +
               nClusters*Distance
             , data = MeltMahalanobis)
etasq(ANOVA)
summary(ANOVA)
anova(ANOVA)
Anova(ANOVA)
model.tables(ANOVA, type = "effects")
model.tables(ANOVA, type = "means")
which(etasq(ANOVA)>.1) #partial eta-squared



## ANOVA Attraction rates ---------------
dimnames(AttractionRates)[[1]] <- c("NonProbModel", "ProbModel")
dimnames(AttractionRates)[[2]] <- c("30Pers", "60Pers", "120Pers") #nPersons
dimnames(AttractionRates)[[3]] <- c("50Obs", "100Obs", "500Obs") #Observations
dimnames(AttractionRates)[[4]] <- c("EqualProportion", "MinorityCluster", "MajorityCluster") #Clustersize
dimnames(AttractionRates)[[5]] <- c("2CL", "4CL") #nClusters
dimnames(AttractionRates)[[6]] <- c("SmallDistance", "MediumDistance", "LargeDistance") #Distance

MeltAttractionRates <- melt.array(AttractionRates)
names(MeltAttractionRates) <- c("ClustModel", "nPers", "nObs", "ClusterSize", "nClusters", "Distance", "Rep", "value" )
str(MeltAttractionRates)

ANOVA <- aov(value ~ ClustModel + nPers + nObs + ClusterSize + nClusters + Distance +
               ClustModel*nPers + ClustModel*nObs + ClustModel*ClusterSize + ClustModel*nClusters + ClustModel*Distance +
               nPers*nObs + nPers*ClusterSize + nPers*nClusters + nPers*Distance +
               nObs*ClusterSize + nObs*nClusters + nObs*Distance +
               ClusterSize*nClusters + ClusterSize*Distance +
               nClusters*Distance
             , data = MeltAttractionRates)

summary(ANOVA)
anova(ANOVA)
Anova(ANOVA)
model.tables(ANOVA, type = "effects")
model.tables(ANOVA, type = "means")
etasq(ANOVA)
which(etasq(ANOVA)>.1)

