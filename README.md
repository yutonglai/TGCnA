# Temporal Gene Coexpression Network Analysis Using A Low-rank plus Sparse Framework

Jinyu Li, Yutong Lai, Chi Zhang, Qi Zhang 

## overview

Various gene network models with distinct physical nature have been widely used in biological studies. For temporal transcriptomic studies, the current dynamic models either ignore the temporal variation in the network structure or fail to scale up to a large number of genes due to severe computational bottlenecks and sample size limitation. On the other hand, correlation-based gene networks are more computationally more affordable, but have not been properly extended to gene expression time-course data.

We propose Temporal Gene Coexpression Network (TGCN) for the transcriptomic time-course data. The mathematical nature of TGCN is the joint modeling of multiple covariance matrices across time points using a “low-rank plus sparse” framework, in which the network similarity across time points is explicitly modeled in the low-rank component. Using both simulations and a real data application, we showed that TGCN improved the covariance estimation loss and identified more robust and interpretable gene modules.

## Related codes

#### 1. Source the function
```
source("functionsNew.R")
```
#### 2. Simulation analysis
```
library(pdfCluster)

load("Data_t5_c5_r15.rdata")
load("Data_t5_c5_r15_modules_True.rdata") # The real modules

## Initialization
tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
numreps = c(15,15,15,15,15) # set the number of replication for each time point

## Return correlation list for each time point
correlationMatrixSim = NULL
correlationMatrixSim = TGCN(tvec, numreps, Data_t5_c5_r15)

## (Opotional) Return module list for each time point
newdataSim = as.data.frame(Data_t5_c5_r15)
modresSim  = TGCN_Cluster(correlationMatrixSim, newdataSim)
Data_t5_c5_r15_modules_Sim  = lapply(1:nt, function(tt) modresSim$modres[[tt]])

## Adjust Rand Index
AdjRandIndex = NULL
for(tt in 1:length(tvec))
{		
	AdjRandIndex[tt] = adj.rand.index(Data_t5_c5_r15_modules_Sim[[tt]], Data_t5_c5_r15_modules_True[[tt]])
}
print(AdjRandIndex)
```

#### 3. Transformed simulation analysis: transform the simulation data to count-type data, then redo the simulation
```
n = nrow(Data_t5_c5_r15)*ncol(Data_t5_c5_r15)
a = 20
b = 0.2
lambda = a*exp(b*Data_t5_c5_r15)
y = matrix(rpois(n,lambda),nrow(Data_t5_c5_r15),ncol(Data_t5_c5_r15))
y = log(y+1)

## Return correlation list for each time point
correlationMatrixTran = NULL
correlationMatrixTran = TGCN(tvec, numreps, y)

## (Opotional) Return module list for each time point
newdataTran = as.data.frame(y)
modresTran  = TGCN_Cluster(correlationMatrixTran, newdataTran)
Data_t5_c5_r15_modules_Tran  = lapply(1:nt, function(tt) modresTran$modres[[tt]])

## Adjust Rand Index
AdjRandIndexTran = NULL
for(tt in 1:length(tvec))
{		
	AdjRandIndexTran[tt] = adj.rand.index(Data_t5_c5_r15_modules_Tran[[tt]], Data_t5_c5_r15_modules_True[[tt]])
}
print(AdjRandIndexTran)
```
Time |T=1        |T=2        |T=3        |T=4        |T=5
----:|:---------:|:---------:|:---------:|:---------:|:---------
Sim  |0.5166732  |0.5875409  |0.5554822  |0.4629107  |0.8039664
Tran |0.6288305  |0.6471858  |0.6928232  |0.5391654  |0.4025173
![corrlistT=3](https://user-images.githubusercontent.com/46899273/56168094-c4573880-5f9f-11e9-93c0-fe95de7430af.PNG)

#### 4. Real data analysis
```
load("braindata10193A.rds")
load("braindata10193B.rds")

## data
braindata10193 = rbind(braindata10193A, braindata10193B)

## initialization
tvec    = c(1.98, 2.80, 3.03, 3.79, 5.04, 8.11, 13.33, 20.56, 43.97, 117.63, 206.37, 384.30) #developmental stages (weeks)
numreps = c(  30,   45,   44,   53,   43,   22,    33,    26,    44,     41,     50,     93) #sample size at each time point
res.mtx = ResidMtx(braindata10193, numreps)  #residual matrix (subtract mean)
res.mtx = Normalize.res.mtx(res.mtx,numreps) #normalize (devide d1 and multiply sqrt(p))

## Return correlation list for each time point
correlationMatrix10193 = NULL
correlationMatrix10193 = TGCN(tvec, numreps, res.mtx)

## (Opotional) Return module list for each time point
newdata10193 = as.data.frame(res.mtx)
modres10193  = TGCN_Cluster(correlationMatrix10193, newdata10193)
```



