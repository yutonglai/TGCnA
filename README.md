# Temporal Gene Coexpression Network Analysis Using A Low-rank plus Sparse Framework

Code for 

Li, Jinyu, Yutong Lai, Chi Zhang, and Qi Zhang. "Temporal Gene Coexpression Network Analysis Using A Low-rank plus Sparse Framework." BioRxiv (2018): 359612.

## Overview


Various gene network models with distinct physical nature have been widely used in biological studies. For temporal transcriptomic studies, the current dynamic models either ignore the temporal variation in the network structure or fail to scale up to a large number of genes due to severe computational bottlenecks and sample size limitation. Although the correlation-based gene networks are computationally affordable, they have limitations after being applied to gene expression time-course data. 

We proposed Temporal Gene Coexpression Network Analysis (TGCnA) framework for the transcriptomic time-course data. The mathematical nature of TGCnA is the joint modeling of multiple covariance matrices across time points using a "low-rank plus sparse" framework, in which the network similarity across time points is explicitly modeled in the low-rank component. 
We demonstrated the advantage of TGCnA in covariance matrix estimation and gene module discovery using both simulation data and real transcriptomic data.

The following is a brief tutorial of the proposed method.

## Tutorial

#### 1. Source the function
```
source("functionsNew.R")
```
#### 2. Analysis of simulated continuous gene expression data (Microarray-like)
```
library(pdfCluster)

load("Data_t5_c5_r15.rdata") # load a dataset simulated as specified in the paper
load("Data_t5_c5_r15_modules_True.rdata") # Load the simulated true modules, used as the ground truth for evaluating module discovery accuracy.

tvec    = c( 1, 2, 3, 4, 5) # The simulated data contain five time points 
numreps = c(15,15,15,15,15) # The simulated data hs 15 replicates at each time point

# Perform TGCnA analysis
# It returns the list of correlation matrices for each time point
correlationMatrixSim = NULL
correlationMatrixSim = TGCN(tvec, numreps, Data_t5_c5_r15)

# (Opotional) Identify modules for each time point
newdataSim = as.data.frame(Data_t5_c5_r15)
modresSim  = TGCN_Cluster(correlationMatrixSim, newdataSim)
Data_t5_c5_r15_modules_Sim  = lapply(1:nt, function(tt) modresSim$modres[[tt]])

# Evalute the module discovery performance by Adjust Rand Index between the simulated true modules and the modules identified from the simulated data
AdjRandIndex = NULL
for(tt in 1:length(tvec))	
    AdjRandIndex[tt] = adj.rand.index(Data_t5_c5_r15_modules_Sim[[tt]], Data_t5_c5_r15_modules_True[[tt]])
```

#### 3. The analysis of simulated count data (RNA-Seq)
```
## Simulate RNA-Seq data with the same group structure as the simulated Microarray data
n = nrow(Data_t5_c5_r15)*ncol(Data_t5_c5_r15)
a = 20
b = 0.2
lambda = a*exp(b*Data_t5_c5_r15) ## define the expectations of the simulation models for the count data by transforming the simulated continuous data
y = matrix(rpois(n,lambda),nrow(Data_t5_c5_r15),ncol(Data_t5_c5_r15)) ## simulate count data.
y = log(y+1) ## use log(count + 1) as the input of TGCnA.

# Perform TGCnA analysis
# It returns the list of correlation matrices for each time point
correlationMatrixTran = NULL
correlationMatrixTran = TGCN(tvec, numreps, y)

# (Opotional) Identify modules for each time point
newdataTran = as.data.frame(y)
modresTran  = TGCN_Cluster(correlationMatrixTran, newdataTran)
Data_t5_c5_r15_modules_Tran  = lapply(1:nt, function(tt) modresTran$modres[[tt]])

# Evalute the module discovery performance by Adjust Rand Index between the simulated true modules and the modules identified from the simulated data
AdjRandIndexTran = NULL
for(tt in 1:length(tvec))	
    AdjRandIndexTran[tt] = adj.rand.index(Data_t5_c5_r15_modules_Tran[[tt]], Data_t5_c5_r15_modules_True[[tt]])
```


![corrlistT=3](https://user-images.githubusercontent.com/46899273/56169175-ebfbd000-5fa2-11e9-8aa0-5e11080ba864.PNG)
##### Figure 1. The true correlation matrix (True), the recovered correlation matrix by TGCnA with simulated continuous data data (Continuous) and counts data (Counts)


##### Table 1. The Adjusted Rand Index for modules recovery by TGCnA with continous data (Continous) and counts data (Counts)
Data / Time |Time=1        |Time=2        |Time=3        |Time=4        |Time=5
-----------:|:------------:|:------------:|:------------:|:------------:|:------------
Continous         |0.5166732     |0.5875409     |0.5554822     |0.4629107     |0.8039664
Counts        |0.6288305     |0.6471858     |0.6928232     |0.5391654     |0.4025173



#### 4. Real data analysis
We load the data, and provide a protocol for analyzing it.

Load the data.
```
readRDS("braindata10193A.rds")
readRDS("braindata10193B.rds")

## data
braindata10193 = rbind(braindata10193A, braindata10193B)
```
##### Table 2. The human BrainSpan-RNA sequence data
GeneID / Sample |1        |2        |3        |4        |5        |6        |......  |524        
---------------:|:-------:|:-------:|:-------:|:-------:|:-------:|:-------:|:------:|:--------
ENSG00000000003 |5.226783 |4.658283 |4.345572 |4.841400 |4.392196 |3.970916 |......  |2.029566     
ENSG00000000419 |5.144586 |4.443982 |4.302681 |4.546363 |4.338313 |3.587409 |......  |4.900434     
ENSG00000001036 |3.471608 |3.072748 |2.971666 |3.208103 |3.200653 |3.031687 |......  |2.757560
......          |......   |......   |......   |......   |......   |......   |......  |......
ENSG00000259770 |3.555986 |3.458692 |3.094457 |3.282903 |3.154652 |2.900298 |......  |3.080229

```
## Meta information of the brain data.
tvec    = c(1.98, 2.80, 3.03, 3.79, 5.04, 8.11, 13.33, 20.56, 43.97, 117.63, 206.37, 384.30) #developmental stages (weeks)
numreps = c(  30,   45,   44,   53,   43,   22,    33,    26,    44,     41,     50,     93) #sample size at each time point
res.mtx = ResidMtx(braindata10193, numreps)  #residual matrix (subtract mean)
res.mtx = Normalize.res.mtx(res.mtx,numreps) #normalize (devide d1 and multiply sqrt(p))

# Perform TGCnA analysis
# It returns the list of correlation matrices for each time point
correlationMatrix10193 = NULL
correlationMatrix10193 = TGCN(tvec, numreps, res.mtx)

# (Opotional) Identify modules for each time point
newdata10193 = as.data.frame(res.mtx)
modres10193  = TGCN_Cluster(correlationMatrix10193, newdata10193)
```



