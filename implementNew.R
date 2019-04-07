source("functionsNew.R")

#### Example

#### #### 1. Simulation Analysis #### ####

#### 1.a. The Simulated Data Analysis. 
##### The simulated data set is a 2000 by 75 matrix of 2000 genes, 5 time points and 15 replications for each time point with 5 clusters 
load("Data_t5_c5_r15.rdata")
tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
numreps = c(15,15,15,15,15) # set the number of replication for each time point

##### Returen correlation list for each time point
correlationMatrix = NULL
correlationMatrix = TGCN(tvec, numreps, Data_t5_c5_r15)

##### (Opotional) Return cluster list for each time point
newdata = as.data.frame(Data_t5_c5_r15)
modres  = TGCN_Cluster(correlationMatrix, newdata)


#### 1.b. The Transformed Simulated Data Analysis.
load("Data_t5_c5_r15.rdata")
tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
numreps = c(15,15,15,15,15) # set the number of replication for each time point

##### Transform the oringal simulated data set Data_t5_c5_r15.rdata to a count data set.
n = nrow(Data_t5_c5_r15)*ncol(Data_t5_c5_r15) # obtain the data set size
a = 20
b = 0.2
lambda = a*exp(b*Data_t5_c5_r15) # generate the Poisson distribution parameters for each data point 
DataTran_t5_c5_r15 = matrix(rpois(n,lambda),nrow(Data_t5_c5_r15),ncol(Data_t5_c5_r15)) # simulate count data from the Poisson distribution
DataTran_t5_c5_r15 = log(y+1) # supress the count data set to the original scale

##### Returen the transformed data correlation list for each time point
correlationMatrixTran = NULL
correlationMatrixTran = TGCN(tvec, numreps, DataTran_t5_c5_r15)

##### (Opotional) Return the transformed data cluster list for each time point
newdataTran = as.data.frame(DataTran_t5_c5_r15)
modresTran  = TGCN_Cluster(correlationMatrixTran, newdataTran)


#### 1.c. Frobenius Loss






