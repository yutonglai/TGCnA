source("functionsNew.R")

#### Example

#### #########################################################################################################################################################
#### #### 1. Simulation Analysis #### ########################################################################################################################
#### #########################################################################################################################################################

#### 1.a. The Simulated Data Analysis. 
##### The simulated data set is a 2000 by 75 matrix of 2000 genes, 5 time points and 15 replications for each time point with 5 clusters 
load("Data_t5_c5_r15.rdata")
tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
numreps = c(15,15,15,15,15) # set the number of replication for each time point

##### Returen correlation matrix for each time point
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

##### Returen the transformed data correlation matrix for each time point
correlationMatrixTran = NULL
correlationMatrixTran = TGCN(tvec, numreps, DataTran_t5_c5_r15)

##### (Opotional) Return the transformed data cluster list for each time point
newdataTran = as.data.frame(DataTran_t5_c5_r15)
modresTran  = TGCN_Cluster(correlationMatrixTran, newdataTran)


#### 1.c. Result plots
load("Pure_t5_c5_r15_s40.rdata") # 40 times simulated data (5 time point, 5 clusters, 15 samples\time point) analysis
# Pure_t15_c5_r15_s40[[1]]: The Frobenius Loss from TGCnA at each time point
# Pure_t15_c5_r15_s40[[2]]: The Frobenius Loss from Naive method at each time point
# Pure_t15_c5_r15_s40[[3]]: The Adjusted Rand Index from TGCnA at each time point
# Pure_t15_c5_r15_s40[[4]]: The Adjusted Rand Index from Naive method at each time point


##### Adjusted Rand Index plot for the simulated data analysis
testaARI = PureTran_t5_c5_r15_s40[[3]]
testbARI = PureTran_t5_c5_r15_s40[[4]]

dataframePro   = data.frame("TGCnA",testaARI)
dataframeNaive = data.frame("Naive",testbARI)

colnames(dataframePro)   = c("label","1","2","3","4","5")
colnames(dataframeNaive) = c("label","1","2","3","4","5")
dataframe=rbind(dataframePro,dataframeNaive)
datamelt=melt(dataframe,id.var="label")

dataframeProMedian=data.frame("ProMedians",colMedians(testaARI))
dataframeNaiveMedian=data.frame("NaiveMedians",colMedians(testbARI))
ggplot(data = datamelt, aes(x=variable, y=value)) + theme_bw() + geom_boxplot(aes(fill = label), group=1) + 
      geom_line(data = dataframeProMedian, aes(linetype = " TGCnA", x = 1:ncol(testaARI), y = colMedians(testaARI), group = 1,size=0.1), size=0.5) + 
      geom_line(data = dataframeNaiveMedian, aes(linetype = "Naive",x = 1:ncol(testaARI), y = colMedians(testbARI), group = 1,size=0.1), size=0.5) +
      labs(x="",y="Adusted Rand Index",title="T=5 C=5 n=15",linetype="") + guides(fill=guide_legend(title=""), shape = guide_legend(override.aes = list(size = 15))) + 
      expand_limits(y=0) + theme(legend.position = "none",  axis.title.y=element_text(size=10), legend.key.width = unit(4,"line"), legend.key.height = unit(2,"line"), legend.text = element_text(size = 15), plot.title=element_text(hjust=0,size=10))



#### #########################################################################################################################################################
#### #### 2. Real Data Analysis #### #########################################################################################################################
#### #########################################################################################################################################################
                                  
#### Real data analysis, 3141 genes
load("braindata3141.rdata") # 3141 genes
nt = length(numreps)
                                    
# return the correlation matrix at each time point     
correlationMatrix3141 = NULL
correlationMatrix3141 = TGCN(tvec, numreps, res.mtx_Original)
                                    
# (Optional) Return cluster list for each time point
newdata3141 = as.data.frame(res.mtx_Original)
modresTimePro3141 = TGCN_Cluster(correlationMatrix3141, newdata3141)

  
