source("functionsNew.R")

#### True modules
  load("Data_t5_c5_r15_corrlistTrue.rdata")
  load("Data_t5_c5_r15.rdata")
	
  ## (Opotional) Return cluster list for each time point
	newdataSim = as.data.frame(Data_t5_c5_r15)
	modresTrue = TGCN_Cluster(corrlistTrue, newdataSim)
	
	## Heatmap for True correlation matrix
	heatmap(Data_t5_c5_r15_corrlistTrue[[1]],col = rev(rainbow(40)))	

	

#### Simulated data (2000 by 75 matrix of 2000 genes, 5 time points and 15 replications for each time point with 5 clusters) 
	load("Data_t5_c5_r15.rdata")
	tvec    = c( 1, 2, 3, 4, 5) # set the time-lapses from the beginning
	numreps = c(15,15,15,15,15) # set the number of replication for each time point

	## Return correlation list for each time point
	correlationMatrixSim = NULL
	correlationMatrixSim = TGCN(tvec, numreps, Data_t5_c5_r15)

	## Heatmap for correlation matrix from simulated data
	heatmap(correlationMatrixSim[[1]],col = rev(rainbow(40)))
	
	## (Opotional) Return cluster list for each time point
	newdataSim = as.data.frame(Data_t5_c5_r15)
	modresSim  = TGCN_Cluster(correlationMatrixSim, newdataSim)
	


#### Transformed simulated data
	load("Data_t5_c5_r15.rdata")
	n = nrow(Data_t5_c5_r15)*ncol(Data_t5_c5_r15)
	a = 20
	b = 0.2
	lambda = a*exp(b*Data_t5_c5_r15)
	y = matrix(rpois(n,lambda),nrow(Data_t5_c5_r15),ncol(Data_t5_c5_r15))
	max(y)
	min(y)
	length(which(y==0))/length(y)
	y = log(y+1)
	
	## Initialization
	tvec    = c( 1, 2, 3, 4, 5) # time points
	numreps = c(15,15,15,15,15) # replicates

	## Return correlation list for each time point
	correlationMatrixTran = NULL
	correlationMatrixTran = TGCN(tvec, numreps, y)

	## Heatmap for correlation matrix from transformed data
	heatmap(correlationMatrixTran[[1]], col = rev(rainbow(40)), name='level')
		
	## (Opotional) Return cluster list for each time point
	newdataTran = as.data.frame(y)
	modresTran  = TGCN_Cluster(correlationMatrixTran, newdataTran)
