# Bootstrap sCCA to identify relevant structural and functional connections that play an important role to the prediction. Prediction of functional connectivity from structural connectivity (one way prediction). Only the functional connectivity matrices are expected to be SPD. 



#DWIData: A list with the structural connectivity matrices. The length of the list will be equal to the number of subjects
#fMRIData: A list with the functional connectivity matrices. The length of the list will be equal to the number of subjects
#cores: Number of cpu units that would be used in PMA.NEW package. Cores are used only in CCA.permute to optimise the penalty values
#log_flag: If TRUE the functional connectivity matrices will be projected to a common average and operations will be performed in the projected space. The prediction will be projected back to SPD space. If FALSE sCCA will be performed on the original values of functional connections
#removeAllZeros: Zeros in structural connectivity matrices are treated as missing values (NAs). If TRUE any connection with a missing value to one or more subjects will be excluded. If FALSE zeros will be used in place of the missing values.
#penaltyX: Shrinkage parameter for the structural data. It takes values from 0 (sparse) to 1 (full). If NULL this parameter it will be optimised. 
#penaltyY: Shrinkage parameter for the functional data. It takes values from 0 (sparse) to 1 (full). If NULL this parameter it will be optimised. 
identificationX2Yspd <- function( DWIData, fMRIData, info=info, cores=1, log_flag=TRUE, removeAllZeros=TRUE, flag_correct=TRUE, startIter=1, finIter=1, penaltyX=NULL, penaltyY=NULL, randomised=TRUE, standardizeC=TRUE ){
	
	outputDirMain <- info$outputDirMain
	typeOfGraph <- info$typeOfGraph
	
	namesRoisX <- rownames(DWIData[[1]]) #just use first subject here assuming similar names for the rest!
	namesRoisY <- rownames(fMRIData[[1]]) #just use first subject here assuming similar names for the rest!	
	
	subjNum <- length(DWIData)
	if( subjNum != length(fMRIData) )
		report(OL$Error,"Number of subjects must be equal across the two variables")
	
	numRCX <- dim(DWIData[[1]])
	numRCY <- dim(fMRIData[[1]])
	if( numRCX[1] != numRCX[2] || numRCY[1] != numRCY[2] )
		report(OL$Error,"Connectivity matrices should be square")
	numRoisX <- numRCX[1]
	numRoisY <- numRCY[1]
	
	if(log_flag){
		numConsY <- ( numRoisY*(numRoisY+1) )/2                                    	#THIS POINT HAS BEEN CHANGED TO SUPPORT THE DIAGONAL
		#This point is for indexing only.It will be used later
		refer_arrayY <- array( 0, c(numRoisY,numRoisY) )
		refer_indUpperY <- which( upper.tri(refer_arrayY,diag=TRUE), arr.ind=TRUE )  #THIS POINT HAS BEEN CHANGED TO SUPPORT THE DIAGONAL
		refer_indLowerY <- cbind(refer_indUpperY[,2],refer_indUpperY[,1])
		refer_indUpperY0 <- which( upper.tri(refer_arrayY,diag=TRUE) )
	} else{
		numConsY <- ( numRoisY*(numRoisY-1) )/2                                    	#THIS POINT DOES NOT SUPPORT THE DIAGONAL
		#This point is for indexing only.It will be used later
		refer_arrayY <- array( 0, c(numRoisY,numRoisY) )
		refer_indUpperY <- which( upper.tri(refer_arrayY,diag=FALSE), arr.ind=TRUE )  #THIS POINT DOES NOT SUPPORT THE DIAGONAL
		refer_indLowerY <- cbind(refer_indUpperY[,2],refer_indUpperY[,1])
		refer_indUpperY0 <- which( upper.tri(refer_arrayY,diag=FALSE) )			
	}

	#indexing for structural connections -- DIAGONAL WILL BE EXCLUDED
	numConsX <- ( numRoisX*(numRoisX-1) )/2                                    	#THIS POINT DOES NOT SUPPORT THE DIAGONAL
	#This point is for indexing only.It will be used later
	refer_arrayX <- array( 0, c(numRoisX,numRoisX) )
	refer_indUpperX <- which( upper.tri(refer_arrayX,diag=FALSE), arr.ind=TRUE )  #THIS POINT DOES NOT SUPPORT THE DIAGONAL
	refer_indLowerX <- cbind(refer_indUpperX[,2],refer_indUpperX[,1])
	refer_indUpperX0 <- which( upper.tri(refer_arrayX,diag=FALSE) )			
	
	#pre-estimate the matrix-logs for each connectivity matrix
	fMRIDataLog <- list()
	for(si in 1:subjNum){
		fMRIDataLog[[si]] <- logm( fMRIData[[si]] )
	}
	
	#estimate the log-mean for all the subjects
	fMRImeanLog <- array(0,c(numRoisY,numRoisY))
	for(si in 1:subjNum){
		fMRImeanLog <- fMRImeanLog + fMRIDataLog[[si]]
	}
	fMRImeanLog <- fMRImeanLog/subjNum
	fMRImean <- as.matrix( expm(fMRImeanLog) )	#project back to the spd matrices
	
	tmpSqm <- sqrtm(fMRImean)
	fMRImeansqm <- tmpSqm$B			#estimate the square root of the mean
	fMRImeansqmI <- tmpSqm$Binv 	#invert the square root

	Xm <- array( 0, c(subjNum, numConsX ) )
	Ym <- array( 0, c(subjNum, numConsY ) )
	count <- 1
	for(si in 1:subjNum){
		#structural data do not need projection
		Xs <- DWIData[[si]]
		Xm[si,] <- Xs[refer_indUpperX0]

		#project each subject to the average mean
		if(log_flag){
			Ys <- fMRImeanCVsqmI %*% fMRIData[[si]] %*% fMRImeanCVsqmI
			Yproj <- logm(Ys)
			Ym[si,] <- Yproj[refer_indUpperY0]	
		} else{
			Ys <- fMRIData[[si]]
			Ym[si,] <- Ys[refer_indUpperY0]
		}

	} #for(si in 1:subjNum)

	print(penaltyX)
	print(penaltyY)
	if( is.null(penaltyX) || is.null(penaltyY) ){
		
		#try load shrinkage parameter values from file
		filepathIn <- file.path(outputDirMain,"penalty.txt")
		if( file.exists(filepathIn) ){
			penalty <- as.matrix( read.table(filepathIn) )
			penaltyX <- penalty[1]
			penaltyY <- penalty[2]
		} else{
			Xmp <- Xm
			Ymp <- Ym
			if(!removeAllZeros){  	#just remove connections that are zero across all subjects -- necessary!
				tmp <- apply(abs(Xmp),2,sum)
				indexXp <- which(tmp!=0)
				Xmp <- Xmp[,indexXp]
			} else{    				#remove any connection with zeros (zeros => NAs)
				tmp <- apply( (Xmp!=0), 2, all )
				indexXp <- which(tmp)
				Xmp <- Xmp[,indexXp]	
			}
		
			res0 <- CCA.permuteCluster(Xmp,Ymp,cores=1,typex="standard",typez="standard",niter=100)
			penaltyX <- res0$bestpenaltyx
			penaltyY <- res0$bestpenaltyz
		
			penalty <- c(penaltyX,penaltyY)
			write.table(penalty,file=filepathIn,row.names=FALSE,col.names=FALSE)
		}  #if( file.exists(filepathIn) )
		
	}  #if( penaltyX==NULL || penaltyY==NULL )

	#bootstrap iterations - resample subjects with replacement
	numVariables <- array(NA,c((finIter-startIter+1),4))
	countv <- 1
	for(ni in startIter:finIter){
		indSample <- sample.int(subjNum,replace=TRUE)
		Xboot <- Xm[indSample,]
		Yboot <- Ym[indSample,]
		#need to exclude connections with zeros across all subjects
		if(!removeAllZeros){  	#just remove connections that are zero across all subjects -- necessary!
			tmp <- apply(abs(Xboot),2,sum)
			indexX <- which(tmp!=0)
			Xboot <- Xboot[,indexX]
			colnames(Xboot) <- paste("I",indexX,sep="")  #keep index information
		} else{    				#remove any connection with zeros (zeros => NAs)
			tmp <- apply((Xboot!=0),2,all)
			indexX <- which(tmp)
			Xboot <- Xboot[,indexX]	
		}
		
		####  apply CCA  ####
		# if(penaltyX==NULL || penaltyY==NULL){
		# 	res0 <- CCA.permuteCluster(Xboot,Yboot,cores=1,typex="standard",typez="standard",niter=100)
		# 	penaltyX <- res0$bestpenaltyx
		# 	penaltyY <- res0$bestpenaltyz
		# }
		res <- CCA( Xboot, Yboot, typex="standard", typez="standard", K=1, niter=3000, penaltyx=penaltyX, penaltyz=penaltyY,randomised=randomised)  #sparce canonical correlation
		
		numVariables[countv,] <- c(dim(Xboot)[2],dim(Yboot)[2], length(which(res$u!=0)), length(which(res$v!=0)))
		countv <- countv+1
		
		Xcoef <- array( 0, c(1, numConsX ) )
		Xcoef[1,indexX] <- res$u
		Ycoef <- res$v

		Xrec <- array( 0, c(numRoisX,numRoisX))
		Xrec[refer_indUpperX] <- Xcoef
		Xrec[refer_indLowerX] <- Xcoef
		Yrec <- array( 0, c(numRoisY,numRoisY))
		Yrec[refer_indUpperY] <- Ycoef
		Yrec[refer_indLowerY] <- Ycoef
		#save results for each iteration
		rownames(Xrec) <- namesRoisX
		colnames(Xrec) <- namesRoisX
		rownames(Yrec) <- namesRoisY
		colnames(Yrec) <- namesRoisY
		write.table(Xrec,file=file.path(outputDirMain,paste("Xrec_",typeOfGraph,"_",as.character(ni),".txt",sep="") ), quote=FALSE )
		write.table(Yrec,file=file.path(outputDirMain,paste("Yrec_",typeOfGraph,"_",as.character(ni),".txt",sep="") ), quote=FALSE )
		
		if(flag_correct){
			#normalise Xboot and Yboot
		    if(standardizeC){
		      sdxboot <- apply(Xboot,2,sd)
		      sdyboot <- apply(Yboot,2,sd)
		      Xboot <- scale(Xboot,TRUE,sdxboot)
		      Yboot <- scale(Yboot,TRUE,sdyboot)
		    }
			Xcov <- cov(Xboot)
			Ycov <- cov(Yboot)
			XcoefCorrect <- array( 0, c(1, numConsX ) )
			XcoefCorrect[1,indexX] <- Xcov %*% res$u 
			YcoefCorrect <- Ycov %*% res$v

			XrecCor <- array( 0, c(numRoisX,numRoisX))
			XrecCor[refer_indUpperX] <- XcoefCorrect
			XrecCor[refer_indLowerX] <- XcoefCorrect
			YrecCor <- array( 0, c(numRoisY,numRoisY))
			YrecCor[refer_indUpperY] <- YcoefCorrect
			YrecCor[refer_indLowerY] <- YcoefCorrect
			
			rownames(XrecCor) <- namesRoisX
			colnames(XrecCor) <- namesRoisX
			rownames(YrecCor) <- namesRoisY
			colnames(YrecCor) <- namesRoisY
			
			#save results for each iteration
			write.table(XrecCor,file=file.path(outputDirMain,paste("XrecCor_",typeOfGraph,"_",as.character(ni),".txt",sep="") ), quote=FALSE )
			write.table(YrecCor,file=file.path(outputDirMain,paste("YrecCor_",typeOfGraph,"_",as.character(ni),".txt",sep="") ), quote=FALSE )
		}  #if(flag_correct)
	} #for(ni in startIter:finIter)
	
	write.table(numVariables,file=file.path(outputDirMain,paste("numVariables_",typeOfGraph,"_",as.character(startIter),"_",as.character(finIter),".txt",sep="") ), quote=FALSE, row.names=FALSE, col.names=FALSE)
}


