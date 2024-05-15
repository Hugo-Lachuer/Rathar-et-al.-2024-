#Center of Mass
CenterMass <- function(Image){
	Xmatrix <- matrix(rep(1:ncol(Image), nrow(Image)), ncol=ncol(Image), nrow=nrow(Image), byrow=TRUE)
	Ymatrix <- matrix(rep(1:nrow(Image), ncol(Image)), ncol=ncol(Image), nrow=nrow(Image), byrow=FALSE)
	Xmass <- sum(Xmatrix * Image)/sum(Image)
	Ymass <- sum(Ymatrix * Image)/sum(Image)
	return(c(Xmass, Ymass))
}

#Remove holes from a mask
RemoveHoles <- function(I){
	
	I <- rbind(rep(0,length=ncol(I)), I, rep(0,length=ncol(I)))
	I <- cbind(rep(0,length=nrow(I)), I, rep(0,length=nrow(I)))

	I <- rbind(rep(1,length=ncol(I)), I, rep(1,length=ncol(I)))
	I <- cbind(rep(1,length=nrow(I)), I, rep(1,length=nrow(I)))
	
	Row <- 2
	Col <- 2
	
	#FloodFill
	P <- matrix(,ncol=2)
	P[1,] <- c(Row,Col)
	I[Row,Col] <- 2
	Pfinal <- P
	
	while(nrow(P) > 0){
	
		P1 <- P
		P <- matrix(,nrow=0,ncol=2)
		
		#4-connexe
		for(i in 1:nrow(P1)){
			if(I[P1[i,1]+1,P1[i,2]] == 0){
				I[P1[i,1]+1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]+1,P1[i,2]))
			}
			if(I[P1[i,1]-1,P1[i,2]] == 0){
				I[P1[i,1]-1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]-1,P1[i,2]))
			}
			if(I[P1[i,1],P1[i,2]+1] == 0){
				I[P1[i,1],P1[i,2]+1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]+1))
			}
			if(I[P1[i,1],P1[i,2]-1] == 0){
				I[P1[i,1],P1[i,2]-1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]-1))
			}
		}
		if(nrow(P) > 0)	{Pfinal <- rbind(Pfinal,P)}
	}
	
	I2 <- matrix(1, ncol=ncol(I), nrow=nrow(I))
	for(i in 1:nrow(Pfinal)){
		I2[Pfinal[i,1], Pfinal[i,2]] <- 0
	}

	I2 <- I2[-c(nrow(I2)-1, nrow(I2)),]
	I2 <- I2[-c(1,2),]
	I2 <- I2[,-c(ncol(I2)-1, ncol(I2))]
	I2 <- I2[,-c(1,2)]
	
	return(I2)
}


#FloodFill
FloodFill <- function(I,Row,Col){

	I <- rbind(rep(0,length=ncol(I)), I, rep(0,length=ncol(I)))
	I <- cbind(rep(0,length=nrow(I)), I, rep(0,length=nrow(I)))
	Row <- Row + 1
	Col <- Col + 1

	if(I[Row,Col] == 1){
		P <- matrix(,ncol=2)
		P[1,] <- c(Row,Col)
		I[Row,Col] <- 2
	} else {
		A1 <- which(I==1, arr.ind=TRUE)
		Row <- A1[1,1]
		Col <- A1[1,2]
		P <- matrix(,ncol=2)
		P[1,] <- c(Row,Col)
		I[Row,Col] <- 2
	}
	
	Pfinal <- P
	
	while(nrow(P) > 0){
	
		P1 <- P
		P <- matrix(,nrow=0,ncol=2)
		
		#4-connexe
		for(i in 1:nrow(P1)){
			if(I[P1[i,1]+1,P1[i,2]] == 1){
				I[P1[i,1]+1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]+1,P1[i,2]))
			}
			if(I[P1[i,1]-1,P1[i,2]] == 1){
				I[P1[i,1]-1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]-1,P1[i,2]))
			}
			if(I[P1[i,1],P1[i,2]+1] == 1){
				I[P1[i,1],P1[i,2]+1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]+1))
			}
			if(I[P1[i,1],P1[i,2]-1] == 1){
				I[P1[i,1],P1[i,2]-1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]-1))
			}
		}
		if(nrow(P) > 0)	{Pfinal <- rbind(Pfinal,P)}
	}

	I2 <- matrix(0, ncol=ncol(I), nrow=nrow(I))
	for(i in 1:nrow(Pfinal)){
		I2[Pfinal[i,1], Pfinal[i,2]] <- 1
	}
	
	I2 <- I2[-nrow(I2),]
	I2 <- I2[-1,]
	I2 <- I2[,-ncol(I2)]
	I2 <- I2[,-1]
	
	return(I2)
}


#Generate CSR cell
CSRCell <- function(n, Mask){

	#n: Number giving the number of points
	#Mask: a binary matrix (image) giving the shape of the cell
	
	random_cell <- matrix(, ncol=2, nrow=n)
	
	if(is.owin(Mask)){
		Mask <- Mask[[10]]*1
	}
	for(i in 1:n){
		First <- TRUE
		x <- y <- 1
		while(Mask[y,x]==0 | First == TRUE){
			First <- FALSE
			x <- sample(1:ncol(Mask),1)
			y <- sample(1:nrow(Mask),1)
		}
		random_cell[i,] <- c(x,y)
	}

	return(random_cell)
}


#AverageCurve
AverageCurve <- function(DataList){

	#DataList: List of matrices with two columns (x and y)

	#Remove NULL
	Index <- which(sapply(DataList, is.null))
	if(length(Index)>0){
		DataList <- DataList[-Index]
	}

	Rmax <- vector(, length=length(DataList))
	for(i in 1:length(Rmax)){
		Rmax[i] <- max(DataList[[i]][,1])
	}
	Rmin <- vector(, length=length(DataList))
	for(i in 1:length(Rmin)){
		Rmin[i] <- min(DataList[[i]][,1])
	}
	Rindex <- seq(median(Rmin), median(Rmax), (median(Rmax)-median(Rmin))/1000)
	
	#Linear interpolation
	DataMatrix <- matrix(,ncol=length(DataList), nrow=length(Rindex))
	for(i in 1:length(DataList)){
		DataMatrix[,i] <- approx(x=DataList[[i]][,1], y=DataList[[i]][,2], xout=Rindex, method="linear")$y
	}
	
	#Average and SEM
	DataAverage  <- DataSEM <- vector(, length=length(Rindex))
	for(i in 1:length(Rindex)){
		DataAverage[i] <- mean(DataMatrix[i,], na.rm=TRUE)
		DataSEM[i] <- sd(DataMatrix[i,],na.rm=TRUE)/sqrt(length(which(is.na(DataMatrix[i,])==FALSE)))
	}
	
	#Output
	Output <- list()
	Output[[1]] <- Rindex
	Output[[2]] <- DataAverage
	Output[[3]] <- DataSEM
	return(Output)
}

#Plot circle
PlotCircle <- function(X, Y, R, Nx, Ny){
	Map <- matrix(0, ncol=Nx, nrow=Ny)
	R <- round(R)
	X <- round(X)
	Y <- round(Y)
	for(i in (X-R):(X+R)){
		for(j in (Y-R):(Y+R)){
			if( (i>0) & (i<=Nx) & (j>0) & (j<=Ny) & ((i-X)^2+(j-Y)^2 <= R^2)){
				Map[j,i] <- 1
			}
		}
	}
	return(Map)
}