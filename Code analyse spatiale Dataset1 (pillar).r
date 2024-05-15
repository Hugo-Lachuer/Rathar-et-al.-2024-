#######################################################################################
##Spatial analysis of actin structures (Dataset1 = Pillar)
##Hugo Lachuer (last update 15/05/2024)
##R version 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting"
#######################################################################################

#######################################################################################
##1/ Functions and packages

#Packages
library(raster)
library(spatstat)

#Import functions
setwd("")	#Path for Functions.R
source("Functions.R")

#######################################################################################
##2/ Import data

#Parameters
Nsim <- 100	#Number of Monte-Carlo simulations

#Get images dimensions
Dimensions <- read.delim("/Dataset1/ImageDimensions.txt")	#Path


#List of files
setwd("/Dataset1")	#Path for data (Dataset1 pillar)
V1 <- list.files(, pattern = ".tif")
V1 <- substr(V1, 1, nchar(V1)-9)	#List of files names
K_List <- Kcross_List <- list()	#List of 2D matrices (x,y)
K_CSR_List <- Kcross_CSR_List <- list()	#List of a list of 2D matrices (x,y)
PPP_list <- list() #list of point pattern

#Offset for FloodFill if center of mass is in a 0 area
offsetX <- c(rep(0,6), 50, rep(0,8))

#Save plot
pdf("Individual plot Dataset1.pdf")

#Loop over all the cells
for(i in 1:length(V1)){

	print(paste0("Analyze of cell ", i, "/", length(V1)))

	#load cell mask
	Mask <- as.matrix(raster(paste0(V1[i],"_mask.tif")))
	Mask <- Mask/max(Mask)
	Center <- round(CenterMass(I=Mask))
	Mask <- FloodFill(I=Mask,Row=Center[2],Col=Center[1]+offsetX[i])	#clear outside
	Mask <- RemoveHoles(I=Mask)
	#plot(raster(Mask), col=c("white","black"))

	#Get conversion factor
	Ind <- which(Dimensions[,2]==paste0(V1[i],"_mask.tif"))
	ConversionFactor <- Dimensions[Ind,3]/Dimensions[Ind,4]	#Pixel per µm

	#Import coordinates
	Actin <- read.csv(paste0(V1[i],".csv"), sep=";")
	Actin <- Actin[,-1]
	Actin <- Actin*ConversionFactor	#Conversion

	#Plot
	#plot(raster(Mask), col=c("white","black"))
	#points(x=Actin[,1]/ncol(Mask), y=1-Actin[,2]/nrow(Mask), pch=16, col="red")

	#Load pillars
	Grid <- read.csv(paste0(V1[i], "_grid.csv"))	#Coordinates of 4 pillars
	L <- mean(Grid$Feret)*ConversionFactor	#Diameter of a pillar
	Grid <- cbind(Grid$X, Grid$Y)*ConversionFactor
	DeltaXY <- 1.8*ConversionFactor	#Spacing between pilar centers is 1.8µm

	#Generate Pillar array
	SeqX <- seq( from=Grid[1,1] - ceiling(Grid[1,1]/DeltaXY)*DeltaXY, to= Grid[1,1]  + ceiling((ncol(Mask)-Grid[1,1])/DeltaXY)*DeltaXY, by=DeltaXY)
	SeqY <- seq( from=Grid[1,2] - ceiling(Grid[1,2]/DeltaXY)*DeltaXY, to= Grid[1,2]  + ceiling((ncol(Mask)-Grid[1,2])/DeltaXY)*DeltaXY, by=DeltaXY)
	Array <- matrix(0, ncol=ncol(Mask), nrow=nrow(Mask))
	for(j in SeqX){
		for(k in SeqY){
			A1 <- PlotCircle(X=j, Y=k, R=L/2, Nx=ncol(Mask), Ny=nrow(Mask))
			Array <- Array + A1
		}
	}

	#Delete actin structures out of the mask
	for(j in nrow(Actin):1){
		if(Mask[round(Actin[j,2]), round(Actin[j,1])] == 0){
			Actin <- Actin[-j,]
		}
	}

	#Final Mask
	MaskFinal <- Mask - (Array & Mask)
	plot(raster(MaskFinal), col=c("white","black"), main=V1[i])
	points(x=Actin[,1]/ncol(MaskFinal), y=1-Actin[,2]/nrow(MaskFinal), pch=16, col="red")


	#######################################################################################
	##3/ Basic Ripley's K function

	#Observed Ripley's K function
	w <- owin(xrange=c(1,ncol(Mask))*1/ConversionFactor, yrange=c(1,nrow(Mask))*1/ConversionFactor, mask=(Mask==1))	#µm units
	p <- ppp(x=Actin[,1]*1/ConversionFactor, y=Actin[,2]*1/ConversionFactor, window=w)
	PPP_list[[i]] <- p
	A2 <- Kest(p, correction="best")
	#plot(x=A2$r, y=A2$trans - pi*(A2$r)^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K function - pir²", xlim=c(0,3), ylim=c(-5,5))


	#Monte-Carlo simulations
	R_CSR <- Ripley_CSR <- list()
	for(j in 1:Nsim){
		ActinCSR <- CSRCell(n=nrow(Actin), Mask=MaskFinal)
		p <- ppp(x=ActinCSR[,1]*1/ConversionFactor, y=ActinCSR[,2]*1/ConversionFactor, window=w)
		A3 <- Kest(p, correction="best")
		Ripley_CSR[[j]] <- cbind(A3$r, A3$trans)
		R_CSR[[j]] <- A3$r
	}
	
	#Save Ripley
	K_List[[i]] <- cbind(A2$r, A2$trans)
	K_CSR_List[[i]] <- Ripley_CSR

	#Interpolation
	R <- unique(unlist(R_CSR))
	K_CSR <- matrix(,ncol=length(R), nrow=Nsim)
	for(j in 1:length(Ripley_CSR)){
		K_CSR[j,] <- approx(x=R_CSR[[j]], y=Ripley_CSR[[j]][,2], xout=R, method="linear")$y
	}

	#Median, Up, Down
	Median <- Up <- Down <- vector()
	for(j in 1:length(R)){
		Median[j] <- median(K_CSR[,j])
		Up[j] <- max(K_CSR[,j])
		Down[j] <- min(K_CSR[,j])
	}


	#Plot
	plot(x=A2$r, y=A2$trans - pi*(A2$r)^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K function - pir²", main=V1[i]) #, xlim=c(0,3), ylim=c(-5,5))
	lines(x=R, y=Median - pi*R^2, lwd=2, col="blue")
	xshade <- c(R[which(is.na(Up)==FALSE)],rev(R[which(is.na(Up)==FALSE)]))
	yshade <- c(Up[which(is.na(Up)==FALSE)]-pi*R[which(is.na(Up)==FALSE)]^2,rev(Down[which(is.na(Up)==FALSE)]-pi*R[which(is.na(Up)==FALSE)]^2))
	polygon(xshade,yshade, border = NA, col=adjustcolor("blue", 0.4))
	legend(x="topleft", lwd=2, col=c("Red", "Blue"), legend=c("Observed", "CSR simulation"))


	#######################################################################################
	##4/ Cross Ripley's K function


	#Pillar coordinates
	SeqX <- seq( from=Grid[1,1] - ceiling(Grid[1,1]/DeltaXY)*DeltaXY, to= Grid[1,1]  + ceiling((ncol(Mask)-Grid[1,1])/DeltaXY)*DeltaXY, by=DeltaXY)
	SeqY <- seq( from=Grid[1,2] - ceiling(Grid[1,2]/DeltaXY)*DeltaXY, to= Grid[1,2]  + ceiling((ncol(Mask)-Grid[1,2])/DeltaXY)*DeltaXY, by=DeltaXY)
	Grid2 <- expand.grid(SeqX, SeqY)
	Grid2 <- Grid2[-which(Grid2[,1]<0 | Grid2[,1]>ncol(MaskFinal)),]
	Grid2 <- Grid2[-which(Grid2[,2]<0 | Grid2[,2]>nrow(MaskFinal)),]
	for(j in nrow(Grid2):1){
		if(round(Grid2[j,1])==0 | round(Grid2[j,2])==0 | round(Grid2[j,1])>ncol(Mask) | round(Grid2[j,2])>nrow(Mask)){
			Grid2 <- Grid2[-j,]
		} else {
			if(Mask[round(Grid2[j,2]), round(Grid2[j,1])] == 0){
				Grid2 <- Grid2[-j,]
			}
		}
	}


	#Plot
	#plot(raster(Mask), col=c("white","black"))
	#points(x=Actin[,1]/ncol(Mask), y=1-Actin[,2]/nrow(Mask), pch=16, col="red")
	#points(x=Grid2[,1]/ncol(Mask), y=1-Grid2[,2]/nrow(Mask), pch=16, col="blue")

	#Create marked point pattern
	X <- c(Actin[,1]*1/ConversionFactor, Grid2[,1]*1/ConversionFactor)
	Y <- c(Actin[,2]*1/ConversionFactor, Grid2[,2]*1/ConversionFactor)
	Marks <- c(rep("Actin", nrow(Actin)), rep("Pillar", nrow(Grid2)))
	Marks <- factor(Marks)
	p <- ppp(x=X, y=Y, window=w, marks=Marks)


	#Plot Kcross
	K01 <- Kcross(p, "Actin", "Pillar", correction="best")
	#plot(x=K01$r, y=K01$trans - pi*(K01$r)^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K_actin,pillar - pir²")	#, xlim=c(0,3), ylim=c(-5,5))
	K10 <- Kcross(p, "Pillar", "Actin", correction="best")
	#plot(x=K10$r, y=K10$trans - pi*(K10$r)^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K_pillar,actin - pir²")	#, xlim=c(0,3), ylim=c(-5,5))

	#Kcross = (K01+K10)/2
	if(all(K01$r!=K10$r)){print("Important problem!")}
	K_Cross <- (K01$trans+K10$trans)/2
	R_Cross <- K01$r

	#Monte-Carlo simulation
	R_Cross_CSR <- Ripley_Cross_CSR <- list()
	for(j in 1:Nsim){

		#CSR
		ActinCSR <- CSRCell(n=nrow(Actin), Mask=MaskFinal)
		X <- c(ActinCSR[,1]*1/ConversionFactor, Grid2[,1]*1/ConversionFactor)
		Y <- c(ActinCSR[,2]*1/ConversionFactor, Grid2[,2]*1/ConversionFactor)
		p <- ppp(x=X, y=Y, window=w, marks=Marks)
		
		#Kcross
		A4 <- Kcross(p, "Actin", "Pillar", correction="best")
		A5 <- Kcross(p, "Pillar", "Actin", correction="best")
		
		#Save
		if(all(A4$r!=A5$r)){print("Important problem!")}
		Ripley_Cross_CSR[[j]] <- cbind(A4$r, (A4$trans+A5$trans)/2)
		R_Cross_CSR[[j]] <- A4$r
	}
	
	
	#Save Ripley Cross
	Kcross_List[[i]] <- cbind(R_Cross, K_Cross)
	Kcross_CSR_List[[i]] <- Ripley_Cross_CSR
	

	#Interpolation
	R_Cross_CSR2 <- unique(unlist(R_Cross_CSR))
	Ripley_Cross_CSR2 <- matrix(,ncol=length(R_Cross_CSR2), nrow=Nsim)
	for(j in 1:Nsim){
		Ripley_Cross_CSR2[j,] <- approx(x=R_Cross_CSR[[j]], y=Ripley_Cross_CSR[[j]][,2], xout=R_Cross_CSR2, method="linear")$y
	}

	#Median, Up, Down
	Median <- Up <- Down <- vector()
	for(j in 1:length(R_Cross_CSR2)){
		Median[j] <- median(Ripley_Cross_CSR2[,j])
		Up[j] <- max(Ripley_Cross_CSR2[,j])
		Down[j] <- min(Ripley_Cross_CSR2[,j])
	}
	
	#Plot Kcross
	plot(x=R_Cross, y=K_Cross - pi*(R_Cross)^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K_actin,pillar - pir²", main=V1[i])	#, xlim=c(0,3), ylim=c(-5,5))
	lines(x=R_Cross_CSR2, y=Median - pi*R_Cross_CSR2^2, lwd=2, col="blue")
	xshade <- c(R[which(is.na(Up)==FALSE)],rev(R[which(is.na(Up)==FALSE)]))
	yshade <- c(Up[which(is.na(Up)==FALSE)]-pi*R[which(is.na(Up)==FALSE)]^2,rev(Down[which(is.na(Up)==FALSE)]-pi*R[which(is.na(Up)==FALSE)]^2))
	polygon(xshade,yshade, border = NA, col=adjustcolor("blue", 0.4))
	legend(x="topleft", lwd=2, col=c("Red", "Blue"), legend=c("Observed", "CSR simulation"))

}

dev.off()

#Saving
setwd("/Dataset1 (Pillar)")	#Path to "Dataset1 (Pillar)" folder
save(K_List, file="K_List.RData")
save(K_CSR_List, file="K_CSR_List.RData")
save(Kcross_List, file="Kcross_List.RData")
save(Kcross_CSR_List, file="Kcross_CSR_List.RData")
save(PPP_list, file="PPP_List.RData")


#######################################################################################
##5/ Average Ripley's K function

#Load data
setwd("/Dataset1")	#Path to "Dataset1 (Pillar)" folder
load(file="K_List.RData")
load(file="K_CSR_List.RData")

#Average Ripley K function
A1 <- AverageCurve(DataList = K_List)
R <- A1[[1]]
Mu <- A1[[2]]
SEM <- A1[[3]]

plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K - pir²", main="Average Ripley function", xlim=c(0,15), ylim=c(-100,800))
#plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K - pir²", main="Average Ripley Kcross function", xlim=c(0,3), ylim=c(-20,100))
xshade <- c(R,rev(R))
yshade <- c(Mu-pi*R^2+SEM,rev(Mu-pi*R^2-SEM))
polygon(xshade,yshade, border = NA, col=adjustcolor("red", 0.4))

#Generate CSR curves averaged over cells
R_Matrix <- Mu_Matrix <- matrix(,ncol=1001,nrow=Nsim)	#1001 because "AverageCurve" function always generates 1001 points
for(i in 1:Nsim){
	A1 <- list()
	for(j in 1:length(K_CSR_List)){
		A1[[j]] <- K_CSR_List[[j]][[i]]
	}
	A2 <- AverageCurve(DataList = A1)
	R_Matrix[i,] <- A2[[1]]
	Mu_Matrix[i,] <- A2[[2]]
}

#Check that R values are always identical
for(i in 1:ncol(R_Matrix)){
	if(length(unique(R_Matrix[,i]))>1){
		print("Problem in x axis indices!")
	}
}

#Generate empirical 95% CI
Average <- Up <- Down <- vector()
for(i in 1:ncol(Mu_Matrix)){
	Average[i] <-  mean(Mu_Matrix[,i])
	Down[i] <- quantile(Mu_Matrix[,i],probs=0.025)
	Up[i] <- quantile(Mu_Matrix[,i],probs=0.975)
}

#Plot
R <- R_Matrix[1,]
lines(x=R, y=Average - pi*R^2, lwd=2, col="blue")
xshade <- c(R,rev(R))
yshade <- c(Up-pi*R^2,rev(Down-pi*R^2))
polygon(xshade,yshade, border = NA, col=adjustcolor("blue", 0.4))

#Legend
legend(x="topleft", lwd=2, col=c("Red", "Blue"), legend=c("Observed", "CSR simulation"))


#######################################################################################
##6/ Average Ripley's K cross function

setwd("/Dataset1")	#Path to "Dataset1 (Pillar)" folder
load(file="Kcross_List.RData")
load(file="Kcross_CSR_List.RData")


#Average Ripley Kcross function
A1 <- AverageCurve(DataList = Kcross_List)
R <- A1[[1]]
Mu <- A1[[2]]
SEM <- A1[[3]]

plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K_actin,pillar - pir²", main="Average Ripley Kcross function", xlim=c(0,0.5), ylim=c(-2,2))
#plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K_actin,pillar - pir²", main="Average Ripley Kcross function", xlim=c(0,3), ylim=c(-20,20))
xshade <- c(R,rev(R))
yshade <- c(Mu-pi*R^2+SEM,rev(Mu-pi*R^2-SEM))
polygon(xshade,yshade, border = NA, col=adjustcolor("red", 0.4))

#Generate CSR curves averaged over cells
R_Matrix <- Mu_Matrix <- matrix(,ncol=1001,nrow=Nsim)	#1001 because "AverageCurve" function always generates 1001 points
for(i in 1:Nsim){
	A1 <- list()
	for(j in 1:length(Kcross_CSR_List)){
		A1[[j]] <- Kcross_CSR_List[[j]][[i]]
	}
	A2 <- AverageCurve(DataList = A1)
	R_Matrix[i,] <- A2[[1]]
	Mu_Matrix[i,] <- A2[[2]]
}

#Check that R values are always identical
for(i in 1:ncol(R_Matrix)){
	if(length(unique(R_Matrix[,i]))>1){
		print("Problem in x axis indices!")
	}
}

#Generate empirical 95% CI
Average <- Up <- Down <- vector()
for(i in 1:ncol(Mu_Matrix)){
	Average[i] <-  mean(Mu_Matrix[,i])
	Down[i] <- quantile(Mu_Matrix[,i],probs=0.025)
	Up[i] <- quantile(Mu_Matrix[,i],probs=0.975)
}

#Plot
R <- R_Matrix[1,]
lines(x=R, y=Average - pi*R^2, lwd=2, col="blue")
xshade <- c(R,rev(R))
yshade <- c(Up-pi*R^2,rev(Down-pi*R^2))
polygon(xshade,yshade, border = NA, col=adjustcolor("blue", 0.4))

#Legend
legend(x="topleft", lwd=2, col=c("Red", "Blue"), legend=c("Observed", "CSR simulation"))
