#######################################################################################
##Spatial analysis of actin structures (Dataset2 = Flat)
##Hugo Lachuer (last update 15/05/2024)
##R version 4.2.2 (2022-10-31 ucrt) -- "Innocent and Trusting"
#######################################################################################

#######################################################################################
##1/ Functions and packages

#Packages
library(raster)
library(spatstat)
library(readxl)

#Import functions
setwd("")	#Path for Functions.R
source("Functions.R")

#######################################################################################
##2/ Import data

#Parameters
Nsim <- 100	#Number of Monte-Carlo simulations

#Get images dimensions
Dimensions <- read.delim("Dataset2/ImageDimensions.txt")	#Path


#List of files
setwd("/Dataset2")	#Path for data (Dataset2 flat)
V1 <- list.files(, pattern = ".tif")
V1 <- substr(V1, 1, nchar(V1)-9)	#List of files names
K_List <- list()	#List of 2D matrices (x,y)
K_CSR_List <- list()	#List of a list of 2D matrices (x,y)
PPP_list <- list() #list of point pattern

#Offset for FloodFill if center of mass is in a 0 area
offsetX <- c(0, 350, rep(0,12))

#Save plot
pdf("Individual plot DataSet2.pdf")

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
	Actin <- read_excel(paste0(V1[i],".xlsx"))
	Actin <- Actin[,-1]
	Actin <- as.matrix(Actin)
	Actin <- Actin*ConversionFactor	#Conversion

	#Plot
	#plot(raster(Mask), col=c("white","black"))
	#points(x=Actin[,1]/ncol(Mask), y=1-Actin[,2]/nrow(Mask), pch=16, col="red")

	#Final Mask
	plot(raster(Mask), col=c("white","black"), main=V1[i])
	points(x=Actin[,1]/ncol(Mask), y=1-Actin[,2]/nrow(Mask), pch=16, col="red")


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
		ActinCSR <- CSRCell(n=nrow(Actin), Mask=Mask)
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

}

dev.off()

#Saving
setwd("/Dataset2")	#Path to "Dataset2 (Flat)" folder
save(K_List, file="K_List.RData")
save(K_CSR_List, file="K_CSR_List.RData")
save(PPP_list, file="PPP_List.RData")


#######################################################################################
##5/ Average Ripley's K function


setwd("/Dataset2")	#Path to "Dataset2 (Flat)" folder
load(file="K_List.RData")
load(file="K_CSR_List.RData")

#Average Ripley K function
A1 <- AverageCurve(DataList = K_List)
R <- A1[[1]]
Mu <- A1[[2]]
SEM <- A1[[3]]

plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K - pir²", main="Average Ripley K function", xlim=c(0,15), ylim=c(-100,800))
#plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K function - pir²", main="Average Ripley K function", xlim=c(0,3), ylim=c(-20,100))
#plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K function - pir²", main="Average Ripley K function", xlim=c(0,0.5), ylim=c(-1,1))
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
##6/ Pillar vs non-pillar


######################################
##Plot pillar data (Dataset 1)

#Load data (code can be started here)
setwd("/Dataset1")		#Path to "Dataset1 (Pillar)" folder
load(file="K_List.RData")
load(file="K_CSR_List.RData")

#Average Ripley K function
A1 <- AverageCurve(DataList = K_List)
R <- A1[[1]]
Mu <- A1[[2]]
SEM <- A1[[3]]

#plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K function - pir²", main="Average Ripley K function", xlim=c(0,3), ylim=c(-20,100))
plot(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="red", xlab="Distance (µm)", ylab="Ripley's K function - pir²", main="Average Ripley K function", xlim=c(0,15), ylim=c(-100,800))
xshade <- c(R,rev(R))
yshade <- c(Mu-pi*R^2+SEM,rev(Mu-pi*R^2-SEM))
polygon(xshade,yshade, border = NA, col=adjustcolor("red", 0.4))
lines(x=R, y=rep(0, length=length(R)), lwd=2, col="black", lty=2)

######################################
##Plot flat data (Dataset 2)

#Load data (code can be started here)
setwd("/Dataset2")
load(file="K_List.RData")
load(file="K_CSR_List.RData")

#Average Ripley K function
A1 <- AverageCurve(DataList = K_List)
R <- A1[[1]]
Mu <- A1[[2]]
SEM <- A1[[3]]

lines(x=R, y=Mu - pi*R^2, lwd=2, type="l", col="blue")
xshade <- c(R,rev(R))
yshade <- c(Mu-pi*R^2+SEM,rev(Mu-pi*R^2-SEM))
polygon(xshade,yshade, border = NA, col=adjustcolor("blue", 0.4))


#Legend
legend(x="topleft", lwd=2, col=c("Red", "Blue"), legend=c("Pillar", "No pillar"))


######################################
##Statistical test

#Load dataset 1
setwd("/Dataset1")		#Path to "Dataset1 (Pillar)" folder
load("PPP_List.RData")
PPP_DataSet1 <- PPP_list

#Load dataset 2
setwd("/Dataset2")		#Path to "Dataset2 (Flat)" folder
load("PPP_List.RData")
PPP_DataSet2 <- PPP_list

#Merge
X <- list()
X[[1]] <- PPP_DataSet1
X[[2]] <- PPP_DataSet2

#Test
studpermu.test(X, summaryfunction = Kest, rinterval = c(0,3), nperm = 999)	#limited to 3µm

studpermu.test(X, summaryfunction = Kest, nperm = 999)	#Unlimited