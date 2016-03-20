### Simulate a spatial, inhomogenous Poisson point process (e.g., animal
### telemetry locations), and fit model using a
### Poisson GLM and a Poisson point process model. 

rm(list=ls())

##########################################################
### Libraries and subroutines
##########################################################

library(mvtnorm)
library(msm)
library(MCMCpack)
library(MASS)
library(sp)
library(rgeos)
library(raster)
# library(data.table)

get.cell.count <- function(s,S){  # tabulate number of locations s per grid cell in S
	# browser()
	idx <- extract(S,s,cellnumbers=TRUE)[,1]  # cell idx
	tab <- data.frame(table(idx))
	tab$cell.idx <- as.numeric(levels(tab$idx))[tab$idx]
	z <- numeric(ncell(S))
	z[tab$cell.idx] <- tab$Freq

	# library(data.table)  # tabulation with data.table faster for T > 1000
	# cell.idx <- data.table(cell.idx=cell.idx)
	# z <- cell.idx[,.N,by=cell.idx]
	z
}


##########################################################
### Simulate environmental covariate
##########################################################

# Spatial support of point process
S <- raster(xmn=0,xmx=1.5,ymn=0,ymx=1,res=0.02)

# Random resource covariate
n <- ncell(S)
xy <- xyFromCell(S, 1:n)
d <- as.matrix(dist(xy))		
phi <- 1
Sigma <- exp(-d/phi)  
Sigma <- t(chol(Sigma))
values(S) <- Sigma %*% rnorm(n)
plot(S)

# Smooth resource covariate
# w <- focalWeight(S,xres(S)*1,"circle")
# S <- focal(S,w,fun=mean,na.rm=TRUE,pad=TRUE)
S <- disaggregate(S,fact=2,method="bilinear")
plot(S)

values(S) <- scale(values(S))  # center and scale covariate

##########################################################
### Simulate animal locations
##########################################################

# Simulate resource selection
beta <- matrix(c(-1,1),2,1)  # resource selection covariates
X <- cbind(1,values(S))  # design matrix
qX <- ncol(X)
lambda <- exp(X%*%beta)

# Simulate observed locations via acceptance-rejection method
n <- 5000  # number of proposed locations
idx <- sample(1:ncell(S),n,replace=TRUE)  # uniform sample of cells
keep <- rbinom(n, 1, prob=lambda[idx]/max(lambda))  # keepers
s <- xyFromCell(S,idx[keep==1])  # USED LOCATIONS
n <- nrow(s)
plot(S)
points(s,col=1,pch=19,cex=0.5)

# Simulate observed locations by sampling cells with probability portional to lambda.
# Is this correct? Same as above?
# n <- 100  # number of locations to simulate
# s <- sample(1:ncell(rsf),n,prob=values(rsf))
# s <- xyFromCell(rsf,s)
# points(s,col=2,pch=19,cex=0.5)

##########################################################
### Fit model via Poisson point process model
##########################################################

source('~/Documents/git/PPM/spatial/spatial.ppp.mcmc.R')
priors <- list(sigma.beta=5)
tune <- list(beta=0.25)
start <- list(beta=rep(0,qX))
out1 <- spatial.ppp.mcmc(s,S,priors,start,tune,n.mcmc=5000)

matplot(out1$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)


##########################################################
### Fit model via Poisson GLM
##########################################################

source('~/Documents/git/GLM/poisson.glm.mcmc.R')
z <- get.cell.count(s,S)  # tabulate the number of locations per grid cell
priors <- list(sigma.beta=5)
tune <- list(beta=0.25)
start <- list(beta=coef(glm(z ~ 0+X, family=poisson())))
out2 <- poisson.glm.mcmc(z,X,priors,start,tune,n.mcmc=5000)

matplot(out2$beta, type="l", lty=1);abline(h=beta,col=1:qX,lty=3)


