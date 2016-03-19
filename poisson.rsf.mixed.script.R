### Simulate locations from multiple animals and environmental covariates
### for the purposes of examining resource selection. Model
### fitting proceeds using a Poisson GLMM

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


##########################################################
### Population-level process model parameters
##########################################################

mu.beta <- matrix(c(-0.5,1),,1)  # mean of betas
qX <- nrow(mu.beta)
rho <- -0.15  # correlation between betas
Sigma <- diag(qX)*0.5  # variance-covariance of betas
Sigma[1,2] <- Sigma[2,1] <- Sigma[1,1]*Sigma[2,2]*rho


##########################################################
### Simulate group-level process model parameters
##########################################################

beta <- t(rmvnorm(J,mu.beta,Sigma))  # betas for each group
plot(t(beta))


##########################################################
### Simulate animal locations
##########################################################

J <- 10  # number of animals
X <- cbind(1,scale(values(S)))  # design matrix
n <- 500  # number of proposed locations

s <- matrix(ncol=2)  # animal locations
g <- numeric()  # grouping variable

for(i in 1:J){  # loop through animals
	idx <- sample(1:ncell(S),n,replace=TRUE)  # uniform sample of cells
	lambda <- exp(X[idx,]%*%beta[,i])	
	keep <- rbinom(n, 1, prob=lambda/max(lambda))  # keepers
	s <- rbind(s,xyFromCell(S,idx[keep==1]))  # USED LOCATIONS
	g <- c(g,rep(i,sum(keep)))
}
s <- s[-1,]

plot(S)
points(s)
table(g)


##########################################################
### Fit model
##########################################################

source('~/Documents/git/GLMM/poisson.glmm.mcmc.R')

z <- c(sapply(1:J,function(x) get.cell.count(s[g==x,],S)))
g.tmp <- rep(1:J,each=ncell(S))
X.tmp <- do.call("rbind",rep(list(X),J))
priors <- list(sigma.beta=5,S0=diag(qX),nu=qX+1)
tune <- list(beta=rep(2.0,J))
start <- list(beta=matrix(0,qX,J),mu.beta=mu.beta,Sigma=Sigma)
out1 <- poisson.glmm.mcmc(z,X.tmp,g.tmp,priors,start,tune,adapt=TRUE,n.mcmc=10000)

# Examine estimates for beta_j
g.idx <- 10  # group idx for plotting beta_j
matplot(out1$beta[,,g.idx],type="l",lty=1);abline(h=beta[,g.idx],col=1:qX,lty=2)

# Examine estimates for mu.beta
matplot(out1$mu.beta,type="l");abline(h=mu.beta,col=1:qX,lty=2)

# Examine estimates for Lambda
matplot(cbind(out1$Sigma[1,1,],out1$Sigma[1,2,]),type="l")
abline(h=c(Sigma[1,1],Sigma[1,2]),lty=2,col=1:qX)
