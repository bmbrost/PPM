#
#
# Bayesian spatial inhomegenous Poisson point process model
#
# Function name: spatial.ppp.MCMC
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 19 MAR 2016
#
# Model statement:
#	s(t) ~ exp(x(s(t))%*%beta)/\int(exp(x(s)%*%beta)ds)
#	beta ~ N(0,sigma.beta^2*I)
#
# Reference:
#
# Required R packages: mvtnorm (if using multivariate normal prior on beta)
#
# Inputs:
#
# s - Tx2 matrix of observed locations (x,y)
# S - rasterized spatial support of the point process with values 
#	corresponding to the covariates for which inference is desired. A
#	raster stack is used for multiple covariates, with one raster layer
# 	per covariate.
# priors - list of priors containing the following elements:
#	1. sigma.beta - Standard deviation of normal prior on beta
# start - list of starting values containing the following elements:
#	1. beta - vector of starting values for resource selection coefficients
# tune - list of tuning parameters containing the following elements:
#	1. beta - Tuning parameter for Metropolis-Hasting update on beta
# adapt - Switch to enable adapative tuning (TRUE/FALSE)
# n.mcmc - Number of desired MCMC iterations
#
#

spatial.ppp.mcmc <- function(s,S,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	# library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	
	###
	###  Setup variable, starting values, and priors
	###
# browser()
	X <- cbind(1,values(S))  # design matrix
	qX <- ncol(X)  # number of betas
	cell.idx <- extract(S,s,cellnumbers=TRUE)[,1]	
	X.used <- X[cell.idx,]
	
	beta <- start$beta
	int <- log(sum(exp(X%*%beta)))  # normalizing constant in denominator of ppp likelihood

	lambda <- exp(X%*%beta)  # intensity of Poisson process

	mu.beta <- rep(0,qX)  # mean of normal prior on beta
	sigma.beta <- priors$sigma.beta
	
	# Sigma.beta <- solve(crossprod(X))
	# c <- sigma.beta/min(diag(Sigma.beta))
	# Sigma.beta <- c*Sigma.beta  # variance-covariance matrix of
		# # normal prior on beta

	###
	### Create receptacles for output
	###
  
	beta.save <- matrix(0,n.mcmc,qX)  # rsf coefficients
	
	keep <- list(beta=0)
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){

		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune$beta <- get.tune(tune$beta,keep.tmp$beta,k)
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update beta
		### 
		
		beta.star <- rnorm(qX,beta,tune$beta)  # proposal for beta
		int.star <- log(sum(exp(X%*%beta.star)))  # normalized constant in ppp likelihood
  		mh.0 <- sum(X.used%*%beta-int)+sum(dnorm(beta,mu.beta,sigma.beta,log=TRUE))
		mh.star <- sum(X.used%*%beta.star-int.star)+sum(dnorm(beta.star,mu.beta,sigma.beta,log=TRUE))
		if(exp(mh.star-mh.0)>runif(1)){
			beta <- beta.star
			int <- int.star
			keep$beta <- keep$beta+1
			keep.tmp$beta <- keep.tmp$beta+1
		}

		###
		###  Save samples 
	    ###

		beta.save[k,] <- beta
	}
	
	keep$beta <- keep$beta/n.mcmc
	cat(paste("\nbeta acceptance rate:",round(keep$beta,2))) 

	###
	### Write output
	###

	list(beta=beta.save,keep=keep,
		s=s,S=S,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}