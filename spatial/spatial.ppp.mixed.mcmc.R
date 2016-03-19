#
#
# Bayesian spatial, inhomogenous Poisson point process model
#
# Function name: spatial.ppp.mixed.mcmc
#
# Author: Brian M. Brost
# Contact: bmbrost@gmail.com
#
# Last updated: 19 MAR 2016
#
# Model statement:
#	s_ij ~ exp(x(s_ij)%*%beta_j)/\int(exp(x(s_j)%*%beta_j)d(s_j))
# 	beta_j ~ N(mu_beta,Sigma.beta)
#	mu_beta ~ N(0,sigma.beta^2*I)
#	Sigma.beta ~ Wish(S_0,nu)
#
# Reference:
#
# Required R packages: mvtnorm
#
# Inputs:
#

# s - Tx2 matrix of observed locations (x,y)
# S - rasterized spatial support of the point process with values 
#	corresponding to the covariates for which inference is desired. A
#	raster stack is used for multiple covariates, with one raster layer
# 	per covariate.


# z - vector of length n containing the count of observations corresponding to 
# 	each row in the design matrix X. Note that the value of z[1] corresponds to
#	X[1,], z[2] corresponds to X[2,], etc.
# X - design matrix of dimension n x qX containing covariates (plus
#	intercept) for which inference is desired
# g - variable that defines groups of observations in z
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

spatial.ppp.mixed.mcmc <- function(s,S,g,priors,start,tune,adapt=TRUE,n.mcmc=1000){

	###
	###  Libraries and Subroutines
	###

	library(mvtnorm)

	get.tune <- function(tune,keep,k,target=0.44){  # adaptive tuning
		# a <- min(0.01,1/sqrt(k))
		a <- min(0.025,1/sqrt(k))
		exp(ifelse(keep<target,log(tune)-a,log(tune)+a))
	}
	
	
	###
	###  Setup variable
	###
# browser()
	J <- length(unique(g))  # number of groups
	g <- as.numeric(g)
	g.idx <- sapply(sort(unique(g)),function(x) which(g==x),simplify=FALSE)
		# indexes of observations in y by group
	n.j <- unlist(lapply(g.idx,length))  # number of observations per group
	n <- length(s)  # total number of observations
	X <- cbind(1,values(S))  # design matrix
	qX <- ncol(X)

	cell.idx <- extract(S,s,cellnumbers=TRUE)[,1]
	X.used <- X[cell.idx,]
		
	###
	###  Starting values, and priors
	###

	beta <- matrix(start$beta,qX,J)	
	int <- sapply(1:J,function(x) log(sum(exp(X%*%beta[,x]))))
		  # normalizing constant in denominator of ppp likelihood
	cp <- lapply(1:J,function(x) solve(crossprod(X.used[g.idx[[x]],])))  # cross-product for propsal dist.

	mu.beta <- matrix(start$mu.beta,qX,1)
	Sigma <- start$Sigma
	Sigma.inv <- solve(Sigma)
	
	Sigma.beta <- priors$sigma.beta^2*diag(qX)
	Sigma.beta.inv <- solve(Sigma.beta)
	S0 <- priors$S0
	nu <- priors$nu

	###
	### Create receptacles for output
	###
  
	beta.save <- array(0,dim=c(n.mcmc,qX,J))
	mu.beta.save <- matrix(0,n.mcmc,qX)
	Sigma.save <- array(0,dim=c(qX,qX,n.mcmc))
	
	keep <- list(beta=rep(0,J))
	keep.tmp <- keep  # track MH accpetance rate for adaptive tuning
	Tb <- 50  # frequency of adaptive tuning

	###
	###  Begin MCMC loop
	###

	for(k in 1:n.mcmc){
# browser()
		if(k%%1000==0) cat(k,"");flush.console()	

		if(adapt==TRUE & k%%Tb==0) {  # Adaptive tuning
			# browser()
			keep.tmp <- lapply(keep.tmp,function(x) x/Tb)
			tune$beta <- sapply(1:J,function(x) get.tune(tune$beta[i],keep.tmp$beta[i],k))
			keep.tmp <- lapply(keep.tmp,function(x) x*0)
	   	} 	

		###
		### Update beta_j
		### 
		
		for(i in 1:J){
			idx <- g.idx[[i]]
			# beta.star <- rnorm(qX,beta[,i],tune$beta)
			# tune.tmp <- tune$beta[i]*solve(crossprod(X.used[idx,]))
			tune.tmp <- tune$beta[i]*cp[[i]]
			beta.star <- c(rmvnorm(1,beta[,i],tune.tmp))
			int.star <- log(sum(exp(X%*%beta.star))) # normalizing constant in ppp likelihood
	  		mh.0 <- sum(X.used[idx,]%*%beta[,i]-int[i])+
				sum(dmvnorm(beta[,i],mu.beta,Sigma,log=TRUE))
			mh.star <- sum(X.used[idx,]%*%beta.star-int.star)+
				sum(dmvnorm(beta.star,mu.beta,Sigma,log=TRUE))
			if(exp(mh.star-mh.0)>runif(1)){
				beta[,i] <- beta.star
				int[i] <- int.star
				keep$beta[i] <- keep$beta[i]+1
				keep.tmp$beta[i] <- keep.tmp$beta[i]+1
			}
		}
		
	  	###
	  	### Sample mu_beta
	  	###

		beta.mean <- apply(beta,1,sum)
		A.inv <- solve(J*Sigma.inv+Sigma.beta.inv)
		b <- Sigma.inv%*%beta.mean
	    mu.beta <- A.inv%*%b+t(chol(A.inv))%*%matrix(rnorm(qX),qX,1)

	  	###
	  	### Sample Sigma
	  	###
# browser()		
	  	Sn <- S0+crossprod(t(beta)-matrix(mu.beta,J,qX,byrow=TRUE))
		Sigma <- solve(rWishart(1,nu+J,solve(Sn))[,,1])
		Sigma.inv <- solve(Sigma)

		###
		###  Save samples 
	    ###

	  	beta.save[k,,] <- beta
		mu.beta.save[k,] <- mu.beta
		Sigma.save[,,k] <- Sigma
	}
	cat("\n")
	keep$beta <- keep$beta/n.mcmc
	cat("beta acceptance rate:",round(keep$beta,2)) 

	###
	### Write output
	###

	list(beta=beta.save,mu.beta=mu.beta.save,Sigma=Sigma.save,keep=keep,
		s=s,S=S,g,priors=priors,start=start,tune=tune,n.mcmc=n.mcmc)
}