require(stringr)
require(reshape2)
require(ggplot2)
require(gridExtra)
require(rjags)
require(foreach)
require(doParallel)
require(mcmcplots)

#read-in time series
y <- read.csv("...", row.names=1)

#------------------------ environmental covariates --------------------
temp <- c(12.38,12.66,16.84,17.77,20.3,21.1,20.41,17.64,17.29,13.81,12.92,12.32,12.39,13.65,15.73,19.24,23,24.11,24.43,17.61,15.76,14.45,12.89,12.36,12.25,17.6,21,22.4,22.88,23,20.64,16.76,15.46,14.39,12.24,13.18)
sal <- c(38.18,38.15,37.6,38.17,38.13,38.17,38.14,38.07,38.16,37.98,37.93,37.92,37.61,37.81,37.59,37.49,37.76,37.75,38.03,37.86,38.17,38.17,38.00,37.82,37.82,NA,37.7,37.83,37.82,NA,38.09,37.21,38.05,38.09,38.24,38.28)
PO4 <- c(0.049,0.044,0.054,0.032,0.045,0.029,0.096,0.102,0.117,0.104,0.108,0.107,0.131,0.112,0.095,0.106,0.078,0.151,0.118,0.139,0.076,0.184,0.112,0.107,0.110,0.096,0.182,0.076,0.105,0.092,0.083,0.101,0.074,0.078,0.130,0.103)
NH4 <- c(0.327,0.363,0.233,0.472,0.449,0.099,0.166,0.156,0.101,0.083,0.067,0.453,0.145,0.144,0.044,0.268,0.055,0.122,0.201,0.225,0.144,0.116,0.616,0.559,0.408,0.037,0.629,0.029,0.687,0.095,0.026,0.052,0.027,1.389,0.390,0.666)
NO2 <- c(0.033,0.031,0.015,0.007,0.011,0.028,0.010,0.034,0.045,0.172,0.231,0.239,0.107,0.084,0.064,0.062,0.029,0.065,0.025,0.055,0.150,0.251,0.275,0.238,0.283,0.016,0.102,0.010,0.038,0.004,0.006,0.060,0.106,0.263,0.222,0.259)
bact <- c(7.75E+05,1.11E+06,9.39E+05,8.37E+05,9.53E+05,6.77E+05,7.05E+05,1.14E+06,9.26E+05,6.93E+05,6.97E+05,4.12E+05,6.90E+05,9.10E+05,9.12E+05,1.06E+06,9.05E+05,7.64E+05,6.11E+05,1.08E+06,4.83E+05,7.68E+05,3.15E+05,6.30E+05,NA,8.75E+05,NA,1.07E+06,NA,8.22E+05,1.09E+06,1.33E+06,8.36E+05,5.91E+05,4.37E+05,4.92E+05)
chl <- c(NA,1.16,0.58,0.12,0.13,0.34,0.34,0.74,0.76,1.13,1,0.96,1.95,0.49,0.44,0.56,0.25,0.44,0.31,0.78,0.53,0.51,0.81,0.96,1.05,0.49,0.54,0.2,0.6,0.32,0.37,2.88,0.89,0.69,0.69,0.67)
#--------------------------------------------------------------------------------------


# data list for rjags model
dat.list <- function(time_series) {
  time_series_t <- as.matrix(t(time_series))
  time_series_Nt <- as.vector(apply(time_series_t, 1, sum))
  
  # Add month as a covariate
  sp_date <- rownames(time_series_t)
  months_raw <- str_split_fixed(sp_date, "_", 3)[,3]
  months <- gsub(x=months_raw, pattern=13, replacement=01)
  months <- gsub(x=months, pattern=14, replacement=02)
  months <- gsub(x=months, pattern=15, replacement=03)
  months <- as.numeric(months)
  
  dataList <- list( y = structure(time_series_t, dim=c( nrow(time_series_t), ncol(time_series_t) ) ) ,
                    # N should sum over time t NOT species i
                    N = time_series_Nt ,
                    NMonths = nrow(time_series_t) ,
                    NSpecies = ncol(time_series_t) ,
                    NLatent = 2 ,
                    # covariates
                    months = months ,
                    temp = temp ,
                    sal = sal ,
                    bact = bact ,
                    chl = chl ,
                    PO4 = PO4 ,
                    NH4 = NH4 ,
                    NO2 = NO2
  )
  return(dataList)
}

#inits list for multi-threading
inits.list2 <- function(dataList) {
  return( c( 
    inits.list(dataList), 
    .RNG.name = "lecuyer::RngStream" , 
    .RNG.seed = randomNumbers(n=1, min=1, max=1e+06) 
  ) )
}

# new inits for rjags model
inits.list <- function(dataList) {
 initsList <- list(
   r=r.new.inits, #from previous runs
   k=k.new.inits, #from previous runs
   alpha.tmp=alpha.tmp.new.inits, #from previous runs
   # the covariates need to have initis otherwise divide-by-zero error in variance partitioning
   tempCov=rep(0.001,nrow(y)),
   salCov=rep(0.001,nrow(y)),
   chlCov=rep(0.001,nrow(y)),
   bactCov=rep(0.001,nrow(y)),
   PO4Cov=rep(0.001,nrow(y)),
   NO2Cov=rep(0.001,nrow(y)),
   NH4Cov=rep(0.001,nrow(y))
  )
 return(initsList)
}

#generate n random number from a sequence of N numbers
randomNumbers <- function(n, min, max) {
  vectorOfnumbers <- seq(min:max)
  shuffleList <- sample(vectorOfnumbers)
  n_number <- shuffleList[n]
  return(n_number)
}

mcmc.combine <- function( ... ) {
  return( as.mcmc.list( sapply( list( ... ), mcmc ) ) )
}

multi_species_model_new <- "
  model {

  ### Multi-species model ###

    # Standardization and centering of covariates
    
    mTemp <- mean( temp[] )
    sdTemp <- sd( temp[] )
    mSal <- mean( sal[] )
    sdSal <- sd( sal[] )
    mChl <- mean( chl[] )
    sdChl <- sd( chl[] )
    mBact <- mean( bact[] )
    sdBact <- sd( bact[] )
    mPO4 <- mean( PO4[] )
    sdPO4 <- sd( PO4[] )
    mNH4 <- mean( NH4[] )
    sdNH4 <- sd( NH4[] )
    mNO2 <- mean( NO2[] )
    sdNO2 <- sd( NO2[] )

    muN[1] <- 0
    for( t in 2:NMonths ) {
      for ( i in 1:NSpecies ) {

        # Stochastic component
        mu[t,i] ~ dnorm( mu.star[t,i], tau[i] )
        y[t,i] ~ dpois( lambda[t,i] )

        # Linear predictor
        log( lambda[t,i] ) <- mu[t,i] + log( N[t] ) + muN[t]

        mu.star[t,i] <- mu[t-1, i] + r[i] * ( 1 - inprod( alpha[i, ], mu[t-1, ] ) / k[i] ) +

        # Covariates standardized to unit variance
        monthCov[i, months[t]] + 
        tempCov[i]*((temp[t]-mTemp)/sdTemp) + 
        salCov[i]*((sal[t]-mSal)/sdSal) + 
        chlCov[i]*((chl[t]-mChl)/sdChl) + 
        bactCov[i]*((bact[t]-mBact)/sdBact) + 
        PO4Cov[i]*((PO4[t]-mPO4)/sdPO4) + 
        NH4Cov[i]*((NH4[t]-mNH4)/sdNH4) + 
        NO2Cov[i]*((NO2[t]-mNO2)/sdNO2) +

        # LVs
        inprod(LV[t,], Loading[i,])
      
      }
      muN[t] ~ dnorm(0,0.01)
    }

    # Priors LVs

    for( t in 1:NMonths ) {
		  for( l in 1:NLatent ) {
			  LV[t,l] ~ dnorm(0,1)  #LVs
		  }
	  }

    for( i in 1:NSpecies ) {
      for( l in 1:NLatent ) {
        Loading[i,l] ~ dnorm(0,1)   #Loadings
    	  lamSp[i,l] <- Loading[i,l] * sd(LV[,l])
      }
    }

    # Covariance & correlation matrices

    for( i in 1:NSpecies ) {
      for( j in 1:NSpecies ) {
        Cov[i,j] <- inprod( lamSp[i,], lamSp[j,] )
        Corr[i,j] <- Cov[i,j] / sqrt( ( Cov[i,i]+1 / tau[i] ) * ( Cov[j,j]+1 / tau[j] ) )   		  
 		  }
    }

    for( i in 1:NSpecies ) {

      # Stationary variance of mu for variance partitioning calculations

      mu.sd[i] <- sd( mu[,i] )
      svar[i] <- mu.sd[i] * mu.sd[i]

      # Interaction coefficient, alpha and the indicator variable g.alpha

      alpha[i,i] <- 1                             
      g.alpha[i,i] <- 1
      alpha.square[i,i] <- 1

      for( j in (i + 1):NSpecies ) {
        
        # Used to calculate total variance for each species

        # Gibbs variable selection (GVS)

        alpha.tmp[i,j] ~ dnorm( 0, alpha.tau[i,j] )
        alpha.tmp[j,i] ~ dnorm( 0, alpha.tau[j,i] )

        alpha[i,j] <- alpha.tmp[i,j] * g.alpha[i,j]
        alpha[j,i] <- alpha.tmp[j,i] * g.alpha[j,i]

        alpha.square[i,j] <- pow( alpha[i,j], 2 )
        alpha.square[j,i] <- pow( alpha[j,i], 2 )

        g.alpha[i,j] ~ dbern( 0.1 )
        g.alpha[j,i] ~ dbern( 0.1 )

        alpha.tau[i,j] <- 1/alpha.var[i,j]
        alpha.tau[j,i] <- 1/alpha.var[j,i]

        alpha.var[i,j] <- (1 - g.alpha[i,j]) * 1 + g.alpha[i,j] * 1
        alpha.var[j,i] <- (1 - g.alpha[j,i]) * 1 + g.alpha[j,i] * 1

      }
    }

    # Variance components

    for( i in 1:NSpecies ) {

    # Variance for each covariate 

      bactVar[i] <- bactCov[i] * bactCov[i]
      tempVar[i] <- tempCov[i] * tempCov[i]
      salVar[i] <- salCov[i] * salCov[i]
      chlVar[i] <- chlCov[i] * chlCov[i]
      NH4Var[i] <- NH4Cov[i] * NH4Cov[i]
      NO2Var[i] <- NO2Cov[i] * NO2Cov[i]
      PO4Var[i] <- PO4Cov[i] * PO4Cov[i]
 
     # Unexplained variance for each species == Environmental forcing
       
       n.sigma2[i] <- Cov[i,i]
 
     # Variance explained by the environmental covariates
       
      coVar[i] <- NO2Var[i] + NH4Var[i] + PO4Var[i] + tempVar[i] + salVar[i] + chlVar[i] + bactVar[i]
 
      # Total variance for every species
  
      totVar[i] <- pow( ( r[i] / k[i]), 2 ) * inprod(alpha.square[i,], svar[]) + n.sigma2[i] + coVar[i]
  
      # Intra/inter-spec. interactions and env. stochasticity
  
      env[i] <- n.sigma2[i] + coVar[i]
      inter[i] <- totInter[i] - intra[i]
      totInter[i] <- pow( ( r[i] / k[i] ), 2 ) * inprod( alpha.square[i,], svar[] )  # Same as totVar except the additional variance terms
      intra[i] <- pow( ( r[i] / k[i] ), 2 ) * svar[i]  # the formula should include alpha.square[i,i] == 1 but since it zero we skip it in the calc

     # Proportions of total variance explained 
 
      varIntra[i] <- intra[i] / totVar[i]
      varInter[i] <- inter[i] / totVar[i]
      varEnv[i] <- n.sigma2[i] / totVar[i]
      varCov[i] <- coVar[i] / totVar[i]
      varEC[i] <- env[i] / totVar[i]
 
     # Prop of env. variance explained by the covariates
 
      bactVarEnv[i] <- bactVar[i] /  env[i]
      tempVarEnv[i] <- tempVar[i] /  env[i]
      salVarEnv[i] <-  salVar[i] /  env[i]
      chlVarEnv[i] <-  chlVar[i] /  env[i]
      NH4VarEnv[i] <-  NH4Var[i] /  env[i]
      NO2VarEnv[i] <-  NO2Var[i] /  env[i]
      PO4VarEnv[i] <-  PO4Var[i] /  env[i]
      propVarEnv[i] <- coVar[i] /  env[i]
 
      bactProp[i] <- bactVar[i] / coVar[i]
      tempProp[i] <- tempVar[i] / coVar[i]
      salProp[i] <-  salVar[i] / coVar[i]
      chlProp[i] <-  chlVar[i] / coVar[i]
      NH4Prop[i] <-  NH4Var[i] / coVar[i]
      NO2Prop[i] <-  NO2Var[i] / coVar[i]
      PO4Prop[i] <-  PO4Var[i] / coVar[i]

    }

    # Set the parameter model (priors) for stoachistic nodes

    for( i in 1:NSpecies ) {

      y[1,i] ~ dpois( lambda[1,i] )
      mu[1,i] ~ dnorm( 0, 0.001 )
      log( lambda[1,i] ) <- mu[1, i] + log( N[1] ) + muN[1]

    # Prior for k on log scale. Large variance = Small precision
      k[i] ~ dexp( 1 )
      r[i] ~ dnorm( 0, 0.1 )

    # Correlation matrix & std deviation
      tau[i] <- pow( sd.corr[i], -2 )
      sd.corr[i] ~ dunif( 0, 1000 )
    }

    # Priors for covariates

    for( i in 1:NSpecies ) {
      for( m in 1:12 ) {
        monthCov[i, m] ~ dnorm( 0, monthCov.tau )
      }
    }
    monthCov.tau <- pow( monthCov.sigma, -2 )
    monthCov.sigma ~ dunif( 0, 100 )

    for( i in 1:NSpecies ) {
      tempCov[i] ~ dnorm( 0, tempCov.tau )
      salCov[i] ~ dnorm( 0, salCov.tau )
      chlCov[i] ~ dnorm( 0, chlCov.tau )
      bactCov[i] ~ dnorm( 0, bactCov.tau )
      PO4Cov[i] ~ dnorm( 0, PO4Cov.tau )
      NH4Cov[i] ~ dnorm( 0, NH4Cov.tau )
      NO2Cov[i] ~ dnorm( 0, NO2Cov.tau )
    }
     
    tempCov.tau <- pow( tempCov.sigma, -2 )
    tempCov.sigma ~ dunif( 0, 100 )
 
    salCov.tau <- pow( salCov.sigma, -2 )
    salCov.sigma ~ dunif( 0, 100 )

    chlCov.tau <- pow( chlCov.sigma, -2 )
    chlCov.sigma ~ dunif( 0, 100 )

    bactCov.tau <- pow( bactCov.sigma, -2 )
    bactCov.sigma ~ dunif( 0, 100 )

    PO4Cov.tau <- pow( PO4Cov.sigma, -2 )
    PO4Cov.sigma ~ dunif( 0, 100 )

    NH4Cov.tau <- pow( NH4Cov.sigma, -2 )
    NH4Cov.sigma ~ dunif( 0, 100 )
 
    NO2Cov.tau <- pow( NO2Cov.sigma, -2 )
    NO2Cov.sigma ~ dunif( 0, 10 )

  # Assign a distribution to the covariates (or other vectors) to handle NAs
  
    for( t in 1:NMonths ) {
      temp[t] ~ dnorm( mu.temp, tau.temp )
      sal[t] ~ dnorm( mu.sal, tau.sal )
      bact[t] ~ dnorm( mu.bact, tau.bact )
      chl[t] ~ dnorm( mu.chl, tau.chl )
      PO4[t] ~ dnorm( mu.PO4, tau.PO4 )
      NH4[t] ~ dnorm( mu.NH4, tau.NH4 )
      NO2[t] ~ dnorm( mu.NO2, tau.NO2 ) 
    }
  
    mu.temp ~ dnorm( 0, 0.0001 )
    tau.temp <- pow( sigma.temp, -2 )
    sigma.temp ~ dunif( 0, 100 )
 
    mu.sal ~ dnorm( 0, 0.0001 )
    tau.sal <- pow( sigma.sal, -2 )
    sigma.sal ~ dunif( 0, 100 )
 
    mu.bact ~ dnorm( 0, 0.0001 )
    tau.bact <- pow( sigma.bact, -2 )
    sigma.bact ~ dunif( 0, 100 )
     
    mu.chl ~ dnorm( 0, 0.0001 )
    tau.chl <- pow( sigma.chl, -2 )
    sigma.chl ~ dunif( 0, 100 )

    mu.PO4 ~ dnorm( 0, 0.0001 )
    tau.PO4 <- pow( sigma.PO4, -2 )
    sigma.PO4 ~ dunif( 0, 100 )

    mu.NH4 ~ dnorm( 0, 0.0001 )
    tau.NH4 <- pow( sigma.NH4, -2 )
    sigma.NH4 ~ dunif( 0, 100 )

    mu.NO2 ~ dnorm( 0, 0.0001 )
    tau.NO2 <- pow( sigma.NO2, -2 )
    sigma.NO2 ~ dunif( 0, 100 )

  }  
"

#------------------------------------------------------------------------------
# RUN THE CHAINS

# The parameters to be monitored
# Output everything!
parameters <- c( "r","k","g.alpha","alpha","mu","mu.star","n.sigma2","tau","sd.corr","LatentSTAR","lamSpSTAR","lamSp","Cov","Corr","totVar","totInter","env","inter","intra","varEnv","varInter","varIntra",
                   "monthCov","bactCov","tempCov","salCov","chlCov","NH4Cov","NO2Cov","PO4Cov",
                   "bactVarEnv","tempVarEnv","salVarEnv","chlVarEnv","NH4VarEnv","NO2VarEnv","PO4VarEnv",
                   "bactProp","tempProp","salProp","chlProp","NH4Prop","NO2Prop","PO4Prop",
                   "bactVar","tempVar","salVar","chlVar","NH4Var","NO2Var","PO4Var",
                   "propVarEnv","coVar","varEC","varCov","varEnvCov")


registerDoParallel(10)

# Number of steps to "tune" the samplers
adaptSteps = 5000              
nChains = 10

numSavedSteps=100000
thinSteps=10
burnInSteps = 50000

nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

run.time <- system.time( {
  
  codaSample <- foreach( i = 1:getDoParWorkers( ) , 
                             .inorder=FALSE , 
                             .combine = "mcmc.combine" ,
                             .multicombine=TRUE ,
                             .packages= "rjags",
                             .verbose=TRUE
  ) %dopar% {
    load.module( "lecuyer" )
    
    # Create, initialize, and adapt the model:
    jagsModel <- jags.model( textConnection(multi_species_model_new) , 
                             data=dat.list(y),
                             inits=inits.list2(dat.list(y)),
                             n.adapt=adaptSteps )
    
    # Burn-in:
    update( jagsModel , n.iter=burnInSteps )
    
    # The saved MCMC chain:
    result <- coda.samples( jagsModel , 
                            variable.names=parameters ,
                            n.iter=nIter , 
                            thin=thinSteps)
    return( result )
  }
  
} )[3]

#name output
save.image("...")
q("no")
