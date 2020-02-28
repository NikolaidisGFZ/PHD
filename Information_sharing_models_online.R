
#####  These synthesis models have been developed by Georgios Nikolaidis #############################################
#####  for the purposes of a simulation experiment that forms part of my PhD Thesis ##################################
#####  The models below build on the code developed by Dias et. al for the TSD series ################################
#####  The code has been reduced to the pairwise meta-analysis case without adjustments for multi-arm trials #########
#####  All the models below consider a RANDOM-EFFECTS approach, but can easily be adapted for a Fixed-effect #########


#### Splitting model (No sharing) ####
# THis model does not impose any information sharing and is only shown here for comparison purposes
# Note that both sources of evidence are analysed in the same model, but in separate loops. There are no common parameters and hence no informatio-sharing

Splitting_T2_model <- function(){
  
  for(i in 1:ns.dir){                       # LOOP THROUGH DIRECT STUDIES
    
    delta[i,1] <- 0                         # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                  # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                    # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])          # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]   # model for linear predictor
      
    }                                       # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2],prec.dir)
    
  }                                         # DIRECT studies Loop closes   
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){       # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                         # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                  # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                    # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])          # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]   # model for linear predictor
      
    }                                       # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2],prec.indir)
    
  }                                         # Indirect studies Loop closes  
  
  d.dir[1]<-0                               # treatment effect is zero for reference treatment in Direct
  
  for (k in 2:nt){  d.indir[k] ~ dnorm(0,.0001) }   # vague priors for treatment effects
  
  tau.dir ~ dunif(0.000001, 5)              # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)               # between-trial precision = (1/between-trial variance)
  tau.indir ~ dunif(0.000001, 5)            # vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)           # between-trial precision = (1/between-trial variance)
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, d.indir, tau.dir, tau.ind, mu


#### Lumping only on relative effects ####

# This model defines separate loops for Direct and Indirect evidence
# But assumes that they have the same random effect mean
# So it lumps only on the relative effect and allows all other parameters to be source-specific 


Lumping_T2_model <- function(){
  
  for(i in 1:ns.dir){                       # LOOP THROUGH DIRECT STUDIES
    
    delta[i,1] <- 0                         # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                  # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                    # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])          # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]   # model for linear predictor
      
    }                                       # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2],prec.dir)
    
  }                                         # DIRECT studies Loop closes   
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){       # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                         # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                  # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                    # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])          # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]   # model for linear predictor
      
    }                                       # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.dir[2],prec.indir)
    
  }                                         # Indirect studies Loop closes  
  
  d.dir[1]<-0                               # treatment effect is zero for reference treatment in Direct
  
  for (k in 2:nt){  d.dir[k] ~ dnorm(0,.0001) }   # vague priors for treatment effects
  
  tau.dir ~ dunif(0.000001, 5)              # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)               # between-trial precision = (1/between-trial variance)
  tau.indir ~ dunif(0.000001, 5)            # vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)           # between-trial precision = (1/between-trial variance)
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, tau.dir, tau.ind, mu

#### Lumping under a common Random-effect  ####

# This model is just a simple Random-effects model which is used to Lump all sources of evidence on all parameters
# So this model does not lump only on relative effectiveness but on heterogeneity as well. 
# It also assumes exchangeability across all studies, regardless of the source they pertain to


RE_T2_model <- function(){
  
  for(i in 1:ns){                      # LOOP THROUGH STUDIES
    w[i,1] <- 0                        # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0                    # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)             # vague priors for all trial baselines
    
    for (k in 1:na[i]) {               # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])     # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]  # model for linear predictor
    }
    
    delta[i,2] ~ dnorm(d[2], prec)
    
  }   
  d[1]<-0                              # treatment effect is zero for reference treatment
  
  
  for (k in 2:nt){  d[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  tau ~ dunif(0.000001, 5)                  # vague prior for between-trial SD
  prec <- pow(tau,-2)                       # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns = ns.dir+ns.indir, nt, na, r, n

## Starting values ##
# Need to include d, tau, mu


#### Lumping on the between-studies heterogeneity ONLY  ####

# This model lumps the two sources of evidence on the between-sources heterogeneity ONLY.
# The heterogeneity is simultaneously estimates in both sources


LUMP_HET_T2_model <- function(){
  
  for(i in 1:ns.dir){                       # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                         # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                  # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                    # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])          # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]   # model for linear predictor
      
    }                                       # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2],prec.dir)
    
  }                                         # Direct studies Loop closes   
  
  d.dir[1]<-0                               # treatment effect is zero for reference treatment in Direct
  
  
  for (k in 2:nt){  d.dir[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  tau.dir ~ dunif(0.000001, 5)                  # vague prior for between-trial SD
  prec.dir <- pow(tau,-2)                   # between-trial precision = (1/between-trial variance)
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){           # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
      
    }                                           # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2],prec.dir)
    
  }                                             # Indirect studies Loop closes  
  
  d.indir[1]<-0                                 # treatment effect is zero for reference treatment
  
 
  for (k in 2:nt){  d.indir[k] ~ dnorm(0,.0001) }  # vague priors for treatment effects
  prec.indir <- pow(tau,-2)                  # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, d.ind, tau, mu


#### Random-walk ####

# This model assumes that that the relative effect of the direct evidence (i.e. d.dir) 
# is centered around that of the indirect evidence (i.e. d.indir)
# The variance is estimates within the model and controls the extent to which the two parameters are forced to be the same


RW_T2_model <- function(){
  
  for(i in 1:ns.dir){                           # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
      
    }                                           # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2],prec.dir)
    
  }                                             # Direct studies Loop closes   
  
  d.dir[1]<-0                                   # treatment effect is zero for reference treatment in Direct
  
  

  for(k in 2:nt){  d.dir[k]~dnorm(d.indir[k], prec.pop) }  # d.dir is centered around d.indir. This substitutes the prior of d.dir 
  prec.pop<-pow(tau.pop, -2)
  tau.pop~dunif(0.000001, 5)
  
  
  tau.dir ~ dunif(0.000001, 5)                  # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)                   # between-trial precision = (1/between-trial variance)
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){           # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
    }                                           # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2],prec.dir)
    
  }                                             # Indirect studies Loop closes  
  
  
  d.indir[1]<-0                                 # treatment effect is zero for reference treatment
   
  
  for (k in 2:nt){  d.indir[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  tau.indir ~ dunif(0.000001, 5)                #  vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)               # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include n.dir, n.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, d.ind, tau.dir, tau.indir, tau.pop, mu


#### Multi-level model ####

# This model assumes that the relative effects of the direct and the indirect evidence (i.e. d.dir & d.indir)
# are drawn from a random effect. I.e. the 3rd level.
# d.dir and d.indir are hence shrunk towards the 3rd level randome effect mean
# The variance of the 3rd level random effect is indicative of the across sources heterogeneity on the treatment effects

MLM_T2_model <- function(){
  
  
  for(i in 1:ns.dir){                           # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
      
    }                                           # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                             # Direct studies Loop closes   
  
  d.dir[1]<-0                                   # treatment effect is zero for reference treatment in Direct
  
  tau.dir ~ dunif(0.000001, 5)                  # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)                   # between-trial precision = (1/between-trial variance)
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){           # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
    }                                           # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2], prec.indir)
    
  }                                             # Indirect studies Loop closes  
  
  d.indir[1]<-0                                 # treatment effect is zero for reference treatment
  
  tau.indir ~ dunif(0.000001, 5)                # vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)               # between-trial precision = (1/between-trial variance)
  
  
  for(k in 2:nt){  # Imposing the additional level. Basic-parameter-specific multi-level model
    d.dir[k]~dnorm(d[k], prec.pop)              # d[k] is hyper-param. I.e. Across sources relative effect 
    d.indir[k]~dnorm(d[k], prec.pop)
  }
  
  for(k in 2:nt){ d[k]~dnorm(0, .001)}          # prior for hyper-param (i.e. hyper-prior)
  prec.pop<-pow(tau.pop, -2)
  tau.pop~dunif(0.000001, 5)                    # heterogeneity across sources' relative treatment effects 
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, d.ind, d, tau.dir, tau.indir, tau.pop, mu

#### Commensurate prior ####

# This model is similar to the random walk in the sense that d.dir is centered around d.indir
# The difference is that the variance of this, effectively prior, is also assumed to follow a spike-and-slab hyperprior
# The spike and slab hyperprior is a mixture of high variance values, which effectivelly do not force commensurability, and low variance values that force commensurability
# The mixture weights for the spike and the slab can either be assumed fixed or estimated within the model (here 50/50 weights are assumed)

PR_COMMENSURATE_T2_model <- function(){
  
  
  for(i in 1:ns.dir){                           # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
      
    }                                           # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                             # Direct studies Loop closes   
  
  d.dir[1]<-0       # treatment effect is zero for reference treatment in Direct
  
  
  # Commensurate prior for d.dir
  
  for(k in 2:nt){  d.dir[k]~dnorm(d.indir[k], ssprec[k]) }
  
  ## ssprec = 1000 forces commensurability
  ## ssprec = 0.001 disconnects Direct and Indirect
  ##  ssprec ~ dgamma(.01, .01) # Standard WinBUGS vague hyperprior
  
  for (k in 2:nt) {
    tee[k,1] ~ dnorm(20,1)                      # R_tau is (essentially) 20
    tee[k,2] ~ dgamma(0.1, 0.1) %_%I(0.1, 5)    # replacement for the true slab
    flip[k] ~ dbern(0.5)                        # p_tau is 0.5. FIXED HERE but can be uncertain (i.e. estimated in the model)
    pick[k] <- flip[k] + 1
    ssprec[k] <- tee[k,pick[k]]
    sstau[k] <- sqrt(1/ssprec[k])
  }
  # for (k in 2:nt){sstau[k]<-sstau.com} # may be used in NMAs 
  
  tau.dir ~ dunif(0.000001, 5)                      # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)                       # between-trial precision = (1/between-trial variance)
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){               # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                                 # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                          # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                            # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])                  # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]           # model for linear predictor
    }                                               # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2], prec.indir)
    
  }                                                 # Indirect studies Loop closes  
  
  
  d.indir[1]<-0       # treatment effect is zero for reference treatment
  
  # vague priors for treatment effects
  for (k in 2:nt){  d.indir[k] ~ dnorm(0,.0001) }
  tau.indir ~ dunif(0.000001, 5)     # vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)   # between-trial precision = (1/between-trial variance)
  
  
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, d.ind, tau.dir, tau.indir, mu
# If the mixture weights of the spike-and-slab components is not fixed, we ll need a prior and a starting value for the probability that goes into the dbern() too



#### Constraint (Direct more negative) #####

# This model assumes that the mean d.dir is smaller than the mean d.indir
# This does imply that the two estimates will not overlap. It only implies that the desired direction is preserved within each MCMC simulation

CONSTRAINT_T2_model <- function(){
  
  for(i in 1:ns.dir){                             # LOOP THROUGH Direct STUDIES
    
    
    delta[i,1] <- 0                              # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                       # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                         # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])               # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]        # model for linear predictor
      
    }                                            # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2],prec.dir)
    
  }                                             # Direct studies Loop closes   
  
  d.dir[1]<-0       # treatment effect is zero for reference treatment in Direct
  
  
  # vague priors for treatment effects
  for (k in 2:nt){  d.dir[k] ~ dnorm(0,.0001) }
  tau.dir ~ dunif(0.000001, 5)     # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)   # between-trial precision = (1/between-trial variance)
  
  
  
  for(i in ns.dir+1:ns.dir+ns.indir){           # LOOP THROUGH Indirect STUDIES
    
    delta[i,1] <- 0                             # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                      # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                        # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])              # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]       # model for linear predictor
      
    }                                           # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2],prec.indir)
    
  }                                             # Indirect studies Loop closes  
  
  
  d.indir[1]<-0       # treatment effect is zero for reference treatment
  
  # vague priors for treatment effects
  for (k in 2:nt){  d.indir[k] ~ dnorm(0,.0001) }
  tau.indir ~ dunif(0.000001, 5)     # vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)   # between-trial precision = (1/between-trial variance)
  
  
  # We want to impose d.dir to be more negative than d.indir 
  # i.e. d.dir < d.indir , d.dir - d.indir < 0,  d.indir - d.dir >= 0
  # i.e. step(d.indir - d.dir) = 1
  
  constraint <- step(d.indir[2] - d.dir[2]) # would need to do it for each [k] in an NMA as long as we think the assumpion holds for each
  b ~ dbern(constraint)
  b <- 1
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n

## Starting values ##
# Need to include d.dir, d.ind, tau.dir, tau.indir, mu


#### Power prior with a= 0.5 #####

# This model discounts the likelihood of the indirect evidence
# To achieve that we have to use the zeros trick in order to program the Custom likelihood function

PR_POWER_T2_model <- function(){
  
  for(i in 1:ns.dir){                             # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                               # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                        # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                          # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])                # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]         # model for linear predictor
      
    }                                             # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
    
  }                                               # Direct studies Loop closes   
  
  d.dir[1]<-0       # treatment effect is zero for reference treatment in Direct
  
  # vague priors for treatment effects
  for (k in 2:nt){  d.dir[k] ~ dnorm(0,.0001) }
  tau.dir ~ dunif(0.000001, 5)                    # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)                     # between-trial precision = (1/between-trial variance)
  
  
  
  const<-10000 # define an arbitrary large constrant
  for(i in ns.dir+1:ns.dir+ns.indir){             # LOOP THROUGH Indirect STUDIES
    
    
    delta[i,1] <- 0                               # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                        # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                          # LOOP THROUGH ARMS
      
      # Implement zeros trick to define custom likelihood
      zeros[i,k]<-0
      zeros[i,k] ~ dpois(phi[i,k])
      phi[i,k] <- alpha * neg.LL[i,k] + const  # alpha comes in as data
      neg.LL[i,k] <- - logfact(n[i,k]) + logfact(r[i,k]) + logfact(n[i,k] - r[i,k]) - r[i,k]*log(p[i,k]) - (n[i,k] - r[i,k])*log(1-p[i,k]) 
      
      logit(p[i,k]) <- mu[i] + delta[i,k]  # model for linear predictor
    }                                      # Arms loop closes  
    
    delta[i,2] ~ dnorm(d.indir[2], prec.indir)
    
  }                                        # Indirect studies Loop closes  
  
  
  # vague priors for treatment effects
  for (k in 1:nt){  d.indir[k] <- d.dir[k] } # imposing common rte
  tau.indir ~ dunif(0.000001, 5)     # vague prior for between-trial SD
  prec.indir <- pow(tau.indir,-2)   # between-trial precision = (1/between-trial variance)
  
  
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n
# AND ADDITIONALLY the desired alpha

## Starting values ##
# Need to include d.dir, d.ind, tau.dir, tau.indir, mu


#### Analyse only indirect data : Intermediate Necessary step ####

# This step is only necessary for two-step models that analyse initially the indirect evidence
# and subsequently use them as data in the analysis of the direct evidence

RE_T2_model <- function(){
  
  for(i in 1:ns){                              # LOOP THROUGH STUDIES
    w[i,1] <- 0                                # adjustment for multi-arm trials is zero for control arm
    delta[i,1] <- 0                            # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                     # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                       # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])             # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]      # model for linear predictor
    }
    
    delta[i,2] ~ dnorm(d[2], prec)
    
  }   
  d[1]<-0                                      # treatment effect is zero for reference treatment
  

  for (k in 2:nt){  d[k] ~ dnorm(0,.0001) }    # vague priors for treatment effects
  tau ~ dunif(0.000001, 5)                     # vague prior for between-trial SD
  prec <- pow(tau,-2)                          # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns=ns.ind, nt, na, r, n


## Starting values ##
# Need to include d (i.e d.ind), tau (i.e. tau.indir), mu


#### Informative-prior using the postrerior mean indirect relative effect ####

# This model is the second step, and requires a previous step where the indirect evidence have been analyses alone  
# It uses an informative priors on d.dir that is based on d.indir and sd.indir
# I.e. the relative effect mean of the indirect evidence and the uncertainty of the MEAN relative effect d.ind 
# It is therefore not appropriate for a future rol-out of a new trial, but only for the MEAN relative effect

PR_INFORMATIVE_RTE_T2_model <- function(){
  
  for(i in 1:ns.dir){                          # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                            # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                     # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                       # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])             # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]      # model for linear predictor
    }                                          # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                            # Direct studies Loop closes   
  
  d.dir[1]<-0                                  # treatment effect is zero for reference treatment in Direct
  
  
  
  # Informative priors for adult treatment effects from the RTEs of Indirect
  # d.indir[k], sd.d.indir[k] come in as data from the previous 1st step where indirect evidence have been analysed alone
  for (k in 2:nt){  d.dir[k] ~ dnorm(d.indir[k], prec.d.indir[k]) } 
  
  for(k in 2:nt){prec.d.indir[k]<-pow(sd.pred.indir[k],-2)}
  sd.pred.indir[2]<- sd.d.indir # comes in as data. is the sd of the mean relative tr. effect
  d.indir[2]<- d.indir          # comes in as data 
  
  # vague prior for between-studies heterogeneity
  tau.dir ~ dunif(0.000001, 5)               # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)                # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns.dir, nt, na, r, n
# AND ADDITIONALLY d.indir, sd.d.indir (from the 1st step that analyses only the indirect data)

## Starting values ##
# Need to include d.dir, tau.dir, mu




#### Informative-prior using the predictive distribution of the indirect relative effect ####

# This model is the second step, and requires a previous step where the indirect evidence have been analyses alone  
# It uses an informative priors on d.dir that is based on d.indir and sd.indir and tau.indir
# I.e. the predictive distribution 
# It may be excessively vague if applied on the mean relative effect of the direct evidence

PR_INF_PREDICTIVE_T2_model <- function(){
  
  for(i in 1:ns.dir){                       # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                         # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                  # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                    # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])          # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]   # model for linear predictor
    }                                       # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                         # Direct studies Loop closes   
  
  d.dir[1]<-0                               # treatment effect is zero for reference treatment in Direct
  
  # Informative priors for adult treatment effects from the RTEs of Indirect
  
  for (k in 2:nt){  d.dir[k] ~ dnorm(d.indir[k], prec.d.indir[k]) } 
  
  for(k in 2:nt){prec.d.indir[k]<-pow(sd.pred.indir[k],-2)}
  sd.pred.indir[2]<- sd.d.indir # comes in as data and is the sd of the predictive distribution 
  d.indir[2]<- d.indir          # comes in as data
  
  
  tau.dir ~ dunif(0.000001, 5)            # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)             # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns.dir, nt, na, r, n
# AND ADDITIONALLY d.indir, sd.d.indir (from the 1st step that analyses only the indirect data). (Note: The sd.d.indir is the sd of the predictive distr)

## Starting values ##
# Need to include d.dir, tau.dir, mu


#### Mixture-prior using the postrerior mean indirect relative effect ####

# This model is the second step, and requires a previous step where the indirect evidence have been analyses alone  
# It uses a mixture prior on d.dir that is based on d.indir and sd.indir AND a vague component
# I.e. the relative effect mean of the indirect evidence and the uncertainty of the MEAN relative effect d.ind 
# It is therefore not appropriate for a future rol-out of a new trial, but only for the MEAN relative effect
# The weight of each component is determined within the model and hence it is an `ADAPTIVE-SHARING` model

PR_MIXTURE_RTE_T2_model <- function(){
  
  for(i in 1:ns.dir){                     # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                       # treatment effect is zero for control arm
    mu[i] ~ dnorm(-0.5,.01)               # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                  # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])        # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      
    }                                     # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                       # Direct studies Loop closes   
  
  d.dir[1]<-0                             # treatment effect is zero for reference treatment in Direct
  
  
  # Informative priors for Direct treatment effects
  # Two components : 1. the RTE of Indirect 2. a vague
  # Mix the two components according to probability vector P
  # essentially p, (1-p)
  # The probability vector is estimated in the model
  # The component distributions are the lambdas
  
  # Needs to be specified for every k in an NMA model
  d.dir[2]<-lambda[T]            
  T~dcat(P[])                 
  P[1:2]~ddirch(alpha[])       
  lambda[1]~dnorm(d.indir[2], prec.d.indir[2]) 
  lambda[2]~dnorm(0, .001)
  
  for(k in 2:nt){prec.d.indir[k]<-pow(sd.rte.indir[k],-2)}
  sd.rte.indir[2]<- sd.d.indir # comes in as data. is the sd of the mean relative tr. effect
  d.indir[2]<- d.indir         # comes in as data
  
  tau.dir ~ dunif(0.000001, 5)     # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)      # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns.dir, nt, na, r, n
# AND ADDITIONALLY d.indir, sd.d.indir (from the 1st step that analyses only the indirect data)
# AND alpha = c(1,1) 

## Starting values ##
# Need to include d.dir, tau.dir, mu
# AND ADDITIONALLY lambda = c(1, 1) and P=c(0.5, 0.5). These are just examples of starting values that can be changed



#### Mixture-prior using the predictive distrubiton of the indirect relative effect #####

# This model is the second step, and requires a previous step where the indirect evidence have been analyses alone  
# It uses a mixture prior on d.dir that is based on d.indir and sd.indir AND a vague component
# I.e. the predictive distribution
# It may be excessively vague if applied on the mean relative effect of the direct evidence
# The weight of each component is determined within the model and hence it is an `ADAPTIVE-SHARING` model


PR_MIXTURE_PREDICTIVE_T2_model <- function(){
  
  for(i in 1:ns.dir){                     # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                       # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                  # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])        # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k] # model for linear predictor
      
    }                                     # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                       # Direct studies Loop closes   
  
  d.dir[1]<-0                             # treatment effect is zero for reference treatment in Direct
  
  # Informative priors for Direct treatment effects
  # Two components : 1. the RTE of Indirect 2. a vague
  # Mix the two components according to probability vector P
  # essentially p, (1-p)
  # The probability vector is estimated in the model
  # The component distributions are the lambdas
  
  # For comparison ALB.vs.PLA
  d.dir[2]<-lambda[T]            
  T~dcat(P[])                 
  P[1:2]~ddirch(alpha[])       
  lambda[1]~dnorm(d.indir[2], prec.d.indir[2]) 
  lambda[2]~dnorm(0, .001)
  
  for(k in 2:nt){prec.d.indir[k]<-pow(sd.rte.indir[k],-2)}
  sd.rte.indir[2]<- sd.d.indir # comes in as data and is the sd of the predictive distribution 
  d.indir[2]<- d.indir         # comes in as data
  
  tau.dir ~ dunif(0.000001, 5)  # vague prior for between-trial SD
  prec.dir <- pow(tau.dir,-2)   # between-trial precision = (1/between-trial variance)
  
  
} 

## Data ##
# Need to include ns.dir, nt, na, r, n
# AND ADDITIONALLY d.indir, sd.d.indir (from the 1st step that analyses only the indirect data) (Note: The sd.d.indir is the sd of the predictive distr)
# AND alpha = c(1,1) 

## Starting values ##
# Need to include d.dir, tau.dir, mu
# AND ADDITIONALLY lambda = c(1, 1) and P=c(0.5, 0.5). These are just examples of starting values that can be changed

#### Log-normal prior on the between-studies heterogeneity ####

# This model requires a previous step where the indirect evidence have been analysed alone
# Importantly, the coda of the tau parameter needs to be saved (e.g. dat <- out.coda[[1]][,"tau"]) when the indirect evidence are initially analysed
# Then, using library(MASS), we can fit a log-normal distribution to the tau coda values and save its estimated parameters. For example,
#           library(MASS)
#           fit_params <- fitdistr(dat,"lognormal")
#           tau.fit.meanlog <- fit_params$estimate[1]
#           tau.fit.sdlog <- fit_params$estimate[2]
# The two parameters tau.fit.meanlog and tau.fit.sdlog will then need to come in to the analysis of the direct evidence as data

LOGNORMAL_PR_HET_T2_model <- function(){
  
  for(i in 1:ns.dir){                             # LOOP THROUGH Direct STUDIES
    
    delta[i,1] <- 0                               # treatment effect is zero for control arm
    mu[i] ~ dnorm(0,.0001)                        # vague priors for all trial baselines
    
    for (k in 1:na[i]) {                          # LOOP THROUGH ARMS
      r[i,k] ~ dbin(p[i,k],n[i,k])                # binomial likelihood
      logit(p[i,k]) <- mu[i] + delta[i,k]         # model for linear predictor
      
    }                                             # Arms loop closes
    
    delta[i,2] ~ dnorm(d.dir[2], prec.dir)
    
  }                                               # Direct studies Loop closes   
  
  d.dir[1]<-0        # treatment effect is zero for reference treatment in Direct
  
  
  for (k in 2:nt){  d.dir[k] ~ dnorm(0,.0001) } # vague priors for treatment effects
  
  # Informative prior on the heterogeneity from the neonate heterogeneity
  
  tau ~ dlnorm(pr.tau.mean, pr.tau.prec) # log-normal prior on tau
  pr.tau.mean <- tau.fit.meanlog         # 1st parameter of log-normal prior (i.e. mean-log). Comes in as data 
  pr.tau.prec <- pow(tau.fit.sdlog, -2)  # 2nd parameter of log-normal prior (i.e. sd-log). Comes in as data
  
  prec.dir <- pow(tau,-2)   # between-trial precision = (1/between-trial variance)
  
} 

## Data ##
# Need to include ns.dir, ns.ind, nt, na, r, n
# AND ADDITIONALLY tau.fit.meanlog, tau.fit.sdlog

## Starting values ##
# Need to include d.dir, tau.dir, mu

