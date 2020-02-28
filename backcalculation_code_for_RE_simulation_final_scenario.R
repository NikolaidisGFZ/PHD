# Back-calculation of small studies probability of an event in treatment arm 
# under either FE or RE
# Based in inverse-variance weights 
# George Nikolaidis


rm(list=ls()) # clear workspace

d_all <- -0.51     # The overall log-odds ratio of the indirect studies
TAU_ALL <- 0.56    # OVerall heterogeneity among all studies

d_big <- qnorm(0.95, mean=d_all, sd=TAU_ALL) 

p_ctl <- 0.366      # Assumes probability of event in control arm

f <- function(x , d_big, p_ctl) ( exp(d_big) - ( ( x*(1-p_ctl) ) 
                                                 / ( p_ctl*(1-x) )  ) )

# find the probability of an event in the treatment arm of the big study
p_big <- uniroot(f, d_big=d_big, p_ctl=p_ctl, lower=0.00000001, upper=1)$root 

N_studies_all <- 10   # total number of indirect studies
N_studies_big <- 1    # number of big studies (probably only 1)
N_studies_small <- N_studies_all - N_studies_big  #number of small studies

SS_ALL_IND <- 2300*8    # total number of indirect patients
pat_share <- 0.85       # proportion of indirect patients included in the biggest study

# number of patients in big study (in both arms)
ss_big <- round(SS_ALL_IND * pat_share, digits=0)  

# number of patients in small studieS (in both arms)
ss_sm <- round((SS_ALL_IND - ss_big) / N_studies_small, digits=0) 

# derivation of within BIG-STUDY LOR variance 1/a + 1/b + 1/c + 1/d
sigma_big_sq <- (1 / (p_ctl * (ss_big/2) ))  +    
  (1 / ((1 - p_ctl) * (ss_big/2) ) )  + 
  (1 / (p_big * (ss_big/2) ))  + 
  (1 / ((1 - p_big) * (ss_big/2) ) )

w_big <- 1 / sigma_big_sq                     # Fixed-effect weight of big study

w_big_star <- 1 / (sigma_big_sq + TAU_ALL^2)  # Random-effect weight of big study


# This bit calculates the probability for a RANDOM-EFFECTS model and converts it to d
f_RE <- function(p_sm ,d_all, d_big, p_ctl, TAU_ALL, ss_sm, w_big_star) (  
  d_all - (   (  (9 * 
                    ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) / 
                                 ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )  
                            + TAU_ALL^2 )  ) * 
                    (log(  ( p_sm*(1-p_ctl) )  /  ( p_ctl*(1-p_sm) )  ) ) ) + 
                   (w_big_star*d_big) )  
              /   ( (9 * ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) / 
                                      ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )  
                                 + TAU_ALL^2 )  )) + w_big_star )  )  )


# probability of an event in the treatment arm of a small study
p_sm_RE <- uniroot(f_RE, d_all=d_all ,d_big=d_big, p_ctl=p_ctl, 
                   TAU_ALL=TAU_ALL, ss_sm=ss_sm, w_big_star=w_big_star, 
                   lower=0.05, upper=0.95, extendInt = "yes")$root  

d_sm_RE <- log(  ( p_sm_RE*(1-p_ctl) )  /  ( p_ctl*(1-p_sm_RE) )  ) 

print(c(p_sm_RE, p_big, d_sm_RE, d_big))








### Heterogeneity recalculation ###

# This is very important !!!
# the probabilities that were derived above, do not necessarily preserve the TAU_ALL
# This bit re-calculates the TAU_ALL

w_sm_est <- ( ( p_sm_RE*(1-p_sm_RE)*p_ctl*(1-p_ctl)*ss_sm  )  
              /  (  ( 2*p_sm_RE*(1-p_sm_RE) ) + ( 2*p_ctl*(1-p_ctl) )  )  )


tau_all_est_sq <- (( 9*w_sm_est*d_sm_RE^2 + w_big*d_big - 
                       9*w_sm_est*d_all^2 - w_big*d_all^2 - 9)    
                   /  (  9*w_sm_est + w_big  -  (  ( 9*w_sm_est^2 + w_big^2 )  
                                                   /  ( 9*w_sm_est + w_big ) ) ) )

tau_all_est <- sqrt(tau_all_est_sq)

print(tau_all_est)


















###### Calculates the probability for a FIXED-EFFECTS model and converts it to d
# Be careful, the Fixed Effects model `breaks` much easier, 
# i.e there is no value for the small studies 
# preserving the original predictive distribution

f_FE <- function(p_sm ,d_all, d_big, p_ctl, ss_sm, w_big) (  
  d_all - (   (  (9 * ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) 
                                 / ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )   )  ) 
                  * (log(  ( p_sm*(1-p_ctl) )  /  
                             ( p_ctl*(1-p_sm) )  ) ) ) + (w_big*d_big) ) 
              /   ( (9 * ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) 
                                    / ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )   )  )) 
                    + w_big )  )  )

# probability of an event in the treatment arm of a small study
p_sm_FE <- uniroot(f_FE, d_all=d_all ,d_big=d_big, 
                   p_ctl=p_ctl, ss_sm=ss_sm, w_big=w_big, 
                   lower=0.05, upper=0.75, extendInt = "yes")$root  

d_sm_FE <- log(  ( p_sm_FE*(1-p_ctl) )  /  ( p_ctl*(1-p_sm_FE) )  ) 

print(c(p_sm_FE, p_sm_RE, p_big, d_sm_FE, d_sm_RE, d_big))















###### Formulas validation code #### 

p_sm <- 0.22 # Assume a probability for the small studies

# derive the variance of the LOR for a small study given p_sm
sigma_sm_sq <- (1 / (p_ctl * (ss_sm/2) ))  +   
  (1 / ((1 - p_ctl) * (ss_sm/2) ) )  + 
  (1 / (p_sm * (ss_sm/2) ))  + 
  (1 / ((1 - p_sm) * (ss_sm/2) ) )

w_sm <- 1 / sigma_sm_sq # derive in FE weight

w_sm_star <- 1 / (sigma_sm_sq + TAU_ALL^2) # derive its RE weight

# check this is equivalent to w_sm
w_sm_2 <- ( ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm  )  
            /  (  ( 2*p_sm*(1-p_sm) ) + ( 2*p_ctl*(1-p_ctl) )  )  ) 

# check this is equivalent to w_sm_star
w_sm_star_2 <-  ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) 
                           / ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )  + TAU_ALL^2 )  ) 

d_sm <- (log(  ( p_sm*(1-p_ctl) )  /  ( p_ctl*(1-p_sm) )  ) )


## If the following formula gives something very close to zero then the formula is correct
zero <-  d_all - (   (  (9 * ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) 
                                        / ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )  + TAU_ALL^2 )  ) 
                         * (log(  ( p_sm*(1-p_ctl) )  /  ( p_ctl*(1-p_sm) )  ) ) ) + (w_big_star*d_big) )  
                     /   ( (9 * ( 1 / ( (  ( 2*p_sm*(1-p_sm) + 2*p_ctl*(1-p_ctl) ) 
                                           / ( p_sm*(1-p_sm)*p_ctl*(1-p_ctl)*ss_sm )  )  + TAU_ALL^2 )  )) + w_big_star )  )

#### End of formulas validation code
