# Kullback-Leibler Divergence Calculation

# The following R function calculates the KL divergence of an estimated normal distribution, given by Norm(mean, sd)
# from a standard normal distribution given by an MCMC sample y

# It uses a kernel density estimate of the true distribution and numerical integration

## INTEGRATE ACROSS THE WHOLE RANGE OF THE MASS OF THE PARAMETER
# i.e. if it is a probability from 0 to 1
# Here we consider relative effects which are defined from -infty to + infty

# I INTEGRATE IN THE FOLLOWING RANGE 
# The lowest value is the lowest p=0.0001 percentile 
# The highest value is the highest p=0.9999 percentile 


mean.true <- 0
sd.true <- 1

y <- rnorm(1e5, mean.true, sd.true) # In an application the CODA is used here

min <- qnorm(p=0.001, mean=mean.true, sd=sd.true)
max <- qnorm(p=0.999, mean=mean.true, sd=sd.true)

#Define KL function

kl <- function(y, mean, sd){
  d <- density(y, from=min, to=max)
  dfn <- with(d, approxfun(x=x, y=y))
  fn <- function(x){
    p <- dfn(x)
    q <- dnorm(x, mean, sd)
    p*(log(p) - log(q))
  }
  
  integrate(fn,min,max)
  
}


# Imagine a plethora of estimated distributions with different (biased) means between -1 and 1

mean.est <- seq(-1, 1, by=0.05)
x <- seq(-4, 4, by=0.005)

plot(x, dnorm(x, mean.true, sd.true), type="l", lwd=4, ylim=c(0,0.41), xlab="Mean", ylab="pdf")

for(i in seq(along=mean.est)){
  lines(x, dnorm(x, mean.est[i], sd.true), col="gray")
}

legend("topleft", col=c("black", "gray"), lwd=c(3,1), c("True", "Simulated"))



# Imagine a plethora of estimated distributions with different (over-precise/under-precise) sd

sd.est <- seq(0.4, 1.6, by=0.03)

plot(x, dnorm(x, mean.true, sd.true), type="l", lwd=3, ylim=c(0,1), xlab="Mean", ylab="pdf")

for(i in seq(along=mean.est)){
  lines(x, dnorm(x, mean.true, sd.est[i]), col="gray")
}

legend("topleft", col=c("black", "gray"), lwd=c(3,1), c("True", "Simulated"))


# Calculating KL

kl.est <- matrix(nrow=length(mean.est), ncol=length(sd.est))

for(i in seq(along=mean.est)){
  for(j in seq(along=sd.est)){
    kl.est[i,j] <- kl(y, mean.est[i], sd.est[j])$value
  }
}

# Relation of Kullback-Leibler to sd, given correct mean

itrue <- which(mean.est==mean.true)

plot(sd.est, kl.est[itrue, ], type="l", xlim=c(0.4,1.8), xlab="Varying sd", ylab="KL div")
abline(v=sd.true, col="gray"); text(x=sd.true, y=0.06, pos=4, "True sd", col="gray")

# Relation of Kullback-Leibler to mean, given correct sd

jtrue <- which(sd.est==sd.true)

plot(mean.est, kl.est[,jtrue], type="l", xlim=c(-1,1), xlab="Varying mean", ylab="KL div")
abline(v=mean.true, col="gray"); text(x=mean.true, y=0.1, pos=4, "True mean", col="gray")


# Contour (Heat map)

library(RColorBrewer)
plot(mean.est, sd.est, type="n", xlab="Simulated mean", ylab="Simulated sd")
levs <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1, 1.5, 2, 2.5, 3, 4, 5, 6)
getPalette = colorRampPalette(brewer.pal(9, "Greens"))
.filled.contour(x=mean.est, y=sd.est, z=kl.est, levels=levs, col=getPalette(length(levs)-1))
contour(x=mean.est, y=sd.est, z=kl.est, levels=levs, add=TRUE)
abline(h=sd.true, col="gray", lty="dashed")
abline(v=mean.true, col="gray", lty="dashed")
