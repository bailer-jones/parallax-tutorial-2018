rootDir <- "/Users/coryn/local/astrostats_2017_ESAC/" # directory containing notebook, Rcode, and data directories

setwd(paste(rootDir, "Rcode", sep=""))
conv <- pi/180 
source("general_functions.R")
source("distance_functions.R")
setwd("../")
library(MASS) # for truehist
library(fields) # for image.plot
library(RColorBrewer) # for colorRampPalette
mypalette <- colorRampPalette(brewer.pal(9, "Greys"), space="rgb", interpolate="linear", bias=2.5)
mycols <- mypalette(64)

# Plot the posteriors to get feel for what they look like
par(mfrow=c(2,2), mar=c(5,5,0.5,1), oma=c(0.1,0.1,0.5,0.1), mgp=c(2.2,0.8,0), cex=1.0) 
# Draw from distance priors (as a check of the code) and plot 
# 1,2,3 are uniform in distance, uniform in space density, and exponentially decreasing space density, respectively.
truehist(r.distprior1(Nsamp=1e5, rmax=1e3), h=20)
truehist(r.distprior2(Nsamp=1e5, rmax=1e3), h=20)
truehist(r.distprior3(Nsamp=1e5, rlen=1e3), h=20)

# Plot unnormalized posteriors over distance for some parallax measurement
# and overplot mode (red), quantiles (blue), 
# and naive parallax inversion and first order uncertainty propagation (green)
w    <- 1e-4
wsd  <- 1e-5
rlen <- 1e3
r <- seq(from=0, to=2e4, length.out=1e3)
plot(r, ud.distpost3(r, w, wsd, rlen), type="l", xlab="distance r [pc]", ylab="P*(r|w,wsd,rlen)")
mode.distpost3(w, wsd, rlen, retall=TRUE)
rMode  <- mode.distpost3(w, wsd, rlen)
rFWHM  <- fwhm.distpost3(w, wsd, rMode=rMode, rlen, rmax=1e6) 
rQuant <- quantiles.distpost3(w, wsd, rlen, probs=c(0.05, 0.5, 0.95), Nsamp=1e5)
abline(v=c(rMode, rQuant), col=c("red", rep.int("blue", length(rQuant))))
M <- ud.distpost3(r=rMode, w, wsd, rlen)
lines(x=c(rFWHM), y=c(M,M)/2, col="red")
abline(v=c(1/w, (1/w)+c(+1,-1)*wsd/w^2), col="green")
c(mode=rMode, rFWHM, rQuant)

# parallaxes are in mas, angles in degrees
dat <- read.csv(paste(rootDir, "data/gdr1set01.csv", sep=""), sep=",")
# convert parallaxes to arcsec so that distances come out as pc
dat$parallax <- 1e-3*dat$parallax
dat$parallax_error <- 1e-3*dat$parallax_error

dim(dat)
# To reduce data size for testing purposes
dat <- dat[1:10,]
dat

rMode  <- Vectorize(mode.distpost3, c("w", "wsd"))(w=dat$parallax, wsd=dat$parallax_error, rlen)
rMode

# The root-finding algorithm in fwhm.distpost3 may give warnings
rFWHM <- Vectorize(fwhm.distpost3, c("w", "wsd", "rMode"))(w=dat$parallax, wsd=dat$parallax_error, rMode=rMode, rlen=rlen, rmax=1e6)
rFWHM <- t(rFWHM)
cbind(rMode, rFWHM)

# My code for computing the quantiles computes the normalization using MCMC, which is very slow.
# A grid would be faster, but its size/spacing would generally need to be optimized for every source separately.
rQuant <- Vectorize(quantiles.distpost3, c("w", "wsd"))(w=dat$parallax, wsd=dat$parallax_error, rlen, probs=c(0.05, 0.5, 0.95), Nsamp=1e5)
rQuant <- t(rQuant)

cbind(w=dat$parallax, fobs=dat$parallax_error/dat$parallax, rMode, rFWHM, rQuant)

par(mfrow=c(2,1), mar=c(5,5,0.5,1), oma=c(0.1,0.1,0.5,0.1), mgp=c(2.2,0.8,0), cex=1.0) 
plot(rMode, rQuant[,2], xlab="mode distance [pc]", ylab="median distance [pc]")
abline(a=0, b=1)
plot(1/dat$parallax, rMode, xlab="inverse parallax [pc]", ylab="mode distance [pc]")
abline(a=0, b=1)

# Use exponentially decreasing space density prior, but marginalize over length scale parameter, rlen. 
# To do this adopt a hyperprior P(rlen), with non-explicit hyperparameters. The new posterior is
# P(r|w) = int P(r,rlen|w) d(rlen)
# P*(r|w) = int P(w|r)P(r|rlen)P(rlen)
# Define the above integrand (to within normalization factor 1/P(w)) as a function, whereby first argument is rlen:
# - hyperprior is a scaled/shifted beta function (defined over rlenRange) with hard-coded hyperparameters.
# - lower limit (min(rlenRange)) should be >0, otherwise prior P(r|rlen) evaluates to NaN (prior undefined)
# This function is P*(r|w), and "hier" in name indicates it is from a "hierarchical" model
ud.distpost3hier <- function(rlen, r, w, wsd, rlenRange) {
  d.like(w=w, r=r, wsd=wsd)*d.distprior3(r=r, rlen=rlen)*d.betagen(x=rlen, xrange=rlenRange, shape1=1, shape2=1) 
}
# Test integration using this integrand
w <- 5e-3
wsd <- 1e-3
r <- 4e2
rlenRange <- c(1e1, 1e3)
integrate.func(f=ud.distpost3hier, lower=rlenRange[1], upper=rlenRange[2], r=r, w=w, wsd=wsd, 
               rlenRange=rlenRange, subdivisions=1e3, rel.tol=1e-6, abs.tol=1e-6) 

# Compute P*(r|w) on a grid
r <- seq(from=0, to=1e3, length.out=1e3)
udDistpost3hier <- double(length=length(r))
for(i in 1:length(r)) {
  udDistpost3hier[i] <- integrate.func(f=ud.distpost3hier, lower=rlenRange[1], upper=rlenRange[2], r=r[i], 
                                       w=w, wsd=wsd, rlenRange=rlenRange, 
                                       subdivisions=1e3, rel.tol=1e-6, abs.tol=1e-6) 
}
# TODO: Some integrals don't converge so return NA. Just set to zero for now
udDistpost3hier[is.na(udDistpost3hier)] <- 0
# Plot and overplot what we get with two different fixed values of rlen
# Scale plots to have maximum of 1 (would be better to normalize them)
mom <- pdfmom(udDistpost3hier, r)
cat("Posterior mean, SD for r =", mom$mean, mom$sd, "kpc\n")
udDistpost3hier <- udDistpost3hier/mom$Z # is now normalized
plot(r, udDistpost3hier, type="l", xlab="distance r [pc]", ylab="P(r|w,wsd)")
udDistpost3 <- ud.distpost3(r, w, wsd, rlen=min(rlenRange))
lines(r, udDistpost3/pdfmom(udDistpost3, r)$Z, lty=2, col="red") # posterior with rlen fixed to minimum of hyperprior range
udDistpost3 <- ud.distpost3(r, w, wsd, rlen=max(rlenRange))
lines(r, udDistpost3/pdfmom(udDistpost3, r)$Z, lty=2, col="blue") # posterior with rlen fixed to maximum of hyperprior range

