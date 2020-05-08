# Lab 04 - Bootstrap
# Author: Alberto Nieto Sandino
# Date: 2018-10-09

# Clean environment 
closeAllConnections()
rm(list=ls())
# Clear plots
graphics.off()
# Number of cores to use for bootstrap
library("parallel")
no_cores <- detectCores() - 1 
library("boot")

########################################################################
## Assignment 1. The gamma distribution and the estimator sample mean ##
########################################################################

n <- 100 # number of samples
k <- 2 # shape parameter
gamma <- 2 # scale parameter
dataset <- rgamma (n, shape = k, scale = gamma)

# Histogram of the data
hist(dataset)

# Theoretical mean and variance
mean_theo <- k*gamma
var_theo <- k*gamma^2

# Our parameter of interest theta = mean
# thus theta_hat = sample mean

# Data-generating experiment
N = 100 # repeated 100 times
sample_mean <- numeric(N) # initiate vector for sample mean values
sample_var <- numeric(N) # initiate vector for sample variance values
for (i in 1:N){
  dataset = rgamma (n, shape = k, scale = gamma)
  sample_mean[i] <- mean(dataset)
  sample_var[i] <- var(dataset)
}

# Theoretical mean and variance of the estimator (sample mean)
mean_hat <- mean(sample_mean)
var_hat <- mean(sample_var)

# Observe the distribution for the sample mean
hist(sample_mean)
abline(v=mean_theo,lwd=2, col="red")

############################################
## Assignment 2. Non-parametric bootstrap ##
############################################

## Assignment 2.1.

# 1) Fix the number of bootstrap re-samples N
N <- 1000
# Generate the empirical cdf 
dataset.big <- rgamma (2000, shape = k, scale = gamma)
# 2) Sample a new data set x* set of size n from x with replacement (this is equivalent
# to sampling from the empirical cdf F^)

BootstrapMean <- function (X=dataset.big){
  x.boot <- sample(X, size = 100, replace = TRUE)
  mean(x.boot)
}

# Use replicate to perform N iterations
boot.replicate <- replicate(N, BootstrapMean())

# Plot of the histogram with the true value of the mean
hist(boot.replicate,breaks=10)
abline(v=mean_theo,lwd=2, col="red")

# Calculate estimted bias and variance of estimator
boot.estimate <- mean(boot.replicate)
orig.estimate <- mean(dataset.big)
bias.nonparam <- boot.estimate - orig.estimate
bias.nonparam
var.nonparam <- var(boot.replicate)
var.nonparam

## Assignment 2.2

# Load library
library(boot)

# First, generate a vector of means for bootstrapping
fun <- function(y = dataset.big, id) {mean(y[id])}

# Boot function for N replicates
boot.out <- boot(dataset.big, fun, 1000, parallel="multicore", ncpus = no_cores)

# Plot of the histogram with the true value of the mean
hist(boot.out$t,breaks=13)
abline(v=mean_theo,lwd=2, col="red")

# Calculate estimted bias and variance of estimator
boot.estimate.boot <- mean(boot.out$t)
bias.nonparam.boot <- boot.estimate.boot - orig.estimate
bias.nonparam.boot
var.nonparam.boot <- var(boot.out$t)
var.nonparam.boot

## Assignment 2.3

# Construc CIs for the estimator

# Basic CI
alpha <- 0.05
l.low <- 2*orig.estimate - sort(boot.out$t)[(1-alpha/2)*1000] # lower limit of CI
l.upper <- 2*orig.estimate - sort(boot.out$t)[(alpha/2)*1000] # upper limit of CI
ci.basic <- c(l.low,l.upper)

hist(boot.out$t,20)
lines(rep(ci.basic[1],2),c(0,850),lty=2)
lines(rep(ci.basic[2],2),c(0,850),lty=2)

# Normal CI
se.nonparam <- sd(boot.out$t)
l.low <- orig.estimate + qnorm(alpha/2, mean = 0, sd = 1)*se.nonparam/sqrt(100)
l.upper <- orig.estimate - qnorm(alpha/2, mean = 0, sd = 1)*se.nonparam/sqrt(100)
ci.normal<- c(l.low,l.upper)

hist(boot.out$t,20)
lines(rep(ci.normal[1],2),c(0,350),lty=2)
lines(rep(ci.normal[2],2),c(0,350),lty=2)

# Percentile CI
l.low <- quantile(boot.out$t, probs = 0.025)
l.upper <- quantile(boot.out$t, probs = 0.975)
ci.percentile <- c(l.low,l.upper)

hist(boot.out$t,20)
lines(rep(ci.percentile[1],2),c(0,350),lty=2)
lines(rep(ci.percentile[2],2),c(0,350),lty=2)

# Pre-programmed confidence intervals 
ci.nonparam <- boot.ci(boot.out = boot.out, conf = 0.95, type = c("basic","norm","perc"), parallel="multicore", ncpus = no_cores)

error.ci <- c(ci.basic,ci.normal,ci.percentile) - 
              c(ci.nonparam$normal[2:3],ci.nonparam$basic[4:5],ci.nonparam$percent[4:5])

## Assignment 2.4

N <- 1000 # number of samples to draw
counter <- rep(0,10)
for (i in 1:N){
  for (j in 1:10){
    dataset.small <- rgamma(j*10, shape = k, scale = gamma)
    boot.out <- boot(dataset.small, fun, N)
    ci <- boot.ci(boot.out, conf = 0.95, type = "norm")$norm
    if (ci[2] <= mean_theo && ci[3] >= mean_theo){
      counter[j] <- counter[j] + 1 # true value within limits
    }
  }
}

plot(counter,xaxt = "n",)
axis(1,at=seq(1,10),labels = seq(10,100,by=10))

########################################
## Assignment 3. Parametric bootstrap ##
########################################

## Assigment 3.1

# Create the distribution test for a gamma distribution
x <- rgamma (length(dataset.big), shape = 2, scale = 2)


mlogl <- function(theta, x) {
  if (length(theta) != 2)
    stop("theta must be vector of length 2")
  k.guess <- theta[1]
  gamma.guess <- theta[2]
  if (k.guess <= 0) stop("The guess for the shape parameter, k, must be positive")
  if (gamma.guess <= 0) stop("The guess for the scale parameter, gamma, must be positive")
  return(- sum(dgamma(x, shape = k.guess, scale = gamma.guess, log = TRUE)))
}

k.guess.start <- mean(x)^2 / var(x)
gamma.guess.start <- var(x) / mean(x)
theta.start <- c(k.guess.start, gamma.guess.start)

fit <- nlm (mlogl, theta.start, x = x, hessian = TRUE, fscale = length(x))
print(fit$estimate)

k_MLE <- fit$estimate[1]
gamma_MLE <- fit$estimate[2]
## Assignment 3.2
n <- 100
N <- 2000
# Generation of random samples gamma distributed with parameters from MLE
gamma_MLE <- rgamma(n, shape = k_MLE, scale = gamma_MLE)

# Boot function for N replicates
mean.param <- rep(0,N)
for (i in 1:N){
  resample_MLE <- sample(gamma_MLE, n, replace = TRUE, prob = NULL)
  mean.param[i] <- mean(resample_MLE)
}

# Plot the sample means 
hist(mean.param)

# Calculate the bias and the variance
bias.param <- mean(mean.param) - mean(gamma_MLE)
bias.param
var.param <- var(mean.param)
var.param


## Assignment 3.3

# Description of how random numbers are generated in parametric bootstrap
gengamma <- function(gendata, theta){
  rgamma(100, shape = theta[1], scale = theta[2])
}

# Generate a vector of means for bootstrapping
fun_param <- function(y = dataset, id) {mean(y[id])}

boot.out.param <- boot(dataset, fun_param, 2000, sim = "parametric", ran.gen = gengamma,
                       mle = fit$estimate, parallel="multicore", ncpus = no_cores)
boot.out.param

plot(boot.out.param)
hist(boot.out.param$t)
abline(v=mean_theo,lwd=2, col="red")

# Calculate estimted bias and variance of estimator
boot.estimate.param <- mean(boot.out.param$t)
bias.param <- boot.estimate.param - orig.estimate
bias.param
var.param <- var(boot.out.param$t)
var.param


# Assignment 3.4

N <- 1000 # number of samples to draw
counter <- rep(0,10)
for (i in 1:N){
  for (j in 1:10){
    dataset.small <- rgamma(j*10, shape = k, scale = gamma)
    # re-estimate parameters
    param_MLE <- fit$estimate
    boot.out <- boot(dataset.small, fun, N, sim = "parametric",
                     ran.gen = gengamma, mle = param_MLE)
    ci <- boot.ci(boot.out, conf = 0.95, type = "norm")$norm
    if (ci[2] <= mean_theo && ci[3] >= mean_theo){
      counter[j] <- counter[j] + 1 # true value within limits
    }
  }
}

plot(counter,xaxt = "n",)
axis(1,at=seq(1,10),labels = seq(10,100,by=10))

#################################
## Assignment 4.Studentized CI ##
#################################


## BASIC CONFIDENCE INTERVAL

fun_var <- function(y, id) {var(y[id])}

N <- 1000 # number of samples to draw
counter <- rep(0,10)
for (i in 1:N){
  for (j in 1:10){
    dataset.small <- rgamma(j*10, shape = k, scale = gamma)
    boot.out <- boot(dataset.small, fun_var, N)
    ci <- boot.ci(boot.out, conf = 0.95, type = "basic")$basic
    if (ci[4] <= var_theo && ci[5] >= var_theo){
      counter[j] <- counter[j] + 1 # true value within limits
    }
  }
}

plot(counter,xaxt = "n",)
axis(1,at=seq(1,10),labels = seq(10,100,by=10))


## STUDENTIZED CONFIDENCE INTERVAL

fun_stud <- function(y, id) {
  boot.out <- boot(y, fun_var, 50, 
                   parallel="multicore", ncpus = no_cores)
  return (c(var(y[id]), var(boot.out$t)))
}

start_time <- Sys.time() # Start measuring time
N <- 100 # number of samples to draw (small bc of computation time)
counter <- rep(0,10)
for (i in 1:N){
  for (j in 1:10){
    dataset.small <- rgamma(j*10, shape = k, scale = gamma)
    boot.out <- boot(dataset.small, fun_stud, 1000, 
                     parallel="multicore", ncpus = no_cores)
    ci <- boot.ci(boot.out, conf = 0.95, type = "stud")$stud
    if (ci[4] <= var_theo && ci[5] >= var_theo){
      counter[j] <- counter[j] + 1 # true value within limits
    }
  }
}

plot(counter,xaxt = "n",)
axis(1,at=seq(1,10),labels = seq(10,100,by=10))


end_time <- Sys.time() # End of measured time
cpu_time <- end_time - start_time
cpu_time









