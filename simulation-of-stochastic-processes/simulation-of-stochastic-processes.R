# Lab 6 - Simulation of stochastic processes
# Author: Alberto Nieto Sandino
# Date: 2018-10-23

# Clean environment 
closeAllConnections()
rm(list=ls())
# Clear plots
graphics.off()
# Clear history
clearhistory <- function() {
  write("", file=".blank")
  loadhistory(".blank")
  unlink(".blank")
}
clearhistory()

# Example 6.3 - Simulation of a Poisson process
x <- cumsum(rexp(50))
y <- cumsum(c(0,rep(1,50)))
plot(stepfun(x,y),xlim = c(0,10), do.points = F)

# Example 6.4 - Simulation of a sample path of a Poisson process
n <- 100
x <- seq(0,10,length = 1000)
y <- cumsum(rpois(1000,1/n))
plot(x,y)

# Example 6.2 - Simulation of a sample path for a standard Wiener process
h <- 1/1000 # step size
t <- seq(0,1,by = h)
inc <- rnorm(1000, mean = 0, sd = sqrt(h))
W <- c(0,cumsum(inc))
plot(t,W,type='l')

##################
## Assignment 1 ##
##################

# Clean environment 
closeAllConnections()
rm(list=ls())
# Clear plots
graphics.off()

lmd = 1 # rate of exponential distribution
t <- seq(0,10, length = 100)

# Define g(t)
g <- function(a){
  if ((a >= 0) && (a <= (1/2))){ # values of g(t) go from 0 to 0.5
    return(1)
  } else{
    return(0)
  }
}

# Because of the form of function g(t), we need to define xi and eta

# form xi so that cdf(xi) is equal or smaller than 10
xi <- rep(0,20) 
i <- 1
while (cumsum(xi)[i] < 10){
  xi[i] = rexp(1, rate = lmd)
  i = i + 1
}
xi[i] <- 0 

# form eta so that cdf(eta) is equal or smaller than 0.5
eta <- rep(0,10)
j <- 1
while (cumsum(eta)[j] < 0.5){
  eta[j] = rexp(1,rate = lmd)
  j = j + 1
}
eta[j] <- 0

# Construct the function to calculate the sample path
shot_noise <- function(t){
  X <- 0 # initialize X
  for (k in 1:(i-1)){
    X = X + g(t-cumsum(xi)[k]) # sumation of xi components
  }
  for (n in 1:(j-1)){
    X = X + g(t+cumsum(eta)[n]) # sumation of eta components
  }
  return (X) # return the computation of shot noise
}

y <- rep(0,length(t))
for (l in 1:length(t)){
  y[l] <- shot_noise(t[l]) # Sample path for the specified time, t
}

plot(stepfun(t, c(0,y)), xlim = c(0,10), do.points = F,
     xlab = "t", ylab = "X(t)", main = "Simulation of the sample path of X(t)")


##################
## Assignment 2 ##
##################

# Clean environment 
closeAllConnections()
rm(list=ls())
# Clear plots
graphics.off()


####################
## Assignment 2.1 ##
####################

lmd <- 1
n <- 10000

zeta <- rep(0,n) 
for (i in 1:n){ # Perform 10000 simulations
  
  # Constructing xi
  xi <- rep(0,100)
  j <- 1
  while (cumsum(xi)[j]<10){
    xi[j] <- rexp(1, rate = lmd)
    j <- j + 1
  }
  xi[j] <- 0
  
  # Constructing eta
  eta <- rep(0,100)
  k <- 1
  while (cumsum(eta)[k]<0.5){
    eta[k] <- rexp(1,rate = lmd)
    k <- k + 1
  }
  eta[j] <- 0
  
  # Constructing M (eq. 6.3)
  M_1 <- -cumsum (eta)[cumsum(eta) < 0.5] 
  M_2 <-  cumsum (xi) [cumsum(xi)  < 10]
  M   <- sort(c(M_1, M_2))
  
  # Check that condition in eq (6.4) holds for M (TRUE-> zeta = 1)
  if (length(M) > 2){
    for (l in 3:length(M)){
      if (M[l] - M[l-2] <= 0.5){
        zeta[i] <- 1
      }
    }
  }
}

# Expected value of zeta (mean)
p_hat <- mean (zeta)
p_hat
# Standard error of zeta 
sigma_p <- sd(zeta)
se_hat <- sigma_p/sqrt(n)
se_hat
# Calculate the normal CI
alpha <- 0.05
z.lower <- qnorm (alpha/2) # lower limit
z.upper <- qnorm (1-(alpha/2)) # upper limit
ci <- c(p_hat + z.lower*se_hat, 
        p_hat + z.upper*se_hat)
ci

####################
## Assignment 2.2 ##
####################

# Clean environment 
closeAllConnections()
rm(list=ls())
# Clear plots
graphics.off()

hits <- rep(0,10000) # initialize counter for hits
lmd = 1 # rate of exponential distribution
h <- 1/100 # distance between time points
t <- seq(0,10, by = h)

# Define g(t)
g <- function(a){
  if ((a >= 0) && (a <= (1/2))){ # values of g(t) go from 0 to 0.5
    return(1)
  } else{
    return(0)
  }
}

start_time <- Sys.time() # Start measuring time
for (n in 1:10000){
  
  # Because of the form of function g(t), we need to define xi and eta
  # form xi so that cdf(xi) is equal or smaller than 10
  xi <- rep(0,200) 
  i <- 1
  while (cumsum(xi)[i] < 10){
    xi[i] = rexp(1, rate = lmd)
    i = i + 1
  }
  xi[i] <- 0 
  
  # form eta so that cdf(eta) is equal or smaller than 0.5
  eta <- rep(0,100)
  j <- 1
  while (cumsum(eta)[j] < 0.5){
    eta[j] = rexp(1,rate = lmd)
    j = j + 1
  }
  eta[j] <- 0
  
  # Construct the function to calculate the sample path
  shot_noise <- function(t){
    X <- 0 # initialize X
    for (k in 1:(i-1)){
      X = X + g(t-cumsum(xi)[k]) # sumation of xi components
    }
    for (n in 1:(j-1)){
      X = X + g(t+cumsum(eta)[n]) # sumation of eta components
    }
    return (X) # return the computation of shot noise
  }
  
  y <- rep(0,length(t))
  for (l in 1:length(t)){
    y[l] <- shot_noise(t[l]) # Sample path for the specified time, t
    if (y[l] > 2){
      hits[n] <- 1
    }
  }
}

# Expected value of zeta (mean)
p_hat <- mean (hits)
p_hat
# Standard error of zeta 
sigma_p <- sd(hits)
se_hat <- sigma_p/sqrt(n)
se_hat
# Calculate the normal CI
alpha <- 0.05
z.lower <- qnorm (alpha/2) # lower limit
z.upper <- qnorm (1-(alpha/2)) # upper limit
ci <- c(p_hat + z.lower*se_hat, 
        p_hat + z.upper*se_hat)
ci

end_time <- Sys.time() # End of measured time
cpu_time <- end_time - start_time
cpu_time


##################
## Assignment 3 ##
##################

# Clean environment 
closeAllConnections()
rm(list=ls())
# Clear plots
graphics.off()

# We start by defining f(t)
f <- function(t){
  fun <- 1 - t^2
  fun [abs(t) > 1] <- 0
  return (fun)
}

n <- 1000 # number of increments
t <- 10
h <- t/n # distance between time points
x <- seq(0,t, by = h)
s_n <- n^2 # Define the summation limits

# Define function to apporximate a stationary Gaussian process X(t)
gaussian_process <- function (t){
  res <- 0
  for (k in (-s_n):s_n){
    inc <- rnorm (1, mean = 0, sd = sqrt(h))
    res = res + f(t+k/n)*inc
  }
  return (res)
}

y = gaussian_process(x)
plot(x, y, xlim = c(0,10), type = 'l',
     xlab = "t", ylab = "X(t)", main = "Simulation of the sample path of X(t)")





