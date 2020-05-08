# Lab 03 - Reliability and survival
# Author: Alberto Nieto Sandino
# Date: 2018-11-02

library(numDeriv)
library(ggplot2)
library(alabama)

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

# Assigment 1 (4p)

R_T <- function(t){
  f1 <- ((1-(1-exp(-sqrt(t)))^3)*(exp(-t/2)))
  return (f1)
}

# Expected life length
E_T = integrate(R_T, lower = 0, upper = Inf)$value
E_T

loga <- function(t){
  f1 <- log((1-(1-exp(-sqrt(t)))^3)*(exp(-t/2)))
  return (f1)
}

# Death rate, r_T
x <- seq(0.001,10,0.001)
r_T <- -grad(loga,x)
plot(r_T, type="l", xaxt = "n", xlab='Time')
axis(1, at=seq(1,10000,1000), labels=1:10)


# Assigment 1.2

integrand <- function(t){
  f1 <- (1-(1-exp(-sqrt(t)))^3)*(exp(-t/2)/2)
  return (f1)
}
P_4 <- integrate(integrand, lower = 0, upper = Inf)$value
P_4


# Assignment 2 (3p)

# Assignment 2.1

############################################################
## WARM COMPONENT ##########################################
############################################################

# Survival function
R_T_warm <- function(t){
  f1 <- ((1-(1-exp(-sqrt(t)))^3)*(1-(1-(exp(-t/2)))^2))
  return (f1)
}

# Expected life length with warm component
E_T_warm = integrate(R_T_warm, lower = 0, upper = Inf)$value
E_T_warm

# Death rate, r_T_warm
loga_warm <- function(t){
  f1 <- log((1-(1-exp(-sqrt(t)))^3)*(1-(1-(exp(-t/2)))^2))
  return (f1)
}

x <- seq(0.001,10,0.001)
r_T_warm <- -grad(loga_warm,x)
plot(r_T_warm, type="l", xaxt = "n", xlab='Time')
axis(1, at=seq(1,10000,1000), labels=1:10)

############################################################
## COLD COMPONENT ##########################################
############################################################

# Calculate the doble integral inside R_T_cold
g <- function(t){
  f <- function(x) {(1-exp(-(1/2)*(t-x)))*(1/2)*exp(-(1/2)*x)}
  integrate(f,0,t)$value
}
int <- Vectorize(g)

# Survival function
R_T_cold<- function(t){
  f1 <- ((1-(1-exp(-sqrt(t)))^3)*(1-int(t)))
  return (f1)
}

# Expected life length with cold component
E_T_cold = integrate(R_T_cold, lower = 0, upper = Inf)$value
E_T_cold

# Death rate, r_T_cold
loga_cold <- function(t){
  f1 <- log((1-(1-exp(-sqrt(t)))^3)*(1-int(t)))
  return (f1)
}

x <- seq(0.001,10,0.001)
r_T_cold <- -grad(loga_cold,x)
plot(r_T_cold, type="l", xaxt = "n", xlab='Time')
axis(1, at=seq(1,10000,1000), labels=1:10)


# Plot of the death rates for both cases and the one from task 1
lim_range <- range(0, r_T, r_T_warm, r_T_cold) # calculate range of values
plot(r_T, type="l", col ="blue", ylim = lim_range)
lines(r_T_warm, type = "l", pch = 22, lty = 2, col = "red")
lines(r_T_cold, type = "l", pch = 23, lty = 3, col = "black")
legend("bottomright", legend = c("r_T", "r_T_warm", "r_T_cold"), cex=0.8, 
       col = c("blue","red","black"), lty = 1:3)



## Assignment 2.2

######################
## WARM COMPONENT ####
######################

# x is rho
E_T_rho <- function(x){
  integrate(
    function(t){
      ((1-(1-exp(-sqrt(t)))^3)*(exp(-t*x)))
    }
    ,lower = 0, upper = Inf)$value
}
E_T_rho_v <- Vectorize(E_T_rho)

# Find the root for the warm redundancy
equation <- function(x){
  E_T_rho_v(x)-E_T_warm
}

rho_warm <- uniroot(equation,interval = c(0.0,0.5))$root
rho_warm

######################
## COLD COMPONENT ####
######################

# x is rho
E_T_rho <- function(x){
  integrate(
    function(t){
      ((1-(1-exp(-sqrt(t)))^3)*(exp(-t*x)))
    }
    ,lower = 0, upper = Inf)$value
}
E_T_rho_v <- Vectorize(E_T_rho)

# Find the root for the warm redundancy
equation <- function(x){
  E_T_rho_v(x)-E_T_cold
}

rho_cold <- uniroot(equation,interval = c(0.0,0.5))$root
rho_cold

########################################################
## Assignment 3 ########################################
########################################################

# Remember that we do minimization -> times -1

# Function for the expected life length
E_T_var <- function(x){
  integrate(
    function(t){
      (-1*(exp(-(x[1]*t)^1/3)*(1-(1-exp(-(x[2]*t)^(1/3)))^2)))
    }
    ,lower = 0, upper = Inf, rel.tol = 1e-10)$value
}

# Perform non-linear optimization with mu,lambda > 0 as constraints
solution <- rep(NA,10)
fval <- rep(NA,10)
mu <- rep(NA,10)
lambda <- rep(NA,10)

for (c in 1:10){
  p0 <- c(1,1)
  # Inequalities: mu, lambda
  hin <- function(x) {
    h <- rep(NA, 1)
    h[1] <- x[1]
    h[2] <- x[2]
    h
  }
  
  # Equalities: cost
  heq <- function(x) {
    h <- rep(NA, 1)
    h[1] <- ((1/5) + (1/x[1]) + 2*((1/5) + (1/x[2]))) - c
    h
  }
  
  solution <- constrOptim.nl(par=p0, fn=E_T_var, heq=heq, hin=hin) 
  fval[c] <- -1*solution$value # reverse the change for minimization
  mu[c] <- solution$par[1] # x[1] = mu
  lambda[c] <- solution$par[2] # x[2] = lambda
}


# Plot of mu and lambda
lim_range <- range(0, mu, lambda) # calculate range of values
plot(mu, type="l", col ="blue", ylim = lim_range, ylab ="Parameters",
     xlab = "Values for the total cost, c")
lines(lambda, type = "l", pch = 22, lty = 2, col = "red")
legend("topright", legend = c("mu", "lambda"), cex=0.8, 
       col = c("blue","red"), lty = 1:2)


# Plot of the function value
lim_range <- range(0, fval) # calculate range of values
plot(fval, type="l", col ="blue", ylim = lim_range, ylab ="Expected life length",
     xlab = "Values for the total cost, c")

