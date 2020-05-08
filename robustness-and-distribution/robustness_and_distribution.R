# Lab 1 by Alberto Nieto Sandino

################################################
# Assignment 1.1
################################################

# Clear environment 
rm(list=ls())
# Close graphics
graphics.off()
# Packages to use
library('zoo')
library('nortest')

# Sample 100 observations from a gamma distribution
num_of_samples = 100
alpha = 2
beta = 2
samples <- rgamma (num_of_samples, shape = alpha, scale = beta)
mean_samp = alpha*beta
var_samp = alpha*beta^2
# Plot a histogram to get an overview of the data
dev.new()
p1 <- hist(samples, breaks = 50, right = FALSE)

# qq-plot
dev.new()
qqnorm(samples, main = "Normal Q-Q Plot",
       xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
       plot.it = TRUE, datax = FALSE)
qqline(samples, col = 2)

# Form the reference distribution, normal distribution
breaks_cdf <- pnorm(p1$breaks, mean = mean_samp, sd = sqrt(var_samp)) # since we assume it's normal 
# distribution. Do we use the normal mean and sd, or do we put the real mean and variance from the
# gamma distribution
null.probs <- rollapply(breaks_cdf, 2, function(x) x[2]-x[1])

## Formal goodness of fit tests
# Chi square test
a <- chisq.test(p1$counts, p = null.probs, rescale.p = TRUE, simulate.p.value = TRUE)

# Lilliefors test
lillie.test(samples)

################################################
# Assignment 1.2
################################################

# Sample 100 observations from a gamma distribution
num_of_samples = 100
alpha = 2
beta = 2
samples <- rgamma (num_of_samples, shape = alpha, scale = beta)
mean_samp = alpha*beta
var_samp = alpha*beta^2

lower <- vector()
upper <- vector()
width <- vector()
counter <- as.integer(0)

for (i in 1:100){
  # 95% confidence interval for the variance
  df <- length(samples) - 1
  var_samp_temp = var(samples) # we don't know it's a gamma distribution
  lower[i] = var_samp_temp * df / qchisq(0.05/2, df, lower.tail = FALSE)
  upper[i] = var_samp_temp * df / qchisq(1 - 0.05/2, df, lower.tail = FALSE)
  width[i] = upper[i] - lower[i]
  if ((lower[i] < var_samp) & (upper > var_samp)){
    counter = counter + 1
  } else {
    counter = counter
  }
  # Generate new gamma-distributed samples for next loop
  samples <- rgamma (num_of_samples, shape = alpha, scale = beta)
}















