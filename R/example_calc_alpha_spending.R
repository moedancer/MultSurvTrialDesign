library(rpact)
library(mvtnorm)

# Roughly corresponds to the disjunctive power in Scenario 1 with w=1
target_power <- 0.9
alpha <- 0.025

# Non-centrality parameter for the normal distribution
ncp <- qnorm(1 - alpha) + qnorm(target_power)

# Get critical values of alpha-spending approach
design_spending <- getDesignGroupSequential(
  kMax = 2,
  alpha = alpha,
  sided = 1,
  typeOfDesign = "asUser",
  informationRates = c(0.5, 1),
  userAlphaSpending = c(0.005, 0.025)
)

# Define mean vector and covariance matrix for the bivariate normal distribution under alternative
corr_mat <- matrix(c(1, sqrt(0.5), sqrt(0.5), 1), ncol = 2, nrow = 2)
drift_vec <- c(ncp / sqrt(2), ncp)

total_power_spending <- 1 -
  pmvnorm(
    lower = -Inf,
    upper = design_spending$criticalValues,
    mean = drift_vec,
    corr = corr_mat
  )[[1]]

total_power_bon <- 1 -
  pmvnorm(
    lower = -Inf,
    upper = c(qnorm(1 - 0.005), qnorm(1 - 0.02)),
    mean = drift_vec,
    corr = corr_mat
  )[[1]]

total_power_spending - total_power_bon
