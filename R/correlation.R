require(simIDM)

# Compute correlation coefficient for PFS and OS in time-homogeneous Markovian illness-death models
# Formula according to Fleischer et al. (2009)
# We work on the basis of the simIDM packages
thm_correlation <- function(simIDM_thm_model) {
  stopifnot(all(
    class(simIDM_thm_model) ==
      c("ExponentialTransition", "TransitionParameters")
  ))
  numerator <- simIDM_thm_model$hazards$h12
  denominator_sq <- simIDM_thm_model$hazards$h01^2 +
    2 * simIDM_thm_model$hazards$h01 * simIDM_thm_model$hazards$h02 +
    simIDM_thm_model$hazards$h12^2
  return(numerator / sqrt(denominator_sq))
}

### Scenario 1
transitionTrt <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
transitionPl <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
thm_correlation(transitionTrt)
thm_correlation(transitionPl)

### Scenario 2
transitionTrt <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
transitionPl <- exponential_transition(h01 = 0.5, h02 = 0.3, h12 = 0.6)
thm_correlation(transitionTrt)
thm_correlation(transitionPl)

### Scenario 3
transitionTrt <- exponential_transition(h01 = 0.14, h02 = 0.112, h12 = 0.25)
transitionPl <- exponential_transition(h01 = 0.18, h02 = 0.15, h12 = 0.255)
thm_correlation(transitionTrt)
thm_correlation(transitionPl)

### Scenario 4
transitionTrt <- exponential_transition(h01 = 0.18, h02 = 0.06, h12 = 0.17)
transitionPl <- exponential_transition(h01 = 0.23, h02 = 0.07, h12 = 0.19)
thm_correlation(transitionTrt)
thm_correlation(transitionPl)
