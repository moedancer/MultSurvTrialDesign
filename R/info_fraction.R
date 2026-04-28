# This code is based on the code provided to the publication of Erdmann et al. (2024), referenced in the main paper
# We use it to get a simulation-based estimate of the expected information fraction for OS at the interim analysis.

library(simIDM)
library(survival)

# Extract key parameters from a set of simulated trials
getEmpPower <- function(
  simStudy,
  nEventPFS,
  nEventOS,
  CriticalPFS = NULL,
  CriticalOS = NULL
) {
  nRep <- length(simStudy)
  ## restrict data first event - PFS
  studyCensored1 <- lapply(
    simStudy,
    censoringByNumberEvents,
    eventNum = nEventPFS,
    typeEvent = "PFS"
  )

  ## restrict data second event -OS
  studyCensored2 <- lapply(
    simStudy,
    censoringByNumberEvents,
    eventNum = nEventOS,
    typeEvent = "OS"
  )
  ## median time
  TimePoints1 <- lapply(
    simStudy,
    getTimePoint,
    eventNum = nEventPFS,
    typeEvent = "PFS",
    byArm = FALSE
  )
  median_time1 <- median(unlist(TimePoints1))

  TimePoints2 <- lapply(
    simStudy,
    getTimePoint,
    eventNum = nEventOS,
    typeEvent = "OS",
    byArm = FALSE
  )
  median_time2 <- median(unlist(TimePoints2))

  events1 <- lapply(
    seq_along(TimePoints1),
    function(t) {
      return(sum(simStudy[[t]]$OSevent[
        (simStudy[[t]]$OStime + simStudy[[t]]$recruitTime) <= TimePoints1[[t]]
      ]))
    }
  )

  events2 <- lapply(
    seq_along(TimePoints2),
    function(t) {
      return(sum(simStudy[[t]]$PFSevent[
        (simStudy[[t]]$PFStime + simStudy[[t]]$recruitTime) <= TimePoints2[[t]]
      ]))
    }
  )

  nOther1 <- mean(unlist(events1))
  nOther2 <- mean(unlist(events2))
  ## significant log-rank test?

  logrank1 <- lapply(studyCensored1, LogRankTest, endpoint = "PFS", CriticalPFS)
  logrank2 <- lapply(studyCensored2, LogRankTest, endpoint = "OS", CriticalOS)

  TestPassed1 <- lapply(logrank1, `[[`, 1)
  TestPassed2 <- lapply(logrank2, `[[`, 1)

  power1 <- 100 * (sum(unlist(TestPassed1)) / nRep)
  power2 <- 100 * (sum(unlist(TestPassed2)) / nRep)

  TestBoth <- (unlist(TestPassed1) + unlist(TestPassed2) == 2)

  powerBoth <- 100 * (sum(TestBoth) / nRep)

  return(c(
    nEventPFS,
    nEventOS,
    power1,
    power2,
    powerBoth,
    nOther1,
    nOther2,
    median_time1,
    median_time2
  ))
}

# Conduct log-rank tests in simulated data
LogRankTest <- function(data, endpoint, Critical) {
  time <- if (endpoint == "OS") {
    data$OStime
  } else if (endpoint == "PFS") {
    data$PFStime
  }
  event <- if (endpoint == "OS") {
    data$OSevent
  } else if (endpoint == "PFS") {
    data$PFSevent
  }
  LogRank <- survdiff(Surv(time, event) ~ trt, data)
  Pvalue <- pchisq(LogRank$chisq, length(LogRank$n) - 1, lower.tail = FALSE)
  Passed <- sqrt(LogRank$chisq) > Critical
  return(list(Passed, Pvalue))
}

## Specify fixed parameters over scenarios
alphaOS <- 0.04
alphaPFS <- 0.01

CriticalOS <- abs(qnorm(alphaOS / 2))
CriticalPFS <- abs(qnorm(alphaPFS / 2))

CriticalOSstart <- CriticalOS
CriticalPFSstart <- CriticalPFS

# Dropout
drprate <- 0.1
drptime <- 12

# Loop over all four main scenarios with w=1 (see main paper for further information)
for (i in 1:4) {
  sec <- i

  if (sec == 1) {
    transitionTrt <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
    transitionPl <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
    nEventOS <- 630
    nEventPFS <- 433
  }

  if (sec == 2) {
    transitionTrt <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
    transitionPl <- exponential_transition(h01 = 0.5, h02 = 0.3, h12 = 0.6)
    nEventOS <- 747
    nEventPFS <- 452
  }

  if (sec == 3) {
    transitionTrt <- exponential_transition(h01 = 0.14, h02 = 0.112, h12 = 0.25)
    transitionPl <- exponential_transition(h01 = 0.18, h02 = 0.15, h12 = 0.255)
    nEventOS <- 742
    nEventPFS <- 644
  }

  if (sec == 4) {
    transitionTrt <- exponential_transition(h01 = 0.18, h02 = 0.06, h12 = 0.17)
    transitionPl <- exponential_transition(h01 = 0.23, h02 = 0.07, h12 = 0.19)
    nEventOS <- 963
    nEventPFS <- 940
  }

  transitionList <- list(transitionPl, transitionTrt)

  SimH1 <- getClinicalTrials(
    nRep = 10000,
    nPat = c(800, 800),
    seed = 1238,
    datType = "1rowPatient",
    transitionByArm = transitionList,
    dropout = list(rate = drprate, time = drptime),
    accrual = list(param = "intensity", value = 25)
  )

  Power <- getEmpPower(SimH1, nEventPFS, nEventOS, CriticalPFS, CriticalOS)
  # Compute information fraction for OS at interim analysis
  cat(Power[6] / Power[2], "\n")
}
