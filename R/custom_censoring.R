### We need this addditional function to simulate censoring because we cannot combine
###   the censoring from the original simIDM functions with the new frailty
custom_censoring <- function(simIDM_widedata, drprate, drptime) {
  n <- dim(simIDM_widedata)[1]

  std_rate <- -log(1 - drprate) / drptime
  sim_dropout <- rexp(n, std_rate)

  simIDM_widedata$CensoredPFS <- ifelse(
    sim_dropout < simIDM_widedata$PFStime,
    1,
    0
  )
  simIDM_widedata$PFSevent <- 1 - simIDM_widedata$CensoredPFS
  simIDM_widedata$PFStime <- pmin(simIDM_widedata$PFStime, sim_dropout)
  simIDM_widedata$PFStimeCal <- simIDM_widedata$PFStime +
    simIDM_widedata$recruitTime

  simIDM_widedata$CensoredOS <- ifelse(
    sim_dropout < simIDM_widedata$OStime,
    1,
    0
  )
  simIDM_widedata$OSevent <- 1 - simIDM_widedata$CensoredOS
  simIDM_widedata$OStime <- pmin(simIDM_widedata$OStime, sim_dropout)
  simIDM_widedata$OStimeCal <- simIDM_widedata$OStime +
    simIDM_widedata$recruitTime

  return(simIDM_widedata)
}
