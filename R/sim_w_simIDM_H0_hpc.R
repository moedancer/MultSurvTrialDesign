hpc_version <- FALSE

if (hpc_version) {
  .libPaths("/home/d/danzerm/R/library")
}

require(foreach)
require(doParallel)

# Setup for parallel simulation
cores_to_use <- 20
my_cluster <- makeCluster(cores_to_use)
registerDoParallel(my_cluster)
if (hpc_version) {
  clusterEvalQ(my_cluster, .libPaths("/home/d/danzerm/R/library"))
}

secs <- 1:4

# Set recruitment rates and set factor to determine maximum sample size per arm
rr_vec <- c(2, 4, 10, 15, 25)
rr_to_n <- 32

# Determine event rates of PFS and OS after which the analyses shall be conducted
r_PFS <- 25 / 64
r_OS <- 38 / 64

combinations <- expand.grid(secs, rr_vec)

# Shape of gamma frailty distribution
# Scale will be chosen to give expected value of 1
shape_gamma_frailty <- 10
scale_gamma_frailty <- 1 / shape_gamma_frailty

sim_num <- 100000

power_metrics <- foreach(c = 1:cores_to_use, .errorhandling = "pass") %dopar%
  {
    if (hpc_version) {
      .libPaths("/home/d/danzerm/R/library")
    }

    source("stat_and_var_calc.R")
    source("custom_censoring.R")

    require(mvtnorm)
    require(simIDM)
    require(dplyr)
    require(rpact)

    sec <- combinations[c, 1]
    # Choose recruitment rate for this core
    rr <- combinations[c, 2]

    test_names <- c("PFS_a1", "OS_a1", "PFS_a2", "OS_a2")
    num_tests <- length(test_names)
    strategy_names <- c(
      "standard",
      "pass_on",
      "shift_to_last",
      "shift_at_interim",
      "standard_gs",
      "gs_pass_on",
      "gs_shift_to_last",
      "gs_shift_at_interim",
      "final_os"
    )
    num_strategies <- length(strategy_names)
    metrics_names <- c("Rej_PFS", "Rej_OS", "Rej_One", "Rej_Both", "Early_Stop")
    num_metrics <- length(metrics_names)

    # Results without frailty
    test_stats <- matrix(NA, ncol = num_tests, nrow = sim_num)
    var_ests <- array(NA, dim = c(sim_num, num_tests, num_tests))
    p_values <- matrix(NA, ncol = num_tests, nrow = sim_num)
    decisions <- array(0, dim = c(sim_num, num_strategies, num_tests))
    results <- array(0, dim = c(sim_num, num_strategies, num_metrics))
    power_metrics <- matrix(NA, nrow = num_strategies, ncol = num_metrics)

    # Results with frailty
    test_stats_frailty <- matrix(NA, ncol = num_tests, nrow = sim_num)
    var_ests_frailty <- array(NA, dim = c(sim_num, num_tests, num_tests))
    p_values_frailty <- matrix(NA, ncol = num_tests, nrow = sim_num)
    decisions_frailty <- array(0, dim = c(sim_num, num_strategies, num_tests))
    results_frailty <- array(0, dim = c(sim_num, num_strategies, num_metrics))
    power_metrics_frailty <- matrix(
      NA,
      nrow = num_strategies,
      ncol = num_metrics
    )

    colnames(test_stats) <-
      dimnames(var_ests)[[2]] <- dimnames(var_ests)[[3]] <-
        colnames(p_values) <-
          dimnames(decisions)[[3]] <- test_names
    dimnames(decisions)[[2]] <- dimnames(results)[[2]] <-
      rownames(power_metrics) <- strategy_names
    dimnames(results)[[3]] <- colnames(power_metrics) <- metrics_names

    colnames(test_stats_frailty) <-
      dimnames(var_ests_frailty)[[2]] <- dimnames(var_ests_frailty)[[3]] <-
        colnames(p_values_frailty) <-
          dimnames(decisions_frailty)[[3]] <- test_names
    dimnames(decisions_frailty)[[2]] <- dimnames(results_frailty)[[2]] <-
      rownames(power_metrics_frailty) <- strategy_names
    dimnames(results_frailty)[[3]] <- colnames(
      power_metrics_frailty
    ) <- metrics_names

    if (sec == 1) {
      transitionTrt <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
      transitionPl <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
      nEventPFS <- rr * rr_to_n * 2 * r_PFS
      nEventOS <- rr * rr_to_n * 2 * r_OS
    }

    if (sec == 2) {
      transitionTrt <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
      transitionPl <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
      nEventPFS <- rr * rr_to_n * 2 * r_PFS
      nEventOS <- rr * rr_to_n * 2 * r_OS
    }

    if (sec == 3) {
      transitionTrt <- exponential_transition(
        h01 = 0.14,
        h02 = 0.112,
        h12 = 0.25
      )
      transitionPl <- exponential_transition(
        h01 = 0.14,
        h02 = 0.112,
        h12 = 0.25
      )
      nEventPFS <- rr * rr_to_n * 2 * r_PFS
      nEventOS <- rr * rr_to_n * 2 * r_OS
    }

    if (sec == 4) {
      transitionTrt <- exponential_transition(
        h01 = 0.18,
        h02 = 0.06,
        h12 = 0.17
      )
      transitionPl <- exponential_transition(h01 = 0.18, h02 = 0.06, h12 = 0.17)
      nEventPFS <- as.integer(rr * rr_to_n * 2 * r_PFS)
      nEventOS <- as.integer(rr * rr_to_n * 2 * r_OS)
    }
    transitionList <- list(transitionPl, transitionTrt)

    drprate <- 0.1
    drptime <- 12

    ### SET UP CLOSED TESTING AND GROUP-SEQUENTIAL PROCEDURES
    bretz_weight_pfs <- 1 / 5
    bretz_weight_os <- 4 / 5

    alpha <- 0.025

    # Group-sequential design for OS at uncorrected level
    os_gs_design <- getDesignGroupSequential(
      kMax = 2,
      alpha = alpha,
      sided = 1,
      typeOfDesign = "asOF",
      informationRates = c(0.4, 1)
    )
    os_stageLevels <- os_gs_design$stageLevels
    os_alphaSpent <- os_gs_design$alphaSpent
    os_alphaIncr <- (c(os_alphaSpent, alpha * bretz_weight_os) -
      c(0, os_alphaSpent))[1:2]
    os_criticalValues <- os_gs_design$criticalValues

    # Group-sequential design for OS at corrected level
    os_gs_design_corrected <- getDesignGroupSequential(
      kMax = 2,
      alpha = alpha * bretz_weight_os,
      sided = 1,
      typeOfDesign = "asOF",
      informationRates = c(0.4, 1)
    )
    os_stageLevels_corrected <- os_gs_design_corrected$stageLevels
    os_alphaSpent_corrected <- os_gs_design_corrected$alphaSpent
    os_alphaIncr_corrected <- (c(
      os_alphaSpent_corrected,
      alpha * bretz_weight_os
    ) -
      c(0, os_alphaSpent_corrected))[1:2]
    os_criticalValues_corrected <- os_gs_design_corrected$criticalValues

    # Group-sequential design for OS with additional weight shifted to final analysis
    os_gs_design_shiftLast <- getDesignGroupSequential(
      kMax = 2,
      alpha = alpha,
      sided = 1,
      typeOfDesign = "asUser",
      userAlphaSpending = c(
        os_alphaSpent_corrected[1],
        os_alphaSpent_corrected[2] + alpha * bretz_weight_pfs
      ),
      informationRates = c(0.4, 1)
    )
    os_stageLevels_shiftLast <- os_gs_design_shiftLast$stageLevels
    os_alphaSpent_shiftLast <- os_gs_design_shiftLast$alphaSpent
    os_alphaIncr_shiftLast <- (c(os_alphaSpent_shiftLast, alpha) -
      c(0, os_alphaSpent_shiftLast))[1:2]
    os_criticalValues_shiftLast <- os_gs_design_shiftLast$criticalValues

    # Group-sequential design for OS with additional weight shifted to interim analysis
    os_gs_design_shiftInterim <- getDesignGroupSequential(
      kMax = 2,
      alpha = alpha,
      sided = 1,
      typeOfDesign = "asUser",
      userAlphaSpending = os_alphaSpent_corrected + alpha * bretz_weight_pfs,
      informationRates = c(0.4, 1)
    )
    os_stageLevels_shiftInterim <- os_gs_design_shiftInterim$stageLevels
    os_alphaSpent_shiftInterim <- os_gs_design_shiftInterim$alphaSpent
    os_alphaIncr_shiftInterim <- (c(os_alphaSpent_shiftInterim, alpha) -
      c(0, os_alphaSpent_shiftInterim))[1:2]
    os_criticalValues_shiftInterim <- os_gs_design_shiftInterim$criticalValues

    for (i in seq(sim_num)) {
      Sim <- getDatasetWideFormat(getOneClinicalTrial(
        nPat = c(rr * rr_to_n, rr * rr_to_n), #seed = 181993,
        transitionByArm = transitionList,
        dropout = list(rate = drprate, time = drptime),
        accrual = list(param = "intensity", value = rr)
      ))

      # Use my own censoring function for same workflow as for frailty data
      # (Frailty data not compatible with censoring from simIDM)
      Sim <- custom_censoring(Sim, drprate = drprate, drptime = drptime)

      studyCensored1 <- censoringByNumberEvents(
        data = Sim,
        eventNum = nEventPFS,
        typeEvent = "PFS"
      )
      studyCensored2 <- censoringByNumberEvents(
        data = Sim,
        eventNum = nEventOS,
        typeEvent = "OS"
      )

      calDatePFS <- getTimePoint(
        data = Sim,
        eventNum = nEventPFS,
        typeEvent = "PFS"
      )
      calDateOS <- getTimePoint(
        data = Sim,
        eventNum = nEventOS,
        typeEvent = "OS"
      )

      complete_data_temp <- merge(
        studyCensored1,
        studyCensored2,
        by = "id",
        all = TRUE,
        suffixes = c("_a1", "_a2")
      )
      complete_data_temp <- complete_data_temp[, c(
        "id",
        "trt_a2",
        "PFStime_a1",
        "PFSevent_a1",
        "OStime_a1",
        "OSevent_a1",
        "PFStime_a2",
        "PFSevent_a2",
        "OStime_a2",
        "OSevent_a2"
      )]
      complete_data_temp[is.na(complete_data_temp)] <- 0
      names(complete_data_temp)[2] <- "group"
      complete_data_temp$group <- complete_data_temp$group - 1

      # Create data with patient-specific frailty (by individual rescaling of time)
      Sim_frailty <- Sim
      frailty <- rgamma(
        n = rr * rr_to_n * 2,
        scale = scale_gamma_frailty,
        shape = shape_gamma_frailty
      )
      Sim_frailty$PFStime <- Sim_frailty$PFStime / frailty
      Sim_frailty$OStime <- Sim_frailty$OStime / frailty
      Sim_frailty$PFStimeCal <- Sim_frailty$recruitTime + Sim_frailty$PFStime
      Sim_frailty$OStimeCal <- Sim_frailty$recruitTime + Sim_frailty$OStime

      # Use my own censoring function for same workflow as for frailty data
      # (Frailty data not compatible with censoring from simIDM)
      Sim_frailty <- custom_censoring(
        Sim_frailty,
        drprate = drprate,
        drptime = drptime
      )

      studyCensored1_frailty <- censoringByNumberEvents(
        data = Sim_frailty,
        eventNum = nEventPFS,
        typeEvent = "PFS"
      )
      studyCensored2_frailty <- censoringByNumberEvents(
        data = Sim_frailty,
        eventNum = nEventOS,
        typeEvent = "OS"
      )

      calDatePFS_frailty <- getTimePoint(
        data = Sim_frailty,
        eventNum = nEventPFS,
        typeEvent = "PFS"
      )
      calDateOS_frailty <- getTimePoint(
        data = Sim_frailty,
        eventNum = nEventOS,
        typeEvent = "OS"
      )

      complete_data_temp_frailty <- merge(
        studyCensored1_frailty,
        studyCensored2_frailty,
        by = "id",
        all = TRUE,
        suffixes = c("_a1", "_a2")
      )
      complete_data_temp_frailty <- complete_data_temp_frailty[, c(
        "id",
        "trt_a2",
        "PFStime_a1",
        "PFSevent_a1",
        "OStime_a1",
        "OSevent_a1",
        "PFStime_a2",
        "PFSevent_a2",
        "OStime_a2",
        "OSevent_a2"
      )]
      complete_data_temp_frailty[is.na(complete_data_temp_frailty)] <- 0
      names(complete_data_temp_frailty)[2] <- "group"
      complete_data_temp_frailty$group <- complete_data_temp_frailty$group - 1

      # Calculate test statistics and covariance estimates for data with and without frailty
      study_result <- test_stat_calc(complete_data_temp)
      test_stats[i, ] <- study_result[[1]]
      var_ests[i, , ] <- study_result[[2]]
      p_values[i, ] <- pnorm(study_result[[1]] / sqrt(diag(study_result[[2]])))

      study_result_frailty <- test_stat_calc(complete_data_temp_frailty)
      test_stats_frailty[i, ] <- study_result_frailty[[1]]
      var_ests_frailty[i, , ] <- study_result_frailty[[2]]
      p_values_frailty[i, ] <- pnorm(
        study_result_frailty[[1]] / sqrt(diag(study_result_frailty[[2]]))
      )

      ##### DETERMINATION OF DECISIONS FOR DIFFERENT TESTING STRATEGIES #####

      ### TESTING STRATEGY 1 ###
      # Test of PFS at interim analysis (level 0.005)
      # Test of OS at final analysis (level 0.02)
      # No propagation
      decisions[i, 1, "PFS_a1"] <- (p_values[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      decisions[i, 1, "OS_a2"] <- (p_values[i, "OS_a2"] <=
        alpha * bretz_weight_os)

      decisions_frailty[i, 1, "PFS_a1"] <- (p_values_frailty[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      decisions_frailty[i, 1, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
        alpha * bretz_weight_os)

      ### TESTING STRATEGY 2 ###
      # Test of PFS at interim analysis (level 0.005)
      # Test of OS at final analysis (level 0.02)
      # with propagation
      decisions[i, 2, "PFS_a1"] <- (p_values[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      if (p_values[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions[i, 2, "OS_a2"] <- (p_values[i, "OS_a2"] <= alpha)
      } else {
        decisions[i, 2, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          alpha * bretz_weight_os)
      }

      decisions_frailty[i, 2, "PFS_a1"] <- (p_values_frailty[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      if (p_values_frailty[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions_frailty[i, 2, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          alpha)
      } else {
        decisions_frailty[i, 2, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          alpha * bretz_weight_os)
      }

      ### TESTING STRATEGY 3 ###
      # Test of PFS at interim analysis (level 0.005)
      # Test of OS at final analysis:
      #  - Exploit dependence with PFS if PFS not rejected
      #  - Propagate level of PFS if PFS rejected
      if (p_values[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions[i, 3, "PFS_a1"] <- 1
        decisions[i, 3, "OS_a2"] <- (p_values[i, "OS_a2"] <= alpha)
      } else {
        temp_factor <- inflation_factor(
          current_levels = alpha * bretz_weight_os,
          previous_levels = alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("PFS_a1", "OS_a2"),
            c("PFS_a1", "OS_a2")
          ]
        )
        decisions[i, 3, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          alpha * bretz_weight_os * temp_factor)
      }

      if (p_values_frailty[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions_frailty[i, 3, "PFS_a1"] <- 1
        decisions_frailty[i, 3, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          alpha)
      } else {
        temp_factor <- inflation_factor(
          current_levels = alpha * bretz_weight_os,
          previous_levels = alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("PFS_a1", "OS_a2"),
            c("PFS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 3, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          alpha * bretz_weight_os * temp_factor)
      }

      ### TESTING STRATEGY 4 ###
      # Test of PFS at interim analysis (level 0.005)
      # Test of OS:
      #  - at interim and final analysis if PFS successful
      #    (spend 0.005 at interim, rest at final)
      #  - only at final if PFS not successful
      #    (exploit dependence)
      if (p_values[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions[i, 4, "PFS_a1"] <- 1
        decisions[i, 4, "OS_a1"] <- (p_values[i, "OS_a1"] <=
          alpha * bretz_weight_pfs)
        temp_factor <- inflation_factor(
          current_levels = alpha * bretz_weight_os,
          previous_levels = alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("OS_a1", "OS_a2"),
            c("OS_a1", "OS_a2")
          ]
        )
        decisions[i, 4, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          alpha * bretz_weight_os * temp_factor)
      } else {
        temp_factor <- inflation_factor(
          current_levels = alpha * bretz_weight_os,
          previous_levels = alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("PFS_a1", "OS_a2"),
            c("PFS_a1", "OS_a2")
          ]
        )
        decisions[i, 4, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          alpha * bretz_weight_os * temp_factor)
      }

      if (p_values_frailty[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions_frailty[i, 4, "PFS_a1"] <- 1
        decisions_frailty[i, 4, "OS_a1"] <- (p_values_frailty[i, "OS_a1"] <=
          alpha * bretz_weight_pfs)
        temp_factor <- inflation_factor(
          current_levels = alpha * bretz_weight_os,
          previous_levels = alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("OS_a1", "OS_a2"),
            c("OS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 4, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          alpha * bretz_weight_os * temp_factor)
      } else {
        temp_factor <- inflation_factor(
          current_levels = alpha * bretz_weight_os,
          previous_levels = alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("PFS_a1", "OS_a2"),
            c("PFS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 4, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          alpha * bretz_weight_os * temp_factor)
      }

      ### TESTING STRATEGY 5 ###
      # Test of PFS at interim analysis (level  0.005)
      # Group-sequential test of OS:
      #   Bounds as specified in Erdmann et al.
      decisions[i, 5, "PFS_a1"] <- (p_values[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      decisions[i, 5, "OS_a1"] <- (p_values[i, "OS_a1"] <=
        os_stageLevels_corrected[1])
      decisions[i, 5, "OS_a2"] <- (p_values[i, "OS_a2"] <=
        os_stageLevels_corrected[2])

      decisions_frailty[i, 5, "PFS_a1"] <- (p_values_frailty[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      decisions_frailty[i, 5, "OS_a1"] <- (p_values_frailty[i, "OS_a1"] <=
        os_stageLevels_corrected[1])
      decisions_frailty[i, 5, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
        os_stageLevels_corrected[2])

      ### TESTING STRATEGY 6 ###
      # Test of PFS at interim analysis (level  0.005)
      # Group-sequential test of OS:
      #   Bounds as specified in Erdmann et al.
      # Propagation from PFS to OS and OS to PFS possible
      if (p_values[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions[i, 6, "PFS_a1"] <- TRUE
        decisions[i, 6, "OS_a1"] <- (p_values[i, "OS_a1"] <= os_stageLevels[1])
        decisions[i, 6, "OS_a2"] <- (p_values[i, "OS_a2"] <= os_stageLevels[2])
      } else if (p_values[i, "OS_a1"] <= os_stageLevels_corrected[1]) {
        decisions[i, 6, "OS_a1"] <- TRUE
        decisions[i, 6, "PFS_a1"] <- (p_values[i, "PFS_a1"] <= alpha)
      } else {
        decisions[i, 6, "OS_a1"] <- (p_values[i, "OS_a1"] <=
          os_stageLevels_corrected[1])
        decisions[i, 6, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          os_stageLevels_corrected[2])
      }

      if (p_values_frailty[i, "PFS_a1"] <= alpha * bretz_weight_pfs) {
        decisions_frailty[i, 6, "PFS_a1"] <- TRUE
        decisions_frailty[i, 6, "OS_a1"] <- (p_values_frailty[i, "OS_a1"] <=
          os_stageLevels[1])
        decisions_frailty[i, 6, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          os_stageLevels[2])
      } else if (p_values_frailty[i, "OS_a1"] <= os_stageLevels_corrected[1]) {
        decisions_frailty[i, 6, "OS_a1"] <- TRUE
        decisions_frailty[i, 6, "PFS_a1"] <- (p_values_frailty[i, "PFS_a1"] <=
          alpha)
      } else {
        decisions_frailty[i, 6, "OS_a1"] <- (p_values_frailty[i, "OS_a1"] <=
          os_stageLevels_corrected[1])
        decisions_frailty[i, 6, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          os_stageLevels_corrected[2])
      }

      ### TESTING STRATEGY 7 ###
      # Test of PFS at interim analysis (level  0.005)
      # Group-sequential test of OS:
      #   Increase alpha to be spent at final analysis by 0.005 if PFS rejected
      temp_factor_interim <- inflation_factor(
        current_levels = c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
        available_level = alpha * bretz_weight_pfs + os_alphaIncr_corrected[1],
        covariance_matrix = var_ests[
          i,
          c("PFS_a1", "OS_a1"),
          c("PFS_a1", "OS_a1")
        ]
      )
      if (
        p_values[i, "PFS_a1"] <= alpha * bretz_weight_pfs * temp_factor_interim
      ) {
        decisions[i, 7, "PFS_a1"] <- 1
        decisions[i, 7, "OS_a1"] <- (p_values[i, "OS_a1"] <=
          os_stageLevels_corrected[1])
        temp_factor <- inflation_factor(
          current_levels = alpha - os_stageLevels_corrected[1],
          previous_levels = os_stageLevels_corrected[1],
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("OS_a1", "OS_a2"),
            c("OS_a1", "OS_a2")
          ]
        )
        decisions[i, 7, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          temp_factor * (alpha - os_stageLevels_corrected[1]))
      } else if (
        p_values[i, "OS_a1"] <= os_alphaIncr_corrected[1] * temp_factor_interim
      ) {
        decisions[i, 7, "PFS_a1"] <- (p_values[i, "PFS_a1"] <= alpha)
      } else {
        temp_factor <- inflation_factor(
          current_levels = os_alphaIncr_corrected[2],
          previous_levels = temp_factor_interim *
            c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("PFS_a1", "OS_a1", "OS_a2"),
            c("PFS_a1", "OS_a1", "OS_a2")
          ]
        )
        decisions[i, 7, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          temp_factor * os_alphaIncr_corrected[2])
      }

      temp_factor_interim <- inflation_factor(
        current_levels = c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
        available_level = alpha * bretz_weight_pfs + os_alphaIncr_corrected[1],
        covariance_matrix = var_ests_frailty[
          i,
          c("PFS_a1", "OS_a1"),
          c("PFS_a1", "OS_a1")
        ]
      )
      if (
        p_values_frailty[i, "PFS_a1"] <=
          alpha * bretz_weight_pfs * temp_factor_interim
      ) {
        decisions_frailty[i, 7, "PFS_a1"] <- 1
        decisions_frailty[i, 7, "OS_a1"] <- (p_values_frailty[i, "OS_a1"] <=
          os_stageLevels_corrected[1])
        temp_factor <- inflation_factor(
          current_levels = alpha - os_stageLevels_corrected[1],
          previous_levels = os_stageLevels_corrected[1],
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("OS_a1", "OS_a2"),
            c("OS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 7, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          temp_factor * (alpha - os_stageLevels_corrected[1]))
      } else if (
        p_values_frailty[i, "OS_a1"] <=
          os_alphaIncr_corrected[1] * temp_factor_interim
      ) {
        decisions_frailty[i, 7, "PFS_a1"] <- (p_values_frailty[i, "PFS_a1"] <=
          alpha)
      } else {
        temp_factor <- inflation_factor(
          current_levels = os_alphaIncr_corrected[2],
          previous_levels = temp_factor_interim *
            c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("PFS_a1", "OS_a1", "OS_a2"),
            c("PFS_a1", "OS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 7, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          temp_factor * os_alphaIncr_corrected[2])
      }

      ### TESTING STRATEGY 8 ###
      # Test of PFS at interim analysis (level  0.005)
      # Group-sequential test of OS:
      #   Increase alpha to be spent at interim and final analysis by 0.005 if PFS rejected
      temp_factor_interim <- inflation_factor(
        current_levels = c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
        available_level = alpha * bretz_weight_pfs + os_alphaIncr_corrected[1],
        covariance_matrix = var_ests[
          i,
          c("PFS_a1", "OS_a1"),
          c("PFS_a1", "OS_a1")
        ]
      )
      if (
        p_values[i, "PFS_a1"] <= alpha * bretz_weight_pfs * temp_factor_interim
      ) {
        decisions[i, 8, "PFS_a1"] <- 1
        decisions[i, 8, "OS_a1"] <- (p_values[i, "OS_a1"] <=
          os_stageLevels_corrected[1] + alpha * bretz_weight_pfs)
        temp_factor <- inflation_factor(
          current_levels = alpha -
            (os_stageLevels_corrected[1] + alpha * bretz_weight_pfs),
          previous_levels = os_stageLevels_corrected[1] +
            alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("OS_a1", "OS_a2"),
            c("OS_a1", "OS_a2")
          ]
        )
        decisions[i, 8, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          temp_factor *
            (alpha - (os_stageLevels_corrected[1] + alpha * bretz_weight_pfs)))
      } else if (
        p_values[i, "OS_a1"] <= os_alphaIncr_corrected[1] * temp_factor_interim
      ) {
        decisions[i, 8, "PFS_a1"] <- (p_values[i, "PFS_a1"] <= alpha)
      } else {
        temp_factor <- inflation_factor(
          current_levels = os_alphaIncr_corrected[2],
          previous_levels = temp_factor_interim *
            c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
          available_level = alpha,
          covariance_matrix = var_ests[
            i,
            c("PFS_a1", "OS_a1", "OS_a2"),
            c("PFS_a1", "OS_a1", "OS_a2")
          ]
        )
        decisions[i, 8, "OS_a2"] <- (p_values[i, "OS_a2"] <=
          temp_factor * os_alphaIncr_corrected[2])
      }

      temp_factor_interim <- inflation_factor(
        current_levels = c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
        available_level = alpha * bretz_weight_pfs + os_alphaIncr_corrected[1],
        covariance_matrix = var_ests_frailty[
          i,
          c("PFS_a1", "OS_a1"),
          c("PFS_a1", "OS_a1")
        ]
      )
      if (
        p_values_frailty[i, "PFS_a1"] <=
          alpha * bretz_weight_pfs * temp_factor_interim
      ) {
        decisions_frailty[i, 8, "PFS_a1"] <- 1
        decisions_frailty[i, 8, "OS_a1"] <- (p_values_frailty[i, "OS_a1"] <=
          os_stageLevels_corrected[1] + alpha * bretz_weight_pfs)
        temp_factor <- inflation_factor(
          current_levels = alpha -
            (os_stageLevels_corrected[1] + alpha * bretz_weight_pfs),
          previous_levels = os_stageLevels_corrected[1] +
            alpha * bretz_weight_pfs,
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("OS_a1", "OS_a2"),
            c("OS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 8, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          temp_factor *
            (alpha - (os_stageLevels_corrected[1] + alpha * bretz_weight_pfs)))
      } else if (
        p_values_frailty[i, "OS_a1"] <=
          os_alphaIncr_corrected[1] * temp_factor_interim
      ) {
        decisions_frailty[i, 8, "PFS_a1"] <- (p_values_frailty[i, "PFS_a1"] <=
          alpha)
      } else {
        temp_factor <- inflation_factor(
          current_levels = os_alphaIncr_corrected[2],
          previous_levels = temp_factor_interim *
            c(alpha * bretz_weight_pfs, os_alphaIncr_corrected[1]),
          available_level = alpha,
          covariance_matrix = var_ests_frailty[
            i,
            c("PFS_a1", "OS_a1", "OS_a2"),
            c("PFS_a1", "OS_a1", "OS_a2")
          ]
        )
        decisions_frailty[i, 8, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
          temp_factor * os_alphaIncr_corrected[2])
      }

      ### TESTING STRATEGY 9 ###
      # Only test OS at final analysis at full level
      decisions[i, 9, "OS_a2"] <- (p_values[i, "OS_a2"] <= alpha)

      decisions_frailty[i, 9, "OS_a2"] <- (p_values_frailty[i, "OS_a2"] <=
        alpha)
    }

    metrics_names <- c("Rej_PFS", "Rej_OS", "Rej_One", "Rej_Both", "Early_Stop")
    results[,, "Rej_PFS"] <- pmax(
      decisions[,, "PFS_a1"],
      decisions[,, "PFS_a2"]
    )
    results[,, "Rej_OS"] <- pmax(decisions[,, "OS_a1"], decisions[,, "OS_a2"])
    results[,, "Rej_One"] <- pmax(results[,, "Rej_PFS"], results[,, "Rej_OS"])
    results[,, "Rej_Both"] <- pmin(results[,, "Rej_PFS"], results[,, "Rej_OS"])
    results[,, "Early_Stop"] <- pmin(
      decisions[,, "PFS_a1"],
      decisions[,, "OS_a1"]
    )

    results_frailty[,, "Rej_PFS"] <- pmax(
      decisions_frailty[,, "PFS_a1"],
      decisions_frailty[,, "PFS_a2"]
    )
    results_frailty[,, "Rej_OS"] <- pmax(
      decisions_frailty[,, "OS_a1"],
      decisions_frailty[,, "OS_a2"]
    )
    results_frailty[,, "Rej_One"] <- pmax(
      results_frailty[,, "Rej_PFS"],
      results_frailty[,, "Rej_OS"]
    )
    results_frailty[,, "Rej_Both"] <- pmin(
      results_frailty[,, "Rej_PFS"],
      results_frailty[,, "Rej_OS"]
    )
    results_frailty[,, "Early_Stop"] <- pmin(
      decisions_frailty[,, "PFS_a1"],
      decisions_frailty[,, "OS_a1"]
    )

    # Prepare additional information
    power_metrics_info <- matrix(NA, ncol = 3, nrow = num_strategies)
    colnames(power_metrics_info) <- c("scenario", "recruitment_rate", "frailty")
    power_metrics_info[, 1] <- sec
    power_metrics_info[, 2] <- rr

    # Compute power/error rate characteristics
    power_metrics_info[, 3] <- "NO"
    power_metrics <- apply(results, MARGIN = c(2, 3), FUN = mean)
    power_metrics <- cbind(power_metrics_info, power_metrics)

    power_metrics_info[, 3] <- "YES"
    power_metrics_frailty <- apply(
      results_frailty,
      MARGIN = c(2, 3),
      FUN = mean
    )
    power_metrics_frailty <- cbind(power_metrics_info, power_metrics_frailty)

    list(power_metrics, power_metrics_frailty)
  }

stopCluster(my_cluster)

save(power_metrics, file = "results/data/closed_testing_fwer.Rda")
