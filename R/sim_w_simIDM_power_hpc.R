hpc_version <- FALSE

if (hpc_version) {
  .libPaths("/home/d/danzerm/R/library")
}

require(foreach)
require(doParallel)

# Setup for parallel simulation
cores_to_use <- 32
my_cluster <- makeCluster(cores_to_use)
registerDoParallel(my_cluster)
if (hpc_version) {
  clusterEvalQ(my_cluster, .libPaths("/home/d/danzerm/R/library"))
}

secs <- 1:4
effect_scales <- c(0, seq(0.6, 1.2, 0.1))
combinations <- expand.grid(secs, effect_scales)

sim_num <- 100000

power_metrics <- foreach(c = 1:cores_to_use, .errorhandling = "pass") %dopar%
  {
    if (hpc_version) {
      .libPaths("/home/d/danzerm/R/library")
    }

    source("stat_and_var_calc.R")

    require(mvtnorm)
    require(simIDM)
    require(dplyr)
    require(rpact)

    sec <- combinations[c, 1]
    #Scaling the effect
    #(0 = no effect,
    # 1 = effect from scenarios of Erdmann et al.,
    # values between 0 and 1 = linear interpolation of hazard rates)
    effect_scale <- combinations[c, 2]

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

    test_stats <- matrix(NA, ncol = num_tests, nrow = sim_num)
    var_ests <- array(NA, dim = c(sim_num, num_tests, num_tests))
    p_values <- matrix(NA, ncol = num_tests, nrow = sim_num)
    decisions <- array(0, dim = c(sim_num, num_strategies, num_tests))
    results <- array(0, dim = c(sim_num, num_strategies, num_metrics))
    power_metrics <- matrix(NA, nrow = num_strategies, ncol = num_metrics)

    colnames(test_stats) <-
      dimnames(var_ests)[[2]] <- dimnames(var_ests)[[3]] <-
        colnames(p_values) <-
          dimnames(decisions)[[3]] <- test_names
    dimnames(decisions)[[2]] <- dimnames(results)[[2]] <-
      rownames(power_metrics) <- strategy_names
    dimnames(results)[[3]] <- colnames(power_metrics) <- metrics_names

    if (sec == 1) {
      transitionTrt <- exponential_transition(h01 = 0.06, h02 = 0.3, h12 = 0.3)
      transitionPl <- exponential_transition(h01 = 0.1, h02 = 0.4, h12 = 0.3)
      nEventPFS <- 433
      nEventOS <- 630
    }

    if (sec == 2) {
      transitionTrt <- exponential_transition(h01 = 0.3, h02 = 0.28, h12 = 0.5)
      transitionPl <- exponential_transition(h01 = 0.5, h02 = 0.3, h12 = 0.6)
      nEventPFS <- 452
      nEventOS <- 747
    }

    if (sec == 3) {
      transitionTrt <- exponential_transition(
        h01 = 0.14,
        h02 = 0.112,
        h12 = 0.25
      )
      transitionPl <- exponential_transition(
        h01 = 0.18,
        h02 = 0.15,
        h12 = 0.255
      )
      nEventPFS <- 644
      nEventOS <- 742
    }

    if (sec == 4) {
      transitionTrt <- exponential_transition(
        h01 = 0.18,
        h02 = 0.06,
        h12 = 0.17
      )
      transitionPl <- exponential_transition(h01 = 0.23, h02 = 0.07, h12 = 0.19)
      nEventPFS <- 940
      nEventOS <- 963
    }

    # Compute interpolation between control group parameters and treatment group parameters
    # according to variable effect_scale
    transitionTrt$hazards <- Map(
      "+",
      transitionPl$hazards,
      lapply(
        Map("-", transitionTrt$hazards, transitionPl$hazards),
        FUN = function(x) x * effect_scale
      )
    )
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
        nPat = c(800, 800), #seed = 181993,
        transitionByArm = transitionList,
        dropout = list(rate = drprate, time = drptime),
        accrual = list(param = "intensity", value = 25)
      ))

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

      study_result <- test_stat_calc(complete_data_temp)
      test_stats[i, ] <- study_result[[1]]
      var_ests[i, , ] <- study_result[[2]]
      p_values[i, ] <- pnorm(study_result[[1]] / sqrt(diag(study_result[[2]])))

      ##### DETERMINATION OF DECISIONS FOR DIFFERENT TESTING STRATEGIES #####

      ### TESTING STRATEGY 1 ###
      # Test of PFS at interim analysis (level 0.005)
      # Test of OS at final analysis (level 0.02)
      # No propagation
      decisions[i, 1, "PFS_a1"] <- (p_values[i, "PFS_a1"] <=
        alpha * bretz_weight_pfs)
      decisions[i, 1, "OS_a2"] <- (p_values[i, "OS_a2"] <=
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

      ### TESTING STRATEGY 6 ###
      # Test of PFS at interim analysis (level  0.005)
      # Group-sequential test of OS:
      #   Bounds as specified in Erdmann et al.
      # Propagation from PFS to OS possible
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

      ### TESTING STRATEGY 9 ###
      # Only test OS at final analysis at full level
      decisions[i, 9, "OS_a2"] <- (p_values[i, "OS_a2"] <= alpha)
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

    power_metrics <- apply(results, MARGIN = c(2, 3), FUN = mean)
    power_metrics_info <- matrix(NA, ncol = 2, nrow = num_strategies)
    colnames(power_metrics_info) <- c("scenario", "effect_scale")
    power_metrics_info[, 1] <- sec
    power_metrics_info[, 2] <- effect_scale
    power_metrics <- cbind(power_metrics_info, power_metrics)

    power_metrics
  }

stopCluster(my_cluster)

save(power_metrics, file = "results/data/closed_testing_power_metrics.Rda")
