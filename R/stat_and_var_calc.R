require(dplyr)
require(mvtnorm)
require(matrixcalc)

### Computation of test statistics and covariance matrix
# complete_study_data has to be a data frame with columns
# group: treatment group indicator (0/1)
# PFStime_a1, PFSevent_a1: censored PFS time and event indicator for analysis 1
# OStime_a1, OSevent_a1: censored OS time and event indicator for analysis 1
# PFStime_a2, PFSevent_a2: censored PFS time and event indicator for analysis 2
# OStime_a2, OSevent_a2: censored OS time and event indicator for analysis 2
test_stat_calc <- function(complete_study_data) {
  n <- dim(complete_study_data)[1]

  test_stats <- rep(NA, 4)
  cov_mat <- matrix(NA, ncol = 4, nrow = 4)

  name_vec <- c("PFS_a1", "OS_a1", "PFS_a2", "OS_a2")
  names(test_stats) <- rownames(cov_mat) <- colnames(cov_mat) <- name_vec

  pfs_a1_obj <- w_calculator(complete_study_data[, c(
    "group",
    "PFStime_a1",
    "PFSevent_a1"
  )])
  test_stats["PFS_a1"] <- 1 / sqrt(n) * pfs_a1_obj[[1]]
  w_pfs_a1 <- pfs_a1_obj[[3]]
  cov_mat["PFS_a1", "PFS_a1"] <- cov_mat["PFS_a1", "PFS_a2"] <- cov_mat[
    "PFS_a2",
    "PFS_a1"
  ] <- pfs_a1_obj[[2]]

  os_a1_obj <- w_calculator(complete_study_data[, c(
    "group",
    "OStime_a1",
    "OSevent_a1"
  )])
  test_stats["OS_a1"] <- 1 / sqrt(n) * os_a1_obj[[1]]
  w_os_a1 <- os_a1_obj[[3]]
  cov_mat["OS_a1", "OS_a1"] <- cov_mat["OS_a1", "OS_a2"] <- cov_mat[
    "OS_a2",
    "OS_a1"
  ] <- os_a1_obj[[2]]

  pfs_a2_obj <- w_calculator(complete_study_data[, c(
    "group",
    "PFStime_a2",
    "PFSevent_a2"
  )])
  test_stats["PFS_a2"] <- 1 / sqrt(n) * pfs_a2_obj[[1]]
  w_pfs_a2 <- pfs_a2_obj[[3]]
  cov_mat["PFS_a2", "PFS_a2"] <- pfs_a2_obj[[2]]

  os_a2_obj <- w_calculator(complete_study_data[, c(
    "group",
    "OStime_a2",
    "OSevent_a2"
  )])
  test_stats["OS_a2"] <- 1 / sqrt(n) * os_a2_obj[[1]]
  w_os_a2 <- os_a2_obj[[3]]
  cov_mat["OS_a2", "OS_a2"] <- os_a2_obj[[2]]

  cov_mat["PFS_a1", "OS_a1"] <- cov_mat["OS_a1", "PFS_a1"] <- (w_pfs_a1 %*%
    w_os_a1) /
    n
  cov_mat["PFS_a1", "OS_a2"] <- cov_mat["OS_a2", "PFS_a1"] <- (w_pfs_a1 %*%
    w_os_a2) /
    n
  cov_mat["PFS_a2", "OS_a1"] <- cov_mat["OS_a1", "PFS_a2"] <- (w_pfs_a2 %*%
    w_os_a1) /
    n
  cov_mat["PFS_a2", "OS_a2"] <- cov_mat["OS_a2", "PFS_a2"] <- (w_pfs_a2 %*%
    w_os_a2) /
    n

  return(list(test_stats, cov_mat))
}

# univ_data has to be a data frame with columns group, time and status
w_calculator <- function(univ_data) {
  univ_data <-
    univ_data %>%
    ungroup() %>%
    rename(group = 1, time = 2, status = 3) %>%
    mutate(totalrank = rank(desc(time), ties.method = "max")) %>%
    group_by(group) %>%
    mutate(grouprank = rank(desc(time), ties.method = "max"))

  new_group_prop <- ifelse(
    univ_data$group == 1,
    univ_data$grouprank / univ_data$totalrank,
    1 - univ_data$grouprank / univ_data$totalrank
  )
  w_count <- univ_data$status * (univ_data$group - new_group_prop)

  haz_incr <- univ_data$status / univ_data$totalrank
  weighted_haz_incr <- new_group_prop * haz_incr
  lr_var_incr <- univ_data$status * new_group_prop * (1 - new_group_prop)

  haz <- cumsum(haz_incr[order(univ_data$time)])[rank(
    univ_data$time,
    ties.method = "first"
  )]
  weighted_haz <- cumsum(weighted_haz_incr[order(univ_data$time)])[rank(
    univ_data$time,
    ties.method = "first"
  )]

  w_comp <- univ_data$group * haz - weighted_haz

  lr_stat <- sum(w_count)
  lr_var <- mean(lr_var_incr)

  return(list(lr_stat, lr_var, w_count - w_comp))
}

### Computation of inflation factors for current analysis
# First argument: Local levels that were actually used in previous analyses
#                 for intersection hypothesis to be investigated
# Second argument: Local levels for current analysis according
#                  to pre-specified alpha-spending functions
# Third argument: Level that shall be exhausted for
inflation_factor <- function(
  current_levels,
  available_level,
  covariance_matrix,
  previous_levels = NULL
) {
  if (!is.positive.semi.definite(covariance_matrix)) {
    return(1)
  } else {
    total_length <- length(current_levels) + length(previous_levels)
    inflation_helper <- function(factor) {
      1 -
        pmvnorm(
          lower = rep(-Inf, total_length),
          upper = c(
            qnorm(1 - previous_levels),
            qnorm(1 - factor * current_levels)
          ) *
            sqrt(diag(covariance_matrix)),
          sigma = covariance_matrix
        ) -
        available_level
    }
    return(uniroot(inflation_helper, lower = 0.5, upper = 20)$root)
  }
}
