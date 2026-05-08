require(ggplot2)
require(ggpubr)

override_r_PFS <- 33 / 64
override_r_OS <- 46 / 64

file_name_string <- "closed_testing_fwer"
if (!is.na(override_r_PFS) & !is.na(override_r_OS)) {
  file_name_string <- paste0(
    file_name_string,
    "_",
    override_r_PFS * 64,
    "_",
    override_r_OS * 64
  )
}

load(paste0("results/data/", file_name_string, ".Rda"))
power_metrics <- unlist(power_metrics, recursive = FALSE)

power_metrics <- lapply(power_metrics, function(mat) {
  mat[which(mat == "YES")] <- "1"
  mat[which(mat == "NO")] <- "0"
  mat
})

power_metrics <- lapply(power_metrics, function(mat) apply(mat, 2, as.numeric))

strategies_old <- c(
  "BON",
  "EX/LAST",
  "EX/FIRST",
  "BON/GS",
  "EX/GS/LAST",
  "EX/GS/FIRST",
  "OS"
)
strategies <- c(
  "BON",
  "REC",
  "EX/LAST",
  "EX/FIRST",
  "BON/GS",
  "REC/GS",
  "EX/GS/LAST",
  "EX/GS/FIRST",
  "OS"
)

power_metrics <- lapply(power_metrics, as.data.frame)
power_metrics <- lapply(power_metrics, FUN = function(x) {
  x$Strategy <- strategies
  rownames(x) <- NULL
  return(x)
})

save_file_string <- "fwer_restr"
if (!is.na(override_r_PFS) & !is.na(override_r_OS)) {
  save_file_string <- paste0(
    save_file_string,
    "_PFS",
    override_r_PFS * 64,
    "_OS",
    override_r_OS * 64
  )
}

### CHOOSE SUBSET OF ALL COMBINATIONS FOR PLOT IN MAIN MANUSCRIPT
plot_list <- list()
for (i in 1:4) {
  save_file_string_scenario <- paste0(
    save_file_string,
    "_scenario",
    i
  )

  strategy_subset <- c("BON", "EX/LAST", "OS")

  power_metrics_matrix <- do.call(rbind, power_metrics)
  plot_matrix <- power_metrics_matrix[
    which(
      power_metrics_matrix$scenario == i &
        power_metrics_matrix$Strategy %in% strategy_subset
    ),
  ]
  plot_matrix$n <- plot_matrix$recruitment_rate * 32 * 2
  plot_matrix$Frailty <- factor(plot_matrix$frailty, labels = c("No", "Yes"))

  # Compute bounds for confidence interval
  alpha <- 0.025
  runs <- 100000
  lb <- alpha - qnorm(0.975) * sqrt(alpha * (1 - alpha) / runs)
  ub <- alpha + qnorm(0.975) * sqrt(alpha * (1 - alpha) / runs)

  fwer_plot <-
    ggplot(
      data = plot_matrix,
      aes(x = n, y = Rej_One, linetype = Frailty, shape = Strategy)
    ) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    ylim(0.022, 0.028) +
    geom_hline(yintercept = 0.025) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = lb,
      ymax = ub,
      alpha = 0.15
    ) +
    ylab("FWER") +
    theme(text = element_text(size = 15))
  ggsave(
    paste0("results/plots/", save_file_string_scenario, ".pdf"),
    fwer_plot,
    device = "pdf",
    width = 9,
    height = 6
  )
  plot_list[[i]] <- fwer_plot
}

overall_fwer_plot <- do.call(
  ggarrange,
  c(plot_list, common.legend = TRUE, legend = "right")
)
ggsave(
  paste0("results/plots/", save_file_string, ".pdf"),
  overall_fwer_plot,
  device = "pdf",
  width = 12,
  height = 8
)

save_file_string <- "fwer_remaining"
if (!is.na(override_r_PFS) & !is.na(override_r_OS)) {
  save_file_string <- paste0(
    save_file_string,
    "_PFS",
    override_r_PFS * 64,
    "_OS",
    override_r_OS * 64
  )
}

### SHOW REMAINING METHODS IN PLOT IN SUPPLEMENT
plot_list <- list()
for (i in 1:4) {
  save_file_string_scenario <- paste0(
    save_file_string,
    "_scenario",
    i
  )

  strategy_subset <- c(
    "REC",
    "EX/FIRST",
    "BON/GS",
    "REC/GS",
    "EX/GS/LAST",
    "EX/GS/FIRST"
  )

  power_metrics_matrix <- do.call(rbind, power_metrics)
  plot_matrix <- power_metrics_matrix[
    which(
      power_metrics_matrix$scenario == i &
        power_metrics_matrix$Strategy %in% strategy_subset
    ),
  ]
  plot_matrix$n <- plot_matrix$recruitment_rate * 32 * 2
  plot_matrix$Frailty <- factor(plot_matrix$frailty, labels = c("No", "Yes"))

  # Compute bounds for confidence interval
  alpha <- 0.025
  runs <- 100000
  lb <- alpha - qnorm(0.975) * sqrt(alpha * (1 - alpha) / runs)
  ub <- alpha + qnorm(0.975) * sqrt(alpha * (1 - alpha) / runs)

  fwer_plot <-
    ggplot(
      data = plot_matrix,
      aes(x = n, y = Rej_One, linetype = Frailty, shape = Strategy)
    ) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    ylim(0.0225, 0.0275) +
    geom_hline(yintercept = 0.025) +
    annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = lb,
      ymax = ub,
      alpha = 0.15
    ) +
    ylab("FWER") +
    theme(text = element_text(size = 15))
  ggsave(
    paste0("results/plots/", save_file_string_scenario, ".pdf"),
    fwer_plot,
    device = "pdf",
    width = 9,
    height = 6
  )
  plot_list[[i]] <- fwer_plot
}

overall_fwer_plot <- do.call(
  ggarrange,
  c(plot_list, common.legend = TRUE, legend = "right")
)
ggsave(
  paste0("results/plots/", save_file_string, ".pdf"),
  overall_fwer_plot,
  device = "pdf",
  width = 12,
  height = 8
)
