require(rio)
require(ggplot2)
require(grid)
require(gridExtra)
require(ggpubr)

load("results/data/closed_testing_power_metrics.Rda")

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

  # Absolute gain in disj. power w.r.t. corresponding Bonferroni-corrected procedure
  x$power_gain_disj <- ifelse(
    grepl("GS", strategies),
    x$Rej_One - x$Rej_One[x$Strategy == "BON/GS"],
    x$Rej_One - x$Rej_One[x$Strategy == "BON"]
  )
  # Relative gain in disj. power w.r.t. corresponding Bonferroni-corrected procedure
  x$rel_power_gain_disj <- ifelse(
    grepl("GS", strategies),
    1 + x$power_gain_disj / x$Rej_One[x$Strategy == "BON/GS"],
    1 + x$power_gain_disj / x$Rej_One[x$Strategy == "BON"]
  )

  # Absolute power loss w.r.t. pure OS testing
  x$power_lost_os <- x$Rej_OS - x$Rej_OS[x$Strategy == "OS"]
  # Absolute power gain w.r.t. corresponding Bonferroni-corrected procedure
  x$power_gain_os <- ifelse(
    grepl("GS", strategies),
    x$Rej_OS - x$Rej_OS[x$Strategy == "BON/GS"],
    x$Rej_OS - x$Rej_OS[x$Strategy == "BON"]
  )
  # Relative power loss w.r.t. pure OS testing
  x$rel_power_os_loss <- 1 + x$power_lost_os / x$Rej_OS[x$Strategy == "OS"]
  # Relative power gain w.r.t. corresponding Bonferroni-corrected procedure
  x$rel_power_os_gain <- ifelse(
    grepl("GS", strategies),
    1 + x$power_gain_os / x$Rej_OS[x$Strategy == "BON/GS"],
    1 + x$power_gain_os / x$Rej_OS[x$Strategy == "BON"]
  )

  return(x)
})

power_metrics_matrix <- do.call(rbind, power_metrics)
power_metrics_list <- lapply(X = 1:4, FUN = function(x) {
  power_metrics_matrix[power_metrics_matrix[, "scenario"] == x, ]
})

power_metrics_list <- lapply(power_metrics_list, as.data.frame)

for (i in 1:4) {
  scen_df <- power_metrics_list[[i]]
  scen_df <- scen_df[scen_df$effect_scale != 0, ]
  scen_df$Strategy <- factor(scen_df$Strategy)

  plot_df <- scen_df[
    scen_df$Strategy %in%
      c("REC", "EX/LAST", "EX/GS/LAST", "EX/FIRST", "OS"),
  ]
  plot_df$rel_power_os_gain[plot_df$Strategy == "EX/GS/LAST"] <- NA
  rel_power_gain_os_plot <-
    ggplot(
      data = plot_df,
      aes(x = effect_scale, y = rel_power_os_gain, group = Strategy)
    ) +
    geom_hline(yintercept = 1, colour = "red", linewidth = 1.5) +
    geom_line(aes(linetype = Strategy), size = 1) +
    scale_linetype_manual(
      values = c(
        "REC" = "twodash",
        "EX/LAST" = "solid",
        "EX/FIRST" = "dashed",
        "EX/GS/LAST" = "dotdash",
        "OS" = "dotted"
      )
    ) +
    geom_point(size = 3) +
    scale_y_continuous(breaks = seq(1, 1.125, 0.025), limits = c(1, 1.125)) +
    xlab("w") +
    ylab("Rel. power compared to Bonferroni correction") +
    theme(text = element_text(size = 15))

  plot_df <- scen_df[
    scen_df$Strategy %in%
      c("REC", "EX/LAST", "EX/GS/LAST", "EX/FIRST", "OS"),
  ]
  plot_df$rel_power_gain_disj[plot_df$Strategy %in% c("REC", "EX/FIRST")] <- NA
  rel_power_gain_disj_plot <-
    ggplot(
      data = plot_df,
      aes(x = effect_scale, y = rel_power_gain_disj, group = Strategy)
    ) +
    geom_hline(yintercept = 1, colour = "red", linewidth = 1.5) +
    geom_line(aes(linetype = Strategy), size = 1) +
    scale_linetype_manual(
      values = c(
        "REC" = "twodash",
        "EX/LAST" = "solid",
        "EX/FIRST" = "dashed",
        "EX/GS/LAST" = "dotdash",
        "OS" = "dotted"
      )
    ) +
    geom_point(size = 3) +
    scale_y_continuous(breaks = seq(0.85, 1.05, 0.05), limits = c(0.85, 1.05)) +
    xlab("w") +
    ylab("Rel. power compared to Bonferroni correction") +
    theme(text = element_text(size = 15))

  rel_power_coll <- ggarrange(
    rel_power_gain_os_plot,
    rel_power_gain_disj_plot,
    common.legend = TRUE,
    legend = "right"
  )
  ggsave(
    paste(
      "results/plots/rel_power_",
      i,
      ".pdf",
      sep = ""
    ),
    rel_power_coll,
    device = "pdf",
    width = 12,
    height = 6
  )
}

scen_df <- power_metrics_list[[1]]
scen_df <- scen_df[scen_df$effect_scale != 0, ]
scen_df$Strategy <- factor(scen_df$Strategy)

plot_df <- scen_df[
  scen_df$Strategy %in%
    c("EX/LAST", "EX/GS/LAST", "EX/FIRST", "OS"),
]
plot_df$rel_power_os_gain[plot_df$Strategy == "EX/GS/LAST"] <- NA
rel_power_gain_os_plot <-
  ggplot(
    data = plot_df,
    aes(x = effect_scale, y = rel_power_os_gain, group = Strategy)
  ) +
  geom_line(aes(linetype = Strategy), size = 1) +
  scale_linetype_manual(
    values = c(
      "EX/LAST" = "solid",
      "EX/FIRST" = "dashed",
      "EX/GS/LAST" = "longdash",
      "OS" = "dotted"
    )
  ) +
  geom_point(size = 3) +
  ylim(1, NA) +
  xlab("w") +
  ylab("Rel. power compared to Bonferroni correction") +
  theme(text = element_text(size = 15))

plot_df <- scen_df[
  scen_df$Strategy %in%
    c("EX/LAST", "EX/GS/LAST", "EX/FIRST", "OS"),
]
plot_df$rel_power_gain_disj[plot_df$Strategy == "EX/FIRST"] <- NA
rel_power_gain_disj_plot <-
  ggplot(
    data = plot_df,
    aes(x = effect_scale, y = rel_power_gain_disj, group = Strategy)
  ) +
  geom_line(aes(linetype = Strategy), size = 1) +
  scale_linetype_manual(
    values = c(
      "EX/LAST" = "solid",
      "EX/FIRST" = "dashed",
      "EX/GS/LAST" = "longdash",
      "OS" = "dotted"
    )
  ) +
  geom_point(size = 3) +
  xlab("w") +
  ylab("Rel. power compared to Bonferroni correction") +
  theme(text = element_text(size = 15))

ggarrange(
  rel_power_gain_os_plot,
  rel_power_gain_disj_plot,
  common.legend = TRUE,
  legend = "right"
)
ggsave(
  "results/plots/rel_power.pdf",
  overall_fwer_plot,
  device = "pdf",
  width = 12,
  height = 8
)
