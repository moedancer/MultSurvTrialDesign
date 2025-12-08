# MultSurvTrialDesign
Exhausting the type I error level in event-driven group-sequential designs with a closed testing procedure for progression-free and overall survival

This repository contains the following R scripts:
- `stat_and var_calc.R`: Compute test statistics for survival data for the endpoints progression-free survival (PFS) and overall survival (OS) in a design with two analyses, as well as the covariance matrix for those four test statistics. Compute inflation factors for local significance levels to exhaust the type I error level in a closed testing procedure.
- `custom_censoring.R`: Generate custom censoring for uncensored data on PFS and OS
- `sim_w_sim_IDM_H0_hpc.R`: Simulate event-driven group-sequential designs with a closed testing procedure for PFS and OS in several scenarios. Estimate the family-wise type I error level under the global null hypothesis for multiple designs.
- `sim_w_sim_IDM_power_hpc.R`: Simulate event-driven group-sequential designs with a closed testing procedure for PFS and OS in several scenarios. Estimate several performance measures (e.g. power for PFS, power for OS, disjunctive power, conjunctive power) for multiple designs with different effects of the experimental therapy.
- `create_figures_fwer.R` and `create_figures.R`. Summarize the estimated type I error levels and performance measures in plots.

This repository contains the following data:
- In the folder `results/data`: Estimated type I error levels under the global null hypothesis (results of `sim_w_sim_IDM_H0_hpc.R`) and estimated performance measures under alternatives (results of `sim_w_sim_IDM_power_hpc.R`).
- In the folder `results/plots`: Figures created with the scripts `create_figures_fwer.R` and `create_figures.R`.