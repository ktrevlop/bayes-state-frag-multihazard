# K. Trevlopoulos
# March 2025

# Bayesian modelling and model comparison for fragility models including
# state-dependent fragility models
#
# Note that the simulation steps are numbered 0 to 3 in the python scripts,
# and in the .csv files, but in the figures here they are numbered 1 to 4.
#
# This script calls the following scripts:
# ./model_1_model_DeltaT_HE_vs_DR.R


library(rethinking)
library(tidyverse)
library(ordinal)
library(ggplot2)
library(here)
library(latex2exp)
library(ggpubr)
library(posterior)
library(bayesplot)
library(loo)
here::here()
set_cmdstan_path("/gpfs/scratch/trevlopoulos/MyLibsR2/cmdstan/cmdstan-2.36.0")
# For loo:
options(mc.cores = 4)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Input ----

print("Processing input...")

# The data for the fragility curves
# comp_thist_results <- read_csv(
#   './data_store/composite_time_history_assembly_EDPs_IMs_CR_LFINF_CDM_11_HEX_2_20240903.csv',
#     show_col_types = FALSE)
orig_dataset <- read_csv(
  './data_store/dataset_CR_LFINF_CDM_11_HEX_2_20240903.csv',
  show_col_types = FALSE)

# Enter the names of the damage states
ds_names <- c('DS0', 'DS1', 'DS2', 'DS3', 'DS4')

# The thresholds of the damage states (maximum inter-story drift)
thresholds <- c(0.0021, 0.0035, 0.0049, 0.0063)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Preparation of the data for the calculations ----

# comp_thist_results will eventually have the data for part_2
comp_thist_results <- select(
  orig_dataset,
  simulation, part_number, load_type, max_int_storey_drift, PGA_g,
  cum_hyst_loop_area, per_of_max_amp)

# Create a column for the damage states based on the thresholds
# ATTENTION: This tags the damage state in a part of the simulation without
# considering the damage state in the previous parts of the simulation.
# This is taken into account by the code that follows.
comp_thist_results <- mutate(
  comp_thist_results,
  DS = case_when(
    max_int_storey_drift < thresholds[1] ~ 0,
    max_int_storey_drift >= thresholds[1] & max_int_storey_drift < thresholds[2] ~ 1,
    max_int_storey_drift >= thresholds[2] & max_int_storey_drift < thresholds[3] ~ 2,
    max_int_storey_drift >= thresholds[3] & max_int_storey_drift < thresholds[4] ~ 3,
    max_int_storey_drift >= thresholds[4] ~ 4,
  )
)


# The normalized period
# The fist egenperiod of the undamanged building (s)
T_0 <- min(comp_thist_results$per_of_max_amp, na.rm = TRUE)

# Shift values one row up to use the period at step 1 (numbering 0-3)
comp_thist_results <- comp_thist_results %>%
  mutate(per_of_max_amp = lead(per_of_max_amp))

comp_thist_results <- mutate( comp_thist_results,
                              lnTnorm = log(per_of_max_amp/T_0),
                              lnDeltaT = log((per_of_max_amp-T_0)/T_0),
)

# The numbers of the simulations with load_type == "strong_gm" at steps 1 and 3
# (for steps numbered 1-4). Reminder: in the python scripts and in the .csv
# files the parts are numbered 0-3.
sim_step_0 <- filter(comp_thist_results, part_number == 0)
# sim_step_0 <- filter(sim_step_0, load_type == "strong_gm")
sim_step_0 <- select(sim_step_0, simulation)

sim_step_2 <- filter(comp_thist_results, part_number == 2)
sim_step_2 <- filter(sim_step_2, load_type == "strong_gm")
sim_step_2 <- select(sim_step_2, simulation)

# The indices of the simulations for the selected hazard types in the
# simulation steps 0 and 2
n_of_sim <- inner_join(sim_step_2, sim_step_0)



# Pivot the DS to separate columns per part of the simulation
# ds_per_part_wide will eventually have the data for part_0
ds_per_part_wide <- pivot_wider(
  comp_thist_results,
  names_from = part_number, values_from = DS, id_cols = simulation)
# Parts numbered 1-4
ds_per_part_wide <- rename(
  ds_per_part_wide,
  DS_part_0 = `0`, DS_part_1 = `1`, DS_part_2 = `2`, DS_part_3 = `3`)

# Keep the lines for the simulations for the selected hazard types in the
# simulation steps 0 and 2
ds_per_part_wide <- left_join(n_of_sim, ds_per_part_wide)
ds_per_part_wide <- select(ds_per_part_wide, simulation, DS_part_0, DS_part_2)

# Make sure that DS_part_2 > DS_part_0
ds_per_part_wide <- rename(ds_per_part_wide, DS_part_2_temp = DS_part_2)
ds_per_part_wide <- ds_per_part_wide %>% rowwise() %>%
  mutate(DS_part_2 = max(DS_part_0, DS_part_2_temp))
ds_per_part_wide <- select(ds_per_part_wide, -DS_part_2_temp)

# The data for part_0
comp_thist_results_part_0 <- filter(comp_thist_results, part_number==0 )
comp_thist_results_part_0 <- left_join(n_of_sim, comp_thist_results_part_0)
comp_thist_results_part_0 <- rename(
  comp_thist_results_part_0,
  DS_part_0 = DS,
  max_int_storey_drift_part_0 = max_int_storey_drift,
  IM_part_0 = PGA_g,
  lnTnorm_part_0 = lnTnorm,
  lnDeltaT_part_1 = lnDeltaT,
  ht_part_0 = load_type
)
comp_thist_results_part_0 <- mutate(
  comp_thist_results_part_0,
  lnIM_part_0 = log(IM_part_0))
comp_thist_results_part_0 <- select(
  comp_thist_results_part_0,
  simulation, IM_part_0, lnIM_part_0, max_int_storey_drift_part_0,
  cum_hyst_loop_area, lnTnorm_part_0, lnDeltaT_part_1, ht_part_0)

# The data for part_2
comp_thist_results_part_2 <- filter(comp_thist_results, part_number==2 )
comp_thist_results_part_2 <- left_join(n_of_sim, comp_thist_results_part_2)
comp_thist_results_part_2 <- rename(
  comp_thist_results_part_2,
  DS_part_2 = DS,
  max_int_storey_drift_part_2 = max_int_storey_drift,
  IM_part_2 = PGA_g )
comp_thist_results_part_2 <- mutate(comp_thist_results_part_2,
                                    lnIM_part_2 = log(IM_part_2))
comp_thist_results_part_2 <- select(comp_thist_results_part_2,
                                    simulation, IM_part_2, lnIM_part_2)


comp_thist_results <- left_join(comp_thist_results_part_0,
                                comp_thist_results_part_2)



# The data for the model_2.x models
model3.x_data_tib <- left_join(ds_per_part_wide, comp_thist_results)
model3.x_data_tib <- select(model3.x_data_tib, -simulation)
model3.x_data_tib$DS_part_0 <- ordered(model3.x_data_tib$DS_part_0)
model3.x_data_tib$DS_part_2 <- ordered(model3.x_data_tib$DS_part_2)

# One-hot encoding for the hazard type
model3.x_data_tib <- mutate( model3.x_data_tib,
                             ht_part_0 = case_when(
                               ht_part_0 == "strong_gm" ~ 0,
                               ht_part_0 == "flood" ~ 1 ) )

model3.x_data_tib <- filter(model3.x_data_tib, lnIM_part_2 != -Inf)

model3.x_data <- list(
  N = nrow(model3.x_data_tib),
  lnIM_part_0 = model3.x_data_tib$lnIM_part_0,
  lnIM_part_2 = model3.x_data_tib$lnIM_part_2,
  DS_part_0 = as.integer( model3.x_data_tib$DS_part_0 ),
  DS_part_2 = as.integer( model3.x_data_tib$DS_part_2 ),
  lndrift_part_0 = log( model3.x_data_tib$max_int_storey_drift_part_0 ),
  lnhen_part_0 = log( model3.x_data_tib$cum_hyst_loop_area ),
  lnTnorm_part_0 = model3.x_data_tib$lnTnorm_part_0,
  lnDeltaT_part_1 = model3.x_data_tib$lnDeltaT_part_1,
  ht_part_0 = as.integer( model3.x_data_tib$ht_part_0 ),
  alpha = rep( 2 , 4 ) ) # delta prior for four damage states

# ATTENTION: DS_part_0 and DS_part_2 now take values 1-5 instead of 0-4.
# This helps with indexing such as 1:DS_part_0, where we would have 0:1.




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Functions ----

# The function waic_calc is not recommended for publication.
# Use the loo package instead.

# Log-mean-exp function (used by waic_calc)
log_mean_exp <- function(x) {
  max_x <- max(x)
  max_x + log(mean(exp(x - max_x)))
}

# This function calculates the WAIC. It takes as input a stan model
waic_calc <- function(stanfit_m) {
  # Extract log_lik matrix (samples × observations)
  # base R matrix
  log_lik_arr <- stanfit_m$draws("log_lik", format = "matrix")
  
  # Compute the log pointwise predictive density (lppd)
  lppd_pw <- apply(log_lik_arr, 2, log_mean_exp)
  
  # Compute effective number of parameters (pointwise var of log_lik)
  n_par_waic <- apply(log_lik_arr, 2, var)
  
  # Final WAIC computation
  waic_out <- -2 * (sum(lppd_pw) - sum(n_par_waic))
  
  model_summary <- list(
    WAIC = waic_out,
    lppd = sum(lppd_pw),
    p_waic = sum(n_par_waic),
    pointwise = data.frame(lppd = lppd, p_waic = p_waic)
  )
  
  return(model_summary)
}




# This function calculates model weights for a series of WAIC values
# Example: round(weight_calc(c(4000, 4001, 4010, 4100)),3)

weight_calc <- function(model_WAICs) {
  model_weight <- tibble(WAIC = model_WAICs)
  model_weight['dWAIC'] <- model_weight$WAIC - min(model_weight$WAIC)
  model_weight['weight'] <- exp(-0.5*model_weight$dWAIC)
  model_weight['weight'] <- model_weight$weight / sum(model_weight$weight)
  return(model_weight$weight)
}




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Auxiliary calculations ----

# source("./model_1_model_DeltaT_HE_vs_DR.R")
load("./DT_thres.RData")
load("./hen_thresholds.RData")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 3.1 ----

print("model_3.1gq")

model_3.1gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear model for eta
    eta <- b0[DS_part_0] * lnIM_part_2 + 
      b1 * sum(append_row(0, delta)[1:DS_part_0]),
    
    # Priors
    vector[max(DS_part_0)]: b0 ~ normal(4, 2),
    b1 ~ normal(0, 2.5),
    cutpoints ~ normal(-4, 2),
    simplex[4]: delta ~ dirichlet(alpha),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] | 
        b0[DS_part_0[i]] * lnIM_part_2[i] +
        b1 * sum(append_row(0, delta)[1:DS_part_0[i]]), 
      cutpoints
    )
  ),
  data = model3.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract and write Stan code to file. Comment if already exported and edited
# writeLines( stancode(model_3.1gq), "model_3_1gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_3_1gq_traceplot.png",
  width = 21.0,
  height = 29.7*0.5,
  units = "cm",
  res = 300
)
traceplot( model_3.1gq )
dev.off()

png(
  filename = "./data_store/checks/model_3_1gq_trankplot.png",
  width = 21.0,
  height = 29.7*0.5,
  units = "cm",
  res = 300
)
trankplot( model_3.1gq )
dev.off()


# Information criteria
model_3.1gq_waic <- WAIC( model_3.1gq )
model_3.1gq_psis <- PSIS( model_3.1gq )

# Summarize results
model_3.1gq_precis <- precis(model_3.1gq, depth = 2)

# Export the parameters
write_csv(
  model_3.1gq_precis |>
    as.matrix() |>
    as_tibble(rownames = "Parameter"),
  "./data_store/fragility_models/model_3.1gq_precis.csv")




## Compute fragility curves ----

# Model parameters
b0 <- model_3.1gq_precis[1:5,'mean']
b1 <- model_3.1gq_precis['b1','mean']
cutpoints <- model_3.1gq_precis[7:10,'mean']
delta <- model_3.1gq_precis[11:14,'mean']


# # Define IM range for plotting
# # Define IM range for plotting
# IM_seq <- seq(from=0.02, to=1.0, by=0.02)
# IM_seq <- c(1e-9, IM_seq)
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Number of damage states
n_ds <- length(grep("^cutpoints\\[", rownames(model_3.1gq_precis) )) + 1


# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               b0[j] * lnIM_seq - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               b0[j] * lnIM_seq - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    temp <- bind_rows(temp)
    temp <- mutate(temp, DS_part_0=j)
    
    # Remove fragility curves for DS<DSj
    for (k in 1:j) {
      temp <- filter(temp, DS != k)
    }
    
    fragility_curves <- bind_rows(fragility_curves,temp)
    
  }
  
}

fragility_model_3.1 <- fragility_curves |> mutate(DS=factor(DS)) |>
  mutate(Hazard = "All", Model = "Model 3.1")




figure_m3.1 <- ggplot(
  fragility_model_3.1, aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(
    ~ DS_part_0, ncol = 2,
    labeller = labeller(DS_part_0 = function(x) paste("DS part 1:", x))) +
  scale_color_viridis_d(name = "Damage State") +
  labs(
    x = "PGA (g)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_m3.1)

ggsave(
  plot = figure_m3.1,
  paste0("./data_store/fragility_models/model_3_1.png"),
  width = 6.0*2.54,
  height = 6.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 3.2 ----

print("model_3.2gq")

model_3.2gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear model for eta
    eta <- b0[DS_part_0] * lnIM_part_2 + 
      b1 * sum(append_row(0, delta)[1:DS_part_0]) +
      b2 * ht_part_0,
    
    # Priors
    vector[max(DS_part_0)]: b0 ~ normal(4, 2),
    b1 ~ normal(0, 2.5),
    b2 ~ normal(0, 2.5),
    cutpoints ~ normal(-4, 2),
    simplex[4]: delta ~ dirichlet(alpha),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] | 
        b0[DS_part_0[i]] * lnIM_part_2[i] +
        b1 * sum(append_row(0, delta)[1:DS_part_0[i]]) +
        b2 * ht_part_0[i], 
      cutpoints
    )
  ),
  data = model3.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# Extract and write Stan code to file. Comment if already exported and edited
# writeLines( stancode(model_3.2gq), "model_3_2gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_3_2gq_traceplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
traceplot( model_3.2gq )
dev.off()

png(
  filename = "./data_store/checks/model_3_2gq_trankplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
trankplot( model_3.2gq )
dev.off()


# Information criteria
model_3.2gq_waic <- WAIC( model_3.2gq )
model_3.2gq_psis <- PSIS( model_3.2gq )

# Summarize results
model_3.2gq_precis <- precis(model_3.2gq, depth = 2)

# Export the parameters
write_csv(
  model_3.2gq_precis |>
    as.matrix() |>
    as_tibble(rownames = "Parameter"),
  "./data_store/fragility_models/model_3.2gq_precis.csv")




## Compute fragility curves ----

load_type = c("Strong GM", "Flood")

# Model parameters
b0 <- model_3.2gq_precis[1:5,'mean']
b1 <- model_3.2gq_precis['b1','mean']
b2 <- model_3.2gq_precis['b2','mean']
cutpoints <- model_3.2gq_precis[8:11,'mean']
delta <- model_3.2gq_precis[12:15,'mean']


# # Define IM range for plotting
# IM_seq <- seq(from=0.02, to=1.0, by=0.02)
# IM_seq <- c(1e-9, IM_seq)
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Number of damage states
n_ds <- length(grep("^cutpoints\\[", rownames(model_3.2gq_precis) )) + 1


# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  for (ht in 0:1) {
    
    if (j == 1 && ht == 0) {
      
      # Compute fragility curves
      fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
        tibble(DS = k+1,
               Probability = plogis(
                 b0[j] * lnIM_seq - cutpoints[k] +
                   b1 * sum(c(0, delta)[1:j]) + b2 * ht), IM = IM_seq)
      })
      
      fragility_curves <- bind_rows(fragility_list)
      fragility_curves <- mutate(fragility_curves, DS_part_0=j)
      fragility_curves["Hazard"] = load_type[ht+1]
      
    } else {
      
      # Compute fragility curves
      temp <- purrr::map(seq_along(cutpoints), function(k) {
        tibble(DS = k+1,
               Probability = plogis(
                 b0[j] * lnIM_seq - cutpoints[k] +
                   b1 * sum(c(0, delta)[1:j]) + b2 * ht), IM = IM_seq)
      })
      
      temp <- bind_rows(temp)
      temp <- mutate(temp, DS_part_0=j)
      temp["Hazard"] = load_type[ht+1]
      
      # Remove fragility curves for DS<DSj
      if (j>1) {
        for (k in 1:j) {
          temp <- filter(temp, DS != k)
        }
      }
      
      
      fragility_curves <- bind_rows(fragility_curves,temp)
      
    }
    
  }
  
}

fragility_model_3.2 <- fragility_curves |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 3.2")




figure_m3.2<- ggplot(
  fragility_model_3.2, aes(x = IM, y = Probability,
                           color = DS, linetype = Hazard)) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ DS_part_0, ncol = 2,
    labeller = labeller(DS_part_0 = function(x) paste("DS part 1:", x))) +
  scale_color_viridis_d(name = "Damage State") +
  scale_linetype_discrete(name = "Hazard in part 1") +
  labs(
    x = "PGA (g)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_m3.2)

ggsave(
  plot = figure_m3.2,
  paste0("./data_store/fragility_models/model_3_2.png"),
  width = 21.0,
  height = 29.7*0.67,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 3.3 ----

print("model_3.3gq")

model_3.3gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear model for eta
    eta <- b0[DS_part_0] * lnIM_part_2 + 
      b1 * sum(append_row(0, delta)[1:DS_part_0]) +
      b2 * ht_part_0 + b3[DS_part_0] * exp(lnDeltaT_part_1),
    
    # Priors
    vector[max(DS_part_0)]: b0 ~ normal(4, 2),
    b1 ~ normal(0, 2.5),
    b2 ~ normal(0, 2.5),
    vector[max(DS_part_0)]: b3 ~ normal(0, 4),
    cutpoints ~ normal(-4, 2),
    simplex[4]: delta ~ dirichlet(alpha),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] | 
        b0[DS_part_0[i]] * lnIM_part_2[i] +
        b1 * sum(append_row(0, delta)[1:DS_part_0[i]]) +
        b2 * ht_part_0[i] +
        b3[DS_part_0[i]] * exp(lnDeltaT_part_1[i]), 
      cutpoints
    )
  ),
  data = model3.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract and write Stan code to file. Comment if already exported and edited
# writeLines( stancode(model_3.3gq), "model_3_3gq.stan")




# Diagnostics ----

# model_3.3a <- cmdstan_model("model_3_3gq.stan")
# 
# model_3.3a_fit <- model_3.3a$sample(
# data = model3.x_data,
# seed = 1234,
# chains = 4,
# parallel_chains = 4,
# iter_sampling = 1000,
# iter_warmup = 1000
# )
# 
# loo( model_3.3a_fit$draws("log_lik", format = "matrix") )

# Computed from 4000 by 4993 log-likelihood matrix.
# 
#          Estimate    SE
# elpd_loo  -3697.8  57.4
# p_loo        11.3   0.3
# looic      7395.6 114.8
# -
# MCSE of elpd_loo is NA.
# MCSE and ESS estimates assume independent draws (r_eff=1).
# 
# Pareto k diagnostic values:
#   Count Pct.    Min. ESS
# (-Inf, 0.7]   (good)     4990  99.9%   3482    
#    (0.7, 1]   (bad)         3   0.1%   <NA>    
#    (1, Inf)   (very bad)    0   0.0%   <NA>    
# See help('pareto-k-diagnostic') for details.




## Checking chains ----

png(
  filename = "./data_store/checks/model_3_3gq_traceplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
traceplot( model_3.3gq, n_cols = 4, max_rows = 10 )
dev.off()

png(
  filename = "./data_store/checks/model_3_3gq_trankplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
trankplot( model_3.3gq, n_cols = 4, max_rows = 10 )
dev.off()


# Information criteria
model_3.3gq_waic <- WAIC( model_3.3gq )
model_3.3gq_psis <- PSIS( model_3.3gq )

# Summarize results
model_3.3gq_precis <- precis(model_3.3gq, depth = 2)

# Export the parameters
write_csv(
  model_3.3gq_precis |>
    as.matrix() |>
    as_tibble(rownames = "Parameter"),
  "./data_store/fragility_models/model_3.3gq_precis.csv")




## Compute fragility curves ----

load_type = c("Strong GM", "Flood")

# Model parameters
b0 <- model_3.3gq_precis[1:5,'mean']
b1 <- model_3.3gq_precis['b1','mean']
b2 <- model_3.3gq_precis['b2','mean']
b3 <- model_3.3gq_precis[8:12,'mean']
cutpoints <- model_3.3gq_precis[13:16,'mean']
delta <- model_3.3gq_precis[17:20,'mean']


# # Define IM range for plotting
# IM_seq <- seq(from=0.02, to=1.0, by=0.02)
# IM_seq <- c(1e-9, IM_seq)
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Number of damage states
n_ds <- length(grep("^cutpoints\\[", rownames(model_3.3gq_precis) )) + 1

# Fix one predictor for the upper and lower bound
thr_lo <- c(1e-9, DT_thres[1:n_ds-1])
thr_up <- DT_thres


# Model x lo

fixed_im <- thr_lo # do not take the log (see Model 2.4.2)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  for (ht in 0:1) {
    
    if (j == 1 && ht == 0) {
      
      # Compute fragility curves
      fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
        tibble(DS = k+1,
               Probability = plogis(
                 b0[j] * lnIM_seq - cutpoints[k] +
                   b1 * sum(c(0, delta)[1:j]) + b2 * ht +
                   b3[j] * fixed_im[j]), IM = IM_seq)
      })
      
      fragility_curves <- bind_rows(fragility_list)
      fragility_curves <- mutate(fragility_curves, DS_part_0=j)
      fragility_curves["Hazard"] = load_type[ht+1]
      
    } else {
      
      # Compute fragility curves
      temp <- purrr::map(seq_along(cutpoints), function(k) {
        tibble(DS = k+1,
               Probability = plogis(
                 b0[j] * lnIM_seq - cutpoints[k] +
                   b1 * sum(c(0, delta)[1:j]) + b2 * ht +
                   b3[j] * fixed_im[j]), IM = IM_seq)
      })
      
      temp <- bind_rows(temp)
      temp <- mutate(temp, DS_part_0=j)
      temp["Hazard"] = load_type[ht+1]
      
      # Remove fragility curves for DS<DSj
      if (j>1) {
        for (k in 1:j) {
          temp <- filter(temp, DS != k)
        }
      }
      
      
      fragility_curves <- bind_rows(fragility_curves,temp)
      
    }
    
  }
  
}


fragility_model_3.3_lo <- fragility_curves |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 3.3 lo")




# Model x up

fixed_im <- thr_up # do not take the log (see Model 2.4.2)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  for (ht in 0:1) {
    
    if (j == 1 && ht == 0) {
      
      # Compute fragility curves
      fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
        tibble(DS = k+1,
               Probability = plogis(
                 b0[j] * lnIM_seq - cutpoints[k] +
                   b1 * sum(c(0, delta)[1:j]) + b2 * ht +
                   b3[j] * fixed_im[j]), IM = IM_seq)
      })
      
      fragility_curves <- bind_rows(fragility_list)
      fragility_curves <- mutate(fragility_curves, DS_part_0=j)
      fragility_curves["Hazard"] = load_type[ht+1]
      
    } else {
      
      # Compute fragility curves
      temp <- purrr::map(seq_along(cutpoints), function(k) {
        tibble(DS = k+1,
               Probability = plogis(
                 b0[j] * lnIM_seq - cutpoints[k] +
                   b1 * sum(c(0, delta)[1:j]) + b2 * ht +
                   b3[j] * fixed_im[j]), IM = IM_seq)
      })
      
      temp <- bind_rows(temp)
      temp <- mutate(temp, DS_part_0=j)
      temp["Hazard"] = load_type[ht+1]
      
      # Remove fragility curves for DS<DSj
      if (j>1) {
        for (k in 1:j) {
          temp <- filter(temp, DS != k)
        }
      }
      
      
      fragility_curves <- bind_rows(fragility_curves,temp)
      
    }
    
  }
  
}

fragility_model_3.3_up <- fragility_curves |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 3.3 up")




figure_m3.3<- ggplot(
  fragility_model_3.3_lo, aes(x = IM, y = Probability,
                              color = DS, linetype = Hazard)) +
  geom_line(linewidth = 1) +
  facet_wrap(
    ~ DS_part_0, ncol = 2,
    labeller = labeller(DS_part_0 = function(x) paste("DS part 1:", x))) +
  scale_color_viridis_d(name = "Damage State") +
  scale_linetype_discrete(name = "Hazard in part 1") +
  labs(
    x = "PGA (g)",
    y = "Probability"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_m3.3)

ggsave(
  plot = figure_m3.3,
  paste0("./data_store/fragility_models/model_3_3_lo.png"),
  width = 21.0,
  height = 29.7*0.67,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot ----

# The state-dependent fragility curves (and some others)
# sdfc <- bind_rows(
#   fragility_model_3.1,
#   fragility_model_3.2,
#   fragility_model_3.3_up,
#   fragility_model_3.3_lo,
# )

sdfc <- bind_rows(
  fragility_model_3.1,
  fragility_model_3.2
)


# The title of the y-axis
y_axis_t <- paste0('$P(DSpt3 \\geq j | DSpt1)$')

# Plot the probabilities of exceeding the thresholds of the damage states
figure_compare <- ggplot(
  sdfc,
  aes( x=IM, y=Probability, color=Model, linetype = Hazard) )+
  geom_line(linewidth = 1)+
  xlab( 'PGA (g)' )+
  ylab( TeX( y_axis_t ) )+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  ) +
  scale_color_viridis_d(end = 0.75) +
  scale_linetype_discrete(name = "Hazard in part 1") +
  # scale_linetype_manual(
  #   values =
  #     setNames( c("solid", "solid", "dotted", "dotted", "dashed", "dashed"),
  #               unique(sdfc$Model) ) )+
  facet_grid(
    cols = vars(DS),
    rows = vars(DS_part_0),
    labeller = labeller(
      DS = function(x) paste0("j = ", x),
      DS_part_0 = function(x) paste0("DSpt1 = ",x)
    )
  )
print(figure_compare)

ggsave(
  plot = figure_compare,
  paste0("./data_store/fragility_models/compare_models_3.1_3.2.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




sdfc <- bind_rows(
  fragility_model_3.1,
  fragility_model_3.3_up,
  fragility_model_3.3_lo,
)


# The title of the y-axis
y_axis_t <- paste0('$P(DSpt3 \\geq j | DSpt1)$')

# Plot the probabilities of exceeding the thresholds of the damage states
figure_compare <- ggplot(
  sdfc,
  aes( x=IM, y=Probability, color=Model, linetype = Hazard) )+
  geom_line(linewidth = 1)+
  xlab( 'PGA (g)' )+
  ylab( TeX( y_axis_t ) )+
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.position = "bottom",
    panel.spacing = unit(1, "lines")
  ) +
  scale_color_viridis_d(end = 0.9) +
  scale_linetype_discrete(name = "Hazard in part 1") +
  # scale_linetype_manual(
  #   values =
  #     setNames( c("solid", "solid", "dotted", "dotted", "dashed", "dashed"),
  #               unique(sdfc$Model) ) )+
  facet_grid(
    cols = vars(DS),
    rows = vars(DS_part_0),
    labeller = labeller(
      DS = function(x) paste0("j = ", x),
      DS_part_0 = function(x) paste0("DSpt1 = ",x)
    )
  )
print(figure_compare)

ggsave(
  plot = figure_compare,
  paste0("./data_store/fragility_models/compare_models_3.1_3.3.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model Comparison ----

reth_comp <- rethinking::compare(model_3.1gq, model_3.2gq, model_3.3gq)
write_csv( reth_comp |>
             as.matrix() |>
             as_tibble(rownames = "Parameter"),
           "./data_store/fragility_models/model_3_comparison_rethinking_logit.csv")

reth_comp_PSIS <- rethinking::compare(
  model_3.1gq, model_3.2gq, model_3.3gq, func = 'PSIS')
write_csv( reth_comp_PSIS |>
             as.matrix() |>
             as_tibble(rownames = "Parameter"),
           "./data_store/fragility_models/model_3_comparison_rethinking_logit_PSIS.csv")


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save Environment ----

# save.image("/gpfs/scratch/trevlopoulos/simulations_CR_LFINF_CDM_11_HEX_2_20240903/frag_curves_pga_haz_eq_totVar_20240903.RData")
save.image("/gpfs/scratch/trevlopoulos/simulations_CR_LFINF_CDM_11_HEX_2_20240903/models_3.x_temp.RData")
