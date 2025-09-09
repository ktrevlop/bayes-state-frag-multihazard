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
#
# If you load an environment from .RData file, execute the following lines up
# until the next Markdown-style comment header


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
# Preprocessing ----

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
sim_step_0 <- filter(sim_step_0, load_type == "strong_gm")
sim_step_0 <- select(sim_step_0, simulation)

sim_step_2 <- filter(comp_thist_results, part_number == 2)
sim_step_2 <- filter(sim_step_2, load_type == "strong_gm")
sim_step_2 <- select(sim_step_2, simulation)

# The indices of the simulations for the selected hazard types in the
# simulation steps 0 and 2
n_of_sim <- inner_join(sim_step_2, sim_step_0)



# Pivot the DS to separate columns per part of the simulation
# ds_per_part_wide will eventually have the data for part_0
ds_per_part_wide <- pivot_wider(comp_thist_results,
                                names_from = part_number, values_from = DS, id_cols = simulation)
# Parts numbered 1-4
ds_per_part_wide <- rename(ds_per_part_wide,
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
)
comp_thist_results_part_0 <- mutate(
  comp_thist_results_part_0,
  lnIM_part_0 = log(IM_part_0))
comp_thist_results_part_0 <- select(
  comp_thist_results_part_0,
  simulation, IM_part_0, lnIM_part_0, max_int_storey_drift_part_0,
  cum_hyst_loop_area, lnTnorm_part_0, lnDeltaT_part_1)

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
model2.x_data_tib <- left_join(ds_per_part_wide, comp_thist_results)
model2.x_data_tib <- select(model2.x_data_tib, -simulation)
model2.x_data_tib$DS_part_0 <- ordered(model2.x_data_tib$DS_part_0)
model2.x_data_tib$DS_part_2 <- ordered(model2.x_data_tib$DS_part_2)
model2.x_data_tib <- filter(model2.x_data_tib, lnIM_part_2 != -Inf)

model2.x_data <- list(
  N = nrow(model2.x_data_tib),
  lnIM_part_0 = model2.x_data_tib$lnIM_part_0,
  lnIM_part_2 = model2.x_data_tib$lnIM_part_2,
  DS_part_0 = as.integer( model2.x_data_tib$DS_part_0 ),
  DS_part_2 = as.integer( model2.x_data_tib$DS_part_2 ),
  lndrift_part_0 = log( model2.x_data_tib$max_int_storey_drift_part_0 ),
  lnhen_part_0 = log( model2.x_data_tib$cum_hyst_loop_area ),
  lnTnorm_part_0 = model2.x_data_tib$lnTnorm_part_0,
  lnDeltaT_part_1 = model2.x_data_tib$lnDeltaT_part_1,
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

source("./model_frag_model_DeltaT_HE_vs_DR.R")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.1 ----

print("Model 1.1")

model_1.1gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    eta <- b0 * lnIM_part_2,
    
    b0 ~ normal(4.0, 2.0),
    cutpoints ~ normal(-4.0, 2.0),
    
    gq> vector[N]: log_lik <-  ordered_logistic_lpmf(
      DS_part_2[i] | b0 * lnIM_part_2[i], cutpoints
    )
    
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_1.1gq <- stancode(model_1.1gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_1.1gq, "model_1_1gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_1_1gq_traceplot.png",
  width = 21.0,
  height = 29.7*0.33,
  units = "cm",
  res = 300
)
traceplot(model_1.1gq)
dev.off()

png(
  filename = "./data_store/checks/model_1_1gq_trankplot.png",
  width = 21.0,
  height = 29.7*0.33,
  units = "cm",
  res = 300
)
trankplot( model_1.1gq )
dev.off()


# Cross-validation and information criteria
model_1.1gq_waic <- WAIC( model_1.1gq )
model_1.1gq_psis <- PSIS( model_1.1gq )

# Summarize results
model_1.1gq_precis <- precis(model_1.1gq, depth = 2)

# Mean parameters
b0 <- model_1.1gq_precis["b0", "mean"]
crow_names <- paste0("cutpoints[", 1:4, "]")
cutpoints <- model_1.1gq_precis[crow_names, "mean"]




## Compute fragility curves ----

# Define IM range for plotting
IM_seq <- seq(from=0.02, to=1.0, by=0.02)
IM_seq <- c(1e-9, IM_seq)
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1


fragility_curves<- purrr::map(seq_along(cutpoints), function(k) {
  tibble(DS = k+1,
         Probability = plogis(
           lnIM_seq * b0 - cutpoints[k]
         ), IM = IM_seq)
})

fragility_model_1.1 <- bind_rows(fragility_curves) |> mutate(DS=factor(DS)) |> 
  mutate(Model = "Model 1.1")

# Plot fragility curves
figure_model_1.1 <- ggplot(fragility_model_1.1,
                           aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DS \\geq j)$") ) +
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_color_viridis_d(name = "j", end = 0.9)
plot(figure_model_1.1)

ggsave(
  paste0("./data_store/fragility_models/model_1_1gq.png"),
  width = 5.0*2.54,
  height = 3.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.1a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 1.1a")

# Compile Stan model with cmdstanr
model_1.1a <- cmdstan_model("model_1_1gq.stan")

# Fit the model
model_1.1a_fit <- model_1.1a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_1.1a_post <- model_1.1a_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_1.1a_post <- select(model_1.1a_post, 2:6)
# Calculate the parameters
model_1.1a_par <- model_1.1a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_1.1a_par, "./data_store/fragility_models/model_1.1a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_1.1a_waic <- waic( model_1.1a_fit$draws("log_lik", format = "matrix") )
# print( model_1.1a_waic$estimates['waic','Estimate'] )
model_1.1a_psis <- loo( model_1.1a_fit$draws("log_lik", format = "matrix") )
# print( model_1.1a_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_1_1a_traceplot <- mcmc_trace(
  model_1.1a_fit$draws(format = "draws_array"),
  pars = c("b0", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 1.1a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1a_traceplot,
  paste0("./data_store/checks/model_1_1a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_1_1a_trankplot <- mcmc_rank_overlay(
  model_1.1a_fit$draws(format = "draws_array"),
  pars = c("b0", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.1a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1a_trankplot,
  paste0("./data_store/checks/model_1_1a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




## Compute fragility curves ----

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

# Posterior mean of b_lnIM_part_2
b0 <- model_1.1a_par$mean[1]

# Posterior mean of cutpoints
cutpoints <- model_1.1a_par$mean[2:5]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1


fragility_curves <- purrr::map(seq_along(cutpoints), function(k) {
  tibble(DS = k+1,
         Probability = plogis(
           lnIM_seq * b0 - cutpoints[k]
         ), IM = IM_seq)
})

fragility_model_1.1a <- bind_rows(fragility_curves) |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 1.1a")




figure_m2.1a <- ggplot(fragility_model_1.1a, aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1)+
  labs(
    x = "PGA (g)",
    y = TeX("$P(DS \\geq j)$")
  )+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_color_viridis_d(name = "j", end = 0.9)
plot(figure_m2.1a)

ggsave(
  plot = figure_m2.1a,
  paste0("./data_store/fragility_models/model_1_1a.png"),
  width = 6.0*2.54,
  height = 6.0*2.54,
  units = c("cm"),
  dpi = 300,
)





#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.1b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 1.1b")

# Compile Stan model with cmdstanr
model_1.1b <- cmdstan_model("model_1_1b.stan")

# Fit the model
model_1.1b_fit <- model_1.1b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_1.1b_post <- model_1.1b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_1.1b_post <- select(model_1.1b_post, 2:6)
# Calculate the parameters
model_1.1b_par <- model_1.1b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_1.1b_par, "./data_store/fragility_models/model_1.1b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_1.1b_waic <- waic( model_1.1b_fit$draws("log_lik", format = "matrix") )
# print( waic_model_1.2b$estimates['waic','Estimate'] )
model_1.1b_psis <- loo( model_1.1b_fit$draws("log_lik", format = "matrix") )
# print( psis_model_1.2b$estimates['looic','Estimate'] )




## Checking chains ----

model_1_1b_traceplot <- mcmc_trace(
  model_1.1b_fit$draws(format = "draws_array"),
  pars = c("b0", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 1.1b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1b_traceplot,
  paste0("./data_store/checks/model_1_1b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_1_1b_trankplot <- mcmc_rank_overlay(
  model_1.1b_fit$draws(format = "draws_array"),
  pars = c("b0", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.1b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1b_trankplot,
  paste0("./data_store/checks/model_1_1b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




## Compute fragility curves ----

# Posterior mean of b_lnIM_part_2
b0 <- model_1.1b_par$mean[1]

# Posterior mean of cutpoints
cutpoints <- model_1.1b_par$mean[2:5]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1


fragility_curves <- purrr::map(seq_along(cutpoints), function(k) {
  tibble(DS = k+1,
         Probability = pnorm(
           lnIM_seq * b0 - cutpoints[k]
         ), IM = IM_seq)
})

fragility_model_1.1b <- bind_rows(fragility_curves) |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 1.1b")




figure_model_1.1b <- ggplot(fragility_model_1.1b, aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1)+
  labs(
    x = "PGA (g)",
    y = TeX("$P(DS \\geq j)$")
  )+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_color_viridis_d(name = "j", end = 0.9)
plot(figure_model_1.1b)

ggsave(
  plot = figure_model_1.1b,
  paste0("./data_store/fragility_models/model_1_1b.png"),
  width = 6.0*2.54,
  height = 6.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.2 ----

print("Model 1.2")

model_1.2gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    eta <- b0 * lnIM_part_0 + b1 * lnIM_part_2,
    
    b0 ~ normal(4, 2),
    b1 ~ normal(4, 2),
    cutpoints ~ normal(-4, 2),
    
    gq> vector[N]: log_lik <-  ordered_logistic_lpmf(
      DS_part_2[i] | b0 * lnIM_part_0[i] + b1 * lnIM_part_2[i],
      cutpoints
    )
    
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_1.2gq <- stancode(model_1.2gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_1.2gq, "model_1_2gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_1_2gq_traceplot.png",
  width = 9.0 * 2.54,
  height = 4.0 * 2.54,
  units = "cm",
  res = 300
)
traceplot( model_1.2gq )
dev.off()

png(
  filename = "./data_store/checks/model_1_2gq_trankplot.png",
  width = 9.0 * 2.54,    # convert cm to inches if needed
  height = 4.0 * 2.54,
  units = "cm",
  res = 300
)
trankplot( model_1.2gq )
dev.off()


# Cross-validation and information criteria
model_1.2gq_waic <- WAIC( model_1.2gq )
model_1.2gq_psis <- PSIS( model_1.2gq )

# Summarize results
model_1.2gq_precis <- precis(model_1.2gq, depth = 2)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.2a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 1.2a")

# Compile Stan model with cmdstanr
model_1.2a <- cmdstan_model("model_1_2gq.stan")

# Fit the model
model_1.2a_fit <- model_1.2a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_1.2a_post <- model_1.2a_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_1.2a_post <- select(model_1.2a_post, 2:7)
# Calculate the parameters
model_1.2a_par <- model_1.2a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_1.2a_par, "./data_store/fragility_models/model_1.2a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_1.2a_waic <- waic( model_1.2a_fit$draws("log_lik", format = "matrix") )
# print( waic_model_2.1b$estimates['waic','Estimate'] )
model_1.2a_psis <- loo( model_1.2a_fit$draws("log_lik", format = "matrix") )
# print( psis_model_1.2a$estimates['looic','Estimate'] )




## Checking chains ----

model_1_2a_traceplot <- mcmc_trace(
  model_1.2a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 1.2a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2a_traceplot,
  paste0("./data_store/checks/model_1_2a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_1_2a_trankplot <- mcmc_rank_overlay(
  model_1.2a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.2a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2a_trankplot,
  paste0("./data_store/checks/model_1_2a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.2b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 1.2b")

# Compile Stan model with cmdstanr
model_1.2b <- cmdstan_model("model_1_2b.stan")

# Fit the model
model_1.2b_fit <- model_1.2b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_1.2b_post <- model_1.2b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_1.2b_post <- select(model_1.2b_post, 2:7)
# Calculate the parameters
model_1.2b_par <- model_1.2b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_1.2b_par, "./data_store/fragility_models/model_1.2b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_1.2b_waic <- waic( model_1.2b_fit$draws("log_lik", format = "matrix") )
# print( waic_model_1.2b$estimates['waic','Estimate'] )
model_1.2b_psis <- loo( model_1.2b_fit$draws("log_lik", format = "matrix") )
# print( psis_model_1.2b$estimates['looic','Estimate'] )




## Checking chains ----

model_1_2b_traceplot <- mcmc_trace(
  model_1.2b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 1.2b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2b_traceplot,
  paste0("./data_store/checks/model_1_2b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_1_2b_trankplot <- mcmc_rank_overlay(
  model_1.2b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.2b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2b_trankplot,
  paste0("./data_store/checks/model_1_2b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.1 ----

print("model_2.1gq")

model_2.1gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    eta <- b0 * lnIM_part_2 + b1 * sum( append_row(0, delta)[1:DS_part_0] ),
    
    b0 ~ normal(4, 2),
    b1 ~ normal(0, 2.5),
    cutpoints ~ normal(-4, 2),
    simplex[4]: delta ~ dirichlet( alpha ),
    
    gq> vector[N]: log_lik <-  ordered_logistic_lpmf(
      DS_part_2[i] | b0 * lnIM_part_2[i] +
        b1 * sum(append_row(0, delta)[1:DS_part_0[i]]), cutpoints
    )
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.1gq <- stancode(model_2.1gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.1gq, "model_2_1gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_2_1gq_traceplot.png",
  width = 9.0 * 2.54,
  height = 6.0 * 2.54,
  units = "cm",
  res = 300
)
traceplot( model_2.1gq )
dev.off()

png(
  filename = "./data_store/checks/model_2_1gq_trankplot.png",
  width = 9.0 * 2.54,
  height = 6.0 * 2.54,
  units = "cm",
  res = 300
)
trankplot( model_2.1gq )
dev.off()


# Information criteria
model_2.1gq_waic <- WAIC( model_2.1gq )
model_2.1gq_psis <- PSIS( model_2.1gq )

# Summarize results
model_2.1gq_precis <- precis(model_2.1gq, depth = 2)



## Compute fragility curves ----

# # Extract posterior samples
# post_model_2.1 <- extract.samples(model_2.1gq)

# Posterior mean of b_lnIM_part_2
# b_lnIM_part_2_mean <- mean(post_model_2.1$b_lnIM_part_2)
b0 <- model_2.1gq_precis[1,'mean']

# Posterior mean of b_DS
# b_DS_mean <- mean(post_model_2.1$b_DS)
b1 <- model_2.1gq_precis[2,'mean']

# Posterior mean of cutpoints
# cutpoints_model_2.1 <- apply(post_model_2.1$cutpoints, 2, mean)
cutpoints <- model_2.1gq_precis[3:6,'mean']

# Posterior mean of delta
# delta_mean_model_2.1 <- colMeans(post_model_2.1$delta)
delta <- model_2.1gq_precis[7:10,'mean']




# Define IM range for plotting
IM_seq <- seq(from=0.02, to=1.0, by=0.02)
IM_seq <- c(1e-9, IM_seq)
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:n_ds) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0 - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0 - cutpoints[k] +
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

fragility_model_2.1a <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.1")




figure_m2.3 <- ggplot(fragility_model_2.1a, aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(~ DS_part_0, ncol = 2,
             labeller = labeller(DS_part_0 = function(x) paste("DSpt1 =", x))) +
  scale_color_viridis_d(name = "j", end = 0.9) +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DSpt3 \\geq j|DSpt1)$")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_m2.3)

ggsave(
  plot = figure_m2.3,
  paste0("./data_store/fragility_models/model_2_1.png"),
  width = 6.0*2.54,
  height = 6.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.1a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.1a")

# Compile Stan model with cmdstanr
model_2.1a <- cmdstan_model("model_2_1gq.stan")

# Fit the model
model_2.1a_fit <- model_2.1a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.1a_post <- model_2.1a_fit$draws(variables = NULL,
                                        format = "data.frame")
# Select the parameters
model_2.1a_post <- select(model_2.1a_post, 2:11)
# Calculate the parameters
model_2.1a_par <- model_2.1a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.1a_par, "./data_store/fragility_models/model_2.1a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.1a_waic <- waic( model_2.1a_fit$draws("log_lik", format = "matrix") )
# print( waic_model_2.1b$estimates['waic','Estimate'] )
model_2.1a_psis <- loo( model_2.1a_fit$draws("log_lik", format = "matrix") )
# print( psis_model_1.2a$estimates['looic','Estimate'] )




## Checking chains ----

model_2_1a_traceplot <- mcmc_trace(
  model_2.1a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.1a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_1a_traceplot,
  paste0("./data_store/checks/model_2_1a_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_1a_trankplot <- mcmc_rank_overlay(
  model_2.1a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.1a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_1a_trankplot,
  paste0("./data_store/checks/model_2_1a_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




## Compute fragility curves ----

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

# Posterior mean of b_lnIM_part_2
b0 <- model_2.1a_par$mean[1]

# Posterior mean of b_DS
b1 <- model_2.1a_par$mean[2]

# Posterior mean of cutpoints
cutpoints <- model_2.1a_par$mean[3:6]

# Posterior mean of delta
delta <- model_2.1a_par$mean[7:10]




# Define IM range for plotting
IM_seq <- seq(from=0.02, to=1.0, by=0.02)
IM_seq <- c(1e-9, IM_seq)
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0 - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0 - cutpoints[k] +
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


fragility_model_2.1a <- fragility_curves |> mutate(DS=factor(DS)) |> 
  mutate(Model = "Model 2.1a")




figure_model_2.1a <- ggplot(
  fragility_model_2.1a,
  aes(x = IM, y = Probability, color = DS))+
  geom_line(linewidth = 1)+
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(~ DS_part_0, ncol = 2,
             labeller = labeller(DS_part_0 = function(x) paste("DSpt1 = ", x))) +
  scale_color_viridis_d(name = "j", end = 0.9) +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DS \\geq j|DSpt1)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_model_2.1a)

ggsave(
  plot = figure_model_2.1a,
  paste0("./data_store/fragility_models/model_2_1a.png"),
  width = 6.0*2.54,
  height = 6.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.1b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 2.1b")

# Compile Stan model with cmdstanr
model_2.1b <- cmdstan_model("model_2_1b.stan")

# Fit the model
model_2.1b_fit <- model_2.1b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.1b_post <- model_2.1b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.1b_post <- select(model_2.1b_post, 2:11)
# Calculate the parameters
model_2.1b_par <- model_2.1b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.1b_par, "./data_store/fragility_models/model_2.1b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.1b_waic <- waic( model_2.1b_fit$draws("log_lik", format = "matrix") )
# print( waic_model_2.1b$estimates['waic','Estimate'] )
model_2.1b_psis <- loo( model_2.1b_fit$draws("log_lik", format = "matrix") )
# print( psis_model_2.1b$estimates['looic','Estimate'] )




## Checking chains ----

model_2_1b_traceplot <- mcmc_trace(
  model_2.1b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.1b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_1b_traceplot,
  paste0("./data_store/checks/model_2_1b_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_1b_trankplot <- mcmc_rank_overlay(
  model_2.1b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.1b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_1b_trankplot,
  paste0("./data_store/checks/model_2_1b_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)



## Compute fragility curves ----

# Posterior mean of b_lnIM_part_2
b0 <- model_2.1b_par$mean[1]

# Posterior mean of b_DS
b1 <- model_2.1b_par$mean[2]

# Posterior mean of cutpoints
cutpoints <- model_2.1b_par$mean[3:6]

# Posterior mean of delta
delta <- model_2.1b_par$mean[7:10]




# Define IM range for plotting
IM_seq <- seq(from=0.02, to=1.0, by=0.02)
IM_seq <- c(1e-9, IM_seq)
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = pnorm(
               lnIM_seq * b0 - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = pnorm(
               lnIM_seq * b0 - cutpoints[k] +
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


fragility_model_2.1b <- fragility_curves |> mutate(DS=factor(DS)) |> 
  mutate(Model = "Model 2.1b")




figure_model_2.1b <- ggplot(
  fragility_model_2.1b,
  aes(x = IM, y = Probability, color = DS))+
  geom_line(linewidth = 1)+
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(~ DS_part_0, ncol = 2,
             labeller = labeller(DS_part_0 = function(x) paste("DSpt1 = ", x))) +
  scale_color_viridis_d(name = "j", end = 0.9) +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DS \\geq j|DSpt1)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_model_2.1b)

ggsave(
  plot = figure_model_2.1b,
  paste0("./data_store/fragility_models/model_2_1b.png"),
  width = 6.0*2.54,
  height = 6.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.2 ----

print("model_2.2gq")

model_2.2gq <- ulam(
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
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.2gq <- stancode(model_2.2gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.2gq, "model_2_2gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_2_2gq_traceplot.png",
  width = 9.0 * 2.54,
  height = 6.0 * 2.54,
  units = "cm",
  res = 300
)
traceplot( model_2.2gq )
dev.off()

png(
  filename = "./data_store/checks/model_2_2gq_trankplot.png",
  width = 9.0 * 2.54,
  height = 6.0 * 2.54,
  units = "cm",
  res = 300
)
trankplot( model_2.2gq )
dev.off()


# Information criteria
model_2.2gq_waic <- WAIC( model_2.2gq )
model_2.2gq_psis <- PSIS( model_2.2gq )

# Summarize results
model_2.2gq_precis <- precis(model_2.2gq, depth = 2)

# Export the parameters
write_csv(
  model_2.2gq_precis |>
    as.matrix() |>
    as_tibble(rownames = "Parameter"),
  "./data_store/fragility_models/model_2.2gq_precis.csv")



## Compute fragility curves ----

# Posterior mean of b_lnIM_part_2
b0 <- model_2.2gq_precis[1:5,'mean']

# Posterior mean of b_DS
b1 <- model_2.2gq_precis[6,'mean']

# Posterior mean of cutpoints
cutpoints <- model_2.2gq_precis[7:10,'mean']

# Posterior mean of delta
delta <- model_2.2gq_precis[11:14,'mean']



# # Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0[j] - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0[j] - cutpoints[k] +
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


fragility_model_2.2gq <- fragility_curves |> mutate(DS = factor(DS)) |>
  mutate(Model = "Model 2.2")




figure_model_2.2gq <- ggplot(fragility_model_2.2gq, aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(~ DS_part_0, ncol = 2,
             labeller = labeller(DS_part_0 = function(x) paste("DSpt1 = ", x))) +
  scale_color_viridis_d(name = "j", end = 0.9) +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DSpt3 \\geq j|DSpt1)$")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_model_2.2gq)

ggsave(
  plot = figure_model_2.2gq,
  paste0("./data_store/fragility_models/model_2_2gq.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.2a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.2a")

# Compile Stan model with cmdstanr
model_2.2a <- cmdstan_model("model_2_2gq.stan")

# Fit the model
model_2.2a_fit <- model_2.2a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.2a_post <- model_2.2a_fit$draws(variables = NULL,
                                        format = "data.frame")
# Select the parameters
model_2.2a_post <- select(model_2.2a_post, 2:15)
# Calculate the parameters
model_2.2a_par <- model_2.2a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.2a_par, "./data_store/fragility_models/model_2.2a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.2a_waic <- waic( model_2.2a_fit$draws("log_lik", format = "matrix") )
# print( model_2.2a_waic$estimates['waic','Estimate'] )
model_2.2a_psis <- loo( model_2.2a_fit$draws("log_lik", format = "matrix") )
# print( model_2.2a_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_2a_traceplot <- mcmc_trace(
  model_2.2a_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.2a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_2a_traceplot,
  paste0("./data_store/checks/model_2_2a_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_2a_trankplot <- mcmc_rank_overlay(
  model_2.2a_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]",  "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.2a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_2a_trankplot,
  paste0("./data_store/checks/model_2_2a_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




## Compute fragility curves ----

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

# Posterior mean of b_lnIM_part_2
b0 <- model_2.2a_par$mean[1:5]

# Posterior mean of b_DS
b1 <- model_2.2a_par$mean[6]

# Posterior mean of cutpoints
cutpoints <- model_2.2a_par$mean[7:10]

# Posterior mean of delta
delta <- model_2.2a_par$mean[11:14]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1


# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0[j] - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               lnIM_seq * b0[j] - cutpoints[k] +
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


fragility_model_2.2a <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.2a")




figure_model_2.2a <- ggplot(fragility_model_2.2a,
                       aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1)+
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(~ DS_part_0, ncol = 2,
             labeller = labeller(DS_part_0 = function(x) paste("DSpt1: = ", x)))+
  scale_color_viridis_d(name = "j") +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DSpt3 \\geq j|DSpt1)$")
  )+
  theme_minimal()+
  theme(legend.position = "bottom")
plot(figure_model_2.2a)

ggsave(
  plot = figure_model_2.2a,
  paste0("./data_store/fragility_models/model_2_2a.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.2b (stan; probit) ----

# This model uses an ordered probit link function

# Compile Stan model with cmdstanr
model_2.2b <- cmdstan_model("model_2_2b.stan")

# Fit the model
model_2.2b_fit <- model_2.2b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.2b_post <- model_2.2b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.2b_post <- select(model_2.2b_post, 2:15)
# Calculate the parameters
model_2.2b_par <- model_2.2b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.2b_par, "./data_store/fragility_models/model_2.2b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.2b_waic <- waic( model_2.2b_fit$draws("log_lik", format = "matrix") )
# print( model_2.2b_waic$estimates['waic','Estimate'] )
model_2.2b_psis <- loo( model_2.2b_fit$draws("log_lik", format = "matrix") )
# print( model_2.2b_psis$estimates['looic','Estimate'] )


## Checking chains ----

model_2_2b_traceplot <- mcmc_trace(
  model_2.2b_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.2b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_2b_traceplot,
  paste0("./data_store/checks/model_2_2b_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_2b_trankplot <- mcmc_rank_overlay(
  model_2.2b_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.2b")+
  ylim(20, 80)+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_2b_trankplot,
  paste0("./data_store/checks/model_2_2b_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




## Compute fragility curves ----

# Posterior mean of b_lnIM_part_2
b0 <- model_2.2b_par$mean[1:5]

# Posterior mean of b_DS
b1 <- model_2.2b_par$mean[6]

# Posterior mean of cutpoints
cutpoints <- model_2.2b_par$mean[7:10]

# Posterior mean of delta
delta <- model_2.2b_par$mean[11:14]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1


# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = pnorm(
               lnIM_seq * b0[j] - cutpoints[k] +
                 b1 * sum(c(0, delta)[1:j])), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = pnorm(
               lnIM_seq * b0[j] - cutpoints[k] +
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


fragility_model_2.2b <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.2b")




figure_model_2.2b <- ggplot(fragility_model_2.2a,
                            aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1)+
  # facet_wrap(~ DS_part_0, ncol = 2, labeller = label_both) +
  facet_wrap(~ DS_part_0, ncol = 2,
             labeller = labeller(DS_part_0 = function(x) paste("DSpt1: = ", x)))+
  scale_color_viridis_d(name = "j") +
  labs(
    x = "PGA (g)",
    y = TeX("$P(DSpt3 \\geq j|DSpt1)$")
  )+
  theme_minimal()+
  theme(legend.position = "bottom")
plot(figure_model_2.2b)

ggsave(
  plot = figure_model_2.2b,
  paste0("./data_store/fragility_models/model_2_2b.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.3 ----

print("model_2.3gq")

model_2.3gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear predictor
    eta <- b0[DS_part_0] * lnIM_part_2 +
      b1 * sum(append_row(0, delta)[1:DS_part_0]) +
      b2[DS_part_0] * lnhen_part_0,
    
    # Priors
    vector[max(DS_part_0)]: b0 ~ normal(4, 2),
    b1 ~ normal(0, 2.5),
    vector[max(DS_part_0)]: b2 ~ normal(0, 2.5),
    cutpoints ~ normal(-4, 2),
    simplex[4]: delta ~ dirichlet(alpha),

    # Log-likelihood for WAIC/LOO
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] |
        b0[DS_part_0[i]] * lnIM_part_2[i] +
        b1 * sum(append_row(0, delta)[1:DS_part_0[i]]) +
        b2[DS_part_0[i]] * lnhen_part_0[i],
      cutpoints
    )
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.3gq <- stancode(model_2.3gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.3gq, "model_2_3gq.stan")
 

## Checking chains ----

png(
  filename = "./data_store/checks/model_2_3gq_traceplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
traceplot(
  model_2.3gq, n_cols = 4, max_rows = 10,
  pars = c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b1",
           "b2[1]", "b2[2]", "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]", "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]", "delta[3]", "delta[4]"))
dev.off()

png(
  filename = "./data_store/checks/model_2_3gq_trankplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
trankplot(
  model_2.3gq, n_cols = 4, max_rows = 10,
  pars = c("b0[1]", "b0[2]", "b0[3]", "b0[4]", "b1",
           "b2[1]", "b2[2]", "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]", "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]", "delta[3]", "delta[4]"))
dev.off()


# Information criteria
model_2.3gq_waic <- WAIC( model_2.3gq )
model_2.3gq_psis <- PSIS( model_2.3gq )

# Summarize results
model_2.3gq_precis <- precis(model_2.3gq, depth = 2)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.3a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.3a")

# Compile Stan model with cmdstanr
model_2.3a <- cmdstan_model("model_2_3gq.stan")

# Fit the model
model_2.3a_fit <- model_2.3a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.3a_post <- model_2.3a_fit$draws(variables = NULL,
                                            format = "data.frame")

# Select the parameters
model_2.3a_post <- select(model_2.3a_post, 2:20)
# Calculate the parameters
model_2.3a_par <- model_2.3a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.3a_par, "./data_store/fragility_models/model_2.3a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.3a_waic <- waic( model_2.3a_fit$draws("log_lik", format = "matrix") )
# print( model_2.3a_waic$estimates['waic','Estimate'] )
model_2.3a_psis <- loo( model_2.3a_fit$draws("log_lik", format = "matrix") )
# print( model_2.3a_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_3a_traceplot <- mcmc_trace(
  model_2.3a_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.3a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_3a_traceplot,
  paste0("./data_store/checks/model_2_3a_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_3a_trankplot <- mcmc_rank_overlay(
  model_2.3a_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]",  "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.3a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_3a_trankplot,
  paste0("./data_store/checks/model_2_3a_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.3b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 2.3b")

# Compile Stan model with cmdstanr
model_2.3b <- cmdstan_model("model_2_3b.stan")

# Fit the model
model_2.3b_fit <- model_2.3b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.3b_post <- model_2.3b_fit$draws(variables = NULL,
                                            format = "data.frame")

# Select the parameters
model_2.3b_post <- select(model_2.3b_post, 2:20)
# Calculate the parameters
model_2.3b_par <- model_2.3b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.3b_par, "./data_store/fragility_models/model_2.3b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.3b_waic <- waic( model_2.3b_fit$draws("log_lik", format = "matrix") )
# print( model_2.3b_waic$estimates['waic','Estimate'] )
model_2.3b_psis <- loo( model_2.3b_fit$draws("log_lik", format = "matrix") )
# print( model_2.3b_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_3b_traceplot <- mcmc_trace(
  model_2.3b_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.3b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_3b_traceplot,
  paste0("./data_store/checks/model_2_3b_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_3b_trankplot <- mcmc_rank_overlay(
  model_2.3b_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]",  "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.3b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_3b_trankplot,
  paste0("./data_store/checks/model_2_3b_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.4 ----

print("model_2.4gq")

model_2.4gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear predictor
    eta <- b0[DS_part_0] * lnIM_part_2 +
      b1 * sum(append_row(0, delta)[1:DS_part_0]) +
      b2[DS_part_0] * exp(lnDeltaT_part_1),
    
    # Priors
    vector[max(DS_part_0)]: b0 ~ normal(4, 2),
    b1 ~ normal(0, 2.5),
    vector[max(DS_part_0)]: b2 ~ normal(0, 2.5),
    cutpoints ~ normal(-4, 2),
    simplex[4]: delta ~ dirichlet(alpha),
    
    # Log-likelihood for WAIC/LOO
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] |
        b0[DS_part_0[i]] * lnIM_part_2[i] +
        b1 * sum(append_row(0, delta)[1:DS_part_0[i]]) +
        b2[DS_part_0[i]] * exp(lnDeltaT_part_1[i]),
      cutpoints
    )
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.4gq <- stancode(model_2.4gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.4gq, "model_2_4gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_2_4gq_traceplot.png",
  width = 29.7,
  height = 21.0,
  units = "cm",
  res = 300
)
traceplot( model_2.4gq, n_cols = 4, max_rows = 10,
           pars = c("b0[1]", "b0[2]", "b0[3]",
                    "b0[4]", "b0[5]", "b0[1]",
                    "b2[1]", "b2[2]", "b2[3]",
                    "b2[4]", "b2[5]", "cutpoints[1]",
                    "cutpoints[2]", "cutpoints[3]", "cutpoints[4]",
                    "delta[1]", "delta[2]", "delta[3]",
                    "delta[4]"))
dev.off()

png(
  filename = "./data_store/checks/model_2_4gq_trankplot.png",
  width = 29.7,    # convert cm to inches if needed
  height = 21.0,
  units = "cm",
  res = 300
)
trankplot( model_2.4gq, n_cols = 4, max_rows = 10,
           pars = c("b0[1]", "b0[2]", "b0[3]",
                    "b0[4]", "b0[5]", "b2[1]",
                    "b2[1]", "b2[2]", "b2[3]",
                    "b2[4]", "b2[5]", "cutpoints[1]",
                    "cutpoints[2]", "cutpoints[3]", "cutpoints[4]",
                    "delta[1]", "delta[2]", "delta[3]",
                    "delta[4]"))
dev.off()


# Information criteria
model_2.4gq_waic <- WAIC( model_2.4gq )
model_2.4gq_psis <- PSIS( model_2.4gq )

# Summarize results
model_2.4gq_precis <- precis(model_2.4gq, depth = 2)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.4a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.4a")

# Compile Stan model with cmdstanr
model_2.4a <- cmdstan_model("model_2_4gq.stan")

# Fit the model
model_2.4a_fit <- model_2.4a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.4a_post <- model_2.4a_fit$draws(variables = NULL,
                                            format = "data.frame")

# Select the parameters
model_2.4a_post <- select(model_2.4a_post, 2:20)
# Calculate the parameters
model_2.4a_par <- model_2.4a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.4a_par, "./data_store/fragility_models/model_2.4a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.4a_waic <- waic( model_2.4a_fit$draws("log_lik", format = "matrix") )
# print( model_2.4a_waic$estimates['waic','Estimate'] )
model_2.4a_psis <- loo( model_2.4a_fit$draws("log_lik", format = "matrix") )
# print( model_2.4a_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_4a_traceplot <- mcmc_trace(
  model_2.4a_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.3a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_4a_traceplot,
  paste0("./data_store/checks/model_2_4a_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_4a_trankplot <- mcmc_rank_overlay(
  model_2.4a_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]",  "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.4a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_4a_trankplot,
  paste0("./data_store/checks/model_2_4a_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.4b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 2.4b")

# Compile Stan model with cmdstanr
model_2.4b <- cmdstan_model("model_2_4b.stan")

# Fit the model
model_2.4b_fit <- model_2.4b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.4b_post <- model_2.4b_fit$draws(variables = NULL,
                                            format = "data.frame")

# Select the parameters
model_2.4b_post <- select(model_2.4b_post, 2:20)
# Calculate the parameters
model_2.4b_par <- model_2.4b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.4b_par, "./data_store/fragility_models/model_2.4b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.4b_waic <- waic( model_2.4b_fit$draws("log_lik", format = "matrix") )
# print( model_2.4b_waic$estimates['waic','Estimate'] )
model_2.4b_psis <- loo( model_2.4b_fit$draws("log_lik", format = "matrix") )
# print( model_2.4b_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_4b_traceplot <- mcmc_trace(
  model_2.4b_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]", "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+
  labs(title = "Model 2.3b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_4b_traceplot,
  paste0("./data_store/checks/model_2_4b_traceplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)


model_2_4b_trankplot <- mcmc_rank_overlay(
  model_2.4b_fit$draws(format = "draws_array"),
  pars = c("b0[1]", "b0[2]",
           "b0[3]", "b0[4]",  "b1",
           "b2[1]", "b2[2]",
           "b2[3]", "b2[4]",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]",
           "delta[1]", "delta[2]",
           "delta[3]", "delta[4]"),
  facet_args = list(ncol = 4, nrow = 5) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.4b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_4b_trankplot,
  paste0("./data_store/checks/model_2_4b_trankplot.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.5 ----

print("model_2.5gq")

model_2.5gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear model using the centered coefficients
    eta <- b0 * lnIM_part_2 + b1 * lndrift_part_0,
    
    # Priors
    b0 ~ lognormal(4, 2),
    b1 ~ lognormal(4, 2),
    cutpoints ~ normal(-4, 2),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] | b0 * lnIM_part_2[i] + b1 * lndrift_part_0[i], cutpoints
    )
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.5gq <- stancode(model_2.5gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.5gq, "model_2_5gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_2_5gq_traceplot.png",
  width = 9.0 * 2.54,
  height = 6.0 * 2.54,
  units = "cm",
  res = 300
)
traceplot( model_2.5gq )
dev.off()

png(
  filename = "./data_store/checks/model_2_5gq_trankplot.png",
  width = 9.0 * 2.54,    # convert cm to inches if needed
  height = 6.0 * 2.54,
  units = "cm",
  res = 300
)
trankplot( model_2.5gq )
dev.off()


# Cross-validation and information criteria
model_2.5gq_waic <- WAIC( model_2.5gq )
model_2.5gq_psis <- PSIS( model_2.5gq )

# Summarize results
model_2.5gq_precis <- precis(model_2.5gq, depth = 2)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.5a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.5a")

# Compile Stan model with cmdstanr
model_2.5a <- cmdstan_model("model_2_5gq.stan")

# Fit the model
model_2.5a_fit <- model_2.5a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.5a_post <- model_2.5a_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.5a_post <- select(model_2.5a_post, 2:7)
# Calculate the parameters
model_2.5a_par <- model_2.5a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.5a_par, "./data_store/fragility_models/model_2.5a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.5a_waic <- waic( model_2.5a_fit$draws("log_lik", format = "matrix") )
# print( model_2.5a_waic$estimates['waic','Estimate'] )
model_2.5a_psis <- loo( model_2.5a_fit$draws("log_lik", format = "matrix") )
# print( model_2.5a_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_5a_traceplot <- mcmc_trace(
  model_2.5a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 2.5a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_5a_traceplot,
  paste0("./data_store/checks/model_2_5a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_2_5a_trankplot <- mcmc_rank_overlay(
  model_2.5a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.5a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_5a_trankplot,
  paste0("./data_store/checks/model_2_5a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.5b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 2.5b")

# Compile Stan model with cmdstanr
model_2.5b <- cmdstan_model("model_2_5b.stan")

# Fit the model
model_2.5b_fit <- model_2.5b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.5b_post <- model_2.5b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.5b_post <- select(model_2.5b_post, 2:7)
# Calculate the parameters
model_2.5b_par <- model_2.5b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.5b_par, "./data_store/fragility_models/model_2.5b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.5b_waic <- waic( model_2.5b_fit$draws("log_lik", format = "matrix") )
# print( model_2.5b_waic$estimates['waic','Estimate'] )
model_2.5b_psis <- loo( model_2.5b_fit$draws("log_lik", format = "matrix") )
# print( model_2.5b_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_5b_traceplot <- mcmc_trace(
  model_2.5b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 2.5b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_5b_traceplot,
  paste0("./data_store/checks/model_2_5b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_2_5b_trankplot <- mcmc_rank_overlay(
  model_2.5b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.5b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_5b_trankplot,
  paste0("./data_store/checks/model_2_5b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.6 ----

print("model_2.6gq")

model_2.6gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear model using the centered coefficients
    eta <- b0 * lnIM_part_2 + b1 * lnhen_part_0,
    
    # Priors
    b0 ~ lognormal(4, 2),
    b1 ~ lognormal(4, 2),
    cutpoints ~ normal(-4, 2),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] | b0 * lnIM_part_2[i] + b1 * lnhen_part_0[i], cutpoints
    )
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.6gq <- stancode(model_2.6gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.6gq, "model_2_6gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_2_6gq_traceplot.png",
  width = 21.0,
  height = 29.7 * 0.25,
  units = "cm",
  res = 300
)
traceplot( model_2.6gq )
dev.off()

png(
  filename = "./data_store/checks/model_2_6gq_trankplot.png",
  width = 21.0,
  height = 29.7 * 0.25,
  units = "cm",
  res = 300
)
trankplot( model_2.6gq )
dev.off()


# Cross-validation and information criteria
model_2.6gq_waic <- WAIC( model_2.6gq )
model_2.6gq_psis <- PSIS( model_2.6gq )

# Summarize results
model_2.6gq_precis <- precis(model_2.6gq, depth = 2)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.6a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.6a")

# Compile Stan model with cmdstanr
model_2.6a <- cmdstan_model("model_2_6gq.stan")

# Fit the model
model_2.6a_fit <- model_2.6a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.6a_post <- model_2.6a_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.6a_post <- select(model_2.6a_post, 2:7)
# Calculate the parameters
model_2.6a_par <- model_2.6a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.6a_par, "./data_store/fragility_models/model_2.6a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.6a_waic <- waic( model_2.6a_fit$draws("log_lik", format = "matrix") )
# print( model_2.6a_waic$estimates['waic','Estimate'] )
model_2.6a_psis <- loo( model_2.6a_fit$draws("log_lik", format = "matrix") )
# print( model_2.6a_psis$estimates['looic','Estimate'] )





## Checking chains ----

model_2_6a_traceplot <- mcmc_trace(
  model_2.6a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 2.6a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_6a_traceplot,
  paste0("./data_store/checks/model_2_6a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_2_6a_trankplot <- mcmc_rank_overlay(
  model_2.6a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.6a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_6a_trankplot,
  paste0("./data_store/checks/model_2_6a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.6b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 2.6b")

# Compile Stan model with cmdstanr
model_2.6b <- cmdstan_model("model_2_6b.stan")

# Fit the model
model_2.6b_fit <- model_2.6b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.6b_post <- model_2.6b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.6b_post <- select(model_2.6b_post, 2:7)
# Calculate the parameters
model_2.6b_par <- model_2.6b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.6b_par, "./data_store/fragility_models/model_2.6b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.6b_waic <- waic( model_2.6b_fit$draws("log_lik", format = "matrix") )
# print( model_2.6b_waic$estimates['waic','Estimate'] )
model_2.6b_psis <- loo( model_2.6b_fit$draws("log_lik", format = "matrix") )
# print( model_2.6b_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_6b_traceplot <- mcmc_trace(
  model_2.6b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 2.6b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_6b_traceplot,
  paste0("./data_store/checks/model_2_6b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_2_6b_trankplot <- mcmc_rank_overlay(
  model_2.6b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.6b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_6b_trankplot,
  paste0("./data_store/checks/model_2_6b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.7 ----

print("model_2.7gq")

model_2.7gq <- ulam(
  alist(
    DS_part_2 ~ ordered_logistic(eta, cutpoints),
    
    # Linear model using the centered coefficients
    eta <- b0 * lnIM_part_2 + b1 * exp(lnDeltaT_part_1),
    
    # Priors
    b0 ~ lognormal(4, 2),
    b1 ~ lognormal(4, 2),
    cutpoints ~ normal(-4, 2),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_2[i] | b0 * lnIM_part_2[i] + b1 * exp(lnDeltaT_part_1[i]),
      cutpoints
    )
  ),
  data = model2.x_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_2.7gq <- stancode(model_2.7gq)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_2.7gq, "model_2_7gq.stan")


## Checking chains ----

png(
  filename = "./data_store/checks/model_2_7gq_traceplot.png",
  width = 21.0,
  height = 29.7 * 0.25,
  units = "cm",
  res = 300
)
traceplot( model_2.7gq )
dev.off()

png(
  filename = "./data_store/checks/model_2_7gq_trankplot.png",
  width = 21.0,
  height = 29.7 * 0.25,
  units = "cm",
  res = 300
)
trankplot( model_2.7gq )
dev.off()


# Cross-validation and information criteria
model_2.7gq_waic <- WAIC( model_2.7gq )
model_2.7gq_psis <- PSIS( model_2.7gq )

# Summarize results
model_2.7gq_precis <- precis(model_2.7gq, depth = 2)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.7a (stan; logit) ----

# This model uses an ordered logit link function

print("Model 2.7a")

# Compile Stan model with cmdstanr
model_2.7a <- cmdstan_model("model_2_7gq.stan")

# Fit the model
model_2.7a_fit <- model_2.7a$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.7a_post <- model_2.7a_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.7a_post <- select(model_2.7a_post, 2:7)
# Calculate the parameters
model_2.7a_par <- model_2.7a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.7a_par, "./data_store/fragility_models/model_2.7a_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.7a_waic <- waic( model_2.7a_fit$draws("log_lik", format = "matrix") )
# print( model_2.7a_waic$estimates['waic','Estimate'] )
model_2.7a_psis <- loo( model_2.7a_fit$draws("log_lik", format = "matrix") )
# print( model_2.7a_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_7a_traceplot <- mcmc_trace(
  model_2.7a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 2.7a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_7a_traceplot,
  paste0("./data_store/checks/model_2_7a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_2_7a_trankplot <- mcmc_rank_overlay(
  model_2.7a_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.7a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_7a_trankplot,
  paste0("./data_store/checks/model_2_7a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 2.7b (stan; probit) ----

# This model uses an ordered logit link function

print("Model 2.7b")

# Compile Stan model with cmdstanr
model_2.7b <- cmdstan_model("model_2_7b.stan")

# Fit the model
model_2.7b_fit <- model_2.7b$sample(
  data = model2.x_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_2.7b_post <- model_2.7b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_2.7b_post <- select(model_2.7b_post, 2:7)
# Calculate the parameters
model_2.7b_par <- model_2.7b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_2.7b_par, "./data_store/fragility_models/model_2.7b_par.csv")


# Cross-validation and information criteria
# See loo output for the Pareto k diagnostic values
model_2.7b_waic <- waic( model_2.7b_fit$draws("log_lik", format = "matrix") )
# print( model_2.7b_waic$estimates['waic','Estimate'] )
model_2.7b_psis <- loo( model_2.7b_fit$draws("log_lik", format = "matrix") )
# print( model_2.7b_psis$estimates['looic','Estimate'] )




## Checking chains ----

model_2_7b_traceplot <- mcmc_trace(
  model_2.7b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+
  labs(title = "Model 2.7a")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_7b_traceplot,
  paste0("./data_store/checks/model_2_7b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)


model_2_7b_trankplot <- mcmc_rank_overlay(
  model_2.7b_fit$draws(format = "draws_array"),
  pars = c("b0", "b1", "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2) )+ 
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 2.7b")+
  theme_minimal()+
  theme(legend.position = "bottom")+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_2_7b_trankplot,
  paste0("./data_store/checks/model_2_7b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.5,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model Comparison ----

reth_comp <- rethinking::compare(model_1.1gq, model_1.2gq, model_2.1gq,
                                 model_2.2gq, model_2.3gq, model_2.4gq,
                                 model_2.5gq, model_2.6gq, model_2.7gq)
# Use write.csv to write the row names
write.csv(reth_comp, "./data_store/fragility_models/model_comparison_rethinking_logit.csv")

# loo_compare(model_2.4a_waic, model_2.2b_waic)

# weight_calc( c(
# model_2.4a_waic$estimates["waic","Estimate"],
# model_2.2b_waic$estimates["waic","Estimate"]
# ) )


comparison_tb = tibble(Model = c("Model 1.1", "Model 1.2", "Model 2.1", "Model 2.2",
                                 "Model 2.3", "Model 2.4",
                                 "Model 2.5", "Model 2.6", "Model 2.7"))
# comparison_tb <- comparison_tb |> mutate(model_n=seq(1,nrow(comparison_tb))) |> 
#   select(model_n, Model)

# Add WAIC and PSIS to the table
comparison_tb['WAIC_logit'] = c(
  model_1.1a_waic$estimates["waic","Estimate"],
  model_1.2a_waic$estimates["waic","Estimate"],
  model_2.1a_waic$estimates["waic","Estimate"],
  model_2.4a_waic$estimates["waic","Estimate"],
  model_2.3a_waic$estimates["waic","Estimate"],
  model_2.4a_waic$estimates["waic","Estimate"],
  model_2.5a_waic$estimates["waic","Estimate"],
  model_2.6a_waic$estimates["waic","Estimate"],
  model_2.7a_waic$estimates["waic","Estimate"]
)

comparison_tb['WAIC_probit'] = c(
  model_1.1b_waic$estimates["waic","Estimate"],
  model_1.2b_waic$estimates["waic","Estimate"],
  model_2.1b_waic$estimates["waic","Estimate"],
  model_2.2b_waic$estimates["waic","Estimate"],
  model_2.3b_waic$estimates["waic","Estimate"],
  model_2.4b_waic$estimates["waic","Estimate"],
  model_2.5b_waic$estimates["waic","Estimate"],
  model_2.6b_waic$estimates["waic","Estimate"],
  model_2.7b_waic$estimates["waic","Estimate"]
)

comparison_tb['PSIS_logit'] = c(
  model_1.1a_psis$estimates["looic","Estimate"],
  model_1.2a_psis$estimates["looic","Estimate"],
  model_2.1a_psis$estimates["looic","Estimate"],
  model_2.4a_psis$estimates["looic","Estimate"],
  model_2.3a_psis$estimates["looic","Estimate"],
  model_2.4a_psis$estimates["looic","Estimate"],
  model_2.5a_psis$estimates["looic","Estimate"],
  model_2.6a_psis$estimates["looic","Estimate"],
  model_2.7a_psis$estimates["looic","Estimate"]
)

comparison_tb['PSIS_probit'] = c(
  model_1.1b_psis$estimates["looic","Estimate"],
  model_1.2b_psis$estimates["looic","Estimate"],
  model_2.1b_psis$estimates["looic","Estimate"],
  model_2.2b_psis$estimates["looic","Estimate"],
  model_2.3b_psis$estimates["looic","Estimate"],
  model_2.4b_psis$estimates["looic","Estimate"],
  model_2.5b_psis$estimates["looic","Estimate"],
  model_2.6b_psis$estimates["looic","Estimate"],
  model_2.7b_psis$estimates["looic","Estimate"]
)

comparison_tb['Weight'] = weight_calc( comparison_tb$WAIC_probit )
comparison_tb <- arrange(comparison_tb, by=WAIC_probit)

write_csv(comparison_tb, "./data_store/fragility_models/model_comparison_2_x.csv")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Make plots ----

source("./plot_curves_1.1_1.2_2.1_2.2.R")
source("./plot_curves_2.2_2.3_2.4.R")
source("./plot_curves_2.5_2.6_2.7.R")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Save Environment ----

save.image("./data_store/frag_curves_pga_DS1_eq_eq_totVar_20240903.RData")

