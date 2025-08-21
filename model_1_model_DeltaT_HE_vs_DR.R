# K. Trevlopoulos
# March 2025

# This script does the following:
#
# 1. Calculates fragility curves for the undamaged building model as a function
# of the lnIM_part_1 (lnPGA - typical fragility curves)
#
# 2. calculates fragility curves for the undamaged building model as a function
# of the period elongation (DeltaT).
#
# 3. Calculates the values of the values of the area of the hysteretic loops
# That corresponds to the thresholds of the interstorey drift.
#
# The calculated values are intended to make fragility curves out of fragility
# surfaces by fixing one of the two predictors.




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Data for part 1 (part_0) ----
# These data will be used for the model for DS as a function of DeltaT

print("Processing part_0 data")

# Shift values one row up to use the period at step 1 (numbering 0-3)
part0_dat <- mutate(orig_dataset, per_of_max_amp = lead(per_of_max_amp))

part0_dat <- filter(part0_dat,
                    part_number==0, load_type == 'strong_gm')

part0_dat <- mutate(
  part0_dat,
  DS_part_0 = case_when(
    max_int_storey_drift < thresholds[1] ~ 1,
    max_int_storey_drift >= thresholds[1] & max_int_storey_drift < thresholds[2] ~ 2,
    max_int_storey_drift >= thresholds[2] & max_int_storey_drift < thresholds[3] ~ 3,
    max_int_storey_drift >= thresholds[3] & max_int_storey_drift < thresholds[4] ~ 4,
    max_int_storey_drift >= thresholds[4] ~ 5,
  )
)


# Remove duplicates. The same input ground motions may have been used in
# multiple simulations
# part0_dat <- select(part0_dat, file_name, gm_component, max_int_storey_drift, PGA_g)
part0_dat <- select(part0_dat, n_of_gm, gm_component,
                    PGA_g, DS_part_0, per_of_max_amp,
                    max_int_storey_drift, hyst_loop_area)
part0_dat <- unique(part0_dat)


# Keep the EDPs for the responses for ONE component from each ground motion
# The numbers of the ground motions in the simulations
part0_dat <- arrange(part0_dat, n_of_gm, gm_component)

# Before filtering out
gm_and_components_unfiltered <- select(part0_dat, n_of_gm, gm_component)

n_of_gms_sim <- select(part0_dat, n_of_gm)
n_of_gms_sim <- unique(n_of_gms_sim)

part0_dat_unfiltered <- part0_dat

for ( j in 1:nrow(n_of_gms_sim) ) {
  # Data for a single record
  temp <- filter(part0_dat, n_of_gm == n_of_gms_sim$n_of_gm[j])
  # If more than one components were used one or more times, keep the first
  # to be found
  if (nrow(temp) > 1) {
    component1 <- temp$gm_component[1]
    temp <- filter(temp, gm_component == component1)
    part0_dat <- filter(part0_dat, n_of_gm != temp$n_of_gm[1] )
    part0_dat <- bind_rows(part0_dat, temp)
  }
}


# Transform parameters
# Avoid calculating the logarithm of a period elongation equal to zero
part0_dat <- filter(part0_dat, per_of_max_amp>T_0)
part0_dat <- mutate( part0_dat,
                     lnDeltaT_part_1 = log((per_of_max_amp-T_0)/T_0),
                     lnIM_part_0 = log(PGA_g))


# For testing
# The gms and components used
gm_and_components <- select(part0_dat, n_of_gm, gm_component)
gm_and_components <- arrange(gm_and_components, n_of_gm, gm_component)

# n_of_gm, gm_component are no longer needed
part0_dat <- select(part0_dat, -n_of_gm, -gm_component)
part0_dat_tib <- part0_dat
part0_dat <- list(
  lnDeltaT_part_1 = part0_dat_tib$lnDeltaT_part_1,
  DS_part_0 = as.integer(part0_dat_tib$DS_part_0),
  N = nrow(part0_dat_tib)
)

# # Export if needed
# write_csv(part0_dat_tib, './part0_dat.csv')


part0_dat_tib <- mutate(
  part0_dat_tib,
  DS_part_0_text = factor(
    case_match(
      DS_part_0, 
      1 ~ "DS0", 2 ~ "DS1", 3 ~ "DS2", 4 ~ "DS3", 5 ~ "DS4",)))




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Boxplot for DeltaT ----

# The title of the y-axis
y_axis_t <- paste0('$\\Delta T$')

boxplot_DeltaT <- ggplot(
  mutate(part0_dat_tib, DS_part_0=factor(DS_part_0)),
  aes(x = DS_part_0,
      y = exp(lnDeltaT_part_1),
      fill = DS_part_0)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_viridis_d( name = "Damage State", end = 0.9 ) +
  xlab('Damage State') +
  ylab( TeX( y_axis_t ) ) +
  theme_minimal()
plot(boxplot_DeltaT)

ggsave(
  plot = boxplot_DeltaT,
  paste0("./data_store/fragility_models/boxplot_DSpart0_DeltaT.png"),
  width = 5.0*2.54,
  height = 3.0*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model DeltaT ----

print("model_DS_DeltaT")

model_DS_DeltaT <- ulam(
  alist(
    DS_part_0 ~ ordered_logistic(eta, cutpoints),
    eta <- b0 * lnDeltaT_part_1,
    
    b0 ~ normal(0, 1),
    cutpoints ~ normal(0, 1.5),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS_part_0[i] |
        b0 * lnDeltaT_part_1[i],
      cutpoints
    )
  ),
  data = part0_dat,
  chains = 4,
  cores = 4,
  iter = 2000
)

print("Calculation complete.")

# # Extract Stan code
# stan_model_DS_DeltaT <- stancode(model_DS_DeltaT)
# # Write Stan code to file. Comment if already exported and edited
# writeLines(stan_model_DS_DeltaT, "model_DS_DeltaT.stan")


### Checking chains ----

png(
  filename = "./data_store/checks/model_DS_DeltaT_traceplot.png",
  width = 9.0 * 2.54,
  height = 3.0 * 2.54,
  units = "cm",
  res = 300
)
traceplot( model_DS_DeltaT )
dev.off()

png(
  filename = "./data_store/checks/model_DS_DeltaT_trankplot.png",
  width = 9.0 * 2.54,    # convert cm to inches if needed
  height = 3.0 * 2.54,
  units = "cm",
  res = 300
)
trankplot( model_DS_DeltaT )
dev.off()


# The parameters
model_DS_DeltaT_precis <- precis(model_DS_DeltaT, depth=2)
# print(model_DS_DeltaT_precis)
#               mean   sd  5.5% 94.5% rhat ess_bulk
# b0            2.75 0.22  2.42  3.11    1  1538.89
# cutpoints[1] -4.61 0.41 -5.28 -3.99    1  1499.89
# cutpoints[2] -3.00 0.32 -3.52 -2.49    1  1876.80
# cutpoints[3] -1.82 0.28 -2.28 -1.39    1  2615.88
# cutpoints[4] -0.83 0.27 -1.27 -0.40    1  3351.69



### Calculate fragility curves ----

# The names of the rows for the cutpoints in precis_model_DS_Tnorm
crow_names <- paste0("cutpoints[", 1:4, "]")
cut_DeltaT <- model_DS_DeltaT_precis[crow_names,'mean']
b_DeltaT <- model_DS_DeltaT_precis['b0','mean']

n_ds <- length(grep("^cutpoints\\[", rownames(model_DS_DeltaT_precis) )) + 1

# The DeltaT thresholds
lnDeltaT_thres <- cut_DeltaT/b_DeltaT
DT_thres <- exp(lnDeltaT_thres)

# save(DT_thres, file="./DT_thres.RData")




# The IM range
DeltaT_seq <- seq(from=0.02, to=2.0, by=0.02)
DeltaT_seq <- c(1e-9, DeltaT_seq)
lnDeltaT_seq <- log(DeltaT_seq)

fragility_DeltaT <- purrr::map(
  seq_along(cut_DeltaT),
  function(k) {
    tibble(DS = k+1,
           Probability = plogis(
             lnDeltaT_seq * b_DeltaT - cut_DeltaT[k]
           ),
           IM = DeltaT_seq)
  }
)
fragility_DeltaT <- bind_rows(fragility_DeltaT)
fragility_DeltaT <- mutate(fragility_DeltaT, DS=factor(DS))

# The titles of the axes
x_axis_t <- paste0('$\\Delta T$')
y_axis_t <- paste0('$P( DS \\geq j )$')

# Plot fragility curves
figure_DeltaT <- ggplot(fragility_DeltaT,
                        aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  xlab( TeX(x_axis_t) )+
  ylab( TeX(y_axis_t) )+
  theme_minimal()+
  scale_color_viridis_d(name = "j", end = 0.9)
plot(figure_DeltaT)

ggsave(
  plot = figure_DeltaT,
  paste0("./data_store/fragility_models/model_DeltaT.png"),
  width = 5.0*2.54,
  height = 3.0*2.54,
  units = c("cm"),
  dpi = 300,
)




## Probabilities of the Damage States ----

# Pivot to calculate the differences between fragility curves
prob_DS_DeltaT <- pivot_wider(fragility_DeltaT,
                              names_from = DS, values_from = Probability)

# Rename columns (probex = probability of exceedance of the DS threshold) 
# prob_DS_DeltaT <- prob_DS_DeltaT %>%
#   rename_with( ~ paste0("probex_DS", seq_along(.)), starts_with("DS ") )
for (j in 2:n_ds) {
  prob_DS_DeltaT <- rename(prob_DS_DeltaT,
                           !!paste0("probex_DS", j) := all_of(j))
}


# Names of the columns for the differences
prob_cols <- prob_DS_DeltaT %>% select(starts_with("probex_DS")) %>% names()
# The differences
prob_DS_DeltaT <- map2_dfc(prob_cols[-length(prob_cols)], prob_cols[-1], ~ {
  new_col <- paste0("diff_", .x, "_", .y)
  tibble(!!new_col := prob_DS_DeltaT[[.x]] - prob_DS_DeltaT[[.y]])
}) %>%
  bind_cols(prob_DS_DeltaT, .)

# Rename the columns. The differences are the probabilities of the DS
prob_DS_DeltaT <- prob_DS_DeltaT %>%
  rename_with(.fn = ~ paste0("prob_DS", seq_along(.x)+1),
              .cols = starts_with("diff_probex_DS"))

# The probabilities of DS 1 and 5
prob_DS_DeltaT <- rename(prob_DS_DeltaT, prob_DS5 = probex_DS5) %>%
  mutate(prob_DS1 = 1-probex_DS2)


# Select the IM (column 1) and the probabilities of the DS
prob_DS_DeltaT <- select(
  prob_DS_DeltaT,
  c(1, (ncol(prob_DS_DeltaT)-n_ds+1) : ncol(prob_DS_DeltaT)) )

# Reorder the columns so that the renaming is correct
prob_DS_DeltaT <- select(prob_DS_DeltaT, order(colnames(prob_DS_DeltaT)))

# Rename the columns. Subtract 1 for naming DS0 to DS4
prob_DS_DeltaT <- prob_DS_DeltaT %>%
  rename_with(.fn = ~ paste0(seq_along(.x)),
              .cols = starts_with("prob_DS"))

prob_DS_DeltaT <- pivot_longer(prob_DS_DeltaT, cols = 2:ncol(prob_DS_DeltaT),
                               names_to = "DS", values_to = "Probability")

prob_DS_DeltaT <- mutate(prob_DS_DeltaT, DS=factor(DS))


# The titles of the axes
x_axis_t <- paste0('$\\Delta T$')
y_axis_t <- paste0('$P( DS=j )$')

# Plot the probability of the DS based on DeltaT
figure_P_DS_DeltaT <- ggplot(
  prob_DS_DeltaT,
  aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  xlab( TeX(x_axis_t) )+
  ylab( TeX(y_axis_t) )+
  theme_minimal()+
  scale_color_viridis_d(name = "j", end = 0.9)
plot(figure_P_DS_DeltaT)

ggsave(
  plot = figure_P_DS_DeltaT,
  paste0("./data_store/fragility_models/model_DeltaT_prob_DS.png"),
  width = 5.0*2.54,
  height = 3.0*2.54,
  units = c("cm"),
  dpi = 300,
)




# Hyst. loop area vs Drift ----

## Part 1 ----

# The title of the y-axis
x_axis_t <- paste0('Hyst. loop area')
y_axis_t <- paste0('Max. inter-storey drift')
figure_title <- paste0('Part 1 - Unique cases')

fig_hen_drift1 <- ggplot(
  # filter(part0_dat_tib, DS_part_0>1),
  part0_dat_tib,
  aes( x=hyst_loop_area, y=max_int_storey_drift ) )+ 
  geom_point()+
  scale_x_log10(
    breaks = 10^seq(-3,2),
    minor_breaks = unlist(lapply(seq(-3,2), function(e) (2:9) * 10^e)),
    limits = c(1e-3,2e2),
  ) +
  scale_y_log10(
    minor_breaks = unlist(lapply(seq(-4,-2), function(e) (2:9) * 10^e)),
    limits = c(1e-4,1e-1),
  ) +
  ggtitle( TeX( figure_title ) )+
  xlab( x_axis_t )+
  ylab( y_axis_t )+
  theme_light()+                                     
  theme(plot.title = element_text(size = 11))
print(fig_hen_drift1)


# Add to the plot lines for thresholds
threshold_plots <- tibble(drift = thresholds)
threshold_plots <- threshold_plots |> mutate(DS = factor(seq(2,5)))
threshold_plots <- threshold_plots |> mutate(pga = 1e-3)
threshold_plots <- threshold_plots |> mutate(pga = 2e2) |> bind_rows(threshold_plots)

fig_hen_drift1b <- fig_hen_drift1+
  geom_path(data=threshold_plots,
            aes(x=pga, y=drift, group = DS, color = DS),
            linewidth = 1.0, linetype = 2)+
  scale_color_viridis_d(name = "Damage State", begin=0.4, end=0.9)+
  theme(legend.position="bottom")
plot(fig_hen_drift1b)


ggsave(
  plot = fig_hen_drift1b,
  paste0("./data_store/fragility_models/drift_vs_hen_DS0_20240903", ".png"),
  width = 5.0*2.54,
  height = 4.5*2.54,
  units = c("cm"),
  dpi = 300,
)




## Hyst. loop area thresholds ----

# Add the log of the area of hysteretic loops and the log of the drift
part0_dat_tib <- mutate(part0_dat_tib, lnHE = log(hyst_loop_area),
                        lnDR = log(max_int_storey_drift))

# Linear regression
regr_drift_hen <- lm(lnHE ~ lnDR,
                     data = filter(part0_dat_tib, DS_part_0>1))

# The values of the area of the hysteretic loops that correspond to the
# thresholds of the max. interstorey drift.
hen_thresholds <- exp( regr_drift_hen$coefficients[1] +
                         regr_drift_hen$coefficients[2]*log(thresholds) )

# # Write thresholds to disk
# save(hen_thresholds, file = "./hen_thresholds.RData")




## Part 4 ----

part2_dat_tib <- orig_dataset

part2_dat_tib <- mutate(
  part2_dat_tib,
  DS = case_when(
    max_int_storey_drift < thresholds[1] ~ 0,
    max_int_storey_drift >= thresholds[1] & max_int_storey_drift < thresholds[2] ~ 1,
    max_int_storey_drift >= thresholds[2] & max_int_storey_drift < thresholds[3] ~ 2,
    max_int_storey_drift >= thresholds[3] & max_int_storey_drift < thresholds[4] ~ 3,
    max_int_storey_drift >= thresholds[4] ~ 4,
  )
)

part2_dat_tib <- part2_dat_tib |> rename(DS_part_0 = DS) |>
  mutate(DS_part_2 = lag(DS_part_0, n=2)) |> 
  filter(part_number==2) |> mutate(DS_part_2 = max(DS_part_0, DS_part_2)) 




# The title of the y-axis
x_axis_t <- paste0('Hyst. loop area')
y_axis_t <- paste0('Max. inter-storey drift')
figure_title <- paste0('Part 3')

fig_hen_drift2 <- ggplot(
  # filter(part2_dat_tib, DS_part_0>1),
  part2_dat_tib,
  aes( x=hyst_loop_area, y=max_int_storey_drift ) )+ 
  geom_point()+
  scale_x_log10(
    breaks = 10^seq(-3,2),
    minor_breaks = unlist(lapply(seq(-3,2), function(e) (2:9) * 10^e)),
    limits = c(1e-3,2e2),
  ) +
  scale_y_log10(
    minor_breaks = unlist(lapply(seq(-4,-2), function(e) (2:9) * 10^e)),
    limits = c(1e-4,1e-1),
  ) +
  ggtitle( TeX( figure_title ) )+
  xlab( x_axis_t )+
  ylab( y_axis_t )+
  theme_light()+                                     
  theme(plot.title = element_text(size = 11))
print(fig_hen_drift2)


fig_hen_drift2b <- fig_hen_drift2+
  geom_path(data=threshold_plots,
            aes(x=pga, y=drift, group = DS, color = DS),
            linewidth = 1.0, linetype = 2)+
  scale_color_viridis_d(name = "Damage State", begin=0.4, end=0.9)+
  theme(legend.position="bottom")
plot(fig_hen_drift2b)


ggsave(
  plot = fig_hen_drift2b,
  paste0("./data_store/fragility_models/drift_vs_hen_part3_20240903", ".png"),
  width = 5.0*2.54,
  height = 4.5*2.54,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Data for Models 1.x ----

model_1_data <- list(
  lnIM = log(part0_dat_tib$PGA_g),
  DS = as.integer(part0_dat_tib$DS_part_0),
  N = nrow(part0_dat_tib)
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.1 ----

print("Model 1.1")

model_1.1gq <- ulam(
  alist(
    DS ~ ordered_logistic(eta, cutpoints),
    
    # Linear model for eta
    eta <- b0 * lnIM,
    
    # Priors
    b0 ~ normal(4, 2),
    cutpoints ~ normal(-4, 2),
    
    # Log-likelihood for model comparison
    gq> vector[N]: log_lik <- ordered_logistic_lpmf(
      DS[i] | b0 * lnIM[i], cutpoints
    )
  ),
  data = model_1_data,
  chains = 4,
  cores = 4,
  iter = 2000
)

# # Extract and write Stan code to file. Comment if already exported and edited
# writeLines( stancode(model_1.1gq) , "model_1_1gq.stan")




## Checking chains ----

png(
  filename = "./data_store/checks/model_1_1gq_traceplot.png",
  width = 21.0,
  height = 29.7*0.33,
  units = "cm",
  res = 300
)
traceplot( model_1.1gq )
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

# Export the parameters
write_csv(
  model_1.1gq_precis |>
    as.matrix() |>
    as_tibble(rownames = "Parameter"),
  "./data_store/fragility_models/model_1_1gq_precis.csv")










## Compute fragility curves ----

# The number of cutpoints and damage states in the model
n_ds <- length(grep("^cutpoints\\[", rownames(model_1.1gq_precis)))

# Mean parameters
b0 <- model_1.1gq_precis["b0", "mean"]
crow_names <- paste0("cutpoints[", 1:n_ds, "]")
cutpoints <- model_1.1gq_precis[crow_names, "mean"]

# Define IM range for plotting
IM_seq <- seq(from=0.02, to=1.0, by=0.02)
IM_seq <- c(1e-9, IM_seq)
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1


fragility_curves <- purrr::map(seq_along(cutpoints), function(k) {
  tibble(DS = k+1,
         Probability = plogis(
           lnIM_seq * b0 - cutpoints[k]
         ), IM = IM_seq)
})
fragility_curves <- bind_rows(fragility_curves)

fragility_curves['DS'] <- factor(fragility_curves$DS)

fragility_model_1a <- fragility_curves
fragility_model_1a["Model"] = "Model 1"




figure_m1 <- ggplot(fragility_model_1a, aes(x = IM, y = Probability, color = DS)) +
  geom_line(linewidth = 1) +
  scale_color_viridis_d(name = "Damage State j", end=0.9) +
  labs(
    title = "Part 1",
    x = "PGA (g)",
    y = TeX("$P(DS \\geq j|DS=1)$")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
plot(figure_m1)

ggsave(
  plot = figure_m1,
  paste0("./data_store/fragility_models/model_1_1.png"),
  width = 5.0*2.54,
  height = 5.0*2.54*0.75,
  units = c("cm"),
  dpi = 300,
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.1a (stan; logit) ----

# quantile( exp(rnorm(n=1000,mean=1,sd=1.0)), probs = c(0.055,0.5,0.955) )
# 1/quantile( exp(rnorm(n=1000,mean=1,sd=1.0)), probs = c(0.055,0.5,0.955) )

print("Model 1.1a")

# Compile Stan model with cmdstanr
model_1.1a <- cmdstan_model("model_1_1gq.stan")

# Fit the model
model_1.1a_fit <- model_1.1a$sample(
  data = model_1_data,
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
model_1.1a_waic <- waic( model_1.1a_fit$draws("log_lik", format = "matrix") )
# print( model_1.1a_waic )
model_1.1a_psis <- loo( model_1.1a_fit$draws("log_lik", format = "matrix") )
# print( model_1.1a_psis )




## Checking chains ----

model_1_1a_traceplot <- mcmc_trace(
  model_1.1a_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  labs(title = "Model 1.1a")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1a_traceplot,
  paste0("./data_store/checks/model_1_1a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)


model_1_1a_trankplot <- mcmc_rank_overlay(
  model_1.1a_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.1a")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1a_trankplot,
  paste0("./data_store/checks/model_1_1a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.2a (stan; logit) ----

print("Model 1.2a")

# Compile Stan model with cmdstanr
model_1.2a <- cmdstan_model("model_1_2gq.stan")

# Fit the model
model_1.2a_fit <- model_1.2a$sample(
  data = model_1_data,
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
model_1.2a_post <- select(model_1.2a_post, 2:6)
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
model_1.2a_waic <- waic( model_1.2a_fit$draws("log_lik", format = "matrix") )
# print( model_1.2a_waic )
model_1.2a_psis <- loo( model_1.2a_fit$draws("log_lik", format = "matrix") )
# print( model_1.2a_psis )




## Checking chains ----

model_1_2a_traceplot <- mcmc_trace(
  model_1.2a_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  labs(title = "Model 1.2a")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2a_traceplot,
  paste0("./data_store/checks/model_1_2a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)


model_1_2a_trankplot <- mcmc_rank_overlay(
  model_1.2a_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.2a")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2a_trankplot,
  paste0("./data_store/checks/model_1_2a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.3a (stan; logit) ----

print("Model 1.3a")

# Compile Stan model with cmdstanr
model_1.3a <- cmdstan_model("model_1_3gq.stan")

# Fit the model
model_1.3a_fit <- model_1.3a$sample(
  data = model_1_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_1.3a_post <- model_1.3a_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_1.3a_post <- select(model_1.3a_post, 2:6)
# Calculate the parameters
model_1.3a_par <- model_1.3a_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_1.3a_par, "./data_store/fragility_models/model_1.3a_par.csv")


# Cross-validation and information criteria
model_1.3a_waic <- waic( model_1.3a_fit$draws("log_lik", format = "matrix") )
# print( model_1.3a_waic )
model_1.3a_psis <- loo( model_1.3a_fit$draws("log_lik", format = "matrix") )
# print( model_1.3a_psis )




## Checking chains ----

model_1_3a_traceplot <- mcmc_trace(
  model_1.3a_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  labs(title = "Model 1.3a")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_3a_traceplot,
  paste0("./data_store/checks/model_1_3a_traceplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)


model_1_3a_trankplot <- mcmc_rank_overlay(
  model_1.3a_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.3a")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_3a_trankplot,
  paste0("./data_store/checks/model_1_3a_trankplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.1b (stan; probit) ----

print("Model 1.1b")

# Compile Stan model with cmdstanr
model_1.1b <- cmdstan_model("model_1_1b.stan")

# Fit the model
model_1.1b_fit <- model_1.1b$sample(
  data = model_1_data,
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
model_1.1b_waic <- waic( model_1.1b_fit$draws("log_lik", format = "matrix") )
# print( model_1.1b_waic )
model_1.1b_psis <- loo( model_1.1b_fit$draws("log_lik", format = "matrix") )
# print( model_1.1b_psis )




## Checking chains ----

model_1_1b_traceplot <- mcmc_trace(
  model_1.1b_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  labs(title = "Model 1.1b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1b_traceplot,
  paste0("./data_store/checks/model_1_1b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)


model_1_1b_trankplot <- mcmc_rank_overlay(
  model_1.1b_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.1b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_1b_trankplot,
  paste0("./data_store/checks/model_1_1b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.2b (stan; probit) ----

print("Model 1.2b")

# Compile Stan model with cmdstanr
model_1.2b <- cmdstan_model("model_1_2b.stan")

# Fit the model
model_1.2b_fit <- model_1.2b$sample(
  data = model_1_data,
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
model_1.2b_post <- select(model_1.2b_post, 2:6)
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
model_1.2b_waic <- waic( model_1.2b_fit$draws("log_lik", format = "matrix") )
# print( model_1.2b_waic )
model_1.2b_psis <- loo( model_1.2b_fit$draws("log_lik", format = "matrix") )
# print( model_1.2b_psis )




## Checking chains ----

model_1_2b_traceplot <- mcmc_trace(
  model_1.2b_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  labs(title = "Model 1.2b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2b_traceplot,
  paste0("./data_store/checks/model_1_2b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)


model_1_2b_trankplot <- mcmc_rank_overlay(
  model_1.2b_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.2b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_2b_trankplot,
  paste0("./data_store/checks/model_1_2b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Model 1.3b (stan; probit) ----

print("Model 1.3b")

# Compile Stan model with cmdstanr
model_1.3b <- cmdstan_model("model_1_3b.stan")

# Fit the model
model_1.3b_fit <- model_1.3b$sample(
  data = model_1_data,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_sampling = 1000,
  iter_warmup = 1000
)

print("Calculation complete.")

# Extract the draws into a data frame
model_1.3b_post <- model_1.3b_fit$draws(variables = NULL,
                                        format = "data.frame")

# Select the parameters
model_1.3b_post <- select(model_1.3b_post, 2:6)
# Calculate the parameters
model_1.3b_par <- model_1.3b_post %>%
  summarise_draws(
    mean,
    sd,
    ~quantile(.x, probs = 0.055),
    ~quantile(.x, probs = 0.955)
  )

# Export the parameters
write_csv(model_1.3b_par, "./data_store/fragility_models/model_1.3b_par.csv")


# Cross-validation and information criteria
model_1.3b_waic <- waic( model_1.3b_fit$draws("log_lik", format = "matrix") )
# print( model_1.3b_waic )
model_1.3b_psis <- loo( model_1.3b_fit$draws("log_lik", format = "matrix") )
# print( model_1.3b_psis )




## Checking chains ----

model_1_3b_traceplot <- mcmc_trace(
  model_1.3b_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  labs(title = "Model 1.3b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_3b_traceplot,
  paste0("./data_store/checks/model_1_3b_traceplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)


model_1_3b_trankplot <- mcmc_rank_overlay(
  model_1.3b_fit$draws(format = "draws_array"),
  pars = c("b0",
           "cutpoints[1]", "cutpoints[2]",
           "cutpoints[3]", "cutpoints[4]"),
  facet_args = list(ncol = 3, nrow = 2))+
  ylim(20, 80)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank() )+
  labs(title = "Model 1.3b")+
  theme_minimal()+
  scale_colour_viridis_d(end = 0.9)

ggsave(
  plot = model_1_3b_trankplot,
  paste0("./data_store/checks/model_1_3b_trankplot.png"),
  width = 21.0,
  height = 29.7*0.33,
  units = c("cm"),
  dpi = 300
)




# Model comparison ----

compar_model_1 = tibble(Model = c("Model 1.1", "Model 1.2", "Model 1.3"))

# Add WAIC and PSIS to the table
compar_model_1['WAIC_logit'] = c(
  model_1.1a_waic$estimates["waic","Estimate"],
  model_1.2a_waic$estimates["waic","Estimate"],
  model_1.3a_waic$estimates["waic","Estimate"]
)

compar_model_1['WAIC_probit'] = c(
  model_1.1b_waic$estimates["waic","Estimate"],
  model_1.2b_waic$estimates["waic","Estimate"],
  model_1.3b_waic$estimates["waic","Estimate"]
)

compar_model_1['PSIS_logit'] = c(
  model_1.1a_psis$estimates["looic","Estimate"],
  model_1.2a_psis$estimates["looic","Estimate"],
  model_1.3a_psis$estimates["looic","Estimate"]
)

compar_model_1['PSIS_probit'] = c(
  model_1.1b_psis$estimates["looic","Estimate"],
  model_1.2b_psis$estimates["looic","Estimate"],
  model_1.3b_psis$estimates["looic","Estimate"]
)

compar_model_1['Weight'] = weight_calc( compar_model_1$WAIC_probit )
compar_model_1 <- arrange(compar_model_1, by=WAIC_probit)

# Export the comparison for the Models 1.x
write_csv(compar_model_1, "./data_store/fragility_models/compar_model_1.csv")


# > sqrt(0.323^2 + 0.25^2 + 0.4^2)
# [1] 0.5716896
# > sqrt(0.382^2 + 0.25^2 + 0.4^2)
# [1] 0.6069794




# Parameter table ----

# Tibbles for the model parameters and for the parameters of the
# fragility curves
maxlik_Bayes_param <- tibble(parameter = c("b0", "κ1", "κ2", "κ3", "κ4"))
maxlik_Bayes_fc <- tibble(parameter = c("σ", "μ2", "μ3", "μ4", "μ5"))


# Max. lik. calculation

model_max_lik <- clm(
  DS ~ lnIM,
  data = part0_dat_tib |>
    select(DS_part_0, lnIM_part_0) |>
    rename(DS=DS_part_0, lnIM=lnIM_part_0) |>
    mutate(DS=factor(DS)),
  link = "probit")

maxlik_Bayes_param["Max. Lik."] <- c(
  as.numeric(model_max_lik$beta), as.numeric(model_max_lik$alpha) )
maxlik_Bayes_fc["Max. Lik."] <- c(
  as.numeric(1/as.numeric(model_max_lik$beta)),
  as.numeric(exp(model_max_lik$alpha/model_max_lik$beta)) )

maxlik_Bayes_param["Model 1.1b"] <- model_1.1b_par$mean
maxlik_Bayes_fc["Model 1.1b"] <- c(
  1/model_1.1b_par$mean[1],
  exp(model_1.1b_par$mean[2:5]/model_1.1b_par$mean[1]) )

maxlik_Bayes_param["Model 1.2b"] <- model_1.2b_par$mean
maxlik_Bayes_fc["Model 1.2b"] <- c(
  1/model_1.2b_par$mean[1],
  exp(model_1.2b_par$mean[2:5]/model_1.2b_par$mean[1]) )

maxlik_Bayes_param["Model 1.3b"] <- model_1.3b_par$mean
maxlik_Bayes_fc["Model 1.3b"] <- c(
  1/model_1.3b_par$mean[1],
  exp(model_1.3b_par$mean[2:5]/model_1.3b_par$mean[1]) )

# Export the model parameters and for the parameters of the fragility curves
write_csv(maxlik_Bayes_param, "./data_store/fragility_models/maxlik_Bayes_param.csv")
write_csv(maxlik_Bayes_fc, "./data_store/fragility_models/maxlik_Bayes_fc.csv")


# Compare max. lik. & Bayesian ----

# Max. lik. calculation

# model_max_lik <- clm(
#   DS ~ lnIM,
#   data = part0_dat_tib |>
#     select(DS_part_0, lnIM_part_0) |>
#     rename(DS=DS_part_0, lnIM=lnIM_part_0) |>
#     mutate(DS=factor(DS)),
#   link = "logit")
# 
# print('Maximum Likelihood calculation (logit): ')
# print(summary(model_max_lik))
# 
# print('Bayesian model (logit): ')
# print(model_1gq_precis)