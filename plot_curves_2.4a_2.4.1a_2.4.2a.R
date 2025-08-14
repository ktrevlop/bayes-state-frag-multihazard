# K. Trevlopoulos
# May 2025

# This script plots fragility models whose parameters are calculated by
# ./frag_models_pga_eq_eq_20240903.R




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.4.1a ----

# This Model has one surface for each damage state.
# Fix one predictor for the upper and lower bound.

# thr_lo <- c(1e-9, hen_thresholds[1:length(hen_thresholds)-1])
thr_lo <- c(5e-3, hen_thresholds[1:length(hen_thresholds)-1])
thr_up <- hen_thresholds




# Posterior means of the parameters

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

b0 <- model_2.4.1a_par$mean[1:5]
b1 <- model_2.4.1a_par$mean[6]
b2 <- model_2.4.1a_par$mean[7:11]
cutpoints <- model_2.4.1a_par$mean[12:15]
delta <- model_2.4.1a_par$mean[16:19]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1




# Model 2.4.1a lo

fixed_im <- log(thr_lo)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b2[j] +
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
               fixed_im[j] * b2[j] +
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

fragility_model_2.4.1a_lo <- fragility_curves |> mutate(DS = factor(DS)) |>
  mutate(Model = "Model 2.4.1a lo")




# Model 2.4.1a up

fixed_im <- log(thr_up)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b2[j] +
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
               fixed_im[j] * b2[j] +
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

fragility_model_2.4.1a_up <- fragility_curves |> mutate(DS = factor(DS)) |>
  mutate(Model = "Model 2.4.1a up")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.4.2a ----

# This Model has one surface for each damage state.
# Fix one predictor for the upper and lower bound.

# The number of cutpoints and damage states in the model
n_ds <- length(grep("^cutpoints\\[", rownames(model_DS_DeltaT_precis) )) + 1

thr_lo <- c(1e-9, DT_thres[1:length(DT_thres)-1])
thr_up <- DT_thres




# Posterior means of the parameters

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

b0 <- model_2.4.2a_par$mean[1:5]
b1 <- model_2.4.2a_par$mean[6]
b2 <- model_2.4.2a_par$mean[7:11]
cutpoints <- model_2.4.2a_par$mean[12:15]
delta <- model_2.4.2a_par$mean[16:19]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)

# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1




# Model 2.4.2a lo

# fixed_im <- log(thr_lo)
fixed_im <- thr_lo # do not take the log, see Model 2.4.2

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b2[j] +
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
               fixed_im[j] * b2[j] +
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

fragility_model_2.4.2a_lo <- fragility_curves |> mutate(DS = factor(DS)) |>
  mutate(Model = "Model 2.4.2a lo")




# Model 2.4.2a up

# fixed_im <- log(thr_up)
fixed_im <- thr_up # do not take the log, see Model 2.4.2

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b2[j] +
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
               fixed_im[j] * b2[j] +
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

fragility_model_2.4.2a_up <- fragility_curves |> mutate(DS = factor(DS)) |>
  mutate(Model = "Model 2.4.2a up")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot ----

# The state-dependent fragility curves (and some others)
sdfc <- bind_rows(
  fragility_model_2.4.1a_lo,
  fragility_model_2.4.1a_up,
  fragility_model_2.4.2a_lo,
  fragility_model_2.4.2a_up,
  fragility_model_2.4a
)




# The title of the y-axis
y_axis_t <- paste0('$P(DSpt3 \\geq j | DSpt1)$')

# Plot the probabilities of exceeding the thresholds of the damage states
figure_compare <- ggplot(
  sdfc,
  aes( x=IM, y=Probability, color=Model, linetype = Model) )+
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
  scale_linetype_manual(
    values =
      setNames( c("dotted", "dotted", "dashed", "dashed", "solid"),
                unique(sdfc$Model) ) )+
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
  paste0("./compare_models_2.4.1a_2.4.2a_2.3a_2.4a.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)
