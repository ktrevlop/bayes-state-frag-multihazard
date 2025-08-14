# K. Trevlopoulos
# May 2025

# This script plots fragility models whose parameters are calculated by
# ./frag_models_pga_eq_eq_20240903.R




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.5a ----

# This Model has one surface for each damage state.
# Fix one predictor for the upper and lower bound.

# The number of cutpoints and damage states in the model
n_ds <- length(grep("^cutpoints\\[", model_2.5a_par$variable)) + 1

# thr_lo <- c(1e-9, thresholds[1:n_ds-1])
thr_lo <- c(1e-4, thresholds[1:n_ds-1])
thr_up <- thresholds





# Posterior means of the parameters

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

b0 <- model_2.5a_par$mean[1]
b1 <- model_2.5a_par$mean[2]
cutpoints <- model_2.5a_par$mean[3:6]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)



# Model x lo

fixed_im <- log(thr_lo)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
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

fragility_model_2.5a_lo <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.5a lo")




# Model x up

fixed_im <- log(thr_up)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
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

fragility_model_2.5a_up <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.5a up")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.6a ----

load("./hen_thresholds.RData")

# This Model has one surface for each damage state.
# Fix one predictor for the upper and lower bound.

# The number of cutpoints and damage states in the model
n_ds <- length(grep("^cutpoints\\[", model_2.6a_par$variable)) + 1

# thr_lo <- c(1e-9, hen_thresholds[1:n_ds-1])
thr_lo <- c(5e-3, hen_thresholds[1:n_ds-1])
thr_up <- hen_thresholds




# Posterior means of the parameters

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

b0 <- model_2.6a_par$mean[1]
b1 <- model_2.6a_par$mean[2]
cutpoints <- model_2.6a_par$mean[3:6]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)




# Model x lo

fixed_im <- log(thr_lo)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
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

fragility_model_2.6a_lo <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.6a lo")




# Model x up

fixed_im <- log(thr_up)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
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

fragility_model_2.6a_up <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.6a up")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.7a ----

# This Model has one surface for each damage state.
# Fix one predictor for the upper and lower bound.

# The number of cutpoints and damage states in the model
n_ds <- length(grep("^cutpoints\\[", model_2.7a_par$variable)) + 1

thr_lo <- c(1e-9, DT_thres[1:n_ds-1])
thr_up <- DT_thres




# Posterior means of the parameters

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

b0 <- model_2.7a_par$mean[1]
b1 <- model_2.7a_par$mean[2]
cutpoints <- model_2.7a_par$mean[3:6]




# Define IM range for plotting
IM_seq <- seq(from=0.02, to=1.0, by=0.02)
IM_seq <- c(1e-9, IM_seq)
lnIM_seq <- log(IM_seq)




# Model x lo

# fixed_im <- log(thr_lo)
fixed_im <- thr_lo # do not take the log, see Model 2.4.2

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
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

fragility_model_2.7a_lo <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.7a lo")



# Model x up

# fixed_im <- log(thr_up)
fixed_im <- thr_up # do not take the log, see Model 2.4.2

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b1 +
                 lnIM_seq * b0 - cutpoints[k]
             ), IM = IM_seq)
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

fragility_model_2.7a_up <- fragility_curves |> mutate(DS = factor(DS)) |> 
  mutate(Model = "Model 2.7a up")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot ----

# The state-dependent fragility curves (and some others)
sdfc <- bind_rows(
  fragility_model_2.5a_lo,
  fragility_model_2.5a_up,
  fragility_model_2.6a_lo,
  fragility_model_2.6a_up,
  fragility_model_2.7a_lo,
  fragility_model_2.7a_up
)




# Plot the probabilities of exceeding the thresholds of the damage states
figure_compare <- ggplot(
  sdfc,
  aes( x=IM, y=Probability, color=Model, linetype=Model) )+
  geom_line(linewidth = 1)+
  xlab( 'PGA (g)' )+
  ylab( TeX( '$P(DSpt3 \\geq j | DSpt1)$' ) )+
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
      setNames( c("solid", "solid", "dotted", "dotted", "dashed", "dashed"),
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
  paste0("./compare_models_2.5a_2.6a_2.7a.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)
