# K. Trevlopoulos
# May 2025

# This script plots fragility models calculated with the script
# ./frag_curves_hen_DS0_totVar_20240903.R




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.1a ----

# Add dummy lines for Model 2.1a for all initial damage states
# Convert DS to numeric to be able to filter next.
# Add 1 to keep the same numbers
fragility_model_2.1a_plot <- mutate(fragility_model_2.1a, DS=as.numeric(DS)+1)
fragility_model_2.1a_plot <- bind_rows(
  mutate(fragility_model_2.1a_plot, DS_part_0 = 1),
  filter(fragility_model_2.1a_plot, DS > 2) %>% mutate(DS_part_0 = 2),
  filter(fragility_model_2.1a_plot, DS > 3) %>% mutate(DS_part_0 = 3),
  filter(fragility_model_2.1a_plot, DS > 4) %>% mutate(DS_part_0 = 4),
) %>% mutate(DS = factor(DS))

# fragility_model_2.1a_plot <- mutate(fragility_model_2.1a_plot, DS=factor(DS))




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Fragility Model 2.2a ----

# This Model has one surface for each damage state.
# Fix one predictor for the upper and lower bound.

# The number of cutpoints and damage states in the model
n_ds <- length(grep("^cutpoints\\[", model_2.1a_par$variable)) + 1

model_1_medians <- exp(
  model_1.1gq_precis[2:5,"mean"]/model_1.1gq_precis[1,"mean"] )

# thr_lo <- c(1e-9, model_1_medians[1:n_ds-1])
thr_lo <- c(5e-2, model_1_medians[1:n_ds-1])
thr_up <- model_1_medians





# Posterior means of the parameters

# # Use the posterior draws, or use precis for rethinking models, e.g.:
# b_DS_mean <- precis_model_2.4['b_DS','mean']

b_lnIM_part_0 <- model_2.2a_par$mean[1]
b_lnIM_part_2 <- model_2.2a_par$mean[2]
cutpoints <- model_2.2a_par$mean[3:6]




# Define IM range for plotting
IM_seq <- c(1e-9,
            seq(from=0.002, to=0.048, by=0.002),
            seq(from=0.05, to=1.0, by=0.02))
lnIM_seq <- log(IM_seq)


# Damage state categories (from cutpoints)
n_ds <- length(cutpoints) + 1




# Model 2.2a lo

fixed_im <- log(thr_lo)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b_lnIM_part_0 +
                 lnIM_seq * b_lnIM_part_2 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b_lnIM_part_0 +
                 lnIM_seq * b_lnIM_part_2 - cutpoints[k]
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

fragility_model_2.2a_lo <- fragility_curves |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 2.2a lo")




# Model 2.2a up

fixed_im <- log(thr_up)

# For each DS_part_0 (j = DS_part_0 + 1)
for (j in 1:(n_ds-1)) {
  
  if (j == 1) {
    
    # Compute fragility curves
    fragility_list <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b_lnIM_part_0 +
                 lnIM_seq * b_lnIM_part_2 - cutpoints[k]
             ), IM = IM_seq)
    })
    
    fragility_curves <- bind_rows(fragility_list)
    fragility_curves <- mutate(fragility_curves, DS_part_0=j)
    
  } else {
    
    # Compute fragility curves
    temp <- purrr::map(seq_along(cutpoints), function(k) {
      tibble(DS = k+1,
             Probability = plogis(
               fixed_im[j] * b_lnIM_part_0 +
                 lnIM_seq * b_lnIM_part_2 - cutpoints[k]
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

fragility_model_2.2a_up <- fragility_curves |> mutate(DS=factor(DS)) |>
  mutate(Model = "Model 2.2a up")




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Plot ----

# The state-dependent fragility curves (and some others)
sdfc <- bind_rows(
  fragility_model_2.1a_plot,
  fragility_model_2.2a_lo,
  fragility_model_2.2a_up,
  fragility_model_2.3a,
  fragility_model_2.4a
)

# sdfc <- bind_rows(
#   fragility_model_2.1a_plot,
#   fragility_model_2.3a,
#   fragility_model_2.4a
# )
# 
# fragility_model_2.2a <- bind_rows(
#   fragility_model_2.2a_lo,
#   fragility_model_2.2a_up,
# 
# )


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
      setNames( c("solid", "dashed", "dashed", "solid", "solid" ),
                unique(sdfc$Model) ) )+
  facet_grid(
    cols = vars(DS),
    rows = vars(DS_part_0),
    labeller = labeller(
      DS = function(x) paste0("j = ", x),
      DS_part_0 = function(x) paste0("DSpt1 = ",x)
    )
  )
plot(figure_compare)

ggsave(
  plot = figure_compare,
  paste0("./data_store/fragility_models/compare_models_2.1a_2.2a_2.3a_2.4a.png"),
  width = 29.7,
  height = 21.0,
  units = c("cm"),
  dpi = 300,
)
