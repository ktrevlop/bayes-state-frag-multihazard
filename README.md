# bayes-state-frag-multihazard
Bayesian state-dependent fragility modeling for multi-hazard sequences (earthquake/flood)

This repository contains the code and the data for our paper.

To calculate the 1.x and 2.x models execute `frag_models_pga_eq_eq_20240903.R`. Once it has been executed, launch a new session and execute the script `frag_models_pga_haz_eq_20240903.R` to calculate the 3.x models.
Use the document outline in RStudio to navigate these scripts, because they contain a large number of code lines. They may be split to smaller files in future versions.

The scripts `plot_curves_1.1_1.2_2.1_2.2.R`, `plot_curves_2.2_2.3_2.4.R`, `plot_curves_2.5_2.6_2.7.R` create the corresponding figures in the paper.

The folders `fragility_models` and `checsks` in `data_store` contain files written by the execution of the code, which include the results in our paper, other results omitted from the paper, because it was already too long.

Information about the R session that may help reproduce the results:
```
> sessionInfo()
R version 4.4.3 (2025-02-28 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 10 x64 (build 19045)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8 
[2] LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8
[4] LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: Europe/Paris
tzcode source: internal

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] loo_2.8.0           bayesplot_1.13.0    ggpubr_0.6.1       
 [4] latex2exp_0.9.6     here_1.0.1          ordinal_2023.12-4.1
 [7] lubridate_1.9.4     forcats_1.0.0       stringr_1.5.1      
[10] dplyr_1.1.4         purrr_1.0.4         readr_2.1.5        
[13] tidyr_1.3.1         tibble_3.3.0        ggplot2_3.5.2      
[16] tidyverse_2.0.0     rethinking_2.42     posterior_1.6.1    
[19] cmdstanr_0.9.0     

loaded via a namespace (and not attached):
 [1] gtable_0.3.6         shape_1.4.6.1        tensorA_0.36.2.1    
 [4] xfun_0.52            processx_3.8.6       rstatix_0.7.2       
 [7] lattice_0.22-6       tzdb_0.5.0           numDeriv_2016.8-1.1 
[10] vctrs_0.6.5          tools_4.4.3          ps_1.9.1            
[13] generics_0.1.4       ucminf_1.2.2         pkgconfig_2.0.3     
[16] Matrix_1.7-2         checkmate_2.3.2      RColorBrewer_1.1-3  
[19] distributional_0.5.0 lifecycle_1.0.4      compiler_4.4.3      
[22] farver_2.1.2         carData_3.0-5        Formula_1.2-5       
[25] pillar_1.11.0        car_3.1-3            MASS_7.3-64         
[28] abind_1.4-8          nlme_3.1-167         tidyselect_1.2.1    
[31] mvtnorm_1.3-3        stringi_1.8.7        rprojroot_2.0.4     
[34] grid_4.4.3           cli_3.6.5            magrittr_2.0.3      
[37] broom_1.0.8          withr_3.0.2          scales_1.4.0        
[40] backports_1.5.0      timechange_0.3.0     matrixStats_1.5.0   
[43] ggsignif_0.6.4       hms_1.1.3            coda_0.19-4.1       
[46] evaluate_1.0.4       knitr_1.50           rlang_1.1.6         
[49] glue_1.8.0           rstudioapi_0.17.1    jsonlite_2.0.0
```