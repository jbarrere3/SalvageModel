#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name _targets.R  
#' @description R script to launch the target pipeline
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options and packages ----------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Load targets
library(targets)
# Load functions
lapply(grep("R$", list.files("R"), value = TRUE), function(x) source(file.path("R", x)))
# install if needed and load packages
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", "sp", "R2jags", "rstan", "cowplot",
                 "ggmcmc", "taxize", "ggspatial", "sf", "ggnewscale", "readxl", "scales", 
                 "FactoMineR", "ade4", "factoextra", "xtable", "MASS", "vegan", "Taxonstand", "WorldFlora", 
                 "lme4", "car", "GGally", "rnaturalearth", "ggspatial", "grid", "betareg")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, 
        clustermq.scheduler = "multicore", 
        dplyr.summarise.inform = FALSE)
tar_option_set(packages = packages.in)
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 1 - Load data ----
  
  # Raw data
  tar_target(datafiles, paste0("data/", list.files("data")), format = "file"),
  tar_target(FUNDIV_tree_FR, fread(grep("FrenchNFI_tree", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_FR, fread(grep("FrenchNFI_plot", datafiles, value = TRUE))), 
  tar_target(FUNDIV_tree_SP, fread(grep("SpanishNFI_tree", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_SP, fread(grep("SpanishNFI_plot.csv", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_SP_bis, fread(grep("SpanishNFI_plot_bis.csv", datafiles, value = TRUE))), 
  tar_target(FUNDIV_tree_FI, fread(grep("FinnishNFI_tree", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_FI, fread(grep("FinnishNFI_plot.csv", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_FI_bis, fread(grep("FinnishNFI_plot_bis", datafiles, value = TRUE))), 
  tar_target(Climate_noFI, fread(grep("NFI_climate", datafiles, value = TRUE))), 
  tar_target(Climate_FI, fread(grep("climate_FI", datafiles, value = TRUE))),
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP, FUNDIV_tree_FI)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP, FUNDIV_plot_FI)),
  tar_target(FUNDIV_plot_bis, subset(rbind(FUNDIV_plot_FR, FUNDIV_plot_SP_bis, FUNDIV_plot_FI_bis), 
                                     disturbance.nature %in% c("biotic", "snow"))),
  tar_target(Climate, rbind(Climate_FI, Climate_noFI)),
  
  # Extract species information (genus, family, order)
  tar_target(species, get_species_info(FUNDIV_tree)),
  
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Step 2 - Prepare data for the model ----
  # 
  # # Keep variables relevant to the model
  # tar_target(data_model, format_data_model(FUNDIV_tree, FUNDIV_plot, Climate, species)), 
  # tar_target(data_model_bis, format_data_model(FUNDIV_tree, FUNDIV_plot_bis, Climate, species)), 
  # 
  # # Prepare data as input for the model with logratio effect
  # tar_target(data_jags, generate_data_jags(data_model)), 
  # tar_target(data_jags_bis, generate_data_jags(data_model_bis)), 
  # 
  # # Prepare data as input for the model with stock effect
  # tar_target(data_jags_stock, generate_data_jags_stock(data_model)), 
  # tar_target(data_jags_stock_bis, generate_data_jags_stock(data_model_bis)), 
  # 
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Step 3 - Model fit ----
  # 
  # # Fit the reference model
  # tar_target(jags.model, fit_mortality(data_jags, n.chains = 3, n.iter = 5000, 
  #                                      n.burn = 1000, n.thin = 20)), 
  # tar_target(jags.model_bis, fit_mortality(data_jags_bis, n.chains = 3, 
  #                                          n.iter = 5000, n.burn = 1000, n.thin = 20)), 
  # 
  # # Fit the model with stocking
  # tar_target(jags.model_stock, fit_mortality_stock(data_jags_stock, n.chains = 3, n.iter = 5000, n.burn = 1000, n.thin = 20)), 
  # tar_target(jags.model_stock_bis, fit_mortality_stock(data_jags_stock_bis, n.chains = 3, n.iter = 5000, n.burn = 1000, n.thin = 20)), 
  # 
  # # Extract model predictions with a mean over all iterations
  # tar_target(disturbance_sensitivity, get_disturbance_sensivity(
  #   jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.75)), 
  # tar_target(disturbance_sensitivity_bis, get_disturbance_sensivity(
  #   jags.model_bis, data_jags_bis, data_model_bis, dbh.ref = 250, I.ref = 0.75)), 
  # 
  # # Extract model prediction for every iteration
  # tar_target(disturbance_sensitivity_full, get_disturbance_sensivity_full(
  #   jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.75)), 
  # tar_target(disturbance_sensitivity_full_bis, get_disturbance_sensivity_full(
  #   jags.model_bis, data_jags_bis, data_model_bis, dbh.ref = 250, I.ref = 0.75)),
  # 
  # 
  # 
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Step 4 - Export plots ----
  # 
  # # Make a map of disturbance intensity for model with all disturbances
  # tar_target(fig_disturbance_intensity, map_disturbance_intensity_ms(
  #   c(jags.model, jags.model_bis), c(data_jags, data_jags_bis), 
  #   rbind(FUNDIV_plot_bis, FUNDIV_plot), "output/fig/ms/map_intensity.png")), 
  # 
  # # Observation vs prediction
  # tar_target(fig_predict_vs_obs_meanProba_ms, plot_predicted_vs_observed(
  #   c(jags.model, jags.model_bis), c(data_jags, data_jags_bis), c(data_model, data_model_bis), 
  #   method = "mean.proba", "output/fig/ms/validation.png"), 
  #   format = "file"),
  # 
  # # Plot trends in disturbance occurrence and severity over time
  # tar_target(fig_disturbance_trends, plot_disturbance_trends(
  #   FUNDIV_plot, FUNDIV_tree, "output/fig/exploratory/disturbance_trends.jpg"), 
  #            format = "file"), 
  # tar_target(fig_disturbance_trends_bis, plot_disturbance_trends(
  #   FUNDIV_plot_bis, FUNDIV_tree, "output/fig/exploratory/disturbance_trends_bis.jpg"), 
  #            format = "file"), 
  # 
  # # Correlation matrix between disturbance sensitivity
  # tar_target(fig_correlation_disturbance, plot_correlation_disturbance(
  #   disturbance_sensitivity, disturbance_sensitivity_bis, 
  #   "output/fig/ms/correlation_disturbance.png")), 
  # 
  # # Verify that species sensitivity is not related to the mean intensity or the number of trees exposed
  # tar_target(fig_sensitivity_vs_intensity_and_ntrees, plot_sensitivity_vs_intensity_and_ntrees(
  #   c(disturbance_sensitivity, disturbance_sensitivity_bis), c(jags.model, jags.model_bis), 
  #   c(data_jags, data_jags_bis), "output/fig/ms/sensitivity_vs_intensity_and_ntrees.png"), format = "file"), 
  # 
  # 
  # # Rhat for the reference model
  # tar_target(fig_rhat_reference, plot_rhat(
  #   c(jags.model, jags.model_bis), "output/fig/ms/rhat.png"), format = "file"), 
  # 
  # 
  # 
  # 
  # 
  # 
  # #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # # -- Step 5 - Export tables and rdata ----
  # 
  # # Extract a table with statistics about disturbance per country
  # tar_target(table_disurbance_stat, export_table_disturbance_stats(
  #   FUNDIV_tree, rbind(FUNDIV_plot_bis, FUNDIV_plot), 
  #   "output/table/disturbance_stat_bis.tex"), format = "file"),
  # 
  # # Export jags objects and correspondence tables
  # tar_target(rdata_dominance, export_jags(
  #   jags.model.in = c(jags.model, jags.model_bis), data_jags.in = c(data_jags, data_jags_bis), 
  #   data_model.in = c(data_model, data_model_bis), file.in = "output/rdata/jags_dominance.Rdata"), format = "file"), 
  # tar_target(rdata_sensitivity_full, export_sensitivity(
  #   disturbance_sensitivity_full, disturbance_sensitivity_full_bis, "output/rdata/sensitivity_full.Rdata"), format = "file"),
  # 
  # 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 - Fit a model with a mixture effect ----
  
  # Format data
  tar_target(data_model_mix, format_data_model_mix(
    FUNDIV_tree, subset(FUNDIV_plot_bis, disturbance.nature != "snow"), 
    Climate, species)), 
  
  # Prepare as input for jags model
  tar_target(data_jags_mix, generate_data_jags_mix(data_model_mix)), 
  
  # Fit mortality model
  tar_target(jags.model_mix, fit_mortality_mix(
    data_jags_mix, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 10)), 
  
  # Fit mortality model with a binary mixture effect
  tar_target(jags.model_mix_bin, fit_mortality_mix_bin(
    data_jags_mix, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 10)), 
  
  # Extract parameters per iteration
  tar_target(param_mix, get_param_mix(jags.model_mix, data_jags_mix, data_model_mix)),
  tar_target(param_mix_bin, get_param_mix(jags.model_mix_bin, data_jags_mix, data_model_mix)),
  
  # Show rhat values
  tar_target(fig_rhat_mixture, plot_rhat(jags.model_mix, "output/fig_mix/rhat.jpg"), 
             format = "file"), 
  tar_target(fig_rhat_mixture_bin, plot_rhat(jags.model_mix_bin, "output/fig_mix/rhat_bin.jpg"), 
             format = "file"), 
  
  # show param estimation per species
  tar_target(fig_param_mix, plot_mix(
    param_mix, file.in = "output/fig_mix/param_per_species.jpg"), format = "file"), 
  tar_target(fig_param_mix_bin, plot_mix_bin(
    param_mix_bin, file.in = "output/fig_mix/param_per_species_bin.jpg"), format = "file"), 
  
  # Plot the share of conifer per plot
  tar_target(fig_share_conifer, plot_distribution_share_conifer(
    FUNDIV_tree, FUNDIV_plot, FUNDIV_plot_bis, species, 
    "output/fig_mix/share_conifer_per_species.jpg"))
  
  
)