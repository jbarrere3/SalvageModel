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
                 "ggmcmc", "taxize", "rnaturalearth", "ggspatial", "sf", "ggnewscale", "readxl", "scales", 
                 "FactoMineR", "ade4", "factoextra", "xtable", "MASS", "vegan", "lme4", "car", "GGally", "grid", 
                 "betareg", "Taxonstand", "WorldFlora")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, 
        clustermq.scheduler = "multiprocess", 
        dplyr.summarise.inform = FALSE)
tar_option_set(packages = packages.in,
               memory = "transient")
future::plan(future::multisession, workers = 6)
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
  tar_target(Climate, fread(grep("NFI_climate", datafiles, value = TRUE))), 
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP)),
  tar_target(FUNDIV_plot_bis, 
             subset(rbind(FUNDIV_plot_FR, FUNDIV_plot_SP_bis), 
                    disturbance.nature %in% c("biotic", "snow"))),
  
  # Extract species information (genus, family, order)
  tar_target(species, get_species_info(FUNDIV_tree)),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 2 - Prepare data for the model ----
  
  # Keep variables relevant to the model
  tar_target(data_model, format_data_model(FUNDIV_tree, FUNDIV_plot, Climate, species)), 
  tar_target(data_model_bis, format_data_model(FUNDIV_tree, FUNDIV_plot_bis, Climate, species)), 
  
  # Prepare data as input for the model with logratio effect
  tar_target(data_jags, generate_data_jags(data_model)), 
  tar_target(data_jags_bis, generate_data_jags(data_model_bis)), 
  
  # Prepare data as input for the model with stock effect
  tar_target(data_jags_stock, generate_data_jags_stock(data_model)), 
  tar_target(data_jags_stock_bis, generate_data_jags_stock(data_model_bis)), 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit ----
  
  # Fit the reference model
  tar_target(jags.model, fit_mortality(
    data_jags, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  tar_target(jags.model_bis, fit_mortality(
    data_jags_bis, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  
  # Fit the model with stocking
  tar_target(jags.model_stock, fit_mortality_stock(
    data_jags_stock, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  tar_target(jags.model_stock_bis, fit_mortality_stock(
    data_jags_stock_bis, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  
  # Extract model predictions as a mean over all iterations
  tar_target(disturbance_sensitivity, get_disturbance_sensivity(
    jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.75)), 
  tar_target(disturbance_sensitivity_bis, get_disturbance_sensivity(
    jags.model_bis, data_jags_bis, data_model_bis, dbh.ref = 250, I.ref = 0.75)), 
  
  # Extract model prediction for every iteration
  tar_target(disturbance_sensitivity_full, get_disturbance_sensivity_full(
    jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.75)), 
  tar_target(disturbance_sensitivity_full_bis, get_disturbance_sensivity_full(
    jags.model_bis, data_jags_bis, data_model_bis, dbh.ref = 250, I.ref = 0.75)),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Export plots ----
  
  # Correlation matrix between disturbance sensitivity
  tar_target(fig_correlation_disturbance, plot_correlation_disturbance(
    disturbance_sensitivity, disturbance_sensitivity_bis, "output/fig/correlation_disturbance.png")), 
  
  # Verify that species sensitivity is not related to the mean intensity or the number of trees exposed
  tar_target(fig_sensitivity_vs_intensity_and_ntrees, plot_sensitivity_vs_intensity_and_ntrees(
    c(disturbance_sensitivity, disturbance_sensitivity_bis), c(jags.model, jags.model_bis), 
    c(data_jags, data_jags_bis), "output/fig/sensitivity_vs_intensity_and_ntrees.png"), format = "file"), 
  
  # Rhat for the reference model
  tar_target(fig_rhat_reference, plot_rhat(c(jags.model, jags.model_bis), "output/fig/rhat.png"), format = "file"), 
  
  # Validation (observed vs predicted probabilities)
  tar_target(fig_predict_vs_obs_meanProba_ms, plot_predicted_vs_observed(
    c(jags.model, jags.model_bis), c(data_jags, data_jags_bis), c(data_model, data_model_bis), 
    method = "mean.proba", "output/fig/validation.png"), 
    format = "file"),
  
  # Trends in disturbance frequency and severity
  tar_target(fig_trend_frequency, plot_trend_disturbance_frequency_ms(
    FUNDIV_plot,FUNDIV_plot_bis, "output/fig/dist_freq_trend.jpg"), format = "file"), 
  tar_target(fig_trend_severity, plot_trend_disturbance_severity_ms(
    FUNDIV_tree, FUNDIV_plot, FUNDIV_plot_bis, "output/fig/dist_sever_trend.jpg"), format = "file"),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 5 - Export tables and rdata----
  
  # Extract a table with statistics about disturbance per country
  tar_target(table_disurbance_stat, export_table_disturbance_stats(
    FUNDIV_tree, rbind(FUNDIV_plot, FUNDIV_plot_bis), 
    "output/table/disturbance_stat_bis.tex"), format = "file"),
  
  
  # Export jags objects and correspondence tables
  tar_target(rdata_dominance, export_jags(
    jags.model.in = c(jags.model, jags.model_bis), 
    data_jags.in = c(data_jags, data_jags_bis),
    data_model.in = c(data_model, data_model_bis), 
    file.in = "output/rdata/jags_dominance.Rdata"), format = "file"), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 Run model with corrected weight ----
  
  tar_target(FUNDIV_tree.cor, correct_weight(FUNDIV_tree)), 
  tar_target(data_model.cor, format_data_model(FUNDIV_tree.cor, FUNDIV_plot, Climate, species)), 
  tar_target(data_model_bis.cor, format_data_model(
    FUNDIV_tree.cor, FUNDIV_plot_bis, Climate, species)), 
  tar_target(data_jags.cor, generate_data_jags(data_model.cor)), 
  tar_target(data_jags_bis.cor, generate_data_jags(data_model_bis.cor)), 
  tar_target(jags.model.cor, fit_mortality(
    data_jags.cor, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  tar_target(jags.model_bis.cor, fit_mortality(
    data_jags_bis.cor, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  tar_target(rdata_dominance.cor, export_jags(
    jags.model.in = c(jags.model.cor, jags.model_bis.cor), 
    data_jags.in = c(data_jags.cor, data_jags_bis.cor),
    data_model.in = c(data_model.cor, data_model_bis.cor), 
    file.in = "output/rdata/jags_dominance_cor.Rdata"), format = "file"), 
  tar_target(disturbance_sensitivity_full.cor, get_disturbance_sensivity_full(
    jags.model.cor, data_jags.cor, data_model.cor, dbh.ref = 250, I.ref = 0.75)), 
  tar_target(disturbance_sensitivity_full_bis.cor, get_disturbance_sensivity_full(
    jags.model_bis.cor, data_jags_bis.cor, data_model_bis.cor, dbh.ref = 250, I.ref = 0.75))
)