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
  
  # Raw files
  tar_target(FUNDIV_tree_FR_file, "data/FrenchNFI_tree.csv", format = "file"), 
  tar_target(FUNDIV_plot_FR_file, "data/FrenchNFI_plot.csv", format = "file"), 
  tar_target(NFI_mountains_file, "data/NFI_mountains.csv", format = "file"), 
  tar_target(Climate_file, "data/NFI_climate.csv", format = "file"), 
  tar_target(NFI_bm_file, "data/background_morta_fnfi.csv", format = "file"),
  
  # Raw data
  tar_target(FUNDIV_tree_FR, fread(FUNDIV_tree_FR_file)), 
  tar_target(FUNDIV_plot_FR, fread(FUNDIV_plot_FR_file)), 
  tar_target(NFI_mountains, fread(NFI_mountains_file)), 
  tar_target(Climate, fread(Climate_file)), 
  tar_target(NFI_bm, fread(NFI_bm_file)), 
  
  # Subset plots in Alps, Vosges and Jura
  tar_target(FUNDIV_tree, filter(FUNDIV_tree_FR, plotcode %in% NFI_mountains$idp)),
  tar_target(FUNDIV_plot, filter(FUNDIV_plot_FR, plotcode %in% NFI_mountains$idp)),
  
  # Extract species information (genus, family, order)
  tar_target(species, get_species_info(FUNDIV_tree)),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 2 - Prepare data for the model ----

  # Keep variables relevant to the model
  tar_target(data_model, format_data_model(
    FUNDIV_tree, FUNDIV_plot, Climate, species, NFI_bm)),

  # Prepare data as input for the model with logratio effect
  tar_target(data_jags, generate_data_jags(data_model)),
  

  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit ----

  # Fit the reference model
  tar_target(jags.model, fit_mortality(
    data_jags, n.chains = 3, n.iter = 12000, n.burn = 2000, n.thin = 10)),
  

  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Export plots ----

  # Rhat for the reference model
  tar_target(fig_rhat_reference, plot_rhat(
    jags.model, "output/Alps/fig/rhat.png"), format = "file"),

  # Validation (observed vs predicted probabilities)
  tar_target(fig_predict_vs_obs_meanProba_ms, plot_predicted_vs_observed(
    jags.model, data_jags, data_model, method = "mean.proba", 
    "output/Alps/fig/validation.png"), format = "file"),
  
  # Prediciton for several dbh and intensity
  tar_target(fig_predictions, plot_predictions(
    rdata_dominance, "output/Alps/fig/prediction.jpg"), format = "file"),
   
   
   
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 5 - Export tables and rdata----
   
  # Export jags objects and correspondence tables
  tar_target(rdata_dominance, export_jags(
    jags.model.in = jags.model, data_jags.in = data_jags, data_model.in = data_model,
    file.in = "output/Alps/rdata/jags_dominance.Rdata"), format = "file")
)