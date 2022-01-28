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
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", "sp", "R2jags", "ggmcmc")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE)
tar_option_set(packages = packages.in)



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  # Download files
  tar_target(files, get_FrenchNFI(), format = "file"), 
  
  # Load raw data
  tar_target(FrenchNFI_tree_raw, fread("data/FrenchNFI/ARBRE.csv")), 
  tar_target(FrenchNFI_plot_raw, fread("data/FrenchNFI/PLACETTE.csv")), 
  tar_target(FrenchNFI_species, fread("data/FrenchNFI_species.csv")), 
  tar_target(Climate, fread("data/Climate.csv")), 
  tar_target(Disturbance, fread("data/Disturbance.csv")), 
  
  # Format raw data 
  tar_target(FrenchNFI, format_FrenchNFI_raw(FrenchNFI_tree_raw, FrenchNFI_plot_raw)),
  
  # Format to FUNDIV template
  tar_target(FUNDIV_tree_FR, format_FrenchNFI_tree_to_FUNDIV(FrenchNFI, FrenchNFI_species)), 
  tar_target(FUNDIV_plot_FR, format_FrenchNFI_plot_to_FUNDIV(FrenchNFI, FUNDIV_tree_FR)), 
  
  # Keep and scale variables relevant to the model
  tar_target(data_model, format_data_model(FUNDIV_tree_FR, FUNDIV_plot_FR, Disturbance, Climate)), 
  tar_target(data_model_scaled, scale_data_model(data_model)), 
  
  # Get parameters for harvest conditional probabilities
  tar_target(param_harvest_proba, get_param_harvest_proba(data_model_scaled)),
  
  # Generate input data for the model
  tar_target(data_jags, generate_data_jags(subset(data_model_scaled, species == "Abies alba"), 
                                           param_harvest_proba, Dj_latent = TRUE, 
                                           data.type = "simulated")),
  
  # Fit jags model
  tar_target(jags.model, fit_mortality(data_jags, n.chains = 3, n.iter = 100, 
                                       n.burn = 20, n.thin = 1, Dj_latent = T)),
  
  # Diagnostic plots
  tar_target(fig_jags.model_chains, 
             plot_convergence(jags.model, file.in = "fig/model_diagnostic/fig_convergence.png"), 
             format = "file"), 
  tar_target(fig_jags.model_true.vs.fitted.param, 
             plot_fitted_vs_true_parameters(data_jags, jags.model, 
                                            "fig/model_diagnostic/fig_true_vs_fitted_param.png"), 
             format = "file")
)