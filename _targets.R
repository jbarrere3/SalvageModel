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
set.seed(2)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Targets workflow --------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list(
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 1 - Load data
  
  # Raw data
  tar_target(datafiles, paste0("data/", list.files("data")), format = "file"),
  tar_target(FUNDIV_tree_FR, fread(datafiles[3])), 
  tar_target(FUNDIV_plot_FR, fread(datafiles[2])), 
  tar_target(FUNDIV_tree_SP, fread(datafiles[6])), 
  tar_target(FUNDIV_plot_SP, fread(datafiles[5])), 
  tar_target(Climate, fread(datafiles[4])), 
  tar_target(BM_equations, fread(datafiles[1])), 
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP)),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 2 - Prepare data for the model
  
  # Keep variables relevant to the model
  tar_target(data_model, format_data_model(FUNDIV_tree_FR, FUNDIV_plot_FR, Climate, BM_equations)), 
  
  # Scale numeric variables
  tar_target(data_model_scaled, scale_data_model(data_model, var = c("dbh", "comp", "sgdd", "wai", "DS"))), 
  
  # Get parameters for harvest probability
  tar_target(param_harvest_proba, get_param_harvest_proba(data_model_scaled, Dj_latent = FALSE)),
  
  # Prepare data as model input
  tar_target(data_jags_sub, generate_data_jags_sub(data_model_scaled, p = param_harvest_proba)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 3 - Model fit and output
  
  # Fit the jags model
  tar_target(jags.model_sub, fit_mortality_sub(data_jags_sub$data_jags, n.chains = 3, n.iter = 2000, n.burn = 1000, n.thin = 1)), 
  
  # Plot convergence 
  tar_target(fig_convergence_sub, plot_convergence(jags.model_sub, data_jags_sub, BM_equations, 
                                                 "fig/real_data/multispecies_submodel/convergence"), 
             format = "file"),
  
  # Plot parameters per species
  tar_target(fig_param_per_species, plot_parameters_per_species(jags.model_sub, data_jags_sub, 
                                                                "fig/real_data/multispecies_submodel/parameters_per_sp.png"), 
             format = "file")
  
)