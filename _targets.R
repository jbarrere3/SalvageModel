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
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", "sp", "R2jags", "ggmcmc", "taxize")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multiprocess")
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
  tar_target(FUNDIV_tree_FR, fread(datafiles[3])), 
  tar_target(FUNDIV_plot_FR, fread(datafiles[2])), 
  tar_target(FUNDIV_tree_SP, fread(datafiles[6])), 
  tar_target(FUNDIV_plot_SP, fread(datafiles[5])), 
  tar_target(Climate, fread(datafiles[4])), 
  tar_target(BM_equations, fread(datafiles[1])), 
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP)),
  
  # Extract species information (genus, family, order)
  tar_target(species, get_species_info(FUNDIV_tree)),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 2 - Prepare data for the model ----
  
  # Keep variables relevant to the model
  tar_target(data_model, format_data_model(FUNDIV_tree_FR, FUNDIV_plot_FR, Climate, BM_equations)), 
  tar_target(data_model_full, format_data_model_full(FUNDIV_tree, FUNDIV_plot, Climate, BM_equations, species)),
  
  # Scale numeric variables
  tar_target(data_model_scaled, scale_data_model(data_model, var = c("dbh", "comp", "sgdd", "wai", "DS"))), 
  tar_target(data_model_full_scaled, scale_data_model(data_model_full, var = c("dbh", "comp", "sgdd", "wai", "DS"))),
  
  # Get parameters for harvest probability
  tar_target(param_harvest_proba, get_param_harvest_proba(data_model_scaled, Dj_latent = FALSE)),
  tar_target(param_harvest_proba_full,  get_param_harvest_proba(data_model_full_scaled, Dj_latent = FALSE)),
  
  # Prepare data as model input
  tar_target(data_jags_sub, generate_data_jags_sub(data_model_scaled, p = param_harvest_proba)),
  tar_target(data_jags_full_sub, 
             generate_data_jags_full_sub(data_model_full_scaled, p = param_harvest_proba_full)),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit and output ----
  
  # Fit the jags model
  # - With France only
  tar_target(jags.model_sub, fit_mortality_sub(
    data_jags_sub$data_jags, n.chains = 3, n.iter = 1000, n.burn = 500, n.thin = 1)), 
  # - With France and Spain
  tar_target(jags.model_full_sub, fit_mortality_full_sub(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 2000, n.burn = 1000, n.thin = 1)), 
  # - With France and Spain and uniform distribution for disturbance intensity
  tar_target(jags.model_full_unif_sub, fit_mortality_full_unif_sub(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 2000, n.burn = 1000, n.thin = 1)), 
  # - Test with latent variable extracted and exponential distribution for power parameter
  tar_target(jags.model_full_sub_test, fit_mortality_full_sub_test(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1)), 
  tar_target(jags.model_full_unif_sub_test, fit_mortality_full_unif_sub_test(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1)), 
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Plot convergence 
  tar_target(fig_convergence_sub, plot_convergence(jags.model_sub, data_jags_sub, BM_equations, 
                                                 "fig/real_data/multispecies_submodel/convergence"), 
             format = "file"),
  tar_target(fig_convergence_full_sub, plot_convergence(jags.model_full_sub, data_jags_full_sub, BM_equations, 
                                                   "fig/real_data/multispecies_submodel_full/convergence"), 
             format = "file"),
  tar_target(fig_convergence_full_unif_sub, plot_convergence(jags.model_full_unif_sub, data_jags_full_sub, BM_equations, 
                                                        "fig/real_data/multispecies_submodel_full_unif/convergence"), 
             format = "file"),
  tar_target(fig_convergence_full_sub_test, plot_convergence(jags.model_full_sub_test, data_jags_full_sub, BM_equations, 
                                                        "fig/real_data/multispecies_submodel_full_test/convergence"), 
             format = "file"),
  tar_target(fig_convergence_full_unif_sub_test, plot_convergence(jags.model_full_unif_sub_test, data_jags_full_sub, BM_equations, 
                                                             "fig/real_data/multispecies_submodel_full_unif_test/convergence"), 
             format = "file"),
  
  # Plot parameters per species
  tar_target(fig_param_per_species, plot_parameters_per_species(jags.model_sub, data_jags_sub, 
                                                                "fig/real_data/multispecies_submodel/parameters_per_sp.png"), 
             format = "file"), 
  tar_target(fig_param_per_species_full, plot_parameters_per_species_full(jags.model_full_sub, data_jags_full_sub, 
                                                                "fig/real_data/multispecies_submodel_full/parameters_per_sp.png"), 
             format = "file"), 
  tar_target(fig_param_per_species_full_unif, plot_parameters_per_species_full(jags.model_full_unif_sub, data_jags_full_sub, 
                                                                          "fig/real_data/multispecies_submodel_full_unif/parameters_per_sp.png"), 
             format = "file"), 
  
  # Plot the predictions of the model
  tar_target(fig_prediction_full, plot_prediction_full(jags.model_full_sub, data_jags_full_sub, data_model_full_scaled, data_model_full,
                                                       "fig/real_data/multispecies_submodel_full/predictions"), 
             format = "file"),
  tar_target(fig_prediction_full_unif, plot_prediction_full(jags.model_full_unif_sub, data_jags_full_sub, data_model_full_scaled, data_model_full,
                                                       "fig/real_data/multispecies_submodel_full_unif/predictions"), 
             format = "file"),
  tar_target(fig_prediction_full_test, plot_prediction_full_test(jags.model_full_sub_test, data_jags_full_sub, data_model_full_scaled, data_model_full,
                                                       "fig/real_data/multispecies_submodel_full_test/predictions"), 
             format = "file"),
  tar_target(fig_prediction_full_unif_test, plot_prediction_full_test(jags.model_full_unif_sub_test, data_jags_full_sub, data_model_full_scaled, data_model_full,
                                                            "fig/real_data/multispecies_submodel_full_unif_test/predictions"), 
             format = "file"),
  
  # Plot disturbance severity vs intensity
  tar_target(fig_intensity_vs_severity_test, 
             plot_intensity_vs_severity(jags.model_full_sub_test, data_jags_full_sub, 
                                        "fig/real_data/testmodel/intensity_vs_severity.png"), 
             format = "file"),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 5 - Exploratory plots ----
  
  tar_target(fig_disturbed_trees_per_species, 
             plot_disturbed_trees_per_species(data_model_full, "fig/exploratory/disturbed_trees_per_species.png"), 
             format = "file")
  
  
  
)