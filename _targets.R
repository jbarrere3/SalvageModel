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
packages.in <- c("dplyr", "ggplot2", "tidyr", "data.table", "sp", "R2jags", "sf",
                 "ggmcmc", "taxize", "rnaturalearth", "ggspatial", "ggnewscale")
for(i in 1:length(packages.in)) if(!(packages.in[i] %in% rownames(installed.packages()))) install.packages(packages.in[i])
# Targets options
options(tidyverse.quiet = TRUE, clustermq.scheduler = "multicore")
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
  tar_target(FUNDIV_plot_SP, fread(grep("SpanishNFI_plot", datafiles, value = TRUE))), 
  tar_target(FUNDIV_tree_FI, fread(grep("FinnishNFI_tree", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_FI, fread(grep("FinnishNFI_plot", datafiles, value = TRUE))), 
  tar_target(Climate_noFI, fread(grep("NFI_climate", datafiles, value = TRUE))), 
  tar_target(Climate_FI, fread(grep("climate_FI", datafiles, value = TRUE))), 
  tar_target(BM_equations, fread(grep("BM_equations", datafiles, value = TRUE))), 
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP, FUNDIV_tree_FI)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP, FUNDIV_plot_FI)),
  tar_target(Climate, rbind(Climate_FI, Climate_noFI)),
  
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
  
  # Prepare data as model input
  tar_target(data_jags_full_sub, generate_data_jags_full_sub(data_model_full_scaled)),
  
  # Simulate some data
  tar_target(parameters_sp, generate_parameters_sp(data_jags_full_sub)),
  tar_target(data_jags_full_sub_simulated, simulate_status_full_sub(data_jags_full_sub, parameters_sp)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit and output ----
  
  # Fit the jags model
  # - With France and Spain
  tar_target(jags.model_full_sub, fit_mortality_full_sub(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1, param.in = paste0("c", c(0:8)))), 
  tar_target(jags.model_full_sub_I, fit_mortality_full_sub(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 10, param.in = c("Istorm", "Ifire", "Iother"))), 
  # - With France and Spain and simulated data
  tar_target(jags.model_full_sub_simulated, fit_mortality_full_sub_simulated(
    data_jags_full_sub_simulated$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1)), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Extract information on how disturbances affect each species (so that parameters that can't be estimated are not shown)
  tar_target(disturbance_species_info, get_disturbance_species_info(data_model_full)),
  
  # Plot convergence 
  tar_target(fig_convergence_full_sub, plot_convergence(jags.model_full_sub, data_jags_full_sub, BM_equations, 
                                                        disturbance_species_info, "fig/real_data/multispecies_submodel_full/convergence"), 
             format = "file"),
  tar_target(fig_convergence_full_sub_simulated, plot_convergence(jags.model_full_sub_simulated, data_jags_full_sub_simulated, BM_equations, 
                                                                  disturbance_species_info, "fig/simulated_data/multispecies_submodel_full/convergence"), 
             format = "file"),
  
  # Plot parameters per species for real data
  tar_target(fig_param_per_species_full, plot_parameters_per_species_full(jags.model_full_sub, data_jags_full_sub, disturbance_species_info, 
                                                                          "fig/real_data/multispecies_submodel_full/parameters_per_sp.png"), 
             format = "file"),
  
  # Plot true vs estimated parameters for simulated data
  tar_target(fig_param_true_vs_estimated, plot_parameters_true_vs_estimated(
    jags.model_full_sub_simulated, data_jags_full_sub_simulated, parameters_sp, 
    disturbance_species_info, "fig/simulated_data/multispecies_submodel_full/true_vs_estimated.png"), 
    format = "file"),
  tar_target(fig_param_true_vs_estimated2, plot_parameters_true_vs_estimated2(
    jags.model_full_sub_simulated, data_jags_full_sub_simulated, parameters_sp, 
    disturbance_species_info, "fig/simulated_data/multispecies_submodel_full/true_vs_estimated2.png"), 
    format = "file"),
  
  # Plot the predictions of the model
  tar_target(fig_prediction_full, plot_prediction_full(jags.model_full_sub, data_jags_full_sub, data_model_full_scaled, data_model_full,
                                                       disturbance_species_info, "fig/real_data/multispecies_submodel_full/predictions"), 
             format = "file"),
  
  # Estimated intensity vs severity
  tar_target(fig_intensity_vs_severity, plot_intensity_vs_severity(
    jags.model_full_sub_I, data_jags_full_sub, "fig/real_data/multispecies_submodel_full/intensity_vs_severity.png"), 
    format = "file"),
  #tar_target(fig_map_intensity_and_severity, map_intensity_and_severity(
   # jags.model_full_sub_I, data_jags_full_sub, FUNDIV_plot, 
   # "fig/real_data/multispecies_submodel_full/map_intensity_and_severity.png"), 
   # format = "file"),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 5 - Exploratory plots ----
  
  tar_target(fig_disturbed_trees_per_species, 
             plot_disturbed_trees_per_species(data_model_full, "fig/exploratory/disturbed_trees_per_species.png"), 
             format = "file"), 
  tar_target(fig_sgdd_species_disturbance, plot_climate_per_species_per_disturbance(
    data_model_full, "sgdd", "fig/exploratory/sgdd_species_disturbance.png"), 
    format = "file"), 
  tar_target(fig_wai_species_disturbance, plot_climate_per_species_per_disturbance(
    data_model_full, "wai", "fig/exploratory/wai_species_disturbance.png"), 
    format = "file")
  
  
  
)