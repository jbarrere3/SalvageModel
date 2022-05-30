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
packages.in <- c("dplyr", "ggplot2", "RCurl", "httr", "tidyr", "data.table", "sp", "R2jags", 
                 "ggmcmc", "taxize", "rnaturalearth", "ggspatial", "sf", "ggnewscale", "readxl", "scales")
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
  tar_target(data_model_full, format_data_model_full(FUNDIV_tree, FUNDIV_plot, Climate, BM_equations, species)),
  
  # Scale numeric variables
  tar_target(data_model_full_scaled, scale_data_model(data_model_full, var = c("dbh", "comp", "sgdd", "wai", "DS"))),
  
  # Prepare data as model input
  tar_target(data_jags_full_sub, generate_data_jags_full_sub(data_model_full_scaled)),
  
  # Simulate some data
  tar_target(parameters_sp, generate_parameters_sp(data_jags_full_sub)),
  tar_target(data_jags_full_sub_simulated, simulate_status_full_sub(data_jags_full_sub, parameters_sp)),
  
  # Prepare data for the model with quadratic diameter
  tar_target(data_model_dqm, format_data_model_dqm(FUNDIV_tree, FUNDIV_plot, Climate, species)), 
  tar_target(data_model_dqm_scaled, scale_data_model_dqm(data_model_dqm)), 
  tar_target(data_jags_dqm, generate_data_jags_dqm(data_model_dqm_scaled)),
  
  
  # Simulate data for a model with a climate effect
  # tar_target(data_jags_full_sub_climate, generate_data_jags_full_sub_climate(data_model_full_scaled)),
  # tar_target(parameters_sp_climate, generate_parameters_sp_climate(data_jags_full_sub_climate)),
  # tar_target(data_jags_full_sub_climate_simulated,
  #            simulate_status_full_sub_climate(data_jags_full_sub_climate, parameters_sp_climate)),
  # 
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit and output ----
  
  # Fit the jags model
  # - With France and Spain
  tar_target(jags.model_full_sub, fit_mortality_full_sub(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 2000, n.burn = 500, n.thin = 1, param.in = paste0("c", c(0:8)))), 
  tar_target(jags.model_full_sub_I, fit_mortality_full_sub(
    data_jags_full_sub$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 10, param.in = c("Istorm", "Ifire", "Iother"))), 
  # tar_target(jags.model_full_sub_climate, fit_mortality_full_sub_climate(
  #   data_jags_full_sub_climate$data_jags, n.chains = 3, n.iter = 2000, n.burn = 500, n.thin = 1)),
  # tar_target(jags.model_full_sub_climate2, fit_mortality_full_sub_climate2(
  #   data_jags_full_sub_climate$data_jags, n.chains = 3, n.iter = 2000, n.burn = 500, n.thin = 1)),
  # # - With France and Spain and simulated data
  tar_target(jags.model_full_sub_simulated, fit_mortality_full_sub_simulated(
    data_jags_full_sub_simulated$data_jags, n.chains = 3, n.iter = 2000, n.burn = 500, n.thin = 1)), 
  # - With simulated data and climate included
  # tar_target(jags.model_full_sub_climate_simulated, fit_mortality_full_sub_climate_simulated(
  #   data_jags_full_sub_climate_simulated$data_jags, n.chains = 3, n.iter = 2000, n.burn = 500, n.thin = 1)),
  # 
  
  # Fit models with quadratic diameter included
  tar_target(jags.model_dqm1, fit_mortality_dqm1(
    data_jags_dqm$data_jags, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 1)),
  tar_target(jags.model_dqm2, fit_mortality_dqm2(
    data_jags_dqm$data_jags, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 1)),
  tar_target(jags.model_dqm3, fit_mortality_dqm3(
    data_jags_dqm$data_jags, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 1)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Extract information on how disturbances affect each species (so that parameters that can't be estimated are not shown)
  tar_target(disturbance_species_info, get_disturbance_species_info(data_model_full)),
  tar_target(disturbance_species_info_dqm, get_disturbance_species_info_climate(data_model_dqm)),
  # tar_target(disturbance_species_info_climate, get_disturbance_species_info_climate(data_model_full)),
  
  # Plot convergence 
  tar_target(fig_convergence_full_sub, plot_convergence(jags.model_full_sub, data_jags_full_sub, BM_equations, 
                                                        disturbance_species_info, "fig/real_data/multispecies_submodel_full/convergence"), 
             format = "file"),
  tar_target(fig_convergence_full_sub_simulated, plot_convergence(jags.model_full_sub_simulated, data_jags_full_sub_simulated, BM_equations, 
                                                                  disturbance_species_info, "fig/simulated_data/multispecies_submodel_full/convergence"), 
             format = "file"),
  tar_target(fig_convergence_dqm1, plot_convergence_climate(jags.model_dqm1, data_jags_dqm, disturbance_species_info_dqm, 
                                                            "fig/real_data/multispecies_submodel_dqm1/convergence"), 
             format = "file"),
  tar_target(fig_convergence_dqm2, plot_convergence_climate(jags.model_dqm2, data_jags_dqm, disturbance_species_info_dqm, 
                                                            "fig/real_data/multispecies_submodel_dqm2/convergence"), 
             format = "file"),
  tar_target(fig_convergence_dqm3, plot_convergence_climate(jags.model_dqm3, data_jags_dqm, disturbance_species_info_dqm, 
                                                            "fig/real_data/multispecies_submodel_dqm3/convergence"), 
             format = "file"),
  # tar_target(fig_convergence_full_sub_climate, plot_convergence_climate(jags.model_full_sub_climate, data_jags_full_sub_climate, 
  #                                                                                 disturbance_species_info_climate, "fig/real_data/multispecies_submodel_full_climate/convergence"), 
  #            format = "file"),
  # tar_target(fig_convergence_full_sub_climate2, plot_convergence_climate(jags.model_full_sub_climate2, data_jags_full_sub_climate, 
  #                                                                       disturbance_species_info_climate, "fig/real_data/multispecies_submodel_full_climate2/convergence"), 
  #            format = "file"),
  # tar_target(fig_convergence_full_sub_climate_simulated, plot_convergence_climate(jags.model_full_sub_climate_simulated, data_jags_full_sub_climate_simulated, 
  #                                                                 disturbance_species_info_climate, "fig/simulated_data/multispecies_submodel_full_climate/convergence"), 
  #            format = "file"),
  # 
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
  tar_target(fig_prediction_full2, plot_prediction_full2(jags.model_full_sub, data_jags_full_sub, data_model_full_scaled, data_model_full,
                                                       disturbance_species_info, "fig/real_data/multispecies_submodel_full/predictions/alltypes.png"), 
             format = "file"),
  tar_target(fig_prediction_dqm1, plot_prediction_dqm1(jags.model_dqm1, data_jags_dqm, data_model_dqm_scaled, data_model_dqm,
                                                          disturbance_species_info_dqm, "fig/real_data/multispecies_submodel_dqm1/predictions_alltypes.png"), 
              format = "file"),
  tar_target(fig_prediction_dqm2, plot_prediction_dqm2(jags.model_dqm2, data_jags_dqm, data_model_dqm_scaled, data_model_dqm,
                                                         disturbance_species_info_dqm, "fig/real_data/multispecies_submodel_dqm2/predictions_alltypes.png"), 
            format = "file"),
  tar_target(fig_prediction_dqm3, plot_prediction_dqm3(jags.model_dqm3, data_jags_dqm, data_model_dqm_scaled, data_model_dqm,
                                                       disturbance_species_info_dqm, "fig/real_data/multispecies_submodel_dqm3/predictions_alltypes.png"), 
            format = "file"),
  # 
  # Estimated intensity vs severity
  tar_target(fig_intensity_vs_severity, plot_intensity_vs_severity(
    jags.model_full_sub_I, data_jags_full_sub, "fig/real_data/multispecies_submodel_full/intensity_vs_severity.png"), 
    format = "file"),
  tar_target(fig_map_intensity_and_severity, map_intensity_and_severity(
    jags.model_full_sub_I, data_jags_full_sub, FUNDIV_plot, 
    "fig/real_data/multispecies_submodel_full/map_intensity_and_severity.png"), 
    format = "file"),
  
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
    format = "file"), 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 - Traits analysis ----
  
  # data files
  tar_target(wood.density_file, "data/traits/GlobalWoodDensityDatabase.xls", format = "file"),
  tar_target(shade.tolerance_file, "data/traits/shade_tolerance_FrenchNFI.csv", format = "file"),
  tar_target(root.depth_file, "data/traits/GRooTAggregateSpeciesVersion.csv", format = "file"),
  tar_target(TRY_file, "data/traits/TRY_data_request_21092.txt", format = "file"),
  tar_target(gbif_file, "data/traits/sp_gbif_climate.csv", format = "file"),
  
  # Compile traits data
  tar_target(traits, compile_traits(wood.density_file, shade.tolerance_file, root.depth_file, 
                                    data_jags_full_sub$species_table$species)),
  tar_target(traits_TRY, compile_traits_TRY(TRY_file, data_jags_full_sub$species_table$species)),
  tar_target(traits_climate, compile_traits_climate(data_model_full, gbif_file)),
  
  # Get disturbance sensitivity
  tar_target(disturbance_sensitivity, get_disturbance_sensivity(
    jags.model_full_sub, data_jags_full_sub, data_model_full_scaled, data_model_full, disturbance_species_info)), 
  
  # Plot disturbance sensitivity against trait values
  tar_target(fig_sensitivity_vs_traits, plot_sensitivity_vs_traits(
    traits, disturbance_sensitivity, "fig/real_data/multispecies_submodel_full/traits/traits"), 
    format = "file"), 
  tar_target(fig_sensitivity_vs_traits_TRY, plot_sensitivity_vs_traits(
    traits_TRY, disturbance_sensitivity, "fig/real_data/multispecies_submodel_full/traits/traits_TRY"), 
    format = "file"), 
  tar_target(fig_sensitivity_vs_traits_climate, plot_sensitivity_vs_traits(
    traits_climate, disturbance_sensitivity, "fig/real_data/multispecies_submodel_full/traits/traits_climate"), 
    format = "file")
  
  
  
  
)