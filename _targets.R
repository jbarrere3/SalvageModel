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
                 "FactoMineR", "ade4", "factoextra", "xtable")
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
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP)),
  
  # Extract species information (genus, family, order)
  tar_target(species, get_species_info(FUNDIV_tree)),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 2 - Prepare data for the model ----
  
  # Keep variables relevant to the model
  tar_target(data_model_dqm, format_data_model_dqm(FUNDIV_tree, FUNDIV_plot, Climate, species)), 
  
  # Scale numeric variables
  tar_target(data_model_dqm_scaled, scale_data_model_dqm(data_model_dqm)), 
  
  # Prepare data as model input
  tar_target(data_jags_dqm, generate_data_jags_dqm(data_model_dqm_scaled)),
  
  # Add weight or climate variables
  tar_target(data_jags_dqm_w, generate_data_jags_dqm2(data_model_dqm, data_model_dqm_scaled)),
  tar_target(data_jags_dqm_w2, generate_data_jags_dqm3(data_model_dqm, data_model_dqm_scaled)),
  tar_target(data_jags_dqm_climate, generate_data_jags_dqm3_climate(data_model_dqm, data_model_dqm_scaled)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit ----
  
  # Fit model with dominance and extract parameters
  tar_target(jags.model_dqm3, fit_mortality_dqm3(
    data_jags_dqm$data_jags, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 1)),
  
  # Fit model with dominance and extract disturbance intensity
  tar_target(jags.model_dqm3_I, fit_mortality_dqm3(
    data_jags_dqm$data_jags, n.chains = 3, n.iter = 50, n.burn = 10, n.thin = 1, 
    param.in = c("Istorm", "Ifire", "Iother"))),
  
  # Fit model with dominance and weight taken into account
  tar_target(jags.model_dqm4, fit_mortality_dqm4(
    data_jags_dqm_w2$data_jags, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 1)),
  
  # Fit model with dominance and a climate effect
  tar_target(jags.model_dqm3_climate, fit_mortality_dqm3_climate(
    data_jags_dqm_climate$data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 1)),
  
  # Fit stan models
  tar_target(stan_file_dqm, "stan/model_dqm.stan", format = "file"),
  tar_target(stan_file_dqm_noweight, "stan/model_dqm_noweight.stan", format = "file"),
  tar_target(data_stan_dqm, data_jags_dqm_w$data_jags[c('Ntrees', 'Nspecies', 'Nplot', 'plot', 'sp', 'w', 'time', 
                                                          'dbh', 'domi', 'Dstorm', 'Dother', 'Dfire', 'd')]), 
  tar_target(stan.model_dqm, stan(file=stan_file_dqm, data=data_stan_dqm, warmup = 50, iter=500, chains = 3, cores = 3,
                                  pars = c(paste0("st", c(0:5)), paste0("ot", c(0:5)), paste0("fi", c(0:5))))),
  tar_target(stan.model_dqm_noweight, stan(file=stan_file_dqm_noweight, data=data_stan_dqm, warmup = 50, iter=200, chains = 3, cores = 3,
                                    pars = c(paste0("st", c(0:5)), paste0("ot", c(0:5)), paste0("fi", c(0:5))))),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Extract information on how disturbances affect each species (so that parameters that can't be estimated are not shown)
  tar_target(disturbance_species_info_dqm, get_disturbance_species_info_climate(data_model_dqm)),

  # Plot convergence 
  tar_target(fig_convergence_dqm3, plot_convergence_climate(jags.model_dqm3, data_jags_dqm, disturbance_species_info_dqm, 
                                                            "fig/real_data/multispecies_submodel_dqm3/convergence"), 
             format = "file"),
  tar_target(fig_convergence_dqm4, plot_convergence_climate(jags.model_dqm4, data_jags_dqm_w2, disturbance_species_info_dqm, 
                                                            "fig/real_data/multispecies_submodel_dqm4/convergence"), 
             format = "file"),
  tar_target(fig_convergence_dqm3_climate, plot_convergence_climate(jags.model_dqm3_climate, data_jags_dqm_climate, disturbance_species_info_dqm, 
                                                            "fig/real_data/multispecies_submodel_dqm3_climate/convergence"), 
             format = "file"),
  tar_target(fig_convergence_stan_dqm, plot_convergence_stan(stan.model_dqm, data_jags_dqm$species_table, disturbance_species_info_dqm, 
                                                            "fig/real_data/stan_model_dqm/convergence"), 
             format = "file"),
  tar_target(fig_convergence_stan_dqm_noweight, plot_convergence_stan(stan.model_dqm_noweight, data_jags_dqm$species_table, disturbance_species_info_dqm, 
                                                            "fig/real_data/stan_model_dqm_noweight/convergence"), 
             format = "file"),
  
  
  # Plot the predictions of the model
  tar_target(fig_prediction_dqm3, plot_prediction_dqm3(jags.model_dqm3, data_jags_dqm, data_model_dqm_scaled, data_model_dqm,
                                                       disturbance_species_info_dqm, "fig/real_data/multispecies_submodel_dqm3/predictions_alltypes.png"), 
            format = "file"),
  tar_target(fig_prediction_dqm4, plot_prediction_dqm3(jags.model_dqm4, data_jags_dqm_w2, data_model_dqm_scaled, data_model_dqm,
                                                       disturbance_species_info_dqm, "fig/real_data/multispecies_submodel_dqm4/predictions_alltypes.png"), 
             format = "file"),
  tar_target(fig_prediction_stan_dqm, plot_prediction_stan_dqm(
    stan.model_dqm, data_jags_dqm$species_table, disturbance_species_info_dqm,
    data_model_dqm_scaled, data_model_dqm, "fig/real_data/stan_model_dqm/predictions.png"), 
    format = "file"),
  tar_target(fig_prediction_stan_dqm_noweight, plot_prediction_stan_dqm(
    stan.model_dqm_noweight, data_jags_dqm$species_table, disturbance_species_info_dqm,
    data_model_dqm_scaled, data_model_dqm, "fig/real_data/stan_model_dqm_noweight/predictions.png"), 
    format = "file"),
  
  # Observation vs prediction
  tar_target(fig_prediction_vs_observation_dqm3, plot_prediction_vs_observation_dqm3(
    jags.model_dqm3, jags.model_dqm3_I, data_jags_dqm, data_model_dqm_scaled, data_model_dqm, 
    "fig/real_data/multispecies_submodel_dqm3/observation_vs_prediction.png")),
  
  # Distribution of disturbance intensity
  tar_target(fig_disturbance_distribution_dqm, plot_disturbance_intensity_distribution(
    jags.model_dqm3_I, data_jags_dqm, "fig/real_data/multispecies_submodel_dqm3/disturbance_distribution.png")),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 5 - Exploratory plots ----
  
  tar_target(fig_disturbed_trees_per_species, 
             plot_disturbed_trees_per_species(data_model_dqm, "fig/exploratory/disturbed_trees_per_species.png"), 
             format = "file"), 
  tar_target(fig_sgdd_species_disturbance, plot_climate_per_species_per_disturbance(
    data_model_dqm, "sgdd", "fig/exploratory/sgdd_species_disturbance.png"), 
    format = "file"), 
  tar_target(fig_wai_species_disturbance, plot_climate_per_species_per_disturbance(
    data_model_dqm, "wai", "fig/exploratory/wai_species_disturbance.png"), 
    format = "file"), 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 - Traits analysis ----
  
  # data files
  tar_target(wood.density_file, "data/traits/GlobalWoodDensityDatabase.xls", format = "file"),
  tar_target(shade.tolerance_file, "data/traits/shade_tolerance_FrenchNFI.csv", format = "file"),
  tar_target(root.depth_file, "data/traits/GRooTAggregateSpeciesVersion.csv", format = "file"),
  tar_target(bark.thickness_file, "data/traits/bark_thickness_NFI.csv", format = "file"),
  tar_target(TRY_file, "data/traits/TRY_data_request_21092.txt", format = "file"),
  tar_target(gbif_file, "data/traits/sp_gbif_climate.csv", format = "file"),
  
  # Compile traits data
  tar_target(traits_TRY, compile_traits_TRY(TRY_file, data_jags_dqm$species_table$species)),
  tar_target(traits_dqm, compile_traits(wood.density_file, shade.tolerance_file, root.depth_file, FUNDIV_tree, 
                                        bark.thickness_file, data_jags_dqm$species_table$species)),
  
  # Get disturbance sensitivity
  tar_target(disturbance_sensitivity_dqm, get_disturbance_sensivity_dqm(
    jags.model_dqm3, data_jags_dqm, data_model_dqm_scaled, data_model_dqm, disturbance_species_info_dqm)),
  
  # Plot PCA of traits against sensitivity with dominance model
  # - storm
  tar_target(pca_storm_sensitivity, plot_pca_traits_vs_sensitivity(
    traits = left_join(traits_dqm, traits_TRY, by = "species")[, c("species", "wood.density_g.cm3", "Root_mass_fraction", "Rooting_depth", "height.dbh.ratio")],
    disturbance_sensitivity = disturbance_sensitivity_dqm[, c("species", "storm.dbh250")],
    file.in = "fig/real_data/multispecies_submodel_dqm3/PCAtraits/storm.jpg", disturbance.in = "storm"), 
    format = "file"),
  # - other
  tar_target(pca_other_sensitivity, plot_pca_traits_vs_sensitivity(
    traits = traits_TRY[, c("species", "TRY_leaf.Chl.content_gm-2", "TRY_leaf.CN.ratio_g.cm3", "TRY_leaf.N.mass_mg.g", 
                            "TRY_leaf.P.mass_mg.g", "TRY_leaf.sla_mm2mg-1", "TRY_leaf.thickness_mm")],
    disturbance_sensitivity = disturbance_sensitivity_dqm[, c("species", "other.dbh250")],
    file.in = "fig/real_data/multispecies_submodel_dqm3/PCAtraits/other.jpg", disturbance.in = "other"), 
    format = "file"), 
  # - fire
  tar_target(pca_fire_sensitivity, plot_pca_traits_vs_sensitivity(
    traits = left_join(traits_dqm, traits_TRY, by = "species")[, c("species", "wood.density_g.cm3", "TRY_leaf.thickness_mm", 
                                                                  "TRY_stomata.conductance_millimolm-2s-1")],
    disturbance_sensitivity = disturbance_sensitivity_dqm[, c("species", "fire.dbh250")],
    file.in = "fig/real_data/multispecies_submodel_dqm3/PCAtraits/fire.jpg", disturbance.in = "fire"), 
    format = "file"), 
  
  # Plot disturbance sensitivity against climate optimum
  tar_target(fig_climate_gbif_vs_sensitivity, plot_climate_gbif_vs_sensitivity(
    disturbance_sensitivity_dqm, gbif_file, "fig/real_data/multispecies_submodel_dqm3/gbif_vs_sensitivity.jpg"), 
    format = "file"), 
  
  # Export some results as latex files
  tar_target(trait_results_latex, 
             export_trait_result_latex(traits = left_join(traits_dqm, traits_TRY, by = "species"), 
                                       disturbance_sensitivity = dplyr::select(disturbance_sensitivity_dqm, "species", "fire" = "fire.dbh250", 
                                                                               "storm" = "storm.dbh250", "other" = "other.dbh250"), 
                                       file.in = "output/model_traits_dqm.tex"), 
             format = "file")
  
  
)