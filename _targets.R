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
                 "FactoMineR", "ade4", "factoextra", "xtable", "MASS", "vegan")
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
  tar_target(FUNDIV_plot_SP, fread(grep("SpanishNFI_plot", datafiles, value = TRUE))), 
  tar_target(FUNDIV_tree_FI, fread(grep("FinnishNFI_tree", datafiles, value = TRUE))), 
  tar_target(FUNDIV_plot_FI, fread(grep("FinnishNFI_plot", datafiles, value = TRUE))), 
  tar_target(Climate_noFI, fread(grep("NFI_climate", datafiles, value = TRUE))), 
  tar_target(Climate_FI, fread(grep("climate_FI", datafiles, value = TRUE))),
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP, FUNDIV_tree_FI)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP, FUNDIV_plot_FI)),
  tar_target(Climate, rbind(Climate_FI, Climate_noFI)),
  
  # Extract species information (genus, family, order)
  tar_target(species, get_species_info(FUNDIV_tree)),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 2 - Prepare data for the model ----
  
  # Keep variables relevant to the model
  tar_target(data_model, format_data_model(FUNDIV_tree, FUNDIV_plot, Climate, species)), 
  
  # Prepare data as model input
  tar_target(data_jags, generate_data_jags(data_model)), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 3 - Model fit ----
  
  # Fit the reference model
  tar_target(jags.model, fit_mortality(data_jags, n.chains = 3, n.iter = 500, n.burn = 100, n.thin = 10)), 
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Plot convergence 
  tar_target(fig_convergence, plot_convergence(jags.model, data_jags, "fig/real_data/reference/convergence"), 
             format = "file"),
  
  # Plot the predictions of the model
  tar_target(fig_prediction_all, plot_prediction_all(
    jags.model, data_jags, data_model, "fig/real_data/reference/prediction/all.png"), 
    format = "file"),
  tar_target(fig_prediction_all2, plot_prediction2(
    jags.model, data_jags, data_model, "fig/real_data/reference/prediction/all2.png"), 
    format = "file"),
  tar_target(fig_prediction, plot_prediction(
    jags.model, data_jags, data_model, "fig/real_data/reference/prediction"), 
    format = "file"),
  
  
  # Observation vs prediction
  tar_target(fig_predict_vs_obs_meanProba, plot_predicted_vs_observed(
    jags.model, data_jags, data_model, method = "mean.proba", "fig/real_data/reference/validation/meanProba.png"), 
    format = "file"),
  tar_target(fig_predict_vs_obs_meanParam, plot_predicted_vs_observed(
    jags.model, data_jags, data_model, method = "mean.param", "fig/real_data/reference/validation/meanParam.png"), 
    format = "file"),
  
  # Distribution of disturbance intensity
  tar_target(fig_intensity_distribution, plot_intensity_distribution(
    jags.model, data_jags, "fig/real_data/reference/intensity_distribution.png"), 
    format = "file"),
  
  # Plot parameter values
  tar_target(fig_param_per_species, plot_param_per_species(
    jags.model, data_jags, data_model, 
    file.in = "fig/real_data/reference/param_per_species.png"), 
    format = "file"),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 5 - Traits analysis ----
  
  ## - data files
  tar_target(wood.density_file, "data/traits/GlobalWoodDensityDatabase.xls", format = "file"),
  tar_target(shade.tolerance_file, "data/traits/shade_tolerance_FrenchNFI.csv", format = "file"),
  tar_target(root.depth_file, "data/traits/GRooTAggregateSpeciesVersion.csv", format = "file"),
  tar_target(bark.thickness_file, "data/traits/bark_thickness_NFI.csv", format = "file"),
  tar_target(TRY_file, "data/traits/TRY_data_request_21092.txt", format = "file"),
  tar_target(gbif_file, "data/traits/sp_gbif_climate.csv", format = "file"),
  
  ## - Compile traits data
  tar_target(traits_TRY, compile_traits_TRY(TRY_file, get_species_list(data_model))),
  tar_target(traits, compile_traits(wood.density_file, shade.tolerance_file, root.depth_file, 
                                    FUNDIV_tree, bark.thickness_file, get_species_list(data_model))), 
  
  
  ##  Make trait by trait regressions
  # -- Storm
  tar_target(fig_trait_by_trait_storm, plot_traits_vs_sensitivity(
    traits = traits, 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                        dbh.ref = 250, I.ref = 0.6), "storm",
    dir.in = "fig/real_data/reference/traits/trait_by_trait_storm"), 
    format = "file"), 
  tar_target(fig_trait_TRY_storm, plot_traits_vs_sensitivity(
    traits = traits_TRY, 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                        dbh.ref = 250, I.ref = 0.6), "storm",
    dir.in = "fig/real_data/reference/traits/TRY/storm"), 
    format = "file"), 
  # -- Fire
  tar_target(fig_trait_by_trait_fire, plot_traits_vs_sensitivity(
    traits = traits, 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                        dbh.ref = 250, I.ref = 0.6), "fire",
    dir.in = "fig/real_data/reference/traits/trait_by_trait_fire"), 
    format = "file"), 
  tar_target(fig_trait_TRY_fire, plot_traits_vs_sensitivity(
    traits = traits_TRY, 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                        dbh.ref = 250, I.ref = 0.6), "fire",
    dir.in = "fig/real_data/reference/traits/TRY/fire"), 
    format = "file"), 
  # -- Other
  tar_target(fig_trait_by_trait_other, plot_traits_vs_sensitivity(
    traits = traits, 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                        dbh.ref = 250, I.ref = 0.6), "other", 
    dir.in = "fig/real_data/reference/traits/trait_by_trait_other"), 
    format = "file"), 
  tar_target(fig_trait_TRY_other, plot_traits_vs_sensitivity(
    traits = traits_TRY, 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                        dbh.ref = 250, I.ref = 0.6), "other",
    dir.in = "fig/real_data/reference/traits/TRY/other"), 
    format = "file"), 
  
  ##  Climate vs sensitivity regressions
  # -- Storm
  tar_target(fig_rda_gbif_storm, plot_rda_traits_vs_sensitivity(
    traits = fread(gbif_file) %>% dplyr::select(species, mean_mat, mean_tmin, mean_map), 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.6)$storm, 
    file.in = "fig/real_data/reference/gbif/rda_storm.jpg"), 
    format = "file"),
  # -- Other
  tar_target(fig_rda_gbif_other, plot_rda_traits_vs_sensitivity(
    traits = fread(gbif_file) %>% dplyr::select(species, mean_mat, mean_tmin, mean_map), 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.6)$other, 
    file.in = "fig/real_data/reference/gbif/rda_other.jpg"), 
    format = "file"), 
  # -- Fire
  tar_target(fig_rda_gbif_fire, plot_rda_traits_vs_sensitivity(
    traits = fread(gbif_file) %>% dplyr::select(species, mean_mat, mean_tmin, mean_map), 
    disturbance_sensitivity = get_disturbance_sensivity(jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.6)$fire, 
    file.in = "fig/real_data/reference/gbif/rda_fire.jpg"), 
    format = "file"), 
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 - Exploratory plots ----
  
  # Rate of dead and harvested trees
  tar_target(fig_harvest_death_rate, plot_harvest_and_death_rate(
    data_model, "fig/exploratory/death_harvest_rates.png"), format = "file"), 
  
  # Make a map of disturbance intensity
  tar_target(fig_disturbance_intensity, map_disturbance_intensity(
    jags.model, data_jags, FUNDIV_plot, "fig/real_data/reference/map_intensity.png"))
  
)