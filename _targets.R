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
        clustermq.scheduler = "multiprocess", 
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
  tar_target(Climate, fread(grep("NFI_climate", datafiles, value = TRUE))), 
  
  # Merge data from the different NFI
  tar_target(FUNDIV_tree, rbind(FUNDIV_tree_FR, FUNDIV_tree_SP)),
  tar_target(FUNDIV_plot, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP)),
  tar_target(FUNDIV_plot_bis, rbind(FUNDIV_plot_FR, FUNDIV_plot_SP_bis)),
  
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
  tar_target(jags.model, fit_mortality(data_jags, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  tar_target(jags.model_bis, fit_mortality(data_jags_bis, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  
  # Fit the model with stocking
  tar_target(jags.model_stock, fit_mortality_stock(data_jags_stock, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  tar_target(jags.model_stock_bis, fit_mortality_stock(data_jags_stock_bis, n.chains = 3, n.iter = 1000, n.burn = 200, n.thin = 20)), 
  
  # Extract model predictions
  tar_target(disturbance_sensitivity, get_disturbance_sensivity(jags.model, data_jags, data_model, 
                                                                dbh.ref = 250, I.ref = 0.6)), 
  tar_target(disturbance_sensitivity_bis, get_disturbance_sensivity(jags.model_bis, data_jags_bis, data_model_bis, 
                                                                dbh.ref = 250, I.ref = 0.6)), 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Plot convergence 
  tar_target(fig_convergence, plot_convergence(jags.model, data_jags, "fig/real_data/reference/convergence"), 
             format = "file"),
  tar_target(fig_convergence_bis, plot_convergence(jags.model_bis, data_jags_bis, "fig/real_data/reference_bis/convergence"), 
             format = "file"),
  tar_target(fig_convergence_stock, plot_convergence(jags.model_stock, data_jags_stock, "fig/real_data/reference_stock/convergence"), 
             format = "file"),
  tar_target(fig_convergence_stock_bis, plot_convergence(jags.model_stock_bis, data_jags_stock_bis, "fig/real_data/reference_stock_bis/convergence"), 
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
  tar_target(fig_param_per_species_bis, plot_param_per_species(
    jags.model_bis, data_jags_bis, data_model_bis, 
    file.in = "fig/real_data/reference_bis/param_per_species.png"), 
    format = "file"),
  tar_target(fig_param_per_species_stock, plot_param_per_species_stock(
    jags.model_stock, data_jags_stock, data_model, 
    file.in = "fig/real_data/reference_stock/param_per_species.png"), 
    format = "file"),
  tar_target(fig_param_per_species_stock_bis, plot_param_per_species_stock(
    jags.model_stock_bis, data_jags_stock_bis, data_model_bis, 
    file.in = "fig/real_data/reference_stock_bis/param_per_species.png"), 
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
  tar_target(gbif_disturbance_file, "data/traits/sp_gbif_disturbance.csv", format = "file"),
  
  ## - Compile traits data
  tar_target(traits_TRY, compile_traits_TRY(TRY_file, get_species_list(data_model))),
  tar_target(traits, compile_traits(wood.density_file, shade.tolerance_file, root.depth_file, 
                                    FUNDIV_tree, bark.thickness_file, get_species_list(data_model))), 
  
  
  ##  Make trait by trait regressions
  # -- For reference model and regular traits
  tar_target(fig_trait_vs_sensitivity, plot_traits_vs_sensitivity_allDist(
    traits, disturbance_sensitivity, "fig/real_data/reference/traits"), format = "file"), 
  # -- For reference model and TRY traits
  tar_target(fig_traitTRY_vs_sensitivity, plot_traits_vs_sensitivity_allDist(
    traits_TRY, disturbance_sensitivity, "fig/real_data/reference/traits_TRY"), format = "file"), 
  # -- For reference bis model and regular traits
  tar_target(fig_trait_vs_sensitivity_bis, plot_traits_vs_sensitivity_allDist(
    traits, disturbance_sensitivity_bis, "fig/real_data/reference_bis/traits"), format = "file"), 
  # -- For reference bis model and TRY traits
  tar_target(fig_traitTRY_vs_sensitivity_bis, plot_traits_vs_sensitivity_allDist(
    traits_TRY, disturbance_sensitivity_bis, "fig/real_data/reference_bis/traits_TRY"), format = "file"), 
  
  ##  Climate vs sensitivity regressions
  # -- For reference model
  tar_target(fig_rda_gbif, plot_rda_climate(disturbance_sensitivity, gbif_file, gbif_disturbance_file, 
                                            "fig/real_data/reference/rda_climate_vs_sensitivity.jpg"), 
             format = "file"),
  # -- For reference model bis
  tar_target(fig_rda_gbif_bis, plot_rda_climate(disturbance_sensitivity_bis, gbif_file, gbif_disturbance_file, 
                                            "fig/real_data/reference_bis/rda_climate_vs_sensitivity.jpg"), 
             format = "file"),
  
  ## Extract trait analysis in a table
  tar_target(table_result_trait, export_trait_result_latex(traits, traits_TRY, disturbance_sensitivity,
                                                           "table/reference_vs_traits.tex"), 
             format = "file"),
  tar_target(table_result_trait_bis, export_trait_result_latex(traits, traits_TRY, disturbance_sensitivity_bis,
                                                           "table/reference_vs_traits_bis.tex"), 
             format = "file"),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 - Exploratory plots ----
  
  # Rate of dead and harvested trees
  tar_target(fig_harvest_death_rate, plot_harvest_and_death_rate(
    data_model, "fig/exploratory/death_harvest_rates.png"), format = "file"), 
  
  # Make a map of disturbance intensity for the reference model
  tar_target(fig_disturbance_intensity, map_disturbance_intensity(
    jags.model, data_jags, FUNDIV_plot, "fig/real_data/reference/map_intensity.png")), 
  tar_target(fig_disturbance_intensity.2, map_disturbance_intensity_ter(
    jags.model, data_jags, FUNDIV_plot, "fig/real_data/reference/map_intensity2.png")), 
  
  # Make a map of disturbance intensity for model with all disturbances
  tar_target(fig_disturbance_intensity_bis, map_disturbance_intensity_bis(
    jags.model_bis, data_jags_bis, FUNDIV_plot_bis, "fig/real_data/reference_bis/map_intensity.png")), 
  tar_target(fig_disturbance_intensity_bis.2, map_disturbance_intensity_ter(
    jags.model_bis, data_jags_bis, FUNDIV_plot_bis, "fig/real_data/reference_bis/map_intensity2.png")), 
  
  # Plot trends in disturbance occurrence and severity over time
  tar_target(fig_disturbance_trends, plot_disturbance_trends(FUNDIV_plot, FUNDIV_tree, "fig/exploratory/disturbance_trends.jpg"), 
             format = "file"), 
  tar_target(fig_disturbance_trends_bis, plot_disturbance_trends(FUNDIV_plot_bis, FUNDIV_tree, "fig/exploratory/disturbance_trends_bis.jpg"), 
             format = "file")
  
)