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
                 "betareg")
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
  tar_target(FUNDIV_plot_bis, 
             subset(rbind(FUNDIV_plot_FR, FUNDIV_plot_SP_bis), disturbance.nature %in% c("biotic", "snow"))),
  
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
                                                                dbh.ref = 250, I.ref = 0.75)), 
  tar_target(disturbance_sensitivity_bis, get_disturbance_sensivity(jags.model_bis, data_jags_bis, data_model_bis, 
                                                                dbh.ref = 250, I.ref = 0.75)), 
  
  # Extract model prediction for every iteration
  tar_target(disturbance_sensitivity_full, get_disturbance_sensivity_full(
    jags.model, data_jags, data_model, dbh.ref = 250, I.ref = 0.75)), 
  tar_target(disturbance_sensitivity_full_bis, get_disturbance_sensivity_full(
    jags.model_bis, data_jags_bis, data_model_bis, dbh.ref = 250, I.ref = 0.75)),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 4 - Plot model outputs ----
  
  # Plot convergence 
  tar_target(fig_convergence, plot_convergence(jags.model, data_jags, "output/fig/real_data/reference/convergence"), 
             format = "file"),
  tar_target(fig_convergence_bis, plot_convergence(jags.model_bis, data_jags_bis, "output/fig/real_data/reference_bis/convergence"), 
             format = "file"),
  tar_target(fig_convergence_stock, plot_convergence(jags.model_stock, data_jags_stock, "output/fig/real_data/reference_stock/convergence"), 
             format = "file"),
  tar_target(fig_convergence_stock_bis, plot_convergence(jags.model_stock_bis, data_jags_stock_bis, "output/fig/real_data/reference_stock_bis/convergence"), 
             format = "file"),
  
  
  # Observation vs prediction
  tar_target(fig_predict_vs_obs_meanProba, plot_predicted_vs_observed(
    jags.model, data_jags, data_model, method = "mean.proba", "output/fig/real_data/reference/validation/meanProba.png"), 
    format = "file"),
  tar_target(fig_predict_vs_obs_meanParam, plot_predicted_vs_observed(
    jags.model, data_jags, data_model, method = "mean.param", "output/fig/real_data/reference/validation/meanParam.png"), 
    format = "file"),
  
  # Distribution of disturbance intensity
  tar_target(fig_intensity_distribution, plot_intensity_distribution(
    jags.model, data_jags, "output/fig/real_data/reference/intensity_distribution.png"), 
    format = "file"),
  
  # Plot parameter values
  tar_target(fig_param_per_species, plot_param_per_species(
    jags.model, data_jags, data_model, 
    file.in = "output/fig/real_data/reference/param_per_species.png"), 
    format = "file"),
  tar_target(fig_param_per_species_bis, plot_param_per_species(
    jags.model_bis, data_jags_bis, data_model_bis, 
    file.in = "output/fig/real_data/reference_bis/param_per_species.jpg"), 
    format = "file"),
  tar_target(fig_param_per_species_stock, plot_param_per_species_stock(
    jags.model_stock, data_jags_stock, data_model, 
    file.in = "output/fig/real_data/reference_stock/param_per_species.png"), 
    format = "file"),
  tar_target(fig_param_per_species_stock_bis, plot_param_per_species_stock(
    jags.model_stock_bis, data_jags_stock_bis, data_model_bis, 
    file.in = "output/fig/real_data/reference_stock_bis/param_per_species.png"), 
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
  tar_target(traits, compile_traits(wood.density_file, shade.tolerance_file, root.depth_file, FUNDIV_tree, 
                                    FUNDIV_plot, bark.thickness_file, get_species_list(data_model))), 
  
  
  ##  Make trait by trait regressions
  # -- For reference model and regular traits
  tar_target(fig_trait_vs_sensitivity, plot_traits_vs_sensitivity_allDist(
    traits, disturbance_sensitivity, "output/fig/real_data/reference/traits"), format = "file"), 
  # -- For reference model and TRY traits
  tar_target(fig_traitTRY_vs_sensitivity, plot_traits_vs_sensitivity_allDist(
    traits_TRY, disturbance_sensitivity, "output/fig/real_data/reference/traits_TRY"), format = "file"), 
  # -- For reference bis model and regular traits
  tar_target(fig_trait_vs_sensitivity_bis, plot_traits_vs_sensitivity_allDist(
    traits, disturbance_sensitivity_bis, "output/fig/real_data/reference_bis/traits"), format = "file"), 
  # -- For reference bis model and TRY traits
  tar_target(fig_traitTRY_vs_sensitivity_bis, plot_traits_vs_sensitivity_allDist(
    traits_TRY, disturbance_sensitivity_bis, "output/fig/real_data/reference_bis/traits_TRY"), format = "file"), 
  
  ##  Climate vs sensitivity regressions
  # -- For reference model
  tar_target(fig_rda_gbif, plot_rda_climate(disturbance_sensitivity, gbif_file, 
                                            "output/fig/real_data/reference/rda_climate_vs_sensitivity.jpg"), 
             format = "file"),
  # -- For reference model bis
  tar_target(fig_rda_gbif_bis, plot_rda_climate(disturbance_sensitivity_bis, gbif_file, 
                                            "output/fig/real_data/reference_bis/rda_climate_vs_sensitivity.jpg"), 
             format = "file"),
  
  ## Extract trait analysis in a table
  tar_target(table_result_trait_all, export_trait_result_latex(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    species, group.in = "all", "output/table/reference_vs_traits_all.tex"), format = "file"),
  tar_target(table_result_trait_conifer, export_trait_result_latex(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    species, group.in = "conifer", "output/table/reference_vs_traits_conifer.tex"), format = "file"),
  tar_target(table_result_trait_broadleaf, export_trait_result_latex(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    species, group.in = "broadleaf", "output/table/reference_vs_traits_broadleaf.tex"), format = "file"),
  tar_target(table_result_trait_allDist, export_trait_allDist_latex(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, "output/table/traits_allDist.tex"), 
    format = "file"),
  
  ## Extract trait analysis with betareg and all mcmc iterations in a table
  tar_target(table_result_trait_full_all, export_trait_result_full_latex(
    traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
    species, group.in = "all", "output/table/reference_vs_traits_full_all.tex"), format = "file"),
  tar_target(table_result_trait_full_broadleaf, export_trait_result_full_latex(
    traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
    species, group.in = "broadleaf", "output/table/reference_vs_traits_full_broadleaf.tex"), format = "file"),
  tar_target(table_result_trait_full_conifer, export_trait_result_full_latex(
    traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
    species, group.in = "conifer", "output/table/reference_vs_traits_full_conifer.tex"), format = "file"),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 6 - Exploratory plots ----
  
  # Rate of dead and harvested trees
  tar_target(fig_harvest_death_rate, plot_harvest_and_death_rate(
    data_model, "output/fig/exploratory/death_harvest_rates.png"), format = "file"), 
  
  # Make a map of disturbance intensity for the reference model
  tar_target(fig_disturbance_intensity, map_disturbance_intensity(
    jags.model, data_jags, FUNDIV_plot, "output/fig/real_data/reference/map_intensity.png")), 
  tar_target(fig_disturbance_intensity.2, map_disturbance_intensity_ter(
    jags.model, data_jags, FUNDIV_plot, "output/fig/real_data/reference/map_intensity2.png")), 
  
  # Plot trends in disturbance occurrence and severity over time
  tar_target(fig_disturbance_trends, plot_disturbance_trends(FUNDIV_plot, FUNDIV_tree, "output/fig/exploratory/disturbance_trends.jpg"), 
             format = "file"), 
  tar_target(fig_disturbance_trends_bis, plot_disturbance_trends(FUNDIV_plot_bis, FUNDIV_tree, "output/fig/exploratory/disturbance_trends_bis.jpg"), 
             format = "file"), 
  
  # Extract a table with statistics about disturbance per country
  tar_target(table_disurbance_stat, export_table_disturbance_stats(
    FUNDIV_tree, rbind(FUNDIV_plot, FUNDIV_plot_bis), "output/table/disturbance_stat_bis.tex"), format = "file"),
  
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 7 - Plots for the manuscript ----
  
  # Climate vs disturbance-related climatic indices
  tar_target(fig_disturbance_climate_ms, plot_disturbance_climate_ms(
    disturbance_sensitivity_full, disturbance_sensitivity_full_bis, gbif_disturbance_file, 
    "output/fig/ms/sensity_vs_climate_disturbance.jpg"), format = "file"),
  
  # PCA of the climatic variables
  tar_target(fig_pca_ms, plot_pca_ms(gbif_file, "output/fig/ms/pca.jpg"), format = "file"),
  
  # RDA of sensitivity vs climate
  tar_target(fig_rda_climate_ms, plot_rda_climate_ms(
    disturbance_sensitivity, disturbance_sensitivity_bis, gbif_file, 
    "output/fig/ms/sensitivity_vs_climate.jpg"), format = "file"), 
  
  # Make a map of disturbance intensity for model with all disturbances
  tar_target(fig_disturbance_intensity_bis.2, map_disturbance_intensity_ms(
    c(jags.model, jags.model_bis[c("snow", "biotic")]), c(data_jags, data_jags_bis[c("snow", "biotic")]), 
    rbind(FUNDIV_plot_bis, FUNDIV_plot), "output/fig/ms/map_intensity.png")), 
  
  # Traits vs sensitivity all disturbance together
  tar_target(fig_traits_vs_sensitivity_ms, plot_traits_vs_sensitivity_ms(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    "output/fig/ms/traits_vs_sensitivity.jpg"), format = "file"), 
  
  # Correlation matrix between disturbance sensitivity
  tar_target(fig_correlation_disturbance, plot_correlation_disturbance(
    disturbance_sensitivity, disturbance_sensitivity_bis, "output/fig/ms/correlation_disturbance.png")), 
  
  # Verify that species sensitivity is not related to the mean intensity or the number of trees exposed
  tar_target(fig_sensitivity_vs_intensity_and_ntrees, plot_sensitivity_vs_intensity_and_ntrees(
    c(disturbance_sensitivity, disturbance_sensitivity_bis), c(jags.model, jags.model_bis), 
    c(data_jags, data_jags_bis), "output/fig/ms/sensitivity_vs_intensity_and_ntrees.png"), format = "file"), 
  
  # Rhat for the reference model
  tar_target(fig_rhat_reference, plot_rhat(c(jags.model, jags.model_bis), "output/fig/ms/rhat.png"), format = "file"), 
  
  # Parameter per species "clean"
  tar_target(fig_param_per_species_ms, plot_param_per_species_ms(
    c(jags.model, jags.model_bis), c(data_jags, data_jags_bis), c(data_model, data_model_bis), 
    file.in = "output/fig/ms/param_per_species.jpg"), format = "file"),
  
  # Map of the disturbances
  tar_target(fig_map_disturbances_ms, map_disturbances_ms(
    FUNDIV_plot, FUNDIV_plot_bis, "output/fig/ms/map_disturbances.jpg"), format = "file"),
  
  # Validation (observed vs predicted probabilities)
  tar_target(fig_predict_vs_obs_meanProba_ms, plot_predicted_vs_observed(
    c(jags.model, jags.model_bis), c(data_jags, data_jags_bis), c(data_model, data_model_bis), 
    method = "mean.proba", "output/fig/ms/validation.png"), 
    format = "file"),
  
  # Trends in disturbance frequency and severity
  tar_target(fig_trend_frequency, plot_trend_disturbance_frequency_ms(
    FUNDIV_plot,FUNDIV_plot_bis, "output/fig/ms/dist_freq_trend.jpg"), format = "file"), 
  tar_target(fig_trend_severity, plot_trend_disturbance_severity_ms(
    FUNDIV_tree, FUNDIV_plot, FUNDIV_plot_bis, "output/fig/ms/dist_sever_trend.jpg"), format = "file"),
  
  # Effect of traits on sensitivity
  tar_target(fig_trait_effect_allsp_ms, plot_trait_effect_ms(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    species, group.in = "all", "output/fig/ms/trait_effect_all.jpg"), format = "file"),
  tar_target(fig_trait_effect_broadleaf_ms, plot_trait_effect_ms(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    species, group.in = "broadleaf", "output/fig/ms/trait_effect_broadleaf.jpg"), format = "file"),
  tar_target(fig_trait_effect_conifer_ms, plot_trait_effect_ms(
    traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
    species, group.in = "conifer", "output/fig/ms/trait_effect_conifer.jpg"), format = "file"),
  
  # Effect of traits on sensitivity using betareg and all iterations of MCMC
  tar_target(fig_trait_effect_allsp_full_ms, plot_trait_effect_full_ms(
    traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, species, group.in = "all", 
    "output/fig/ms/trait_effect_full_all.jpg"), format = "file"),
  tar_target(fig_trait_effect_broadleaf_full_ms, plot_trait_effect_full_ms(
    traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, species, group.in = "broadleaf", 
    "output/fig/ms/trait_effect_full_broadleaf.jpg"), format = "file"),
  tar_target(fig_trait_effect_conifer_full_ms, plot_trait_effect_full_ms(
    traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, species, group.in = "conifer", 
    "output/fig/ms/trait_effect_full_conifer.jpg"), format = "file"),
  
  # Effect of climate on sensitivity
  tar_target(fig_climate_effect_ms, plot_climate_effect_ms(
    gbif_file, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
    file.in = "output/fig/ms/climate_effect.jpg"), format = "file"),
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # -- Step 8 - Export outputs of the model ----
  
  # Export jags objects and correspondence tables
  tar_target(rdata_dominance, export_jags(
    jags.model.in = c(jags.model, jags.model_bis), data_jags.in = c(data_jags, data_jags_bis), 
    file.in = "output/rdata/jags_dominance.Rdata"), format = "file"), 
  tar_target(rdata_stock, export_jags(
    jags.model.in = c(jags.model_stock, jags.model_stock_bis), data_jags.in = c(data_jags_stock, data_jags_stock_bis), 
    file.in = "output/rdata/jags_stock.Rdata"), format = "file")
)