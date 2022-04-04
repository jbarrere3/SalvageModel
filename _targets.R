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
  # - Step 1 - Download and format raw data
  
  # Download files
  tar_target(files, get_FrenchNFI(), format = "file"), 
  
  # Load raw data
  tar_target(FrenchNFI_tree_raw, fread(files[1])), 
  tar_target(FrenchNFI_plot_raw, fread(files[7])), 
  tar_target(FrenchNFI_species, fread("data/FrenchNFI_species.csv")), 
  tar_target(Climate, fread("data/Climate.csv")), 
  tar_target(Disturbance, fread("data/Disturbance.csv")), 
  
  # Format raw data 
  tar_target(FrenchNFI, format_FrenchNFI_raw(FrenchNFI_tree_raw, FrenchNFI_plot_raw)),
  
  # Format to FUNDIV template
  tar_target(FUNDIV_tree_FR, format_FrenchNFI_tree_to_FUNDIV(FrenchNFI, FrenchNFI_species)), 
  tar_target(FUNDIV_plot_FR, format_FrenchNFI_plot_to_FUNDIV(FrenchNFI, FUNDIV_tree_FR)), 
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 2 - Prepare data for the model
  
  # Keep and scale variables relevant to the model
  tar_target(data_model, format_data_model(FUNDIV_tree_FR, FUNDIV_plot_FR, Disturbance, Climate)), 
  tar_target(data_model_scaled, scale_data_model(data_model)), 
  tar_target(data_model_scaled_BM_D, data_model_scaled %>% mutate(dbh = dbh - min(dbh) + 0.01)), 
  
  # Get parameters for harvest conditional probabilities
  tar_target(param_harvest_proba, get_param_harvest_proba(data_model_scaled, Dj_latent = FALSE)),
  tar_target(param_harvest_proba_Djlatent, get_param_harvest_proba(data_model_scaled, Dj_latent = TRUE)),
  
  # Generate input data for the model with disturbance latent and parameters a0 and a1 given
  tar_target(data_simulated, simulate_status(subset(data_model_scaled, species == "Abies alba"), 
                                             param_harvest_proba_Djlatent, Dj_latent = TRUE)),
  tar_target(data_jags_simulated, generate_data_jags_from_simulated(data_simulated, Dj_latent = TRUE)),
  
  # Generate input data for the model with disturbance given
  tar_target(data_simulated_D, simulate_status(subset(data_model_scaled, species == "Abies alba"), 
                                             param_harvest_proba, Dj_latent = FALSE)),
  tar_target(data_jags_simulated_D, generate_data_jags_from_simulated(data_simulated_D, Dj_latent = FALSE)),
  
  # Generate input data for the model with disturbance given for new BM functions
  # - Background mortality equation 1
  tar_target(data_simulated_BM1_D, simulate_status_BM_D(subset(data_model_scaled_BM_D, species == "Abies alba"), 
                                               param_harvest_proba, BM.equation = 1)),
  tar_target(data_jags_simulated_BM1_D, generate_data_jags_from_simulated(data_simulated_BM1_D, Dj_latent = FALSE)),
  # - Background mortality equation 2
  tar_target(data_simulated_BM2_D, simulate_status_BM_D(subset(data_model_scaled_BM_D, species == "Abies alba"), 
                                                        param_harvest_proba, BM.equation = 2)),
  tar_target(data_jags_simulated_BM2_D, generate_data_jags_from_simulated(data_simulated_BM2_D, Dj_latent = FALSE)),
  # - Background mortality equation 3
  tar_target(data_simulated_BM3_D, simulate_status_BM_D(subset(data_model_scaled_BM_D, species == "Abies alba"), 
                                                        param_harvest_proba, BM.equation = 3)),
  tar_target(data_jags_simulated_BM3_D, generate_data_jags_from_simulated(data_simulated_BM3_D, Dj_latent = FALSE)),
  # - Background mortality equation 4
  tar_target(data_simulated_BM4_D, simulate_status_BM_D(subset(data_model_scaled_BM_D, species == "Abies alba"), 
                                                        param_harvest_proba, BM.equation = 4)),
  tar_target(data_jags_simulated_BM4_D, generate_data_jags_from_simulated(data_simulated_BM4_D, Dj_latent = FALSE)),
  # - Background mortality equation 5
  tar_target(data_simulated_BM5_D, simulate_status_BM_D(subset(data_model_scaled_BM_D, species == "Abies alba"), 
                                                        param_harvest_proba, BM.equation = 5)),
  tar_target(data_jags_simulated_BM5_D, generate_data_jags_from_simulated(data_simulated_BM5_D, Dj_latent = FALSE)),
  # - Background mortality equation 6
  tar_target(data_simulated_BM6_D, simulate_status_BM_D(subset(data_model_scaled_BM_D, species == "Abies alba"), 
                                                        param_harvest_proba, BM.equation = 6)),
  tar_target(data_jags_simulated_BM6_D, generate_data_jags_from_simulated(data_simulated_BM6_D, Dj_latent = FALSE)),
  
  # Generate input data for jags with real data 
  tar_target(param_disturbance_proba, get_param_disturbance_proba(data_model_scaled)), 
  # -- With silver fir
  tar_target(data_jags_A.alba, generate_data_jags(subset(data_model_scaled, species == "Abies alba"),
                                                  param_harvest_proba)), 
  tar_target(data_jags_A.alba_D, generate_data_jags_D(subset(data_model_scaled, species == "Abies alba"),
                                                  param_harvest_proba)), 
  
  # -- With hornbeam
  tar_target(data_jags_C.betulus, generate_data_jags(subset(data_model_scaled, species == "Carpinus betulus"),
                                                  param_harvest_proba)), 
  tar_target(data_jags_C.betulus_D, generate_data_jags_D(subset(data_model_scaled, species == "Carpinus betulus"),
                                                      param_harvest_proba)), 
  
  # -- With oak
  tar_target(data_jags_Q.robur, generate_data_jags(subset(data_model_scaled, species == "Quercus robur"),
                                                     param_harvest_proba)), 
  tar_target(data_jags_Q.robur_D, generate_data_jags_D(subset(data_model_scaled, species == "Quercus robur"),
                                                         param_harvest_proba)),
  # -- With spruce
  tar_target(data_jags_P.abies, generate_data_jags(subset(data_model_scaled, 
                                                          (species == "Picea abies") & (comp < 7)),
                                                   param_harvest_proba)), 
  tar_target(data_jags_P.abies_D, generate_data_jags_D(subset(data_model_scaled, species == "Picea abies"),
                                                         param_harvest_proba)),
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 3 - Fit the jags model
  
  # Fit jags models
  tar_target(jags.model_D, fit_mortality_D(data_jags_simulated_D, n.chains = 3, n.iter = 5000, 
                                       n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_BM1_D, fit_mortality_BM1_D(data_jags_simulated_BM1_D, n.chains = 3, n.iter = 5000, 
                                           n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_BM2_D, fit_mortality_BM2_D(data_jags_simulated_BM2_D, n.chains = 3, n.iter = 5000, 
                                                   n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_BM3_D, fit_mortality_BM3_D(data_jags_simulated_BM3_D, n.chains = 3, n.iter = 5000, 
                                                   n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_BM4_D, fit_mortality_BM4_D(data_jags_simulated_BM4_D, n.chains = 3, n.iter = 5000, 
                                                   n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_BM5_D, fit_mortality_BM5_D(data_jags_simulated_BM5_D, n.chains = 3, n.iter = 5000, 
                                                   n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_BM6_D, fit_mortality_BM6_D(data_jags_simulated_BM6_D, n.chains = 3, n.iter = 5000, 
                                                   n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_A.alba_D, fit_mortality_D(data_jags_A.alba_D, n.chains = 3, 
                                                  n.iter = 5000, n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_C.betulus_D, fit_mortality_D(data_jags_C.betulus_D, n.chains = 3, 
                                                     n.iter = 5000, n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_Q.robur_D, fit_mortality_D(data_jags_Q.robur_D, n.chains = 3, 
                                                   n.iter = 5000, n.burn = 1000, n.thin = 1)),
  tar_target(jags.model_P.abies_D, fit_mortality_D(data_jags_P.abies_D, n.chains = 3, 
                                                   n.iter = 5000, n.burn = 1000, n.thin = 1)),
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 4 - Plot the outputs of the model
  
  # Convergence
  tar_target(fig_jags.model_chains_D, 
             plot_convergence(jags.model_D, file.in = "fig/simulated_data/fig_convergence_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_chains_BM1_D, 
             plot_convergence(jags.model_BM1_D, file.in = "fig/simulated_data/fig_convergence_BM1_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_chains_BM2_D, 
             plot_convergence(jags.model_BM2_D, file.in = "fig/simulated_data/fig_convergence_BM2_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_chains_BM3_D, 
             plot_convergence(jags.model_BM3_D, file.in = "fig/simulated_data/fig_convergence_BM3_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_chains_BM4_D, 
             plot_convergence(jags.model_BM4_D, file.in = "fig/simulated_data/fig_convergence_BM4_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_chains_BM5_D, 
             plot_convergence(jags.model_BM5_D, file.in = "fig/simulated_data/fig_convergence_BM5_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_chains_BM6_D, 
             plot_convergence(jags.model_BM6_D, file.in = "fig/simulated_data/fig_convergence_BM6_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_A.alba_chains_D, 
             plot_convergence(jags.model_A.alba_D, file.in = "fig/real_data/fig_convergence_Aalba_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_C.betulus_chains_D, 
             plot_convergence(jags.model_C.betulus_D, file.in = "fig/real_data/fig_convergence_Cbetulus_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_Q.robur_chains_D, 
             plot_convergence(jags.model_Q.robur_D, file.in = "fig/real_data/fig_convergence_Qrobur_D.png"), 
             format = "file"), 
  tar_target(fig_jags.model_P.abies_chains_D, 
             plot_convergence(jags.model_P.abies_D, file.in = "fig/real_data/fig_convergence_P.abies_D.png"), 
             format = "file"), 
  
  # Parameters value
  tar_target(fig_fitted_vs_true_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -3, b2 = 3, b3 = -3, 
                                                              b4 = 3, c0 = 0, c1 = 3, c2 = -3)), 
                                            jags.model_D, "fig/simulated_data/fitted_vs_true_D.png"), 
             format = "file"), 
  tar_target(fig_fitted_vs_true_BM1_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -2, b2 = 1, b3 = 1.5, 
                                                              b4 = 3, b5 = -1, c0 = 0, c1 = 3, c2 = -3)), 
                                            jags.model_BM1_D, "fig/simulated_data/fitted_vs_true_BM1_D.png"), 
             format = "file"), 
  tar_target(fig_fitted_vs_true_BM2_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -2, b2 = 1, b3 = 1.5, b4 = 3,
                                                              b5 = -1, b6 = 0, b7 = -2, b8 = 1, b9 = 1.5, 
                                                              c0 = 0, c1 = 3, c2 = -3)), 
                                            jags.model_BM2_D, "fig/simulated_data/fitted_vs_true_BM2_D.png"), 
             format = "file"), 
  tar_target(fig_fitted_vs_true_BM3_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -2, b2 = 1, b3 = 1.5, b4 = 3,
                                                              b5 = -1, b6 = 0, b7 = -2, b8 = 1, b9 = 1.5, 
                                                              b10 = 3, b11 = -1, c0 = 0, c1 = 3, c2 = -3)), 
                                            jags.model_BM3_D, "fig/simulated_data/fitted_vs_true_BM3_D.png"), 
             format = "file"), 
  tar_target(fig_fitted_vs_true_BM4_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -2, b2 = 1, b3 = 1.5, b4 = 3,
                                                              b5 = -1, b6 = 0, b7 = -2, c0 = 0, c1 = 3, 
                                                              c2 = -3)), 
                                            jags.model_BM4_D, "fig/simulated_data/fitted_vs_true_BM4_D.png"), 
             format = "file"), 
  tar_target(fig_fitted_vs_true_BM5_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -2, b2 = 1, b3 = 1.5, b4 = 3,
                                                              b5 = -1, b6 = 0, b7 = -2, b8 = 1, b9 = 1.5, 
                                                              b10 = 3, b11 = -1, c0 = 0, 
                                                              c1 = 3, c2 = -3)), 
                                            jags.model_BM5_D, "fig/simulated_data/fitted_vs_true_BM5_D.png"), 
             format = "file"), 
  tar_target(fig_fitted_vs_true_BM6_D, 
             plot_fitted_vs_true_parameters(list(param = list(b0 = 0, b1 = -2, b2 = 1, b3 = 1.5, b4 = 3,
                                                              b5 = -1, b6 = 0, b7 = -2, b8 = 1, b9 = 1.5, 
                                                              b10 = 3, b11 = -1, b12 = -1, b13 = 0, c0 = 0, 
                                                              c1 = 3, c2 = -3)), 
                                            jags.model_BM6_D, "fig/simulated_data/fitted_vs_true_BM6_D.png"), 
             format = "file"), 
  tar_target(fig_param_species_D, 
             plot_parameters_species(list.in = list(A.alba = jags.model_A.alba_D, 
                                                    C.betulus = jags.model_C.betulus_D, 
                                                    Q.robur = jags.model_Q.robur_D, 
                                                    P.abies = jags.model_P.abies_D), 
                                     file.in = "fig/real_data/param_per_species_D.png"), 
             format = "file"),
  
  # Exploratory plots
  tar_target(fig_harvest_proba, plot_harvest_probability(data_model, data_model_scaled, 
                                                         "fig/exploratory/harvest_proba.png"), 
             format = "file"), 
  
  
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # - Step 5 - New version of the model
  tar_target(BM_equations, fread("data/BM_equations.csv")), 
  tar_target(data_model2, format_data_model2(FUNDIV_tree_FR, FUNDIV_plot_FR, Climate, BM_equations)), 
  tar_target(data_model2_scaled, scale_data_model2(data_model2, var = c("dbh", "comp", "sgdd", "wai", "DS"))), 
  tar_target(param_harvest_proba2, get_param_harvest_proba2(data_model2_scaled, Dj_latent = FALSE)),
  tar_target(data_jags2, generate_data_jags2(data_model2_scaled, param_harvest_proba2)),
  tar_target(jags.model2, fit_mortality2(data_jags2$data_jags, n.chains = 3, n.iter = 1500, n.burn = 500, n.thin = 1)), 
  tar_target(fig_convergence2, plot_convergence2(jags.model2, data_jags2, BM_equations, 
                                                 "fig/real_data/multispecies_model/convergence"), 
             format = "file")
  
  
)