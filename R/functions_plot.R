#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_plot.R  
#' @description R script containing all functions relative to data
#               visualisation
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 1. Generic functions ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Function to get the path of a file, and create directories if they don't exist
#' @param file.in character: path of the file, filename included (ex: "plot/plot.png")
create_dir_if_needed <- function(file.in){

  path.in <- strsplit(file.in, "/")[[1]]
  if(length(path.in) > 1){
    for(i in 1:(length(path.in)-1)){
      if(i == 1) path.in_i <- path.in[i]
      else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
      if(!dir.exists(path.in_i)) dir.create(path.in_i)
    }
  }
}









#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 2. Exploratory plots ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot harvest probability depending on dbh, disturbance and land property
#' @param data_model_scaled Input data for the model centered and scaled
#' @param data_model Input data for the model
#' @param file.in Path and file where to save the plot
plot_harvest_probability <- function(data_model, data_model_scaled, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  # Get property data
  data_property <- fread("data/FrenchNFI_property.csv") %>%
    mutate(property = case_when(propriete %in% c(1:2) ~ "public", 
                                propriete == 3 ~ "private")) %>%
    dplyr::select(plotcode = idp, property) %>%
    filter(plotcode %in% data_model_scaled$plotcode)
  # Format data
  data_priorharvest <- data_model_scaled %>%
    dplyr::select(plotcode, dbh, a, h, D, DS) %>%
    filter(!(a == 1 & D == 1)) %>% # Remove alive trees in disturbed plots
    dplyr::select(-a) %>%
    left_join(data_property, by = "plotcode")
  # Fit model for undisturbed plots
  mod.undist <- glm(h ~ dbh+property, 
                    data = subset(data_priorharvest, D == 0), 
                    family = binomial(link = "logit"))
  # Fit model for disturbed plots
  mod.dist <- glm(h ~ DS*property + dbh, 
                  data = subset(data_priorharvest, D == 1), 
                  family = binomial(link = "logit"))
  # Predict harvest probability
  fit <- data_model_scaled %>% 
    dplyr::select(treecode, plotcode, dbh, DS) %>%
    left_join(data_property, by = "plotcode")
  fit$undisturbed <- predict(mod.undist, newdata = fit, type = "response")
  fit$disturbed <- predict(mod.dist, newdata = fit, type = "response")
  fit <- fit %>%
    dplyr::select(-dbh, -DS) %>%
    left_join((data_model %>% dplyr::select(treecode, dbh, DS, D)), 
              by = "treecode") %>%
    gather(key = disturbance, value = harvest, undisturbed, disturbed) %>%
    filter(!(disturbance == "undisturbed" & D == 1)) %>%
    filter(!(disturbance == "disturbed" & D == 0)) 
  # Plot the fit
  plot.out <- fit %>%
    mutate(status = case_when(disturbance == "undisturbed" ~ "all trees in undisturbed plots", 
                              TRUE ~ "dead trees in disturbed plots")) %>% 
    ggplot(aes(x = dbh, y = harvest, group = interaction(DS, property), 
               linetype = property, colour = DS)) + 
    geom_line(size = 1) + 
    scale_color_gradient(low = "black", high = "red") + 
    facet_wrap(~ status) + 
    theme_bw() + 
    ylab("harvest probability") +
    xlab("DBH (mm)") +
    xlim(0, 1000)
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 17, height = 8, units = "cm", dpi = 600)
  return(file.in)
  
}





#' Plot harvest probability depending on dbh, disturbance and land property
#' @param data_model tree level dataset with the variables necessary for the model
#' @param file.in Path and file where to save the plot
plot_disturbed_trees_per_species <- function(data_model, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  plot.out <- data_model %>%
    mutate(Dnone = ifelse(D == 0, 1, 0)) %>%
    gather(key = "disturbance.type", value = "dist.present", "Dfire", "Dother", "Dstorm", "Dnone") %>%
    filter(dist.present == 1) %>%
    mutate(disturbance.type = gsub("D", "", disturbance.type)) %>%
    group_by(disturbance.type, species) %>%
    summarize(n = n()) %>%
    ggplot(aes(x = species, y = n)) + 
    geom_bar(stat = "identity", colour = "black") + 
    facet_wrap(~ disturbance.type, scales = "free_x") + 
    coord_flip() + 
    theme_bw() + 
    ylab("Number of trees impacted")
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 18, height = 16, units = "cm", dpi = 600)
  return(file.in)
  
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 3. Diagnostic plots ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Plot Markov chain convergence of a rjags object
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param BM_equations dataframe indicating which background mortality to use per species
#' @param disturbance_species_info information on which disturbance affected which species
#' @param dir.in Path where to save the plot
plot_convergence <- function(jags.model, data_jags, BM_equations, disturbance_species_info, 
                             dir.in){
  
  # Initialize output
  out <- c()
  
  # Species included in the simulations
  species.in <- data_jags$species_table$species
  
  # Species for which we have enough disturbed trees
  species.to.keep.in <- (data.frame(ID = disturbance_species_info$species_parameter_to_keep) %>%
                           mutate(species = gsub("\\..+", "", ID)) %>%
                           dplyr::select(-ID) %>%
                           distinct())$species
  
  # Vector storing all parameters extracted
  params <- c()
    
  # Loop on all species
  for(i in 1:length(species.in)){
    
    # Only make the figure if the species have been sufficiently exposed to 
    # to at least one disturbance
    if(species.in[i] %in% species.to.keep.in){
      # Create directories if needed
      file.in.i <- paste0(dir.in, "/fig_convergence_", gsub(" ", "-", species.in[i]), ".png")
      create_dir_if_needed(file.in.i)
      # Add to the output
      out <- c(out, file.in.i)
      
      ## - Identify the columns to extract
      # Get the species code
      species.code.in <- data_jags$species_table$sp[i]
      # Get the list of parameters based on the bm equation of the species
      bm.in <- ifelse((species.in[i] %in% BM_equations$species), 
                      BM_equations$BM_eq[which(BM_equations$species == species.in[i])], 1) 
      if(bm.in == 1) params.in <- c(paste0("b", c(0:5)), paste0("c", c(0:9)))
      if(bm.in == 2) params.in <- c(paste0("b", c(0:9)), paste0("c", c(0:9)))
      if(bm.in == 3) params.in <- c(paste0("b", c(0:11)), paste0("c", c(0:9)))
      if(bm.in == 4) params.in <- c(paste0("b", c(0:7)), paste0("c", c(0:9)))
      if(bm.in == 5) params.in <- c(paste0("b", c(0:11)), paste0("c", c(0:9)))
      if(bm.in == 6) params.in <- c(paste0("b", c(0:13)), paste0("c", c(0:9)))
      # Add the species code to have the complete list of columns to sample
      params.in <- paste0(params.in, "[", species.code.in, "]")
      params <- c(params, params.in)
      
      ## - Make the plot
      plot.i <- ggs(as.mcmc(jags.model)) %>%
        filter(Parameter %in% params.in) %>%
        mutate(Parameter = gsub(paste0("\\[", species.code.in, "\\]"), "", Parameter), 
               Chain = as.factor(Chain), 
               ID = paste(species.in[i], Parameter, sep = "."), 
               value = ifelse(ID %in% disturbance_species_info$species_parameter_to_keep, value, NA_real_)) %>%
        ggplot(aes(x = Iteration, y = value, colour = Chain, group = Chain)) + 
        geom_line() + 
        facet_wrap(~ Parameter, scales = "free") + 
        scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
        theme_bw() + 
        ggtitle(species.in[i])
      
      ## - Save the plot
      ggsave(file.in.i, plot.i, width = 17, height = 12, units = "cm", dpi = 600)
    }
    
    
    
  }
  
  
  # return the name of all the plots made
  return(out)
}



#' Plot harvest probability depending on dbh, disturbance and land property
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param file.in Path and file where to save the plot
plot_parameters_per_species <- function(jags.model, data_jags, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  plot.out <- ggs(as.mcmc(jags.model)) %>%
    filter(Parameter != "deviance") %>%
    mutate(sp = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
           Parameter = gsub("\\[.+\\]", "", Parameter)) %>%
    left_join(data_jags$species_table, by = "sp") %>% 
    dplyr::select(-sp) %>%
    spread(key = species, value = value) %>%
    mutate(.PRIOR = case_when(Parameter %in% paste0("c", c(0, 2, 3, 5)) ~ rnorm(nrow(.), 0, 1/(0.5^2)), 
                              Parameter %in% c("c1", "c4") ~ rnorm(nrow(.), 6, 1), 
                              TRUE ~ NA_real_)) %>%
    gather(key = "species", value = "value", c(unique(data_jags$species_table$species), ".PRIOR")) %>%
    group_by(species, Parameter) %>%
    summarize(mean = mean(value), 
              sd = sd(value)) %>%
    mutate(type = ifelse(species == ".PRIOR", "prior", "species")) %>%
    ggplot(aes(x = species, y = mean, color = type)) + 
    geom_errorbar(aes(ymin = mean - sd, ymax = mean+sd), 
                  width = 0)  + 
    geom_point(size = 0.5) + 
    coord_flip() + 
    scale_color_manual(values = c("red", "black")) +
    facet_wrap(~ Parameter, scales = "free_x", nrow = 1) + 
    ylab("Parameter value") + xlab("") +
    theme(legend.position = "none", 
          axis.line=element_line(), 
          panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank())
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 21, height = 7, units = "cm", dpi = 600)
  return(file.in)
  
}



#' Plot harvest probability depending on dbh, disturbance and land property
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param disturbance_species_info information on which disturbance affected which species
#' @param file.in Path and file where to save the plot
plot_parameters_per_species_full <- function(jags.model, data_jags, disturbance_species_info, 
                                             file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  plot.out <- ggs(as.mcmc(jags.model)) %>%
    mutate(sp = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
           Parameter = gsub("\\[.+\\]", "", Parameter)) %>%
    filter(Parameter %in% paste0("c", c(0:10))) %>%
    left_join(data_jags$species_table, by = "sp") %>% 
    dplyr::select(-sp) %>%
    spread(key = species, value = value) %>%
    mutate(.PRIOR = case_when(Parameter %in% c("c0", "c3", "c6", "c2", "c5", "c8") ~ rnorm(nrow(.), 0, 1), 
                              Parameter %in% c("c1", "c4", "c7") ~ rexp(nrow(.), 0.5), 
                              TRUE ~ NA_real_)) %>%
    gather(key = "species", value = "value", c(unique(data_jags$species_table$species), ".PRIOR")) %>%
    group_by(species, Parameter) %>%
    summarize(mean = mean(value), 
              sd = sd(value)) %>%
    mutate(ID = paste(species, Parameter, sep = ".")) %>%
    filter(ID %in% c(disturbance_species_info$species_parameter_to_keep, 
                     paste0(".PRIOR.c", c(0:8)))) %>%
    mutate(type = ifelse(species == ".PRIOR", "prior", "species")) %>%
    ggplot(aes(x = species, y = mean, color = type)) + 
    geom_errorbar(aes(ymin = mean - sd, ymax = mean+sd), 
                  width = 0)  + 
    geom_point(size = 0.5) + 
    coord_flip() + 
    scale_color_manual(values = c("red", "black")) +
    facet_wrap(~ Parameter, scales = "free_x", nrow = 2) + 
    ylab("Parameter value") + xlab("") +
    theme(legend.position = "none", 
          axis.line=element_line(), 
          panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank())
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 21, height = 20, units = "cm", dpi = 600)
  return(file.in)
  
}

#' Plot estimated vs true parameter value (+ prior) in the case of a model ran with simulated data
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param parameters_sp table containing the true parameters value by species
#' @param disturbance_species_info information on which disturbance affected which species
#' @param file.in Path and file where to save the plot
plot_parameters_true_vs_estimated <- function(jags.model, data_jags, parameters_sp, 
                                              disturbance_species_info, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  plot.out <- ggs(as.mcmc(jags.model)) %>%
    filter(Parameter != "deviance") %>%
    mutate(sp = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
           Parameter = gsub("\\[.+\\]", "", Parameter)) %>%
    left_join(data_jags$species_table, by = "sp") %>% 
    dplyr::select(-sp) %>%
    spread(key = species, value = value) %>%
    mutate(.PRIOR = case_when(Parameter %in% c("c0", "c3", "c6", "c2", "c5", "c8") ~ rnorm(nrow(.), 0, 1), 
                              Parameter %in% c("c1", "c4", "c7") ~ rexp(nrow(.), 0.5), 
                              TRUE ~ NA_real_)) %>%
    gather(key = "species", value = "value", c(unique(data_jags$species_table$species), ".PRIOR")) %>%
    group_by(species, Parameter) %>%
    summarize(mean = mean(value), 
              sd = sd(value)) %>%
    mutate(type = ifelse(species == ".PRIOR", "prior", "species")) %>%
    mutate(ID = paste(species, Parameter, sep = ".")) %>%
    left_join((data_jags$species_table %>%
                 left_join(parameters_sp, by = "sp") %>%
                 gather(key = "Parameter", value = "value.true", paste0("c", c(0:8))) %>%
                 mutate(ID = paste(species, Parameter, sep = ".")) %>%
                 dplyr::select(ID, value.true)), 
              by = "ID") %>%
    filter(ID %in% c(disturbance_species_info$species_parameter_to_keep, 
                     paste0(".PRIOR.c", c(0:8)))) %>%
    ggplot(aes(x = species, y = mean, color = type)) + 
    geom_errorbar(aes(ymin = mean - sd, ymax = mean+sd), 
                  width = 0)  + 
    geom_point(size = 1) + 
    geom_point(aes(y = value.true), inherit.aes = TRUE, color = "blue", size = 1) +
    coord_flip() + 
    scale_color_manual(values = c("red", "black")) +
    facet_wrap(~ Parameter, nrow = 1) + 
    ylab("Parameter value") + xlab("") +
    theme(legend.position = "none", 
          axis.line=element_line(), 
          panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank())
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 21, height = 9, units = "cm", dpi = 600)
  return(file.in)
  
}



#' Plot estimated vs true parameter value (+ prior) in the case of a model ran with simulated data
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param parameters_sp table containing the true parameters value by species
#' @param disturbance_species_info information on which disturbance affected which species
#' @param file.in Path and file where to save the plot
plot_parameters_true_vs_estimated2 <- function(jags.model, data_jags, parameters_sp, 
                                              disturbance_species_info, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  plot.out <- ggs(as.mcmc(jags.model)) %>%
    filter(Parameter != "deviance") %>%
    mutate(sp = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
           Parameter = gsub("\\[.+\\]", "", Parameter)) %>%
    left_join(data_jags$species_table, by = "sp") %>% 
    dplyr::select(-sp) %>%
    spread(key = species, value = value) %>%
    mutate(.PRIOR = case_when(Parameter %in% c("c0", "c3", "c6", "c2", "c5", "c8", "c9", "c10") ~ rnorm(nrow(.), 0, 1), 
                              Parameter %in% c("c1", "c4", "c7") ~ rexp(nrow(.), 0.5), 
                              TRUE ~ NA_real_)) %>%
    gather(key = "species", value = "value", c(unique(data_jags$species_table$species), ".PRIOR")) %>%
    group_by(species, Parameter) %>%
    summarize(mean = mean(value), 
              sd = sd(value)) %>%
    mutate(type = ifelse(species == ".PRIOR", "prior", "species")) %>%
    mutate(ID = paste(species, Parameter, sep = ".")) %>%
    left_join((data_jags$species_table %>%
                 left_join(parameters_sp, by = "sp") %>%
                 gather(key = "Parameter", value = "value.true", 
                        colnames(.)[which(colnames(.) %in% paste0("c", c(0:11)))]) %>%
                 mutate(ID = paste(species, Parameter, sep = ".")) %>%
                 dplyr::select(ID, value.true)), 
              by = "ID") %>%
    filter(ID %in% c(disturbance_species_info$species_parameter_to_keep, 
                     paste0(".PRIOR.c", c(0:11)))) %>%
    filter(type != "prior") %>%
    mutate(Parameter = factor(
      Parameter, levels = paste0("c", c(0:11))[which(paste0("c", c(0:11)) %in% unique(.$Parameter))])) %>%
    ggplot(aes(x = value.true, y = mean, fill = species)) + 
    geom_errorbar(aes(ymin = mean - sd, ymax = mean+sd), 
                  width = 0)  + 
    geom_point(size = 1, shape = 21, color = "black") + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~ Parameter, scales = "free") + 
    ylab("Estimated parameter value") + xlab("True parameter value") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.title = element_blank(),
          legend.text = element_text(size = 7, face = "italic"),
          legend.key = element_rect(fill = alpha("white", 0.0))) 
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 23, height = 14, units = "cm", dpi = 600)
  return(file.in)
  
}







#' Plot death probability predicted by the model
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param data_model_scaled Tree data formatted for the IPM and scaled 
#' @param dir.in Path where to save the plot
plot_prediction_full <- function(jags.model, data_jags, data_model_scaled, data_model, 
                                 disturbance_species_info, dir.in){
  
  # Initialize output
  out <- paste0(dir.in, c("/storm", "/other", "/fire"), ".png")
  
  # Create the diretories that are needed
  for(i in 1:length(out)) create_dir_if_needed(out[i])
  
  # Identify parameters per species
  param_per_species <- ggs(as.mcmc(jags.model)) %>%
    filter(Parameter != "deviance") %>%
    mutate(sp = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
           Parameter = gsub("\\[.+\\]", "", Parameter)) %>%
    left_join(data_jags$species_table, by = "sp") %>% 
    group_by(species, Parameter) %>%
    summarize(mean = mean(value)) %>%
    spread(key = "Parameter", value = "mean")
  
  # Model to scale dbh
  scale_dbh <- lm(dbh.scaled ~ dbh, 
                  data = data.frame(dbh = data_model$dbh, 
                                    dbh.scaled = data_model_scaled$dbh))
  
  # - Preformat the data before plotting
  data_for_plot <- expand.grid(species = unique(param_per_species$species), 
                               dbh = c(100:1200), 
                               Intensity = c(0.2, 0.4, 0.6, 0.8, 1)) %>%
    mutate(dbh.scaled = predict(scale_dbh, newdata = .)) %>%
    left_join(param_per_species, by = "species") %>%
    gather(key = "Parameter", value = "value", paste0("c", c(0:8))) %>%
    mutate(ID = paste(species, Parameter, sep = "."), 
           value = ifelse(ID %in% disturbance_species_info$species_parameter_to_keep, value, NA_real_)) %>%
    dplyr::select(-ID) %>%
    spread(key = Parameter, value = value) %>%
    mutate(pdstorm = plogis(c0 + c1*Intensity*dbh.scaled^c2), 
           pdother = plogis(c3 + c4*Intensity*dbh.scaled^c5), 
           pdfire = plogis(c6 + c7*Intensity*dbh.scaled^c8)) 
  
  # - Plot for storm
  plot.storm <- data_for_plot %>%
    ggplot(aes(x = dbh, y = pdstorm, group = Intensity, color = Intensity)) + 
    geom_line() + 
    scale_color_gradient(low = "blue", high = "black") + 
    facet_wrap(~ species) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank()) + 
    ylab("Probability to die from a storm") + 
    xlab("dbh (mm)")
  
  # - Plot for other
  plot.other <- data_for_plot %>%
    ggplot(aes(x = dbh, y = pdother, group = Intensity, color = Intensity)) + 
    geom_line() + 
    scale_color_gradient(low = "green", high = "black") + 
    facet_wrap(~ species) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank()) + 
    ylab("Probability to die from an other disturbance") + 
    xlab("dbh (mm)")
  
  # - Plot for fire
  plot.fire <- data_for_plot %>%
    ggplot(aes(x = dbh, y = pdfire, group = Intensity, color = Intensity)) + 
    geom_line() + 
    scale_color_gradient(low = "orange", high = "black") + 
    facet_wrap(~ species) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank()) + 
    ylab("Probability to die from a fire") + 
    xlab("dbh (mm)")
  
  
  # - Save the three plots
  ggsave(out[1], plot.storm, width = 25, height = 20, units = "cm", dpi = 600)
  ggsave(out[2], plot.other, width = 25, height = 20, units = "cm", dpi = 600)
  ggsave(out[3], plot.fire, width = 25, height = 20, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(out)
}



#' Plot the relation between disturbance severity (observed) and intensity (estimated)
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param file.in Path and file where to save the plot
plot_intensity_vs_severity <- function(jags.model, data_jags, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - make the plot
  plot.out <- data.frame(plot = data_jags$data_jags$plot, 
                         Ifire = data_jags$data_jags$Dfire, 
                         Istorm = data_jags$data_jags$Dstorm, 
                         Iother = data_jags$data_jags$Dother) %>%
    distinct() %>%
    gather(key = "variable", value = "value", "Ifire", "Istorm", "Iother") %>%
    filter(value == 1) %>%
    mutate(Parameter = paste0(variable, "[", plot, "]")) %>%
    left_join((ggs(as.mcmc(jags.model)) %>%
                 group_by(Parameter) %>%
                 summarize(intensity.mean = mean(value, na.rm = T), 
                           intensity.sd = sd(value, na.rm = T))), 
              by = "Parameter") %>%
    left_join((data.frame(plot = data_jags$data_jags$plot, 
                          state = data_jags$data_jags$state) %>%
                 mutate(dead = ifelse(state == 3, 0, 1)) %>%
                 group_by(plot) %>%
                 summarize(severity = sum(dead)/n())), 
              by = "plot") %>%
    mutate(disturbance.type = gsub("I", "", variable)) %>%
    dplyr::select(plot, disturbance.type, severity, intensity.mean, intensity.sd) %>%
    filter(!is.na(intensity.mean)) %>%
    ggplot(aes(x = severity, y = intensity.mean)) + 
    geom_point(aes(size = 1/intensity.sd), shape = 21, 
               alpha = 0.3, fill = "black", color = "black") + 
    theme_bw() + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_smooth(se = FALSE, method = "loess", color = "red") + 
    facet_wrap(~ disturbance.type) + 
    xlab("Observed disturbance severity") + 
    ylab("Estimated disturbance intensity") + 
    theme(strip.background = element_blank(), 
          legend.position = "none")
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 23, height = 8, units = "cm", dpi = 600)
  return(file.in)
  
}

