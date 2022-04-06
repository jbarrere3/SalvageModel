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



#' Plot Markov chain convergence of a rjags object
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param BM_equations dataframe indicating which background mortality to use per species
#' @param dir.in Path where to save the plot
plot_convergence <- function(jags.model, data_jags, BM_equations, dir.in){
  
  # Initialize output
  out <- c()
  
  # Species included in the simulations
  species.in <- data_jags$species_table$species
  
  # Vector storing all parameters extracted
  params <- c()
    
  # Loop on all species
  for(i in 1:length(species.in)){
    
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
    if(bm.in == 1) params.in <- c(paste0("b", c(0:5)), paste0("c", c(0:6)))
    if(bm.in == 2) params.in <- c(paste0("b", c(0:9)), paste0("c", c(0:6)))
    if(bm.in == 3) params.in <- c(paste0("b", c(0:11)), paste0("c", c(0:6)))
    if(bm.in == 4) params.in <- c(paste0("b", c(0:7)), paste0("c", c(0:6)))
    if(bm.in == 5) params.in <- c(paste0("b", c(0:11)), paste0("c", c(0:6)))
    if(bm.in == 6) params.in <- c(paste0("b", c(0:13)), paste0("c", c(0:6)))
    # Add the species code to have the complete list of columns to sample
    params.in <- paste0(params.in, "[", species.code.in, "]")
    params <- c(params, params.in)
    
    ## - Make the plot
    plot.i <- ggs(as.mcmc(jags.model)) %>%
      filter(Parameter %in% params.in) %>%
      mutate(Parameter = gsub(paste0("\\[", species.code.in, "\\]"), "", Parameter), 
             Chain = as.factor(Chain)) %>%
      ggplot(aes(x = Iteration, y = value, colour = Chain, group = Chain)) + 
      geom_line() + 
      facet_wrap(~ Parameter, scales = "free") + 
      scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
      theme_bw() + 
      ggtitle(species.in[i])
    
    ## - Save the plot
    ggsave(file.in.i, plot.i, width = 17, height = 12, units = "cm", dpi = 600)
    
  }
  
  ## - Plot of all the other parameters
  ## - Make the plot
  plot.other <- ggs(as.mcmc(jags.model)) %>%
    filter(!(Parameter %in% params)) %>%
    mutate(Chain = as.factor(Chain)) %>%
    ggplot(aes(x = Iteration, y = value, colour = Chain, group = Chain)) + 
    geom_line() + 
    facet_wrap(~ Parameter, scales = "free") + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    theme_bw()
  # Create directories if needed
  file.other <- paste0(dir.in, "/fig_other_param.png")
  create_dir_if_needed(file.other)
  # Add to the output
  out <- c(out, file.other)
  ## - Save the plot
  ggsave(file.other, plot.other, width = 17, height = 8, units = "cm", dpi = 600)
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