#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Prepare data for model with real data ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Get family and order for each species present in the dataset
#' @param FUNDIV_tree Tree table with a column entitled "species"
get_species_info <- function(FUNDIV_tree){
  
  # Extract genus and family of all species present
  extraction.original <- cbind(data.frame(species.original = unique(FUNDIV_tree$species)), 
                               TPL(splist = unique(FUNDIV_tree$species)))
  
  # Table linking family to order
  data(vascular.families)
  
  # Format the output
  out <- extraction.original %>%
    dplyr::select(species = species.original, 
                  genus = New.Genus) %>% 
    left_join((extraction.original %>%
                 dplyr::select(genus = New.Genus, family = Family) %>%
                 filter(!is.na(family)) %>%
                 filter(family != "") %>%
                 distinct()), 
              by = "genus") %>%
    left_join((vascular.families %>%
                 dplyr::select(family = Family, group = Group)), 
              by = "family") %>%
    filter(!is.na(genus)) %>%
    mutate(group = ifelse(is.na(group), group, paste(toupper(substr(group, 1, 1)), 
                                                     substr(group, 2, nchar(group)), sep="")))
  
  # return output
  return(out)
  
}
#' Format the data before fitting the reference mortality model
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plot FUNDIV plot table
#' @param Climate dataframe containing climatic data per plot
#' @param species dataframe containing species information
format_data_model <- function(FUNDIV_tree, FUNDIV_plot, Climate, species){
  
  ## - Begin formatting final dataset
  data <- FUNDIV_tree %>%
    # Add plot level data 
    left_join(Climate, by = "plotcode") %>% 
    left_join((FUNDIV_plot %>% 
                 dplyr::select(plotcode, time = yearsbetweensurveys, disturbance.nature)), 
              by = "plotcode") %>% 
    # Only keep disturbed plots since the model is only fitted with disturbed plots
    filter(disturbance.nature != "none") %>%
    # Remove ingrowth
    filter(treestatus != 1) %>%
    # Remove recruited trees that died
    filter(dbh1 > 100) %>%
    # Remove plots with null severity
    mutate(dead = ifelse(treestatus == 2, 0, 1)) %>%
    group_by(plotcode) %>%
    mutate(dead.sum = sum(dead)) %>%
    filter(dead.sum > 0) %>%
    ungroup() %>%
    dplyr::select(-dead, -dead.sum) %>%
    # Add species information
    left_join((species %>% dplyr::select(species, group)), by = "species")
  
  ## - Aggregated species that are not represented well enough
  species.aggregated <- data %>%
    group_by(species, group, plotcode, disturbance.nature) %>% 
    summarize(n.indiv.per.plotcode = n()) %>%
    ungroup() %>%
    group_by(species, group, disturbance.nature) %>%
    summarize(n.per.species = sum(n.indiv.per.plotcode), 
              n.plot.per.species = n()) %>%
    ungroup() %>%
    mutate(enough.individuals = ifelse((n.per.species > 150 & n.plot.per.species > 15), 1, 0), 
           species.ag = case_when((enough.individuals == 0 & group == "Angiosperms") ~ "Other broadleaf", 
                                  (enough.individuals == 0 & group == "Gymnosperms") ~ "Other conifer",
                                  enough.individuals == 1  ~ species)) %>%
    dplyr::select(group, species, disturbance.nature, species.ag) %>%
    filter(!(species %in% c("", " "))) %>%
    spread(key = "disturbance.nature", value = "species.ag") %>%
    # Replace NA by other
    ungroup() %>%
    gather(key = "disturbance.nature", value = "species.ag", unique(data$disturbance.nature)) %>%
    mutate(species.ag = case_when((is.na(species.ag) & group == "Angiosperms") ~ "Other broadleaf", 
                                  (is.na(species.ag) & group == "Gymnosperms") ~ "Other conifer", 
                                  TRUE ~ species.ag)) %>%
    group_by(species, disturbance.nature) %>%
    mutate(species.ag = ifelse(length(grep("\\ ", species.ag)) == 0, paste0(species.ag, " sp"), species.ag)) %>%
    ungroup() %>%
    drop_na() %>%
    spread(key = "disturbance.nature", value = "species.ag") %>%
    dplyr::select(-group)
  
  ## - Initialize the list output and loop on all types of disturbance
  out <- list()
  disturbances.in <- unique(data$disturbance.nature)
  for(i in 1:length(disturbances.in)){
    # Output dataset for disturbance i
    data.i <- data %>%
      # Filter by disturbance nature
      filter(disturbance.nature == disturbances.in[i]) %>%
      # Join species aggregated table
      left_join((species.aggregated %>% 
                   dplyr::select("species", "species.ag" = disturbances.in[i])), 
                by = "species") %>%
      # Compute basal area of competitors
      group_by(plotcode) %>%
      mutate(comp=sum(ba_ha1) - ba_ha1) %>%
      # Compute stand quadratic diameter
      filter(!is.na(weight1)) %>%
      mutate(dbh2.W = (dbh1^2)*weight1, 
             dqm = sqrt(sum(dbh2.W, na.rm = TRUE)/sum(weight1, na.rm  = TRUE)), 
             stock = sum(ba_ha1)) %>%
      ungroup() %>%
      # compute logratio dbh/dqm
      mutate(log_dbh_dqm = log(dbh1/dqm)) %>%
      # Reformat status columns
      mutate(h = case_when(treestatus == 3 ~ 1, T ~ 0), 
             d = case_when(treestatus %in% c(4, 5) ~ 1, T ~ 0), 
             a = case_when(treestatus == 2 ~ 1, T ~ 0)) %>%
      # Reformat disturbance columns
      filter(!is.na(sgdd) & !is.na(wai) & !is.na(comp) & !is.na(species.ag)) %>%
      # Select variables and remove NAs
      dplyr::select(plotcode, treecode, country, species = species.ag, time, h, d, a, dbh = dbh1, comp, 
                    sgdd, wai, dqm, log_dbh_dqm, stock)
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
    
  }
  
  # return output list
  return(out)
}



#' Center and Scale variables before fitting one fo the models that includes dqm
#' @param data_model Tree data formated for the IPM. 
#' @param var Variables to scale (character vector)
scale_data_model <- function(data_model,  var = c("dbh", "dqm", "comp", "sgdd", "wai", "log_dbh_dqm", "stock")){
  id_var <- which(colnames(data_model) %in% var) 
  out <- data_model
  unscaled = as.matrix(out[, id_var])
  scaled = scale(unscaled)
  out[, id_var] <- as.data.frame(scaled)
  # dbh and dqm has a power effect, so need to make them positive
  if("dbh" %in% var) if(min(out$dbh) < 0) out$dbh <- out$dbh - min(out$dbh) + 0.01 
  if("dqm" %in% var) if(min(out$dqm) < 0) out$dqm <- out$dqm - min(out$dqm) + 0.01 
  colnames(out) <- colnames(data_model)
  return(out)
}



#' generate data for the jags model
#' @param data_model dataset where the tree status (dead, alive or harvested) is not simulated
generate_data_jags <- function(data_model){
  
  # Initialize output list
  out <- list()
  
  # Loop on all disturbances
  disturbances.in <- names(data_model)
  
  for(i in 1:length(disturbances.in)){
    # Store dataframe for disturbance i in object data.i 
    eval(parse(text = paste0("data.i <- data_model$", disturbances.in[i])))
    # Scale the dataset
    data.i_scaled <- scale_data_model(data.i)
    # Replace plotcode by a number
    plotcode.table.i <- data.frame(plotcode = unique(data.i_scaled$plotcode), 
                                   plot = c(1:length(unique(data.i_scaled$plotcode))))
    # Replace species by a number 
    species.table.i <- data.frame(species = unique(data.i_scaled$species), 
                                  sp = c(1:length(unique(data.i_scaled$species))))
    # Replace country by a number
    country.table.i <- data.frame(country = unique(data.i_scaled$country), 
                                  co = c(1:length(unique(data.i_scaled$country))))
    # Format final data
    data.out.i <- data.i_scaled %>%
      # Determine if tree is dead or not
      mutate(d = ifelse(a == 1, 0, 1)) %>%
      # Add plot code
      left_join(plotcode.table.i, by = "plotcode") %>%
      # Add species code
      left_join(species.table.i, by = "species") %>%
      # Add country code
      left_join(country.table.i, by = "country")
    
    # Create the data list
    data_jags.i <- list(
      Ntrees = NROW(data.out.i),
      Nplot = NROW(plotcode.table.i), 
      Nspecies = NROW(species.table.i),
      Ncountry = NROW(country.table.i),
      plot = data.out.i$plot, 
      sp = data.out.i$sp,
      co = data.out.i$co,
      d = data.out.i$d,
      time = data.out.i$time,
      dbh = data.out.i$dbh,
      logratio = data.out.i$log_dbh_dqm)
    # Final list for disturbance i
    out.i <- list(data_jags = data_jags.i, 
                  plotcode.table = plotcode.table.i, 
                  species.table = species.table.i, 
                  country.table = country.table.i)
    # Add list i to the global list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- out.i")))
    
  }
  
  return(out)
}


#' generate data for the jags model with stocking as explanatory variable
#' @param data_model dataset where the tree status (dead, alive or harvested) is not simulated
generate_data_jags_stock <- function(data_model){
  
  # Initialize output list
  out <- list()
  
  # Loop on all disturbances
  disturbances.in <- names(data_model)
  
  for(i in 1:length(disturbances.in)){
    # Store dataframe for disturbance i in object data.i 
    eval(parse(text = paste0("data.i <- data_model$", disturbances.in[i])))
    # Scale the dataset
    data.i_scaled <- scale_data_model(data.i)
    # Replace plotcode by a number
    plotcode.table.i <- data.frame(plotcode = unique(data.i_scaled$plotcode), 
                                   plot = c(1:length(unique(data.i_scaled$plotcode))))
    # Replace species by a number 
    species.table.i <- data.frame(species = unique(data.i_scaled$species), 
                                  sp = c(1:length(unique(data.i_scaled$species))))
    # Replace country by a number
    country.table.i <- data.frame(country = unique(data.i_scaled$country), 
                                  co = c(1:length(unique(data.i_scaled$country))))
    # Format final data
    data.out.i <- data.i_scaled %>%
      # Determine if tree is dead or not
      mutate(d = ifelse(a == 1, 0, 1)) %>%
      # Add plot code
      left_join(plotcode.table.i, by = "plotcode") %>%
      # Add species code
      left_join(species.table.i, by = "species") %>%
      # Add country code
      left_join(country.table.i, by = "country")
    
    # Create the data list
    data_jags.i <- list(
      Ntrees = NROW(data.out.i),
      Nplot = NROW(plotcode.table.i), 
      Nspecies = NROW(species.table.i),
      Ncountry = NROW(country.table.i),
      plot = data.out.i$plot, 
      sp = data.out.i$sp,
      co = data.out.i$co,
      d = data.out.i$d,
      time = data.out.i$time,
      dbh = data.out.i$dbh,
      stock = data.out.i$stock)
    # Final list for disturbance i
    out.i <- list(data_jags = data_jags.i, 
                  plotcode.table = plotcode.table.i, 
                  species.table = species.table.i, 
                  country.table = country.table.i)
    # Add list i to the global list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- out.i")))
    
  }
  
  return(out)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Extract model outputs   ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Extract intensity for each iteration in each plot
#' @param jags.model.in a rjags object
#' @param data_jags.in input for jags.model.in
extract_intensity_per_plotcode <- function(jags.model.in, data_jags.in){
  
  ggs(as.mcmc(jags.model.in)) %>%
    # Remove deviance from parameters
    filter(Parameter != "deviance") %>%
    # Separate species code from parameter, and add disturbance
    mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter)) %>%
    # Only keep intensity
    filter(Parameter.true == "I") %>%
    mutate(plot = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter)))) %>%
    # Add species name
    left_join(data_jags.in$plotcode.table, by = "plot") %>%
    # Merge iteration and chain
    mutate(i = paste0("chain", Chain, "_iter", Iteration)) %>%
    # Finish formatting
    dplyr::select(plotcode, i, value) %>%
    spread(key = "i", value = "value")
  
}

#' Extract parameter per species for each iteration
#' @param jags.model.in a rjags object
#' @param data_jags.in input for jags.model.in
extract_param_per_species <- function(jags.model.in, data_jags.in){
  
  out <- ggs(as.mcmc(jags.model.in)) %>%
    # Remove deviance from parameters
    filter(Parameter != "deviance") %>%
    # Separate species code from parameter, and add disturbance
    mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter), 
           sp_country = gsub("\\]", "", gsub(".+\\[", "", Parameter)), 
           sp = as.integer(gsub("\\,.+", "", sp_country)), 
           co = ifelse(sp == sp_country, NA_integer_, as.integer(gsub(".+\\,", "", sp_country)))) %>%
    # Add country information
    left_join(data_jags.in$country.table, by = "co") %>%
    group_by(sp, Parameter.true, Chain, Iteration, country) %>%
    summarise(value = mean(value)) %>%
    ungroup() %>%
    # Merge iteration and chain
    mutate(iter = paste0("chain", Chain, "_iter", Iteration)) %>%
    # Remove intensity
    filter(Parameter.true != "I") %>%
    # Add species name
    left_join(data_jags.in$species.table, by = "sp") %>%
    # Finish formatting
    dplyr::select(species, iter, parameter = Parameter.true, value, country)
  
  # Fill the blanks for parameters that are not country specific
  out <- rbind.data.frame(
    subset(out, !is.na(country)), 
    merge(out %>% filter(is.na(country)) %>% dplyr::select(-country), 
          data.frame(country = data_jags.in$country.table$country))) %>%
    spread(key = "parameter", value = "value")
  
  return(out)
  
}




#' Predict tree-level death probability by averaging probability per iteration
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param dir.in Directory where to save the plots
predict_jags_meanProba <- function(jags.model, data_jags, data_model){
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Loop on all disturbances to extract data
  for(i in 1:length(disturbances.in)){
    
    # Format data 
    data.i <- data_model[[i]] %>%
      mutate(disturbance = disturbances.in[i],
             dbh.scaled = data_jags[[i]]$data_jags$dbh, 
             logratio.scaled = data_jags[[i]]$data_jags$logratio, 
             dead = ifelse(a == 1, 0, 1)) %>%
      dplyr::select(treecode, plotcode, country, species, dbh, logratio = log_dbh_dqm, 
                    disturbance, dbh.scaled, logratio.scaled, time, dead) %>%
      # Add intensity per plotcode
      left_join(extract_intensity_per_plotcode(jags.model[[i]], data_jags[[i]]), 
                by = "plotcode") %>%
      gather(key = "iter", value = "I", colnames(.)[grep("chain", colnames(.))]) %>%
      # Add parameter per species
      left_join(extract_param_per_species(jags.model[[i]], data_jags[[i]]), 
                by = c("species", "iter", "country")) %>%
      # Add missing parameter
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, NA_real_)) %>%
      # Calculate death probability depending on the model
      mutate(pd = case_when(disturbance %in% c("storm", "snow") ~ 1 - (1 - plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))^time, 
                            TRUE ~ 1 - (1 - plogis(a0 + b*I*dbh.scaled^c))^time)) %>%
      # Average probabilities by treecode
      group_by(treecode) %>%
      summarise(p = mean(pd), 
                p_05 = as.numeric(quantile(pd, probs = 0.05)), 
                p_95 = as.numeric(quantile(pd, probs = 0.95)))
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
  }
  
  # return the name of all the plots made
  return(out)
}



#' Predict tree-level death probability by averaging parameters per iteration
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param dir.in Directory where to save the plots
predict_jags_meanParam <- function(jags.model, data_jags, data_model){
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Loop on all disturbances to extract data
  for(i in 1:length(disturbances.in)){
    
    # Format data 
    data.i <- data_model[[i]] %>%
      mutate(disturbance = disturbances.in[i],
             dbh.scaled = data_jags[[i]]$data_jags$dbh, 
             logratio.scaled = data_jags[[i]]$data_jags$logratio, 
             dead = ifelse(a == 1, 0, 1)) %>%
      dplyr::select(treecode, plotcode, country, species, dbh, logratio = log_dbh_dqm, 
                    disturbance, dbh.scaled, logratio.scaled, time, dead) %>%
      # Add intensity per plotcode
      left_join(extract_intensity_per_plotcode(jags.model[[i]], data_jags[[i]]), 
                by = "plotcode") %>%
      gather(key = "iter", value = "I", colnames(.)[grep("chain", colnames(.))]) %>%
      # Add parameter per species
      left_join(extract_param_per_species(jags.model[[i]], data_jags[[i]]), 
                by = c("species", "iter", "country")) %>%
      # Add missing parameter
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, NA_real_)) %>%
      # Compute mean parameter and intensity
      group_by(treecode, plotcode, disturbance, logratio.scaled, dbh.scaled, time) %>%
      summarize(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c), I = mean(I)) %>% 
      ungroup() %>%
      # Calculate death probability depending on the model
      mutate(p = case_when(disturbance %in% c("storm", "snow") ~ 1 - (1 - plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))^time, 
                           TRUE ~ 1 - (1 - plogis(a0 + b*I*dbh.scaled^c))^time)) %>%
      dplyr::select(treecode, p)
    
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
  }
  
  # return the name of all the plots made
  return(out)
}


#' Get the list of all species involved in the model
#' @param data_model
get_species_list <- function(data_model){
  
  # Initialize output
  out <- c()
  
  # Extract species for each disturbance dataset
  for(i in 1:length(names(data_model))) out <- c(out, data_model[[i]]$species)
  
  # Remove duplicates 
  out <- unique(out)
  
  # return output
  return(out)
}

#' Get disturbance sensitivity by species
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM.
#' @param dbh.ref dbh value for which to predict values
#' @param logratio.ref value of the ratio dbh/dqm for which to predict values
#' @param I.ref disturbance intensity for which to predict values
get_disturbance_sensivity <- function(jags.model, data_jags, data_model, 
                                      dbh.ref = 300, logratio.ref = 0, I.ref = 0.5){
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Loop on all disturbances to extract data
  for(i in 1:length(disturbances.in)){
    
    # Model to scale logratio
    scale_logratio <- lm(logratio.scaled ~ logratio, 
                         data = data.frame(logratio = data_model[[i]]$log_dbh_dqm, 
                                           logratio.scaled = data_jags[[i]]$data_jags$logratio))
    # Model to scale dbh
    scale_dbh <- lm(dbh.scaled ~ dbh, 
                    data = data.frame(dbh = data_model[[i]]$dbh, 
                                      dbh.scaled = data_jags[[i]]$data_jags$dbh))
    
    # Parameter per species
    param_per_species.i <- extract_param_per_species(jags.model[[i]], data_jags[[i]])
    if(!(disturbances.in[[i]] %in% c("storm", "snow"))) param_per_species.i$a1 = NA_real_
    param_per_species.i <- param_per_species.i %>%
      left_join((data_model[[i]] %>%
                   group_by(species, country) %>%
                   summarize(weight = n())), 
                by = c("country", "species")) %>%
      mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
      group_by(species, iter) %>%
      summarise(a0 = sum(a0*weight)/sum(weight), 
                a1 = sum(a1*weight)/sum(weight), 
                b = sum(b*weight)/sum(weight), 
                c = sum(c*weight)/sum(weight))
    
    # Format data 
    data.i <- expand.grid(species = data_jags[[i]]$species.table$species,
                          dbh = dbh.ref, logratio = logratio.ref, I = I.ref, 
                          iter = unique(param_per_species.i$iter), 
                          disturbance = disturbances.in[i]) %>%
      # Scale variables when needed
      mutate(dbh.scaled = predict(scale_dbh, newdata = .), 
             logratio.scaled = predict(scale_logratio, newdata = .)) %>%
      # Add parameters per species
      left_join(param_per_species.i, by = c("species", "iter")) %>%
      # Compute probabilities
      mutate(pd = case_when(disturbance %in% c("storm", "snow") ~ plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c), 
                            TRUE ~ plogis(a0 + b*I*dbh.scaled^c)), 
             pd = 1 - (1 - pd)^5) %>%
      # Average per species
      group_by(species) %>%
      summarise(p = mean(pd), 
                p_025 = as.numeric(quantile(pd, probs = 0.025)), 
                p_975 = as.numeric(quantile(pd, probs = 0.975)))
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
  }
  
  # return the name of all the plots made
  return(out)
}


#' Get disturbance sensitivity by species for each iteration of the model
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM.
#' @param dbh.ref dbh value for which to predict values
#' @param logratio.ref value of the ratio dbh/dqm for which to predict values
#' @param I.ref disturbance intensity for which to predict values
get_disturbance_sensivity_full <- function(jags.model, data_jags, data_model, 
                                           dbh.ref = 300, logratio.ref = 0, I.ref = 0.5){
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Loop on all disturbances to extract data
  for(i in 1:length(disturbances.in)){
    
    # Model to scale logratio
    scale_logratio <- lm(logratio.scaled ~ logratio, 
                         data = data.frame(logratio = data_model[[i]]$log_dbh_dqm, 
                                           logratio.scaled = data_jags[[i]]$data_jags$logratio))
    # Model to scale dbh
    scale_dbh <- lm(dbh.scaled ~ dbh, 
                    data = data.frame(dbh = data_model[[i]]$dbh, 
                                      dbh.scaled = data_jags[[i]]$data_jags$dbh))
    
    # Parameter per species
    param_per_species.i <- extract_param_per_species(jags.model[[i]], data_jags[[i]])
    if(!(disturbances.in[[i]] %in% c("storm", "snow"))) param_per_species.i$a1 = NA_real_
    param_per_species.i <- param_per_species.i %>%
      left_join((data_model[[i]] %>%
                   group_by(species, country) %>%
                   summarize(weight = n())), 
                by = c("country", "species")) %>%
      mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
      group_by(species, iter) %>%
      summarise(a0 = sum(a0*weight)/sum(weight), 
                a1 = sum(a1*weight)/sum(weight), 
                b = sum(b*weight)/sum(weight), 
                c = sum(c*weight)/sum(weight))
    
    # Format data 
    data.i <- expand.grid(species = data_jags[[i]]$species.table$species,
                          dbh = dbh.ref, logratio = logratio.ref, I = I.ref, 
                          iter = unique(param_per_species.i$iter), 
                          disturbance = disturbances.in[i]) %>%
      # Scale variables when needed
      mutate(dbh.scaled = predict(scale_dbh, newdata = .), 
             logratio.scaled = predict(scale_logratio, newdata = .)) %>%
      # Add parameters per species
      left_join(param_per_species.i, by = c("species", "iter")) %>%
      # Compute probabilities
      mutate(pd = case_when(disturbance %in% c("storm", "snow") ~ plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c), 
                            TRUE ~ plogis(a0 + b*I*dbh.scaled^c)), 
             pd = 1 - (1 - pd)^5) %>%
      # Keep variables of interest
      dplyr::select(species, iter, p = pd)
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
  }
  
  # return the name of all the plots made
  return(out)
}




#' Export jags model as R object
#' @param jags.model.in List of fags model to save (one list item per disturbance type)
#' @param data_jags.in list of dajagas object to save (one list item per disturbance type)
#' @param data_jags.in list of data_model object to save (one list item per disturbance type)
#' @param file.in name of the file to save, including path
export_jags <- function(jags.model.in, data_jags.in, data_model.in, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Initialize correspondence table
  corresp.tables <- list()
  # Initialize weight table
  weight.tables <- list()
  # Initialize scaling table
  scale.tables <- list()
  
  # Loop on all disturbances 
  for(i in 1:length(names(data_jags.in))){
    
    # Table containing the parameters to scale dbh and logratio
    # -- First run some models
    mod.dbh.i <- lm(dbh.scaled ~ dbh, data = data.frame(
      dbh = data_model.in[[i]]$dbh, dbh.scaled = data_jags.in[[i]]$data_jags$dbh))
    mod.logratio.i <- lm(logratio.scaled ~ logratio,  data = data.frame(
      logratio = data_model.in[[i]]$log_dbh_dqm, logratio.scaled = data_jags.in[[i]]$data_jags$logratio))
    # -- Store coefficients in a table
    scale.table.i <- list(dbh.intercept = summary(mod.dbh.i)$coefficients[1, 1], 
                          dbh.slope = summary(mod.dbh.i)$coefficients[2, 1], 
                          logratio.intercept = summary(mod.logratio.i)$coefficients[1, 1], 
                          logratio.slope = summary(mod.logratio.i)$coefficients[2, 1])
    # -- Add to the list of scaling tables
    eval(parse(text = paste0("scale.tables$", names(data_jags.in)[i], " <- scale.table.i")))
    
    
    
    # Table containing the weight of each species per disturbance type
    weight.table.i <- data_model.in[[i]] %>%
      group_by(species, country) %>%
      summarize(weight = n())
    eval(parse(text = paste0("weight.tables$", names(data_jags.in)[i], " <- weight.table.i")))
    
    
    
    
    # Correspondence tables
    corresp.table.i <- data_jags.in[[i]][-1]
    eval(parse(text = paste0("corresp.tables$", names(data_jags.in)[i], " <- corresp.table.i")))
  }
  
  # Rename jags object
  jags.list <- jags.model.in
  
  # Save the two objects
  save(corresp.tables, jags.list, scale.tables, weight.tables, file = file.in)
  
  # Return the name of the file
  return(file.in)
}



#' Export sensitivity
#' @param disturbance_sensitivity_full sensitivity for fire, other and storm for all mcmc iterations
#' @param disturbance_sensitivity_full_bis sensitivity for snow and biotic for all mcmc iterations
#' @param file.in name of the file to save, including path
export_sensitivity <- function(disturbance_sensitivity_full, disturbance_sensitivity_full_bis, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Save the two objects
  save(disturbance_sensitivity_full, disturbance_sensitivity_full_bis, file = file.in)
  
  # Return the name of the file
  return(file.in)
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Extract general information ---------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Export a latex table with statistics about disturbance and number of trees and plots
#' @param FUNDIV_tree Tree table formatted for FUNDIV
#' @param FUNDIV_plot Plot table formatted for FUNDIV
#' @param file.in Name of the file to save
export_table_disturbance_stats <- function(FUNDIV_tree, FUNDIV_plot, file.in){
  
  # Initialize final table
  table.out <- data.frame(col1 = "")
  
  # Table with the information on all countries 
  table.allcountries <- FUNDIV_tree %>%
    left_join((FUNDIV_plot %>% dplyr::select(plotcode, disturbance.nature)), 
              by = "plotcode") %>%
    left_join((FUNDIV_plot %>% group_by(country) %>% summarize(n.plots.per.country = n())), 
              by = "country") %>%
    filter(disturbance.nature %in% c("storm", "fire", "other", "snow", "biotic")) %>%
    filter(treestatus != 1) %>%
    mutate(dead = ifelse(treestatus == 2, 0, 1)) %>%
    group_by(country, plotcode, disturbance.nature, n.plots.per.country) %>%
    summarize(n.trees = n(), 
              death.rate = sum(dead, na.rm = TRUE)/n(),
              ba.plot = sum(ba_ha1, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(country, disturbance.nature) %>%
    summarize(n.plots = n(), 
              perc.plots.per.country = round(100*(n()/n.plots.per.country), digits = 2),
              n.trees.per.plot = mean(n.trees, na.rm = TRUE), 
              n.trees = sum(n.trees, na.rm = TRUE), 
              death.rate = mean(death.rate, na.rm = TRUE), 
              basal.area = mean(ba.plot)) %>%
    distinct() %>%
    gather(key = "variable", value = "value", colnames(.)[c(3:dim(.)[2])]) %>%
    mutate(value = ifelse(variable %in% c("death.rate", "perc.plots.per.country"), 
                          as.character(round(value, digits = 2)), 
                          as.character(round(value, digits = 0)))) %>%
    mutate(variable = gsub("\\.", "\\ ", variable)) %>%
    spread(key = "disturbance.nature", value = "value") %>%
    mutate(across(everything(), ~replace_na(.x, "")))
  
  # Finish to format the first line of the reference table
  for(i in 3:dim(table.allcountries)[2]){
    eval(parse(text = paste0("table.out$col", (i-1), " <- colnames(table.allcountries)[", i, "]")))
  }
  
  # Create one table per country, and add it to the final table
  for(j in 1:length(unique(table.allcountries$country))){
    # First line of the table for country j
    head.table.j <- as.data.frame(
      matrix(nrow = 1, ncol = dim(table.out)[2], data = "", 
             dimnames = list(c(""), colnames(table.out)))) %>%
      mutate(col1 = toupper(unique(table.allcountries$country)[j]))
    # Add an empty line on top to separate countries
    if(j > 1) head.table.j <- rbind.data.frame(
      matrix(nrow = 1, ncol = dim(table.out)[2], data = "", 
             dimnames = list(c(""), colnames(table.out))),
      head.table.j)
    # Table with values for country j
    table.country.j <- table.allcountries %>%
      filter(country == unique(table.allcountries$country)[j]) %>%
      ungroup() %>%
      dplyr::select(-country)
    colnames(table.country.j) <- colnames(head.table.j)
    # Add to the final table
    table.out <- rbind.data.frame(table.out, head.table.j, table.country.j)
  }
  
  # Remove column names
  colnames(table.out) <- NULL
  
  # Create output dir if necessary
  create_dir_if_needed(file.in)
  
  # create a tex file
  print(xtable(table.out, type = "latex", label = "dist_stat_table",
               caption = "Mean number of trees, plots, number of trees per plot, 
               basal area, death rate and percentage of disturbed plots in each 
               country and for each disturbance type"), 
        include.rownames=FALSE, hline.after = c(0, 1, dim(table.out)[1]), 
        include.colnames = FALSE, caption.placement = "top", file = file.in)
  
  #return output
  return(file.in)
}
