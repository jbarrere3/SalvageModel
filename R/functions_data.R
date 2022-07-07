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
  # Create a dataset with one line per species
  species <- FUNDIV_tree %>%
    dplyr::select(species) %>%
    distinct() %>%
    # Add the genus
    mutate(genus = gsub(" .+", "", species), 
           sp = NA_character_)
  # Add species
  for(i in 1:dim(species)[1]) species$sp[i] <- strsplit(species$species[i], split = " ")[[1]][2]
  # Set to NA when species name is not known
  species <- species %>%
    mutate(sp = ifelse(sp %in% c("sp.", "sp", "", "x"), NA_character_, sp),
           speciesname = ifelse(is.na(sp), NA_character_, paste(genus, sp, sep = " "))) 
  
  # Create a table connecting species to family
  species_to_family <- species %>%
    filter(!is.na(speciesname)) %>%
    dplyr::select(speciesname) %>%
    distinct()
  species_to_family <- tax_name(sci = species_to_family$speciesname, get = "family", db = "ncbi", messages = FALSE)
  
  # Create a table linking the genus to the family
  genus_to_family <- species_to_family %>%
    dplyr::select(speciesname = query, family) %>%
    mutate(genus = gsub(" .+", "", speciesname)) %>%
    dplyr::select(genus, family) %>%
    filter(!is.na(family)) %>%
    distinct()
  
  # Create a table linking the family to the branch
  family_to_branch <- tpl_families()
  
  # Merge all datasets
  species.final <- species %>%
    mutate(genus = gsub(" ", "", genus)) %>%
    left_join(genus_to_family, by = "genus") %>%
    left_join(family_to_branch, by = "family")
  species.final <- species.final %>%
    left_join((species.final %>% 
                 filter(!is.na(speciesname)) %>%
                 mutate(group2 = ifelse(family %in% c("Fagaceae", "Fabaceae"), "Angiosperms", group)) %>%
                 dplyr::select(genus, group2) %>%
                 filter(!is.na(group2)) %>%
                 distinct()), 
              by = "genus") %>%
    dplyr::select(species, group = group2, family, genus) %>%
    distinct()
  
  return(species.final)
  
}

#' Get family and order for each species present in the dataset
#' @param FUNDIV_tree Tree table with a column entitled "species"
get_species_info2 <- function(FUNDIV_tree){
  
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
    mutate(enough.individuals = ifelse((n.per.species > 200 & n.plot.per.species > 20), 1, 0), 
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
             dqm = sqrt(sum(dbh2.W, na.rm = TRUE)/sum(weight1, na.rm  = TRUE))) %>%
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
                    sgdd, wai, dqm, log_dbh_dqm)
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
    
  }
  
  # return output list
  return(out)
}



#' Center and Scale variables before fitting one fo the models that includes dqm
#' @param data_model Tree data formated for the IPM. 
#' @param var Variables to scale (character vector)
scale_data_model <- function(data_model,  var = c("dbh", "dqm", "comp", "sgdd", "wai", "log_dbh_dqm")){
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
      mutate(a1 = ifelse(disturbance == "storm", a1, NA_real_)) %>%
      # Calculate death probability depending on the model
      mutate(pd = case_when(disturbance == "storm" ~ 1 - (1 - plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))^time, 
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
      mutate(a1 = ifelse(disturbance == "storm", a1, NA_real_)) %>%
      # Compute mean parameter and intensity
      group_by(treecode, plotcode, disturbance, logratio.scaled, dbh.scaled, time) %>%
      summarize(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c), I = mean(I)) %>% 
      ungroup() %>%
      # Calculate death probability depending on the model
      mutate(p = case_when(disturbance == "storm" ~ 1 - (1 - plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))^time, 
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
    if(disturbances.in[[i]] != "storm") param_per_species.i$a1 = NA_real_
    param_per_species.i <- param_per_species.i %>%
      group_by(species, iter) %>%
      summarise(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c))
    
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
      mutate(pd = case_when(disturbance == "storm" ~ plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c), 
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



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Manage traits data      ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#' Compile all traits data
#' @param bark.thickness_file file containing data on bark thickness
#' @param wood.density_file file containing data on wood density
#' @param FUNDIV_tree original tree dataset (to get height diameter ratio)
#' @param species.in character vector of all the species for which to extract the data
compile_traits <- function(wood.density_file, shade.tolerance_file, root.depth_file, 
                           FUNDIV_tree, bark.thickness_file, species.in){
  
  # Wood density (Chave et al. 2008 + Dryad to quote)
  wood.density <- read_xls(wood.density_file, sheet = "Data")
  colnames(wood.density) <- c("n", "family", "species", "wood.density_g.cm3", "region", "reference")
  
  # Shade tolerance (Niimenets)
  shade.tolerance <- fread(shade.tolerance_file) %>%
    mutate(shade.tolerance = as.numeric(gsub("\\,", "\\.", shade.tolerance.mean))) %>%
    dplyr::select(species, shade.tolerance)
  
  # Rooting depth (Guerrero-Ramirez et al. 2021 - Groot database)
  root.depth <- fread(root.depth_file) %>%
    mutate(species = paste(genusTNRS, speciesTNRS, sep = " ")) %>%
    filter(traitName %in% c("Rooting_depth", "Root_mass_fraction")) %>%
    filter(species %in% species.in) %>%
    dplyr::select("species", "traitName", "meanSpecies") %>%
    spread(key = "traitName", value = "meanSpecies")
  
  # Bark thickness from NFI data
  bark.thickness <- fread(bark.thickness_file)
  
  ## - Height diameter ratio from FUNDIV_tree
  # Format data (remove unknown height and species)
  FUNDIV_tree.in <- FUNDIV_tree %>%
    filter(height1 > 0) %>%
    filter(!is.na(species))
  # Initialize output
  height.dbh.ratio.data <- data.frame(species = character(0), 
                                      height.dbh.ratio = numeric(0))
  # Loop on all species
  for(sp in unique(FUNDIV_tree.in$species)){
    model.sp = lm(height1 ~ dbh1, data = subset(FUNDIV_tree.in, species == sp))
    # The ratio is the slope of the regression height vs dbh
    height.dbh.ratio.data <- rbind.data.frame(
      height.dbh.ratio.data, data.frame(species = sp, 
                                        height.dbh.ratio = as.numeric(coefficients(model.sp)[2])))
  }
  
  # Global trait dataset
  traits <- data.frame(species = species.in) %>%
    # Add wood density
    left_join((wood.density %>% 
                 group_by(species) %>%
                 summarize(wood.density_g.cm3 = mean(wood.density_g.cm3))),
              by = "species") %>%
    # Add shade tolerance
    left_join(shade.tolerance, by = "species") %>%
    # Add shade tolerance
    left_join(bark.thickness, by = "species") %>%
    # Add root depth
    left_join(root.depth, by = "species") %>%
    # Add height dbh ratio
    left_join(height.dbh.ratio.data, by = "species")
  
  return(traits)
}

#' Compile traits related to the climatic conditions eahc species is exposed to
#' @param data_model Data frame used to fit the model, with variables unscaled
#' @param gbif_file File containing gbif climatic data
compile_traits_climate <- function(data_model, gbif_file){
  
  gbif_data <- fread(gbif_file) %>%
    dplyr::select(species, mean_mat, mean_tmin, mean_map)
  
  data_model %>%
    group_by(species) %>%
    summarize(sgdd.mean = mean(sgdd, na.rm = TRUE), 
              wai.mean = mean(wai, na.rm = TRUE), 
              comp.mean = mean(comp, na.rm = TRUE)) %>%
    left_join(gbif_data, by = "species")
}


#' Compile all traits data
#' @param TRY_file file containing TRY request
#' @param species.in character vector of all the species for which to extract the data
compile_traits_TRY <- function(TRY_file, species.in){
  
  # -- Translate TRY traits code into abbreviated traits name
  TRY.traits.name <- data.frame(
    TraitID = c(24, 3117, 146, 14, 56, 15, 46, 65, 1111, 2809, 2807, 2808, 159, 30, 318, 
                31, 719, 59, 819, 45, 773, 413, 324, 1229, 153, 865, 837, 3446), 
    trait = c("bark.thickness", "leaf.sla", "leaf.CN.ratio", "leaf.N.mass", "leaf.NP.ratio", 
              "leaf.P.mass", "leaf.thickness", "root.type", "seedbank.density", 
              "seedbank.duration", "seedbank.n.layers", "seedbank.thickness.toplayer", 
              "seedbank.type", "tolerance.drought", "tolerance.fire", "tolerance.frost", 
              "xylem.hydraulic.vulnerability", "plant.lifespan", "plant.resprouting.capacity", 
              "stomata.conductance", "crown.height", "leaf.Chl.content", "crown.length", 
              "wood.Nmass", "budbank.height.distribution", "budbank.seasonality", 
              "bark.structure", "plant.biomass")
  )
  
  
  
  ## -- Compile numeric traits data
  traits.TRY <-data.table::fread(TRY_file) %>%
    filter(!is.na(TraitID)) %>%
    filter(!is.na(StdValue)) %>%
    filter(AccSpeciesName %in% species.in) %>%
    left_join(TRY.traits.name, by = "TraitID") %>%
    mutate(trait = paste("TRY", trait, gsub("\\ ", "", UnitName), sep = "_"), 
           trait = gsub("\\/", ".", trait)) %>%
    rename("species" = "AccSpeciesName") %>%
    group_by(species, trait) %>%
    summarize(trait.value = mean(StdValue, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(trait) %>%
    mutate(n.species.per.trait = n()) %>%
    filter(n.species.per.trait >= 10) %>% 
    dplyr::select(-n.species.per.trait) %>%
    spread(key = trait, value = trait.value)  
  
  return(traits.TRY)
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Outdated functions ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Generate species-specific parameters to simulate data
#' @param data.jags.in dataframe containing plot and tree variable centered and scaled
#' @param param character vector containing the name of all the parameters to initialize
generate_parameters_sp <- function(data.jags.in, param = paste0("c", c(0:8))){
  out <- data.frame(sp = data.jags.in$species_table$sp)
  for(i in 1:length(param)){
    # intercept
    if(param[i] %in% c("c0", "c3", "c6")){
      out$value <- round(runif(dim(out)[1], -6, -4), digits = 3)
    } 
    # multiplicative parameter
    if(param[i] %in% c("c1", "c4", "c7")){
      out$value <- round(runif(dim(out)[1], 3, 6), digits = 3)
    } 
    # power parameter
    if(param[i] %in% c("c2", "c5", "c8")){
      out$value <- round(rnorm(dim(out)[1], mean = 0, sd = 0.5), digits = 1)
    } 
    
    colnames(out)[i+1] <- param[i]
  }
  return(out)
}





#' Simulate tree status data for the jags model with disturbance given and specific background
#' mortality equations (see supplementary of Kunstler et al. 2020 for original equations)
#' @param data.jags.in dataframe containing plot and tree variable centered and scaled
#' @param parameters_sp dataframe containing the parameters value per species
#' @author Björn Reineking, Julien Barrere
simulate_status_full_sub <- function(data.jags.in, parameters_sp){
  
  
  
  # species table extended with parameters value
  species.table.extended <- data.jags.in$species_table %>%
    left_join(parameters_sp, by = "sp")
  
  # plot table extended with disturbance intensity
  plot.table.extended <- data.frame(plot = data.jags.in$data_jags$plot, 
                                    d = data.jags.in$data_jags$d, 
                                    Dfire = data.jags.in$data_jags$Dfire, 
                                    Dstorm = data.jags.in$data_jags$Dstorm, 
                                    Dother = data.jags.in$data_jags$Dother) %>%
    group_by(plot, Dfire, Dstorm, Dother) %>%
    summarise(severity = sum(d)/n()) %>%
    ungroup() %>%
    mutate(Intensity.fire = rbeta(dim(.)[1], 0.66, 0.36),
           Intensity.fire = case_when(Intensity.fire < 0.001 ~ 0.001, 
                                      Intensity.fire > 0.999 ~ 0.999, 
                                      TRUE ~ Intensity.fire), 
           Intensity.storm = rbeta(dim(.)[1], 0.65, 2.66),
           Intensity.storm = case_when(Intensity.storm < 0.001 ~ 0.001, 
                                       Intensity.storm > 0.999 ~ 0.999, 
                                       TRUE ~ Intensity.storm), 
           Intensity.other = rbeta(dim(.)[1], 0.48, 1.77),
           Intensity.other = case_when(Intensity.other < 0.001 ~ 0.001, 
                                       Intensity.other > 0.999 ~ 0.999, 
                                       TRUE ~ Intensity.other), 
           Ifire = Dfire*Intensity.fire, 
           Istorm = Dstorm*Intensity.storm, 
           Iother = Dother*Intensity.other) %>%
    ungroup() %>%
    dplyr::select(plot, Ifire, Istorm, Iother)
  
  # Calculate the state for each tree
  newstatus.table <- data.frame(plot = data.jags.in$data_jags$plot, 
                                sp = data.jags.in$data_jags$sp, 
                                time = data.jags.in$data_jags$time, 
                                dbh = data.jags.in$data_jags$dbh, 
                                Dfire = data.jags.in$data_jags$Dfire, 
                                Dstorm = data.jags.in$data_jags$Dstorm, 
                                Dother = data.jags.in$data_jags$Dother) %>%
    # Add species parameter
    left_join(species.table.extended, by = "sp") %>%
    # Add intensity
    left_join(plot.table.extended, by = "plot") %>%
    # calculate probabilities for the three different states
    mutate(pdstorm = 1 - (1 - plogis(c0 + c1*Istorm*dbh^c2))^time, 
           pdother = 1 - (1 - plogis(c3 + c4*Iother*dbh^c5))^time, 
           pdfire = 1 - (1 - plogis(c6 + c7*Ifire*dbh^c8))^time, 
           pdD = 1 - (1 - Dstorm*pdstorm)*(1 - Dother*pdother)*(1 - Dfire*pdfire), 
           d = NA_integer_)
  
  # death probability
  for(i in 1:dim(newstatus.table)[1]) newstatus.table$d[i] <- rbinom(1, 1, newstatus.table$pdD[i])
  
  # Final output
  out <- data.jags.in
  # Add disturbance intensity
  out$data_jags$Istorm = newstatus.table$Istorm
  out$data_jags$Iother = newstatus.table$Iother
  out$data_jags$Ifire = newstatus.table$Ifire
  # Add tree status
  out$data_jags$d <- newstatus.table$d
  
  return(out)
}








#' Get disturbance sensitivity based on the model with dominance
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param data_model_scaled Tree data formatted for the IPM and scaled 
get_disturbance_sensivity_dqm <- function(jags.model, data_jags, data_model_scaled, data_model, 
                                          disturbance_species_info){
  
  # Identify parameters per species
  param_per_species <- ggs(as.mcmc(jags.model)) %>%
    filter(Parameter != "deviance") %>%
    mutate(sp = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
           Parameter = gsub("\\[.+\\]", "", Parameter)) %>%
    left_join(data_jags$species_table, by = "sp") %>% 
    group_by(species, Parameter) %>%
    dplyr::summarize(mean = mean(value)) %>%
    spread(key = "Parameter", value = "mean")
  
  # Model to scale dbh
  scale_dbh <- lm(dbh.scaled ~ dbh, 
                  data = data.frame(dbh = data_model$dbh, 
                                    dbh.scaled = data_model_scaled$dbh))
  
  # Disturbance - species combination to keep 
  disturbance.species.to.keep <- (disturbance_species_info$species_disturbance_table %>%
                                    mutate(ID = paste(species, disturbance.type, sep = ".")) %>%
                                    filter(n.indiv >= 50 & n.plot >= 5))$ID
  
  # dbh range per species
  dbh.range.per.sp <- data_model %>%
    group_by(species) %>%
    dplyr::summarize(dbh.min = min(dbh, na.rm = TRUE), 
                     dbh.max = max(dbh, na.rm = TRUE))
  
  # - Pre-format the data before plotting
  out <- expand.grid(species = unique(param_per_species$species), 
                     dbh = c(100, 250, 700), 
                     Intensity = 0.7) %>%
    mutate(dbh.scaled = predict(scale_dbh, newdata = .)) %>%
    # Remove dbh outside of dbh range per species
    left_join(dbh.range.per.sp, by = "species") %>%
    filter(dbh >= dbh.min & dbh <= dbh.max) %>%
    dplyr::select(-dbh.min, -dbh.max) %>%
    # Add parameters per species
    left_join(param_per_species, by = "species") %>%
    # Calculate death probabilities by species and by disturbance
    mutate(storm = plogis(st3 + st4*Intensity*dbh.scaled^st5), 
           other = plogis(ot3 + ot4*Intensity*dbh.scaled^ot5), 
           fire = plogis(fi3 + fi4*Intensity*dbh.scaled^fi5)) %>%
    dplyr::select(species, dbh, storm, other, fire) %>%
    gather(key = "disturbance", value = "sensitivity", "storm", "other", "fire") %>%
    # filter predictions for which we don't have enough data
    mutate(ID = paste(species, disturbance, sep = ".")) %>%
    filter(ID %in% disturbance.species.to.keep) %>%
    mutate(trait.name = paste0(disturbance, ".dbh", dbh)) %>%
    dplyr::select(species, trait.name, sensitivity) %>%
    spread(key = trait.name, value = sensitivity)
  
  
  # return the formatted dataset
  return(out)
}




#' Create a latex table with the model results
#' @param traits dataframe contianing trait values per species
#' @param disturbance_sensivity dataframe containing the sensitivity to each disturbance
#' @param file.in Name of the file to save
export_trait_result_latex <- function(traits, disturbance_sensitivity, file.in){
  # Initiate tables for fire, storm and other
  tables <- list()
  
  # Loop on all traits
  for(i in 1:(dim(traits)[2] - 1)){
    # Identify the name of trait i
    trait.i <- colnames(traits)[i+1]
    # Create a table with only species and trait i
    traits.i <- traits %>% dplyr::select("species", "trait" = trait.i)
    # Loop on all type of disturbances
    for(j in 1:(dim(disturbance_sensitivity)[2] - 1)){
      # Identify disturbance j
      disturbance.j <- colnames(disturbance_sensitivity)[j+1]
      # Create a table with trait i and sensitivity to disturbance j
      data.ij <- traits.i %>%
        left_join((disturbance_sensitivity %>%
                     dplyr::select("species", "sensitivity" = disturbance.j) %>%
                     mutate(sensitivity.logit = log(sensitivity/(1 - sensitivity)))), 
                  by = "species") %>%
        drop_na()
      # Fit a model
      model.ij <- lm(sensitivity.logit ~ trait, data = data.ij)
      # Results
      table.ij <- data.frame(
        trait = trait.i, 
        Est. = as.character(round(summary(model.ij)$coefficients[2, 1], digits = 3)), 
        Fval. = as.character(round(anova(model.ij)[1, 4], digits = 2)), 
        p = scales::pvalue(anova(model.ij)[1, 5], accuracy = 0.01)
      )
      # Add to the list containing the final results
      if(i == 1) eval(parse(text = paste0("tables$", disturbance.j, " <- table.ij ")))
      else eval(parse(text = paste0("tables$", disturbance.j, " <- rbind.data.frame(tables$", disturbance.j, ", table.ij)")))
    }
  }
  
  # Build the final dataframe
  out <- rbind.data.frame(
    data.frame(col1 = c("", "traits"), 
               col2 = c("", "Est."), col3 = c("fire", "F value"), col4 = c("", "p"), 
               col5 = c("", "Est."), col6 = c("storm", "F value"), col7 = c("", "p"), 
               col8 = c("", "Est."), col9 = c("other", "F value"), col10 = c("", "p")), 
    tables$fire %>%
      rename(col1 = trait, col2 = Est., col3 = Fval., col4 = p) %>%
      left_join((tables$storm %>%
                   rename(col1 = trait, col5 = Est., col6 = Fval., col7 = p)), 
                by = "col1") %>%
      left_join((tables$other %>%
                   rename(col1 = trait, col8 = Est., col9 = Fval., col10 = p)), 
                by = "col1")
  )
  
  colnames(out) <- NULL
  
  # Create output dir if necessary
  create_dir_if_needed(file.in)
  
  # create a tex file
  print(xtable(out, type = "latex"), file = file.in)
  
  #return output
  return(file.in)
}
