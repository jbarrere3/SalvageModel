#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_data.R  
#' @description R script containing all functions relative to data
#               importation and formatting
#' @author Lukas Heiland (re-arranged by Julien BARRERE)
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Format the data before fiting IPM or mortality model
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plot FUNDIV plot table
#' @param Climate dataframe containing climatic data per plot
#' @param BM_equations dataframe indicating which background mortality to use per species
format_data_model <- function(FUNDIV_tree, FUNDIV_plot, Climate, BM_equations){
  
  ## - Step 1 - Format species information
  species_bm <- FUNDIV_tree %>%
    # Remove ingrowth trees
    filter(treestatus != 1) %>% 
    # Count number of individual per species
    group_by(species) %>% 
    summarize(n = n()) %>%
    # Rem0ve species with less than 4000 individual trees
    filter(n > 4000) %>%
    arrange(desc(n)) %>%
    # Add table with info on background mortality equations
    left_join(BM_equations, by = "species") %>%
    # Turn it into a binary variable
    mutate(BM_eq = ifelse(is.na(BM_eq), 1, BM_eq), # If no info on BM equation, use the simpler
           bm1 = ifelse(BM_eq == 1, 1, 0), 
           bm2 = ifelse(BM_eq == 2, 1, 0), 
           bm3 = ifelse(BM_eq == 3, 1, 0), 
           bm4 = ifelse(BM_eq == 4, 1, 0), 
           bm5 = ifelse(BM_eq == 5, 1, 0), 
           bm6 = ifelse(BM_eq == 6, 1, 0)) %>%
    dplyr::select(-BM_eq)
  
  ## - Ajust climate depending on the input dataset
  if(class(FUNDIV_tree$plotcode) == "integer"){
    Climate <- Climate %>% 
      filter(plotcode %in% as.character(FUNDIV_tree$plotcode)) %>%
      mutate(plotcode = as.integer(plotcode))
  } 
  
  ## - Step 2 - Format final dataset
  out <- FUNDIV_tree %>%
    # Add plot level data 
    left_join(Climate, by = "plotcode") %>% 
    left_join((FUNDIV_plot %>% 
                 dplyr::select(plotcode, time = yearsbetweensurveys, 
                               disturbance.severity, disturbance.nature) %>%
                 mutate(DS = case_when(disturbance.severity == 0 ~ 0, 
                                       disturbance.severity == 1 ~ 0.125, 
                                       disturbance.severity == 2 ~ 0.375, 
                                       disturbance.severity == 3 ~ 0.625, 
                                       disturbance.severity == 4 ~ 0.875))), 
              by = "plotcode") %>% 
    # Keep plots only composed of species present in large number in the dataset
    filter(plotcode %in% (FUNDIV_tree %>% 
                            filter(treestatus != 1) %>%
                            mutate(has_bm = ifelse(species %in% species_bm$species, 1, 0)) %>%
                            group_by(plotcode) %>%
                            summarize(prop.bm = sum(has_bm)/n()) %>%
                            filter(prop.bm == 1))$plotcode) %>% 
    # Add species information
    left_join(species_bm, by = "species") %>%
    # Compute basal area of competitors
    group_by(plotcode) %>%
    mutate(comp=sum(ba_ha1) - ba_ha1) %>%
    ungroup() %>%
    filter(treestatus != 1) %>%
    mutate(h = case_when(treestatus == 3 ~ 1, T ~ 0), 
           d = case_when(treestatus %in% c(4, 5) ~ 1, T ~ 0), 
           a = case_when(treestatus == 2 ~ 1, T ~ 0), 
           disturbance.nature = ifelse(is.na(disturbance.nature), "other", disturbance.nature),
           Dfire = ifelse(disturbance.nature == "fire", 1, 0), 
           Dstorm = ifelse(disturbance.nature == "storm", 1, 0),
           Dother = ifelse(disturbance.nature == "other", 1, 0), 
           D = 1 - (1 - Dfire)*(1 - Dstorm)*(1 - Dother)) %>%
    filter(!is.na(sgdd) & !is.na(wai) & !is.na(comp)) %>%
    # Select variables and remove NAs
    dplyr::select(plotcode, treecode, species, time, h, d, a, dbh = dbh1, comp, sgdd, 
                  wai, D, Dfire, Dstorm, Dother, DS, bm1, bm2, bm3, bm4, bm5, bm6) 
}

#' Format the data before fiting IPM or mortality model
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plot FUNDIV plot table
#' @param Climate dataframe containing climatic data per plot
#' @param BM_equations dataframe indicating which background mortality to use per species
#' @param species dataset containing species name, genus, family and branch of each species
format_data_model_full <- function(FUNDIV_tree, FUNDIV_plot, Climate, BM_equations, species){
  
  ## - Step 1 - Format species information
  species_bm <- FUNDIV_tree %>%
    # Remove ingrowth trees
    filter(treestatus != 1) %>% 
    # Count number of individual per species
    group_by(species) %>% 
    summarize(n.species = n()) %>%
    # Add species information
    left_join((species %>% dplyr::select(species, group)), by = "species") %>%
    # Adjust species name depending on number
    mutate(species2 = case_when((n.species < 5000 & group == "Angiosperms") ~ "Other broadleaf", 
                                (n.species < 5000 & group == "Gymnosperms") ~ "Other conifer",
                                TRUE ~ species)) %>%
    group_by(species2) %>%
    mutate(n.species2 = sum(n.species)) %>%
    filter(n.species2 > 5000) %>%
    # Add table with info on background mortality equations
    left_join(BM_equations, by = c("species2" = "species")) %>%
    # Turn it into a binary variable
    mutate(BM_eq = ifelse(is.na(BM_eq), 1, BM_eq), # If no info on BM equation, use the simpler
           bm1 = ifelse(BM_eq == 1, 1, 0), 
           bm2 = ifelse(BM_eq == 2, 1, 0), 
           bm3 = ifelse(BM_eq == 3, 1, 0), 
           bm4 = ifelse(BM_eq == 4, 1, 0), 
           bm5 = ifelse(BM_eq == 5, 1, 0), 
           bm6 = ifelse(BM_eq == 6, 1, 0)) %>%
    dplyr::select(-BM_eq, -group, -n.species, -n.species2)
  
  ## - Ajust climate depending on the input dataset
  if(class(FUNDIV_tree$plotcode) == "integer"){
    Climate <- Climate %>% 
      filter(plotcode %in% as.character(FUNDIV_tree$plotcode)) %>%
      mutate(plotcode = as.integer(plotcode))
  } 
  
  ## - Step 2 - Format final dataset
  out <- FUNDIV_tree %>%
    # Add plot level data 
    left_join(Climate, by = "plotcode") %>% 
    left_join((FUNDIV_plot %>% 
                 dplyr::select(plotcode, time = yearsbetweensurveys, 
                               disturbance.severity, disturbance.nature) %>%
                 mutate(DS = case_when(disturbance.severity == 0 ~ 0, 
                                       disturbance.severity == 1 ~ 0.125, 
                                       disturbance.severity == 2 ~ 0.375, 
                                       disturbance.severity == 3 ~ 0.625, 
                                       disturbance.severity == 4 ~ 0.875))), 
              by = "plotcode") %>% 
    # Keep plots only composed of species present in large number in the dataset
    filter(plotcode %in% (FUNDIV_tree %>% 
                            filter(treestatus != 1) %>%
                            mutate(has_bm = ifelse(species %in% species_bm$species, 1, 0)) %>%
                            group_by(plotcode) %>%
                            summarize(prop.bm = sum(has_bm)/n()) %>%
                            filter(prop.bm == 1))$plotcode) %>% 
    # Add species information
    left_join(species_bm, by = "species") %>%
    # Compute basal area of competitors
    group_by(plotcode) %>%
    mutate(comp=sum(ba_ha1) - ba_ha1) %>%
    ungroup() %>%
    filter(treestatus != 1) %>%
    mutate(h = case_when(treestatus == 3 ~ 1, T ~ 0), 
           d = case_when(treestatus %in% c(4, 5) ~ 1, T ~ 0), 
           a = case_when(treestatus == 2 ~ 1, T ~ 0), 
           disturbance.nature = ifelse(is.na(disturbance.nature), "other", disturbance.nature),
           Dfire = ifelse(disturbance.nature == "fire", 1, 0), 
           Dstorm = ifelse(disturbance.nature == "storm", 1, 0),
           Dother = ifelse(disturbance.nature == "other", 1, 0), 
           D = 1 - (1 - Dfire)*(1 - Dstorm)*(1 - Dother)) %>%
    filter(!is.na(sgdd) & !is.na(wai) & !is.na(comp)) %>%
    # Select variables and remove NAs
    dplyr::select(plotcode, treecode, species = species2, time, h, d, a, dbh = dbh1, comp, sgdd, 
                  wai, D, Dfire, Dstorm, Dother, DS, bm1, bm2, bm3, bm4, bm5, bm6) 
}


#' Center and Scale variables before fitting a model
#' @param data_model Tree data formated for the IPM. 
#' @param var Variables to scale (character vector)
scale_data_model <- function(data_model, 
                              var = c("dbh", "comp", "sgdd", "wai", "DA.Senf", "DS.Senf", "DS")){
  id_var <- which(colnames(data_model) %in% var) 
  if("DS" %in% var) out <- data_model %>% mutate(DS = ifelse(DS == 0, NA_real_, DS))
  else out <- data_model
  unscaled = as.matrix(out[, id_var])
  scaled = scale(unscaled)
  out[, id_var] <- as.data.frame(scaled)
  if("DS" %in% var) out$DS[which(is.na(out$DS))] <- -99
  # dbh is logged in BM equation, so need to make it positive
  if("dbh" %in% var) if(min(out$dbh) < 0) out$dbh <- out$dbh - min(out$dbh) + 0.01 
  colnames(out) <- colnames(data_model)
  return(out)
}






#' generate data for the jags model using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param p list containing the parameters for the harvest conditional probabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags_full_sub <- function(data){
  
  # Remove undisturbed plots
  data_sub <- data %>% filter(D == 1)
  
  # Replace plotcode by a number
  plotcode.table <- data.frame(plotcode = unique(data_sub$plotcode), 
                               plot = c(1:length(unique(data_sub$plotcode))))
  
  # Replace species by a number 
  species.table <- data.frame(species = unique(data_sub$species), 
                              sp = c(1:length(unique(data_sub$species))))
  
  # Identify other-disturbed spruce-dominated plots with 100% mortality
  plots.reference.other <- (data %>%
                              filter(Dother == 1) %>%
                              group_by(plotcode, species) %>%
                              mutate(dead = ifelse(a == 0, 1, 0)) %>%
                              summarize(n = n(), n.dead = sum(dead)) %>%
                              ungroup() %>%
                              group_by(plotcode) %>%
                              mutate(prop.species = n/sum(n), 
                                     prop.species.dead = n.dead/n) %>%
                              filter(species == "Picea abies" & prop.species == 1 & prop.species.dead == 1))$plotcode 
  
  # Identify storm-disturbed spruce-dominated plots with 100% mortality
  plots.reference.storm <- (data %>%
                              filter(Dstorm == 1) %>%
                              group_by(plotcode, species) %>%
                              mutate(dead = ifelse(a == 0, 1, 0)) %>%
                              summarize(n = n(), n.dead = sum(dead)) %>%
                              ungroup() %>%
                              group_by(plotcode) %>%
                              mutate(prop.species = n/sum(n), 
                                     prop.species.dead = n.dead/n) %>%
                              filter(species == "Picea abies" & prop.species > 0.8 & prop.species.dead == 1))$plotcode
  
  # Identify storm-disturbed spruce-dominated plots with 100% mortality
  plots.reference.fire <- (data %>%
                             filter(Dfire == 1) %>%
                             group_by(plotcode, species) %>%
                             mutate(dead = ifelse(a == 0, 1, 0)) %>%
                             summarize(n = n(), n.dead = sum(dead)) %>%
                             ungroup() %>%
                             group_by(plotcode) %>%
                             mutate(prop.species = n/sum(n), 
                                    prop.species.dead = n.dead/n) %>%
                             filter(species == "Pinus halepensis" & prop.species == 1 & 
                                      prop.species.dead == 1 & n > 25))$plotcode
  
  # Format final data
  data.in <- data_sub %>%
    # Set disturbance intensity for storm and other
    mutate(Ifire = ifelse(plotcode %in% plots.reference.fire, 0.7, NA_real_),
           Istorm = ifelse(plotcode %in% plots.reference.storm, 0.7, NA_real_),
           Iother = ifelse(plotcode %in% plots.reference.other, 0.7, NA_real_)) %>%
    # Determine if tree is dead or not
    mutate(d = ifelse(a == 1, 0, 1)) %>%
    # Add plot code
    left_join(plotcode.table, by = "plotcode") %>%
    # Add species code
    left_join(species.table, by = "species")
  
  
  
  # Create the output list
  data_jags <- list(
    Ntrees = NROW(data.in),
    Nplot = NROW(plotcode.table), 
    Nspecies = NROW(species.table),
    plot = data.in$plot, 
    sp = data.in$sp,
    d = data.in$d,
    time = data.in$time,
    dbh = data.in$dbh,
    Dfire = data.in$Dfire, 
    Dstorm = data.in$Dstorm, 
    Dother = data.in$Dother, 
    Ifire = data.in$Ifire, 
    Istorm = data.in$Istorm, 
    Iother = data.in$Iother)
  
  # Final output
  out <- list(data_jags, plotcode.table, species.table)
  names(out) = c("data_jags", "plotcode_table", "species_table")
  return(out)
}




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




#' Generate species-specific parameters to simulate data
#' @param data.jags.in dataframe containing plot and tree variable centered and scaled
#' @param param character vector containing the name of all the parameters to initialize
generate_parameters_sp <- function(data.jags.in, param = paste0("c", c(0:8))){
  out <- data.frame(sp = data.jags.in$species_table$sp)
  for(i in 1:length(param)){
    # intercept
    if(param[i] %in% c("c0", "c3", "c6")){
      out$value <- round(rnorm(dim(out)[1], mean = -4, sd = 0.5), digits = 1)
    } 
    # multiplicative parameter
    if(param[i] %in% c("c1", "c4", "c7")){
      out$value <- round(rnorm(dim(out)[1], mean = 3.5, sd = 0.5), digits = 1)
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
    mutate(Intensity = severity + rnorm(1, -0.05, 0.05),
           Intensity.corrected = case_when(Intensity < 0 ~ -Intensity, 
                                           Intensity > 1 ~ 1 - 2*(Intensity - 1), 
                                           TRUE ~ Intensity), 
           Ifire = Dfire*Intensity.corrected, 
           Istorm = Dstorm*Intensity.corrected, 
           Iother = Dother*Intensity.corrected) %>%
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