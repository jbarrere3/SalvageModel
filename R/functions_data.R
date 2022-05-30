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
    filter(dbh1 >= 100) %>%
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


#' Format the data before fiting IPM or mortality model, with mean quadratic diameter included
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plot FUNDIV plot table
#' @param Climate dataframe containing climatic data per plot
#' @param BM_equations dataframe indicating which background mortality to use per species
format_data_model_dqm <- function(FUNDIV_tree, FUNDIV_plot, Climate, species){
  
  ## - Begin formatting final dataset
  out <- FUNDIV_tree %>%
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
  species.aggregated <- out %>%
    group_by(species, group, plotcode) %>% 
    summarize(n.species.per.plotcode = n()) %>%
    ungroup() %>%
    group_by(species, group) %>%
    summarize(n.per.species = sum(n.species.per.plotcode), 
              n.plot.per.species = n()) %>%
    mutate(enough.individuals = ifelse((n.per.species > 400 & n.plot.per.species > 30), 1, 0),
           species.ag = case_when((enough.individuals == 0 & group == "Angiosperms") ~ "Other broadleaf", 
                                  (enough.individuals == 0 & group == "Gymnosperms") ~ "Other conifer",
                                  enough.individuals == 1  ~ species)) %>%
    dplyr::select(species, species.ag)
  
  ## - Continue formatting fo the final table
  out <- out %>%
    left_join(species.aggregated, by = "species") %>%
    # Compute basal area of competitors
    group_by(plotcode) %>%
    mutate(comp=sum(ba_ha1) - ba_ha1) %>%
    # Compute stand quadratic diameter
    filter(!is.na(weight1)) %>%
    mutate(dbh2.W = (dbh1^2)*weight1, 
           dqm = sqrt(sum(dbh2.W, na.rm = TRUE)/sum(weight1, na.rm  = TRUE))) %>%
    ungroup() %>%
    mutate(dominated = ifelse(dbh1 < dqm, 1, 0)) %>%
    # Reformat status columns
    mutate(h = case_when(treestatus == 3 ~ 1, T ~ 0), 
           d = case_when(treestatus %in% c(4, 5) ~ 1, T ~ 0), 
           a = case_when(treestatus == 2 ~ 1, T ~ 0)) %>%
    # Reformat disturbance columns
    mutate(disturbance.nature = ifelse(is.na(disturbance.nature), "other", disturbance.nature),
           Dfire = ifelse(disturbance.nature == "fire", 1, 0), 
           Dstorm = ifelse(disturbance.nature == "storm", 1, 0),
           Dother = ifelse(disturbance.nature == "other", 1, 0), 
           D = 1 - (1 - Dfire)*(1 - Dstorm)*(1 - Dother)) %>%
    filter(!is.na(sgdd) & !is.na(wai) & !is.na(comp) & !is.na(species.ag)) %>%
    # Select variables and remove NAs
    dplyr::select(plotcode, treecode, species = species.ag, time, h, d, a, dbh = dbh1, comp, 
                  sgdd, wai, dqm, D, Dfire, Dstorm, Dother, dominated) 
  
  return(out)
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

#' Center and Scale variables before fitting one fo the models that includes dqm
#' @param data_model Tree data formated for the IPM. 
#' @param var Variables to scale (character vector)
scale_data_model_dqm <- function(data_model, 
                             var = c("dbh", "dqm", "comp", "sgdd", "wai")){
  id_var <- which(colnames(data_model) %in% var) 
  if("DS" %in% var) out <- data_model %>% mutate(DS = ifelse(DS == 0, NA_real_, DS))
  else out <- data_model
  unscaled = as.matrix(out[, id_var])
  scaled = scale(unscaled)
  out[, id_var] <- as.data.frame(scaled)
  if("DS" %in% var) out$DS[which(is.na(out$DS))] <- -99
  # dbh and dqm has a power effect, so need to make them positive
  if("dbh" %in% var) if(min(out$dbh) < 0) out$dbh <- out$dbh - min(out$dbh) + 0.01 
  if("dqm" %in% var) if(min(out$dqm) < 0) out$dqm <- out$dqm - min(out$dqm) + 0.01 
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



#' generate data for the jags model using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @author Björn Reineking, Julien Barrere
generate_data_jags_full_sub_climate <- function(data){
  
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
    sgdd = data.in$sgdd, 
    wai = data.in$wai,
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



#' generate data for the jags model including dqm using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param p list containing the parameters for the harvest conditional probabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags_dqm <- function(data){
  
  # Replace plotcode by a number
  plotcode.table <- data.frame(plotcode = unique(data$plotcode), 
                               plot = c(1:length(unique(data$plotcode))))
  
  # Replace species by a number 
  species.table <- data.frame(species = unique(data$species), 
                              sp = c(1:length(unique(data$species))))
  
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
  data.in <- data %>%
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
    dqm = data.in$dqm,
    domi = data.in$dominated,
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


#' Generate species-specific parameters to simulate data
#' @param data.jags.in dataframe containing plot and tree variable centered and scaled
#' @param param character vector containing the name of all the parameters to initialize
generate_parameters_sp_climate <- function(data.jags.in, param = c(paste0("st", c(0:4)), paste0("ot", c(0:4)), paste0("fi", c(0:4)))){
  out <- data.frame(sp = data.jags.in$species_table$sp)
  for(i in 1:length(param)){
    # intercept
    if(param[i] %in% paste0(c("st", "ot", "fi"), 0)){
      out$value <- round(runif(dim(out)[1], -6, -4), digits = 3)
    } 
    # multiplicative parameter
    if(param[i] %in% paste0(c("st", "ot", "fi"), 1)){
      out$value <- round(runif(dim(out)[1], 3, 6), digits = 3)
    } 
    # power parameter
    if(param[i] %in% paste0(c("st", "ot", "fi"), 2)){
      out$value <- round(rnorm(dim(out)[1], mean = 0, sd = 0.5), digits = 1)
    } 
    # sgdd parameter
    if(param[i] %in% paste0(c("st", "ot", "fi"), 3)){
      out$value <- round(rnorm(dim(out)[1], mean = 0, sd = 1), digits = 1)
    } 
    # wai parameter
    if(param[i] %in% paste0(c("st", "ot", "fi"), 4)){
      out$value <- round(rnorm(dim(out)[1], mean = 0, sd = 1), digits = 1)
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


#' Simulate tree status data for the jags model with disturbance given and specific background
#' mortality equations (see supplementary of Kunstler et al. 2020 for original equations)
#' @param data.jags.in dataframe containing plot and tree variable centered and scaled
#' @param parameters_sp dataframe containing the parameters value per species
#' @author Björn Reineking, Julien Barrere
simulate_status_full_sub_climate <- function(data.jags.in, parameters_sp){
  
  
  
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
                                sgdd = data.jags.in$data_jags$sgdd, 
                                wai = data.jags.in$data_jags$wai, 
                                Dfire = data.jags.in$data_jags$Dfire, 
                                Dstorm = data.jags.in$data_jags$Dstorm, 
                                Dother = data.jags.in$data_jags$Dother) %>%
    # Add species parameter
    left_join(species.table.extended, by = "sp") %>%
    # Add intensity
    left_join(plot.table.extended, by = "plot") %>%
    # calculate probabilities for the three different states
    mutate(pdstorm = 1 - (1 - plogis(st0 + st1*Istorm*dbh^st2 + st3*sgdd + st4*wai))^time, 
           pdother = 1 - (1 - plogis(ot0 + ot1*Iother*dbh^ot2 + ot3*sgdd + ot4*wai))^time, 
           pdfire = 1 - (1 - plogis(fi0 + fi1*Ifire*dbh^fi2 + fi3*sgdd + fi4*wai))^time, 
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






#' Get info on which species are affected by which disturbance, 
#' and thus which parameters per species to keep as "robust"
#' @param data_model dataset formatted to fit the model
#' @return a list with two elements: 
#'         - species_disturbance_table: table giving the n of indiv and 
#'           plot per species affected by each type of disturbance
#'         - species_parameter_to_keep is a vector of species-parameter combination with enough
#'           plots and individuals
get_disturbance_species_info <- function(data_model){
  
  # Initialize output
  out <- list()
  
  # Table indicating the number of individuals per species affected by each type of disturbance
  out$species_disturbance_table <- data_model %>%
    gather(key = "disturbance.type", value = "dist.present", "Dfire", "Dother", "Dstorm") %>%
    filter(dist.present == 1) %>%
    mutate(disturbance.type = gsub("D", "", disturbance.type)) %>%
    group_by(disturbance.type, species, plotcode) %>%
    summarize(n.indiv.per.plot = n()) %>%
    group_by(disturbance.type, species) %>%
    summarize(n.indiv = sum(n.indiv.per.plot), 
              n.plot = n()) 
  
  # Vector of species - parameter combination to keep
  out$species_parameter_to_keep <- (out$species_disturbance_table %>%
                                      filter(n.indiv > 150 & n.plot > 15) %>%
                                      mutate(c0 = ifelse(disturbance.type == "storm", 1, 0), 
                                             c1 = ifelse(disturbance.type == "storm", 1, 0), 
                                             c2 = ifelse(disturbance.type == "storm", 1, 0), 
                                             c3 = ifelse(disturbance.type == "other", 1, 0), 
                                             c4 = ifelse(disturbance.type == "other", 1, 0), 
                                             c5 = ifelse(disturbance.type == "other", 1, 0), 
                                             c6 = ifelse(disturbance.type == "fire", 1, 0), 
                                             c7 = ifelse(disturbance.type == "fire", 1, 0), 
                                             c8 = ifelse(disturbance.type == "fire", 1, 0)) %>%
                                      gather(key = "Parameter", value = "Present", paste0("c", c(0:8))) %>%
                                      filter(Present == 1) %>%
                                      mutate(ID = paste(species, Parameter, sep = ".")))$ID
  
  return(out)
}


#' Get info on which species are affected by which disturbance, 
#' and thus which parameters per species to keep as "robust"
#' @param data_model dataset formatted to fit the model
#' @return a list with two elements: 
#'         - species_disturbance_table: table giving the n of indiv and 
#'           plot per species affected by each type of disturbance
#'         - species_parameter_to_keep is a vector of species-parameter combination with enough
#'           plots and individuals
get_disturbance_species_info_climate <- function(data_model, 
                                                 param.per.disturbance = list(storm = paste0("st", c(0:5)), 
                                                                              fire = paste0("fi", c(0:5)), 
                                                                              other = paste0("ot", c(0:5)))){
  
  # Initialize output
  out <- list()
  
  # Table indicating the number of individuals per species affected by each type of disturbance
  out$species_disturbance_table <- data_model %>%
    gather(key = "disturbance.type", value = "dist.present", "Dfire", "Dother", "Dstorm") %>%
    filter(dist.present == 1) %>%
    mutate(disturbance.type = gsub("D", "", disturbance.type)) %>%
    group_by(disturbance.type, species, plotcode) %>%
    summarize(n.indiv.per.plot = n()) %>%
    group_by(disturbance.type, species) %>%
    summarize(n.indiv = sum(n.indiv.per.plot), 
              n.plot = n()) 
  
  # Parameters per species per disturbance
  param.per.sp.per.dist <- expand.grid(disturbance = c("storm", "fire", "other"), 
                                       Parameter = as.character(unlist(param.per.disturbance)), 
                                       species = unique(data_model$species)) %>%
    mutate(value = case_when(disturbance == "storm" & Parameter %in% param.per.disturbance$storm ~ 1, 
                             disturbance == "other" & Parameter %in% param.per.disturbance$other ~ 1, 
                             disturbance == "fire" & Parameter %in% param.per.disturbance$fire ~ 1, 
                             TRUE ~ 0)) %>%
    spread(key = "Parameter", value = "value") %>%
    mutate(ID.spdist = paste(disturbance, species, sep = ".")) %>%
    dplyr::select("ID.spdist", as.character(unlist(param.per.disturbance)))
  
  
  
  
  # Vector of species - parameter combination to keep
  out$species_parameter_to_keep <- (out$species_disturbance_table %>%
                                      filter(n.indiv > 150 & n.plot > 15) %>%
                                      mutate(ID.spdist = paste(disturbance.type, species, sep = ".")) %>%
                                      left_join(param.per.sp.per.dist, by = "ID.spdist") %>%
                                      gather(key = "Parameter", value = "value", as.character(unlist(param.per.disturbance))) %>%
                                      filter(value == 1) %>%
                                      mutate(ID = paste(species, Parameter, sep = ".")))$ID
  
  return(out)
}



#' Extract from the model disturbance sensitivity
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param data_model_scaled Tree data formatted for the IPM and scaled 
get_disturbance_sensivity <- function(jags.model, data_jags, data_model_scaled, data_model, 
                                      disturbance_species_info){
  
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
  # Parameters to keep per species
  species_parameter_to_keep <- (disturbance_species_info$species_disturbance_table %>%
                                  filter(n.indiv > 50 & n.plot > 5) %>%
                                  mutate(c0 = ifelse(disturbance.type == "storm", 1, 0), 
                                         c1 = ifelse(disturbance.type == "storm", 1, 0), 
                                         c2 = ifelse(disturbance.type == "storm", 1, 0), 
                                         c3 = ifelse(disturbance.type == "other", 1, 0), 
                                         c4 = ifelse(disturbance.type == "other", 1, 0), 
                                         c5 = ifelse(disturbance.type == "other", 1, 0), 
                                         c6 = ifelse(disturbance.type == "fire", 1, 0), 
                                         c7 = ifelse(disturbance.type == "fire", 1, 0), 
                                         c8 = ifelse(disturbance.type == "fire", 1, 0)) %>%
                                  gather(key = "Parameter", value = "Present", paste0("c", c(0:8))) %>%
                                  filter(Present == 1) %>%
                                  mutate(ID = paste(species, Parameter, sep = ".")))$ID
  
  # - Pre-format the data before plotting
  out <- expand.grid(species = unique(param_per_species$species), 
                     dbh = c(150, 400, 700), 
                     Intensity = 0.8) %>%
    mutate(dbh.scaled = predict(scale_dbh, newdata = .)) %>%
    left_join(param_per_species, by = "species") %>%
    gather(key = "Parameter", value = "value", paste0("c", c(0:8))) %>%
    mutate(ID = paste(species, Parameter, sep = "."), 
           value = ifelse(ID %in% species_parameter_to_keep, value, NA_real_)) %>%
    dplyr::select(-ID) %>%
    spread(key = Parameter, value = value) %>%
    mutate(storm = plogis(c0 + c1*Intensity*dbh.scaled^c2), 
           other = plogis(c3 + c4*Intensity*dbh.scaled^c5), 
           fire = plogis(c6 + c7*Intensity*dbh.scaled^c8)) %>%
    dplyr::select(species, dbh, storm, other, fire) %>%
    gather(key = "disturbance", value = "sensitivity", "storm", "other", "fire") %>%
    mutate(trait.name = paste0(disturbance, ".dbh", dbh)) %>%
    dplyr::select(-dbh, -disturbance) %>%
    spread(key = trait.name, value = sensitivity)
  
  
  # return the formatted dataset
  return(out)
}





#' Compile all traits data
#' @param bark.thickness_file file containing data on bark thickness
#' @param wood.density_file file containing data on wood density
#' @param species.in character vector of all the species for which to extract the data
compile_traits <- function(wood.density_file, shade.tolerance_file, root.depth_file, species.in){
  
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
  
  # Global trait dataset
  traits <- data.frame(species = species.in) %>%
    # Add wood density
    left_join((wood.density %>% 
                 group_by(species) %>%
                 summarize(wood.density_g.cm3 = mean(wood.density_g.cm3))),
              by = "species") %>%
    # Add shade tolerance
    left_join(shade.tolerance, by = "species") %>%
    # Add root depth
    left_join(root.depth, by = "species")
  
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


