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


#' Get parameters harvest probabilities
#' @param data_model_scaled
#' @param Dj_latent
get_param_harvest_proba <- function(data_model_scaled, Dj_latent){
  
  # Initialize output
  out <- list()
  
  # Use the right DS depending on Dj_latent
  if(Dj_latent) data_priorharvest <- data_model_scaled %>% dplyr::select(plotcode, dbh, a, h, D, DS = DS.Senf)
  else data_priorharvest <- data_model_scaled %>% dplyr::select(plotcode, dbh, a, h, D, DS)
  
  # Remove alive trees in disturbed plots
  data_priorharvest <- data_priorharvest %>%
    filter(!(a == 1 & D == 1)) %>% 
    dplyr::select(-a)
  
  
  # Fit model for undisturbed plots
  mod.undist <- glm(h ~ dbh, 
                    data = subset(data_priorharvest, D == 0), 
                    family = binomial(link = "logit"))
  # Get parameters e0 and e1 (cf. doc Multinomial harvest model)
  out$e0 <- as.numeric(mod.undist$coefficients[1])
  out$e1 <- as.numeric(mod.undist$coefficients[2])
  
  # Fit model for disturbed plots
  mod.dist <- glm(h ~ dbh*DS, 
                  data = subset(data_priorharvest, D == 1), 
                  family = binomial(link = "logit"))
  
  # Get parameters d0, d1, d2 and d3 (cf. doc Multinomial harvest model)
  out$d0 <- as.numeric(mod.dist$coefficients[1])
  out$d1 <- as.numeric(mod.dist$coefficients[2])
  out$d2 <- as.numeric(mod.dist$coefficients[3])
  out$d3 <- as.numeric(mod.dist$coefficients[4])
  
  return(out)
}

#' generate data for the jags model using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param p list containing the parameters for the harvest conditional probabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags <- function(data, p = param_harvest_proba){
  
  # Replace plotcode by a number
  plotcode.table <- data.frame(plotcode = unique(data$plotcode), 
                               plot = c(1:length(unique(data$plotcode))))
  
  # Replace species by a number 
  species.table <- data.frame(species = unique(data$species), 
                              sp = c(1:length(unique(data$species))))
  
  # Format final data
  data.in <- data %>%
    # Calculate harvest conditional probabilities
    mutate(phdD = plogis(p$d0 + p$d1*dbh + p$d2*DS + p$d3*dbh*DS), 
           phadBM = plogis(p$e0 + p$e1*dbh)) %>%
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
    state = apply(dplyr::select(data.in, h, d, a), 1, which.max), 
    dbh = data.in$dbh,
    comp = data.in$comp,
    sgdd = data.in$sgdd,
    wai = data.in$wai,
    phdD = data.in$phdD,
    phadBM = data.in$phadBM, 
    D = data.in$D, 
    Dfire = data.in$Dfire, 
    Dstorm = data.in$Dstorm, 
    Dother = data.in$Dother,
    bm1 = data.in$bm1,
    bm2 = data.in$bm2, 
    bm3 = data.in$bm3, 
    bm4 = data.in$bm4, 
    bm5 = data.in$bm5, 
    bm6 = data.in$bm6
  )
  
  # Final output
  out <- list(data_jags, plotcode.table, species.table)
  names(out) = c("data_jags", "plotcode_table", "species_table")
  return(out)
}





#' generate data for the jags model using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param p list containing the parameters for the harvest conditional probabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags_sub <- function(data, p = param_harvest_proba){
  
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
  
  # Format final data
  data.in <- data_sub %>%
    # Calculate harvest conditional probabilities
    mutate(phdD = plogis(p$d0 + p$d1*dbh + p$d2*DS + p$d3*dbh*DS)) %>%
    # Remove (for now) fire disturbance, and only keep storm and other (fire classified as other)
    mutate(Dother = ifelse(Dfire == 1, 1, Dother)) %>%
    # Set disturbance intensity for storm and other
    mutate(Istorm = ifelse(plotcode %in% plots.reference.storm, 0.7, NA_real_),
           Iother = ifelse(plotcode %in% plots.reference.other, 0.7, NA_real_)) %>%
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
    state = apply(dplyr::select(data.in, h, d, a), 1, which.max),
    time = data.in$time,
    dbh = data.in$dbh,
    phdD = data.in$phdD,
    Dstorm = data.in$Dstorm, 
    Dother = data.in$Dother, 
    Istorm = data.in$Istorm, 
    Iother = data.in$Iother)
  
  # Final output
  out <- list(data_jags, plotcode.table, species.table)
  names(out) = c("data_jags", "plotcode_table", "species_table")
  return(out)
}




#' generate data for the jags model using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param p list containing the parameters for the harvest conditional probabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags_full_sub <- function(data, p = param_harvest_proba){
  
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
                             filter(species == "Pinus halepensis" & prop.species > 0.8 & 
                                      prop.species.dead == 1 & n > 20))$plotcode
  
  # Format final data
  data.in <- data_sub %>%
    # Calculate harvest conditional probabilities
    mutate(phdD = plogis(p$d0 + p$d1*dbh + p$d2*DS + p$d3*dbh*DS)) %>%
    # Set disturbance intensity for storm and other
    mutate(Ifire = ifelse(plotcode %in% plots.reference.fire, 0.7, NA_real_),
           Istorm = ifelse(plotcode %in% plots.reference.storm, 0.7, NA_real_),
           Iother = ifelse(plotcode %in% plots.reference.other, 0.7, NA_real_)) %>%
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
    state = apply(dplyr::select(data.in, h, d, a), 1, which.max),
    time = data.in$time,
    dbh = data.in$dbh,
    phdD = data.in$phdD,
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
