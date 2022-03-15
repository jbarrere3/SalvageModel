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


#' Get file from its url and write it on disk, at a specified location. 
#' @param dir.name.in Directory where the file should be written (ex: "data/BWI")
#' @param url.in URL where to download the file.
get_and_write <- function(dir.name.in, url.in){
  
  # Write directories if they do not exist
  path.in <- strsplit(dir.name.in, "/")[[1]]
  for(i in 1:length(path.in)){
    if(i == 1) path.in_i <- path.in[i]
    else path.in_i <- paste(path.in_i, path.in[i], sep = "/")
    if(!dir.exists(path.in_i)) dir.create(path.in_i)
  }
  
  # Write file on the disk
  url.in_split <- strsplit(url.in, "/")[[1]]
  file.in <- paste(dir.name.in, url.in_split[length(url.in_split)], sep = "/")
  if(!file.exists(file.in)){
    try(GET(url.in, authenticate('guest', ""), write_disk(file.in, overwrite = TRUE)))
    # Specific case of zip file: unzip and delete zip file
    if("zip" %in% strsplit(file.in, split = "\\.")[[1]]){
      unzip(file.in, exdir = dir.name.in, overwrite = T)
      print(paste0("---Getting and unzipping ", file.in))
      unlink(file.in)
    }else{print(paste0("---Getting ", file.in))}
  } 
}


#' Function to write all FRENCH NFI data on disk
#' @return a dataframe containing the directories where each file is stored, and the url to get these files
get_FrenchNFI <- function(){
  
  dir.in = c("data/FrenchNFI")
  url.in = c("https://inventaire-forestier.ign.fr/dataifn/data/export_dataifn_2005_2020.zip")
  get_and_write(dir.in, url.in)
  return(paste(dir.in, list.files(dir.in), sep = "/"))
}


#' Format the raw FrenchNFI data
#' @author Natheo BEAUCHAMP (rearranged by Julien BARRERE)
#' @param FrenchNFI_tree_raw Raw french NFI tree table
#' @param FrenchNFI_plot_raw Raw french NFI plot table
#' @return A list with two elements: formated tree table (FrenchNFI_tree) and plot table (FrenchNFI_plot)
format_FrenchNFI_raw <- function(FrenchNFI_tree_raw, FrenchNFI_plot_raw){
  ## - Format plot level data
  FrenchNFI_plot <- FrenchNFI_plot_raw %>% 
    # FILTER BY NUMBER OF VISITS
    group_by(IDP) %>% 
    mutate(Nvisits = n(),
           lastCampagne = max(CAMPAGNE)) %>% 
    filter(Nvisits == 2) %>% 
    mutate(firstVisit = !(CAMPAGNE == lastCampagne)) %>% 
    # FILTER TO REMOVE FIRST CAMPAIGN BECAUSE NO REMEASURENT OF C13
    mutate(firstCampaign = (CAMPAGNE == 2005 & firstVisit) | (CAMPAGNE == 2010 & !firstVisit)) %>% 
    filter(!firstCampaign) %>% 
    # FILTER BY VEGETATION TYPE
    mutate(isForest = CSA %in% c(1,2,3)) %>% 
    group_by(IDP) %>% 
    mutate(isForestBothVisits = (sum(isForest)==2)) %>%
    filter(isForestBothVisits) %>%
    # FILTER BY INTERCEPTION WITH FOREST EDGE
    # Set modalities of before-2007 campaign to fit those of post-2007 
    mutate(PLISI = case_when(
      CAMPAGNE > 2007  ~ ifelse(is.na(PLISI), 1, PLISI + 1), 
      CAMPAGNE <= 2007 ~ ifelse(is.na(PLISI), 1, PLISI + 0) )) %>%
    mutate(isNotAtEdge = (PLISI==1)) %>% 
    group_by(IDP) %>% 
    mutate(isNotAtEdgeBothVisits = (sum(isNotAtEdge)==2))
  
  
  ## FILTER THE TREE DATASET 
  
  # Filter tree dataset to keep only filtered plots
  FrenchNFI_tree <- FrenchNFI_tree_raw %>% 
    filter(IDP %in% unique(FrenchNFI_plot$IDP))
  
  
  # Get IDP of plots that are visited twice
  plots_visited_twice <- FrenchNFI_tree %>% 
    group_by(IDP) %>% 
    summarize(Nvisits = length(unique((CAMPAGNE)))) %>% 
    filter(Nvisits == 2)
  plots_visited_twice <- plots_visited_twice$IDP
  
  # Filter for keeping only plots visited twice in tree and stand datasets
  FrenchNFI_tree <- FrenchNFI_tree %>% 
    filter(IDP %in% plots_visited_twice)
  FrenchNFI_plot <- FrenchNFI_plot %>% 
    filter(IDP %in% plots_visited_twice)
  
  
  
  
  # FILTER THE TREE DATASET FOR KEEPING ONLY STANDS WHERE TREES HAS BEES REMEASURED AT THE SECOND VISIT
  FrenchNFI_tree <- FrenchNFI_tree %>% 
    left_join(FrenchNFI_plot %>% select(IDP, CAMPAGNE, firstVisit),
              by = c("IDP", "CAMPAGNE")) %>% 
    mutate(isNotRemeasured = ifelse(!firstVisit & VEGET5 == 0, is.na(C13), NA))
  # IDs of remeasured stands
  FrenchNFI_plot_remeasured <- (FrenchNFI_tree %>% 
                                  group_by(IDP) %>% 
                                  summarize(SumNotRemeasured = sum(isNotRemeasured, na.rm = TRUE)) %>% 
                                  mutate(standsRemeasured = (SumNotRemeasured==0)) %>% 
                                  filter(standsRemeasured))$IDP
  # Filter datasets with only remeasured stands
  FrenchNFI_tree <- FrenchNFI_tree %>%
    filter(IDP %in% FrenchNFI_plot_remeasured)
  FrenchNFI_plot <- FrenchNFI_plot %>% 
    filter(IDP %in% FrenchNFI_plot_remeasured)
  
  
  
  # FILTER THE TREE DATASET FOR KEEPING ONLY LIVING TREES AT FIRST VISIT
  # Searching for IDP, A with VEGET != 0
  FrenchNFI_tree <- FrenchNFI_tree %>%
    mutate(isVeget0 = ifelse(firstVisit, VEGET==0, NA)) %>% 
    group_by(IDP, A) %>% 
    mutate(firstVisitIsVeget0 = sum(isVeget0, na.rm = TRUE)==1) %>% 
    filter(firstVisitIsVeget0) 
  FrenchNFI_plot_livingtrees <- unique(FrenchNFI_tree$IDP)
  FrenchNFI_plot <- FrenchNFI_plot %>% 
    filter(IDP %in% FrenchNFI_plot_livingtrees)
  
  out <- list(FrenchNFI_tree, FrenchNFI_plot)
  return(out)
}



#' Format French NFI tree data to FUNIV original format
#' @author Natheo BEAUCHAMP (rearranged by Julien BARRERE)
#' @param FrenchNFI french NFI data formatted. 
#' @param FrenchNFI_species FrenchNFI species table
format_FrenchNFI_tree_to_FUNDIV <- function(FrenchNFI, FrenchNFI_species){
  FrenchNFI[[1]] %>% 
    group_by() %>% 
    select(CAMPAGNE, IDP, firstVisit, A, ESPAR, C13, HTOT, VEGET, VEGET5, W) %>% 
    pivot_wider(names_from = firstVisit,
                names_sep = "_",
                values_from = c(CAMPAGNE, ESPAR, C13, HTOT, VEGET, VEGET5, W)) %>% 
    mutate(country = "France") %>%
    rename(plotcode = IDP) %>%
    mutate(treecode = paste(plotcode, A, sep = "_")) %>%
    dplyr::select(country, plotcode, treecode, VEGET5 = VEGET5_FALSE,
                  c13_1 = C13_TRUE, c13_2 = C13_FALSE, weight1 = W_TRUE, 
                  height1 = HTOT_TRUE, height2 = HTOT_FALSE, sp_code = ESPAR_TRUE) %>% 
    mutate(dbh1 = c13_1/(pi) * 1000,
           dbh2 = c13_2/(pi) * 1000) %>% 
    filter(!(dbh1 < 100 & dbh2 < 100)) %>%
    mutate(treestatus = case_when((VEGET5 == "0" & dbh1 >= 100) ~ 2,
                                  (VEGET5 == "0" & dbh1 < 100) ~ 1,
                                  VEGET5 %in% c("6", "7") ~ 3,
                                  VEGET5 %in% c("1", "2", "5", "A", "C", "M", "T") ~ 4,
                                  VEGET5 == "N" ~ 5), 
      ba1 = (pi*((dbh1/1000)/2)^2),
      ba2 = (pi*((dbh2/1000)/2)^2), 
      weight2 = NA_real_) %>% 
    mutate(ba_ha1 = ba1*weight1,
           ba_ha2 = ba2*weight1,
           bachange_ha_yr = weight1*(ba2 - ba1)/5) %>% 
    left_join(FrenchNFI_species %>%
                dplyr::select(code, latinName), 
              by = c("sp_code" = "code")) %>%
    dplyr::select(treecode, plotcode, species = latinName, treestatus, dbh1, dbh2, 
                  height1, height2, ba1, ba_ha1, ba2, ba_ha2, bachange_ha_yr, 
                  weight1, weight2, country)
}



#' Format French NFI plot data to FUNIV original format
#' @author Natheo BEAUCHAMP (rearranged by Julien BARRERE)
#' @param FrenchNFI french NFI data formatted. 
#' @param FUNDIV_tree_FR_raw FrenchNFI tree table formatted for FUNDIV
format_FrenchNFI_plot_to_FUNDIV <- function(FrenchNFI, FUNDIV_tree_FR_raw){
  # Select interest variables
  out <- FrenchNFI[[2]] %>% 
    group_by() %>% 
    dplyr::select(IDP, CAMPAGNE, firstVisit, XL, YL, INCID, NINCID) %>% 
    pivot_wider(names_from = firstVisit,
                names_sep = "_",
                values_from = c(CAMPAGNE, XL, YL, INCID, NINCID))  %>% 
    mutate(XL_SAME = (XL_TRUE==XL_FALSE),
           YL_SAME = (YL_TRUE==YL_FALSE)) %>% 
    dplyr::select(IDP, surveydate1 = CAMPAGNE_TRUE, surveydate2 = CAMPAGNE_FALSE,
                  XL = XL_TRUE, YL = YL_TRUE, incid = INCID_FALSE, nincid = NINCID_FALSE) %>% 
    mutate(country = "France", 
           yearsbetweensurveys = 5, 
           biome = NA_real_, 
           cluster = NA_real_) %>%
    mutate(disturbance.severity = case_when(nincid == 0 ~ 0, 
                                            TRUE ~ as.numeric(incid)), 
           disturbance.nature = case_when(nincid == 0 ~ "none", 
                                          nincid == 1 ~ "fire", 
                                          nincid == 4 ~ "storm", 
                                          nincid %in% c(2, 3, 5) ~ "other", 
                                          is.na(nincid) ~ if_else(incid == 0, "none", NA_character_))) %>%
    rename(plotcode = IDP)
  
  # Convert coord from Lambert93 French national projection system to WGS84 projection system
  out.coordinates <- out %>% dplyr::select(XL, YL)
  coordinates(out.coordinates) <- c("XL", "YL")
  proj4string(out.coordinates) <- CRS("+init=epsg:2154") 
  out.coordinates <- spTransform(out.coordinates,  CRS("+init=epsg:4326"))
  out$longitude <- coordinates(out.coordinates)[,"XL"]
  out$latitude <- coordinates(out.coordinates)[,"YL"]
  
  # Get management and basal area per hectare from tree level data
  out.management.weight <- FUNDIV_tree_FR_raw %>%
    mutate(harvest = case_when(treestatus == 3 ~ 1, TRUE ~ 0)) %>%
    group_by(plotcode) %>%
    summarise(ba_ha1 = sum(ba_ha1, na.rm = T),
              ba_ha2 = sum(ba_ha2, na.rm = T), 
              n.harvest = sum(harvest, na.rm = T)) %>%
    mutate(management = case_when(n.harvest > 0 ~ 1, TRUE ~ 0)) %>%
    dplyr::select(plotcode, ba_ha1, ba_ha2, management)
  
  
  # Finish formating
  out <- out %>%
    left_join(out.management.weight, by = "plotcode") %>%
    dplyr::select(plotcode, cluster, country, longitude, latitude, 
                  yearsbetweensurveys, surveydate1, surveydate2, biome, 
                  ba_ha1, ba_ha2, management, disturbance.severity, disturbance.nature)
  
  return(out)
}




#' Format the data before fiting IPM or mortality model
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plot FUNDIV plot table
#' @param Disturbance dataframe containing Senf disturbance area and severity per plot
#' @param Climate dataframe containing climatic data per plot
format_data_model <- function(FUNDIV_tree, FUNDIV_plot, Disturbance, Climate){
  FUNDIV_tree %>%
    # Add plot level data 
    merge(Disturbance, by = "plotcode", all.x = T, all.y = F) %>%
    merge(Climate, by = "plotcode", all.x = T, all.y = F) %>% 
    merge((FUNDIV_plot %>% dplyr::select(plotcode, disturbance.severity, disturbance.nature)), 
          by = "plotcode", all.x = T, all.y = F) %>% 
    group_by(plotcode) %>%
    # Compute basal area of competitors
    mutate(ba_ha1_plot=sum(ba_ha1)) %>%
    ungroup() %>%
    mutate(comp = ba_ha1_plot - ba_ha1, 
           h = case_when(treestatus == 3 ~ 1, T ~ 0), 
           d = case_when(treestatus %in% c(4, 5) ~ 1, T ~ 0), 
           a = case_when(treestatus == 2 ~ 1, T ~ 0), 
           D = case_when(disturbance.nature == "none" ~ 0, T ~ 1),
           DS = case_when(disturbance.severity == 0 ~ 0, 
                          disturbance.severity == 1 ~ 0.125, 
                          disturbance.severity == 2 ~ 0.375, 
                          disturbance.severity == 3 ~ 0.625, 
                          disturbance.severity == 4 ~ 0.875)) %>%
    filter(!is.na(sgdd) & !is.na(wai) & !is.na(Senf.DS)& !is.na(disturbance.severity)) %>%
    # Select variables and remove NAs
    dplyr::select(plotcode, treecode, species, h, d, a, dbh = dbh1, comp, sgdd, 
                  wai, DA.Senf = Senf.DA, DS.Senf = Senf.DS, D, DS) 
}


#' Center and Scale variables before fitting a model
#' @param data_model Tree data formated for the IPM. 
#' @param var Variables to scale (character vector)
scale_data_model <- function(data_model, 
                             var = c("dbh", "comp", "sgdd", "wai", "DA.Senf", "DS.Senf", "DS")){
  id_var <- which(colnames(data_model) %in% var) 
  out <- data_model
  unscaled = as.matrix(out[, id_var])
  scaled = scale(unscaled)
  out[, id_var] <- as.data.frame(scaled)
  colnames(out) <- colnames(data_model)
  return(out)
}


#' Get parameters harvest probabilities
#' @param data_model_scaled
get_param_harvest_proba <- function(data_model_scaled){
  
  # Initialize output
  out <- list()
  
  # Format data
  data_priorharvest <- data_model_scaled %>%
    dplyr::select(plotcode, dbh, a, h, D, DS) %>%
    filter(!(a == 1 & D == 1)) %>% # Remove alive trees in disturbed plots
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



#' Get parameters disturbance probabilities
#' @param data_model_scaled
get_param_disturbance_proba <- function(data_model_scaled){
  
  # Initialize output
  out <- list()
  
  # Format data
  data.in <- data_model_scaled %>%
    dplyr::select(plotcode, D, DA = DA.Senf) %>%
    distinct()
  
  # Fit model
  mod <- glm(D ~ DA, data = data.in, family = binomial(link = "logit"))
  # Get parameters e0 and e1 (cf. doc Multinomial harvest model)
  out$a0 <- as.numeric(mod$coefficients[1])
  out$a1 <- as.numeric(mod$coefficients[2])
  
  return(out)
}

#' Simulate tree status data for the jags model
#' @param data.in dataframe containing plot and tree variable centered and scaled
#' @param param_harvest_proba list containing the parameters for the harvest conditional proobabilities
#' @param Dj_latent Boolean specifying if the variable Dj (occurence of a disturbance) comes from 
#'                  the data or should be estimated
#' @param par list containing the parameters value used to simulate the data
#' @author Björn Reineking, Julien Barrere
simulate_status <- function(data.in, param_harvest_proba, Dj_latent, 
                            par = list(a0 = 0, a1 = 3, b0 = 0, b1 = -3, b2 = 3, 
                                       b3 = -3, b4 = 3, c0 = 0, c1 = 3, c2 = -3)) {
  
  # Add simulated disturbance events
  plot_data <- data.in %>%
    dplyr::select(plotcode, DA = DA.Senf, DS = DS.Senf, D) %>% distinct(plotcode, .keep_all = TRUE) 
  if(Dj_latent) {
    plot_data <- plot_data %>% mutate(D = rbinom(NROW(plot_data), size = 1, prob = plogis(par$a0 + par$a1 * DA)))
  }
  tree_data <- left_join(dplyr::select(data.in, -D, -DS), plot_data, by = "plotcode")
  
  d0 = param_harvest_proba$d0
  d1 = param_harvest_proba$d1
  d2 = param_harvest_proba$d2
  d3 = param_harvest_proba$d3
  e0 = param_harvest_proba$e0
  e1 = param_harvest_proba$e1
  
  # Compute probabilities based on parameters and data
  tree_data <- tree_data %>%
    mutate(pdD = plogis(par$c0 + par$c1*DS + par$c2*DS*dbh), 
           pdBM = plogis(par$b0 + par$b1*dbh + par$b2*comp + par$b3*sgdd + par$b4*wai), 
           phdD = plogis(d0 + d1*dbh + d2*DS + d3*dbh*DS), 
           phadBM = plogis(e0 + e1*dbh),
           pdDj = 1 - (1 - pdD)*(1 - pdBM),
           ph = (1 - D)*phadBM + D*pdDj*phdD, 
           pd = (1 - D)*pdBM*(1 - phadBM) + D*pdDj*(1 - phdD),
           pa = 1 - (ph + pd))
  
  # Multinational realization of tree state - we have to call rmultinom for each tree separately
  hda <- t(apply(tree_data[, c("ph", "pd", "pa")], 1, function(x) rmultinom(n = 1, size = 1, prob = x)))
  colnames(hda) <- c("h", "d", "a")
  tree_data <- cbind(dplyr::select(tree_data, -h, -d, -a), hda)
  tree_data
}


#' generate data for the jags model using simulated data
#' @param data_simulated dataset where the tree status (dead, alive or harvested) is simulated
#' @param Dj_latent Boolean specifying if the variable Dj (occurence of a disturbance) comes from 
#'                  the data or should be estimated
#' @author Björn Reineking, Julien Barrere
generate_data_jags_from_simulated <- function(data_simulated, Dj_latent) {
  if (Dj_latent) {
    # separate data at plot and tree level
    plot_data <- data_simulated %>%
      dplyr::select(plotcode, D, DA = DA.Senf, DS = DS.Senf) %>% distinct(plotcode, .keep_all = TRUE)
    
    # Add variable "plot" from 1:n_plots, and sort plot data by this variable.
    # This is important because we will use the variable "plot" of a given tree as row number to access 
    # the corresponding plot data in jags
    plot_data <- plot_data %>% mutate(plot = as.integer(as.factor(plotcode))) %>% arrange(plot)
    # Tree data
    data_simulated <- left_join(dplyr::select(data_simulated, -D, -DS), plot_data, by = "plotcode")
  }
  data_jags <- list(
    Ntrees = NROW(data_simulated),
    dbh = data_simulated$dbh,
    comp = data_simulated$comp,
    sgdd = data_simulated$sgdd,
    wai = data_simulated$wai,
    DS = data_simulated$DS.Senf,
    phdD = data_simulated$phdD,
    phadBM = data_simulated$phadBM,
    state = apply(dplyr::select(data_simulated, h, d, a), 1, which.max)
  )
  
  if (Dj_latent) {
    data_jags$Nplots <- NROW(plot_data)
    data_jags$DA <- plot_data$DA
    data_jags$plot <- data_simulated$plot
  } else {
    data_jags$D <- data_simulated$D
  }
  data_jags
}

#' generate data for the jags model using real data
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param param_harvest_proba list containing the parameters for the harvest conditional proobabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags <- function(data, param_harvest_proba){
  # Add harvest conditional probabilities
  d0 = param_harvest_proba$d0
  d1 = param_harvest_proba$d1
  d2 = param_harvest_proba$d2
  d3 = param_harvest_proba$d3
  e0 = param_harvest_proba$e0
  e1 = param_harvest_proba$e1
  data.in <- data %>%
    mutate(phdD = plogis(d0 + d1*dbh + d2*DS + d3*dbh*DS), 
           phadBM = plogis(e0 + e1*dbh))
  
  # Separate data at the plot and tree level
  plot_data <- data.in %>%
    dplyr::select(plotcode, D, DA = DA.Senf, DS = DS.Senf) %>% distinct(plotcode, .keep_all = TRUE)
  # Add variable "plot" from 1:n_plots, and sort plot data by this variable.
  # This is important because we will use the variable "plot" of a given tree as row number to access 
  # the corresponding plot data in jags
  plot_data <- plot_data %>% mutate(plot = as.integer(as.factor(plotcode))) %>% arrange(plot)
  # Tree data
  data.in <- left_join(dplyr::select(data.in, -D, -DS), plot_data, by = "plotcode")
  
  # Create the output list
  data_jags <- list(
    Ntrees = NROW(data.in),
    dbh = data.in$dbh,
    comp = data.in$comp,
    sgdd = data.in$sgdd,
    wai = data.in$wai,
    DS = data.in$DS.Senf,
    phdD = data.in$phdD,
    phadBM = data.in$phadBM,
    state = apply(dplyr::select(data.in, h, d, a), 1, which.max), 
    Nplots = NROW(plot_data), 
    DA = plot_data$DA, 
    plot = data.in$plot
  )
  
  return(data_jags)
}


#' Add probability that a disturbance occurs to data_jags
#' @param data_jags input data for the jags model (must have been generated with Dj_latent = TRUE)
#' @param a0 numeric value for a0 parameter (intercept for disturbance occurrence)
#' @param a1 numeric value for a1 parameter (effect on DA on disturbance occurrence)
add_pD_to_data_jags <- function(data_jags, a0, a1){
  out <- data_jags
  out$pD <- plogis(a0 + a1*out$DA)
  out <- out[-which(names(out) == "DA")]
  return(out)
}


#' generate data for the jags model using real data where disturbance occurrence and severity is given
#' @param data dataset where the tree status (dead, alive or harvested) is not simulated
#' @param param_harvest_proba list containing the parameters for the harvest conditional proobabilities
#' @author Björn Reineking, Julien Barrere
generate_data_jags_D <- function(data, param_harvest_proba){
  # Add harvest conditional probabilities
  d0 = param_harvest_proba$d0
  d1 = param_harvest_proba$d1
  d2 = param_harvest_proba$d2
  d3 = param_harvest_proba$d3
  e0 = param_harvest_proba$e0
  e1 = param_harvest_proba$e1
  data.in <- data %>%
    mutate(phdD = plogis(d0 + d1*dbh + d2*DS + d3*dbh*DS), 
           phadBM = plogis(e0 + e1*dbh))
  
  # Create the output list
  data_jags <- list(
    Ntrees = NROW(data.in),
    dbh = data.in$dbh,
    comp = data.in$comp,
    sgdd = data.in$sgdd,
    wai = data.in$wai,
    DS = data.in$DS.Senf,
    phdD = data.in$phdD,
    phadBM = data.in$phadBM,
    state = apply(dplyr::select(data.in, h, d, a), 1, which.max), 
    D = data.in$D
  )
  
  return(data_jags)
}