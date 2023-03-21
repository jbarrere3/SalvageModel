#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_mixture.R  
#' @description R script containing additional functions to fit a model 
#'              with a mixture effect 
#' @author Julien BARRERE
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Format the data before fitting the reference mortality model
#' @param FUNDIV_tree FUNDIV tree table
#' @param FUNDIV_plot FUNDIV plot table
#' @param Climate dataframe containing climatic data per plot
#' @param species dataframe containing species information
format_data_model_mix <- function(FUNDIV_tree, FUNDIV_plot, Climate, species){
  
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
    # Remove empty species 
    filter(species != " ") %>%
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
      # Join species aggregated table
      left_join((species.aggregated %>% 
                   dplyr::select("species", "species.ag" = "biotic")), 
                by = "species") %>%
      # Compute basal area of competitors
      group_by(plotcode) %>%
      mutate(comp=sum(ba_ha1) - ba_ha1) %>%
      # Compute share of mixture
      mutate(conifer = ifelse(group == "Angiosperms", 0, 1), 
             broadleaf = ifelse(group == "Gymnosperms", 0, 1), 
             share_conifer = sum(ba_ha1*conifer)/sum(ba_ha1), 
             share_broadleaf = sum(ba_ha1*broadleaf)/sum(ba_ha1), 
             mix = ifelse(conifer == 1, share_broadleaf, share_conifer), 
             mixed = ifelse(mix %in% c(0, 1), 0, 1),
             mix = ifelse(mix < 0.001, 0.001, mix),
             mix_logit = log(mix/(1 - mix))) %>%
      dplyr::select(-conifer, -broadleaf, -share_conifer, -share_broadleaf) %>%
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
      dplyr::select(plotcode, treecode, country, species = species.ag, time, h, d, 
                    a, dbh = dbh1, comp, sgdd, wai, dqm, log_dbh_dqm, stock, mix, 
                    mix_logit, mixed)
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- data.i")))
    
  }
  
  # return output list
  return(out)
}



#' generate data for the jags model
#' @param data_model data set formatted to fit mixture effect model
generate_data_jags_mix <- function(data_model){
  
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
      logratio = data.out.i$log_dbh_dqm, 
      mix_logit = data.out.i$mix_logit, 
      mixed = data.out.i$mixed)
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


#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_biotic <- function(data_jags, n.chains, n.iter, n.burn, n.thin, 
                                 param = c("a0", "a1", "b", "c", "I")){
  
  # Initialize time
  start <- Sys.time()
  
  ## - For prior of intensity, fit beta distribution on severity
  # Prepare numeric vector containing severity values
  DS_distr <- (data.frame(plot = data_jags$plot, 
                          dead = data_jags$d) %>%
                 group_by(plot) %>%
                 summarise(DS = sum(dead)/n()) %>%
                 mutate(DS = case_when(DS < 0.001 ~ 0.001, 
                                       DS > 0.999 ~ 0.999, 
                                       TRUE ~ DS)))$DS
  # Fit distribution based on this vector, and get parameters
  beta1 <- as.numeric(fitdistr(DS_distr, densfun = "beta", 
                               start = list(shape1 = 1, shape2 = 2))[[1]][1])
  beta2 <- as.numeric(fitdistr(DS_distr, densfun = "beta", 
                               start = list(shape1 = 1, shape2 = 2))[[1]][2])
  
  ## - Write the model
  mortality_model <- paste0(
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pd[i])
      # Probability to die from a storm disturbance
      logit(p[i]) = a0[sp[i], co[i]] + a1[sp[i]]*mix_logit[i] + b[sp[i]]*I[plot[i]]*(dbh[i]^c[sp[i]])
      pd[i] = 1 - (1 - p[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      a1[s] ~ dnorm(0, 0.1)
      b[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c[s] ~ dnorm(0, 1) 
      
      for(country in 1:Ncountry){
       a0[s, country] ~ dnorm(0, 0.1)
      }
    }
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      I[k] ~ dbeta(", beta1, ", ", beta2, ") T(0.001,0.999)
    }
  }"
  )
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model, tmp)
  out <- R2jags::jags.parallel(data = data_jags,
                               param = param,
                               model.file = tmp,
                               n.chains = n.chains,
                               n.iter = n.iter,
                               n.burnin = n.burn,
                               n.thin = n.thin,
                               DIC = TRUE)
  
  # Print the computation time
  stop <- Sys.time()
  print(paste0("Computation time: ", 
               round(as.numeric(difftime(stop, start, units = "mins")), digits = 1), 
               " min."))
  
  return(out)
}




#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_biotic_bin <- function(data_jags, n.chains, n.iter, n.burn, n.thin, 
                                     param = c("a0", "a1", "b", "c", "I")){
  
  # Initialize time
  start <- Sys.time()
  
  ## - For prior of intensity, fit beta distribution on severity
  # Prepare numeric vector containing severity values
  DS_distr <- (data.frame(plot = data_jags$plot, 
                          dead = data_jags$d) %>%
                 group_by(plot) %>%
                 summarise(DS = sum(dead)/n()) %>%
                 mutate(DS = case_when(DS < 0.001 ~ 0.001, 
                                       DS > 0.999 ~ 0.999, 
                                       TRUE ~ DS)))$DS
  # Fit distribution based on this vector, and get parameters
  beta1 <- as.numeric(fitdistr(DS_distr, densfun = "beta", 
                               start = list(shape1 = 1, shape2 = 2))[[1]][1])
  beta2 <- as.numeric(fitdistr(DS_distr, densfun = "beta", 
                               start = list(shape1 = 1, shape2 = 2))[[1]][2])
  
  ## - Write the model
  mortality_model <- paste0(
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pd[i])
      # Probability to die from a storm disturbance
      logit(p[i]) = a0[sp[i], co[i]] + a1[sp[i]]*mixed[i] + b[sp[i]]*I[plot[i]]*(dbh[i]^c[sp[i]])
      pd[i] = 1 - (1 - p[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      a1[s] ~ dnorm(0, 0.1)
      b[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c[s] ~ dnorm(0, 1) 
      
      for(country in 1:Ncountry){
       a0[s, country] ~ dnorm(0, 0.1)
      }
    }
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      I[k] ~ dbeta(", beta1, ", ", beta2, ") T(0.001,0.999)
    }
  }"
  )
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model, tmp)
  out <- R2jags::jags.parallel(data = data_jags,
                               param = param,
                               model.file = tmp,
                               n.chains = n.chains,
                               n.iter = n.iter,
                               n.burnin = n.burn,
                               n.thin = n.thin,
                               DIC = TRUE)
  
  # Print the computation time
  stop <- Sys.time()
  print(paste0("Computation time: ", 
               round(as.numeric(difftime(stop, start, units = "mins")), digits = 1), 
               " min."))
  
  return(out)
}


#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A list containing one rjags object per disturbance
fit_mortality_mix <- function(data_jags, n.chains, n.iter, n.burn, n.thin){
  
  # Identify disturbances
  disturbances.in <- names(data_jags)
  
  # Initialize output
  out <- list()
  
  # Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    print(paste0("Fit mortality model for ", disturbances.in[i], " disturbance"))
    # There's a specific model for storm
    if(disturbances.in[i] %in% c("storm", "snow")){
      jags.i <- fit_mortality_storm(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    } 
    # A mixture effect model for biotic
    if(disturbances.in[i] %in% c("biotic")){
      jags.i <- fit_mortality_biotic(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    } 
    # And one "classical model" for fire and other disturbance
    if(disturbances.in[i] %in% c("fire", "other")){
      jags.i <- fit_mortality_other(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    }
    # Add model to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- jags.i")))
  }
  
  # Return output
  return(out)
}

#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A list containing one rjags object per disturbance
fit_mortality_mix_bin <- function(data_jags, n.chains, n.iter, n.burn, n.thin){
  
  # Identify disturbances
  disturbances.in <- names(data_jags)
  
  # Initialize output
  out <- list()
  
  # Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    print(paste0("Fit mortality model for ", disturbances.in[i], " disturbance"))
    # There's a specific model for storm
    if(disturbances.in[i] %in% c("storm", "snow")){
      jags.i <- fit_mortality_storm(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    } 
    # A mixture effect model for biotic
    if(disturbances.in[i] %in% c("biotic")){
      jags.i <- fit_mortality_biotic_bin(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    } 
    # And one "classical model" for fire and other disturbance
    if(disturbances.in[i] %in% c("fire", "other")){
      jags.i <- fit_mortality_other(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    }
    # Add model to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- jags.i")))
  }
  
  # Return output
  return(out)
}


#' Function to get parameter values per iteration for each species and disturbance
#' @param jags.model.in Jags model to save
#' @param data_jags.in dajagas object to save
#' @param data_model.in Tree data formatted for model (data frame)
get_param_mix <- function(jags.model.in, data_jags.in, data_model.in){
  
  # Initialize correspondence table
  corresp.tables <- list()
  # Initialize weight table
  weight.tables <- list()
  # Initialize scaling table
  scale.tables <- list()
  
  # Loop on all disturbances 
  for(i in 1:length(names(data_jags.in))){
    
    # Table containing the parameters to scale dbh, logratio and mix_logit
    # -- First run some models
    mod.dbh.i <- lm(dbh.scaled ~ dbh, data = data.frame(
      dbh = data_model.in[[i]]$dbh, dbh.scaled = data_jags.in[[i]]$data_jags$dbh))
    mod.logratio.i <- lm(logratio.scaled ~ logratio,  data = data.frame(
      logratio = data_model.in[[i]]$log_dbh_dqm, 
      logratio.scaled = data_jags.in[[i]]$data_jags$logratio))
    mod.mix_logit.i <- lm(mix_logit.scaled ~ mix_logit,  data = data.frame(
      mix_logit = data_model.in[[i]]$mix_logit, 
      mix_logit.scaled = data_jags.in[[i]]$data_jags$mix_logit))
    # -- Store coefficients in a table
    scale.table.i <- list(dbh.intercept = summary(mod.dbh.i)$coefficients[1, 1], 
                          dbh.slope = summary(mod.dbh.i)$coefficients[2, 1], 
                          logratio.intercept = summary(mod.logratio.i)$coefficients[1, 1], 
                          logratio.slope = summary(mod.logratio.i)$coefficients[2, 1], 
                          mix_logit.intercept = summary(mod.mix_logit.i)$coefficients[1, 1], 
                          mix_logit.slope = summary(mod.mix_logit.i)$coefficients[2, 1])
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
  
  # Loop on all disturbances
  for(i in 1:length(names(jags.list))){
    
    # Format the table containing parameter value per species, disturbance and iteration
    param.table.i <- ggs(as.mcmc(jags.list[[i]])) %>%
      # Add disturbance
      mutate(disturbance = names(jags.list)[i]) %>%
      # Extract information on parameter, species and country
      mutate(Param = gsub("\\[.+", "", Parameter), 
             sp = as.integer(ifelse(Param == "a0", gsub(".+\\[", "", gsub("\\,.+", "", Parameter)), 
                                    gsub(".+\\[", "", gsub("\\]", "", Parameter)))), 
             co = as.integer(ifelse(Param == "a0", gsub(".+\\,", "", gsub("\\]", "", Parameter)), NA_integer_))) %>%
      # Remove the estimation of intensity and deviance
      filter(Param != "I") %>%
      filter(Param != "deviance") %>%
      # Add name of the country and species, and weight of each species per country
      left_join(corresp.tables[[i]]$country.table, by = "co") %>%
      left_join(corresp.tables[[i]]$species.table, by = "sp") %>%
      left_join(weight.tables[[i]], by = c("species", "country")) %>%
      # No weight for the parameters that do not rely on the country
      mutate(weight = ifelse(Param == "a0", weight, 1)) %>%
      mutate(weight = ifelse(is.na(weight), 1, weight)) %>%
      # Summarize Parameter value per country (only apply to a0)
      group_by(disturbance, Iteration, Chain, Param, species) %>%
      summarize(val = sum(value*weight, na.rm = TRUE)/sum(weight, na.rm = TRUE)) %>%
      # Format to get one column per parameter
      spread(key = Param, value = val) %>%
      # Set a1 to 0 (dominance effect) if disturbance is not storm or snow
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow", "biotic"), a1, 0)) %>%
      # Add parameters to scale dbh and logratio
      mutate(dbh.intercept = scale.tables[[i]]$dbh.intercept, 
             dbh.slope = scale.tables[[i]]$dbh.slope, 
             logratio.intercept = scale.tables[[i]]$logratio.intercept, 
             logratio.slope = scale.tables[[i]]$logratio.slope, 
             mix_logit.intercept = scale.tables[[i]]$mix_logit.intercept, 
             mix_logit.slope = scale.tables[[i]]$mix_logit.slope)
    
    # Store table in the final table
    if(i == 1) param.table <- param.table.i
    else param.table <- rbind(param.table, param.table.i)
  }
  
  # return the list
  return(param.table)
  
}



#' Figure of parameter value per species for the manuscript
#' @param param_per_iteration value of each parameter at each iteration per species
#' @param disturbance.in Disturbance for which to make the plot
#' @param dbh.ref reference dbh (in mm) to use to calculate sensitivity
#' @param I.ref reference disturbance intensity (from 0 to 1) to use to calculate sensitivity
#' @param time.ref reference time between two measurements (in years) to use to calculate sensitivity
#' @param mix.ref reference proportion of mixture
#' @param file.in Name and path of the plot to save
plot_mix <- function(
  param, disturbance.in = "biotic", dbh.ref = 250, I.ref = 0.75,
  time.ref = 5, mix.ref = 0.3, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify the right color for the disturbance plotted
  color.in <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                 filter(disturbance == disturbance.in))$color
  
  # Format param
  param.in <- param %>%
    filter(disturbance == disturbance.in) %>%
    mutate(iter = paste0("chain", Chain, "_iter", Iteration)) %>%
    ungroup() %>%
    dplyr::select(-disturbance, -Chain, -Iteration)
  
  # Calculate sensitivity, dbh effect and mixture effect per species
  data.in <- expand.grid(species = unique(param.in$species),
                         dbh = dbh.ref, mix = mix.ref, I = I.ref, 
                         iter = unique(param.in$iter)) %>%
    # Add parameters per species
    left_join(param.in, by = c("species", "iter")) %>%
    # Scale variables when needed
    mutate(mix_logit = log(mix/(1 - mix)), 
           dbh.scaled = dbh.intercept + dbh.slope*dbh, 
           mix_logit.scaled = mix_logit.intercept + mix_logit.slope*mix_logit) %>%
    # Compute probabilities
    mutate(pd = 1 - (1 - plogis(a0 + a1*mix_logit.scaled + b*I*dbh.scaled^c))^time.ref) %>%
    # Average per species
    dplyr::select(species, iter, sensitivity = pd, dbh_effect = c, mixture_effect = a1) %>%
    filter(!(species %in% c("Other conifer", "Other broadleaf"))) %>%
    gather(key = "parameter", value = "value", "sensitivity", "dbh_effect", "mixture_effect") %>%
    group_by(species, parameter) %>%
    summarize(mean = mean(value, na.rm = TRUE), 
              lwr = quantile(value, 0.025, na.rm = TRUE), 
              upr = quantile(value, 0.975, na.rm = TRUE))
  
  # Make final plot
  plot.out = data.in %>%
    mutate(parameter = gsub("\\_", "\\ ", parameter), 
           species = factor(species, levels = (
             (data.in %>% 
                filter(parameter == "sensitivity") %>% 
                arrange(mean))$species)), 
           parameter = factor(parameter, levels = c(
             "sensitivity", "mixture effect", "dbh effect")), 
           signif = ifelse(lwr > 0 | upr < 0, "yes", "no")) %>%
    ggplot(aes(x = species, y = mean, ymin = lwr, ymax = upr, color = signif)) +
    geom_errorbar(width = 0) +
    geom_point(fill = color.in, shape = 21) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() + 
    facet_wrap(~ parameter, scales = "free_x") + 
    scale_color_manual(values = c(`yes` = "black", `no` = "gray")) +
    xlab("") + ylab("") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          axis.text.y = element_text(size = 9, face = "italic"), 
          axis.text.x = element_text(size = 7), 
          strip.text = element_text(size = 16), 
          legend.position = "none")
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 25, height = 12, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}

#' Figure of parameter value per species for the manuscript
#' @param param_per_iteration value of each parameter at each iteration per species
#' @param disturbance.in Disturbance for which to make the plot
#' @param dbh.ref reference dbh (in mm) to use to calculate sensitivity
#' @param I.ref reference disturbance intensity (from 0 to 1) to use to calculate sensitivity
#' @param time.ref reference time between two measurements (in years) to use to calculate sensitivity
#' @param mixed.ref Mixed plot (1) or mono specific (0) to calculate sensitivity
#' @param file.in Name and path of the plot to save
plot_mix_bin <- function(
  param, disturbance.in = "biotic", dbh.ref = 250, I.ref = 0.75,
  time.ref = 5, mixed.ref = 0, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify the right color for the disturbance plotted
  color.in <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                 filter(disturbance == disturbance.in))$color
  
  # Format param
  param.in <- param %>%
    filter(disturbance == disturbance.in) %>%
    mutate(iter = paste0("chain", Chain, "_iter", Iteration)) %>%
    ungroup() %>%
    dplyr::select(-disturbance, -Chain, -Iteration)
  
  # Calculate sensitivity, dbh effect and mixture effect per species
  data.in <- expand.grid(species = unique(param.in$species),
                         dbh = dbh.ref, mixed = mixed.ref, I = I.ref, 
                         iter = unique(param.in$iter)) %>%
    # Add parameters per species
    left_join(param.in, by = c("species", "iter")) %>%
    # Scale variables when needed
    mutate(dbh.scaled = dbh.intercept + dbh.slope*dbh) %>%
    # Compute probabilities
    mutate(pd = 1 - (1 - plogis(a0 + a1*mixed + b*I*dbh.scaled^c))^time.ref) %>%
    # Average per species
    dplyr::select(species, iter, sensitivity = pd, dbh_effect = c, mixture_effect = a1) %>%
    filter(!(species %in% c("Other conifer", "Other broadleaf"))) %>%
    gather(key = "parameter", value = "value", "sensitivity", "dbh_effect", "mixture_effect") %>%
    group_by(species, parameter) %>%
    summarize(mean = mean(value, na.rm = TRUE), 
              lwr = quantile(value, 0.025, na.rm = TRUE), 
              upr = quantile(value, 0.975, na.rm = TRUE))
  
  # Make final plot
  plot.out = data.in %>%
    mutate(parameter = gsub("\\_", "\\ ", parameter), 
           species = factor(species, levels = (
             (data.in %>% 
                filter(parameter == "sensitivity") %>% 
                arrange(mean))$species)), 
           parameter = factor(parameter, levels = c(
             "sensitivity", "mixture effect", "dbh effect")), 
           signif = ifelse(lwr > 0 | upr < 0, "yes", "no")) %>%
    ggplot(aes(x = species, y = mean, ymin = lwr, ymax = upr, color = signif)) +
    geom_errorbar(width = 0) +
    geom_point(fill = color.in, shape = 21) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() + 
    facet_wrap(~ parameter, scales = "free_x") + 
    scale_color_manual(values = c(`yes` = "black", `no` = "gray")) +
    xlab("") + ylab("") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          axis.text.y = element_text(size = 9, face = "italic"), 
          axis.text.x = element_text(size = 7), 
          strip.text = element_text(size = 16), 
          legend.position = "none")
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 25, height = 12, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}