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



#' Plot observed vs predicted probabilities
#' @param jags.model rjags object where parameters value are extracted
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM.
#' @param method character indicating if proba ("mean.proba") or parameters 
#'               ("mean.param") are averaged per iteration 
#' @param file.in Name of the file to save
plot_predicted_vs_observed <- function(jags.model, data_jags, data_model, 
                                       method = "mean.proba", file.in){
  
  # Create the directories that are needed
  create_dir_if_needed(file.in)
  
  # Predict death probability for each model
  if(method == "mean.proba") jags.prediction <- predict_jags_meanProba(jags.model, data_jags, data_model)
  if(method == "mean.param") jags.prediction <- predict_jags_meanParam(jags.model, data_jags, data_model)
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                              color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                     filter(disturbance %in% disturbances.in) %>%
                     arrange(disturbance))$color
  
  # Loop on all disturbances to extract data
  for(i in 1:length(disturbances.in)){
    
    # Format data 
    data.i <- data_model[[i]] %>%
      # Add disturbance and whether the tree is dead or not
      mutate(disturbance = disturbances.in[i],
             dead = ifelse(a == 1, 0, 1)) %>%
      dplyr::select(treecode, disturbance, species, dead) %>%
      # Add model prediction
      left_join(jags.prediction[[i]], by = "treecode") 
    
    if(i == 1) data <- data.i
    else data <- rbind.data.frame(data, data.i)
  }
  
  # - Make the plot
  plot.out <- data %>%
    mutate(p = round(p*20, digits = 0)/20) %>%
    group_by(p, disturbance) %>%
    summarize(death.rate = sum(dead)/n(), 
              mean.proba = mean(p),
              n = n()) %>%
    filter(n > 5) %>%
    ggplot(aes(x = death.rate, y = mean.proba, fill = n)) + 
    geom_point(shape = 21, color = "black", alpha = 0.7) + 
    geom_abline(intercept = 0, slope = 1) +
    scale_fill_gradient(low = "blue", high = "red") +
    xlab("Observed death rate") + ylab("Average predicted\ndeath probability") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.key = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 13, height = 11, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
}



#' Plot the rhat of the parameters of a jags model
#' @param jags.model list of rjags object, each element being associated to a disturbance
#' @param file.in Name of the file to save, path included
plot_rhat <- function(jags.model, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Name of the disturbances
  disturbances.in <- names(jags.model)
  
  # Colors for plotting
  colors.in <- (data.frame(disturbance = disturbances.in) %>%
                  left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                       color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                            by = "disturbance"))$color
  
  # Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    
    print(paste0("Calculating rhat for parameters related to ", disturbances.in[i], " disturbance"))
    
    # Convert jags model for distrubance i in a mcmcm object 
    mcmc.i <- ggs(as.mcmc(jags.model[[i]]))
    
    # Initialize dataset for disturbance i
    data.i <- mcmc.i %>%
      mutate(Param.cat = gsub("\\[.+\\]", "", Parameter)) %>%
      filter(Param.cat %in% c("a0", "a1", "b", "c", "I")) %>%
      mutate(Param.cat = factor(Param.cat, levels = c("a0", "a1", "b", "c", "I")), 
             disturbance = disturbances.in[i], 
             rhat = NA_real_) %>%
      dplyr::select(Param.cat, Parameter, disturbance, rhat) %>%
      distinct()
    
    # Loop on all parameters of disturbance i
    for(j in 1:dim(data.i)[1]){
      # Convert mcmc i, parameter j in a matrix with row = iter and col = chain
      matrix.ij <- as.matrix(
        mcmc.i %>%
          filter(Parameter == data.i$Parameter[j]) %>%
          mutate(Chain = paste0("chain", Chain)) %>%
          dplyr::select(Iteration, Chain, value) %>%
          spread(key = "Chain", value = "value") %>%
          dplyr::select(-Iteration))
      # calculate rhat from this matrix
      data.i$rhat[j] <- Rhat(matrix.ij)
    }
    
    # Add to final dataset
    if(i == 1) data.out <- data.i
    if(i > 1) data.out <- rbind(data.out, data.i)
    
  }
  
  # Final plot
  plot.out <- data.out %>%
    mutate(disturbance = factor(disturbance, levels = disturbances.in)) %>%
    ggplot(aes(x = rhat))  + 
    geom_histogram(color = "black", aes(fill = disturbance)) +
    scale_fill_manual(values = colors.in) +
    facet_grid(Param.cat ~ disturbance, scales = "free_y") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.position = "none", 
          strip.text.y = element_text(angle = 360))  + 
    geom_vline(xintercept = 1.1, linetype = "dashed")
  
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 25, height = 15, units = "cm", dpi = 600, bg = "white")
  
  # return the name of all the plots made
  return(file.in)
  
  
}


#' Plot the predictions of the model
#' @param rdata.file rdata file containing correspondance table and model fit
#' @param file.out Name of the file to save, including path
plot_predictions = function(rdata.file, file.out){
  
  # Create directory if needed 
  create_dir_if_needed(file.out)
  
  # Load rdata
  load(rdata.file)
  
  # Loop on all disturbances
  for(i in 1:length(names(jags.list))){
    
    # Format the table containing parameter value per species, disturbance and iteration
    param.table.i <- ggs(as.mcmc(jags.list[[i]])) %>%
      # Add disturbance
      mutate(disturbance = names(jags.list)[i]) %>%
      # Extract information on parameter, species and country
      mutate(Param = gsub("\\[.+", "", Parameter), 
             sp = as.integer(gsub(".+\\[", "", gsub("\\]", "", Parameter)))) %>%
      # Remove the estimation of intensity and deviance
      filter(Param != "I") %>%
      filter(Param != "deviance") %>%
      # Add name of the country and species, and weight of each species per country
      left_join(corresp.tables[[i]]$species.table, by = "sp") %>%
      # Remove useless columns
      dplyr::select(-Parameter, -sp) %>%
      # Format to get one column per parameter
      spread(key = Param, value = value) %>%
      # Set a1 to 0 (dominance effect) if disturbance is not storm or snow
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, 0)) %>%
      # Add parameters to scale dbh and logratio
      mutate(dbh.intercept = scale.tables[[i]]$dbh.intercept, 
             dbh.slope = scale.tables[[i]]$dbh.slope, 
             logratio.intercept = scale.tables[[i]]$logratio.intercept, 
             logratio.slope = scale.tables[[i]]$logratio.slope)
    
    # Store table in the final table
    if(i == 1) param.table <- param.table.i
    else param.table <- rbind(param.table, param.table.i)
  }
  
  
  # Summarize parameter values
  data.param = param.table %>%
    gather(key = "param", value = "value", "a0", "a1", "b", "c", "dbh.intercept", 
           "dbh.slope", "logratio.intercept", "logratio.slope") %>%
    group_by(species, param) %>%
    summarize(value.mean = mean(value)) %>%
    spread(key = "param", value = "value.mean")
  
  
  # Prepare data with predicitons for plotting
  data.fit = expand.grid(species = unique(param.table$species), 
                         dqm = 250, dbh = c(110:800), 
                         I = c(0.1, 0.5, 0.9)) %>%
    left_join(data.param, by = "species") %>%
    mutate(logratio = log(dbh/dqm), 
           logratio.scaled = logratio.intercept + logratio.slope*logratio, 
           dbh.scaled = dbh.intercept + dbh.slope*dbh, 
           p = plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))
  
  # Make the plot
  plot.out = data.fit %>%
    filter(species %in% c("Picea abies", "Abies alba", "Fagus sylvatica")) %>%
    ggplot(aes(x = dbh, y = p, color = I, group = I)) + 
    geom_line() + 
    facet_wrap(~ species) + 
    theme_bw()
  
  
  # - Save the plot
  ggsave(file.out, plot.out, width = 18, height = 8, units = "cm", 
         dpi = 600, bg = "white")
  
  # return the list
  return(file.out)
}

