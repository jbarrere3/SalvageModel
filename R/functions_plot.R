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




#' Plot the rate of dead and harvested trees
#' @param data_model Data used to fit the model
#' @param file.in Plot to save
plot_harvest_and_death_rate <- function(data_model, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify disturbances
  disturbances.in <- names(data_model)
  
  # Loop on all disturbances to paste the different dataset
  for(i in 1:length(disturbances.in)){
    if(i == 1) data <- data_model[[i]]
    else data <- rbind(data, data_model[[i]])
  }
  
  # Make the plot
  plot.out <- data %>%
    mutate(dbh = round(dbh/100, digits = 0)*100) %>%
    group_by(dbh, species) %>%
    summarise(dead = sum(d)/n(), 
              harvested = sum(h)/n(), 
              n = n(), 
              label.pos = (n() - sum(a))/n()) %>%
    filter(n > 10) %>%
    gather(key = "status", value = "rate", "dead", "harvested") %>%
    ggplot(aes(x = dbh, y = rate, fill = status)) + 
    geom_bar(stat="identity", color = "black") + 
    geom_text(nudge_y = 0.1, aes(y = label.pos, label = n), size = 2, inherit.aes = TRUE) + 
    facet_wrap(~ species) + 
    scale_fill_manual(values = c("#CCD5A2", "#D4A373")) +
    xlim(0, 1000) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.key = element_blank(), 
          legend.title = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(file.in)
  
}

#' Plot the relation between disturbance severity (observed) and intensity (estimated)
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param FUNDIV_plot Plot table formatted for FUNDIV
#' @param file.in Path and file where to save the plot
map_disturbance_intensity <- function(jags.model, data_jags, FUNDIV_plot, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Initialize output
  data.plot <- list()
  
  ## - Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    # Create spatial object with data for diturbance i
    data.i <- extract_intensity_per_plotcode(jags.model[[i]], data_jags[[i]]) %>%
      gather(key = "iter", value = "I", colnames(.)[which(colnames(.) != "plotcode")]) %>%
      group_by(plotcode) %>%
      summarize(intensity = mean(I)) %>%
      mutate(disturbance = disturbances.in[i]) %>%
      left_join(FUNDIV_plot, by = "plotcode") %>%
      dplyr::select("plotcode", "latitude", "longitude", 
                    "disturbance", "intensity") %>% 
      spread(key = "disturbance", value = "intensity") %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
    
    # Add to the output list
    eval(parse(text = paste0("data.plot$", disturbances.in[i], " <- data.i")))
  }
  
  ## - First plot for the countries only
  plot.out <- ne_countries(scale = "medium", returnclass = "sf") %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(fill = "#343A40", color = "lightgray", show.legend = F) + 
    annotation_scale(location = "br", width_hint = 0.2) +
    geom_rect(aes(xmin = -9.7, xmax = 9.7, ymin = 36, ymax = 51.5), 
              color = "#BA181B", fill = NA) +
    geom_rect(aes(xmin = 20, xmax = 31.5, ymin = 59.5, ymax = 70.5), 
              color = "#BA181B", fill = NA)   + 
    # Add other disturbance
    geom_sf(data = data.plot$other, shape = 16, aes(color = other), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#90A955", high = "black", midpoint = 0.3)  +
    new_scale_colour() +
    # Add storm disturbance
    geom_sf(data = data.plot$storm, shape = 16, aes(color = storm), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#4361EE", high = "black", midpoint = 0.3)  +
    new_scale_colour() +
    # Add fire disturbance
    geom_sf(data = data.plot$fire, shape = 16, aes(color = fire), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#F77F00", high = "black", midpoint = 0.3) +
    # Finish formatting
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.box = "horizontal",
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"), 
          legend.position = c(0.2, 0.9))
  
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 26, height = 33, units = "cm", dpi = 600)
  return(file.in)
  
}



#' Plot the relation between disturbance severity (observed) and intensity (estimated)
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param FUNDIV_plot Plot table formatted for FUNDIV
#' @param file.in Path and file where to save the plot
map_disturbance_intensity_bis <- function(jags.model, data_jags, FUNDIV_plot, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Initialize output
  data.plot <- list()
  
  ## - Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    # Create spatial object with data for diturbance i
    data.i <- extract_intensity_per_plotcode(jags.model[[i]], data_jags[[i]]) %>%
      gather(key = "iter", value = "I", colnames(.)[which(colnames(.) != "plotcode")]) %>%
      group_by(plotcode) %>%
      summarize(intensity = mean(I)) %>%
      mutate(disturbance = disturbances.in[i]) %>%
      left_join(FUNDIV_plot, by = "plotcode") %>%
      dplyr::select("plotcode", "latitude", "longitude", 
                    "disturbance", "intensity") %>% 
      spread(key = "disturbance", value = "intensity") %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
    
    # Add to the output list
    eval(parse(text = paste0("data.plot$", disturbances.in[i], " <- data.i")))
  }
  
  ## - First plot for the countries only
  plot.out <- ne_countries(scale = "medium", returnclass = "sf") %>%
    mutate(keep = ifelse(sovereignt %in% c("France", "Spain", "Finland"), "yes", "no")) %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(aes(fill = keep), color = "white", show.legend = F) + 
    scale_fill_manual(values = c("#8D99AE", "#343A40")) +
    new_scale_colour() +
    # Add other disturbance
    geom_sf(data = data.plot$other, shape = 16, aes(color = other), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#5F0F40", high = "black", midpoint = 0.3)  +
    new_scale_colour() +
    # Add storm disturbance
    geom_sf(data = data.plot$storm, shape = 16, aes(color = storm), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#4361EE", high = "black", midpoint = 0.3)  +
    new_scale_colour() +
    # Add biotic disturbance
    geom_sf(data = data.plot$biotic, shape = 16, aes(color = biotic), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#90A955", high = "black", midpoint = 0.3)  +
    new_scale_colour() +
    # Add snow disturbance
    geom_sf(data = data.plot$snow, shape = 16, aes(color = snow), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#006D77", high = "black", midpoint = 0.3)  +
    new_scale_colour() +
    # Add fire disturbance
    geom_sf(data = data.plot$fire, shape = 16, aes(color = fire), 
            show.legend = "point", size = 1) +
    scale_color_gradient2(low = "white", mid = "#F77F00", high = "black", midpoint = 0.3) +
    # Finish formatting
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    guides(fill = FALSE) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.box = "horizontal",
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"), 
          legend.position = c(0.3, 0.9))
  
  
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 33.8, height = 42.9, units = "cm", dpi = 600)
  return(file.in)
  
}


#' Plot the relation between disturbance severity (observed) and intensity (estimated)
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param FUNDIV_plot Plot table formatted for FUNDIV
#' @param file.in Path and file where to save the plot
map_disturbance_intensity_ter <- function(jags.model, data_jags, FUNDIV_plot, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Initialize output
  data.plot <- list()
  
  ## - Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    # Create spatial object with data for diturbance i
    data.i <- extract_intensity_per_plotcode(jags.model[[i]], data_jags[[i]]) %>%
      gather(key = "iter", value = "I", colnames(.)[which(colnames(.) != "plotcode")]) %>%
      group_by(plotcode) %>%
      summarize(intensity = mean(I)) %>%
      mutate(disturbance = paste0(toupper(substr(disturbances.in[i], 1, 1)), 
                                  substr(disturbances.in[i], 2, nchar(disturbances.in[i])))) %>%
      mutate(disturbance = factor(disturbance)) %>%
      left_join(FUNDIV_plot, by = "plotcode") %>%
      dplyr::select("plotcode", "latitude", "longitude", "intensity", "disturbance") %>% 
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
    
    # Add to the output dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
  }
  
  ## - Final plot
  plot.out <- ne_countries(scale = "medium", returnclass = "sf") %>%
    mutate(keep = ifelse(sovereignt %in% c("France", "Spain", "Finland"), "yes", "no")) %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(aes(fill = keep), color = "white", show.legend = F) + 
    scale_fill_manual(values = c("#8D99AE", "#343A40")) +
    annotation_scale(location = "br", width_hint = 0.2) + 
    geom_sf(data = data, shape = 16, aes(color = intensity), 
            show.legend = "point", size = 0.1) +
    scale_color_gradient2(low = "#EDF2F4", mid = "#EF233C", high = "#D90429", midpoint = 0.4)  +
    facet_wrap(~ disturbance) +
    guides(fill = FALSE) +
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 15)) + 
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71))
  
  ## - save the plot
  if(length(disturbances.in) == 3) ggsave(file.in, plot.out, width = 25, height = 12.5, units = "cm", dpi = 600)
  if(length(disturbances.in) > 3) ggsave(file.in, (plot.out + theme(legend.position = c(0.9, 0.2),
                                                                    legend.justification = c(0.9, 0.2))),
                                         width = 20, height = 18, units = "cm", dpi = 600)
  return(file.in)
  
}



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 3. Diagnostic plots ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Plot Markov chain convergence of a rjags object
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param dir.in Path where to save the plot
plot_convergence <- function(jags.model, data_jags, dir.in){
  
  # Initialize output
  out <- c()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Loop on all disturbances to extract data for plotting
  for(i in 1:length(disturbances.in)){
    # Extract output for disturbance i
    data.i <- ggs(as.mcmc(jags.model[[i]])) %>%
      # Remove deviance from parameters
      filter(Parameter != "deviance") %>%
      # Separate species code from parameter, add disturbance, average country intercept
      mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter), 
             sp_country = gsub("\\]", "", gsub(".+\\[", "", Parameter)), 
             sp = as.integer(gsub("\\,.+", "", sp_country)), 
             country = ifelse(sp == sp_country, NA_integer_, 
                              as.integer(gsub(".+\\,", "", sp_country))), 
             disturbance = disturbances.in[i]) %>%
      group_by(sp, Parameter.true, Chain, Iteration, disturbance) %>%
      summarise(value = mean(value)) %>%
      ungroup() %>%
      # Remove intensity
      filter(Parameter.true != "I") %>%
      # Add species name
      left_join(data_jags[[i]]$species.table, by = "sp") %>%
      dplyr::select(species, Parameter = Parameter.true, Chain, Iteration, 
                    disturbance, value)
    # Add output to a global mutli disturbance dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
  }
  
  # List of all species studied
  species.in <- unique(arrange(data, species)$species)
  
  # Loop on all species
  for(j in 1:length(species.in)){
    
    # Create directories if needed
    file.in.j <- paste0(dir.in, "/convergence_", gsub(" ", "-", species.in[j]), ".png")
    create_dir_if_needed(file.in.j)
    # Add to the output
    out <- c(out, file.in.j)
    ## - Make the plot
    plot.j <- data %>%
      mutate(Chain = as.factor(Chain), 
             Parameter = as.factor(Parameter), 
             disturbance = as.factor(disturbance)) %>%
      filter(species == species.in[j]) %>%
      ggplot(aes(x = Iteration, y = value, colour = Chain, group = Chain)) + 
      geom_line() + 
      facet_grid(Parameter ~ disturbance, scales = "free", drop = FALSE) + 
      scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
      theme(panel.background = element_rect(fill = "white", color = "black"), 
            panel.grid = element_blank(), 
            strip.background = element_blank(), 
            legend.position = "none", 
            strip.text.y = element_text(angle = 360)) + 
      ggtitle(species.in[j])
    
    ## - Save the plot
    ggsave(file.in.j, plot.j, width = 17, height = 12, units = "cm", dpi = 600)
    
  }
  
  # return the name of all the plots made
  return(out)
}



#' Distribution of disturbance intensity predicted by the model
#' @param jags.model jags model 
#' @param data_jags data used for the jags model
#' @param file.in Name and location of the file to save
plot_intensity_distribution <- function(jags.model, data_jags, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Loop on all disturbances to extract data for plotting
  for(i in 1:length(disturbances.in)){
    # Extract output for disturbance i
    data.i <- ggs(as.mcmc(jags.model[[i]])) %>%
      # Remove deviance from parameters
      filter(Parameter != "deviance") %>%
      # Separate species code from parameter, and add disturbance
      mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter)) %>%
      # Only keep intensity
      filter(Parameter.true == "I") %>%
      mutate(plot = as.integer(gsub("\\]", "", gsub(".+\\[", "", Parameter))), 
             disturbance = disturbances.in[i]) %>%
      # Add species name
      left_join(data_jags[[i]]$plotcode.table, by = "plot") %>%
      group_by(plotcode, disturbance) %>%
      summarize(I = mean(value))
    # Add output to a global mutli disturbance dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
  }
  
  ## - Make the plot
  # Intensity per plot
  plot.out <- data %>%
    ggplot(aes(x = I)) + 
    geom_histogram(color = "lightgray", fill = "black", binwidth = 0.1) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank()) + 
    facet_wrap(~ disturbance, nrow = 1, scales = "free_y") + 
    xlab("Estimated disturbance intensity") 
  
  
  # - Save the three plots
  ggsave(file.in, plot.out, width = 20, height = 7, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(file.in)
}


#' Plot death probability predicted by the reference model
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param file.in Name of the file to save
plot_prediction_all <- function(jags.model, data_jags, data_model, file.in){
  
  # Create the diretories that are needed
  create_dir_if_needed(file.in)
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                              color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                     filter(disturbance %in% disturbances.in) %>%
                     arrange(disturbance))$color
  
  # Loop on all disturbances to extract data for plotting
  for(i in 1:length(disturbances.in)){
    
    # Model to scale logratio
    scale_logratio <- lm(logratio.scaled ~ logratio, 
                         data = data.frame(logratio = data_model[[i]]$log_dbh_dqm, 
                                           logratio.scaled = data_jags[[i]]$data_jags$logratio))
    
    # Format data for plotting
    data.i <- data_model[[i]] %>%
      mutate(disturbance = disturbances.in[i], 
             ratio = exp(log_dbh_dqm),
             ratio = round(ratio*2, digits = 0)/2,
             logratio = log(ratio)) %>%
      dplyr::select(species, dbh, ratio, logratio, disturbance) %>%
      mutate(dbh.scaled = data_jags[[i]]$data_jags$dbh, 
             logratio.scaled = predict(scale_logratio, newdata = .), 
             I = 0.5) %>%
      filter(ratio %in% c(0.5, 1, 2)) %>%
      # Add parameters per species
      left_join((ggs(as.mcmc(jags.model[[i]])) %>%
                   # Remove deviance from parameters
                   filter(Parameter != "deviance") %>%
                   # Separate species code from parameter, and add disturbance
                   mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter), 
                          sp_country = gsub("\\]", "", gsub(".+\\[", "", Parameter)), 
                          sp = as.integer(gsub("\\,.+", "", sp_country)), 
                          country = ifelse(sp == sp_country, NA_integer_, 
                                           as.integer(gsub(".+\\,", "", sp_country))), 
                          disturbance = disturbances.in[i]) %>%
                   group_by(sp, Parameter.true, Chain, Iteration, disturbance) %>%
                   summarise(value = mean(value)) %>%
                   ungroup() %>%
                   # Remove intensity
                   filter(Parameter.true != "I") %>%
                   # Add species name
                   left_join(data_jags[[i]]$species.table, by = "sp") %>%
                   group_by(species, Parameter.true) %>%
                   summarise(value = mean(value)) %>%
                   spread(key = "Parameter.true", value = "value")), 
                by = "species") %>%
      mutate(ratio = paste0("dbh/dqm = ", ratio)) %>%
      filter(dbh < 750)
    
    # Calculate death probability, depending on the model
    if(disturbances.in[i] == "storm") data.i <- mutate(data.i, pd = plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))
    if(disturbances.in[i] != "storm") data.i <- mutate(data.i, pd = plogis(a0 + b*I*dbh.scaled^c))
    
    # Only keep columns of interest
    data.i <- data.i %>% 
      dplyr::select(dbh, pd, ratio, disturbance, species)
    
    # Add output to a global mutli disturbance dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
  }
  
  # - Pre-format the data before plotting
  plot.out <- data %>%
    ggplot(aes(x = dbh, y = pd, group = interaction(ratio, disturbance), 
               color = disturbance, linetype = ratio)) + 
    geom_line() + 
    facet_wrap(~ species) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.key = element_blank(),
          legend.title = element_blank()) + 
    scale_color_manual(values = color.vector) +
    ylab("Probability to die from a disturbance") + 
    xlab("dbh (mm)")
  
  # - Save the three plots
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(file.in)
}



#' Plot death probability predicted by the reference model
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param dir.in Directory where to save the plots
plot_prediction <- function(jags.model, data_jags, data_model, dir.in){
  
  # Initialize output
  out <- c()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                              color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                     filter(disturbance %in% disturbances.in) %>%
                     arrange(disturbance))$color
  
  # Loop on all disturbances to extract data for plotting
  for(i in 1:length(disturbances.in)){
    
    # Create the diretories that are needed
    file.in.i <- paste0(dir.in, "/", disturbances.in[i], ".png")
    create_dir_if_needed(file.in.i)
    
    # Add file to the output
    out <- c(out, file.in.i)
    
    # Model to scale logratio
    scale_logratio <- lm(logratio.scaled ~ logratio, 
                         data = data.frame(logratio = data_model[[i]]$log_dbh_dqm, 
                                           logratio.scaled = data_jags[[i]]$data_jags$logratio))
    
    # Format data for plotting
    data.i <- data_model[[i]] %>%
      mutate(disturbance = disturbances.in[i], 
             ratio = exp(log_dbh_dqm),
             ratio = round(ratio*2, digits = 0)/2,
             logratio = log(ratio)) %>%
      dplyr::select(species, dbh, ratio, logratio, disturbance) %>%
      mutate(dbh.scaled = data_jags[[i]]$data_jags$dbh, 
             logratio.scaled = predict(scale_logratio, newdata = .), 
             I_02 = 0.2, I_05 = 0.5, I_08 = 0.8) %>%
      gather(key = "I_code", value = "I", "I_02", "I_05", "I_08") %>%
      filter(ratio %in% c(0.5, 1, 2)) %>%
      # Add parameters per species
      left_join((ggs(as.mcmc(jags.model[[i]])) %>%
                   # Remove deviance from parameters
                   filter(Parameter != "deviance") %>%
                   # Separate species code from parameter, and add disturbance
                   mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter), 
                          sp_country = gsub("\\]", "", gsub(".+\\[", "", Parameter)), 
                          sp = as.integer(gsub("\\,.+", "", sp_country)), 
                          country = ifelse(sp == sp_country, NA_integer_, 
                                           as.integer(gsub(".+\\,", "", sp_country))), 
                          disturbance = disturbances.in[i]) %>%
                   group_by(sp, Parameter.true, Chain, Iteration, disturbance) %>%
                   summarise(value = mean(value)) %>%
                   ungroup() %>%
                   # Remove intensity
                   filter(Parameter.true != "I") %>%
                   # Add species name
                   left_join(data_jags[[i]]$species.table, by = "sp") %>%
                   group_by(species, Parameter.true) %>%
                   summarise(value = mean(value)) %>%
                   spread(key = "Parameter.true", value = "value")), 
                by = "species") %>%
      mutate(ratio = paste0("dbh/dqm = ", ratio)) %>%
      filter(dbh < 750) %>%
      filter(dbh > 150)
    
    # Calculate death probability, depending on the model
    if(disturbances.in[i] == "storm") data.i <- mutate(data.i, pd = plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))
    if(disturbances.in[i] != "storm") data.i <- mutate(data.i, pd = plogis(a0 + b*I*dbh.scaled^c))
    
    # - Pre-format the data before plotting
    plot.out.i <- data.i %>% 
      dplyr::select(dbh, pd, ratio, species, I) %>%
      mutate(intensity = paste0("I = ", as.character(I))) %>%
      ggplot(aes(x = dbh, y = pd, group = interaction(ratio, intensity), 
                 color = intensity, linetype = ratio)) + 
      geom_line() + 
      facet_wrap(~ species) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            strip.background = element_blank(), 
            panel.grid = element_blank(), 
            legend.key = element_blank(),
            legend.title = element_blank()) + 
      scale_color_manual(values = color.vector) +
      ylab(paste0("Probability to die from ", disturbances.in[i], " disturbance")) + 
      xlab("dbh (mm)")
    
    # - Save the three plots
    ggsave(file.in.i, plot.out.i, width = 25, height = 20, units = "cm", dpi = 600)
  }
  
  # return the name of all the plots made
  return(out)
}

#' Plot death probability predicted by the reference model
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM. 
#' @param file.in Name of the file to save
plot_prediction2 <- function(jags.model, data_jags, data_model, file.in){
  
  # Create the diretories that are needed
  create_dir_if_needed(file.in)
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                              color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                     filter(disturbance %in% disturbances.in) %>%
                     arrange(disturbance))$color
  
  # Loop on all disturbances to extract data for plotting
  for(i in 1:length(disturbances.in)){
    
    # Parameter per species and per iteration
    param.per.species.i <- ggs(as.mcmc(jags.model[[i]])) %>%
      # Remove deviance from parameters
      filter(Parameter != "deviance") %>%
      # Separate species code from parameter, and add disturbance
      mutate(Parameter.true = gsub("\\[.+\\]", "", Parameter), 
             sp_country = gsub("\\]", "", gsub(".+\\[", "", Parameter)), 
             sp = as.integer(gsub("\\,.+", "", sp_country)), 
             country = ifelse(sp == sp_country, NA_integer_, 
                              as.integer(gsub(".+\\,", "", sp_country))), 
             disturbance = disturbances.in[i]) %>%
      group_by(sp, Parameter.true, Chain, Iteration, disturbance) %>%
      summarise(value = mean(value)) %>%
      ungroup() %>%
      # Remove intensity
      filter(Parameter.true != "I") %>%
      # Add species name
      left_join(data_jags[[i]]$species.table, by = "sp") %>%
      # Create new column combining chain and iteration
      mutate(iter = paste0("chain", Chain, "_iter", Iteration)) %>%
      # Compute mean parameter value 
      dplyr::select(species, Parameter.true, iter, disturbance, value) %>%
      spread(key = "Parameter.true", value = "value") %>%
      mutate(a1 = ifelse(disturbance == "storm", a1, NA_real_))
    
    
    # Model to scale logratio
    scale_logratio <- lm(logratio.scaled ~ logratio, 
                         data = data.frame(logratio = data_model[[i]]$log_dbh_dqm, 
                                           logratio.scaled = data_jags[[i]]$data_jags$logratio))
    
    # Model to scale dbh
    scale_dbh <- lm(dbh.scaled ~ dbh, 
                    data = data.frame(dbh = data_model[[i]]$dbh, 
                                      dbh.scaled = data_jags[[i]]$data_jags$dbh))
    
    # Format data for plotting
    data.i <- expand.grid(species = unique(data_model[[i]]$species), 
                          dbh = c(100:750), 
                          logratio = 0, 
                          I = 0.5, 
                          iter = unique(param.per.species.i$iter)) %>%
      # Remove dbh out of the range
      left_join((data_model[[i]] %>%
                   mutate(dbh = round(dbh, digits = 0)) %>%
                   group_by(species) %>%
                   summarize(dbh.min = min(dbh), 
                             dbh.max = max(dbh))),
                by = "species") %>%
      filter(dbh > dbh.min & dbh < dbh.max) %>%
      # scale dbh and logratio
      mutate(dbh.scaled = predict(scale_dbh, newdata = .), 
             logratio.scaled = predict(scale_logratio, newdata = .)) %>%
      # Add parameter per species
      left_join(param.per.species.i, by = c("species", "iter")) %>%
      # Compute probabilities
      mutate(pd = ifelse(disturbance == "storm", 
                         plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c), 
                         plogis(a0 + b*I*dbh.scaled^c))) %>%
      # Average per species and dbh
      group_by(species, dbh, disturbance) %>%
      mutate(p = mean(pd), 
             p_975 = as.numeric(quantile(pd, probs = 0.975)), 
             p_025 = as.numeric(quantile(pd, probs = 0.025)))
    
    
    # Add output to a global mutli disturbance dataset
    if(i == 1) data <- data.i
    else data <- rbind(data, data.i)
  }
  
  # - Pre-format the data before plotting
  plot.out <- data %>%
    ggplot(aes(x = dbh, y = p, group = disturbance, fill = disturbance,
               color = disturbance)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = p_025, ymax = p_975), alpha = 0.5) +
    facet_wrap(~ species) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.key = element_blank(),
          legend.title = element_blank()) + 
    scale_color_manual(values = color.vector) +
    scale_fill_manual(values = color.vector) +
    ylab("Probability to die from a disturbance") + 
    xlab("dbh (mm)")
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(file.in)
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
    mutate(p = round(p*10, digits = 0)/10) %>%
    group_by(p, disturbance, species) %>%
    summarize(death.rate = sum(dead)/n(), 
              n = n()) %>%
    filter(n > 10) %>%
    ggplot(aes(x = death.rate, y = p, fill = disturbance)) + 
    geom_point(shape = 21, color = "black", alpha = 0.7) + 
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(~ species) +
    scale_fill_manual(values = color.vector) +
    xlab("Observed death rate") + ylab("Predicted death probability") + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          strip.background = element_blank(), 
          panel.grid = element_blank(), 
          legend.key = element_blank())
  
  # - Save the plot
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(file.in)
}




#' Get disturbance sensitivity by species
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM.
#' @param dbh.ref dbh value for which to predict values
#' @param logratio.ref value of the ratio dbh/dqm for which to predict values
#' @param I.ref disturbance intensity for which to predict values
#' @param file.in Name and path of the plot to save
plot_param_per_species <- function(jags.model, data_jags, data_model, 
                                   dbh.ref = 300, logratio.ref = 0, I.ref = 0.5, 
                                   file.in){
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = disturbances.in) %>%
                     left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                               by = "disturbance"))$color
  
  # Initialize number of species per disturbance
  n.sp.per.dist <- c()
  
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
      summarise(sensitivity_mean = mean(pd), 
                sensitivity_ql = as.numeric(quantile(pd, probs = 0.025)), 
                sensitivity_qh = as.numeric(quantile(pd, probs = 0.975)), 
                dbh.effect_mean = mean(c), 
                dbh.effect_ql = as.numeric(quantile(c, probs = 0.025)), 
                dbh.effect_qh = as.numeric(quantile(c, probs = 0.975)), 
                dominance.effect_mean = mean(a1, na.rm = TRUE), 
                dominance.effect_ql = as.numeric(quantile(a1, probs = 0.025, na.rm = TRUE)), 
                dominance.effect_qh = as.numeric(quantile(a1, probs = 0.975, na.rm = TRUE))) %>%
      # Arrange species by sensitivity value
      ungroup() %>%
      filter(!(species %in% c("Other conifer", "Other broadleaf"))) %>%
      mutate(species = factor(species, levels = .$species[order(.$sensitivity_mean)])) %>%
      # Format for plotting
      gather(key = "variable", value = "value", colnames(.)[which(colnames(.) != "species")]) %>%
      separate(col = "variable", into = c("parameter", "value.type"), sep = "_") %>%
      spread(key = "value.type", value = "value") %>%
      mutate(parameter = gsub("\\.", "\\ ", parameter)) %>%
      mutate(parameter = ifelse(parameter != "sensitivity", parameter, 
                                paste0(disturbances.in[i], " sensitivity"))) %>%
      mutate(parameter = factor(parameter, 
                                levels = c(paste0(disturbances.in[i], " sensitivity"), 
                                           "dbh effect", "dominance effect")))
    
    # Make the plot
    plot.i <- data.i %>%
      ggplot(aes(x = species, y = mean, ymin = ql, ymax = qh)) +
      geom_errorbar(width = 0, color = color.vector[i]) +
      geom_point(color = color.vector[i]) + 
      coord_flip() + 
      facet_wrap(~ parameter, scales = "free_x") + 
      geom_hline(yintercept = 0, linetype = "dashed") + 
      xlab("") + ylab("") +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            strip.background = element_blank())
    
    # Add to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- plot.i")))
    
    # Add also the number of species
    n.sp.per.dist <- c(n.sp.per.dist, length(unique(data.i$species)))
  }
  
  
  # Assemble all plots
  plot.out <- plot_grid(plotlist = out, ncol = 1, align = "v", rel_heights = (n.sp.per.dist+10),
                        labels = paste0("(", letters[c(1:length(disturbances.in))], ")"))
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 18, height = 30, units = "cm", dpi = 600)
  return(file.in)
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Traits analysis ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Make a plot showing disturbance sensitivity against traits with an RDA
#' @param traits dataset containing trait values per species
#' @param disturbance_sensitivity dataset containing disturbance sensitivity per species
#' @param file.in File to export
plot_rda_traits_vs_sensitivity <- function(traits, disturbance_sensitivity, 
                                           file.in){
  
  ##%%%%%%%%%
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  
  ##%%%%%%%%%
  ## - Build plot
  
  # Remove na in trait table
  traits <- traits %>%
    # Add sensitivity
    left_join(disturbance_sensitivity, by = "species") %>%
    drop_na()
  
  # scale trait dataset
  traits.scaled <- scale(as.matrix(traits %>% dplyr::select(colnames(.)[!(colnames(.) %in% c("species", "p", "p_975", "p_025"))])))
  
  # Make pca
  rda <- rda(traits.scaled ~ p, traits,  scale = FALSE)
  
  
  # Results of the rda
  res.rda.ind <- data.frame(species = traits$species, 
                            rda1 = as.numeric(summary(rda)$sites[, 1])) %>%
    mutate(species.short = paste0(substr(species, 1, 2), substr(gsub(".+\\ ", "", species), 1, 2)))
  
  # Identify maximum coordinates of species in rda
  max.rda1 <- max(sqrt(res.rda.ind$rda1^2)) + 0.1
  
  # Data for model
  data.model <- disturbance_sensitivity %>%
    # Logit transformation
    mutate(sensitivity.logit = log(p/(1 - p)), 
           w = 1/(p_975 - p_025)) %>%
    # Add pca
    left_join(res.rda.ind, by = "species") %>%
    drop_na()
  
  
  # Make regression sensitivity vs rda
  model.rda1 <- lm(sensitivity.logit ~ rda1, weights = w, data = data.model)
  
  # Data for the predictions
  data.fit <- data.frame(rda1 = c(round(-max.rda1*1000, digits = 0):round(max.rda1*1000, digits = 0))/1000) %>%
    mutate(fit.logit = predict(model.rda1, newdata = .), 
           fit.se = predict(model.rda1, newdata = ., se.fit = TRUE)$se.fit,
           fit = plogis(fit.logit), 
           fit.inf = plogis(fit.logit - fit.se), 
           fit.sup = plogis(fit.logit + fit.se), 
           p = NA_real_)
  
  # Plot predictions pca1
  plot.model <- data.model %>%
    ggplot(aes(x = rda1, y = p, group = 1)) + 
    geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0) +
    geom_point(size = 2, shape = 21, fill = "lightgray", color = "black") + 
    geom_line(data = data.fit, aes(y = fit), inherit.aes = TRUE) + 
    geom_ribbon(data = data.fit, aes(ymin = fit.inf, ymax = fit.sup), 
                alpha = 0.5, fill = "gray", inherit.aes = TRUE) +
    xlab(paste0("RDA1 (", round(summary(rda)$cont$importance[2, 1]*100, digits = 2), "%)")) + 
    ylab(paste0("Disturbance sensitivity")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank()) + 
    ggtitle(paste0("F = ", round(anova(model.rda1)[1, 4], digits = 1), ", ",
                   scales::pvalue(anova(model.rda1)[1, 5], add_p = TRUE, accuracy = 0.01))) + 
    xlim(-max.rda1, max.rda1)
  
  # Plot rda
  plot.rda <- data.frame(trait = as.character(rownames(summary(rda)$species)), 
                         rda1 = summary(rda)$species[, 1], 
                         position = c(1:(dim(summary(rda)$species)[1]))) %>%
    ggplot(aes(x = 0, xend = rda1, y = position, yend = position)) + 
    geom_segment(arrow = arrow(length = unit(0.2, "cm"))) + 
    scale_y_continuous(breaks = c(1:(dim(summary(rda)$species)[1])), 
                       labels = as.character(rownames(summary(rda)$species)), 
                       expand = c(0.25, 0)) +
    geom_vline(xintercept = 0, color = "darkgray") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank()) + 
    ylab("") + xlim(-max.rda1, max.rda1) + 
    xlab(paste0("RDA1 (", round(summary(rda)$cont$importance[2, 1]*100, digits = 2), "%)"))
  
  
  # Plot all graphs together
  library(cowplot)
  plot.out <- plot_grid(plot.rda, plot.model, nrow = 2, rel_heights = c(0.5, 1), 
                        labels = c("(a)", "(b)"), align = "v")
  
  ##%%%%%%%%%
  ## - Export plot
  ggsave(file.in, plot.out, width = 14, height = 11, units = "cm", dpi = 600)
  return(file.in)
  
}


#' Make a plot showing sensitivity to each disturbance against traits with an RDA
#' @param traits dataset containing trait values per species
#' @param disturbance_sensitivity list of dataset containing disturbance sensitivity per species
#' @param file.in Where to save the plot
plot_rda_traits_vs_sensitivity_allDist <- function(traits, disturbance_sensitivity, file.in){
  
  ##%%%%%%%%%
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  ## - Initialize output
  plots.out <- list()
  
  ## - Loop on all type of disturbances
  for(i in 1:length(disturbances.in)){
    
    # Remove na in trait table
    traits.i <- traits %>%
      # Add sensitivity
      left_join(disturbance_sensitivity[[i]], by = "species") %>%
      drop_na()
    
    # scale trait dataset
    traits.scaled.i <- scale(as.matrix(traits.i %>% dplyr::select(colnames(.)[!(colnames(.) %in% c("species", "p", "p_975", "p_025"))])))
    
    # Make pca
    rda <- rda(traits.scaled.i ~ p, traits.i,  scale = FALSE)
    
    
    # Results of the rda
    res.rda.ind <- data.frame(species = traits.i$species, 
                              rda1 = as.numeric(summary(rda)$sites[, 1])) %>%
      mutate(species.short = paste0(substr(species, 1, 2), substr(gsub(".+\\ ", "", species), 1, 2)))
    
    # Identify maximum coordinates of species in rda
    max.rda1 <- max(sqrt(res.rda.ind$rda1^2)) + 0.1
    
    # Data for model
    data.model <- disturbance_sensitivity[[i]] %>%
      # Logit transformation
      mutate(sensitivity.logit = log(p/(1 - p)), 
             w = 1/(p_975 - p_025)) %>%
      # Add pca
      left_join(res.rda.ind, by = "species") %>%
      drop_na()
    
    
    # Make regression sensitivity vs rda
    model.rda1 <- lm(sensitivity.logit ~ rda1, weights = w, data = data.model)
    
    # Data for the predictions
    data.fit <- data.frame(rda1 = c(round(-max.rda1*1000, digits = 0):round(max.rda1*1000, digits = 0))/1000) %>%
      mutate(fit.logit = predict(model.rda1, newdata = .), 
             fit.se = predict(model.rda1, newdata = ., se.fit = TRUE)$se.fit,
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.logit - fit.se), 
             fit.sup = plogis(fit.logit + fit.se), 
             p = NA_real_)
    
    # Plot predictions pca1
    plot.model <- data.model %>%
      ggplot(aes(x = rda1, y = p, group = 1)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = color.in[i], color = "#343A40") + 
      geom_line(data = data.fit, aes(y = fit), inherit.aes = TRUE, color = color.in[i]) + 
      geom_ribbon(data = data.fit, aes(ymin = fit.inf, ymax = fit.sup), 
                  alpha = 0.5, fill = color.in[i], inherit.aes = TRUE) +
      xlab(paste0("RDA1 (", round(summary(rda)$cont$importance[2, 1]*100, digits = 2), "%)")) + 
      ylab("Disturbance sensitivity") +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank()) + 
      xlim(-max.rda1, max.rda1) + 
      scale_y_continuous(breaks = c(0:5)*0.2)
    
    # Plot rda
    plot.rda <- data.frame(trait = as.character(rownames(summary(rda)$species)), 
                           rda1 = summary(rda)$species[, 1], 
                           position = c(1:(dim(summary(rda)$species)[1]))) %>%
      ggplot(aes(x = 0, xend = rda1, y = position, yend = position)) + 
      geom_segment(arrow = arrow(length = unit(0.1, "cm")), type = "closed") + 
      scale_y_continuous(breaks = c(1:(dim(summary(rda)$species)[1])), 
                         labels = as.character(rownames(summary(rda)$species)), 
                         expand = c(0.25, 0)) +
      geom_vline(xintercept = 0, color = "darkgray") +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank()) + 
      ylab("") + xlim(-max.rda1, max.rda1) + xlab("") +
      labs(title = paste0(toupper(substr(disturbances.in[i], 1, 1)), 
                          substr(disturbances.in[i], 2, nchar(disturbances.in[i]))), 
           subtitle = paste0("F = ", round(anova(model.rda1)[1, 4], digits = 1), ", ",
                             scales::pvalue(anova(model.rda1)[1, 5], add_p = TRUE, accuracy = 0.01)))
    
    # Adjust axis and titles depending on graph position
    if(!((i == 1) | (i == 4 & length(disturbances.in) > 3))){
      plot.rda <- plot.rda + theme(axis.text.y = element_blank(), 
                                   axis.ticks.y = element_blank())
      plot.model <- plot.model + ylab("") 
    } 
    
    
    # Plot all graphs together
    plot.i <- plot_grid(plot.rda, plot.model, nrow = 2, rel_heights = c(0.6, 1), align = "v")
    
    # Add to the output list
    eval(parse(text = paste0("plots.out$", disturbances.in[i], " <- plot.i")))
  }
  
  # Adjust final plot dimensions depending on the number of disturbances
  if(length(disturbances.in) > 3){
    # Draw final plot
    plot.out <- plot_grid(plotlist = plots.out, nrow = 2, align = "hv", 
                          labels = paste0("(", letters[c(1:length(disturbances.in))], ")"))
    
    # Export plot
    ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600)
  }
  
  if(length(disturbances.in) <= 3){
    # Draw final plot
    plot.out <- plot_grid(plotlist = plots.out, nrow = 1, align = "hv", 
                          labels = paste0("(", letters[c(1:length(disturbances.in))], ")"))
    
    # Export plot
    ggsave(file.in, plot.out, width = 30, height = 10, units = "cm", dpi = 600)
  }
  
  return(file.in)
  
}



#' Make one plot per trait of disturbance sensitivity against traits
#' @param traits dataset containing trait values per species
#' @param disturbance_sensitivity dataset containing disturbance sensitivity per species
#' @param disturbance.in Character indicating the type of disturbance
#' @param directory where to export the plots
plot_traits_vs_sensitivity <- function(traits, disturbance_sensitivity, disturbance.in, dir.in){
  
  ## - Initialize output
  out <- c()
  
  # Use the right disturbance sensitivity dataset
  eval(parse(text = paste0("disturbance_sensitivity.in <- disturbance_sensitivity$", disturbance.in)))
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                              color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) %>%
                     filter(disturbance == disturbance.in) %>%
                     arrange(disturbance))$color
  
  # Loop on all traits
  for(trait.i in colnames(traits)[which(colnames(traits) != "species")]){
    
    # Create directory if needed
    file.in.i <- paste0(dir.in, "/", gsub("\\.", "", trait.i), ".png")
    create_dir_if_needed(file.in.i)
    
    # Format data for the model
    data.i <- disturbance_sensitivity.in %>%
      left_join((traits %>% dplyr::select("species", "trait" = trait.i)), 
                by = "species") %>%
      # Logit transformation and weight
      mutate(sensitivity.logit = log(p/(1 - p)), 
             w = 1/(p_975 - p_025)) %>%
      # remove NA
      drop_na()
    
    # Fit a model
    model.i <- lm(sensitivity.logit ~ trait, weights = w, data = data.i)
    
    # Name of the trait formatted for plotting
    trait.label = ifelse(
      substr(trait.i, 1, 3) != "TRY", trait.i,
      paste0(gsub("\\_.+", "", gsub("TRY\\_", "", trait.i)), " (", gsub(".+\\_", "", gsub("TRY\\_", "", trait.i)), ")"))
    
    # Data with the predictions of the model
    data.fit <- data.frame(
      trait = c(round(min(data.i$trait)*1000, digits = 0):round(max(data.i$trait)*1000, digits = 0))/1000) %>%
      mutate(fit.logit = predict(model.i, newdata = .), 
             fit.se = predict(model.i, newdata = ., se.fit = TRUE)$se.fit,
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.logit - fit.se), 
             fit.sup = plogis(fit.logit + fit.se), 
             p = NA_real_)
    
    # Plot predictions 
    plot.i <- data.i %>%
      ggplot(aes(x = trait, y = p)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = color.in, color = "#343A40") + 
      geom_line(data = data.fit, aes(y = fit, group = 1), color = color.in, inherit.aes = TRUE) + 
      geom_ribbon(data = data.fit, aes(ymin = fit.inf, ymax = fit.sup), 
                  alpha = 0.5, fill = color.in, inherit.aes = TRUE) +
      ylab(paste0("Sensitivity to ", disturbance.in, " disturbance")) + xlab(trait.label) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank()) + 
      ggtitle(paste0("F = ", round(anova(model.i)[1, 4], digits = 1), ", ",
                     scales::pvalue(anova(model.i)[1, 5], add_p = TRUE, accuracy = 0.01))) 
    
    # Save the plot
    ggsave(file.in.i, plot.i, width = 14, height = 11, units = "cm", dpi = 600)
    
    # Add filename to the output
    out <- c(out, file.in.i)
    
  }
  
  return(out)
  
}

#' Plot traits vs disturbance sensitivity for all disturbances
#' @param traits dataset containing trait values per species
#' @param disturbance_sensitivity dataset containing disturbance sensitivity per species
#' @param dir.in where to export the plots
plot_traits_vs_sensitivity_allDist <- function(traits, disturbance_sensitivity, dir.in){
  
  # Initialize output
  out <- c()
  
  # All disturbances 
  disturbances.in <- names(disturbance_sensitivity)
  
  # Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    # Apply function for disturbance i
    files.i <- plot_traits_vs_sensitivity(traits, disturbance_sensitivity, disturbances.in[i], 
                                          paste(dir.in, disturbances.in[i], sep = "/"))
    # Add file to the final output
    out <- c(out, files.i)
  }
  
  # Return output
  return(files.i)
}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Outdated functions ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#' Plot the comparison of disturbance intensity by two different models
#' @param jags.model.1 jags model containing intensity value for model 1
#' @param data_jags.1 data used for the model 1
#' @param jags.model.2 jags model containing intensity value for model 2
#' @param data_jags.2 data used for the model 2
#' @param file.in Name and location of the file to save
plot_intensity1_vs_intensity2 <- function(jags.model.1, data_jags.1, jags.model.2, 
                                          data_jags.2, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  
  
  ## - Build the plot
  # Intensity for the first model
  data.intensity.1 <- data.frame(plot = data_jags.1$data_jags$plot, 
                                 Ifire = data_jags.1$data_jags$Dfire, 
                                 Istorm = data_jags.1$data_jags$Dstorm, 
                                 Iother = data_jags.1$data_jags$Dother) %>%
    distinct() %>%
    gather(key = "variable", value = "value", "Ifire", "Istorm", "Iother") %>%
    mutate(Parameter = paste0(variable, "[", plot, "]")) %>%
    filter(value == 1) %>%
    left_join((ggs(as.mcmc(jags.model.1)) %>%
                 filter(substr(Parameter, 1, 1) == "I") %>%
                 group_by(Parameter) %>%
                 dplyr::summarize(intensity.mod1 = mean(value, na.rm = T))), 
              by = "Parameter") %>%
    mutate(disturbance = gsub("I", "", variable)) %>%
    left_join(data_jags.1$plotcode_table, by = "plot") %>%
    dplyr::select(plotcode, disturbance, intensity.mod1)
  # Intensity for the second model
  data.intensity.2 <- data.frame(plot = data_jags.2$data_jags$plot, 
                                 Ifire = data_jags.2$data_jags$Dfire, 
                                 Istorm = data_jags.2$data_jags$Dstorm, 
                                 Iother = data_jags.2$data_jags$Dother) %>%
    distinct() %>%
    gather(key = "variable", value = "value", "Ifire", "Istorm", "Iother") %>%
    mutate(Parameter = paste0(variable, "[", plot, "]")) %>%
    filter(value == 1) %>%
    left_join((ggs(as.mcmc(jags.model.2)) %>%
                 filter(substr(Parameter, 1, 1) == "I") %>%
                 group_by(Parameter) %>%
                 dplyr::summarize(intensity.mod2 = mean(value, na.rm = T))), 
              by = "Parameter") %>%
    mutate(disturbance = gsub("I", "", variable)) %>%
    left_join(data_jags.2$plotcode_table, by = "plot") %>%
    dplyr::select(plotcode, intensity.mod2)
  # Plot the two intensities
  plot.out <- data.intensity.1 %>%
    left_join(data.intensity.2, by = "plotcode") %>%
    ggplot(aes(x = intensity.mod1, y = intensity.mod2, fill = disturbance)) + 
    geom_point(shape = 21, color = "black") + 
    geom_abline(intercept = 0, slope = 1) + 
    facet_wrap(~ disturbance) + 
    scale_fill_manual(values = c("#F77F00", "#90A955", "#4361EE")) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          legend.position = "none")
  
  
  # - Save the three plots
  ggsave(file.in, plot.out, width = 16, height = 6, units = "cm", dpi = 600)
  
  # return the name of all the plots made
  return(file.in)
}
