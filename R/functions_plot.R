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
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600, bg = "white")
  
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
  ggsave(file.in, plot.out, width = 26, height = 33, units = "cm", dpi = 600, bg = "white")
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
  ggsave(file.in, plot.out, width = 33.8, height = 42.9, units = "cm", dpi = 600, bg = "white")
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
  if(length(disturbances.in) <= 3) ggsave(file.in, plot.out, width = 25, height = 12.5, units = "cm", dpi = 600, bg = "white")
  if(length(disturbances.in) > 3) ggsave(file.in, (plot.out + theme(legend.position = c(0.9, 0.2),
                                                                    legend.justification = c(0.9, 0.2))),
                                         width = 20, height = 18, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}



#' Function to plot trends in disturbance occurence and severity over time
#' @param FUNDIV_plot plot table formatted for FUNDIV
#' @param FUNDIV_tree tree table formatted for FUNDIV
#' @param file.in Name of the file too save (including path)
plot_disturbance_trends <- function(FUNDIV_plot, FUNDIV_tree, file.in){
  
  # Create directory for the pplot if needed
  create_dir_if_needed(file.in)
  
  # Severity per plotcode
  sever.per.plot <- FUNDIV_tree %>%
    filter(treestatus != 1) %>%
    mutate(dead = ifelse(treestatus == 2, 0, 1)) %>%
    group_by(plotcode) %>%
    summarize(severity = sum(dead)/n())
  
  # All the countries present in the dataset
  countries.in <- unique(FUNDIV_plot$country)
  
  # Type of disturbances
  disturbances.in <- unique(subset(FUNDIV_plot, disturbance.nature != "none")$disturbance.nature)
  # Add a category for all disturbances (only if all disturbances are included)
  if("storm" %in% disturbances.in & "fire" %in% disturbances.in){
    disturbances.in <- c(disturbances.in, "all")
  }
  
  # Create a vector of colors for plotting
  color.df <- data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm", "all"), 
                         color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE", "#D90429")) %>%
    filter(disturbance %in% disturbances.in) %>%
    arrange(disturbance)
  
  # Loop on all countries
  for(i in 1:length(countries.in)){
    
    for(j in 1:length(disturbances.in)){
      
      # Initialize data for ij by keeping only country i
      data.ij <- FUNDIV_plot %>% filter(country == countries.in[i]) %>% left_join(sever.per.plot, by = "plotcode")
      
      # Identify which plots were disturbed by disturbance j
      if(disturbances.in[j] == "all") data.ij <- mutate(data.ij, disturbed = ifelse(disturbance.nature %in% disturbances.in, 1, 0))
      if(disturbances.in[j] != "all") data.ij <- mutate(data.ij, disturbed = ifelse(disturbance.nature == disturbances.in[j], 1, 0))
      
      # Make sure that there are occurrences of disturbance j in country i
      if(sum(data.ij$disturbed, na.rm = TRUE) > 50){
        
        # One data set for proportion, one for severity
        # -- Proportion
        data.p.ij <- data.ij %>%
          group_by(surveydate2) %>% 
          summarize(n.plots.per.year = n(), 
                    p = sum(disturbed)/n()) %>%
          filter(n.plots.per.year > 50) %>%
          mutate(p = ifelse(p < 0.0001, 0.0001, p), 
                 p.logit = log(p/(1 - p)), 
                 country = countries.in[i], 
                 disturbance = disturbances.in[j])
        # -- severity
        data.s.ij <- data.ij %>%
          filter(disturbed == 1) %>%
          filter(surveydate2 %in% data.p.ij$surveydate2) %>%
          group_by(surveydate2) %>% 
          summarize(s = mean(severity, na.rm = TRUE), 
                    n = n()) %>%
          mutate(s = ifelse(s < 0.0001, 0.0001, s), 
                 s = ifelse(s > 0.9999, 0.9999, s), 
                 s.logit = log(s/(1 - s)), 
                 country = countries.in[i], 
                 disturbance = disturbances.in[j])
        
        # Fit a model for proportion and for severity
        model.p.ij <- lm(p.logit ~ surveydate2, data = data.p.ij, weights = n.plots.per.year)
        model.s.ij <- lm(s.logit ~ surveydate2, data = data.s.ij, weights = n)
        
        # Put statistics in another table
        # -- Model with proportion of disturbed plots
        stat.p.ij <- data.frame(country = countries.in[i], 
                                disturbance = disturbances.in[j], 
                                Est = summary(model.p.ij)$coef[2, 1], 
                                SE = summary(model.p.ij)$coef[2, 2], 
                                Fval = anova(model.p.ij)[1, 4], 
                                pval = scales::pvalue(anova(model.p.ij)[1, 5], add_p = TRUE, accuracy = 0.01))
        # -- Model with severity of disturbed plots
        stat.s.ij <- data.frame(country = countries.in[i], 
                                disturbance = disturbances.in[j], 
                                Est = summary(model.s.ij)$coef[2, 1], 
                                SE = summary(model.s.ij)$coef[2, 2], 
                                Fval = anova(model.s.ij)[1, 4], 
                                pval = scales::pvalue(anova(model.s.ij)[1, 5], add_p = TRUE, accuracy = 0.01))
        
        # Finish formatting
        # -- Proportion
        data.p.ij <- data.p.ij %>%
          mutate(fit.logit = predict(model.p.ij, newdata = .), 
                 fit.se = predict(model.p.ij, newdata = ., se.fit = TRUE)$se.fit,
                 fit = plogis(fit.logit), 
                 fit.inf = plogis(fit.logit - fit.se), 
                 fit.sup = plogis(fit.logit + fit.se)) %>%
          dplyr::select(surveydate2, country, disturbance, p, fit, fit.inf, fit.sup, n = n.plots.per.year)
        # -- Severity
        data.s.ij <- data.s.ij %>%
          mutate(fit.logit = predict(model.s.ij, newdata = .), 
                 fit.se = predict(model.s.ij, newdata = ., se.fit = TRUE)$se.fit,
                 fit = plogis(fit.logit), 
                 fit.inf = plogis(fit.logit - fit.se), 
                 fit.sup = plogis(fit.logit + fit.se)) %>%
          dplyr::select(surveydate2, country, disturbance, s, fit, fit.inf, fit.sup, n)
        
        # Build final dataframe
        if(i == 1 & j == 1){
          data.p <- data.p.ij
          data.s <- data.s.ij
          data.p.stat <- stat.p.ij
          data.s.stat <- stat.s.ij
        } else {
          data.p <- rbind(data.p, data.p.ij)
          data.s <- rbind(data.s, data.s.ij)
          data.p.stat <- rbind(data.p.stat, stat.p.ij)
          data.s.stat <- rbind(data.s.stat, stat.s.ij)
        }
        
      }
      
    }
  }
  
  # Plot for proportions
  plot.proportion <- data.p %>%
    mutate(disturbance = factor(disturbance, levels = color.df$disturbance)) %>%
    ggplot(aes(x = surveydate2, y = p, color = disturbance, fill = disturbance, group = 1)) + 
    geom_line(size = 0.2, inherit.aes = TRUE) + 
    geom_line(aes(y = fit), inherit.aes = TRUE) + 
    geom_ribbon(aes(ymin = fit.inf, ymax = fit.sup), alpha = 0.2, colour = NA) + 
    geom_point() + 
    #geom_text(aes(y = (p+0.03), label = paste0("(", n, ")")), size = 2.5) +
    facet_grid(disturbance ~ country, scales = "free_x") +
    scale_color_manual(values = color.df$color) + 
    scale_fill_manual(values = color.df$color) + 
    geom_text(data = (data.p.stat %>% 
                        left_join((data.p %>% group_by(country) %>% summarize(surveydate2 = min(surveydate2))), 
                                  by = "country") %>%
                        mutate(p = 0.25, lab = paste0("F = ", round(Fval, digits = 1), ", ", pval))), 
              aes(label = lab), inherit.aes = TRUE, hjust = 0, color = "black") + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_text(face = "bold"),
          legend.position = "none") + 
    xlab("Date of inventory 2") + ylab("Disturbance frequency over the study period")
  
  # Plot for severity
  plot.severity <- data.s %>%
    mutate(disturbance = factor(disturbance, levels = color.df$disturbance)) %>%
    ggplot(aes(x = surveydate2, y = s, color = disturbance, fill = disturbance, group = 1)) + 
    geom_line(size = 0.2, inherit.aes = TRUE) + 
    geom_line(aes(y = fit), inherit.aes = TRUE) + 
    geom_ribbon(aes(ymin = fit.inf, ymax = fit.sup), alpha = 0.2, colour = NA) + 
    geom_point() + 
    geom_text(aes(y = (s+0.1), label = paste0("(", n, ")")), size = 2.5) +
    facet_grid(disturbance ~ country, scales = "free_x") +
    scale_color_manual(values = color.df$color) + 
    scale_fill_manual(values = color.df$color) + 
    geom_text(data = (data.s.stat %>% 
                        left_join((data.s %>% group_by(country) %>% summarize(surveydate2 = min(surveydate2))), 
                                  by = "country") %>%
                        mutate(s = 1.1, lab = paste0("F = ", round(Fval, digits = 1), ", ", pval))), 
              aes(label = lab), inherit.aes = TRUE, hjust = 0, color = "black") + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_text(face = "bold"),
          legend.position = "none") + 
    xlab("Date of inventory 2") + ylab("Mean severity of plots disturbed")
  
  # Final plot
  plot.out <- plot_grid(plot.proportion, plot.severity, nrow = 1, labels = c("(a)", "(b)"), align = "h")
  
  # Save the plot
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600, bg = "white")
  
  # return the name of the saved plot
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
    ggsave(file.in.j, plot.j, width = 17, height = 12, units = "cm", dpi = 600, bg = "white")
    
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
  ggsave(file.in, plot.out, width = 20, height = 7, units = "cm", dpi = 600, bg = "white")
  
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
    if(disturbances.in[i] %in% c("storm", "snow")) data.i <- mutate(data.i, pd = plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))
    if(!(disturbances.in[i] %in% c("storm", "snow"))) data.i <- mutate(data.i, pd = plogis(a0 + b*I*dbh.scaled^c))
    
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
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600, bg = "white")
  
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
    if(disturbances.in[i] %in% c("storm", "snow")) data.i <- mutate(data.i, pd = plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c))
    if(!(disturbances.in[i] %in% c("storm", "snow"))) data.i <- mutate(data.i, pd = plogis(a0 + b*I*dbh.scaled^c))
    
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
    ggsave(file.in.i, plot.out.i, width = 25, height = 20, units = "cm", dpi = 600, bg = "white")
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
      mutate(a1 = ifelse(disturbance %in% c("storm", "snow"), a1, NA_real_))
    
    
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
      mutate(pd = ifelse(disturbance %in% c("storm", "snow"), 
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
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600, bg = "white")
  
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
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600, bg = "white")
  
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
                                   dbh.ref = 300, logratio.ref = 0, I.ref = 0.75, 
                                   file.in){
  
  # Initialize output
  out <- list()
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
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
    if(!(disturbances.in[i] %in% c("storm", "snow"))) param_per_species.i$a1 = NA_real_
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
      mutate(pd = case_when(disturbance %in% c("storm", "snow") ~ plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c), 
                            TRUE ~ plogis(a0 + b*I*dbh.scaled^c))) %>%
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
  ggsave(file.in, plot.out, width = 18, height = 30, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}





#' Get disturbance sensitivity by species for the model with stock included
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM.
#' @param dbh.ref dbh value for which to predict values
#' @param stock.ref value of the stock (ba/ha) for which to predict values
#' @param I.ref disturbance intensity for which to predict values
#' @param file.in Name and path of the plot to save
plot_param_per_species_stock <- function(jags.model, data_jags, data_model, 
                                   dbh.ref = 300, stock.ref = 20, I.ref = 0.75, 
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
    
    # Model to scale stock
    scale_stock <- lm(stock.scaled ~ stock, 
                         data = data.frame(stock = data_model[[i]]$stock, 
                                           stock.scaled = data_jags[[i]]$data_jags$stock))
    # Model to scale dbh
    scale_dbh <- lm(dbh.scaled ~ dbh, 
                    data = data.frame(dbh = data_model[[i]]$dbh, 
                                      dbh.scaled = data_jags[[i]]$data_jags$dbh))
    
    # Parameter per species
    param_per_species.i <- extract_param_per_species(jags.model[[i]], data_jags[[i]])
    param_per_species.i <- param_per_species.i %>%
      group_by(species, iter) %>%
      summarise(a0 = mean(a0), a1 = mean(a1), b = mean(b), c = mean(c))
    
    # Format data 
    data.i <- expand.grid(species = data_jags[[i]]$species.table$species,
                          dbh = dbh.ref, stock = stock.ref, I = I.ref, 
                          iter = unique(param_per_species.i$iter), 
                          disturbance = disturbances.in[i]) %>%
      # Scale variables when needed
      mutate(dbh.scaled = predict(scale_dbh, newdata = .), 
             stock.scaled = predict(scale_stock, newdata = .)) %>%
      # Add parameters per species
      left_join(param_per_species.i, by = c("species", "iter")) %>%
      # Compute probabilities
      mutate(pd = plogis(a0 + a1*stock.scaled + b*I*dbh.scaled^c)) %>%
      # Average per species
      group_by(species) %>%
      summarise(sensitivity_mean = mean(pd), 
                sensitivity_ql = as.numeric(quantile(pd, probs = 0.025)), 
                sensitivity_qh = as.numeric(quantile(pd, probs = 0.975)), 
                dbh.effect_mean = mean(c), 
                dbh.effect_ql = as.numeric(quantile(c, probs = 0.025)), 
                dbh.effect_qh = as.numeric(quantile(c, probs = 0.975)), 
                stock.effect_mean = mean(a1, na.rm = TRUE), 
                stock.effect_ql = as.numeric(quantile(a1, probs = 0.025, na.rm = TRUE)), 
                stock.effect_qh = as.numeric(quantile(a1, probs = 0.975, na.rm = TRUE))) %>%
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
                                           "dbh effect", "stock effect")))
    
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
  ggsave(file.in, plot.out, width = 18, height = 30, units = "cm", dpi = 600, bg = "white")
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
  ggsave(file.in, plot.out, width = 14, height = 11, units = "cm", dpi = 600, bg = "white")
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
    ggsave(file.in, plot.out, width = 30, height = 10, units = "cm", dpi = 600, bg = "white")
  }
  
  return(file.in)
  
}


#' Function to plot results of the climate analysis 
#' @param gbif_file Name of the file containing species mean climate
#' @param gbif_disturbance_file Name of the file containing species mean disturbance index
#' @param disturbance_sensitivity list of dataset containing disturbance sensitivity per species
#' @param file.in Where to save the plot
plot_rda_climate <- function(disturbance_sensitivity, gbif_file, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  # Compile trait dataset
  traits <- fread(gbif_file) %>% dplyr::select(species, mat, tmin, map)
  
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
    
    # Make rda
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
    plot.rda <- data.frame(trait = c("map", "mat", "tmin"), 
                           position = c(1:3)) %>%
      left_join(data.frame(trait = as.character(rownames(summary(rda)$species)), 
                           rda1 = summary(rda)$species[, 1]), by = "trait") %>%
      mutate(trait = factor(trait, levels = c("map", "tmin", "mat", "windspeed", "fwi", "swe"))) %>%
      ggplot(aes(x = 0, xend = rda1, y = position, yend = position)) + 
      geom_segment(arrow = arrow(length = unit(0.1, "cm")), type = "closed") + 
      scale_y_continuous(breaks = c(1:3), 
                         labels = c("map", "mat", "tmin"), 
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
    plot.i <- plot_grid(plot.rda, plot.model, nrow = 2, rel_heights = c(0.7, 1), align = "v")
    
    # Add to the output list
    eval(parse(text = paste0("plots.out$", disturbances.in[i], " <- plot.i")))
  }
  
   
  
  
  
  
  
  ## - Plot a pca of species position along the climatic variables tested
  
  # - Climatic values per species
  gbif_climate <- fread(gbif_file) %>%
    dplyr::select(species, mat, tmin, map) %>%
    drop_na()
  
  # - Make PCA 
  pca <- prcomp((gbif_climate %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # - Extract the coordinates of the ndividuals on pca axis
  res.ind <- data.frame(species = gbif_climate$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2]) %>%
    mutate(sp = paste(substr(gsub("\\ .+", "", species), 1, 2), 
                      substr(gsub(".+\\ ", "", species), 1, 2), 
                      sep = "."))
  
  # - Extract the coordinates of the variables on pca axis
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2])
  
  # - Minimum and maximum in each pca axis
  pca.xmin <- -max(abs(res.ind$pca1))
  pca.xmax <- max(abs(res.ind$pca1))
  pca.ymin <- -max(abs(res.ind$pca2))
  pca.ymax <- max(abs(res.ind$pca2))
  
  # Make the plot
  plot.pca <- res.ind %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_point(fill = "#023E8A", color = "black", shape = 21) +
    geom_text(aes(label = sp), nudge_y = 0.1, color = "#023E8A", size = 3) +
    geom_segment(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
                 aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.1, "cm")), 
                 type = "closed", color = "#D90429") + 
    geom_text(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
              aes(label = var), color = "#D90429", size = 5, 
              nudge_x = ifelse(res.var$pca1 < 0, pca.xmin/12, pca.xmax/12)) +
    geom_hline(size = 0.2, yintercept = 0, color = "#6C757D", linetype = "dashed") + 
    geom_vline(size = 0.2, xintercept = 0, color = "#6C757D", linetype = "dashed") + 
    xlim((pca.xmin-0.2), (pca.xmax+0.2)) + 
    ylim((pca.ymin-0.2), (pca.ymax+0.2)) +
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank())
  

  
  ## - Assemble pca and rda plots and adjust dimensions depending on the number of disturbances
  
  # If more than three disturbances
  if(length(disturbances.in) > 3){
    # assemble rda plots
    plot.rda <- plot_grid(plotlist = plots.out, nrow = 2, align = "hv", 
                          labels = paste0("(", letters[c(2:(length(disturbances.in)+1))], ")"))
    # Add pca plot
    plot.out <- plot_grid(plot_grid((ggplot() + theme_void()), plot.pca, (ggplot() + theme_void()), 
                                    ncol = 1, rel_heights = c(0.1, 1, 0.4)), 
                          (ggplot() + theme_void()),
                          plot.rda, nrow = 1, labels = c("(a)", "", ""), rel_widths = c(0.6, 0.1, 1))
    # Save the plot
    ggsave(file.in, plot.out, width = 35, height = 20, units = "cm", dpi = 600)
  }
  
  # If up to three disturbances
  if(length(disturbances.in) <= 3){
    # assemble rda plots
    plot.rda <- plot_grid(plotlist = plots.out, nrow = 1, align = "hv", 
                          labels = paste0("(", letters[c(2:(length(disturbances.in)+1))], ")"))
    # Add pca plot
    plot.out <- plot_grid(plot.pca, plot.rda, nrow = 1, align = "v", 
                          labels = c("(a)", ""), rel_widths = c(0.37, 1))
    # Save the plot
    ggsave(file.in, plot.out, width = 35, height = 13, units = "cm", dpi = 600, bg = "white")
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
    ggsave(file.in.i, plot.i, width = 14, height = 11, units = "cm", dpi = 600, bg = "white")
    
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
## Plots for the manuscript ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




#' Function to plot results of the climate analysis for the manuscript
#' @param gbif_file Name of the file containing species mean climate
#' @param gbif_disturbance_file Name of the file containing species mean disturbance index
#' @param disturbance_sensitivity list of dataset containing disturbance sensitivity per species
#' @param disturbance_sensitivity_bis list of dataset containing disturbance sensitivity to snow and biotic
#' @param file.in Where to save the plot
plot_disturbance_climate_ms <- function(disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
                                        gbif_disturbance_file, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- c(disturbance_sensitivity_full, disturbance_sensitivity_full_bis[c("snow", "biotic")])
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  
  # Initialize plot list
  plots.disturbance <- list()
  
  # Couples disturbance - index
  disturbance.index <- data.frame(disturbance = c("snow", "storm", "fire"), 
                                  index = c("swe", "windspeed", "fwi"), 
                                  name = c("Snow Water \n Equivalent", "Mean windspeed", "Fire Weather \n Index"), 
                                  color = c("#006D77", "#4361EE", "#F77F00"))
  
  # Loop on all couples
  for(j in 1:dim(disturbance.index)[1]){
    # Data for the model
    data.model.j <- as.data.frame(disturbance_sensitivity.in[disturbance.index$disturbance[j]])
    colnames(data.model.j) <- c("species", "iter", "p")
    data.model.j <- data.model.j %>%
      left_join(fread(gbif_disturbance_file) %>% dplyr::select("species", "index" = disturbance.index$index[j]), 
                by = "species") %>%
      drop_na() %>%
      mutate(sensitivity.logit = log(p/(1 - p))) %>%
      group_by(species, index) %>%
      summarize(p_025 = quantile(p, probs = 0.025), 
                p_975 = quantile(p, probs = 0.975), 
                p = mean(p), 
                p.logit = mean(sensitivity.logit), 
                w = 1/var(sensitivity.logit))
    
    # model
    model.j <- lm(p.logit ~ index, weights = w, data = data.model.j)
    # Data with predictions
    data.fit.j <- data.frame(index = c(round(min(data.model.j$index)*100, digits = 0):
                                         round(max(data.model.j$index)*100, digits = 0)/100)) %>%
      mutate(fit.logit = predict(model.j, newdata = .), 
             fit.lwr = predict(model.j, newdata = ., interval = "confidence")[, 2],
             fit.upr = predict(model.j, newdata = ., interval = "confidence")[, 3],
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.lwr), 
             fit.sup = plogis(fit.upr), 
             p = NA_real_)
    
    
    # Plot
    plot.j <- data.model.j %>%
      ggplot(aes(x = index, y = p, group = 1)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = disturbance.index$color[j], color = "#343A40") + 
      xlab(disturbance.index$name[j]) + 
      ylab(paste0("Sensitivity to \n", disturbance.index$disturbance[j])) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            plot.title = element_text(size = 17, face = "italic"), 
            axis.title = element_text(size = 17)) + 
      scale_y_continuous(breaks = c(0:5)*0.2) + 
      ggtitle(paste0("F = ", round(anova(model.j)[1, 4], digits = 1), ", ",
                     scales::pvalue(anova(model.j)[1, 5], add_p = TRUE, accuracy = 0.01)))
    
    # Add line and confidence interval only if the regression is significant
    if(anova(model.j)[1, 5] <= 0.05){
      plot.j <- plot.j  + 
        geom_line(data = data.fit.j, aes(y = fit), inherit.aes = TRUE, color = disturbance.index$color[j])  + 
        geom_ribbon(data = data.fit.j, aes(ymin = fit.inf, ymax = fit.sup), 
                    alpha = 0.5, fill = disturbance.index$color[j], inherit.aes = TRUE)
      
    }
    # Add to the output list
    eval(parse(text = paste0("plots.disturbance$", disturbance.index$index[j], " <- plot.j")))
  }
  
  
  # Final plot
  plot.out <- plot_grid(plotlist = plots.disturbance[c("windspeed", "fwi", "swe")], nrow = 1, align = "h", 
                        labels = paste0("(", letters[c(1:3)], ")"), scale = 0.9)
  
  
  # Save the plot
  ggsave(file.in, plot.out, width = 25, height = 8, units = "cm", dpi = 600, bg = "white")
  
  return(file.in)
  
}



#' Function to plot results of the climate analysis for the manuscript
#' @param gbif_file Name of the file containing species mean climate
#' @param disturbance_sensitivity list of dataset containing disturbance sensitivity per species
#' @param disturbance_sensitivity_bis list of dataset containing disturbance sensitivity to snow and biotic
#' @param file.in Where to save the plot
plot_rda_climate_ms <- function(disturbance_sensitivity, disturbance_sensitivity_bis, 
                                gbif_file, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- c(disturbance_sensitivity, disturbance_sensitivity_bis[c("snow", "biotic")])
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  # Compile trait dataset
  traits <- fread(gbif_file) %>% dplyr::select(species, mat, tmin, map)
  
  
  ## - Initialize output
  plots.out <- list()
  
  ## - Loop on all type of disturbances
  for(i in 1:length(disturbances.in)){
    
    # Remove na in trait table
    traits.i <- traits %>%
      # Add sensitivity
      left_join(disturbance_sensitivity.in[[i]], by = "species") %>%
      mutate(p.logit = log(p/(1 - p))) %>%
      drop_na()
    
    # scale trait dataset
    traits.scaled.i <- scale(as.matrix(traits.i %>% dplyr::select(colnames(.)[!(colnames(.) %in% c("species", "p", "p_975", "p_025"))])))
    
    # Make rda
    rda <- rda(traits.scaled.i ~ p.logit, traits.i,  scale = FALSE)
    
    # Test the significance of the rda model
    rda.test <- anova(rda, permutations = 1000, by = "axis")
    
    
    # Results of the rda
    res.rda.ind <- data.frame(species = traits.i$species, 
                              rda1 = as.numeric(summary(rda)$sites[, 1])) %>%
      mutate(species.short = paste0(substr(species, 1, 2), substr(gsub(".+\\ ", "", species), 1, 2)))
    
    # Identify maximum coordinates of species in rda
    max.rda1 <- max(sqrt(res.rda.ind$rda1^2)) + 0.1
    
    # Data for model
    data.model <- disturbance_sensitivity.in[[i]] %>%
      # Logit transformation
      mutate(sensitivity.logit = log(p/(1 - p)), 
             w = 1/(p_975 - p_025)) %>%
      # Add pca
      left_join(res.rda.ind, by = "species") %>%
      drop_na()
    
    
    # Make regression sensitivity vs rda
    model.rda1 <- lm(sensitivity.logit ~ rda1, data = data.model)
    
    # Data for the predictions
    data.fit <- data.frame(rda1 = c(round(-max.rda1*1000, digits = 0):round(max.rda1*1000, digits = 0))/1000) %>%
      mutate(fit.logit = predict(model.rda1, newdata = .), 
             fit.lwr = predict(model.rda1, newdata = ., interval = "confidence")[, 2],
             fit.upr = predict(model.rda1, newdata = ., interval = "confidence")[, 3],
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.lwr), 
             fit.sup = plogis(fit.upr), 
             p = NA_real_)
    
    # Plot predictions rda1
    plot.model <- data.model %>%
      ggplot(aes(x = rda1, y = p, group = 1)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = color.in[i], color = "#343A40") +
      xlab(paste0("RDA1 (", round(summary(rda)$cont$importance[2, 1]*100, digits = 2), "%)")) + 
      ylab(expression(atop("Disturbance", "sensitivity"))) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            axis.title = element_text(size = 16), 
            axis.text = element_text(size =10)) + 
      xlim(-max.rda1, max.rda1) + 
      scale_y_continuous(breaks = c(0:5)*0.2, limits = c(0, 1))
    
    # Add line and confidence interval only if the regression is significant
    if(rda.test[1, 4] <= 0.05){
      plot.model <- plot.model  + 
        geom_line(data = data.fit, aes(y = fit), inherit.aes = TRUE, color = color.in[i]) + 
        geom_ribbon(data = data.fit, aes(ymin = fit.inf, ymax = fit.sup), 
                    alpha = 0.5, fill = color.in[i], inherit.aes = TRUE)
    }
    
    # Plot rda
    plot.rda <- data.frame(trait = c("map", "mat", "tmin"), 
                           position = c(1:3)) %>%
      left_join(data.frame(trait = as.character(rownames(summary(rda)$species)), 
                           rda1 = summary(rda)$species[, 1]), by = "trait") %>%
      mutate(trait = factor(trait, levels = c("map", "tmin", "mat", "windspeed", "fwi", "swe"))) %>%
      ggplot(aes(x = 0, xend = rda1, y = position, yend = position)) + 
      geom_segment(arrow = arrow(length = unit(0.1, "cm")), type = "closed") + 
      scale_y_continuous(breaks = c(1:3), 
                         labels = c("map", "mat", "tmin"), 
                         expand = c(0.25, 0)) +
      geom_vline(xintercept = 0, color = "darkgray") + 
      ylab("") + xlim(-max.rda1, max.rda1) + xlab("") +
      labs(title = paste0(toupper(substr(disturbances.in[i], 1, 1)), 
                          substr(disturbances.in[i], 2, nchar(disturbances.in[i]))), 
           subtitle = paste0("F = ", round(rda.test[1, 3], digits = 1), ", ",
                             scales::pvalue(rda.test[1, 4], add_p = TRUE, accuracy = 0.01))) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            plot.title = element_text(size = 21), 
            axis.text.y = element_text(size = 15), 
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(),
            plot.subtitle = element_text(size = 17, face = "italic"))
    
    # Adjust axis and titles depending on graph position
    if(i > 1){
      plot.rda <- plot.rda + theme(axis.text.y = element_blank(), 
                                   axis.ticks.y = element_blank())
      plot.model <- plot.model + ylab("") 
    } 
    
    
    # Plot all graphs together
    plot.i <- plot_grid(plot.rda, plot.model, nrow = 2, rel_heights = c(0.7, 1), align = "v")
    
    # Add to the output list
    eval(parse(text = paste0("plots.out$", disturbances.in[i], " <- plot.i")))
  }
  
  # Add an empty plot to manage spacing between graphs
  plots.out$empty <- ggplot() + theme_void()
  
  # Final rda plot
  plot.out <- plot_grid(plotlist = plots.out[c("storm", "fire", "other", "empty", "biotic", "snow")], nrow = 1, 
                        align = "hv", rel_widths = c(1.2, 1, 1, 0.25, 1, 1),
                        labels = c("(a)", "", "", "(b)", "", ""), label_size = 19)
  
  # Save the plot
  ggsave(file.in, plot.out, width = 35, height = 10, units = "cm", dpi = 600, bg = "white")
  
  return(file.in)
  
}



#' Function to plot pca of three climatic variables for each species
#' @param gbif_file Name of the file containing species mean climate
#' @param file.in Where to save the plot
plot_pca_ms <- function(gbif_file, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  
  # - Climatic values per species
  gbif_climate <- fread(gbif_file) %>%
    dplyr::select(species, mat, tmin, map) %>%
    drop_na()
  
  # - Make PCA 
  pca <- prcomp((gbif_climate %>% dplyr::select(-species)), 
                center = T, scale = T)
  
  # - Extract the coordinates of the ndividuals on pca axis
  res.ind <- data.frame(species = gbif_climate$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2]) %>%
    mutate(sp = paste(substr(gsub("\\ .+", "", species), 1, 2), 
                      substr(gsub(".+\\ ", "", species), 1, 2), 
                      sep = "."))
  
  # - Extract the coordinates of the variables on pca axis
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2])
  
  # - Minimum and maximum in each pca axis
  pca.xmin <- -max(abs(res.ind$pca1))
  pca.xmax <- max(abs(res.ind$pca1))
  pca.ymin <- -max(abs(res.ind$pca2))
  pca.ymax <- max(abs(res.ind$pca2))
  
  # Make the plot
  plot.out <- res.ind %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_point(fill = "#023E8A", color = "black", shape = 21) +
    geom_text(aes(label = sp), nudge_y = 0.1, color = "#023E8A", size = 3) +
    geom_segment(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
                 aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.1, "cm")), 
                 type = "closed", color = "#D90429") + 
    geom_text(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
              aes(label = var), color = "#D90429", size = 5, 
              nudge_x = ifelse(res.var$pca1 < 0, pca.xmin/12, pca.xmax/12)) +
    geom_hline(size = 0.2, yintercept = 0, color = "#6C757D", linetype = "dashed") + 
    geom_vline(size = 0.2, xintercept = 0, color = "#6C757D", linetype = "dashed") + 
    xlim((pca.xmin-0.2), (pca.xmax+0.2)) + 
    ylim((pca.ymin-0.2), (pca.ymax+0.2)) +
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 15))
  
  
  # Save the plot
  ggsave(file.in, plot.out, width = 14, height = 14, units = "cm", dpi = 600, bg = "white")
  
  return(file.in)
  
}


#' Plot correlation matrix between the sensitivity to each disturbance
#' @param disturbance_sensitivity dataset containing disturbance sensitivity per species
#' @param disturbance_sensitivity_bis list of dataset containing disturbance sensitivity to snow and biotic
#' @param file.in Name of the file to save (including path)
plot_correlation_disturbance <- function(disturbance_sensitivity, disturbance_sensitivity_bis, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- c(disturbance_sensitivity, disturbance_sensitivity_bis)
  
  ## - Loop on all disturbance types to assemble data sets
  for(i in 1:length(names(disturbance_sensitivity.in))){
    data.in.i <- as.data.frame(disturbance_sensitivity.in[[i]]) %>%
      mutate(disturbance = names(disturbance_sensitivity.in)[i])
    if(i == 1) data.in <- data.in.i
    if(i > 1) data.in <- rbind(data.in, data.in.i)
  }
  
  ## - Format the dataset to get one column per disturbance, and one line per species
  data.in <- data.in %>% 
    dplyr::select(species, disturbance, p) %>%
    spread(key = disturbance, value = p)
  
  ## - Make the plot
  plot.out <- ggpairs(data.in %>% dplyr::select(-species)) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), 
          panel.grid = element_blank(), 
          strip.background = element_blank()) 
  
  ## - Export plot
  ggsave(file.in, plot.out, width = 15, height = 15, units = "cm", dpi = 600, bg = "white")
  
  return(file.in)
}



#' Plot to verify that species sensitivity is not related to the mean intensity or to the number of trees exposed
#' @param disturbance_sensitivity.in dataset containing disturbance sensitivity per species
#' @param jags.model.in object generated by a function form the script functions_analysis.R
#' @param data_jags.in input data for the model jags.model.in
#' @param file.in Name of the file to save (path included)
plot_sensitivity_vs_intensity_and_ntrees <- function(disturbance_sensitivity.in, jags.model.in, data_jags.in, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Name of the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  # Initialize plot lists
  plots.intensity <- list()
  plots.ntrees <- list()
  
  # Loop on all distrubance types
  for(i in 1:length(disturbances.in)){
    
    # intensity per plotcode 
    intensity.i <- extract_intensity_per_plotcode(jags.model.in[[i]], data_jags.in[[i]]) %>%
      gather(key = "iteration", value = "intensity", colnames(.)[c(2:dim(.)[2])]) %>%
      group_by(plotcode) %>%
      summarize(intensity = mean(intensity))
    
    # Calculate mean intensity per species and cross with disturbance sensitivity
    data.i <- data.frame(sp = data_jags.in[[i]]$data_jags$sp, 
                         plot = data_jags.in[[i]]$data_jags$plot) %>%
      left_join(data_jags.in[[i]]$plotcode.table, by = "plot") %>%
      left_join(data_jags.in[[i]]$species.table, by = "sp") %>%
      group_by(plotcode, species) %>%
      summarize(n = n()) %>%
      ungroup() %>%
      group_by(plotcode) %>%
      mutate(rate = n/sum(n)) %>%
      left_join(intensity.i, by = "plotcode") %>%
      ungroup() %>%
      group_by(species) %>%
      summarize(intensity.mean = sum(rate*intensity)/sum(rate), 
                n.tree.exposed = sum(n)) %>%
      left_join(disturbance_sensitivity.in[[i]], by = "species") %>%
      mutate(sensitivity.logit = log(p/(1 - p)), 
             w = 1/(p_975 - p_025))
    
    # Models of disturbance sensitivity against intensity and number of trees
    model.intensity.i <- lm(sensitivity.logit ~ intensity.mean, weights = w, data = data.i)
    model.ntrees.i <- lm(sensitivity.logit ~ n.tree.exposed, weights = w, data = data.i)
    
    # Data for the predictions of the intensity model
    data.fit.intensity <- data.frame(intensity.mean = c(round(min(data.i$intensity.mean)*100, digits = 0):
                                                          round(max(data.i$intensity.mean)*100, digits = 0))/100) %>%
      mutate(fit.intensity.logit = predict(model.intensity.i, newdata = .)) %>%
      mutate(fit.intensity.se = predict(model.intensity.i, newdata = ., se.fit = TRUE)$se.fit) %>%
      mutate(fit.intensity = plogis(fit.intensity.logit)) %>%
      mutate(fit.intensity.inf = plogis(fit.intensity.logit - fit.intensity.se)) %>%
      mutate(fit.intensity.sup = plogis(fit.intensity.logit + fit.intensity.se)) %>%
      mutate(p = NA_real_)
    
    # Data for the predictions of the ntrees model
    data.fit.ntrees <- data.frame(n.tree.exposed = c(min(data.i$n.tree.exposed):max(data.i$n.tree.exposed))) %>%
      mutate(fit.ntrees.logit = predict(model.ntrees.i, newdata = .)) %>%
      mutate(fit.ntrees.se = predict(model.ntrees.i, newdata = ., se.fit = TRUE)$se.fit) %>%
      mutate(fit.ntrees = plogis(fit.ntrees.logit)) %>%
      mutate(fit.ntrees.inf = plogis(fit.ntrees.logit - fit.ntrees.se)) %>%
      mutate(fit.ntrees.sup = plogis(fit.ntrees.logit + fit.ntrees.se)) %>%
      mutate(p = NA_real_)
    
    # Plot for intensity
    plot.intensity.i <- data.i %>%
      ggplot(aes(x = intensity.mean, y = p, group = 1)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = color.in[i], color = "#343A40") + 
      geom_line(data = data.fit.intensity, aes(y = fit.intensity), inherit.aes = TRUE, color = color.in[i]) + 
      geom_ribbon(data = data.fit.intensity, aes(ymin = fit.intensity.inf, ymax = fit.intensity.sup), 
                  alpha = 0.5, fill = color.in[i], inherit.aes = TRUE) +
      xlab("Mean disturbance intensity") + 
      ylab("Disturbance sensitivity") +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank()) + 
      scale_y_continuous(breaks = c(0:5)*0.2) +
      labs(title = paste0(toupper(substr(disturbances.in[i], 1, 1)), 
                          substr(disturbances.in[i], 2, nchar(disturbances.in[i]))), 
           subtitle = paste0("F = ", round(anova(model.intensity.i)[1, 4], digits = 1), ", ",
                             scales::pvalue(anova(model.intensity.i)[1, 5], add_p = TRUE, accuracy = 0.01)))
    
    # Plot for ntrees
    plot.ntrees.i <- data.i %>%
      ggplot(aes(x = n.tree.exposed, y = p, group = 1)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975), width = 0, color = "#343A40") +
      geom_point(size = 2, shape = 21, fill = color.in[i], color = "#343A40") + 
      geom_line(data = data.fit.ntrees, aes(y = fit.ntrees), inherit.aes = TRUE, color = color.in[i]) + 
      geom_ribbon(data = data.fit.ntrees, aes(ymin = fit.ntrees.inf, ymax = fit.ntrees.sup), 
                  alpha = 0.5, fill = color.in[i], inherit.aes = TRUE) +
      xlab("Number of trees exposed") + 
      ylab("Disturbance sensitivity") +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank()) + 
      scale_y_continuous(breaks = c(0:5)*0.2) +
      labs(title = paste0(toupper(substr(disturbances.in[i], 1, 1)), 
                          substr(disturbances.in[i], 2, nchar(disturbances.in[i]))), 
           subtitle = paste0("F = ", round(anova(model.ntrees.i)[1, 4], digits = 1), ", ",
                             scales::pvalue(anova(model.ntrees.i)[1, 5], add_p = TRUE, accuracy = 0.01)))
    
    # Add to the plot lists
    eval(parse(text = paste0("plots.intensity$", disturbances.in[i], " <- plot.intensity.i")))
    eval(parse(text = paste0("plots.ntrees$", disturbances.in[i], " <- plot.ntrees.i")))
    
  }
  
  # Assemble all plots
  plot.out <- plot_grid(plot_grid(plotlist = plots.intensity, nrow = 1, align = "hv", scale = 0.9), 
                        plot_grid(plotlist = plots.ntrees, nrow = 1, align = "hv", scale = 0.9), 
                        nrow = 2, align = "hv", labels = c("(a)", "(b)"))
  
  # Save the plot
  ggsave(file.in, plot.out, width = 30, height = 13, units = "cm", dpi = 600, bg = "white")
  
  return(file.in)
}




#' Figure of parameter value per species for the manuscript
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param data_model Tree data formatted for the IPM.
#' @param dbh.ref dbh value for which to predict values
#' @param logratio.ref value of the ratio dbh/dqm for which to predict values
#' @param I.ref disturbance intensity for which to predict values
#' @param file.in Name and path of the plot to save
plot_param_per_species_ms <- function(jags.model, data_jags, data_model, 
                                      dbh.ref = 300, logratio.ref = 0, I.ref = 0.75, 
                                      file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  # Initialize output
  out <- list()
  
  # Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Create a vector of colors for plotting
  color.vector <- (data.frame(disturbance = disturbances.in) %>%
                     left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                               by = "disturbance"))$color
  
  # Create a vector of colors for the plot legend
  color.legend <- (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                     left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                          color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                               by = "disturbance"))$color
  
  # Create a legend for the plot
  plot.legend <- cowplot::get_legend(
    data.frame(x = c(1:5), y = c(1:5), ymin = c(0:4), ymax = c(2:6), 
               disturbance = factor(disturbances.in, 
                                    levels = c("storm", "fire", "other", "biotic", "snow"))) %>%
      ggplot(aes(x = x, y = y, color = disturbance)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) + 
      scale_color_manual(values = color.legend) + 
      theme(legend.title = element_blank(), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 19))
  )
  
  # Loop on the three factor levels
  for(j in 1:3){
    
    # Level of factor j
    factor.j <- c("species sensitivity", "dbh effect", "dominance effect")[j]
    
    # Initialize the list that will contain the plots of plot i
    plot.list.j <- list()
    
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
      if(!(disturbances.in[i] %in% c("storm", "snow"))) param_per_species.i$a1 = NA_real_
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
        mutate(pd = case_when(disturbance %in% c("storm", "snow") ~ plogis(a0 + a1*logratio.scaled + b*I*dbh.scaled^c), 
                              TRUE ~ plogis(a0 + b*I*dbh.scaled^c))) %>%
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
        mutate(parameter = ifelse(parameter != "sensitivity", parameter, "species sensitivity")) %>%
        mutate(parameter = factor(parameter, levels = c("species sensitivity", "dbh effect", "dominance effect")))
      
      # Case where the plot should be empty
      if(factor.j == "dominance effect" & disturbances.in[i] %in% c("other", "fire", "biotic")){
        
        # Make an empty plot, or use the legend if graph at the top right
        if(j == 3 & disturbances.in[i] == "other") plot.ij <- plot.legend
        else plot.ij <- ggplot() + theme_void() 
        
      }else{
        
        # Build plot ij
        plot.ij <- data.i %>%
          mutate(parameter = as.character(parameter)) %>%
          filter(parameter == factor.j) %>%
          ggplot(aes(x = species, y = mean, ymin = ql, ymax = qh)) +
          geom_errorbar(width = 0, color = color.vector[i]) +
          geom_point(color = color.vector[i]) + 
          coord_flip() + 
          facet_wrap(~ parameter, scales = "free_x") + 
          xlab("") + ylab("") +
          theme(panel.background = element_rect(color = "black", fill = "white"), 
                panel.grid = element_blank(), 
                strip.background = element_blank(), 
                axis.text.y = element_text(size = 9, face = "italic"), 
                axis.text.x = element_text(size = 7), 
                strip.text = element_text(size = 16), 
                plot.margin = unit(c(0, 0, 0, 0), "cm"))
        
        # If dbh or dominance effect, remove axis ticks and add a vertical line for 0
        if(j > 1) plot.ij <- plot.ij + 
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + 
            geom_hline(yintercept = 0, linetype = "dashed") +
            ylim(-2.3, 2.3)
        
        # If sensitivity, scale between 0 and 1
        if(j == 1) plot.ij <- plot.ij + ylim(0, 1)
        
        # If not storm disturbance (positioned on top), remove strip title
        if(disturbances.in[i] != "storm") plot.ij <- plot.ij + theme(strip.text = element_blank())
        
        # If not snow disturbance (positioned at bottom), remove x-axis text and ticks
        if(disturbances.in[i] != "snow") plot.ij <- plot.ij + 
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        
      }
      
      # Add to the list
      eval(parse(text = paste0("plot.list.j$", disturbances.in[i], " <- plot.ij")))
      
      # Add also the number of species
      n.sp.per.dist <- c(n.sp.per.dist, length(unique(data.i$species)))
      
    }
    
    # Assemble plots ij, add label only if first column
    if(j == 1){
      plot.j <- plot_grid(plotlist = plot.list.j[c("storm", "fire", "other", "biotic", "snow")], 
                          ncol = 1, align = "v", labels = c("(a)", "", "", "(b)", ""),
                          rel_heights = (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                                           left_join(data.frame(disturbance = disturbances.in, h = (n.sp.per.dist + 10)), 
                                                     by = "disturbance"))$h)
    }else{
      plot.j <- plot_grid(plotlist = plot.list.j[c("storm", "fire", "other", "biotic", "snow")], 
                          ncol = 1, align = "v", 
                          rel_heights = (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                                           left_join(data.frame(disturbance = disturbances.in, h = (n.sp.per.dist + 10)), 
                                                     by = "disturbance"))$h)
    }
    
    # Add to the output list
    eval(parse(text = paste0("out$", gsub("\\ ", "\\.", factor.j), " <- plot.j")))
    
  }
  
  
  # Assemble all plots
  plot.out <- plot_grid(plotlist = out, nrow = 1, align = "h", rel_widths = c(1.5, 1, 1))
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 27, height = 30, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}




#' Plot the location of each plot and the nature of the disturbance
#' @param FUNDIV_plot Plot table formatted for FUNDIV
#' @param FUNDIV_plot_bis Plot table formatted for FUNDIV with only biotic and snow
#' @param file.in Path and file where to save the plot
map_disturbances_ms <- function(FUNDIV_plot, FUNDIV_plot_bis, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - Format data for plotting
  data <- FUNDIV_plot %>%
    filter(!(plotcode %in% FUNDIV_plot_bis$plotcode)) %>%
    rbind(FUNDIV_plot_bis) %>%
    filter(disturbance.nature %in% c("biotic", "storm", "snow", "fire", "other")) %>%
    mutate(disturbance = factor(disturbance.nature, 
                                levels = c("biotic", "fire", "other", "snow", "storm"))) %>%
    dplyr::select(plotcode, latitude, longitude, disturbance)  %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
  
  
  ## - Final plot
  plot.out <- ne_countries(scale = "medium", returnclass = "sf") %>%
    mutate(keep = ifelse(sovereignt %in% c("France", "Spain", "Finland"), "yes", "no")) %>%
    ggplot(aes(geometry = geometry)) +
    geom_sf(aes(fill = keep), color = "white", show.legend = F) + 
    scale_fill_manual(values = c("#8D99AE", "#343A40")) +
    geom_sf(data = data, shape = 16, aes(color = disturbance), size = 1.3) +
    scale_color_manual(values = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) +
    coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) +
    guides(fill = FALSE, 
           colour = guide_legend(override.aes = list(size=12))) +
    annotation_scale(location = "br", width_hint = 0.2) + 
    theme(panel.background = element_rect(color = 'black', fill = 'white'), 
          panel.grid = element_blank(),
          axis.text = element_text(size = 35),
          legend.text = element_text(size = 45),
          legend.title = element_blank(),
          legend.position = c(0.1, 0.9), 
          legend.key = element_blank())
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 33.8, height = 42.9, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}


#' Plot senstivity to all disturbances against traits on one multipanel
#' @param traits dataset containing trait values per species
#' @param traits_TRY dataset containing trait values per species from TRY
#' @param disturbance_sensitivity dataset containing disturbance sensitivity per species
#' @param disturbance_sensitivity_bis list of dataset containing disturbance sensitivity to snow and biotic
#' @param file.in Name of the file to save (including path)
plot_traits_vs_sensitivity_ms <- function(traits, traits_TRY, disturbance_sensitivity, 
                                          disturbance_sensitivity_bis, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Initialize plot list
  plots.out <- list()
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- c(disturbance_sensitivity, disturbance_sensitivity_bis[c("snow", "biotic")])
  
  ## - Loop on all disturbance types to assemble data sets
  for(i in 1:length(names(disturbance_sensitivity.in))){
    data.in.i <- as.data.frame(disturbance_sensitivity.in[[i]]) %>%
      mutate(disturbance = names(disturbance_sensitivity.in)[i])
    if(i == 1) data.in <- data.in.i
    if(i > 1) data.in <- rbind(data.in, data.in.i)
  }
  
  ## - Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood density" = "wood.density_g.cm3", 
      "Shade tolerance" = "shade.tolerance", 
      "Root mass fraction" = "Root_mass_fraction", 
      "Bark thickness" = "bark.thickness_mm", 
      "H to dbh ratio" = "height.dbh.ratio", 
      "Leaf CN ratio" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf NP ratio" = "TRY_leaf.NP.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf Pmass" = "TRY_leaf.P.mass_mg.g", 
      "SLA" = "TRY_leaf.sla_mm2mg-1", 
      "Leaf thickness" = "TRY_leaf.thickness_mm", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Stomata conductance" = "TRY_stomata.conductance_millimolm-2s-1", 
      "Maximum growth" = "growth.max"
    )
  
  
  # Loop on all traits
  for(trait.i in colnames(traits.in)[which(colnames(traits.in) != "species")]){
    
    # Format data for the model
    data.i <- data.in %>%
      left_join((traits.in %>% dplyr::select("species", "trait" = trait.i)), 
                by = "species") %>%
      # Logit transformation and weight
      mutate(sensitivity.logit = log(p/(1 - p)), 
             w = 1/(p_975 - p_025)) %>%
      # remove NA
      drop_na()
    
    # Fit a model
    model.i <- lmer(sensitivity.logit ~ trait + (1|species), weights = w, data = data.i)
    
    # Data with the predictions of the model
    data.fit <- data.frame(
      trait = c(round(min(data.i$trait)*100, digits = 0):round(max(data.i$trait)*100, digits = 0))/100) %>%
      # Extract model coefficients
      mutate(a = summary(model.i)$coefficients[1,1], 
             a.se = summary(model.i)$coefficients[1,2], 
             b = summary(model.i)$coefficients[2,1], 
             b.se = summary(model.i)$coefficients[2,2]) %>%
      # Predict output
      mutate(fit.logit = a + b*trait, 
             fit.logit.inf = a - a.se + (b - b.se)*trait,
             fit.logit.sup = a + a.se + (b + b.se)*trait,
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.logit.inf), 
             fit.sup = plogis(fit.logit.sup), 
             p = NA_real_)
    
    # Plot predictions 
    plot.i <- data.i %>%
      mutate(disturbance = factor(disturbance, levels = c("biotic", "fire", "other", "snow", "storm"))) %>%
      ggplot(aes(x = trait, y = p)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975, color =disturbance), width = 0, alpha = 0.5) +
      geom_point(size = 1, aes(color = disturbance), alpha = 0.5) + 
      scale_color_manual(values = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) +
      geom_line(data = data.fit, aes(y = fit, group = 1), inherit.aes = TRUE) + 
      geom_line(data = data.fit, aes(y = fit.inf, group = 1), linetype = "dashed", inherit.aes = TRUE) + 
      geom_line(data = data.fit, aes(y = fit.sup, group = 1), linetype = "dashed", inherit.aes = TRUE) + 
      ylab("Disturbance sensitivity") + xlab(trait.i) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.position = "none", 
            plot.title = element_text(size = 11), 
            axis.title = element_text(size = 13)) + 
      ggtitle(paste0("Chisq = ", round(Anova(model.i)[1, 1], digits = 2), ", ",
                     scales::pvalue(Anova(model.i)[1, 3], add_p = TRUE, accuracy = 0.01))) 
    
    
    # If the model is significant, add to the output list
    if(round(Anova(model.i)[1, 3], digits = 2) <= 0.05){
      eval(parse(text = paste0("plots.out$", gsub("\\ ", "", trait.i), " <- plot.i")))
    }
  }
  
  # Create a legend for the plot
  plot.legend <- cowplot::get_legend(
    data.frame(x = c(1:5), y = c(1:5), ymin = c(0:4), ymax = c(2:6), 
               disturbance = factor(c("storm", "fire", "other", "biotic", "snow"), 
                                    levels = c("storm", "fire", "other", "biotic", "snow"))) %>%
      ggplot(aes(x = x, y = y, color = disturbance)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) + 
      scale_color_manual(values = c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77")) + 
      theme(legend.title = element_blank(), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 19)))
    
  # Final plot
  plot.out <- plot_grid(
    plot_grid(plotlist = plots.out, align = "hv", nrow = 2, scale = 0.9), 
    plot.legend, nrow = 1, rel_widths = c(1, 0.25)
  )
  
  # Save the plot
  ggsave(file.in, plot.out, width = 21, height = 14, units = "cm", dpi = 600)
  
  return(file.in)
  
}


#' Plot senstivity to all disturbances against traits on one multipanel
#' @param traits dataset containing trait values per species
#' @param traits_TRY dataset containing trait values per species from TRY
#' @param disturbance_sensitivity_full dataset containing disturbance sensitivity per species for each mcmc iteration
#' @param disturbance_sensitivity_full_bis list of dataset containing disturbance sensitivity to snow and biotic for each mcmc iteration
#' @param file.in Name of the file to save (including path)
plot_traits_vs_sensitivity_varweight_ms <- function(traits, traits_TRY, disturbance_sensitivity_full, 
                                                    disturbance_sensitivity_full_bis, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Initialize plot list
  plots.out <- list()
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- c(disturbance_sensitivity_full, disturbance_sensitivity_full_bis[c("snow", "biotic")])
  
  ## - Loop on all disturbance types to assemble data sets
  for(i in 1:length(names(disturbance_sensitivity.in))){
    data.in.i <- as.data.frame(disturbance_sensitivity.in[[i]]) %>%
      mutate(disturbance = names(disturbance_sensitivity.in)[i])
    if(i == 1) data.in <- data.in.i
    if(i > 1) data.in <- rbind(data.in, data.in.i)
  }
  
  ## - Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood density" = "wood.density_g.cm3", 
      "Shade tolerance" = "shade.tolerance", 
      "Root mass fraction" = "Root_mass_fraction", 
      "Bark thickness" = "bark.thickness_mm", 
      "H to dbh ratio" = "height.dbh.ratio", 
      "Leaf CN ratio" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf NP ratio" = "TRY_leaf.NP.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf Pmass" = "TRY_leaf.P.mass_mg.g", 
      "SLA" = "TRY_leaf.sla_mm2mg-1", 
      "Leaf thickness" = "TRY_leaf.thickness_mm", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Stomata conductance" = "TRY_stomata.conductance_millimolm-2s-1", 
      "Maximum growth" = "growth.max"
    )
  
  
  # Loop on all traits
  for(trait.i in colnames(traits.in)[which(colnames(traits.in) != "species")]){
    
    # Format data for the model
    data.i <- data.in %>%
      # Logit transformation and weight
      mutate(p.logit = log(p/(1 - p))) %>%
      group_by(species, disturbance) %>%
      summarize(w = 1/var(p.logit), 
                p = mean(p),
                p_025.logit = quantile(p.logit, probs = 0.025), 
                p_975.logit = quantile(p.logit, probs = 0.975), 
                p.logit = mean(p.logit)) %>%
      mutate(p_025 = plogis(p_025.logit), 
             p_975 = plogis(p_975.logit)) %>%
      left_join((traits.in %>% dplyr::select("species", "trait" = trait.i)), 
                by = "species") %>%
      # remove NA
      drop_na()
    
    # Scale weight so that the sum equals the number of observations
    data.i$w <- data.i$w*dim(data.i)[1]/sum(data.i$w)
    
    # Fit a model
    model.i <- lmer(p.logit ~ trait + (1|species), weights = w, data = data.i)
    
    # Data with the predictions of the model
    data.fit <- data.frame(
      trait = c(round(min(data.i$trait)*100, digits = 0):round(max(data.i$trait)*100, digits = 0))/100) %>%
      # Extract model coefficients
      mutate(a = summary(model.i)$coefficients[1,1], 
             a.se = summary(model.i)$coefficients[1,2], 
             b = summary(model.i)$coefficients[2,1], 
             b.se = summary(model.i)$coefficients[2,2]) %>%
      # Predict output
      mutate(fit.logit = a + b*trait, 
             fit.logit.inf = a - a.se + (b - b.se)*trait,
             fit.logit.sup = a + a.se + (b + b.se)*trait,
             fit = plogis(fit.logit), 
             fit.inf = plogis(fit.logit.inf), 
             fit.sup = plogis(fit.logit.sup), 
             p = NA_real_)
    
    # Plot predictions 
    plot.i <- data.i %>%
      mutate(disturbance = factor(disturbance, levels = c("biotic", "fire", "other", "snow", "storm"))) %>%
      ggplot(aes(x = trait, y = p)) + 
      geom_errorbar(aes(ymin = p_025, ymax = p_975, color =disturbance), width = 0, alpha = 0.5) +
      geom_point(size = 1, aes(color = disturbance), alpha = 0.5) + 
      scale_color_manual(values = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) +
      geom_line(data = data.fit, aes(y = fit, group = 1), inherit.aes = TRUE) + 
      geom_line(data = data.fit, aes(y = fit.inf, group = 1), linetype = "dashed", inherit.aes = TRUE) + 
      geom_line(data = data.fit, aes(y = fit.sup, group = 1), linetype = "dashed", inherit.aes = TRUE) + 
      ylab("Disturbance \n sensitivity") + xlab(trait.i) +
      theme(panel.background = element_rect(color = "black", fill = "white"), 
            panel.grid = element_blank(), 
            legend.position = "none", 
            plot.title = element_text(size = 11), 
            axis.title = element_text(size = 13)) + 
      ggtitle(paste0("Chisq = ", round(Anova(model.i)[1, 1], digits = 2), ", ",
                     scales::pvalue(Anova(model.i)[1, 3], add_p = TRUE, accuracy = 0.01))) 
    
    
    # If the model is significant, add to the output list
    if(round(Anova(model.i)[1, 3], digits = 2) <= 0.05){
      eval(parse(text = paste0("plots.out$", gsub("\\ ", "", trait.i), " <- plot.i")))
    }
  }
  
  # Create a legend for the plot
  plot.legend <- cowplot::get_legend(
    data.frame(x = c(1:5), y = c(1:5), ymin = c(0:4), ymax = c(2:6), 
               disturbance = factor(c("storm", "fire", "other", "biotic", "snow"), 
                                    levels = c("storm", "fire", "other", "biotic", "snow"))) %>%
      ggplot(aes(x = x, y = y, color = disturbance)) + 
      geom_point() + 
      geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0) + 
      scale_color_manual(values = c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77")) + 
      theme(legend.title = element_blank(), 
            legend.key = element_blank(), 
            legend.text = element_text(size = 19)))
  
  # Final plot
  plot.out <- plot_grid(
    plot_grid(plotlist = plots.out, align = "hv", nrow = 1, scale = 0.9), 
    plot.legend, nrow = 1, rel_widths = c(1, 0.25)
  )
  
  # Save the plot
  if(length(names(plots.out)) == 1) ggsave(file.in, plot.out, width = 13, height = 8, units = "cm", dpi = 600, bg = "white")
  if(length(names(plots.out)) == 2) ggsave(file.in, plot.out, width = 21, height = 8, units = "cm", dpi = 600, bg = "white")
  if(length(names(plots.out)) == 3) ggsave(file.in, plot.out, width = 28, height = 6.5, units = "cm", dpi = 600, bg = "white")
  if(length(names(plots.out)) > 3) ggsave(file.in, plot.out, width = 35, height = 6, units = "cm", dpi = 600, bg = "white")

  return(file.in)
  
}



#' Plot the relation between disturbance severity (observed) and intensity (estimated)
#' along with the distribution of each disturbance intensity
#' @param jags.model rjags object
#' @param data_jags input used for the jags model
#' @param FUNDIV_plot Plot table formatted for FUNDIV
#' @param file.in Path and file where to save the plot
map_disturbance_intensity_ms <- function(jags.model, data_jags, FUNDIV_plot, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - Identify disturbances 
  disturbances.in <- names(jags.model)
  
  # Initialize output
  plots.out <- list()
  
  ## - Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    # Create spatial object with data for disturbance i
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
    
    # Histogram plot
    hist.i <- data.i %>%
      as.data.frame() %>%
      dplyr::select(intensity) %>%
      mutate(cat = cut(intensity, breaks = c(0:10)/10)) %>%
      mutate(cat = ifelse(intensity < 0.1, 0, as.numeric(substr(cat, 2, 4)))) %>%
      group_by(cat) %>%
      summarise(n = n()) %>%
      ggplot(aes(x = cat, y = n, fill = cat)) +
      geom_bar(colour = "black", stat = "identity")  + 
      scale_fill_gradientn(colours = colorRampPalette(c("#EDF2F4", "#EF233C", "#D90429"))(10)) + 
      scale_x_continuous(breaks = c(0:9)/10) +
      xlab("") + ylab("") +
      theme(panel.background = element_rect(color = 'black', fill = 'white'), 
            panel.grid = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank(), 
            legend.position = "none", 
            plot.margin = unit(c(0, 0, 0, 0), "cm")) 
    
    # Plot for disturbance i with histogram inserted
    plot.i <- ne_countries(scale = "medium", returnclass = "sf") %>%
      mutate(disturbance = disturbances.in[i]) %>%
      mutate(keep = case_when(
        disturbance %in% c("fire", "storm", "other") ~ ifelse(sovereignt %in% c("France", "Spain", "Finland"), "yes", "no"), 
        disturbance %in% c("biotic", "snow") ~ ifelse(sovereignt %in% c("Spain", "Finland"), "yes", "no"))) %>%
      ggplot(aes(geometry = geometry)) +
      geom_sf(aes(fill = keep), color = "white", show.legend = F, size = 0.2) + 
      scale_fill_manual(values = c("#8D99AE", "#343A40")) +
      geom_sf(data = data.i, shape = 16, aes(color = intensity), 
              show.legend = "point", size = 0.4) +
      scale_color_gradient2(low = "#EDF2F4", mid = "#EF233C", high = "#D90429", midpoint = 0.4)  +
      guides(fill = FALSE) +
      theme(panel.background = element_rect(color = 'black', fill = 'white'), 
            panel.grid = element_blank(), 
            plot.title = element_text(size = 15, hjust = 0.5), 
            legend.position = "none", 
            axis.text = element_blank(), 
            axis.ticks = element_blank()) + 
      coord_sf(xlim = c(-10, 32), ylim = c(36, 71)) + 
      ggtitle(paste0(toupper(substr(disturbances.in[i], 1, 1)), 
                     substr(disturbances.in[i], 2, nchar(disturbances.in[i])))) +
      annotation_custom(ggplotGrob(hist.i), xmin = -12, xmax = 3, ymin = 62, ymax = 72)
    
    # Add to the output plot list
    eval(parse(text = paste0("plots.out$", disturbances.in[i], " <- plot.i")))
    
  }
  
  ## - Add a legend at the end of the list
  plots.out$legend <- cowplot::get_legend(plot.i + theme(legend.position = "left"))
  
  ## - Final plot
  plot.out <- plot_grid(plotlist = plots.out[c("storm", "fire", "other", "biotic", "snow", "legend")], 
                        nrow = 2, labels = c("(a)", "", "", "(b)", "", "")) 
  
  ## - save the plot
  ggsave(file.in, plot.out, width = 20, height = 18, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}


#' Plot trends in disturbance frequency
#' @param FUNDIV_plot plot table formatted for FUNDIV with biotic and snow in other category
#' @param FUNDIV_plot_bis plot table formatted for FUNDIV with biotic and snow only
#' @param file.in Name of the file too save (including path)
plot_trend_disturbance_frequency_ms <- function(FUNDIV_plot, FUNDIV_plot_bis, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Merge FUNDIV plot table to have biotic and snow as separate disturbance categories
  FUNDIV_plot.in <- FUNDIV_plot %>%
    filter(!(plotcode %in% FUNDIV_plot_bis$plotcode)) %>%
    rbind.data.frame(FUNDIV_plot_bis)
  
  ## - Arrange the data for plotting
  data <- FUNDIV_plot.in %>%
    # Removee unknown disturbances
    filter(!is.na(disturbance.nature)) %>%
    # Only keep france and spain (not enough temporal variation in Finland)
    filter(country %in% c("France", "Spain")) %>%
    # Calculate the median year for the sampling of each plot
    mutate(medianyear = floor((surveydate1 + surveydate2)/2)) %>%
    # Calculate the number of ocurrence of each disturbance per year and country
    group_by(country, medianyear, disturbance.nature) %>%
    summarize(n.plots = n()) %>%
    # Convert into a proportion and remove cases with no disturbances
    ungroup() %>% group_by(country, medianyear) %>%
    mutate(proportion.plots = n.plots/sum(n.plots, na.rm = TRUE), 
           total.plots.per.year = sum(n.plots, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(disturbance.nature != "none") %>%
    dplyr::select(-n.plots) %>%
    filter(total.plots.per.year > 1000) %>%
    # Convert disturbance in factor to control the order in the plot
    mutate(disturbance.nature = factor(
      disturbance.nature, levels = c("biotic", "fire", "other", "snow", "storm"))) %>%
    # Convert in an annual proportion by dividing by the number of year between surveys
    mutate(proportion.plots = ifelse(country == "France", proportion.plots/5, proportion.plots/10))
  
  ## - Data to add the number of plots above each bar in the plot
  data.label <- data %>%
    ungroup() %>% group_by(country, medianyear) %>%
    summarize(n = mean(total.plots.per.year), 
              proportion.plots = sum(proportion.plots) + 0.0003) %>%
    mutate(disturbance.nature = NA, label = paste0("(", n, ")"))
  
  ## - Plot proportion of disturbed plots
  plot.out <- data %>%
    ggplot(aes(x = as.factor(medianyear), y = proportion.plots*100, fill = disturbance.nature))  +
    geom_bar(stat="identity", color = "black", alpha = 0.65) +
    scale_fill_manual(values = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) +
    facet_grid(. ~ country, scales = "free_x", space = "free") + 
    xlab("Median year between census 1 and 2") + ylab("Annual proportion of plots disturbed (%)") +
    geom_text(data = data.label, aes(label = label), inherit.aes = TRUE, size = 2.5) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text = element_text(size = 18),
          axis.title = element_text(size = 15),
          legend.title = element_blank(), 
          legend.text = element_text(size = 15),
          panel.spacing = unit(2, "lines"), 
          axis.text.x = element_text(angle = 360))
  
  ## - Save the plot
  ggsave(file.in, plot.out, width = 19, height = 12, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}



#' Plot trends in disturbance severity
#' @param FUNDIV_tree tree table formatted for FUNDIV
#' @param FUNDIV_plot plot table formatted for FUNDIV with biotic and snow in other category
#' @param FUNDIV_plot_bis plot table formatted for FUNDIV with biotic and snow only
#' @param file.in Name of the file too save (including path)
plot_trend_disturbance_severity_ms <- function(FUNDIV_tree, FUNDIV_plot, 
                                               FUNDIV_plot_bis, file.in){
  
  ## - Create directory if needed
  create_dir_if_needed(file.in)
  
  ## - Merge FUNDIV plot table to have biotic and snow as separate disturbance categories
  FUNDIV_plot.in <- FUNDIV_plot %>%
    filter(!(plotcode %in% FUNDIV_plot_bis$plotcode)) %>%
    rbind.data.frame(FUNDIV_plot_bis)
  
  ## - Arrange data to study changes in disturbance severity
  data.severity <- FUNDIV_plot.in %>%
    # Only keep disturbed plots
    filter(disturbance.nature %in% c("biotic", "fire", "other", "snow", "storm")) %>%
    # Only keep france and spain (not enough temporal variation in Finland)
    filter(country %in% c("France", "Spain")) %>%
    # Calculate the median year for the sampling of each plot
    mutate(medianyear = floor((surveydate1 + surveydate2)/2)) %>%
    # Join severity per plot
    left_join((FUNDIV_tree %>%
                 filter(treestatus != 1) %>%
                 filter(dbh1 > 100) %>%
                 mutate(dead = ifelse(treestatus == 2, 0, 1)) %>%
                 group_by(plotcode) %>%
                 summarize(severity = sum(dead)/n())), 
              by = "plotcode") %>%
    # Remove null severities
    filter(severity > 0) %>%
    # Calculate annual severity by dividing by the time interval between two studies
    mutate(severity.annual = severity/yearsbetweensurveys) %>%
    # Calculate mean and confidence interval per year, disturbance and country
    ungroup() %>%
    group_by(country, medianyear, disturbance.nature) %>%
    summarize(severity.mean = mean(severity.annual), 
              severity.low = quantile(severity.annual, probs = 0.025), 
              severity.high = quantile(severity.annual, probs = 0.975),
              n = n()) %>%
    # Add label and label position
    mutate(label = paste0("(", n, ")")) %>%
    # Convert disturbance in factor to control the order in the plot
    mutate(disturbance.nature = factor(
      disturbance.nature, levels = c("biotic", "fire", "other", "snow", "storm")))
  
  ## - Plot trends in disturbance severity
  plot.out <- data.severity %>%
    ggplot(aes(x = as.factor(medianyear), y = severity.mean, color = disturbance.nature))  +
    geom_point() +
    geom_errorbar(aes(ymin = severity.low, ymax = severity.high), width = 0) + 
    scale_color_manual(values = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")) +
    facet_grid(disturbance.nature ~ country, scales = "free_x", space = "free") + 
    xlab("Median year between census 1 and 2") + ylab("Mean severity") +
    geom_text(aes(label = label, y = severity.high), size = 2.5, nudge_y = 0.05) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.y = element_blank(),
          strip.text = element_text(size = 18),
          axis.title = element_text(size = 15),
          legend.title = element_blank(), 
          legend.key = element_blank(),
          legend.text = element_text(size = 15),
          panel.spacing = unit(1, "lines"), 
          axis.text.x = element_text(angle = 90))
  
  ## - Save the plot
  ggsave(file.in, plot.out, width = 16, height = 16, units = "cm", dpi = 600, bg = "white")
  return(file.in)
}


#' Plot the effect of traits on disturbance sensitivity
#' @param traits dataframe containing trait values per species
#' @param traits_TRY dataframe containing trait values from TRY per species
#' @param disturbance_sensivity dataframe containing the sensitivity to each disturbance
#' @param disturbance_sensivity_bis dataframe containing the sensitivity to biotic and snow
#' @param species Table containing species information
#' @param group.in character indicating which species to include ("all", "conifer" or "broadleaf")
#' @param file.in Name of the file to save
plot_trait_effect_ms <- function(traits, traits_TRY, disturbance_sensitivity, disturbance_sensitivity_bis, 
                                 species, group.in = "all", file.in){
  
  # Create dir if needed
  create_dir_if_needed(file.in)
  
  # merge disturbance sensitivity 
  disturbance_sensitivity.in <- c(disturbance_sensitivity, disturbance_sensitivity_bis[c("biotic", "snow")])
  
  # Identify the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood dens." = "wood.density_g.cm3", 
      "Shade tol." = "shade.tolerance", 
      "Root mass frac." = "Root_mass_fraction", 
      "Bark thick." = "bark.thickness_mm", 
      "H/dbh ratio" = "height.dbh.ratio", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Max. growth" = "growth.max", 
      "Leaf C/N" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf thick." = "TRY_leaf.thickness_mm", 
      "Stomata cond." = "TRY_stomata.conductance_millimolm-2s-1"
    )
  
  # Identify the species to select depending on the group 
  species.to.select <- (species %>%
                          mutate(group.in = group.in) %>%
                          mutate(keep = case_when(group.in == "conifer" ~ ifelse(group == "Gymnosperms", 1, 0), 
                                                  group.in == "broadleaf" ~ ifelse(group == "Angiosperms", 1, 0), 
                                                  group.in == "all" ~ 1)) %>%
                          filter(keep == 1))$species
  
  # And filter trait table
  traits.in <- traits.in %>% filter(species %in% species.to.select)
  
  # Center and scale the trait values
  traits.in <- scale_data_model(traits.in, var = colnames(traits.in)[c(2:dim(traits.in)[2])])
  
  # Loop on all traits
  for(i in 1:(dim(traits.in)[2] - 1)){
    # Identify the name of trait i
    trait.i <- colnames(traits.in)[i+1]
    # Create a table with only species and trait i
    traits.i <- traits.in %>% dplyr::select("species", "trait" = trait.i)
    # Loop on all type of disturbances
    for(j in 1:length(disturbances.in)){
      # Create a table with trait i and sensitivity to disturbance j
      data.ij <- traits.i %>%
        left_join((disturbance_sensitivity.in[[j]] %>%
                     mutate(sensitivity.logit = log(p/(1 - p)), 
                            w = 1/(p_975 - p_025))), 
                  by = "species") %>%
        drop_na()
      # Only perform a test if there is enough data
      if(dim(data.ij)[1] > 3){
        # Fit a model
        model.ij <- lm(sensitivity.logit ~ trait, weights = w, data = data.ij)
        # Results
        table.ij <- data.frame(
          trait = trait.i, 
          disturbance = disturbances.in[j],
          n = dim(data.ij)[1],
          Est = summary(model.ij)$coefficients[2, 1], 
          Est.se = summary(model.ij)$coefficients[2, 2], 
          Est.sup = confint.lm(model.ij)[2, 2], 
          Est.inf = confint.lm(model.ij)[2, 1],
          p = anova(model.ij)[1, 5]
        )
      }else{table.ij <- data.frame(trait = trait.i, disturbance = disturbances.in[j],
                                   n = NA_real_, Est = NA_real_, Est.se = NA_real_, 
                                   Est.sup = NA_real_, Est.inf = NA_real_, p = 1)}
      
      
      # Add to the list containing the final results
      if(i == 1 & j == 1) data <- table.ij
      else data <- rbind.data.frame(data, table.ij)
    }
  }
  
  # Add categories for each trait
  data <- data %>%
    mutate(trait.category = case_when(
      trait %in% c("Wood dens.", "Lifespan", "Max. growth") ~ "Growth vs.\n  survival", 
      trait %in% c("Leaf thick.", "Stomata cond.") ~ "Drought \n traits", 
      trait %in% c("Leaf C/N", "Leaf Nmass") ~ "Growth vs.\n defense", 
      TRUE ~ "Other\ntraits"
    ), 
    significance = ifelse(p <= 0.05, "*", ""), 
    label = ifelse(n > 3, paste0("(", n, ") ", significance), ""))
  
  ## - Make the plot
  plot.out <- data %>%
    mutate(disturbance = factor(disturbance, levels = c("storm", "fire", "other", "biotic",  "snow"))) %>%
    ggplot(aes(x = trait, y = Est, color = disturbance)) + 
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.3) + 
    geom_point(size = 1) + 
    geom_errorbar(aes(ymin = Est.inf, ymax = Est.sup), width = 0) + 
    facet_grid(trait.category ~ disturbance, scales = "free_y", space = "free_y") +
    geom_text(aes(label = label, y = max(data$Est.sup, na.rm = TRUE)), 
              size = 2.5, nudge_y = 3, hjust = "inward", show.legend = F) +
    scale_color_manual(values = c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77")) +
    xlab("") + ylab("Trait effect on disturbance sensitivity") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 13),
          strip.text.y = element_text(size = 13, angle = 360),
          axis.title = element_text(size = 15),
          legend.position = "none",
          axis.text.x = element_text(size = 10, angle = 360),
          axis.text.y = element_text(size = 10, angle = 360)) + 
    coord_flip() + 
    ylim(min(data$Est.inf, na.rm = T), (max(data$Est.sup, na.rm = T) + 3))
  
  ## - Save the plot
  ggsave(file.in, plot.out, width = 23, height = 9, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}



#' Plot the effect of traits on disturbance sensitivity using betareg and all mcmc iteration as input
#' @param traits dataframe containing trait values per species
#' @param traits_TRY dataframe containing trait values from TRY per species
#' @param disturbance_sensivity_full dataframe containing the sensitivity to each disturbance
#' @param disturbance_sensivity_full_bis dataframe containing the sensitivity to biotic and snow
#' @param species Table containing species information
#' @param group.in character indicating which species to include ("all", "conifer" or "broadleaf")
#' @param file.in Name of the file to save
plot_trait_effect_full_ms <- function(traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
                                      species, group.in = "all", file.in){
  
  # Create dir if needed
  create_dir_if_needed(file.in)
  
  # merge disturbance sensitivity 
  disturbance_sensitivity.in <- c(disturbance_sensitivity_full, disturbance_sensitivity_full_bis[c("biotic", "snow")])
  
  # Identify the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood dens." = "wood.density_g.cm3", 
      "Shade tol." = "shade.tolerance", 
      "Root mass frac." = "Root_mass_fraction", 
      "Bark thick." = "bark.thickness_mm", 
      "H/dbh ratio" = "height.dbh.ratio", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Max. growth" = "growth.max", 
      "Leaf C/N" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf thick." = "TRY_leaf.thickness_mm", 
      "Stomata cond." = "TRY_stomata.conductance_millimolm-2s-1"
    )
  
  # Identify the species to select depending on the group 
  species.to.select <- (species %>%
                          mutate(group.in = group.in) %>%
                          mutate(keep = case_when(group.in == "conifer" ~ ifelse(group == "Gymnosperms", 1, 0), 
                                                  group.in == "broadleaf" ~ ifelse(group == "Angiosperms", 1, 0), 
                                                  group.in == "all" ~ 1)) %>%
                          filter(keep == 1))$species
  
  # And filter trait table
  traits.in <- traits.in %>% filter(species %in% species.to.select)
  
  # Center and scale the trait values
  traits.in <- scale_data_model(traits.in, var = colnames(traits.in)[c(2:dim(traits.in)[2])])
  
  # Loop on all traits
  for(i in 1:(dim(traits.in)[2] - 1)){
    # Identify the name of trait i
    trait.i <- colnames(traits.in)[i+1]
    # Create a table with only species and trait i
    traits.i <- traits.in %>% dplyr::select("species", "trait" = trait.i)
    # Loop on all type of disturbances
    for(j in 1:length(disturbances.in)){
      
      # Create a table with trait i and sensitivity to disturbance j
      data.ij <- disturbance_sensitivity.in[[j]] %>%
        mutate(sensitivity.logit = log(p/(1 - p)), 
               w = 1/length(unique(.$iter))) %>%
        left_join((traits.i), 
                  by = "species") %>%
        drop_na()
      
      # Only perform a test if there is enough data
      if(length(unique(data.ij$species)) > 3){
        
        # Fit a model depending on model type chosen
        model.ij <- betareg(p ~ trait, weights = w, data = data.ij, link = "logit")
        
        # Extract results
        table.ij <- data.frame(
          trait = trait.i, 
          disturbance = disturbances.in[j],
          n = length(unique(data.ij$species)),
          Est.sup = confint(model.ij)[2, 2], 
          Est.inf = confint(model.ij)[2, 1], 
          Est = coefficients(summary(model.ij))$mean[2, 1]
        )
        
      }else{table.ij <- data.frame(trait = trait.i, disturbance = disturbances.in[j], n = NA_real_, 
                                   Est.sup = NA_real_, Est.inf = NA_real_, Est = NA_real_)}
      
      
      # Add to the list containing the final results
      if(i == 1 & j == 1) data <- table.ij
      else data <- rbind.data.frame(data, table.ij)
    }
  }
  
  # Add categories for each trait
  data <- data %>%
    mutate(trait.category = case_when(
      trait %in% c("Wood dens.", "Lifespan", "Max. growth") ~ "Growth vs.\n  survival", 
      trait %in% c("Leaf thick.", "Stomata cond.") ~ "Drought \n traits", 
      trait %in% c("Leaf C/N", "Leaf Nmass") ~ "Growth vs.\n defense", 
      TRUE ~ "Other\ntraits"
    ), 
    significance = ifelse((Est.inf > 0 | Est.sup < 0), "*", ""), 
    label = ifelse(n > 3, paste0("(", n, ") ", significance), ""))
  
  ## - Make the plot
  plot.out <- data %>%
    mutate(disturbance = factor(disturbance, levels = c("storm", "fire", "other", "biotic",  "snow"))) %>%
    ggplot(aes(x = trait, y = Est, color = disturbance)) + 
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.3) + 
    geom_point(size = 1) + 
    geom_errorbar(aes(ymin = Est.inf, ymax = Est.sup), width = 0) + 
    facet_grid(trait.category ~ disturbance, scales = "free_y", space = "free_y") +
    geom_text(aes(label = label, y = max(data$Est.sup, na.rm = TRUE)), 
              size = 2.5, nudge_y = 3, hjust = "inward", show.legend = F) +
    scale_color_manual(values = c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77")) +
    xlab("") + ylab("Trait effect on disturbance sensitivity") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 13),
          strip.text.y = element_text(size = 13, angle = 360),
          axis.title = element_text(size = 15),
          legend.position = "none",
          axis.text.x = element_text(size = 10, angle = 360),
          axis.text.y = element_text(size = 10, angle = 360)) + 
    coord_flip() + 
    ylim(min(data$Est.inf, na.rm = T), (max(data$Est.sup, na.rm = T) + 3))
  
  ## - Save the plot
  ggsave(file.in, plot.out, width = 23, height = 9, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}


#' Plot the effect of traits on disturbance sensitivity using betareg and inverse variance as weight
#' @param traits dataframe containing trait values per species
#' @param traits_TRY dataframe containing trait values from TRY per species
#' @param disturbance_sensivity_full dataframe containing the sensitivity to each disturbance
#' @param disturbance_sensivity_full_bis dataframe containing the sensitivity to biotic and snow
#' @param species Table containing species information
#' @param group.in character indicating which species to include ("all", "conifer" or "broadleaf")
#' @param file.in Name of the file to save
plot_trait_effect_varweight_ms <- function(traits, traits_TRY, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, 
                                           species, group.in = "all", file.in){
  
  # Create dir if needed
  create_dir_if_needed(file.in)
  
  # merge disturbance sensitivity 
  disturbance_sensitivity.in <- c(disturbance_sensitivity_full, disturbance_sensitivity_full_bis[c("biotic", "snow")])
  
  # Identify the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Rearrange traits table
  traits.in <- traits %>%
    left_join(traits_TRY, by = "species") %>%
    dplyr::select(
      "species", 
      "Wood dens." = "wood.density_g.cm3", 
      "Shade tol." = "shade.tolerance", 
      "Root mass frac." = "Root_mass_fraction", 
      "Bark thick." = "bark.thickness_mm", 
      "H/dbh ratio" = "height.dbh.ratio", 
      "Lifespan" = "TRY_plant.lifespan_year", 
      "Max. growth" = "growth.max", 
      "Leaf C/N" = "TRY_leaf.CN.ratio_g.cm3", 
      "Leaf Nmass" = "TRY_leaf.N.mass_mg.g", 
      "Leaf thick." = "TRY_leaf.thickness_mm", 
      "Stomata cond." = "TRY_stomata.conductance_millimolm-2s-1"
    )
  
  # Identify the species to select depending on the group 
  species.to.select <- (species %>%
                          mutate(group.in = group.in) %>%
                          mutate(keep = case_when(group.in == "conifer" ~ ifelse(group == "Gymnosperms", 1, 0), 
                                                  group.in == "broadleaf" ~ ifelse(group == "Angiosperms", 1, 0), 
                                                  group.in == "all" ~ 1)) %>%
                          filter(keep == 1))$species
  
  # And filter trait table
  traits.in <- traits.in %>% filter(species %in% species.to.select)
  
  # Center and scale the trait values
  traits.in <- scale_data_model(traits.in, var = colnames(traits.in)[c(2:dim(traits.in)[2])])
  
  # Loop on all traits
  for(i in 1:(dim(traits.in)[2] - 1)){
    # Identify the name of trait i
    trait.i <- colnames(traits.in)[i+1]
    # Create a table with only species and trait i
    traits.i <- traits.in %>% dplyr::select("species", "trait" = trait.i)
    # Loop on all type of disturbances
    for(j in 1:length(disturbances.in)){
      
      # Create a table with trait i and sensitivity to disturbance j
      data.ij <- disturbance_sensitivity.in[[j]] %>%
        mutate(p.logit = log(p/(1 - p))) %>%
        group_by(species) %>%
        summarize(w = 1/var(p.logit), 
                  p.logit = mean(p.logit)) %>%
        mutate(p = plogis(p.logit)) %>%
        left_join((traits.i), 
                  by = "species") %>%
        drop_na()
      
      # Only perform a test if there is enough data
      if(length(unique(data.ij$species)) > 3){
        
        # Fit a model depending on model type chosen
        model.ij <- lm(p.logit ~ trait, weights = w, data = data.ij)
        
        # Extract results
        table.ij <- data.frame(
          trait = trait.i, 
          disturbance = disturbances.in[j],
          n = dim(data.ij)[1],
          Est.sup = confint.lm(model.ij)[2, 2], 
          Est.inf = confint.lm(model.ij)[2, 1],
          Est = summary(model.ij)$coefficients[2, 1]
        )
        
      }else{table.ij <- data.frame(trait = trait.i, disturbance = disturbances.in[j], n = NA_real_, 
                                   Est.sup = NA_real_, Est.inf = NA_real_, Est = NA_real_)}
      
      
      # Add to the list containing the final results
      if(i == 1 & j == 1) data <- table.ij
      else data <- rbind.data.frame(data, table.ij)
    }
  }
  
  # Add categories for each trait
  data <- data %>%
    mutate(trait.category = case_when(
      trait %in% c("Wood dens.", "Lifespan", "Max. growth") ~ "Growth vs.\n  survival", 
      trait %in% c("Leaf thick.", "Stomata cond.") ~ "Drought \n traits", 
      trait %in% c("Leaf C/N", "Leaf Nmass") ~ "Growth vs.\n defense", 
      TRUE ~ "Other\ntraits"
    ), 
    significance = ifelse((Est.inf > 0 | Est.sup < 0), "*", ""), 
    label = ifelse(n > 3, paste0("(", n, ") ", significance), ""))
  
  ## - Make the plot
  plot.out <- data %>%
    mutate(disturbance = factor(disturbance, levels = c("storm", "fire", "other", "biotic",  "snow"))) %>%
    ggplot(aes(x = trait, y = Est, color = disturbance)) + 
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 0.3) + 
    geom_point(size = 1) + 
    geom_errorbar(aes(ymin = Est.inf, ymax = Est.sup), width = 0) + 
    facet_grid(trait.category ~ disturbance, scales = "free_y", space = "free_y") +
    geom_text(aes(label = label, y = max(data$Est.sup, na.rm = TRUE)), 
              size = 2.5, nudge_y = 3, hjust = "inward", show.legend = F) +
    scale_color_manual(values = c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77")) +
    xlab("") + ylab("Trait effect on disturbance sensitivity") +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(), 
          strip.text.x = element_text(size = 13),
          strip.text.y = element_text(size = 13, angle = 360),
          axis.title = element_text(size = 15),
          legend.position = "none",
          axis.text.x = element_text(size = 10, angle = 360),
          axis.text.y = element_text(size = 10, angle = 360)) + 
    coord_flip() + 
    ylim(min(data$Est.inf, na.rm = T), (max(data$Est.sup, na.rm = T) + 3))
  
  ## - Save the plot
  ggsave(file.in, plot.out, width = 23, height = 9, units = "cm", dpi = 600, bg = "white")
  return(file.in)
  
}




#' Function to plot results of the climate analysis for the manuscript
#' @param gbif_file Name of the file containing climate data per species
#' @param disturbance_sensitivity_full list of dataset containing disturbance sensitivity per species
#' @param disturbance_sensitivity_full_bis list of dataset containing disturbance sensitivity to snow and biotic
#' @param file.in Where to save the plot
plot_climate_effect_ms <- function(gbif_file, disturbance_sensitivity_full, disturbance_sensitivity_full_bis, file.in){
  
  # Create directory if needed
  create_dir_if_needed(file.in)
  
  
  # - Make PCA 
  pca <- prcomp((fread(gbif_file) %>%
                   dplyr::select(species, mat, tmin, map) %>%
                   drop_na() %>% 
                   dplyr::select(-species)), 
                center = T, scale = T)
  # - Extract the coordinates of the individuals on pca axis
  data.climate <- data.frame(species = (fread(gbif_file) %>%
                                          dplyr::select(species, mat, tmin, map) %>%
                                          drop_na())$species, 
                             pca1 = get_pca_ind(pca)[[1]][, 1], 
                             pca2 = get_pca_ind(pca)[[1]][, 2]) 
  
  # - Extract the coordinates of the individuals on pca axis
  res.ind <- data.frame(species = (fread(gbif_file) %>%
                                     dplyr::select(species, mat, tmin, map) %>%
                                     drop_na())$species, 
                        pca1 = get_pca_ind(pca)[[1]][, 1], 
                        pca2 = get_pca_ind(pca)[[1]][, 2]) %>%
    mutate(sp = paste(substr(gsub("\\ .+", "", species), 1, 2), 
                      substr(gsub(".+\\ ", "", species), 1, 2), 
                      sep = "."))
  
  # - Extract the coordinates of the variables on pca axis
  res.var <- data.frame(var = rownames(get_pca_var(pca)[[1]]), 
                        pca1 = get_pca_var(pca)[[1]][, 1], 
                        pca2 = get_pca_var(pca)[[1]][, 2])
  
  # - Minimum and maximum in each pca axis
  pca.xmin <- -max(abs(res.ind$pca1))
  pca.xmax <- max(abs(res.ind$pca1))
  pca.ymin <- -max(abs(res.ind$pca2))
  pca.ymax <- max(abs(res.ind$pca2))
  
  # Make the plot
  plot.pca <- res.ind %>%
    ggplot(aes(x = pca1, y = pca2)) + 
    geom_point(fill = "#023E8A", color = "black", shape = 21) +
    geom_text(aes(label = sp), nudge_y = 0.1, color = "#023E8A", size = 3) +
    geom_segment(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
                 aes(x = 0, xend = pca1, y = 0, yend = pca2), 
                 arrow = arrow(length = unit(0.1, "cm")), 
                 type = "closed", color = "#D90429") + 
    geom_text(data = (res.var %>% mutate(pca1 = pca1*1.5, pca2 = pca2*1.5)), 
              aes(label = var), color = "#D90429", size = 5, 
              nudge_x = ifelse(res.var$pca1 < 0, pca.xmin/12, pca.xmax/12)) +
    geom_hline(size = 0.2, yintercept = 0, color = "#6C757D", linetype = "dashed") + 
    geom_vline(size = 0.2, xintercept = 0, color = "#6C757D", linetype = "dashed") + 
    xlim((pca.xmin-0.2), (pca.xmax+0.2)) + 
    ylim((pca.ymin-0.2), (pca.ymax+0.2)) +
    xlab(paste0("PCA1 (", round(summary(pca)$importance[2, 1]*100, digits = 2), "%)")) +
    ylab(paste0("PCA2 (", round(summary(pca)$importance[2, 2]*100, digits = 2), "%)")) +
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 15))
  
  
  ## - Assemble the two disturbance sensitivity files
  disturbance_sensitivity.in <- c(disturbance_sensitivity_full, disturbance_sensitivity_full_bis[c("snow", "biotic")])
  
  ## - Names of the disturbances
  disturbances.in <- names(disturbance_sensitivity.in)
  
  # Create a vector of colors for plotting
  color.in <- (data.frame(disturbance = disturbances.in) %>%
                 left_join(data.frame(disturbance = c("biotic", "fire", "other", "snow", "storm"), 
                                      color = c("#90A955", "#F77F00", "#5F0F40", "#006D77", "#4361EE")), 
                           by = "disturbance"))$color
  
  # Compile trait dataset
  data.in <- data.climate %>% rename("var.1" = colnames(.)[2], "var.2" = colnames(.)[3])
  
  
  ## - Initialize output
  plots.out <- list()
  
  ## - Loop on all type of disturbances
  for(i in 1:length(disturbances.in)){
    
    # Data to fit the model for disturbance i
    data.i <- disturbance_sensitivity.in[[i]]  %>%
      mutate(p.logit = log(p/(1 - p))) %>%
      group_by(species) %>%
      summarize(w = 1/var(p.logit), 
                p.logit = mean(p.logit)) %>%
      left_join(data.in, by = "species") %>%
      drop_na()
    
    # Fit model depending on weighting method
    model.i = lm(p.logit ~ scale(var.1) + scale(var.2), weights = w, data = data.i)
    
    # Build result table to plot results
    results.i <- data.frame(disturbance = disturbances.in[i],
                            n = length(unique(data.i$species)),
                            var = colnames(data.climate)[c(2, 3)])
    if(dim(data.i)[1] > 3) results.i <- results.i %>% 
      mutate(est.low = as.numeric(confint(model.i)[c(2, 3), 1]), 
             est.high = as.numeric(confint(model.i)[c(2, 3), 2]))
    else results.i <- results.i %>% mutate(est.low = NA_real_, est.high = NA_real_)
    
    
    # Different extraction of estimate depending on model type
    results.i <- results.i %>% 
      mutate(est = coefficients(summary(model.i))[c(2, 3), 1]) %>% 
      mutate(est = ifelse(is.na(est.low), NA_real_, est))
    
    # Add to the final table
    if(i == 1) results <- results.i
    else results <- rbind(results, results.i)
    
  }
  
  # Label for the plot
  breaks.label <- (data.frame(disturbance = c("storm", "fire", "other", "biotic", "snow")) %>%
                     left_join(results %>%
                                 mutate(label = paste0(disturbance, "\n (n = ", n, ")")) %>%
                                 dplyr::select(disturbance, label) %>%
                                 distinct(), 
                               by = "disturbance"))$label
  
  # Final plot
  plot.effect <- results %>%
    mutate(disturbance = factor(disturbance, levels = rev(c("storm", "fire", "other", "biotic", "snow"))), 
           var = toupper(var)) %>%
    ggplot(aes(x = disturbance, y = est, color = disturbance)) + 
    geom_point() + 
    geom_errorbar(aes(ymin = est.low, ymax = est.high), width = 0) + 
    scale_color_manual(values = rev(c("#4361EE", "#F77F00", "#5F0F40", "#90A955", "#006D77"))) +
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ylab("Effect on disturbance sensitivity") + 
    xlab("") + 
    scale_x_discrete(labels = rev(breaks.label)) +
    facet_wrap( ~ var, nrow = 1) + 
    coord_flip() + 
    theme(panel.background = element_rect(color = "black", fill = "white"), 
          panel.grid = element_blank(), 
          strip.background = element_blank(),
          legend.position = "none", 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 14), 
          axis.title = element_text(size = 15), 
          strip.text = element_text(size = 15))
  
  
  # Save the plot
  plot.out <- plot_grid(plot.pca, plot.effect, nrow = 1, scale = 0.9, 
                        labels = c("(a)", "(b)"), rel_widths = c(1, 1.3), align = "h")
  ggsave(file.in, plot.out, width = 25, height = 10, units = "cm", dpi = 600, bg = "white")
  
  # Return file name
  return(file.in)
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




