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
  disturbances.in <- c(unique(subset(FUNDIV_plot, disturbance.nature != "none")$disturbance.nature), "all")
  
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
      if(sum(data.ij$disturbed, na.rm = TRUE) > 0){
        
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
  ggsave(file.in, plot.out, width = 25, height = 20, units = "cm", dpi = 600)
  
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
  ggsave(file.in, plot.out, width = 18, height = 30, units = "cm", dpi = 600)
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
                                         dbh.ref = 300, stock.ref = 20, I.ref = 0.5, 
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


#' Function to plot results of the climate analysis 
#' @param gbif_file Name of the file containing species mean climate
#' @param gbif_disturbance_file Name of the file containing species mean disturbance index
#' @param disturbance_sensitivity list of dataset containing disturbance sensitivity per species
#' @param file.in Where to save the plot
plot_rda_climate <- function(disturbance_sensitivity, gbif_file, gbif_disturbance_file, file.in){
  
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
    
    # Add variable depending on the disturbance
    if(disturbances.in[i] == "storm") traits.i <- left_join(traits.i, (fread(gbif_disturbance_file) %>% dplyr::select(species, windspeed)), by = "species")
    if(disturbances.in[i] == "snow") traits.i <- left_join(traits.i, (fread(gbif_disturbance_file) %>% dplyr::select(species, swe)), by = "species")
    if(disturbances.in[i] == "fire") traits.i <- left_join(traits.i, (fread(gbif_disturbance_file) %>% dplyr::select(species, fwi)), by = "species")
    
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
    plot.rda <- data.frame(trait = c("windspeed", "fwi", "swe", "map", "mat", "tmin"), 
                           position = c(1:6)) %>%
      left_join(data.frame(trait = as.character(rownames(summary(rda)$species)), 
                           rda1 = summary(rda)$species[, 1]), by = "trait") %>%
      mutate(trait = factor(trait, levels = c("map", "tmin", "mat", "windspeed", "fwi", "swe"))) %>%
      ggplot(aes(x = 0, xend = rda1, y = position, yend = position)) + 
      geom_segment(arrow = arrow(length = unit(0.1, "cm")), type = "closed") + 
      scale_y_continuous(breaks = c(1:6), 
                         labels = c("windspeed", "fwi", "swe", "map", "mat", "tmin"), 
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
    ggsave(file.in, plot.out, width = 35, height = 13, units = "cm", dpi = 600)
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



