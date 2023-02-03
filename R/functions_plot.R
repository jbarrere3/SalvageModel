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

