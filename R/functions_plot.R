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


