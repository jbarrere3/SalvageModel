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

# Add comment to test branching. 

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


#' Plot Markov chain convergence of a rjags object
#' @param fit.rjags.in rjags object
#' @param file.in Path and file where to save the plot
#' @param title.in plot title
plot_convergence <- function(fit.rjags.in, file.in, title.in = ""){
  
  # Create directories if needed
  create_dir_if_needed(file.in)
  
  # make the plot
  plot.in <- ggs(as.mcmc(fit.rjags.in)) %>%
    filter(Parameter != "deviance") %>%
    mutate(Chain = as.factor(Chain)) %>%
    ggplot(aes(x = Iteration, y = value, colour = Chain, group = Chain)) + 
    geom_line() + 
    facet_wrap(~ Parameter, scales = "free") + 
    scale_color_manual(values = c("#335C67", "#E09F3E", "#9E2A2B")) +
    theme_bw() + 
    ggtitle(title.in) 
  
  # save the plot
  ggsave(file.in, plot.in, width = 17, height = 12, units = "cm", dpi = 600)
  return(file.in)
}


#' Function to compare the parameters true value and estimation by the bayesian model
#' @param data_jags simulated data fot the bayesian model
#' @param jags.model rjags object: output of the model fitted with data_jags
#' @param file.in Path and file where to save the plot
plot_fitted_vs_true_parameters <- function(data_jags, jags.model, file.in){
  
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  ## - Make the plot
  #  -- Get fitted parameters
  fitted.paramater.in <- as.data.frame(jags.model$BUGSoutput$sims.matrix) %>%
    dplyr::select(-deviance) %>%
    gather(key = "parameter", value = "value") 
  # -- Get true parameters
  true.parameter.in <- data.frame(parameter = names(data_jags$param), 
                                  value = as.numeric(unlist(data_jags$param)))
  # -- Plot fitted vs true
  plot.in <- fitted.paramater.in %>%
    ggplot(aes(x = value)) + 
    geom_histogram(aes(y = stat(density)), 
                   colour = "black", bins = 30) + 
    facet_wrap(~ parameter, scales = "free") + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(data = true.parameter.in, 
               aes(xintercept = value), 
               color = "red") + 
    xlab("Parameter value") +
    theme_bw()
  
  ## - save the plot
  ggsave(file.in, plot.in, width = 17, height = 12, units = "cm", dpi = 600)
  return(file.in)
}


#' Compare generated and fitted probabilities
#' @param data_jags input for the jags model
#' @param jags.model output of the jags model
#' @param data_model_scaled data centered and scaled
#' @param file.in Path and file where to save the plot
compare_fitted_generated_probabilities <- function(data_jags, jags.model, data_model_scaled, 
                                                   file.in){
  ## - Create directories if needed
  create_dir_if_needed(file.in)
  
  # Generated probabilities
  proba_generated <- data_model_scaled %>%
    # Compute intermediate probabilities based on parameters and data
    mutate(pdD = data_jags$param$c0 + data_jags$param$c1*DS + data_jags$param$c2*DS*dbh, 
           pdBM = data_jags$param$b0 + data_jags$param$b1*dbh + data_jags$param$b2*comp + data_jags$param$b3*sgdd + data_jags$param$b4*wai, 
           phdD = data_jags$data$d0 + data_jags$data$d1*dbh + data_jags$data$d2*DS + data_jags$data$d3*dbh*DS, 
           phadBM = data_jags$data$e0 + data_jags$data$e1*dbh) %>%
    # Apply inverse logit function to constrain probabilities between 0 and 1
    mutate(pdD = exp(pdD)/(1 + exp(pdD)), 
           pdBM = exp(pdBM)/(1 + exp(pdBM)), 
           phdD = exp(phdD)/(1 + exp(phdD)), 
           phadBM = exp(phadBM)/(1 + exp(phadBM))) %>%
    mutate(pdDj = 1 - (1 - pdD)*(1 - pdBM)) %>%
    # Compute probability to be alive, dead or harvested
    mutate(ph = (1 - D)*phadBM + D*pdDj*phdD, 
           pd = (1 - D)*pdBM*(1 - phadBM) + D*pdDj*(1 - phdD)) %>%
    mutate(pa = 1 - (ph + pd)) %>%
    dplyr::select(ph, pd, pa, pdD, pdBM, pdDj)
  
  # Fitted parameters
  param_fitted <- ggs(as.mcmc(jags.model)) %>%
    filter(Parameter != "deviance") %>%
    spread(key = Parameter, value = value) 
  # Fitted probabilities
  proba_fitted <- array(NA_real_, dim = c(dim(proba_generated), dim(param_fitted)[1]))
  for(i in 1:dim(param_fitted)[1]){
    if(floor(i/1000) == i/1000) print(paste0("Simulation ", i, "/", dim(param_fitted)[1]))
    proba.sim.i <- data_model_scaled %>%
      # Compute intermediate probabilities based on parameters and data
      mutate(pdD = param_fitted$c0[i] + param_fitted$c1[i]*DS + param_fitted$c2[i]*DS*dbh, 
             pdBM = param_fitted$b0[i] + param_fitted$b1[i]*dbh + param_fitted$b2[i]*comp + param_fitted$b3[i]*sgdd + param_fitted$b4[i]*wai, 
             phdD = data_jags$data$d0 + data_jags$data$d1*dbh + data_jags$data$d2*DS + data_jags$data$d3*dbh*DS, 
             phadBM = data_jags$data$e0 + data_jags$data$e1*dbh) %>%
      # Apply inverse logit function to constrain probabilities between 0 and 1
      mutate(pdD = exp(pdD)/(1 + exp(pdD)), 
             pdBM = exp(pdBM)/(1 + exp(pdBM)), 
             phdD = exp(phdD)/(1 + exp(phdD)), 
             phadBM = exp(phadBM)/(1 + exp(phadBM))) %>%
      mutate(pdDj = 1 - (1 - pdD)*(1 - pdBM)) %>%
      # Compute probability to be alive, dead or harvested
      mutate(ph = (1 - D)*phadBM + D*pdDj*phdD, 
             pd = (1 - D)*pdBM*(1 - phadBM) + D*pdDj*(1 - phdD)) %>%
      mutate(pa = 1 - (ph + pd)) %>%
      dplyr::select(ph, pd, pa, pdD, pdBM, pdDj)
    proba_fitted[, ,i] <- as.matrix(proba.sim.i)
  }
  proba_fitted <- as.data.frame(apply(proba_fitted, c(1, 2), mean)) 
  colnames(proba_fitted) <-  colnames(proba_generated)
  
  plot.in <- cbind((proba_fitted %>% gather(key = "probability", value = "fitted_value")), 
        (proba_generated %>% gather(key = "proba", value = "generated_value"))) %>%
    dplyr::select(-proba) %>%
    ggplot(aes(x = generated_value, y = fitted_value)) + 
    geom_point() +
    facet_wrap(~ probability, scales = "free") + 
    theme_bw()
  
  ## - save the plot
  ggsave(file.in, plot.in, width = 17, height = 12, units = "cm", dpi = 600)
  return(file.in)
}
