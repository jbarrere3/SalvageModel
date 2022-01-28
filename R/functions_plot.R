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
