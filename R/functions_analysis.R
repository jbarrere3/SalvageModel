#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#### SCRIPT INTRODUCTION ####
#
#' @name functions_analysis.R  
#' @description R script containing all functions relative to data
#               analysis
#
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 1. Functions for the reference model ------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jagsList containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_storm <- function(data_jags, n.chains, n.iter, n.burn, n.thin, 
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
      logit(p[i]) = a0[sp[i], co[i]] + a1[sp[i]]*logratio[i] + b[sp[i]]*I[plot[i]]*(dbh[i]^c[sp[i]])
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
fit_mortality_other <- function(data_jags, n.chains, n.iter, n.burn, n.thin, 
                                param = c("a0", "b", "c", "I")){
  
  
  
  ## - Initialize time
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
  mortality_model <- paste0("model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pd[i])
      # Probability to die from a storm disturbance
      logit(p[i]) = a0[sp[i], co[i]] + b[sp[i]]*I[plot[i]]*(dbh[i]^c[sp[i]])
      pd[i] = 1 - (1 - p[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
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
  }")
  
  
  
  
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
fit_mortality <- function(data_jags, n.chains, n.iter, n.burn, n.thin){
  
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
    # And one "classical model" for fire and other disturbance
    if(disturbances.in[i] %in% c("fire", "other", "biotic")){
      jags.i <- fit_mortality_other(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    }
    # Add model to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- jags.i")))
  }
  
  # Return output
  return(out)
}






#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## 2. Functions for the model with stocking --------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#' Fit a mortality model including stocking
#' @param data_jagsList containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_storm_stock <- function(data_jags, n.chains, n.iter, n.burn, n.thin, 
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
      logit(p[i]) = a0[sp[i], co[i]] + a1[sp[i]]*stock[i] + b[sp[i]]*I[plot[i]]*(dbh[i]^c[sp[i]])
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
fit_mortality_stock <- function(data_jags, n.chains, n.iter, n.burn, n.thin){
  
  # Identify disturbances
  disturbances.in <- names(data_jags)
  
  # Initialize output
  out <- list()
  
  # Loop on all disturbances
  for(i in 1:length(disturbances.in)){
    
    print(paste0("Fit mortality model for ", disturbances.in[i], " disturbance"))
    
    # There's a specific model for storm
    jags.i <- fit_mortality_storm_stock(data_jags[[i]][[1]], n.chains, n.iter, n.burn, n.thin)
    
    # Add model to the output list
    eval(parse(text = paste0("out$", disturbances.in[i], " <- jags.i")))
  }
  
  # Return output
  return(out)
}


