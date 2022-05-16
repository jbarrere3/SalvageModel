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










#' Fit the mortality model with disturbance event given in the data for several countries
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_full_sub <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, param.in){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pdD[i])
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstorm[i]) = c0[sp[i]] + (c1[sp[i]]*Istorm[plot[i]]*(dbh[i]^c2[sp[i]])) 
      pdstorm[i] = 1 - (1 - pstorm[i])^time[i]
      
      # Probability to die from an other disturbance
      logit(pother[i]) = c3[sp[i]] + (c4[sp[i]]*Iother[plot[i]]*(dbh[i]^c5[sp[i]])) 
      pdother[i] = 1 - (1 - pother[i])^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfire[i]) = c6[sp[i]] + (c7[sp[i]]*Ifire[plot[i]]*(dbh[i]^c8[sp[i]])) 
      pdfire[i] = 1 - (1 - pfire[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      c0[s] ~ dnorm(0, 0.1)
      c1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c2[s] ~ dnorm(0, 1) 
      c3[s] ~ dnorm(0, 0.1)
      c4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c5[s] ~ dnorm(0, 1) 
      c6[s] ~ dnorm(0, 0.1)
      c7[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c8[s] ~ dnorm(0, 1) 
    }
    
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      Ifire[k] ~ dbeta(0.66, 0.36) T(0.001,0.999)
      Istorm[k] ~ dbeta(0.65, 2.66) T(0.001,0.999)
      Iother[k] ~ dbeta(0.48, 1.77) T(0.001,0.999)
    }
    
    
    
    
  }"
  
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model_D, tmp)
  out <- R2jags::jags.parallel(data = data_jags.in,
                               param = param.in,
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


#' Fit the mortality model with disturbance event given in the data for several countries and a climate effect
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_full_sub_climate <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, 
                                           param.in = c(paste0("st", c(0:4)), 
                                                        paste0("ot", c(0:4)), 
                                                        paste0("fi", c(0:4)))){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pdD[i])
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstorm[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]])) + st3[sp[i]]*sgdd[i] + st4[sp[i]]*wai[i]
      pdstorm[i] = 1 - (1 - pstorm[i])^time[i]
      
      # Probability to die from an other disturbance
      logit(pother[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]])) + ot3[sp[i]]*sgdd[i] + ot4[sp[i]]*wai[i]
      pdother[i] = 1 - (1 - pother[i])^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfire[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]])) + fi3[sp[i]]*sgdd[i] + fi4[sp[i]]*wai[i]
      pdfire[i] = 1 - (1 - pfire[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.5)
      st4[s] ~ dnorm(0, 0.5)
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.5)
      ot4[s] ~ dnorm(0, 0.5)
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.5)
      fi4[s] ~ dnorm(0, 0.5)
      
    
    }
    
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      Ifire[k] ~ dbeta(0.66, 0.36) T(0.001,0.999)
      Istorm[k] ~ dbeta(0.65, 2.66) T(0.001,0.999)
      Iother[k] ~ dbeta(0.48, 1.77) T(0.001,0.999)
    }
    
    
    
    
  }"
  
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model_D, tmp)
  out <- R2jags::jags.parallel(data = data_jags.in,
                               param = param.in,
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


#' Fit the mortality model with disturbance event given in the data for several countries and a climate effect
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_full_sub_climate2 <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, 
                                           param.in = c(paste0("st", c(0:3)), 
                                                        paste0("ot", c(0:3)), 
                                                        paste0("fi", c(0:3)))){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pdD[i])
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstorm[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]])) + st3[sp[i]]*sgdd[i]
      pdstorm[i] = 1 - (1 - pstorm[i])^time[i]
      
      # Probability to die from an other disturbance
      logit(pother[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]])) + ot3[sp[i]]*sgdd[i]
      pdother[i] = 1 - (1 - pother[i])^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfire[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]])) + fi3[sp[i]]*sgdd[i]
      pdfire[i] = 1 - (1 - pfire[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.5)
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.5)
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.5)

    
    }
    
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      Ifire[k] ~ dbeta(0.66, 0.36) T(0.001,0.999)
      Istorm[k] ~ dbeta(0.65, 2.66) T(0.001,0.999)
      Iother[k] ~ dbeta(0.48, 1.77) T(0.001,0.999)
    }
    
    
    
    
  }"
  
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model_D, tmp)
  out <- R2jags::jags.parallel(data = data_jags.in,
                               param = param.in,
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



#' Fit the mortality model with disturbance event given in the data for several countries
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_full_sub_simulated <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pdD[i])
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstorm[i]) = c0[sp[i]] + (c1[sp[i]]*Istorm[i]*(dbh[i]^c2[sp[i]])) 
      pdstorm[i] = 1 - (1 - pstorm[i])^time[i]
      
      # Probability to die from an other disturbance
      logit(pother[i]) = c3[sp[i]] + (c4[sp[i]]*Iother[i]*(dbh[i]^c5[sp[i]])) 
      pdother[i] = 1 - (1 - pother[i])^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfire[i]) = c6[sp[i]] + (c7[sp[i]]*Ifire[i]*(dbh[i]^c8[sp[i]])) 
      pdfire[i] = 1 - (1 - pfire[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      c0[s] ~ dnorm(0, 0.1)
      c1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c2[s] ~ dnorm(0, 1) 
      c3[s] ~ dnorm(0, 0.1)
      c4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c5[s] ~ dnorm(0, 1) 
      c6[s] ~ dnorm(0, 0.1)
      c7[s] ~ dnorm(0, 0.1) T(0.01, 100)
      c8[s] ~ dnorm(0, 1) 
    
    }
    
    
    
  }"
  
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model_D, tmp)
  out <- R2jags::jags.parallel(data = data_jags.in,
                               param = paste0("c", c(0:8)),
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


#' Fit the mortality model with disturbance event given in the data for several countries
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_full_sub_climate_simulated <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pdD[i])
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstorm[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[i]*(dbh[i]^st2[sp[i]])) + st3[sp[i]]*sgdd[i] + st4[sp[i]]*wai[i]
      pdstorm[i] = 1 - (1 - pstorm[i])^time[i]
      
      # Probability to die from an other disturbance
      logit(pother[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[i]*(dbh[i]^ot2[sp[i]])) + ot3[sp[i]]*sgdd[i] + ot4[sp[i]]*wai[i]
      pdother[i] = 1 - (1 - pother[i])^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfire[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[i]*(dbh[i]^fi2[sp[i]])) + fi3[sp[i]]*sgdd[i] + fi4[sp[i]]*wai[i]
      pdfire[i] = 1 - (1 - pfire[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.5)
      st4[s] ~ dnorm(0, 0.5)
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.5)
      ot4[s] ~ dnorm(0, 0.5)
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.5)
      fi4[s] ~ dnorm(0, 0.5)
      
    
    }
    
    
    
  }"
  
  
  
  ## - Fit the model in parallel
  tmp <- tempfile()
  writeLines(mortality_model_D, tmp)
  out <- R2jags::jags.parallel(data = data_jags.in,
                               param = c(paste0("st", c(0:4)), 
                                         paste0("ot", c(0:4)), 
                                         paste0("fi", c(0:4))),
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





