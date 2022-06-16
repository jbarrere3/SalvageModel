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








#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_dqm3 <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, 
                               param.in = c(paste0("st", c(0:5)), 
                                            paste0("ot", c(0:5)), 
                                            paste0("fi", c(0:5)))){
  
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
      logit(pstormdomi[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]]))
      logit(pstorm[i]) = st3[sp[i]] + (st4[sp[i]]*Istorm[plot[i]]*(dbh[i]^st5[sp[i]]))
      pdstorm[i] = 1 - (1 - (domi[i]*pstormdomi[i] + (1 - domi[i])*pstorm[i]))^time[i]
      
      # Probability to die from an other disturbance
      logit(potherdomi[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]]))
      logit(pother[i]) = ot3[sp[i]] + (ot4[sp[i]]*Iother[plot[i]]*(dbh[i]^ot5[sp[i]]))
      pdother[i] = 1 - (1 - (domi[i]*potherdomi[i] + (1 - domi[i])*pother[i]))^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfiredomi[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]]))
      logit(pfire[i]) = fi3[sp[i]] + (fi4[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi5[sp[i]]))
      pdfire[i] = 1 - (1 - (domi[i]*pfiredomi[i] + (1 - domi[i])*pfire[i]))^time[i]
      
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.1)
      st4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st5[s] ~ dnorm(0, 1) 
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.1)
      ot4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot5[s] ~ dnorm(0, 1) 
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.1)
      fi4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi5[s] ~ dnorm(0, 1) 
      

    
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


#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_dqm3_climate <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, 
                               param.in = c(paste0("st", c(0:7)), 
                                            paste0("ot", c(0:7)), 
                                            paste0("fi", c(0:7)))){
  
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
      logit(pstormdomi[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]])) + st6[sp[i]]*sgdd[i]
      logit(pstorm[i]) = st3[sp[i]] + (st4[sp[i]]*Istorm[plot[i]]*(dbh[i]^st5[sp[i]])) + st7[sp[i]]*sgdd[i]
      pdstorm[i] = 1 - (1 - (domi[i]*pstormdomi[i] + (1 - domi[i])*pstorm[i]))^time[i]
      
      # Probability to die from an other disturbance
      logit(potherdomi[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]])) + ot6[sp[i]]*sgdd[i]
      logit(pother[i]) = ot3[sp[i]] + (ot4[sp[i]]*Iother[plot[i]]*(dbh[i]^ot5[sp[i]])) + ot7[sp[i]]*sgdd[i]
      pdother[i] = 1 - (1 - (domi[i]*potherdomi[i] + (1 - domi[i])*pother[i]))^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfiredomi[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]])) + fi6[sp[i]]*sgdd[i]
      logit(pfire[i]) = fi3[sp[i]] + (fi4[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi5[sp[i]])) + fi7[sp[i]]*sgdd[i]
      pdfire[i] = 1 - (1 - (domi[i]*pfiredomi[i] + (1 - domi[i])*pfire[i]))^time[i]
      
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.1)
      st4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st5[s] ~ dnorm(0, 1) 
      st6[s] ~ dnorm(0, 0.1) 
      st7[s] ~ dnorm(0, 0.1) 
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.1)
      ot4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot5[s] ~ dnorm(0, 1) 
      ot6[s] ~ dnorm(0, 0.1) 
      ot7[s] ~ dnorm(0, 0.1) 
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.1)
      fi4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi5[s] ~ dnorm(0, 1) 
      fi6[s] ~ dnorm(0, 0.1) 
      fi7[s] ~ dnorm(0, 0.1) 
      

    
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






#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_dqm4 <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, 
                               param.in = c(paste0("st", c(0:5)), 
                                            paste0("ot", c(0:5)), 
                                            paste0("fi", c(0:5)))){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbin(pdD[i], w[i])
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstormdomi[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]]))
      logit(pstorm[i]) = st3[sp[i]] + (st4[sp[i]]*Istorm[plot[i]]*(dbh[i]^st5[sp[i]]))
      pdstorm[i] = 1 - (1 - (domi[i]*pstormdomi[i] + (1 - domi[i])*pstorm[i]))^time[i]
      
      # Probability to die from an other disturbance
      logit(potherdomi[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]]))
      logit(pother[i]) = ot3[sp[i]] + (ot4[sp[i]]*Iother[plot[i]]*(dbh[i]^ot5[sp[i]]))
      pdother[i] = 1 - (1 - (domi[i]*potherdomi[i] + (1 - domi[i])*pother[i]))^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfiredomi[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]]))
      logit(pfire[i]) = fi3[sp[i]] + (fi4[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi5[sp[i]]))
      pdfire[i] = 1 - (1 - (domi[i]*pfiredomi[i] + (1 - domi[i])*pfire[i]))^time[i]
      
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.1)
      st4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st5[s] ~ dnorm(0, 1) 
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.1)
      ot4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot5[s] ~ dnorm(0, 1) 
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.1)
      fi4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi5[s] ~ dnorm(0, 1) 
      

    
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


#' Fit a mortality model including the mean quadratic diameter (dqm)
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param param.in parameters to extract
#' @return A rjags object
fit_mortality_dqm5 <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, 
                               param.in = c(paste0("st", c(0:5)), 
                                            paste0("ot", c(0:5)), 
                                            paste0("fi", c(0:5)))){
  
  # Initialize time
  start <- Sys.time()
  
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      # Probability that the tree died from a disturbance
      d[i] ~ dbern(pdD[i])
      pdD[i] ~ dbeta(mu[i]*(100/w[i]), (1 - mu[i])*(100/w[i]))
      mu[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstormdomi[i]) = st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]]))
      logit(pstorm[i]) = st3[sp[i]] + (st4[sp[i]]*Istorm[plot[i]]*(dbh[i]^st5[sp[i]]))
      pdstorm[i] = 1 - (1 - (domi[i]*pstormdomi[i] + (1 - domi[i])*pstorm[i]))^time[i]
      
      # Probability to die from an other disturbance
      logit(potherdomi[i]) = ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]]))
      logit(pother[i]) = ot3[sp[i]] + (ot4[sp[i]]*Iother[plot[i]]*(dbh[i]^ot5[sp[i]]))
      pdother[i] = 1 - (1 - (domi[i]*potherdomi[i] + (1 - domi[i])*pother[i]))^time[i]
      
      # Probability to die from a fire disturbance
      logit(pfiredomi[i]) = fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]]))
      logit(pfire[i]) = fi3[sp[i]] + (fi4[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi5[sp[i]]))
      pdfire[i] = 1 - (1 - (domi[i]*pfiredomi[i] + (1 - domi[i])*pfire[i]))^time[i]
      
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ dnorm(0, 0.1)
      st1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st2[s] ~ dnorm(0, 1) 
      st3[s] ~ dnorm(0, 0.1)
      st4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      st5[s] ~ dnorm(0, 1) 
      ot0[s] ~ dnorm(0, 0.1)
      ot1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot2[s] ~ dnorm(0, 1) 
      ot3[s] ~ dnorm(0, 0.1)
      ot4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      ot5[s] ~ dnorm(0, 1) 
      fi0[s] ~ dnorm(0, 0.1)
      fi1[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi2[s] ~ dnorm(0, 1) 
      fi3[s] ~ dnorm(0, 0.1)
      fi4[s] ~ dnorm(0, 0.1) T(0.01, 100)
      fi5[s] ~ dnorm(0, 1) 
      

    
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



