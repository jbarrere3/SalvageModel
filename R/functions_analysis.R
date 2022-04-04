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

#' Fit the mortality model with disturbance event as a latent variable
#' @param data_jags.in List contianing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the mortality model 
  mortality_model <- function() {
    # Loop over all plots
    for (j in 1:Nplots) {
      D[j] ~ dbern(pD[j])
      logit(pD[j]) = a0 + a1*DA[j]
    }
    
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[plot[i]])*phadBM[i] + D[plot[i]]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[plot[i]])*pdBM[i]*(1 - phadBM[i]) + D[plot[i]]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[plot[i]])*(1 - pdBM[i])*(1 - phadBM[i]) + D[plot[i]]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(pdBM[i]) = b0 + b1*dbh[i] + b2*comp[i] + b3*sgdd[i] + b4*wai[i]
      # Probability to die from a disturbance
      logit(pdD[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
    }
    
    # Priors
    a0 ~ dnorm(0, 0.01)
    a1 ~ dnorm(0, 0.01)
    b0 ~ dnorm(0, 0.01)
    b1 ~ dnorm(0, 0.01)
    b2 ~ dnorm(0, 0.01)
    b3 ~ dnorm(0, 0.01)
    b4 ~ dnorm(0, 0.01)
    c0 ~ dnorm(0, 0.01)
    c1 ~ dnorm(0, 0.01)
    c2 ~ dnorm(0, 0.01)
  }
  
  ## - Function to initialize priors
  initjags <- function(){
    return(list(a0 = runif(1, -5, 5),
                a1 = runif(1, 0, 5),
                b0 = runif(1, -5, 5),
                b1 = runif(1, -5, 0),
                b2 = runif(1, 2, 5),
                b3 = runif(1, -5, 0),
                b4 = runif(1, 0, 5),
                c0 = runif(1, -5, 5),
                c1 = runif(1, 2, 5),
                c2 = runif(1, -5, 0)))
  }
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = names(initjags()),
                      inits = list(initjags(),
                                   initjags(),
                                   initjags()),
                      model.file = mortality_model,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}







#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(pdBM[i]) = b0 + b1*dbh[i] + b2*comp[i] + b3*sgdd[i] + b4*wai[i]
      # Probability to die from a disturbance
      logit(pdD[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
    }
    
    # Priors
    b0 ~ dnorm(0, 0.01)
    b1 ~ dnorm(0, 0.01)
    b2 ~ dnorm(0, 0.01)
    b3 ~ dnorm(0, 0.01)
    b4 ~ dnorm(0, 0.01)
    c0 ~ dnorm(0, 1)
    c1 ~ dnorm(0, 1)
    c2 ~ dnorm(0, 1)
  }
  
  ## - Function to initialize priors
  initjags_D <- function(){
    return(list(
      b0 = runif(1, -5, 5),
      b1 = runif(1, -5, 0),
      b2 = runif(1, 2, 5),
      b3 = runif(1, -5, 0),
      b4 = runif(1, 0, 5),
      c0 = runif(1, -5, 5),
      c1 = runif(1, 2, 5),
      c2 = runif(1, -5, 0)))
  }
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = names(initjags_D()),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}


#' Fit the mortality model with parameters (a0 and a1) for disturbance occurrence given
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_a <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the mortality model 
  mortality_model_a <- function() {
    
    # Loop over all plots
    for (j in 1:Nplots) {
      D[j] ~ dbern(pD[j])
    }
    
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[plot[i]])*phadBM[i] + D[plot[i]]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[plot[i]])*pdBM[i]*(1 - phadBM[i]) + D[plot[i]]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[plot[i]])*(1 - pdBM[i])*(1 - phadBM[i]) + D[plot[i]]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(pdBM[i]) = b0 + b1*dbh[i] + b2*comp[i] + b3*sgdd[i] + b4*wai[i]
      # Probability to die from a disturbance
      logit(pdD[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
    }
    
    # Priors
    b0 ~ dnorm(0, 0.01)
    b1 ~ dnorm(0, 0.01)
    b2 ~ dnorm(0, 0.01)
    b3 ~ dnorm(0, 0.01)
    b4 ~ dnorm(0, 0.01)
    c0 ~ dnorm(0, 0.01)
    c1 ~ dnorm(0, 0.01)
    c2 ~ dnorm(0, 0.01)
  }
  
  ## - Function to initialize priors
  initjags <- function(){
    return(list(b0 = runif(1, -5, 5),
                b1 = runif(1, -5, 0),
                b2 = runif(1, 2, 5),
                b3 = runif(1, -5, 0),
                b4 = runif(1, 0, 5),
                c0 = runif(1, -5, 5),
                c1 = runif(1, 2, 5),
                c2 = runif(1, -5, 0)))
  }
  
  ## - Fit the model
  out_a <- R2jags::jags(data = data_jags.in,
                        param = names(initjags()),
                        inits = list(initjags(),
                                     initjags(),
                                     initjags()),
                        model.file = mortality_model_a,
                        n.chains = n.chains,
                        n.iter = n.iter,
                        n.burnin = n.burn,
                        n.thin = n.thin,
                        DIC = TRUE, 
                        progress.bar = "text")
  
  return(out_a)
}



#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_BM1_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(mu1[i]) = b0 + b1*dbh[i] + b2*log(dbh[i]) + b3*comp[i] + b4*(1/sgdd[i]) + b5*(1/wai[i])
      pdBM[i] = 1 - (1 - mu1[i])^5
      # Probability to die from a disturbance
      logit(mu2[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    # Priors
    b0 ~ dnorm(0, 0.5)
    b1 ~ dnorm(0, 0.5)
    b2 ~ dnorm(0, 0.5)
    b3 ~ dnorm(0, 0.5)
    b4 ~ dnorm(0, 0.5)
    b5 ~ dnorm(0, 0.5)
    c0 ~ dnorm(0, 0.5)
    c1 ~ dnorm(0, 0.5)
    c2 ~ dnorm(0, 0.5)
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:5)), paste0("c", c(0:2))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}

#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_BM2_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(mu1[i]) = b0 + b1*dbh[i] + b2*log(dbh[i]) + b3*comp[i] + b4*(1/sgdd[i]) + 
        b5*(1/wai[i]) + b6*(1/sgdd[i])*dbh[i] + b7*(1/sgdd[i])*log(dbh[i]) +
        b8*(1/wai[i])*dbh[i] + b9*(1/wai[i])*log(dbh[i])
      pdBM[i] = 1 - (1 - mu1[i])^5
      # Probability to die from a disturbance
      logit(mu2[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    # Priors
    b0 ~ dnorm(0, 0.5)
    b1 ~ dnorm(0, 0.5)
    b2 ~ dnorm(0, 0.5)
    b3 ~ dnorm(0, 0.5)
    b4 ~ dnorm(0, 0.5)
    b5 ~ dnorm(0, 0.5)
    b6 ~ dnorm(0, 0.5)
    b7 ~ dnorm(0, 0.5)
    b8 ~ dnorm(0, 0.5)
    b9 ~ dnorm(0, 0.5)
    c0 ~ dnorm(0, 0.5)
    c1 ~ dnorm(0, 0.5)
    c2 ~ dnorm(0, 0.5)
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:9)), paste0("c", c(0:2))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}

#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_BM3_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(mu1[i]) = b0 + b1*dbh[i] + b2*log(dbh[i]) + b3*comp[i] + b4*(1/sgdd[i]) + 
        b5*(1/wai[i]) + b6*(1/sgdd[i])*dbh[i] + b7*(1/sgdd[i])*log(dbh[i]) + b8*(1/wai[i])*dbh[i] + 
        b9*(1/wai[i])*log(dbh[i]) + b10*(1/wai[i])*comp[i] + b11*(1/sgdd[i])*comp[i]
      pdBM[i] = 1 - (1 - mu1[i])^5
      # Probability to die from a disturbance
      logit(mu2[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    # Priors
    b0 ~ dnorm(0, 0.5)
    b1 ~ dnorm(0, 0.5)
    b2 ~ dnorm(0, 0.5)
    b3 ~ dnorm(0, 0.5)
    b4 ~ dnorm(0, 0.5)
    b5 ~ dnorm(0, 0.5)
    b6 ~ dnorm(0, 0.5)
    b7 ~ dnorm(0, 0.5)
    b8 ~ dnorm(0, 0.5)
    b9 ~ dnorm(0, 0.5)
    b10 ~ dnorm(0, 0.5)
    b11 ~ dnorm(0, 0.5)
    c0 ~ dnorm(0, 0.5)
    c1 ~ dnorm(0, 0.5)
    c2 ~ dnorm(0, 0.5)
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:11)), paste0("c", c(0:2))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}

#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_BM4_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(mu1[i]) = b0 + b1*dbh[i] + b2*log(dbh[i]) + b3*comp[i] + b4*sgdd[i] + 
        b5*sgdd[i]^2 + b6*wai[i] + b7*wai[i]^2
      pdBM[i] = 1 - (1 - mu1[i])^5
      # Probability to die from a disturbance
      logit(mu2[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    # Priors
    b0 ~ dnorm(0, 0.5)
    b1 ~ dnorm(0, 0.5)
    b2 ~ dnorm(0, 0.5)
    b3 ~ dnorm(0, 0.5)
    b4 ~ dnorm(0, 0.5)
    b5 ~ dnorm(0, 0.5)
    b6 ~ dnorm(0, 0.5)
    b7 ~ dnorm(0, 0.5)
    c0 ~ dnorm(0, 0.5)
    c1 ~ dnorm(0, 0.5)
    c2 ~ dnorm(0, 0.5)
    
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:7)), paste0("c", c(0:2))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}


#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_BM5_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(mu1[i]) = b0 + b1*dbh[i] + b2*log(dbh[i]) + b3*comp[i] + b4*sgdd[i] + 
        b5*sgdd[i]^2 + b6*wai[i] + b7*wai[i]^2 + b8*sgdd[i]*dbh[i] + b9*sgdd[i]*log(dbh[i]) + 
        b10*wai[i]*dbh[i] + b11*wai[i]*log(dbh[i])
      pdBM[i] = 1 - (1 - mu1[i])^5
      # Probability to die from a disturbance
      logit(mu2[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    # Priors
    b0 ~ dnorm(0, 0.5)
    b1 ~ dnorm(0, 0.5)
    b2 ~ dnorm(0, 0.5)
    b3 ~ dnorm(0, 0.5)
    b4 ~ dnorm(0, 0.5)
    b5 ~ dnorm(0, 0.5)
    b6 ~ dnorm(0, 0.5)
    b7 ~ dnorm(0, 0.5)
    b8 ~ dnorm(0, 0.5)
    b9 ~ dnorm(0, 0.5)
    b10 ~ dnorm(0, 0.5)
    b11 ~ dnorm(0, 0.5)
    c0 ~ dnorm(0, 0.5)
    c1 ~ dnorm(0, 0.5)
    c2 ~ dnorm(0, 0.5)
    
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:11)), paste0("c", c(0:2))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}


#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_BM6_D <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      # Probability to die from background mortality
      logit(mu1[i]) = b0 + b1*dbh[i] + b2*log(dbh[i]) + b3*comp[i] + b4*sgdd[i] + 
        b5*sgdd[i]^2 + b6*wai[i] + b7*wai[i]^2 + b8*sgdd[i]*dbh[i] + b9*sgdd[i]*log(dbh[i]) + 
        b10*wai[i]*dbh[i] + b11*wai[i]*log(dbh[i]) + b12*sgdd[i]*comp[i] + b13*wai[i]*comp[i]
      pdBM[i] = 1 - (1 - mu1[i])^5
      # Probability to die from a disturbance
      logit(mu2[i]) = c0 + c1*DS[i] + c2*dbh[i]*DS[i]
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    # Priors
    b0 ~ dnorm(0, 0.5)
    b1 ~ dnorm(0, 0.5)
    b2 ~ dnorm(0, 0.5)
    b3 ~ dnorm(0, 0.5)
    b4 ~ dnorm(0, 0.5)
    b5 ~ dnorm(0, 0.5)
    b6 ~ dnorm(0, 0.5)
    b7 ~ dnorm(0, 0.5)
    b8 ~ dnorm(0, 0.5)
    b9 ~ dnorm(0, 0.5)
    b10 ~ dnorm(0, 0.5)
    b11 ~ dnorm(0, 0.5)
    b12 ~ dnorm(0, 0.5)
    b13 ~ dnorm(0, 0.5)
    c0 ~ dnorm(0, 0.5)
    c1 ~ dnorm(0, 0.5)
    c2 ~ dnorm(0, 0.5)
    
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:13)), paste0("c", c(0:2))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}



#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality2 <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- function() {
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = (1 - D[i])*phadBM[i] + D[i]*pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = (1 - D[i])*pdBM[i]*(1 - phadBM[i]) + D[i]*pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - D[i])*(1 - pdBM[i])*(1 - phadBM[i]) + D[i]*(1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      
      ## - Probability to die from background mortality
      logit(mu1[i]) = b0[sp[i]] + b1[sp[i]]*dbh[i] + b2[sp[i]]*log(dbh[i]) + b3[sp[i]]*comp[i] + 
        # 1st equation of Background mortality
        bm1[i]*(b4[sp[i]]*(1/sgdd[i]) + b5[sp[i]]*(1/wai[i])) + 
        # 2nd equation of Background mortality
        bm2[i]*(b4[sp[i]]*(1/sgdd[i]) + b5[sp[i]]*(1/wai[i]) + b6[sp[i]]*(1/sgdd[i])*dbh[i] + 
                  b7[sp[i]]*(1/sgdd[i])*log(dbh[i])+ b8[sp[i]]*(1/wai[i])*dbh[i] + 
                  b9[sp[i]]*(1/wai[i])*log(dbh[i])) + 
        # 3rd equation of Background mortality
        bm3[i]*(b4[sp[i]]*(1/sgdd[i]) + b5[sp[i]]*(1/wai[i]) + b6[sp[i]]*(1/sgdd[i])*dbh[i] + 
                  b7[sp[i]]*(1/sgdd[i])*log(dbh[i])+ b8[sp[i]]*(1/wai[i])*dbh[i] + 
                  b9[sp[i]]*(1/wai[i])*log(dbh[i]) + b10[sp[i]]*(1/wai[i])*comp[i] + 
                  b11[sp[i]]*(1/sgdd[i])*comp[i]) + 
        # 4th equation of Background mortality
        bm4[i]*(b4[sp[i]]*sgdd[i] + b5[sp[i]]*(sgdd[i]^2) + 
                  b6[sp[i]]*wai[i] + b7[sp[i]]*(wai[i]^2)) + 
        # 5th equation of Background mortality
        bm5[i]*(b4[sp[i]]*sgdd[i] + b5[sp[i]]*(sgdd[i]^2) + 
                  b6[sp[i]]*wai[i] + b7[sp[i]]*(wai[i]^2) + 
                  b8[sp[i]]*sgdd[i]*dbh[i] + b9[sp[i]]*sgdd[i]*log(dbh[i]) + 
                  b10[sp[i]]*wai[i]*dbh[i] + b11[sp[i]]*wai[i]*log(dbh[i])) + 
        # 6th equation of Background mortality
        bm6[i]*(b4[sp[i]]*sgdd[i] + b5[sp[i]]*(sgdd[i]^2) + 
                  b6[sp[i]]*wai[i] + b7[sp[i]]*(wai[i]^2) + 
                  b8[sp[i]]*sgdd[i]*dbh[i] + b9[sp[i]]*sgdd[i]*log(dbh[i]) + 
                  b10[sp[i]]*wai[i]*dbh[i] + b11[sp[i]]*wai[i]*log(dbh[i]) + 
                  b12[sp[i]]*sgdd[i]*comp[i] + b13[sp[i]]*wai[i]*comp[i])
      # Add the "year" effect
      pdBM[i] = 1 - (1 - mu1[i])^5
      
      ## - Probability to die from a disturbance
      logit(mu2[i]) = c0[sp[i]] + 
        Dstorm[i]*(c1[sp[i]]*DI[plot[i]]*(dbh[i]^c2[sp[i]])) + 
        Dfire[i]*(c3[sp[i]]*DI[plot[i]]*(dbh[i]^c4[sp[i]])) + 
        Dother[i]*(c5[sp[i]]*DI[plot[i]]*(dbh[i]^c6[sp[i]]))
      pdD[i] = 1 - (1 - mu2[i])^5
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      b0[s] ~ dnorm(0, 0.5)
      b1[s] ~ dnorm(0, 0.5)
      b2[s] ~ dnorm(0, 0.5)
      b3[s] ~ dnorm(0, 0.5)
      b4[s] ~ dnorm(0, 0.5)
      b5[s] ~ dnorm(0, 0.5)
      b6[s] ~ dnorm(0, 0.5)
      b7[s] ~ dnorm(0, 0.5)
      b8[s] ~ dnorm(0, 0.5)
      b9[s] ~ dnorm(0, 0.5)
      b10[s] ~ dnorm(0, 0.5)
      b11[s] ~ dnorm(0, 0.5)
      b12[s] ~ dnorm(0, 0.5)
      b13[s] ~ dnorm(0, 0.5)
      c0[s] ~ dnorm(0, 0.5)
      c1[s] ~ dnorm(0, 0.5)
      c2[s] ~ dnorm(0, 0.5)
      c3[s] ~ dnorm(0, 0.5)
      c4[s] ~ dnorm(0, 0.5)
      c5[s] ~ dnorm(0, 0.5)
      c6[s] ~ dnorm(0, 0.5)
    }
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      DI[k] ~ dunif(0, 1)
    }
    
  }
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = c(paste0("b", c(0:13)), paste0("c", c(0:6))),
                      model.file = mortality_model_D,
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}

