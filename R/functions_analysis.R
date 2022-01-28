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

#' Fit the mortality model
#' @param data_jags.in List contianing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @param Dj_latent Boolean specifying if the variable Dj (occurence of a disturbance) comes from 
#'                  the data or should be estimated
#' @return A rjags object
fit_mortality <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin, Dj_latent){
  
  ## - Write the model
  
  # Loop on all plots
  model.plot.loops <- c("# Loop on all plots
                        for(j in 1:Nplots){", 
                        "
                        }")
  if(Dj_latent){
    model.plot.loops[1] <- paste0(model.plot.loops[1], 
                                  "
                                  D[j] ~ dbern(pD[j])
                                  logit(pD[j]) = a0 + a1*DA[j]
                                  ")
  }
  
  # Loop on all trees
  model.trees.loop <- "
  # Loop on trees
    for(i in beg[j]:end[j]){
      
      # Probability to observe harvested tree
      h[i] ~ dbern(ph[i])
      ph[i] = (1 - D[j])*phadBM[i] + D[j]*pdDj[i]*phdD[i]
      
      # Probability to observe dead tree
      d[i] ~ dbern(pd[i])
      pd[i] = (1 - D[j])*pdBM[i]*(1 - phadBM[i]) + D[j]*pdDj[i]*(1 - phdD[i])
      
      # Probability to observe alive tree
      a[i] ~ dbern(pa[i])
      pa[i] = (1 - D[j])*(1 - pdBM[i])*(1 - phadBM[i]) + D[j]*(1 - pdDj[i])
      
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - pdBM[i])*(1 - pdD[i])
      
      # Probability to die from background mortality
      logit(pdBM[i]) = b0 + b1*dbh[i] + b2*comp[i] + b3*sgdd[j] + b4*wai[j]
      
      # Probability to die from a disturbance
      logit(pdD[i]) = c0 + c1*DS[j] + c2*dbh[i]*DS[j]
      
      # Probability to be harvested knowing that the tree is dead from a disturbance
      logit(phdD[i]) = d0 + d1*dbh[i] + d2*DS[j] + d3*dbh[i]*DS[j]
      
      # Probability to be harvested in undisturbed plots
      logit(phadBM[i]) = e0 + e1*dbh[i]
    }
  "
  
  
  # Priors
  model.priors <- "
  b0 ~ dnorm(0, 0.01)
  b1 ~ dnorm(-1.5, 0.01)
  b2 ~ dnorm(3, 0.01)
  b3 ~ dnorm(-1.5, 0.01)
  b4 ~ dnorm(1.5, 0.01)
  c0 ~ dnorm(0, 0.01)
  c1 ~ dnorm(3, 0.01)
  c2 ~ dnorm(-1.5, 0.01)
  "
  if(Dj_latent){model.priors <- paste0("
                                       a0 ~ dnorm(0, 0.01)
                                       a1 ~ dnorm(1.5, 0.01)
                                       ", model.priors)}
  
  # Full Model
  model <- paste0("model{
                  ", model.plot.loops[1], 
                  model.trees.loop, model.plot.loops[2], 
                  model.priors, "}")
 
  
  ## - Function to initialize priors
  if(Dj_latent){
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
    param.in <- c("a0", "a1", "b0", "b1", "b2", "b3", 
                  "b4", "c0", "c1", "c2")
  } else {
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
    param.in <- c("b0", "b1", "b2", "b3", 
                  "b4", "c0", "c1", "c2")
  }
  
  
  # Run model
  out <- R2jags::jags(data = data_jags.in$data,
                      param = param.in,
                      inits = list(initjags(),
                                   initjags(),
                                   initjags()),
                      model.file = textConnection(model),
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}