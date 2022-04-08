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





#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
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





#' Fit the mortality model with disturbance event given in the data
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_sub <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
      state[i] ~ dcat(proba[i, 1:3])
      # Probability to observe harvested tree
      proba[i, 1] = pdDj[i]*phdD[i]
      # Probability to observe dead tree
      proba[i, 2] = pdDj[i]*(1 - phdD[i])
      # Probability to observe alive tree
      proba[i, 3] = (1 - pdDj[i])
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
      # Probability to die from a storm disturbance
      logit(pstorm[i]) = c0[sp[i]] + (c1[sp[i]]*Istorm[plot[i]]*(dbh[i]^c2[sp[i]])) 
      pdstorm[i] = 1 - (1 - pstorm[i])^5
      # Probability to die from an other disturbance
      logit(pother[i]) = c3[sp[i]] + (c4[sp[i]]*Iother[plot[i]]*(dbh[i]^c5[sp[i]])) 
      pdother[i] = 1 - (1 - pother[i])^time[i]
    }
    
    ## - Priors
    # Priors at species level
    for(s in 1:Nspecies){
      c0[s] ~ dnorm(0, 0.5)
      c1[s] ~ dnorm(6, 1)
      c2[s] ~ dnorm(0, 0.5)
      c3[s] ~ dnorm(0, 0.5)
      c4[s] ~ dnorm(6, 1)
      c5[s] ~ dnorm(0, 0.5)
    }
    # Disturbance intensity at plot level
    #for(k in 1:Nplot){
     # Istorm[k] ~ dunif(0, 1)
      #Iother[k] ~ dunif(0, 1)
    #}
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      Istorm[k] ~ dbeta(0.25*5, (1 - 0.25)*5) T(0.001,0.999)
      Iother[k] ~ dbeta(0.25*5, (1 - 0.25)*5) T(0.001,0.999)
    }
    
    
    
  }"
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = paste0("c", c(0:5)),
                      model.file = textConnection(mortality_model_D),
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}



#' Fit the mortality model with disturbance event given in the data for several countries
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_full_sub <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      state[i] ~ dcat(proba[i, 1:3])
      
      # Probability to observe harvested tree
      proba[i, 1] = pdDj[i]*phdD[i]
      
      # Probability to observe dead tree
      proba[i, 2] = pdDj[i]*(1 - phdD[i])
      
      # Probability to observe alive tree
      proba[i, 3] = (1 - pdDj[i])
      
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
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
      c0[s] ~ dnorm(-2, 1)
      c1[s] ~ dnorm(3, 1)
      c2[s] ~ dnorm(0, 0.5)
      c3[s] ~ dnorm(-2, 1)
      c4[s] ~ dnorm(3, 1)
      c5[s] ~ dnorm(0, 0.5)
      c6[s] ~ dnorm(-2, 1)
      c7[s] ~ dnorm(3, 1)
      c8[s] ~ dnorm(0, 0.5)
    }
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      Ifire[k] ~ dbeta(0.66, 0.36) T(0.001,0.999)
      Istorm[k] ~ dbeta(0.65, 2.66) T(0.001,0.999)
      Iother[k] ~ dbeta(0.48, 1.77) T(0.001,0.999)
    }
    
    
    
  }"
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = paste0("c", c(0:8)),
                      model.file = textConnection(mortality_model_D),
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}




#' Fit the mortality model with disturbance event given in the data for several countries
#' @param data_jags.in List containing all the inputs of the model
#' @param n.chains numeric: Number of MCMC Markov chains
#' @param n.iter numeric: Number of iterations
#' @param n.burn numeric: Burn-in
#' @param n.thin numeric: Thinning rate
#' @return A rjags object
fit_mortality_full_unif_sub <- function(data_jags.in, n.chains, n.iter, n.burn, n.thin){
  
  ## - Write the model
  mortality_model_D <- 
    "model{
    for (i in 1:Ntrees) {
    
      state[i] ~ dcat(proba[i, 1:3])
      
      # Probability to observe harvested tree
      proba[i, 1] = pdDj[i]*phdD[i]
      
      # Probability to observe dead tree
      proba[i, 2] = pdDj[i]*(1 - phdD[i])
      
      # Probability to observe alive tree
      proba[i, 3] = (1 - pdDj[i])
      
      # Probability to die in disturbed plots
      pdDj[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i])
      
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
      c0[s] ~ dnorm(-2, 1)
      c1[s] ~ dnorm(3, 1)
      c2[s] ~ dnorm(0, 0.5)
      c3[s] ~ dnorm(-2, 1)
      c4[s] ~ dnorm(3, 1)
      c5[s] ~ dnorm(0, 0.5)
      c6[s] ~ dnorm(-2, 1)
      c7[s] ~ dnorm(3, 1)
      c8[s] ~ dnorm(0, 0.5)
    }
    
    # Disturbance intensity at plot level
    for(k in 1:Nplot){
      Ifire[k] ~ dunif(0, 1)
      Istorm[k] ~ dunif(0, 1)
      Iother[k] ~ dunif(0, 1)
    }
    
    
    
  }"
  
  
  
  ## - Fit the model
  out <- R2jags::jags(data = data_jags.in,
                      param = paste0("c", c(0:8)),
                      model.file = textConnection(mortality_model_D),
                      n.chains = n.chains,
                      n.iter = n.iter,
                      n.burnin = n.burn,
                      n.thin = n.thin,
                      DIC = TRUE, 
                      progress.bar = "text")
  
  return(out)
}