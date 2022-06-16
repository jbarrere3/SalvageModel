data {
  int<lower=1> Ntrees;
  int<lower=1> Nspecies;
  int<lower=1> Nplot;
  int<lower=1,upper=Nplot> plot[Ntrees];
  int<lower=1,upper=Nspecies> sp[Ntrees];
  vector [Ntrees] w;
  vector <lower=1> [Ntrees] time;
  vector <lower=0> [Ntrees] dbh;
  int <lower=0, upper=1> d[Ntrees];
  int <lower=0, upper=1> domi[Ntrees];
  int <lower=0, upper=1> Dstorm[Ntrees];
  int <lower=0, upper=1> Dother[Ntrees];
  int <lower=0, upper=1> Dfire[Ntrees];
}

parameters {
  // Parameters for storm
  real st0 [Nspecies]; 
  real<lower=0.001> st1 [Nspecies]; 
  real<lower=-5, upper=5> st2 [Nspecies]; 
  real st3 [Nspecies]; 
  real<lower=0.001> st4 [Nspecies]; 
  real<lower=-5, upper=5> st5 [Nspecies]; 
  
  // Parameters for other disturbances
  real ot0 [Nspecies]; 
  real<lower=0.001> ot1 [Nspecies]; 
  real<lower=-5, upper=5> ot2 [Nspecies]; 
  real ot3 [Nspecies]; 
  real<lower=0.001> ot4 [Nspecies]; 
  real<lower=-5, upper=5> ot5 [Nspecies]; 
  
  // Parameters for fire
  real fi0 [Nspecies]; 
  real<lower=0.001> fi1 [Nspecies]; 
  real<lower=-5, upper=5> fi2 [Nspecies]; 
  real fi3 [Nspecies]; 
  real<lower=0.001> fi4 [Nspecies]; 
  real<lower=-5, upper=5> fi5 [Nspecies]; 
  
  // Storm, other and fire intensity
  real<lower=0, upper = 1> Istorm [Nplot]; 
  real<lower=0, upper = 1> Iother [Nplot]; 
  real<lower=0, upper = 1> Ifire [Nplot]; 
  
}


transformed parameters {
  // Fitted probabilities
  real<lower=0, upper = 1> pdstorm [Ntrees];
  real<lower=0, upper = 1> pdother [Ntrees];
  real<lower=0, upper = 1> pdfire [Ntrees];
  real<lower=0, upper = 1> pdD [Ntrees];
  
  // Transformation
  for (i in 1:Ntrees) {
    
      // Probability to die from a storm disturbance
      pdstorm[i] = 1 - (1 - (domi[i]*inv_logit(st0[sp[i]] + (st1[sp[i]]*Istorm[plot[i]]*(dbh[i]^st2[sp[i]]))) + 
      (1 - domi[i])*inv_logit(st3[sp[i]] + (st4[sp[i]]*Istorm[plot[i]]*(dbh[i]^st5[sp[i]])))))^time[i];
      
      // Probability to die from an other disturbance
      pdother[i] = 1 - (1 - (domi[i]*inv_logit(ot0[sp[i]] + (ot1[sp[i]]*Iother[plot[i]]*(dbh[i]^ot2[sp[i]]))) + 
      (1 - domi[i])*inv_logit(ot3[sp[i]] + (ot4[sp[i]]*Iother[plot[i]]*(dbh[i]^ot5[sp[i]])))))^time[i];
      
      // Probability to die from a fire disturbance
      pdfire[i] = 1 - (1 - (domi[i]*inv_logit(fi0[sp[i]] + (fi1[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi2[sp[i]]))) + 
      (1 - domi[i])*inv_logit(fi3[sp[i]] + (fi4[sp[i]]*Ifire[plot[i]]*(dbh[i]^fi5[sp[i]])))))^time[i];
      
      // Probability that the tree died from a disturbance
      pdD[i] = 1 - (1 - Dfire[i]*pdfire[i])*(1 - Dstorm[i]*pdstorm[i])*(1 - Dother[i]*pdother[i]);
    }
}


model{

    // Priors at species level
    for(s in 1:Nspecies){
      st0[s] ~ normal(-5, 10);
      st1[s] ~ normal(0, 10);
      st2[s] ~ normal(0, 1); 
      st3[s] ~ normal(-5, 10);
      st4[s] ~ normal(0, 10); 
      st5[s] ~ normal(0, 1); 
      ot0[s] ~ normal(-5, 10);
      ot1[s] ~ normal(0, 10); 
      ot2[s] ~ normal(0, 1); 
      ot3[s] ~ normal(-5, 10);
      ot4[s] ~ normal(0, 10); 
      ot5[s] ~ normal(0, 1); 
      fi0[s] ~ normal(-5, 10);
      fi1[s] ~ normal(0, 10); 
      fi2[s] ~ normal(0, 1); 
      fi3[s] ~ normal(-5, 10);
      fi4[s] ~ normal(0, 10); 
      fi5[s] ~ normal(0, 1); 
    }
    
    // Disturbance intensity at plot level
    for(k in 1:Nplot){
      Ifire[k] ~ beta(0.66, 0.36); 
      Istorm[k] ~ beta(0.65, 2.66); 
      Iother[k] ~ beta(0.48, 1.77); 
    }
    
    // Likelihood
    for (i in 1:Ntrees) {
      target += w[i]*bernoulli_logit_lpmf(d[i] | inv_logit(pdD[i]));
    }
}
