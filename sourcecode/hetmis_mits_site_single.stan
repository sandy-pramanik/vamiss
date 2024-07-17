
data {
  
  int nCause;
  int nSite;
  int errormat[nSite*nCause,nCause];
  int nObs;
  int idnonzero[nObs];
  int fpindex[nCause,nCause-1];
  
  vector<lower=0>[2] shape_pull;
  vector<lower=0>[2] shape_hetsens;
  vector<lower=0>[2] shape_hetrfp;
  
  real<lower=0> shapeBase;
  vector<lower=0>[2] sensprior_shape;
  vector<lower=0>[nCause] Dirprior_shape;
  
}

parameters {
  
  vector<lower=0,upper=1>[nCause] a_i;
  simplex[nCause] alpha;
  
  real<lower=0> pss_pull;
  vector<lower=0,upper=1>[nCause] sens_i;
  simplex[nCause-1] q_i[nCause];
  
  real<lower=0> pss_hetsens;
  vector<lower=0,upper=1>[nSite*nCause] sens_si;
  
  real<lower=0> pss_hetrfp;
  simplex[nCause-1] q_si[nSite*nCause];
  
}

transformed parameters {
  
  simplex[nCause] mu_si[nSite*nCause];
  vector[nObs] loglik;
  
  for(s in 1:nSite){
    
    for(i in 1:nCause){
      
      mu_si[(s-1)*nCause + i,i] = sens_si[(s-1)*nCause + i];
      mu_si[(s-1)*nCause + i,fpindex[i,]] = (1-sens_si[(s-1)*nCause + i])*q_si[(s-1)*nCause + i];
    
    }
    
  }
  
  for(l in 1:nObs){
    
    loglik[l] = multinomial_lpmf(errormat[idnonzero[l],] | mu_si[idnonzero[l]]);
    
  }
  
}

model {
  
  target += beta_lpdf(a_i | sensprior_shape[1], sensprior_shape[2]);
  target += dirichlet_lpdf(alpha | Dirprior_shape);
  
  target += beta_lpdf(1/(pss_pull+1) | shape_pull[1], shape_pull[2]) - 2*log(pss_pull+1);
  target += beta_lpdf(1/(pss_hetsens+1) | shape_hetsens[1], shape_hetsens[2]) - 2*log(pss_hetsens+1);
  target += beta_lpdf(1/(pss_hetrfp+1) | shape_hetsens[1], shape_hetsens[2]) - 2*log(pss_hetrfp+1);
  
  for(i in 1:nCause){
    
    target += beta_lpdf(sens_i[i] | shapeBase + 2*pss_pull*(1 - (1-a_i[i])*(1-alpha[i])), 
                                    shapeBase + 2*pss_pull*(1-a_i[i])*(1-alpha[i]));
    target += dirichlet_lpdf(q_i[i] | shapeBase + (nCause-1)*pss_pull*(alpha[fpindex[i,]]/(1-alpha[i])));
  
    for(s in 1:nSite){
      
      target += beta_lpdf(sens_si[(s-1)*nCause + i] | shapeBase + 2*pss_hetsens*sens_i[i], 
                                                      shapeBase + 2*pss_hetsens*(1-sens_i[i]));
      target += dirichlet_lpdf(q_si[(s-1)*nCause + i] | shapeBase + (nCause-1)*pss_hetrfp*q_i[i]);
      
    }
    
  }
  
  target += loglik;
  
}

generated quantities {

  vector<lower=0,upper=1>[nCause] sens_i_pull = 1 - (1-a_i) .* (1-alpha);
  simplex[nCause-1] q_i_pull[nCause];
  simplex[nCause] mu_i_pull[nCause];
  
  simplex[nCause] mu_i[nCause];

  vector<lower=0,upper=1>[nCause] sens_si_pred;
  simplex[nCause-1] q_si_pred[nCause];
  simplex[nCause] mu_si_pred[nCause];
    
  for(i in 1:nCause){
    
    q_i_pull[i] = alpha[fpindex[i,]]/(1-alpha[i]);
    mu_i_pull[i,i] = sens_i_pull[i];
    mu_i_pull[i,fpindex[i,]] = (1-sens_i_pull[i])*q_i_pull[i];
  
    mu_i[i,i] = sens_i[i];
    mu_i[i,fpindex[i,]] = (1-sens_i[i])*q_i[i];
    
    sens_si_pred[i] = beta_rng(shapeBase + 2*pss_hetsens*sens_i[i], 
                               shapeBase + 2*pss_hetsens*(1-sens_i[i]));
    q_si_pred[i] = dirichlet_rng(shapeBase + (nCause-1)*pss_hetrfp*q_i[i]);
    mu_si_pred[i,i] = sens_si_pred[i];
    mu_si_pred[i,fpindex[i,]] = (1-sens_si_pred[i])*q_si_pred[i];

  }
  
}
