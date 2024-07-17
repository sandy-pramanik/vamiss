
data {
  
  int nCause;
  int nAlgo;
  int aj[nAlgo,nCause];
  
  simplex[nCause] p_uncalib_obs;
  vector<lower=0>[nCause] Mmatprior_asDirich[nAlgo,nCause];
  real<lower=0> pss;
  
}

parameters {
  
  simplex[nCause] Mmat[nAlgo,nCause];
  simplex[nCause] p_calib;
  // real<lower=0> pss;
  
}

transformed parameters {
  
  matrix[nCause,nCause] Mmat_byalgo[nAlgo];
  simplex[nCause] q[nAlgo];
  vector[nAlgo] loglik;
  
  for(k in 1:nAlgo){
    
    for(i in 1:nCause){
      
      Mmat_byalgo[k][i,] = Mmat[k,i]';
      
    }
  
    q[k] = (Mmat_byalgo[k]')*p_calib;
  
    loglik[k] = multinomial_lpmf(aj[k,] | q[k]);
  
  }
  
}

model {
  
  // target += beta_lpdf(1/(pss+1) | 0.5, 0.5) - 2*log(pss+1);
  target += dirichlet_lpdf(p_calib | 1 + nCause*pss*p_uncalib_obs);
  // target += dirichlet_lpdf(p_calib | nCause*pss*p_uncalib_obs);
  // target += dirichlet_lpdf(p_calib | rep_vector(1, nCause));
  
  for(k in 1:nAlgo){
  
    for(i in 1:nCause){
      
      target += dirichlet_lpdf(Mmat[k,i] | Mmatprior_asDirich[k,i]);
      
    }
    
  }
  
  target += loglik;
  
}
