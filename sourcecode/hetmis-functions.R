

# loading libraries ----
library(reshape2)
library(tidyverse)
library(doParallel)
library(rstan)
rstan_options(auto_write = TRUE)


# stan object of the model ----
hetmis_mits_single.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                                    'hetmis_mits_single.stan'),
                                                   auto_write = rstan_options("auto_write"=TRUE))
hetmis_mits_single_pss.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                                        'hetmis_mits_single_pss.stan'),
                                                       auto_write = rstan_options("auto_write"=TRUE))
hetmis_mits_site_sens_single.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                                              'hetmis_mits_site_sens_single.stan'),
                                                             auto_write = rstan_options("auto_write"=TRUE))
hetmis_mits_site_sens_single_pss.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                                                  'hetmis_mits_site_sens_single_pss.stan'),
                                                                 auto_write = rstan_options("auto_write"=TRUE))
hetmis_mits_site_single.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                                         'hetmis_mits_site_single.stan'),
                                                        auto_write = rstan_options("auto_write"=TRUE))
hetmis_mits_site_single_pss.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                                             'hetmis_mits_site_single_pss.stan'),
                                                            auto_write = rstan_options("auto_write"=TRUE))
seqcalib_mmat.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                          'seqcalib_mmat.stan'),
                                         auto_write = rstan_options("auto_write"=TRUE))
seqcalib.stan_object = rstan::stan_model(file = file.path(sourcecode.path,
                                                          'seqcalib.stan'),
                                         auto_write = rstan_options("auto_write"=TRUE))


# helper function ----
## for combining outputs from replicated study
mysimcombine_insample = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = list.combined[[1]]
  
  for(k in 2:length.list.combined){
    
    list.out$uncalib = rbind(list.out$uncalib, list.combined[[k]]$uncalib)
    
    list.out$calib = abind::abind(list.out$calib,
                                  list.combined[[k]]$calib,
                                  along = 4)
    
  }
  
  return(list.out)
  
}

mysimcombine_outsample = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = list.combined[[1]]
  
  for(k in 2:length.list.combined){
    
    list.out$Mmat = abind::abind(list.out$Mmat,
                                 list.combined[[k]]$Mmat,
                                 along = 5)
    
    list.out$p_calib = abind::abind(list.out$p_calib,
                                    list.combined[[k]]$p_calib,
                                    along = 4)
    
  }
  
  return(list.out)
  
}

mysimcombine = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = list.combined[[1]]
  
  for(k in 2:length.list.combined){
    
    list.out$Mmat.obs = abind::abind(list.out$Mmat.obs,
                                     list.combined[[k]]$Mmat.obs,
                                     along = 4)
    
    list.out$Mmat.insample = abind::abind(list.out$Mmat.insample,
                                          list.combined[[k]]$Mmat.insample,
                                          along = 6)
    
    list.out$ic = abind::abind(list.out$ic,
                               list.combined[[k]]$ic,
                               along = 3)
    
    list.out$csmf.insample_uncalib = abind::abind(list.out$csmf.insample_uncalib,
                                                  list.combined[[k]]$csmf.insample_uncalib,
                                                  along = 3)
    
    list.out$csmf.insample_calib = abind::abind(list.out$csmf.insample_calib,
                                                list.combined[[k]]$csmf.insample_calib,
                                                along = 5)
    
    list.out$Mmat.outsample = abind::abind(list.out$Mmat.outsample,
                                           list.combined[[k]]$Mmat.outsample,
                                           along = 6)
    
    list.out$csmf.outsample_calib = abind::abind(list.out$csmf.outsample_calib,
                                                 list.combined[[k]]$csmf.outsample_calib,
                                                 along = 5)
    
  }
  
  return(list.out)
  
}

mychampscombine_outsample = function(...){
  
  list.combined = list(...)
  length.list.combined = length(list.combined)
  list.out = list.combined[[1]]
  
  for(k in 2:length.list.combined){
    
    list.out = abind::abind(list.out,
                            list.combined[[k]],
                            along = 5)
    
  }
  
  return(list.out)
  
}


# insample ----
hetmis <- function(va_labeled = NULL, gold_standard = NULL, errormat_by_site = NULL, 
                   model.choice = c('hom', 'het_part', 'het',
                                    'separate_empirical', 'pooled_empirical'),
                   cause.type = "single",
                   pss.method = 'learn', shapeBase = .5,
                   shape_pull = c(.5, .5),
                   shape_hetsens = c(.5, .5), shape_hetrfp = c(.5, .5),
                   pss_fixed = NULL,
                   sensprior_shape = c(1,1), Dirprior_shape = NULL,
                   nMCMC = 5000, nBurn = 2000, nThin = 1, adapt_delta_stan = .9,
                   sourcecode.path = file.path(getwd(), "sourcecode"),
                   seed = 1, verbose = T, saveoutput = T,
                   output_dir = NULL, output_filename = NULL){
  
  input.list = as.list(environment())
  
  if(verbose){
    
    refresh.stan = max((nBurn + nMCMC*nThin)/10, 1)
    
  }else{
    
    refresh.stan = 0
    
  }
  
  ## datalist for modeling ----
  if(any(is.null(va_labeled), is.null(gold_standard)) && is.null(errormat_by_site)){
    
    stop("Need to provide either both 'va_labeled' and 'gold_standard', or 'errormat_by_site'")
    
  }else if((!is.null(va_labeled)) && (!is.null(gold_standard))){
    
    # have column names check
    if(any(is.null(colnames(va_labeled)), is.null(colnames(gold_standard)))) stop("'va_labeled' and 'gold_standard' should have column names and they should be same")
    
    # same column names check
    if(!identical(colnames(va_labeled), colnames(gold_standard))) stop("'va_labeled' and 'gold_standard' should have the same column names")
    
    # same row names check if present. making sure rows in 'va_labeled' and 'gold_standard' from the same patient
    if((!is.null(rownames(va_labeled))) && (!is.null(rownames(gold_standard))) && (!identical(rownames(va_labeled), rownames(gold_standard)))){
      
      stop("rownames of 'va_labeled' and 'gold_standard' (patient id's) should be the same")
      
    }else if(any(is.null(rownames(va_labeled)), is.null(rownames(gold_standard)))){
      
      print("Assuming the rows in 'va_labeled' and 'gold_standard' correspond to the same patient")
      
    }
    
    # same site names check
    if(!identical(va_labeled$site, gold_standard$site)) stop("'site' in 'va_labeled' and 'gold_standard' do not match")
    
    # all unique causes
    causes = colnames(gold_standard)[-ncol(gold_standard)]
    
    # all unique sites
    ## gold_standard
    if(!any(is.character(gold_standard$site), is.factor(gold_standard$site))){
      
      stop("'site' in 'gold_standard' must be either character or factor")
      
    }else if(is.character(gold_standard$site)){
      
      gold_standard$site = factor(x = gold_standard$site, 
                                  levels = sort(unique(gold_standard$site)))
      
    }
    
    ## va_labeled
    if(!any(is.character(va_labeled$site), is.factor(va_labeled$site))){
      
      stop("'site' in 'va_labeled' must be either character or factor")
      
    }else if(is.character(va_labeled$site)){
      
      va_labeled$site = factor(x = va_labeled$site, 
                               levels = sort(unique(va_labeled$site)))
      
    }
    
    sites = levels(gold_standard$site)
    
    # error matrix by site
    if(cause.type=='single'){
      
      errormat_by_site = Reduce('rbind',
                                lapply(sites, function(site_s){
                                  
                                  if(verbose) print(site_s)
                                  
                                  mits_sub = gold_standard %>% dplyr::filter(site==site_s)
                                  champs_sub = va_labeled %>% dplyr::filter(site==site_s)
                                  
                                  error_site_s = as.matrix(t(mits_sub[,causes]))%*%as.matrix(champs_sub[,causes])
                                  mits.total_site_s = rowSums(error_site_s)
                                  
                                  colnames(error_site_s) <- causes
                                  rownames(error_site_s) <- paste0(site_s, '_', causes, '(', mits.total_site_s, ')')
                                  
                                  error_site_s
                                  
                                }))
      
      mits.total_by_site = rowSums(errormat_by_site)
      
    }
    
    nCause = length(causes)
    nSite = length(sites)
    
  }else if(!is.null(errormat_by_site)){
    
    if(is.null(names(errormat_by_site))) stop("The components of 'errormat_by_site' should have the respective sites as names")
    sites = names(errormat_by_site)
    
    if(cause.type=='single'){
      
      errormat_by_site_new = do.call('rbind',
                                     lapply(sites,
                                            FUN = function(site_s){
                                              
                                              if(is.null(colnames(errormat_by_site[[site_s]]))) stop(paste0("No column names in the errormat for ", site_s,
                                                                                                            ". Should have VA causes as column names."))
                                              if(is.null(rownames(errormat_by_site[[site_s]]))) stop(paste0("No row names in the errormat for ", site_s,
                                                                                                            ". Should have MITS causes as column names."))
                                              if(!identical(rownames(errormat_by_site[[site_s]]), colnames(errormat_by_site[[site_s]]))) stop(paste0("Row and column names in the errormat for ", site_s,
                                                                                                                                                     " is either not the same or not in the same order"))
                                              mits.total_site_s = rowSums(errormat_by_site[[site_s]])
                                              
                                              rownames(errormat_by_site[[site_s]]) <- paste0(site_s, '_', rownames(errormat_by_site[[site_s]]),
                                                                                             '(', mits.total_site_s, ')')
                                              
                                              errormat_by_site[[site_s]]
                                              
                                            }))
      errormat_by_site = errormat_by_site_new
      rm(errormat_by_site_new)
      
      # all unique causes
      causes = colnames(errormat_by_site)
      
      mits.total_by_site = rowSums(errormat_by_site)
      nCause = length(causes)
      nSite = length(sites)
      
    }
    
  }
  
  mits_by_site = rep(causes, nSite)
  sites_by_site = rep(sites, each = nCause)
  
  mitsid_by_site = match(mits_by_site, causes)
  siteid_by_site = match(sites_by_site, sites)
  
  # identifying sites with at least one mits case
  if(cause.type=='single'){
    
    # error proportions
    errorpropmat_by_site = errormat_by_site/mits.total_by_site
    
    idnonNA_by_site = unname(which(mits.total_by_site!=0))
    
    errormat_by_site_nonNA = errormat_by_site[idnonNA_by_site,]
    mits.total_by_site_nonNA = mits.total_by_site[idnonNA_by_site]
    errorpropmat_by_site_nonNA = errorpropmat_by_site[idnonNA_by_site,]
    
    mits_by_site_nonNA = mits_by_site[idnonNA_by_site]
    sites_by_site_nonNA = sites_by_site[idnonNA_by_site]
    
    mitsid_by_site_nonNA = mitsid_by_site[idnonNA_by_site]
    siteid_by_site_nonNA = siteid_by_site[idnonNA_by_site]
    
  }
  
  # datalist
  datalist = list('causes' = causes, 'sites' = sites,
                  'nCause' = nCause, 'nSite' = nSite,
                  'errormat_by_site' = errormat_by_site,
                  'mits.total_by_site' = mits.total_by_site,
                  'errorpropmat_by_site' = errorpropmat_by_site,
                  'mits_by_site' = mits_by_site, 'sites_by_site' = sites_by_site,
                  'mitsid_by_site' = mitsid_by_site, 'siteid_by_site' = siteid_by_site,
                  'idnonNA_by_site' = idnonNA_by_site,
                  'errormat_by_site_nonNA' = errormat_by_site_nonNA,
                  'mits.total_by_site_nonNA' = mits.total_by_site_nonNA,
                  'errorpropmat_by_site_nonNA' = errorpropmat_by_site_nonNA,
                  'mits_by_site_nonNA' = mits_by_site_nonNA, 'sites_by_site_nonNA' = sites_by_site_nonNA,
                  'mitsid_by_site_nonNA' = mitsid_by_site_nonNA, 'siteid_by_site_nonNA' = siteid_by_site_nonNA,
                  'fpindex' = do.call('rbind',
                                      lapply(1:nCause,
                                             FUN = function(i){
                                               
                                               (1:nCause)[-i]
                                               
                                             })),
                  'va_labeled' = va_labeled, 'gold_standard' = gold_standard,
                  'siteid_labeled' = match(va_labeled$site, sites))
  
  # default Dirichlet shape parameters
  if(is.null(Dirprior_shape)) Dirprior_shape = rep(1, datalist$nCause)
  
  if(model.choice=='hom'){
    
    ## mits ----
    if(pss.method=='learn'){
      
      ### learn pss ----
      #### stan fit ----
      stanfit = rstan::sampling(hetmis_mits_single.stan_object,
                                data = list('nCause' = datalist$nCause,
                                            'nSite' = datalist$nSite,
                                            'errormat' = datalist$errormat_by_site,
                                            'nObs' = length(datalist$idnonNA_by_site),
                                            'idnonzero' = datalist$idnonNA_by_site,
                                            'fpindex' = datalist$fpindex,
                                            'shape_pull' = shape_pull,
                                            'shapeBase' = shapeBase,
                                            'sensprior_shape' = sensprior_shape,
                                            'Dirprior_shape' = Dirprior_shape),
                                chains = 1, 
                                iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                control = list('adapt_delta' = adapt_delta_stan),
                                seed = seed,
                                refresh = refresh.stan)
      
      #### MCMC output ----
      MCMCout = rstan::extract(stanfit)
      
      #### max Rhat ----
      max_Rhat = max(apply(X = MCMCout$a_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$alpha, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     rstan::Rhat(c(MCMCout$pss_pull)),
                     apply(X = MCMCout$sens_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$mu_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }))
      
      #### min bulk ESS ----
      min_ess_bulk = min(apply(X = MCMCout$a_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$alpha, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         rstan::ess_bulk(c(MCMCout$pss_pull)),
                         apply(X = MCMCout$sens_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$mu_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }))/nMCMC
      
      # pss types and all pss in one matrix
      MCMCout$pss = cbind.data.frame('pull' = MCMCout$pss_pull)
      
    }else if(pss.method=='fixed'){
      
      # setting default pss_fixed
      if(is.null(pss_fixed)) pss_fixed = list('pss_pull' = 2)
      
      ### fixed pss ----
      #### stan fit ----
      stanfit = rstan::sampling(hetmis_mits_single_pss.stan_object,
                                data = list('nCause' = datalist$nCause,
                                            'nSite' = datalist$nSite,
                                            'errormat' = datalist$errormat_by_site,
                                            'nObs' = length(datalist$idnonNA_by_site),
                                            'idnonzero' = datalist$idnonNA_by_site,
                                            'fpindex' = datalist$fpindex,
                                            'pss_pull' = pss_fixed$pss_pull,
                                            'shapeBase' = shapeBase,
                                            'sensprior_shape' = sensprior_shape,
                                            'Dirprior_shape' = Dirprior_shape),
                                chains = 1, 
                                iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                control = list('adapt_delta' = adapt_delta_stan),
                                seed = seed,
                                refresh = refresh.stan)
      
      #### MCMC output ----
      MCMCout = rstan::extract(stanfit)
      
      #### max Rhat ----
      max_Rhat = max(apply(X = MCMCout$a_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$alpha, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$sens_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$mu_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }))
      
      #### min bulk ESS ----
      min_ess_bulk = min(apply(X = MCMCout$a_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$alpha, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$sens_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$mu_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }))/nMCMC
      
      # pss types and all pss in one matrix
      MCMCout$pss = cbind.data.frame('pull' = pss_fixed$pss_pull)
      
    }
    
    m_insample = array(dim = c(nMCMC, datalist$nCause*datalist$nSite,
                               datalist$nCause),
                       dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                       colnames(datalist$errormat_by_site)))
    m_insample[,datalist$idnonNA_by_site,] = MCMCout$mu_si[,datalist$idnonNA_by_site,]
    MCMCout$m_insample = m_insample
    
    m_pred = MCMCout$mu_si_pred
    dimnames(m_pred)[[2]] = datalist$causes
    dimnames(m_pred)[[3]] = datalist$causes
    MCMCout$m_pred = m_pred
    
    MCMCout$est.type = 'homogen'
    
  }else if(model.choice=='het_part'){
    
    ## mits+site_sens ----
    if(pss.method=='learn'){
      
      ### learn pss ----
      #### stan fit ----
      stanfit = rstan::sampling(hetmis_mits_site_sens_single.stan_object,
                                data = list('nCause' = datalist$nCause,
                                            'nSite' = datalist$nSite,
                                            'errormat' = datalist$errormat_by_site,
                                            'nObs' = length(datalist$idnonNA_by_site),
                                            'idnonzero' = datalist$idnonNA_by_site,
                                            'fpindex' = datalist$fpindex,
                                            'shape_pull' = shape_pull,
                                            'shape_hetsens' = shape_hetsens,
                                            'shape_hetrfp' = shape_hetrfp,
                                            'shapeBase' = shapeBase,
                                            'sensprior_shape' = sensprior_shape,
                                            'Dirprior_shape' = Dirprior_shape),
                                chains = 1, 
                                iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                control = list('adapt_delta' = adapt_delta_stan),
                                seed = seed,
                                refresh = refresh.stan)
      
      #### MCMC output ----
      MCMCout = rstan::extract(stanfit)
      
      #### max Rhat ----
      max_Rhat = max(apply(X = MCMCout$a_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$alpha, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     rstan::Rhat(c(MCMCout$pss_pull)),
                     apply(X = MCMCout$sens_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     rstan::Rhat(c(MCMCout$pss_hetsens)),
                     apply(X = MCMCout$sens_si, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$mu_si, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }))
      
      #### min bulk ESS ----
      min_ess_bulk = min(apply(X = MCMCout$a_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$alpha, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         rstan::ess_bulk(c(MCMCout$pss_pull)),
                         apply(X = MCMCout$sens_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         rstan::ess_bulk(c(MCMCout$pss_hetsens)),
                         apply(X = MCMCout$sens_si, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$mu_si, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }))/nMCMC
      
      # pss types and all pss in one data frame
      MCMCout$pss = cbind.data.frame('pull' = MCMCout$pss_pull,
                                       'hetsens' = MCMCout$pss_hetsens)
      
    }else if(pss.method=='fixed'){
      
      # setting default pss_fixed
      if(is.null(pss_fixed)) pss_fixed = list('pss_pull' = 2,
                                              'pss_hetsens' = 2)
      
      ### fixed pss ----
      #### stan fit ----
      stanfit = rstan::sampling(hetmis_mits_site_sens_single_pss.stan_object,
                                data = list('nCause' = datalist$nCause,
                                            'nSite' = datalist$nSite,
                                            'errormat' = datalist$errormat_by_site,
                                            'nObs' = length(datalist$idnonNA_by_site),
                                            'idnonzero' = datalist$idnonNA_by_site,
                                            'fpindex' = datalist$fpindex,
                                            'pss_pull' = pss_fixed$pss_pull,
                                            'pss_hetsens' = pss_fixed$pss_hetsens,
                                            'shape_pull' = shape_pull,
                                            'shape_hetsens' = shape_hetsens,
                                            'shapeBase' = shapeBase,
                                            'sensprior_shape' = sensprior_shape,
                                            'Dirprior_shape' = Dirprior_shape),
                                chains = 1, 
                                iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                control = list('adapt_delta' = adapt_delta_stan),
                                seed = seed,
                                refresh = refresh.stan)
      
      #### MCMC output ----
      MCMCout = rstan::extract(stanfit)
      
      #### max Rhat ----
      max_Rhat = max(apply(X = MCMCout$a_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$alpha, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$sens_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$sens_si, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$mu_si, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }))
      
      #### min bulk ESS ----
      min_ess_bulk = min(apply(X = MCMCout$a_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$alpha, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$sens_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$sens_si, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$mu_si, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }))/nMCMC
      
      # pss types and all pss in one data frame
      MCMCout$pss = cbind.data.frame('pull' = pss_fixed$pss_pull,
                                       'hetsens' = pss_fixed$pss_hetsens)
      
    }
    
    m_insample = array(dim = c(nMCMC, datalist$nCause*datalist$nSite,
                               datalist$nCause),
                       dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                       colnames(datalist$errormat_by_site)))
    m_insample[,datalist$idnonNA_by_site,] = MCMCout$mu_si[,datalist$idnonNA_by_site,]
    MCMCout$m_insample = m_insample
    
    m_pred = MCMCout$mu_si_pred
    dimnames(m_pred)[[2]] = datalist$causes
    dimnames(m_pred)[[3]] = datalist$causes
    MCMCout$m_pred = m_pred
    
    MCMCout$est.type = 'heterogen'
    
  }else if(model.choice=='het'){
    
    ## mits+site ----
    if(pss.method=='learn'){
      
      ### learn pss ----
      #### stan fit ----
      stanfit = rstan::sampling(hetmis_mits_site_single.stan_object,
                                data = list('nCause' = datalist$nCause,
                                            'nSite' = datalist$nSite,
                                            'errormat' = datalist$errormat_by_site,
                                            'nObs' = length(datalist$idnonNA_by_site),
                                            'idnonzero' = datalist$idnonNA_by_site,
                                            'fpindex' = datalist$fpindex,
                                            'shape_pull' = shape_pull,
                                            'shape_hetsens' = shape_hetsens,
                                            'shape_hetrfp' = shape_hetrfp,
                                            'shapeBase' = shapeBase,
                                            'sensprior_shape' = sensprior_shape,
                                            'Dirprior_shape' = Dirprior_shape),
                                chains = 1, 
                                iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                control = list('adapt_delta' = adapt_delta_stan),
                                seed = seed,
                                refresh = refresh.stan)
      
      #### MCMC output ----
      MCMCout = rstan::extract(stanfit)
      
      #### max Rhat ----
      max_Rhat = max(apply(X = MCMCout$a_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$alpha, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     rstan::Rhat(c(MCMCout$pss_pull)),
                     apply(X = MCMCout$sens_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     rstan::Rhat(c(MCMCout$pss_hetsens)),
                     apply(X = MCMCout$sens_si, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     rstan::Rhat(c(MCMCout$pss_hetrfp)),
                     apply(X = MCMCout$q_si, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$mu_si, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }))
      
      #### min bulk ESS ----
      min_ess_bulk = min(apply(X = MCMCout$a_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$alpha, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         rstan::ess_bulk(c(MCMCout$pss_pull)),
                         apply(X = MCMCout$sens_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         rstan::ess_bulk(c(MCMCout$pss_hetsens)),
                         apply(X = MCMCout$sens_si, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         rstan::ess_bulk(c(MCMCout$pss_hetrfp)),
                         apply(X = MCMCout$q_si, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$mu_si, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }))/nMCMC
      
      # pss types and all pss in one data frame
      MCMCout$pss = cbind.data.frame('pull' = MCMCout$pss_pull,
                                       'hetsens' = MCMCout$pss_hetsens,
                                       'hetrfp' = MCMCout$pss_hetrfp)
      
    }else if(pss.method=='fixed'){
      
      # setting default pss_fixed
      if(is.null(pss_fixed)) pss_fixed = list('pss_pull' = 2,
                                              'pss_hetsens' = 2,
                                              'pss_hetrfp' = 2)
      
      ### fixed pss ----
      #### stan fit ----
      stanfit = rstan::sampling(hetmis_mits_site_single_pss.stan_object,
                                data = list('nCause' = datalist$nCause,
                                            'nSite' = datalist$nSite,
                                            'errormat' = datalist$errormat_by_site,
                                            'nObs' = length(datalist$idnonNA_by_site),
                                            'idnonzero' = datalist$idnonNA_by_site,
                                            'fpindex' = datalist$fpindex,
                                            'pss_pull' = pss_fixed$pss_pull,
                                            'pss_hetsens' = pss_fixed$pss_hetsens,
                                            'pss_hetrfp' = pss_fixed$pss_hetrfp,
                                            'shape_pull' = shape_pull,
                                            'shape_hetsens' = shape_hetsens,
                                            'shape_hetrfp' = shape_hetrfp,
                                            'shapeBase' = shapeBase,
                                            'sensprior_shape' = sensprior_shape,
                                            'Dirprior_shape' = Dirprior_shape),
                                chains = 1, 
                                iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                                control = list('adapt_delta' = adapt_delta_stan),
                                seed = seed,
                                refresh = refresh.stan)
      
      #### MCMC output ----
      MCMCout = rstan::extract(stanfit)
      
      #### max Rhat ----
      max_Rhat = max(apply(X = MCMCout$a_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$alpha, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$sens_i, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_i, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$sens_si, 2,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$q_si, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }),
                     apply(X = MCMCout$mu_si, 2:3,
                           FUN = function(v){
                             
                             rstan::Rhat(v)
                             
                           }))
      
      #### min bulk ESS ----
      min_ess_bulk = min(apply(X = MCMCout$a_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$alpha, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$sens_i, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_i, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$sens_si, 2,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$q_si, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }),
                         apply(X = MCMCout$mu_si, 2:3,
                               FUN = function(v){
                                 
                                 rstan::ess_bulk(v)
                                 
                               }))/nMCMC
      
      # pss types and all pss in one data frame
      MCMCout$pss = cbind.data.frame('pull' = pss_fixed$pss_pull,
                                       'hetsens' = pss_fixed$pss_hetsens,
                                       'hetrfp' = pss_fixed$pss_hetrfp)
      
    }
    
    m_insample = array(dim = c(nMCMC, datalist$nCause*datalist$nSite,
                               datalist$nCause),
                       dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                       colnames(datalist$errormat_by_site)))
    m_insample[,datalist$idnonNA_by_site,] = MCMCout$mu_si[,datalist$idnonNA_by_site,]
    MCMCout$m_insample = m_insample
    
    m_pred = MCMCout$mu_si_pred
    dimnames(m_pred)[[2]] = datalist$causes
    dimnames(m_pred)[[3]] = datalist$causes
    MCMCout$m_pred = m_pred
    
    MCMCout$est.type = 'heterogen'
    
  }else if(model.choice=='separate_empirical'){
    
    ## separate empirical ----
    if(cause.type=='single'){
      
      ### single ----
      mu_si = array(dim = c(3, datalist$nCause*datalist$nSite,
                            datalist$nCause),
                    dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                    colnames(datalist$errormat_by_site)))
      q_si = array(dim = c(3, datalist$nCause*datalist$nSite,
                           datalist$nCause-1),
                   dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                   NULL))
      for(l in datalist$idnonNA_by_site){
        
        mu_si[,l,] =  t(DescTools::MultinomCI(x = datalist$errormat_by_site[l,],
                                              method = 'goodman', conf.level = .95, sides = 'two.sided'))
        
        q_si[,l,] = t(DescTools::MultinomCI(x = datalist$errormat_by_site[l,-datalist$mitsid_by_site[l]],
                                            method = 'goodman', conf.level = .95, sides = 'two.sided'))
        
      }
      
      MCMCout = list('mu_si' = mu_si, 'q_si' = q_si)
      
    }
    
    m_insample = MCMCout$mu_si
    # dimnames(m_insample)[[2]] = rownames(datalist$errormat_by_site)
    # dimnames(m_insample)[[3]] = colnames(datalist$errormat_by_site)
    MCMCout$m_insample = m_insample
    
    MCMCout$est.type = 'heterogen'
    
  }else if(model.choice=='pooled_empirical'){
    
    ## pooled empirical ----
    if(cause.type=='single'){
      
      ## single ----
      datalist$errormat = Reduce('+',
                                 lapply(1:datalist$nSite,
                                        FUN = function(s){
                                          
                                          datalist$errormat_by_site[((s-1)*datalist$nCause + 1):(s*datalist$nCause),]
                                          
                                        }))
      mu_i = array(dim = c(3, datalist$nCause,
                           datalist$nCause))
      q_i = array(dim = c(3, datalist$nCause,
                          datalist$nCause-1))
      for(i in 1:datalist$nCause){
        
        mu_i[,i,] =  t(DescTools::MultinomCI(x = datalist$errormat[i,],
                                             method = 'goodman', conf.level = .95, sides = 'two.sided'))
        
        q_i[,i,] = t(DescTools::MultinomCI(x = datalist$errormat[i,-i],
                                           method = 'goodman', conf.level = .95, sides = 'two.sided'))
        
      }
      
      mu_si = array(dim = c(3, datalist$nCause*datalist$nSite,
                            datalist$nCause),
                    dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                    colnames(datalist$errormat_by_site)))
      q_si = array(dim = c(3, datalist$nCause*datalist$nSite,
                           datalist$nCause-1),
                   dimnames = list(NULL, rownames(datalist$errormat_by_site),
                                   NULL))
      for(s in 1:datalist$nSite){
        
        mu_si[,((s-1)*datalist$nCause + 1):(s*datalist$nCause),] = mu_i
        q_si[,((s-1)*datalist$nCause + 1):(s*datalist$nCause),] = q_i
        
      }
      
      MCMCout = list('mu_i' = mu_i, 'q_i' = q_i,
                     'mu_si' = mu_si, 'q_si' = q_si)
      
    }
    
    MCMCout$m_insample = MCMCout$mu_si
    
    MCMCout$est.type = 'homogen'
    
  }
  
  ## preparing output list ----
  if(!(model.choice %in% c('separate_empirical', 'pooled_empirical'))){
    
    ### information criterion for model selection performacne ----
    #### loo ic ----
    MCMCout$loo.out = loo::loo(MCMCout$loglik, 
                               r_eff = loo::relative_eff(exp(MCMCout$loglik), 
                                                         chain_id = rep(1, nrow(MCMCout$loglik))), 
                               cores = 1)
    
    #### waic ----
    MCMCout$waic.out = loo::waic(MCMCout$loglik)
    
    ic.df = rbind(MCMCout$waic.out$estimates,
                  MCMCout$loo.out$estimates)
    
    ic.df.melt = as.numeric(t(ic.df))
    names(ic.df.melt) = paste0(rep(rownames(ic.df), each = ncol(ic.df)),'_',rep(colnames(ic.df), nrow(ic.df)))
    
    MCMCout$ic.df = ic.df.melt
    
    #### mcmc diagnostic ----
    MCMCout$mcmc.diagnostic = data.frame('max_Rhat' = max_Rhat, 'min_ess_bulk' = min_ess_bulk,
                                        'num_divergent' = rstan::get_num_divergent(stanfit),
                                        'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit))
    
    if(verbose){
      
      cat('\n')
      print(MCMCout$mcmc.diagnostic)
      cat('\n')
      
    }
    
    ## list of country-specific misclassification matrix posterior/predictive distributions ----
    Mmat = lapply(datalist$sites,
                  FUN = function(site_s){
                    
                    if(cause.type=='single'){
                      
                      Mmat_s = MCMCout$m_insample[,datalist$sites_by_site==site_s,]
                      NA.causeid = which(datalist$mits.total_by_site[datalist$sites_by_site==site_s]==0)
                      
                      if(length(NA.causeid)>0){Mmat_s[,NA.causeid,] = MCMCout$m_pred[,NA.causeid,]}
                      
                    }
                    
                    dimnames(Mmat_s)[[2]] = dimnames(Mmat_s)[[3]] = datalist$causes
                    
                    return(Mmat_s)
                    
                  })
    Mmat = c(Mmat, list(MCMCout$m_pred))
    names(Mmat) = c(datalist$sites, 'other')
    dimnames(Mmat$other)[[2]] = dimnames(Mmat$other)[[3]] = datalist$causes
    
    ## list of country-specific misclassification matrix point estimates (posterior mean) ----
    Mmat.postmean = lapply(1:length(Mmat),
                           FUN = function(s){
                             
                             apply(Mmat[[s]], 2:3, mean)
                             
                           })
    names(Mmat.postmean) = names(Mmat)
    
    ## list of dirichlet parameters approximating country-specific misclassification posterior/predictive distributions ----
    Mmat.asDirich = lapply(1:length(Mmat),
                           FUN = function(s){
                             
                             alpha.Dirich_s = t(apply(Mmat[[s]], 2,
                                                      FUN = function(v){
                                                        
                                                        mle.out = sirt::dirichlet.mle(x = v)
                                                        mle.out$alpha
                                                        
                                                      }))
                             
                             rownames(alpha.Dirich_s) = colnames(alpha.Dirich_s) = datalist$causes
                             
                             return(alpha.Dirich_s)
                             
                           })
    names(Mmat.asDirich) = names(Mmat)
    
    # output list
    fitout = c(list('MCMCout' = MCMCout),
               list('stanobj' = stanfit, 
                    'Mmat' = Mmat,
                    'Mmat.postmean' = Mmat.postmean, 'Mmat.asDirich' = Mmat.asDirich,
                    'datalist' = datalist))
    
  }else if(model.choice %in% c('separate_empirical', 'pooled_empirical')){
    
    # list of country-specific misclassification matrix
    Mmat = lapply(1:datalist$nSite,
                  FUN = function(s){
                    
                    Mmat_s = MCMCout$m_insample[,datalist$siteid_by_site==s,]
                    
                    dimnames(Mmat_s)[[2]] = dimnames(Mmat_s)[[3]] = datalist$causes
                    
                    return(Mmat_s)
                    
                  })
    names(Mmat) = datalist$sites
    
    # output list
    fitout = c(list('MCMCout' = MCMCout),
               list('Mmat' = Mmat, 'datalist' = datalist))
    
  }
  
  if(saveoutput){
    
    if(is.null(output_dir)) output_dir = 'hetmis_output'
    if(is.null(output_filename)) output_filename = paste0(model.choice, '_', cause.type)
    
    input.list.now = as.list(environment())
    input.list = input.list.now[names(input.list)]
    # fitout$input$output_dir = output_dir
    # fitout$input$output_filename = output_filename
    
    fitout = c(fitout, list('input' = input.list))
    
    if(!dir.exists(output_dir)){
      
      dir.create(output_dir)
      
    }
    
    saveRDS(fitout, file.path(output_dir, output_filename))
    
  }else{
    
    input.list.now = as.list(environment())
    input.list = input.list.now[names(input.list)]
    
    fitout = c(fitout, list('input' = input.list))
    
    return(fitout)
    
  }
  
}


# sequential va-calibration ----
calibratedva_v2 <- function(va_unlabeled = NULL, cause.type = "single",
                            model.choice = c("seqpshrink_mfixed", "seqpshrink_up"), 
                            Mmat.fixed = NULL,
                            Mmat.asDirich = NULL, pss = 4,
                            nMCMC = 10000, nBurn = 50000, nThin = 1, adapt_delta_stan = .95,
                            # pseudo_samplesize = 100,
                            sourcecode.path = file.path(getwd(), "sourcecode"),
                            seed = 1, verbose = T, saveoutput = T,
                            output_dir = NULL, output_filename = NULL){
  
  input.list = as.list(environment())
  
  if(verbose){
    
    refresh.stan = max((nBurn + nMCMC*nThin)/10, 1)
    
  }else{
    
    refresh.stan = 0
    
  }
  
  
  # data checks and preparing for modeling ----
  if(!is.list(va_unlabeled)){
    
    K = 1
    if(!any(is.vector(va_unlabeled), is.matrix(va_unlabeled))){
      
      stop("'va_unlabeled' must be either a vector (death counts for each cause), or a matrix (individuals along rows, causes along columns).")
      
    }else if(is.matrix(va_unlabeled)){
      
      if(is.null(colnames(va_unlabeled))) stop("Columns of 'va_unlabeled' must have causes as names.")  
      
      causes = colnames(va_unlabeled)
      nCause = length(causes)
      
      if(is.null(rownames(va_unlabeled))) print("Rows of 'va_unlabeled' doesn't have names. Assuming they correspond to individuals.")
      
      singlecause.check = (rowSums(va_unlabeled)==1)&(apply(X = va_unlabeled, 1, FUN = function(v){sum(v!=0)})==1)
      if(!all(singlecause.check)){
        
        cat("\n")
        print(va_unlabeled[singlecause.check,])
        cat("\n")
        stop("'va_unlabeled' doesn't look like single cause predictions for the above cases.")
        
      }
      
      va_deaths = array(dim = c(K, nCause), 
                        dimnames = list('algo1', causes))
      va_deaths[1,] = colSums(va_unlabeled)
      
    }else if(is.vector(va_unlabeled)){
      
      if(is.null(names(va_unlabeled))) stop("Components of 'va_unlabeled' must have causes as names.")  
      
      causes = names(va_unlabeled)
      nCause = length(causes)
      
      count.check = sum(va_unlabeled)==sum(floor(va_unlabeled))
      if(!count.check) stop("'va_unlabeled' must be a vector of death counts for each cause.")
      
      va_deaths = array(dim = c(K, nCause), 
                        dimnames = list('algo1', causes))
      va_deaths[1,] = va_unlabeled
      
    }
    
  }else{
    
    if(is.null(names(va_unlabeled))) stop("Components of 'va_unlabeled' must have algorithms as names.")
    
    K = length(va_unlabeled)
    for(k in 1:K){
      
      if(is.matrix(va_unlabeled[[k]])){
        
        if(is.null(colnames(va_unlabeled))) stop(paste0("Columns of 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," must have causes as names."))
        
        if(is.null(rownames(va_unlabeled[[k]]))) print(paste0("Rows of 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," doesn't have names. Assuming they correspond to individuals."))
        
        singlecause.check = (rowSums(va_unlabeled[[k]])==1)&(apply(X = va_unlabeled[[k]], 1, FUN = function(v){sum(v!=0)})==1)
        if(!all(singlecause.check)){
          
          cat("\n")
          print(va_unlabeled[[k]][singlecause.check,])
          cat("\n")
          stop(paste0("'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," doesn't look like single cause predictions for the above cases."))
          
        }
        
        causes_temp = colnames(va_unlabeled[[k]])
        
        if(k==1){
          
          causes = causes_temp
          nCause = length(causes)
          
          va_deaths = array(dim = c(K, nCause), 
                            dimnames = list(names(va_unlabeled), causes))
          
        }else if(k>1){
          
          if(!identical(causes_temp, causes)) stop("Causes input through 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," does not match with algorithm ", names(va_unlabeled)[k-1],".")
          
        }
        
        va_deaths[k,] = colSums(va_unlabeled[[k]])
        
        
      }else if(is.vector(va_unlabeled[[k]])){
        
        if(is.null(names(va_unlabeled[[k]]))) stop("Components of 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," must have causes as names.")
        
        causes_temp = names(va_unlabeled[[k]])
        
        count.check = sum(va_unlabeled[[k]])==sum(floor(va_unlabeled[[k]]))
        if(!count.check) stop(paste0("'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," doesn't look like a count vector. It must be a vector of death counts for each cause."))
        
        if(k==1){
          
          causes = causes_temp
          nCause = length(causes)
          
          va_deaths = array(dim = c(K, nCause), 
                            dimnames = list(names(va_unlabeled), causes))
          
        }else if(k>1){
          
          if(!identical(causes_temp, causes)) stop("Causes input through 'va_unlabeled' for algorithm ", names(va_unlabeled)[k]," does not match with algorithm ", names(va_unlabeled)[k-1],".")
          
        }
        
        va_deaths[k,] = va_unlabeled[[k]]
        
      }
      
    }
    
  }
  
  
  # observed uncalibrated csmf to shrink to
  p.uncalib.obs = colSums(va_deaths)/sum(va_deaths)
  if(any(p.uncalib.obs==0)){
    
    p.uncalib.obs[which.max(p.uncalib.obs)] = p.uncalib.obs[which.max(p.uncalib.obs)] - 0.001
    p.uncalib.obs[p.uncalib.obs==0] = 0.001/sum(p.uncalib.obs==0)
    p.uncalib.obs[which.max(p.uncalib.obs)] = 1 - sum(p.uncalib.obs[-which.max(p.uncalib.obs)])
    
  }
  
  
  # calibration modeling ----
  if(!(model.choice %in% c("seqpshrink_mfixed", "seqpshrink_up"))) {
    
    stop("'model.choice' must be either mshrink, pshrink, seqpshrink_mfixed, or seqpshrink_up")
    
  }else if(model.choice=="seqpshrink_mfixed") {
    
    ## fixed Mmat ----
    
    # checks
    if(is.null(Mmat.fixed)) stop("Must specify 'Mmat.fixed' for model.choice 'seqpshrink_mfixed'.")
    if(K==1){
      
      if(!is.matrix(Mmat.fixed)) stop("'Mmat.fixed' must be a matrix with MITS causes along the rows and VA causes along the columns.")
      if(is.null(colnames(Mmat.fixed))) stop("'Mmat.fixed' must have VA causes as column names.")
      if(!identical(colnames(Mmat.fixed), causes)) stop("VA causes specified in 'Mmat.fixed' do not match with that specified in 'va_unlabeled'.")
      if(is.null(rownames(Mmat.fixed))) stop("'Mmat.fixed' must have MITS causes as row names.")
      if(!identical(rownames(Mmat.fixed), causes)) stop("MITS causes specified in 'Mmat.fixed' do not match with that specified in 'va_unlabeled'.")
      if(!all.equal(target = unname(rowSums(Mmat.fixed)), current = rep(1, nrow(Mmat.fixed)))) stop("Each row of 'Mmat.fixed' must sum to 1.")
      # if(any(unname(rowSums(Mmat.fixed))!=1)) stop("Each row of 'Mmat.fixed' must sum to 1.")
      
      Mmat.touse = array(dim = c(K, nCause, nCause),
                         dimnames = list("algo1", causes, causes))
      Mmat.touse[1,,] = Mmat.fixed
      
    }else if(K>1){
      
      if(!is.list(Mmat.fixed)) stop("'Mmat.fixed' must be a list of matrices for each algorithm. Rows and columns of each matrix must correspond to MITS and VA causes.")
      
      if(is.null(names(Mmat.fixed))) stop("Components of 'Mmat.fixed' must have algorithms as names.")
      if(!identical(names(Mmat.fixed), names(va_unlabeled))) stop("Algorithm names specified in 'Mmat.fixed' and 'va_unlabeled' do not match.")
      
      Mmat.touse = array(dim = c(K, nCause, nCause),
                         dimnames = list(names(Mmat.fixed), causes, causes))
      
      for(k in 1:K){
        
        if(!is.matrix(Mmat.fixed[[k]])) stop(paste0("'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k]," must be a matrix with MITS causes along the rows and VA causes along the columns."))
        if(is.null(colnames(Mmat.fixed[[k]]))) stop(paste0("'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must have VA causes as column names."))
        if(!identical(colnames(Mmat.fixed[[k]]), causes)) stop(paste0("VA causes specified in 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " do not match with that specified in 'va_unlabeled'."))
        if(is.null(rownames(Mmat.fixed[[k]]))) stop(paste0("'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must have MITS causes as row names."))
        if(!identical(rownames(Mmat.fixed[[k]]), causes)) stop(paste0("MITS causes specified in 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " do not match with that specified in 'va_unlabeled'."))
        if(!all.equal(target = unname(rowSums(Mmat.fixed[[k]])), current = rep(1, nrow(Mmat.fixed[[k]])))) stop(paste0("Each row of 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must sum to 1."))
        # if(any(rowSums(Mmat.fixed)!=1)) stop(paste0("Each row of 'Mmat.fixed' for algorithm ", names(Mmat.fixed)[k], " must sum to 1."))
        
        Mmat.touse[k,,] = Mmat.fixed[[k]]
        
      }
      
    }
    
    # stan fit ====
    # Mmat.init = Mmat.prior
    # for(k in 1:K) Mmat.init[k,,] = Mmat.prior[k,,]/rowSums(Mmat.prior[k,,])
    
    # print(Mmat.init)
    # print(p.uncalib.obs)
    
    stanfit = rstan::sampling(seqcalib_mmat.stan_object,
                              pars = c('p_calib', 'loglik'),
                              include = T,
                              data = list('nCause' = nCause,
                                          'nAlgo' = K,
                                          'aj' = va_deaths,
                                          'p_uncalib_obs' = p.uncalib.obs,
                                          'Mmat' = Mmat.touse,
                                          'pss' = pss
                              ),
                              chains = 1, 
                              iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                              control = list('adapt_delta' = adapt_delta_stan),
                              seed = seed,
                              refresh = refresh.stan#,
                              # init = list(list('Mmat' = Mmat.init,
                              #                  'p_calib' = p.uncalib.obs))
    )
    
    MCMCout = rstan::extract(stanfit) # MCMC output
    
    # max Rhat
    max_Rhat = max(apply(X = MCMCout$p_calib, 2,
                         FUN = function(v){
                           
                           rstan::Rhat(v)
                           
                         }))
    
    #### min bulk ESS ----
    min_ess_bulk = min(apply(X = MCMCout$p_calib, 2,
                             FUN = function(v){
                               
                               rstan::ess_bulk(v)
                               
                             }))/nMCMC
    
  }else if(model.choice=="seqpshrink_up") {
    
    ## Mmat prior as independent Dirichlet ----
    
    # checks
    if(is.null(Mmat.asDirich)) stop("Must specify 'Mmat.asDirich' for model.choice 'seqpshrink_up'.")
    if(K==1){
      
      if(!is.matrix(Mmat.asDirich)) stop("'Mmat.asDirich' must be a matrix with MITS causes along the rows and VA causes along the columns.")
      if(is.null(colnames(Mmat.asDirich))) stop("'Mmat.asDirich' must have VA causes as column names.")
      if(!identical(colnames(Mmat.asDirich), causes)) stop("VA causes specified in 'Mmat.asDirich' do not match with that specified in 'va_unlabeled'.")
      if(is.null(rownames(Mmat.asDirich))) stop("'Mmat.asDirich' must have MITS causes as row names.")
      if(!identical(rownames(Mmat.asDirich), causes)) stop("MITS causes specified in 'Mmat.asDirich' do not match with that specified in 'va_unlabeled'.")
      if(sum(Mmat.asDirich==0)>0) stop("'Mmat.asDirich' cannot have zeros.")
      
      Mmat.prior = array(dim = c(K, nCause, nCause),
                         dimnames = list("algo1", causes, causes))
      Mmat.prior[1,,] = Mmat.asDirich
      
    }else if(K>1){
      
      if(!is.list(Mmat.asDirich)) stop("'Mmat.asDirich' must be a list of matrices for each algorithm. Rows and columns of each matrix must correspond to MITS and VA causes.")
      
      if(is.null(names(Mmat.asDirich))) stop("Components of 'Mmat.asDirich' must have algorithms as names.")
      if(!identical(names(Mmat.asDirich), names(va_unlabeled))) stop("Algorithm names specified in 'Mmat.asDirich' and 'va_unlabeled' do not match.")
      
      Mmat.prior = array(dim = c(K, nCause, nCause),
                         dimnames = list(names(Mmat.asDirich), causes, causes))
      
      for(k in 1:K){
        
        if(!is.matrix(Mmat.asDirich[[k]])) stop(paste0("'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k]," must be a matrix with MITS causes along the rows and VA causes along the columns."))
        if(is.null(colnames(Mmat.asDirich[[k]]))) stop(paste0("'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " must have VA causes as column names."))
        if(!identical(colnames(Mmat.asDirich[[k]]), causes)) stop(paste0("VA causes specified in 'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " do not match with that specified in 'va_unlabeled'."))
        if(is.null(rownames(Mmat.asDirich[[k]]))) stop(paste0("'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " must have MITS causes as row names."))
        if(!identical(rownames(Mmat.asDirich[[k]]), causes)) stop(paste0("MITS causes specified in 'Mmat.asDirich' for algorithm ", names(Mmat.asDirich)[k], " do not match with that specified in 'va_unlabeled'."))
        
        Mmat.prior[k,,] = Mmat.asDirich[[k]]
        
      }
      
    }
    
    # stan fit ====
    # Mmat.init = Mmat.prior
    # for(k in 1:K) Mmat.init[k,,] = Mmat.prior[k,,]/rowSums(Mmat.prior[k,,])
    
    # print(Mmat.init)
    # print(p.uncalib.obs)
    
    stanfit = rstan::sampling(seqcalib.stan_object,
                              pars = c('Mmat', 'p_calib', 'loglik'),
                              include = T,
                              data = list('nCause' = nCause,
                                          'nAlgo' = K,
                                          'aj' = va_deaths,
                                          'p_uncalib_obs' = p.uncalib.obs,
                                          # 'p_uncalib_obs' = rep(1/nCause, nCause),
                                          'Mmatprior_asDirich' = Mmat.prior,
                                          'pss' = pss
                              ),
                              chains = 1, 
                              iter = nBurn + nMCMC*nThin, warmup = nBurn, thin = nThin,
                              control = list('adapt_delta' = adapt_delta_stan),
                              seed = seed,
                              refresh = refresh.stan#,
                              # init = list(list('Mmat' = Mmat.init,
                              #                  'p_calib' = p.uncalib.obs))
    )
    
    MCMCout = rstan::extract(stanfit) # MCMC output
    
    # max Rhat
    max_Rhat = max(apply(X = MCMCout$Mmat, 2:4,
                         FUN = function(v){
                           
                           rstan::Rhat(v)
                           
                         }),
                   apply(X = MCMCout$p_calib, 2,
                         FUN = function(v){
                           
                           rstan::Rhat(v)
                           
                         }))
    
    #### min bulk ESS ----
    min_ess_bulk = min(apply(X = MCMCout$Mmat, 2:4,
                             FUN = function(v){
                               
                               rstan::ess_bulk(v)
                               
                             }),
                       apply(X = MCMCout$p_calib, 2,
                             FUN = function(v){
                               
                               rstan::ess_bulk(v)
                               
                             }))/nMCMC
    
    dimnames(MCMCout$Mmat) = list(NULL, names(Mmat.asDirich), causes, causes)
    
    MCMCout$Mmat.postmean = apply(MCMCout$Mmat, 2:4, mean)
    
  }
  
  
  ### model diagnostics ----
  # loo ic
  MCMCout$loo.out = loo::loo(MCMCout$loglik, 
                             r_eff = loo::relative_eff(exp(MCMCout$loglik), 
                                                       chain_id = rep(1, nrow(MCMCout$loglik))), 
                             cores = 1)
  
  # waic
  MCMCout$waic.out = loo::waic(MCMCout$loglik)
  
  ic.df = rbind(MCMCout$waic.out$estimates,
                MCMCout$loo.out$estimates)
  
  ic.df.melt = as.numeric(t(ic.df))
  names(ic.df.melt) = paste0(rep(rownames(ic.df), each = ncol(ic.df)),'_',rep(colnames(ic.df), nrow(ic.df)))
  
  MCMCout$ic.df = ic.df.melt
  
  MCMCout$mcmc.diagnostic = data.frame('max_Rhat' = max_Rhat, 'min_ess_bulk' = min_ess_bulk,
                                       'num_divergent' = rstan::get_num_divergent(stanfit),
                                       'num_max_treedepth' = rstan::get_num_max_treedepth(stanfit))
  
  cat("\n")
  print(MCMCout$mcmc.diagnostic)
  cat("\n")
  
  colnames(MCMCout$p_calib) = causes
  MCMCout$p_calib.postmean = colMeans(MCMCout$p_calib)
  
  
  ## output list ----
  output = c(list('MCMCout' = MCMCout, "p_uncalib_obs" = p.uncalib.obs),
             list('input' = input.list))
  if(saveoutput){
    
    if(is.null(output_dir)) output_dir = 'calibratedva_output'
    if(is.null(output_filename)) output_filename = model.choice
    
    if(!dir.exists(output_dir)){
      
      dir.create(output_dir)
      
    }
    
    saveRDS(output, file.path(output_dir, output_filename))
    
  }else{
    
    return(output)
    
  }
  
}


# function for creating figures given a data and model outputs ----
makefigs = function(output.list, output.empirical,
                    empirical.type1, empirical.type2, site.level = NULL, ic.llim = NULL,
                    plot.title, plot.subtitle, xlab.title, ylab.title,
                    cause.label = NULL, site.label = NULL, method.label = NULL, ic.label = NULL,
                    yrange.pss = NULL,
                    pss.label = NULL, add.mean = T, pss.inv, mybreaks,
                    method.colors = NULL,
                    pointshape = NULL, pointcolors = NULL, pointstroke = 2, pt.size = 7, opacity = .3,
                    linetype = NULL, linewidth = 1.5,
                    barcolors,
                    xtext.angle = 30, 
                    xhjust = .5, xvjust = .5,
                    limits.colorgradient = NULL, n.breaks.colorgradient = NULL,
                    plot.title.size = 60, plot.subtitle.size = 50, axis.title.size = 50,
                    axis.text.size = 30, strip.text.size = 40, 
                    legend.title.size = 40, legend.key = 30, legend.text.size = 40,
                    legend.key.width = 2, legend.key.height = 2,
                    legend.key.size = 2, legend.spacing.x = 2, 
                    legend.spacing.y = 2, 
                    remove.01 = F,
                    plot.margin, nrow.legend = 2, legend.position = 'bottom',
                    toplot = 'mu'){
  
  plotlist = NULL
  
  if(toplot=='mu_post'){
    
    ## posterior misclassification rates ----
    
    plotdf = cause.level = #site.level = 
      model.est.type = NULL
    
    # temporary model names for each output
    model.names = paste0('model', 1:length(output.list))
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    for(i in 1:length(output.list)){
      
      if(output.list[[i]]$input$model.choice %in% c('mits',
                                                    'mits+site_sens', 'mits+site',
                                                    'separate_wdecomp', 'separate_wodecomp')){
        
        postsumm = apply(output.list[[i]]$MCMCout$m_insample, 2:3,
                         FUN = function(v){
                           
                           c(mean(v),
                             quantile(x = v, probs = c(.025, .975), na.rm = T))
                           
                         })
        
        plotdf = rbind.data.frame(plotdf,
                                  data.frame('est' = as.numeric(postsumm[1,,]),
                                             'llim' = as.numeric(postsumm[2,,]),
                                             'ulim' = as.numeric(postsumm[3,,]),
                                             'mits' = rep(output.list[[i]]$datalist$mits_by_site,
                                                          output.list[[i]]$datalist$nCause),
                                             'va' = rep(output.list[[i]]$datalist$causes, 
                                                        each = output.list[[i]]$datalist$nCause*
                                                          output.list[[i]]$datalist$nSite),
                                             'site' = rep(output.list[[i]]$datalist$sites_by_site,
                                                          output.list[[i]]$datalist$nCause),
                                             'method' = model.names[i],
                                             'est.type' = output.list[[i]]$MCMCout$est.type))
        
        model.est.type = c(model.est.type, output.list[[i]]$MCMCout$est.type)
        
      }else if(output.list[[i]]$input$model.choice %in% c('separate_empirical', 'pooled_empirical')){
        
        plotdf = rbind.data.frame(plotdf,
                                  data.frame('est' = as.numeric(output.list[[i]]$MCMCout$m_insample[1,,]),
                                             'llim' = as.numeric(output.list[[i]]$MCMCout$m_insample[2,,]),
                                             'ulim' = as.numeric(output.list[[i]]$MCMCout$m_insample[3,,]),
                                             'mits' = rep(output.list[[i]]$datalist$mits_by_site,
                                                          output.list[[i]]$datalist$nCause),
                                             'va' = rep(output.list[[i]]$datalist$causes, 
                                                        each = output.list[[i]]$datalist$nCause*
                                                          output.list[[i]]$datalist$nSite),
                                             'site' = rep(output.list[[i]]$datalist$sites_by_site,
                                                          output.list[[i]]$datalist$nCause),
                                             'method' = model.names[i],
                                             'est.type' = output.list[[i]]$MCMCout$est.type))
        
        model.est.type = c(model.est.type, output.list[[i]]$MCMCout$est.type)
        
      }
      
    }
    
    # print(model.names)
    # print(model.est.type)
    
    plotdf$linestyle = plotdf$shapestyle = 
      paste0(plotdf$method, '_', plotdf$est.type)
    plotdf$linestyle[plotdf$est.type=='heterogen'] = NA
    plotdf$shapestyle[plotdf$est.type=='homogen'] = NA
    
    if(is.null(cause.level)) cause.level = output.list[[1]]$datalist$causes
    if(is.null(site.level)) site.level = output.list[[1]]$datalist$sites
    
    # making sure of the country and method order for plotting
    plotdf$mits = factor(x = plotdf$mits, levels = cause.level)
    plotdf$va = factor(x = plotdf$va, levels = cause.level)
    plotdf$site = factor(x = plotdf$site, levels = site.level)
    plotdf$method = factor(x = plotdf$method, levels = model.names)
    plotdf$linestyle = factor(x = plotdf$linestyle)
    plotdf$shapestyle = factor(x = plotdf$shapestyle)
    
    print(levels(plotdf$linestyle))
    print(levels(plotdf$shapestyle))
    
    if(is.null(cause.label)) cause.label = cause.level
    
    if(output.list[[1]]$input$cause.type=='single'){
      
      mits.freq = unlist(lapply(1:output.list[[1]]$datalist$nCause,
                                FUN = function(i){
                                  
                                  sum(output.list[[1]]$datalist$mits.total_by_site[output.list[[1]]$datalist$mitsid_by_site==i])
                                  
                                }))
      va.freq = unname(colSums(output.list[[1]]$datalist$errormat_by_site))
      
      cause.label.row = paste0(cause.label, ' (', mits.freq, ')')
      cause.label.col = paste0(cause.label, ' (', va.freq, ')')
      
    }
    names(cause.label.row) = names(cause.label.col) = cause.level
    
    if(is.null(site.label)) site.label = site.level
    # names(site.label) = site.level
    
    idmodel.homogen = model.est.type=='homogen'
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = plotdf,
                      ggplot2::aes(x = site, y = est)) +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(mits~va,
                          labeller = ggplot2::labeller(mits = cause.label.row,
                                                       va = cause.label.col)) +
      ggplot2::scale_x_discrete(labels = site.label)
    
    if(length(model.names)==1){
      
      if(model.est.type=='homogen'){
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_rect(data = subset(plotdf, est.type=='homogen'),
                             aes(xmin = -Inf, xmax = Inf, ymin = llim, ymax = ulim),
                             fill = 'grey85', alpha = opacity) +
          ggplot2::geom_hline(data = subset(plotdf, est.type=='homogen'),
                              ggplot2::aes(yintercept = est),
                              linetype = 'solid', color = 'black', linewidth = linewidth)
        # ggplot2::geom_ribbon(data = subset(plotdf, est.type=='homogen'),
        #                      aes(x = as.numeric(site),
        #                          # xmin = as.numeric(site)-1, xmax = as.numeric(site)+1, 
        #                          ymin = llim, ymax = ulim),
        #                      fill = 'grey50', alpha = opacity)
        
      }
      
      if(model.est.type=='heterogen'){
        
        if(is.null(pointcolors)) pointcolors = scales::hue_pal()(length(site.level))
        # names(pointcolors) = site.level
        
        if(is.null(pointshape)) pointshape = 8
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_point(data = subset(plotdf, est.type=='heterogen'),
                              ggplot2::aes(x = site, y = est, group = method, color = site),
                              shape = pointshape, size = pt.size, stroke = pointstroke,
                              position = ggplot2::position_dodge(width = 1)) +
          ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim, color = site),
                                 width=.5, linewidth = linewidth, alpha = opacity,
                                 position = ggplot2::position_dodge(width = 1)) +
          ggplot2::scale_color_manual(values = pointcolors,
                                      labels = site.label)
        
      }
      
    }else{
      
      if(sum(idmodel.homogen)>0){
        
        if(is.null(linetype)) linetype = 1:length(levels(plotdf$linestyle))
        names(linetype) = levels(plotdf$linestyle)
        
        linetype.label = method.label[idmodel.homogen]
        names(linetype.label) = levels(plotdf$linestyle)
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_rect(data = subset(plotdf, est.type=='homogen'),
                             aes(xmin = -Inf, xmax = Inf, ymin = llim, ymax = ulim),
                             fill = 'grey85', alpha = opacity) +
          ggplot2::geom_hline(data = subset(plotdf, est.type=='homogen'),
                              ggplot2::aes(yintercept = est, linetype = linestyle),
                              color = 'black', linewidth = linewidth) +
          # ggplot2::geom_ribbon(data = subset(plotdf, est.type=='homogen'),
          #             aes(x = as.numeric(site), ymin = llim, ymax = ulim),
          #             fill = 'grey50', alpha = opacity) +
          ggplot2::scale_linetype_manual(values = linetype,
                                         labels = linetype.label)
        
      }
      
      if(sum(!idmodel.homogen)>0){
        
        if(is.null(pointcolors)) pointcolors = scales::hue_pal()(length(site.level))
        # names(pointcolors) = site.level
        
        if(is.null(pointshape)) pointshape = 8 + 0:(length(levels(plotdf$shapestyle))-1)
        names(pointshape) = levels(plotdf$shapestyle)
        
        shape.label = method.label[!idmodel.homogen]
        names(shape.label) = levels(plotdf$shapestyle)
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_point(data = subset(plotdf, est.type=='heterogen'),
                              ggplot2::aes(x = site, y = est, 
                                           group = shapestyle, color = site, shape = shapestyle),
                              size = pt.size, stroke = pointstroke,
                              position = ggplot2::position_dodge(width = 1)) +
          ggplot2::geom_errorbar(data = subset(plotdf, est.type=='heterogen'),
                                 ggplot2::aes(x = site, ymin = llim, ymax = ulim, 
                                              group = shapestyle, color = site),
                                 width=.5, linewidth = linewidth, alpha = opacity,
                                 position = ggplot2::position_dodge(width = 1)) +
          ggplot2::scale_color_manual(values = pointcolors,
                                      labels = site.label) +
          ggplot2::scale_shape_manual(values = pointshape,
                                      labels = shape.label)
        
      }
      
    }
    
    plotdf.plot = plotdf.plot +
      ggplot2::geom_rect(data = subset(plotdf, mits==va), 
                         fill = NA, colour = "red3", linewidth = 8,
                         xmin = -Inf, xmax = Inf,
                         ymin = -Inf, ymax = Inf) +
      ggplot2::geom_rect(data = subset(plotdf, mits!=va), 
                         fill = NA, colour = "black", linewidth = 3,
                         xmin = -Inf, xmax = Inf,
                         ymin = -Inf, ymax = Inf) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=plot.title.size, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=plot.subtitle.size),
                     axis.title.x = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size, 
                                                         angle = xtext.angle, hjust = xhjust, vjust = xvjust),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 1.5),
                     axis.ticks.length.x = unit(.5, "cm"),
                     axis.ticks.y = element_line(linewidth = 1.5),
                     axis.ticks.length.y = unit(.5, "cm"),
                     # axis.ticks.x = ggplot2::element_blank(),
                     # panel.border = element_rect(color='black', linetype = "solid", 
                     #                             fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(legend.key.width, "cm"), 
                     legend.key.height = ggplot2::unit(legend.key.height, "cm"), 
                     legend.key.size = ggplot2::unit(legend.key, "cm"),
                     legend.spacing.x = ggplot2::unit(legend.spacing.x, 'cm'), 
                     legend.spacing.y = ggplot2::unit(legend.spacing.y, 'cm'), 
                     legend.text=ggplot2::element_text(size=legend.text.size),
                     legend.position = legend.position)
    
    if(!missing(nrow.legend)){
      
      plotdf.plot = plotdf.plot +
        ggplot2::guides(color = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                                      # override.aes = list(
                                                      #   linetype = "solid",
                                                      #   shape = c(16, rep(NA, length(model.names)-1)))
        ),
        shape = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                      # override.aes = list(
                                      #   linetype = "solid",
                                      #   shape = c(16, rep(NA, length(model.names)-1)))
        ),
        linetype = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                         # override.aes = list(
                                         #   linetype = "solid",
                                         #   shape = c(16, rep(NA, length(model.names)-1)))
        )) #+
      # ggplot2::guides(linetype = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
      #                                               override.aes = list(
      #                                                 linetype = linetype_toplot
      #                                                 )),
      #                 shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
      #                                               override.aes = list(
      #                                                 shape = shape_toplot
      #                                               )),
      #                 color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
      #                                               override.aes = list(
      #                                                 color = color_toplot
      #                                               ))
      #                 )
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='sens'){
    
    ## sensitivities ----
    
    if(!missing(output.list)){
      
      # temporary model names for each output
      model.names = paste('model', 1:length(output.list), sep = '')
      
      plotdf = do.call('rbind.data.frame',
                       lapply(X = 1:length(output.list),
                              FUN = function(X){
                                
                                postsumm = apply(output.list[[X]]$MCMCout$musij_by_mits, 2:3,
                                                 FUN = function(v){
                                                   
                                                   c(mean(v),
                                                     quantile(x = v, probs = c(.025, .975), na.rm = T))
                                                   
                                                 })
                                
                                data.frame('est' = as.numeric(postsumm[1,,]),
                                           'llim' = as.numeric(postsumm[2,,]),
                                           'ulim' = as.numeric(postsumm[3,,]),
                                           'mits' = rep(rep(output.list[[X]]$datalist$causes, 
                                                            each = length(output.list[[X]]$datalist$sites)),
                                                        length(output.list[[X]]$datalist$causes)),
                                           'va' = rep(output.list[[X]]$datalist$causes, 
                                                      each = length(output.list[[X]]$datalist$causes)*
                                                        length(output.list[[X]]$datalist$sites)),
                                           'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                            length(output.list[[X]]$datalist$causes)),
                                                        length(output.list[[X]]$datalist$causes)),
                                           'method' = model.names[X])
                                
                              }))
      
      names(cause.label) = output.list[[1]]$datalist$causes
      
    }else{
      
      # data frame for ggplot
      plotdf = NULL
      
      # temporary model names for each output
      model.names = NULL
      
    }
    
    # adding empirical data
    if(!missing(output.empirical)){
      
      # data frame for ggplot
      plotdf = rbind.data.frame(data.frame('est' = as.numeric(output.empirical$musij_by_mits[1,,]),
                                           'llim' = as.numeric(output.empirical$musij_by_mits[2,,]),
                                           'ulim' = as.numeric(output.empirical$musij_by_mits[3,,]),
                                           'mits' = rep(rep(output.empirical$datalist$causes, 
                                                            each = length(output.empirical$datalist$sites)),
                                                        length(output.empirical$datalist$causes)),
                                           'va' = rep(output.empirical$datalist$causes, 
                                                      each = length(output.empirical$datalist$causes)*
                                                        length(output.empirical$datalist$sites)),
                                           'site' = rep(rep(output.empirical$datalist$sites,
                                                            length(output.empirical$datalist$causes)),
                                                        length(output.empirical$datalist$causes)),
                                           'method' = 'empirical'),
                                plotdf)
      
      if(is.null(names(cause.label))){
        
        names(cause.label) = output.empirical$datalist$causes
        
      }
      
      # temporary model names for each output
      model.names = c('empirical', model.names)
      
    }
    
    # making sure of the country and method order for plotting
    site.level = sort(unique(plotdf$site))
    plotdf$site = factor(x = plotdf$site, levels = site.level)  # alphabetical country order
    plotdf$method = factor(x = plotdf$method, levels = model.names)  # method order empirical (if given) then other outputs in their input order
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = subset(plotdf, mits==va),
                      ggplot2::aes(x = site, y = est, 
                                   group = method, color = site)) +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(.~va,
                          labeller = ggplot2::labeller(va = cause.label)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=plot.title.size, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=plot.subtitle.size),
                     axis.title.x = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size, 
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     # axis.ticks.x = ggplot2::element_blank(),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"), 
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position)
    
    if(length(model.names)==1){
      
      plotdf.plot = plotdf.plot +
        ggplot2::geom_point(shape = pointshape, size = pt.size, stroke = pointstroke,
                            position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim),
                               width=.5, linewidth = linewidth, alpha = opacity,
                               position = ggplot2::position_dodge(width = 1))
      
    }else{
      
      plotdf.plot = plotdf.plot +
        ggplot2::geom_point(ggplot2::aes(shape = method), size = pt.size, 
                            stroke = pointstroke,
                            position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim),
                               width=.5, linewidth = linewidth, alpha = opacity,
                               position = ggplot2::position_dodge(width = 1))
      
    }
    
    if(!missing(nrow.legend)){
      
      plotdf.plot = plotdf.plot +
        ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
                                                      override.aes = list(size = 10)),
                        shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
                                                      override.aes = list(size = 10)))
      
    }
    
    if(length(model.names)>1){
      
      names(pointshape) = model.names
      plotdf.plot = plotdf.plot + 
        ggplot2::scale_shape_manual(values = pointshape,
                                    labels = method.label)
      
    }
    
    if(!missing(pointcolors)){
      
      plotdf.plot = plotdf.plot + 
        ggplot2::scale_color_manual(values = pointcolors,
                                    labels = site.label) +
        ggplot2::scale_x_discrete(labels = site.label)
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='rfp'){
    
    ## rfp ----
    
    if(!missing(output.list)){
      
      # temporary model names for each output
      model.names = paste('model', 1:length(output.list), sep = '')
      
      plotdf = do.call('rbind.data.frame',
                       lapply(X = 1:length(output.list),
                              FUN = function(X){
                                
                                postsumm = apply(output.list[[X]]$MCMCout$musij_by_mits, 2:3,
                                                 FUN = function(v){
                                                   
                                                   c(mean(v),
                                                     quantile(x = v, probs = c(.025, .975), na.rm = T))
                                                   
                                                 })
                                
                                data.frame('est' = as.numeric(postsumm[1,,]),
                                           'llim' = as.numeric(postsumm[2,,]),
                                           'ulim' = as.numeric(postsumm[3,,]),
                                           'mits' = rep(rep(output.list[[X]]$datalist$causes, 
                                                            each = length(output.list[[X]]$datalist$sites)),
                                                        length(output.list[[X]]$datalist$causes)),
                                           'va' = rep(output.list[[X]]$datalist$causes, 
                                                      each = length(output.list[[X]]$datalist$causes)*
                                                        length(output.list[[X]]$datalist$sites)),
                                           'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                            length(output.list[[X]]$datalist$causes)),
                                                        length(output.list[[X]]$datalist$causes)),
                                           'method' = model.names[X])
                                
                              }))
      
      names(cause.label) = output.list[[1]]$datalist$causes
      
    }else{
      
      # data frame for ggplot
      plotdf = NULL
      
      # temporary model names for each output
      model.names = NULL
      
    }
    
    # adding empirical data
    if(!missing(output.empirical)){
      
      # data frame for ggplot
      plotdf = rbind.data.frame(do.call('rbind.data.frame',
                                        lapply(X = 1:length(output.empirical$datalist$causes),
                                               FUN = function(X){
                                                 
                                                 data.frame('est' = as.numeric(output.empirical$qsij_by_mits[1,((X-1)*length(output.empirical$datalist$sites) + 1):(X*length(output.empirical$datalist$sites)),]),
                                                            'llim' = as.numeric(output.empirical$qsij_by_mits[2,((X-1)*length(output.empirical$datalist$sites) + 1):(X*length(output.empirical$datalist$sites)),]),
                                                            'ulim' = as.numeric(output.empirical$qsij_by_mits[3,((X-1)*length(output.empirical$datalist$sites) + 1):(X*length(output.empirical$datalist$sites)),]),
                                                            'mits' = output.empirical$datalist$causes[X],
                                                            'va' = rep(output.empirical$datalist$causes[-X], 
                                                                       each = length(output.empirical$datalist$sites)),
                                                            'site' = rep(output.empirical$datalist$sites,
                                                                         length(output.empirical$datalist$causes)-1),
                                                            'method' = 'empirical')
                                                 
                                                 
                                               })),
                                plotdf)
      
      if(is.null(names(cause.label))){
        
        names(cause.label) = output.empirical$datalist$causes
        
      }
      
      # temporary model names for each output
      model.names = c('empirical', model.names)
      
    }
    
    # making sure of the country and method order for plotting
    site.level = sort(unique(plotdf$site))
    plotdf$site = factor(x = plotdf$site, levels = site.level)  # alphabetical country order
    plotdf$method = factor(x = plotdf$method, levels = model.names)  # method order empirical (if given) then other outputs in their input order
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = plotdf,
                      ggplot2::aes(x = site, y = est, 
                                   group = method, color = site)) +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(mits~va,
                          labeller = ggplot2::labeller(mits = cause.label,
                                                       va = cause.label)) +
      ggplot2::geom_rect(data = subset(plotdf, mits!=va), 
                         fill = NA, colour = "black", linewidth = 3,
                         xmin = -Inf, xmax = Inf,
                         ymin = -Inf, ymax = Inf) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=plot.title.size, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=plot.subtitle.size),
                     axis.title.x = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size, 
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     # axis.ticks.x = ggplot2::element_blank(),
                     # panel.border = element_rect(color='black', linetype = "solid", 
                     #                             fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"), 
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position)
    
    if(length(model.names)==1){
      
      plotdf.plot = plotdf.plot +
        ggplot2::geom_point(shape = pointshape, size = pt.size, stroke = pointstroke,
                            position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim),
                               width=.5, linewidth = linewidth, alpha = opacity,
                               position = ggplot2::position_dodge(width = 1))
      
    }else{
      
      plotdf.plot = plotdf.plot +
        ggplot2::geom_point(ggplot2::aes(shape = method), size = pt.size, 
                            stroke = pointstroke,
                            position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim),
                               width=.5, linewidth = linewidth, alpha = opacity,
                               position = ggplot2::position_dodge(width = 1))
      
    }
    
    if(!missing(nrow.legend)){
      
      plotdf.plot = plotdf.plot +
        ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
                                                      override.aes = list(size = 10)),
                        shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
                                                      override.aes = list(size = 10)))
      
    }
    
    if(length(model.names)>1){
      
      names(pointshape) = model.names
      plotdf.plot = plotdf.plot + 
        ggplot2::scale_shape_manual(values = pointshape,
                                    labels = method.label)
      
    }
    
    if(!missing(pointcolors)){
      
      plotdf.plot = plotdf.plot + 
        ggplot2::scale_color_manual(values = pointcolors,
                                    labels = site.label) +
        ggplot2::scale_x_discrete(labels = site.label)
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='pull'){
    
    ## pull ----
    
    if(!missing(output.list)){
      
      # temporary model names for each output
      model.names = paste('model', 1:length(output.list), sep = '')
      
      plotdf = do.call('rbind.data.frame',
                       lapply(X = 1:length(output.list),
                              FUN = function(X){
                                
                                postsumm = apply(output.list[[X]]$MCMCout$musij_by_mits, 2:3,
                                                 FUN = function(v){
                                                   
                                                   c(mean(v),
                                                     quantile(x = v, probs = c(.025, .975), na.rm = T))
                                                   
                                                 })
                                
                                data.frame('est' = as.numeric(postsumm[1,,]),
                                           'llim' = as.numeric(postsumm[2,,]),
                                           'ulim' = as.numeric(postsumm[3,,]),
                                           'mits' = rep(rep(output.list[[X]]$datalist$causes, 
                                                            each = length(output.list[[X]]$datalist$sites)),
                                                        length(output.list[[X]]$datalist$causes)),
                                           'va' = rep(output.list[[X]]$datalist$causes, 
                                                      each = length(output.list[[X]]$datalist$causes)*
                                                        length(output.list[[X]]$datalist$sites)),
                                           'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                            length(output.list[[X]]$datalist$causes)),
                                                        length(output.list[[X]]$datalist$causes)),
                                           'method' = model.names[X])
                                
                              }))
      
      names(cause.label) = output.list[[1]]$datalist$causes
      
    }else{
      
      # data frame for ggplot
      plotdf = NULL
      
      # temporary model names for each output
      model.names = NULL
      
    }
    
    # adding empirical data
    if(!missing(output.empirical)){
      
      # data frame for ggplot
      plotdf = rbind.data.frame(data.frame('est' = as.numeric(output.empirical$musij_by_mits[1,,]),
                                           'llim' = as.numeric(output.empirical$musij_by_mits[2,,]),
                                           'ulim' = as.numeric(output.empirical$musij_by_mits[3,,]),
                                           'mits' = rep(rep(output.empirical$datalist$causes, 
                                                            each = length(output.empirical$datalist$sites)),
                                                        length(output.empirical$datalist$causes)),
                                           'va' = rep(output.empirical$datalist$causes, 
                                                      each = length(output.empirical$datalist$causes)*
                                                        length(output.empirical$datalist$sites)),
                                           'site' = rep(rep(output.empirical$datalist$sites,
                                                            length(output.empirical$datalist$causes)),
                                                        length(output.empirical$datalist$causes)),
                                           'method' = 'empirical'),
                                plotdf)
      
      if(is.null(names(cause.label))){
        
        names(cause.label) = output.empirical$datalist$causes
        
      }
      
      # temporary model names for each output
      model.names = c('empirical', model.names)
      
    }
    
    # making sure of the country and method order for plotting
    plotdf$mits = factor(x = plotdf$mits, levels = output.empirical$datalist$causes)  # alphabetical country order
    plotdf$va = factor(x = plotdf$va, levels = output.empirical$datalist$causes)  # alphabetical country order
    plotdf$site = factor(x = plotdf$site, levels = output.empirical$datalist$sites)  # alphabetical country order
    plotdf$method = factor(x = plotdf$method, levels = model.names)  # method order empirical (if given) then other outputs in their input order
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = subset(plotdf, mits!=va),
                      ggplot2::aes(x = va, y = est#, 
                                   # group = site, 
                                   # color = site
                      )) +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(.~va,
                          labeller = ggplot2::labeller(va = cause.label)) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=60, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=50),
                     axis.title.x = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = 40, 
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = 40),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     # axis.ticks.x = ggplot2::element_blank(),
                     # panel.border = element_rect(color='black', linetype = "solid", 
                     #                             fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = 40, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = 40, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"), 
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position)
    
    if(length(model.names)==1){
      
      plotdf.plot = plotdf.plot +
        ggplot2::geom_point(shape = pointshape, size = pt.size, stroke = pointstroke,
                            position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim),
                               width=.5, linewidth = linewidth, alpha = opacity,
                               position = ggplot2::position_dodge(width = 1))
      
    }else{
      
      plotdf.plot = plotdf.plot +
        ggplot2::geom_point(ggplot2::aes(shape = method), size = pt.size, 
                            stroke = pointstroke,
                            position = ggplot2::position_dodge(width = 1)) +
        ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim),
                               width=.5, linewidth = linewidth, alpha = opacity,
                               position = ggplot2::position_dodge(width = 1))
      
    }
    
    if(!missing(nrow.legend)){
      
      plotdf.plot = plotdf.plot +
        ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
                                                      override.aes = list(size = 10)),
                        shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
                                                      override.aes = list(size = 10)))
      
    }
    
    if(!missing(pointcolors)){
      
      plotdf.plot = plotdf.plot + 
        ggplot2::scale_color_manual(values = pointcolors,
                                    labels = site.label) +
        ggplot2::scale_x_discrete(labels = cause.label)
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='effsize'){
    
    ## pss ----
    
    # temporary model names for each output
    model.names = names(output.list)
    
    plotdf = do.call('rbind.data.frame',
                     lapply(X = 1:length(output.list),
                            FUN = function(X){
                              
                              pss.mat = matrix(nrow = output.list[[X]]$input$nMCMC,
                                               ncol = 3)
                              colnames(pss.mat) = c("pull", "hetsens", "hetrfp")
                              pss.mat[,1:ncol(output.list[[X]]$MCMCout$pss)] = as.matrix(output.list[[X]]$MCMCout$pss)
                              
                              return(data.frame('pss' = as.numeric(pss.mat),
                                                'pss_type' = rep(colnames(pss.mat),
                                                                 each = output.list[[X]]$input$nMCMC),
                                                'method' = model.names[X]))
                              
                            }))
    
    plotdf = rbind.data.frame(plotdf,
                              data.frame('pss' = as.numeric(output.list[[1]]$datalist$errormat_by_site),
                                         'pss_type' = 'observed',
                                         'method' = 'observed'))
    model.names = c('observed', model.names)
    
    pss_types = c("pull", "hetsens", "hetrfp", 'observed')
    
    # making sure of the method order for plotting
    plotdf$pss_type = factor(x = plotdf$pss_type, pss_types)  # method order empirical (if given) then other outputs in their input order
    plotdf$method = factor(x = plotdf$method, levels = model.names)  # method order empirical (if given) then other outputs in their input order
    
    if(is.null(pss.label)){
      
      pss.label = pss_types
      names(pss.label) = pss_types
      
    }
    
    if(is.null(method.label)){
      
      method.label = model.names
      names(method.label) = model.names
      
    }
    
    plotdf$pss_type.method = paste0(plotdf$pss_type, '_', plotdf$method)
    pss_type.method.unique = unique(plotdf$pss_type.method)
    
    plotdf.mean = do.call('rbind.data.frame',
                          lapply(X = 1:length(pss_type.method.unique),
                                 FUN = function(X){
                                   
                                   plotdf_X = subset(plotdf, pss_type.method==pss_type.method.unique[X])
                                   
                                   return(data.frame('pss_mean' = mean(plotdf_X$pss),
                                                     'pss_type' = plotdf_X$pss_type[1],
                                                     'method' = plotdf_X$method[1]))
                                   
                                 }))
    plotdf.mean$pss_type = factor(x = plotdf.mean$pss_type, levels = levels(plotdf$pss_type))  # method order empirical (if given) then other outputs in their input order
    plotdf.mean$method = factor(x = plotdf.mean$method, levels = levels(plotdf$method))  # method order empirical (if given) then other outputs in their input order
    
    if(missing(pointshape)){
      
      pointshape = c(8, 15 + 0:(length(output.empirical)-1))
      names(pointshape) = levels(model.names)
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = plotdf) +
      # ggplot2::facet_grid(.~pss_type,
      #                     labeller = ggplot2::labeller(pss_type = pss.label)) +
      ggplot2::facet_wrap(.~pss_type, scales = "free_x", nrow = 1, ncol = length(pss_types),
                          labeller = ggplot2::labeller(pss_type = pss.label),
                          # labeller = label_parsed
                          # labeller = ggplot2::as_labeller(x = pss.label,
                          #                                 default = label_parsed
                          #                                 )
                          ) +
      ggplot2::geom_violin(ggplot2::aes(x = method, y = pss, fill = method),
                           trim = T, alpha = opacity, linewidth = linewidth,
                           position = position_dodge()) + 
      ggplot2::scale_fill_manual(values = method.colors,
                                 labels = method.label) +
      ggplot2::scale_x_discrete(labels = method.label) +
      ggplot2::coord_cartesian(ylim = yrange.pss) +
      ggplot2::scale_y_continuous(
        trans = scales::sqrt_trans(),
        breaks = mybreaks,
        labels = mybreaks
      ) +
      ggplot2::geom_point(data = plotdf.mean,
                          ggplot2::aes(x = method, y = pss_mean, shape = method),
                          stroke = pointstroke, size = pt.size, color = "black") +
      ggplot2::scale_shape_manual(values = pointshape,
                                  labels = method.label) +
      ggplot2::guides(
        color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE),
        shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE),
      ) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=plot.title.size, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=plot.subtitle.size),
                     axis.title.x = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     # axis.ticks.x = element_line(linewidth = 1.5),
                     axis.ticks.length.x = unit(.5, "cm"),
                     axis.ticks.y = element_line(linewidth = 1.5),
                     axis.ticks.length.y = unit(.5, "cm"),
                     axis.ticks.x = ggplot2::element_blank(),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 1, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 1, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(3, "cm"), 
                     legend.key.size = ggplot2::unit(legend.key, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=legend.text.size),
                     legend.key.spacing.x = ggplot2::unit(2, 'cm'),
                     legend.position = legend.position) +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='se_heatmap'){
    
    ## se as heatmap ----
    
    # temporary model names for each output
    model.names = paste('model', 1:length(output.list), sep = '')
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    plotdf = do.call('rbind.data.frame',
                     lapply(X = 1:length(output.list),
                            FUN = function(X){
                              
                              postsumm = apply(output.list[[X]]$MCMCout$musij_by_mits, 
                                               2:3, mean)
                              
                              sqpostloss = (postsumm - output.empirical$musij_by_mits[1,,])
                              
                              rownames(sqpostloss) = paste0(rep(output.list[[X]]$datalist$causes,
                                                                each = length(output.list[[X]]$datalist$sites)), '_',
                                                            rep(output.list[[X]]$datalist$sites,
                                                                length(output.list[[X]]$datalist$causes)))
                              colnames(sqpostloss) = output.list[[X]]$datalist$causes
                              
                              sqpostloss_melt = reshape2::melt(sqpostloss)
                              
                              cbind.data.frame(sqpostloss_melt,
                                               'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                                length(output.list[[X]]$datalist$causes)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'mits' = rep(rep(output.list[[X]]$datalist$causes,
                                                                each = length(output.list[[X]]$datalist$sites)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'method' = model.names[X])
                              
                            }))
    
    plotdf$Var1 = factor(x = plotdf$Var1,
                         levels = rev(paste0(rep(output.list[[1]]$datalist$causes,
                                                 each = length(output.list[[1]]$datalist$sites)), '_',
                                             rep(output.list[[1]]$datalist$sites,
                                                 length(output.list[[1]]$datalist$causes)))))
    plotdf$Var2 = factor(x = plotdf$Var2,
                         levels = output.list[[1]]$datalist$causes)
    plotdf$site = factor(x = plotdf$site,
                         levels = output.list[[1]]$datalist$sites)
    plotdf$mits = factor(x = plotdf$mits,
                         levels = rev(output.list[[1]]$datalist$causes))
    plotdf$method = factor(x = plotdf$method, levels = model.names)
    
    # for manual cause label
    if(missing(cause.label)) cause.label = levels(plotdf$mits)
    if(missing(site.label)) site.label = levels(plotdf$site)
    site.freq = unlist(lapply(1:output.list[[1]]$datalist$nSite,
                              FUN = function(s){
                                
                                sum(output.list[[1]]$datalist$mits.total_by_mits[output.list[[1]]$datalist$siteid_by_mits==s])
                                
                              }))
    site.label = paste0(site.label, ' (', site.freq, ')')
    
    names(cause.label) = output.list[[1]]$datalist$causes
    names(site.label) = output.list[[1]]$datalist$sites
    
    plotdf.plot = ggplot2::ggplot(plotdf, aes(Var2, mits, fill = value)) + 
      ggplot2::geom_tile(color="white", linewidth=.5) +
      # ggplot2::scale_fill_viridis_c(option = 'viridis', direction = -1) +
      # ggplot2::scale_fill_gradient(low="white", high="blue4") +
      ggplot2::scale_fill_gradient2(low="red", mid = 'white', high="blue",
                                    limits = c(limits.colorgradient[1], limits.colorgradient[2]),
                                    n.breaks = n.breaks.colorgradient) +
      ggplot2::facet_grid(method~site, #site~method,
                          labeller = ggplot2::labeller(method = method.label,
                                                       site = site.label)
      ) +
      ggplot2::scale_x_discrete(labels = cause.label) +
      ggplot2::scale_y_discrete(labels = cause.label) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=60, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=50),
                     axis.title.x = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size,
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     axis.ticks.y = element_line(linewidth = 2),
                     axis.ticks.length.y = unit(.5, "cm"),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = 40, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = 40, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     # legend.title = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 40, face = "bold"),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"),
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position) #+
    # ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)),
    #                 shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)))
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title, fill = NULL)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='postintscore_heatmap'){
    
    ## interval score of 95% credible interval as heatmap ----
    
    # temporary model names for each output
    model.names = paste('model', 1:length(output.list), sep = '')
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    plotdf = do.call('rbind.data.frame',
                     lapply(X = 1:length(output.list),
                            FUN = function(X){
                              
                              credI = apply(output.list[[X]]$MCMCout$musij_by_mits, 2:3,
                                            FUN = function(v){
                                              
                                              quantile(x = v, probs = c(.025, .975), na.rm = T)
                                              
                                            })
                              
                              postintscore = 
                                matrix(data = scoringutils::interval_score(true_values = as.numeric(output.empirical$musij_by_mits[1,,]),
                                                                           lower = as.numeric(credI[1,,]),
                                                                           upper = as.numeric(credI[2,,]),
                                                                           interval_range = 95),
                                       nrow = length(output.list[[X]]$datalist$causes)*
                                         length(output.list[[X]]$datalist$sites),
                                       ncol = length(output.list[[X]]$datalist$causes))
                              
                              rownames(postintscore) = paste0(rep(output.list[[X]]$datalist$causes,
                                                                  each = length(output.list[[X]]$datalist$sites)), '_',
                                                              rep(output.list[[X]]$datalist$sites,
                                                                  length(output.list[[X]]$datalist$causes)))
                              colnames(postintscore) = output.list[[X]]$datalist$causes
                              
                              postintscore_melt = reshape2::melt(postintscore)
                              
                              cbind.data.frame(postintscore_melt,
                                               'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                                length(output.list[[X]]$datalist$causes)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'mits' = rep(rep(output.list[[X]]$datalist$causes,
                                                                each = length(output.list[[X]]$datalist$sites)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'method' = model.names[X])
                              
                            }))
    
    plotdf$Var1 = factor(x = plotdf$Var1,
                         levels = rev(paste0(rep(output.list[[1]]$datalist$causes,
                                                 each = length(output.list[[1]]$datalist$sites)), '_',
                                             rep(output.list[[1]]$datalist$sites,
                                                 length(output.list[[1]]$datalist$causes)))))
    plotdf$Var2 = factor(x = plotdf$Var2,
                         levels = output.list[[1]]$datalist$causes)
    plotdf$site = factor(x = plotdf$site,
                         levels = output.list[[1]]$datalist$sites)
    plotdf$mits = factor(x = plotdf$mits,
                         levels = rev(output.list[[1]]$datalist$causes))
    plotdf$method = factor(x = plotdf$method, levels = model.names)
    
    # for manual cause label
    if(missing(cause.label)) cause.label = levels(plotdf$mits)
    if(missing(site.label)) site.label = levels(plotdf$site)
    
    names(cause.label) = output.list[[1]]$datalist$causes
    names(site.label) = output.list[[1]]$datalist$sites
    
    plotdf.plot = ggplot2::ggplot(plotdf, aes(Var2, mits, fill = value)) + 
      ggplot2::geom_tile(color="white", linewidth=.5) +
      # ggplot2::scale_fill_viridis_c(option = 'viridis', direction = -1) +
      # ggplot2::scale_fill_gradient(low="white", high="blue4") +
      ggplot2::scale_fill_gradient2(low="red", mid = 'white', high="blue",
                                    limits = c(limits.colorgradient[1], limits.colorgradient[2]),
                                    n.breaks = n.breaks.colorgradient) +
      ggplot2::facet_grid(method~site, #site~method,
                          labeller = ggplot2::labeller(method = method.label,
                                                       site = site.label)
      ) +
      ggplot2::scale_x_discrete(labels = cause.label) +
      ggplot2::scale_y_discrete(labels = cause.label) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=60, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=50),
                     axis.title.x = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size,
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     axis.ticks.y = element_line(linewidth = 2),
                     axis.ticks.length.y = unit(.5, "cm"),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = 40, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = 40, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     # legend.title = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 40, face = "bold"),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"),
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position) #+
    # ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)),
    #                 shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)))
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title, fill = NULL)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='ic'){
    
    ## ic ----
    
    # temporary model names for each output
    model.names = paste('model', 1:length(output.list), sep = '')
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    plotdf = do.call('rbind.data.frame',
                     lapply(X = 1:length(output.list),
                            FUN = function(X){
                              
                              rbind.data.frame(data.frame('ic' = output.list[[X]]$MCMCout$ic.df["waic_Estimate"],
                                                          'se' = output.list[[X]]$MCMCout$ic.df["waic_SE"],
                                                          'ic_type' = 'WAIC',
                                                          'method' = model.names[X]),
                                               data.frame('ic' = output.list[[X]]$MCMCout$ic.df["looic_Estimate"],
                                                          'se' = output.list[[X]]$MCMCout$ic.df["looic_SE"],
                                                          'ic_type' = 'LOO-IC',
                                                          'method' = model.names[X]))
                              
                            }))
    
    # making sure of the method order for plotting
    plotdf$ic_type = factor(x = plotdf$ic_type, levels = c('LOO-IC', 'WAIC'))
    plotdf$method = factor(x = plotdf$method, levels = model.names)  # method order empirical (if given) then other outputs in their input order
    
    if(missing(barcolors)) barcolors = scales::hue_pal()(length(model.names))
    names(barcolors) = model.names
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = plotdf,
                      ggplot2::aes(x = method, y = ic)) +
      ggplot2::facet_grid(cols = ggplot2::vars(ic_type)#,
                          # labeller = ggplot2::labeller(ic_type = ic.label)
      ) +
      ggplot2::scale_x_discrete(labels = method.label) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=plot.title.size, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=plot.subtitle.size),
                     axis.title.x = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     # axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size, 
                     #                                     angle = xtext.angle, hjust = xhjust, vjust = xvjust),
                     # axis.ticks.x = element_line(linewidth = 1),
                     # axis.ticks.length.x = unit(.5, "cm"),
                     axis.text.x = ggplot2::element_blank(),
                     axis.ticks.x = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(1.5, "cm"), 
                     legend.key.height = ggplot2::unit(1.5, "cm"), 
                     legend.key.size = ggplot2::unit(legend.key, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), 
                     legend.key.spacing.x = ggplot2::unit(1, 'cm'), 
                     legend.text=ggplot2::element_text(size=legend.text.size),
                     legend.position = legend.position)
    
    if(length(model.names)==1){
      
      plotdf.plot = plotdf.plot +
        # ggplot2::coord_cartesian(ylim = c(ic.llim, NA)) +
        geom_bar(fill = scales::hue_pal()(1), stat="identity", color="black",
                 width = .5,
                 linewidth = linewidth, alpha = opacity) #+
      # scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) #+
      # ggplot2::geom_point(shape = pointshape, size = pt.size, stroke = pointstroke,
      #                     position = ggplot2::position_dodge(width = 1)) +
      # ggplot2::geom_errorbar(ggplot2::aes(x = method, ymin = ic-se, ymax = ic+se),
      #                        width=.5, linewidth = linewidth, alpha = opacity,
      #                        position = ggplot2::position_dodge(width = 1))
      
    }else{
      
      plotdf.plot = plotdf.plot +
        # ggplot2::coord_cartesian(ylim = c(ic.llim, NA)) +
        geom_bar(ggplot2::aes(fill = method), stat="identity", color="black",
                 width = .5,
                 linewidth = linewidth, alpha = opacity) +
        scale_fill_manual(values = barcolors,
                          labels = method.label) #+
      # ggplot2::geom_point(ggplot2::aes(shape = method), size = pt.size, 
      #                     stroke = pointstroke,
      #                     position = ggplot2::position_dodge(width = 1)) +
      # ggplot2::geom_errorbar(ggplot2::aes(x = method, ymin = ic-se, ymax = ic+se),
      #                        width=.5, linewidth = linewidth, alpha = opacity,
      #                        position = ggplot2::position_dodge(width = 1))
      
    }
    
    if(!missing(nrow.legend)){
      
      plotdf.plot = plotdf.plot +
        ggplot2::guides(fill = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                                     # override.aes = list(
                                                     #   linetype = "solid",
                                                     #   shape = c(16, rep(NA, length(model.names)-1)))
        ))
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='mu_pred'){
    
    ## predictive misclassification rates ----
    
    plotdf = cause.level = #site.level = 
      model.pred.type = NULL
    
    # temporary model names for each output
    model.names = paste0('model', 1:length(output.list))
    
    for(i in 1:length(output.list)){
      
      predsumm = apply(output.list[[i]]$m_pred, 2:3,
                       FUN = function(v){
                         
                         c(mean(v),
                           quantile(x = v, probs = c(.025, .975), na.rm = T))
                         
                       })
      
      plotdf = rbind.data.frame(plotdf,
                                data.frame('value' = as.numeric(predsumm[1,,]),
                                           'llim' = as.numeric(predsumm[2,,]),
                                           'ulim' = as.numeric(predsumm[3,,]),
                                           'mits' = rep(rep(output.list[[i]]$datalist$causes, 
                                                            output.list[[i]]$datalist$nSite),
                                                        output.list[[i]]$datalist$nCause),
                                           'va' = rep(output.list[[i]]$datalist$causes, 
                                                      each = output.list[[i]]$datalist$nCause*
                                                        output.list[[i]]$datalist$nSite),
                                           'site' = rep(rep(output.list[[i]]$datalist$sites,
                                                            each = output.list[[i]]$datalist$nCause),
                                                        output.list[[i]]$datalist$nCause),
                                           'method' = model.names[i],
                                           'pred.type' = output.list[[i]]$pred.type))
      
      model.pred.type = c(model.pred.type, output.list[[i]]$pred.type)
      
    }
    
    if(!missing(output.empirical)){
      
      plotdf = rbind.data.frame(plotdf,
                                data.frame('value' = as.numeric(output.empirical$MCMCout$m_insample[1,,]),
                                           'llim' = as.numeric(output.empirical$MCMCout$m_insample[2,,]),
                                           'ulim' = as.numeric(output.empirical$MCMCout$m_insample[3,,]),
                                           'mits' = rep(rep(output.list[[i]]$datalist$causes, 
                                                            output.list[[i]]$datalist$nSite),
                                                        output.list[[i]]$datalist$nCause),
                                           'va' = rep(output.list[[i]]$datalist$causes, 
                                                      each = output.list[[i]]$datalist$nCause*
                                                        output.list[[i]]$datalist$nSite),
                                           'site' = rep(rep(output.list[[i]]$datalist$sites,
                                                            each = output.list[[i]]$datalist$nCause),
                                                        output.list[[i]]$datalist$nCause),
                                           'method' = 'empirical',
                                           'pred.type' = output.empirical$MCMCout$est.type))
      
      model.names = c(model.names, 'empirical')
      model.pred.type = c(model.pred.type, output.empirical$MCMCout$est.type)
      
    }
    
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    # print(model.names)
    # print(model.est.type)
    
    plotdf$linestyle = plotdf$shapestyle = 
      paste0(plotdf$method, '_', plotdf$pred.type)
    plotdf$linestyle[plotdf$pred.type=='heterogen'] = NA
    plotdf$shapestyle[plotdf$pred.type=='homogen'] = NA
    
    alllabel = data.frame('model.names' = model.names, 'method.label' = method.label,
                          'style' = paste0(model.names, '_', model.pred.type),
                          'pred.type' = model.pred.type)
    
    print(alllabel)
    
    if(is.null(cause.level)) cause.level = output.list[[1]]$datalist$causes
    if(is.null(site.level)) site.level = output.list[[1]]$datalist$sites
    
    # making sure of the country and method order for plotting
    plotdf$mits = factor(x = plotdf$mits, levels = cause.level)
    plotdf$va = factor(x = plotdf$va, levels = cause.level)
    plotdf$site = factor(x = plotdf$site, levels = site.level)
    plotdf$method = factor(x = plotdf$method, levels = model.names)
    plotdf$linestyle = factor(x = plotdf$linestyle)
    plotdf$shapestyle = factor(x = plotdf$shapestyle)
    
    print(levels(plotdf$linestyle))
    print(levels(plotdf$shapestyle))
    
    if(is.null(cause.label)) cause.label = cause.level
    mits.freq = unlist(lapply(1:output.list[[1]]$datalist$nCause,
                              FUN = function(i){
                                
                                sum(output.list[[1]]$datalist$mits.total_by_site[output.list[[1]]$datalist$mitsid_by_site==i])
                                
                              }))
    va.freq = unname(colSums(output.list[[1]]$datalist$errormat_by_site))
    
    cause.label.row = paste0(cause.label, ' (', mits.freq, ')')
    cause.label.col = paste0(cause.label, ' (', va.freq, ')')
    names(cause.label.row) = names(cause.label.col) = cause.level
    
    if(is.null(site.label)) site.label = site.level
    # names(site.label) = site.level
    
    idmodel.homogen = model.pred.type=='homogen'
    
    # mu ggplot object
    plotdf.plot =
      ggplot2::ggplot(data = plotdf,
                      ggplot2::aes(x = site, y = value)) +
      ggplot2::ylim(0,1) +
      ggplot2::facet_grid(mits~va,
                          labeller = ggplot2::labeller(mits = cause.label.row,
                                                       va = cause.label.col)) +
      ggplot2::geom_rect(data = subset(plotdf, mits==va), 
                         fill = NA, colour = "red3", linewidth = 8,
                         xmin = -Inf, xmax = Inf,
                         ymin = -Inf, ymax = Inf) +
      ggplot2::geom_rect(data = subset(plotdf, mits!=va), 
                         fill = NA, colour = "black", linewidth = 3,
                         xmin = -Inf, xmax = Inf,
                         ymin = -Inf, ymax = Inf) +
      ggplot2::scale_x_discrete(labels = site.label) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=plot.title.size, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=plot.subtitle.size),
                     axis.title.x = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=axis.title.size,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size, 
                                                         angle = xtext.angle, hjust = xhjust, vjust = xvjust),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 1),
                     axis.ticks.length.x = unit(.5, "cm"),
                     # axis.ticks.x = ggplot2::element_blank(),
                     # panel.border = element_rect(color='black', linetype = "solid", 
                     #                             fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = strip.text.size, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     legend.title = ggplot2::element_blank(),
                     legend.key.width = ggplot2::unit(legend.key.width, "cm"), 
                     legend.key.height = ggplot2::unit(legend.key.height, "cm"), 
                     legend.key.size = ggplot2::unit(legend.key, "cm"),
                     legend.spacing.x = ggplot2::unit(legend.spacing.x, 'cm'), 
                     legend.spacing.y = ggplot2::unit(legend.spacing.y, 'cm'), 
                     legend.text=ggplot2::element_text(size=legend.text.size),
                     legend.position = legend.position)
    
    if(length(model.names)==1){
      
      if(model.pred.type=='homogen'){
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_hline(data = subset(plotdf, pred.type=='homogen'),
                              ggplot2::aes(yintercept = value),
                              linetype = 'solid', color = 'black', linewidth = linewidth) +
          ggplot2::geom_ribbon(data = subset(plotdf, pred.type=='homogen'),
                               aes(x = as.numeric(site),
                                   # xmin = as.numeric(site)-1, xmax = as.numeric(site)+1, 
                                   ymin = llim, ymax = ulim),
                               fill = 'grey50', alpha = opacity)
        
      }
      
      if(model.pred.type=='heterogen'){
        
        if(is.null(pointcolors)) pointcolors = scales::hue_pal()(length(site.level))
        # names(pointcolors) = site.level
        
        if(is.null(pointshape)) pointshape = 8
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_point(data = subset(plotdf, pred.type=='heterogen'),
                              ggplot2::aes(x = site, y = value, group = method, color = site),
                              shape = pointshape, size = pt.size, stroke = pointstroke,
                              position = ggplot2::position_dodge(width = 1)) +
          ggplot2::geom_errorbar(ggplot2::aes(x = site, ymin = llim, ymax = ulim, color = site),
                                 width=.5, linewidth = linewidth, alpha = opacity,
                                 position = ggplot2::position_dodge(width = 1)) +
          ggplot2::scale_color_manual(values = pointcolors,
                                      labels = site.label)
        
      }
      
    }else{
      
      if(sum(!idmodel.homogen)>0){
        
        if(is.null(pointcolors)) pointcolors = scales::hue_pal()(length(site.level))
        # names(pointcolors) = site.level
        
        # alllabel = data.frame('model.names' = model.names, 'method.label' = method.label,
        #                       'style' = paste0(model.names, '_', model.pred.type))
        
        if(is.null(pointshape)){
          
          pointshape = 8 + 0:(length(levels(plotdf$shapestyle))-1)
          
        }else{
          
          pointshape = pointshape[match(levels(plotdf$shapestyle),
                                        alllabel$style)]
          
        }
        names(pointshape) = levels(plotdf$shapestyle)
        
        shape.label = alllabel$method.label[match(levels(plotdf$shapestyle),
                                                  alllabel$style)]
        names(shape.label) = names(pointshape)
        
        print(pointshape)
        print(shape.label)
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_point(data = subset(plotdf, pred.type=='heterogen'),
                              ggplot2::aes(x = site, y = value, 
                                           group = shapestyle, color = site, shape = shapestyle),
                              size = pt.size, stroke = pointstroke,
                              position = ggplot2::position_dodge(width = 1)) +
          ggplot2::geom_errorbar(data = subset(plotdf, pred.type=='heterogen'),
                                 ggplot2::aes(x = site, ymin = llim, ymax = ulim, 
                                              group = shapestyle, color = site),
                                 width=.5, linewidth = linewidth, alpha = opacity,
                                 position = ggplot2::position_dodge(width = 1)) +
          ggplot2::scale_color_manual(values = pointcolors,
                                      labels = site.label) +
          ggplot2::scale_shape_manual(values = pointshape,
                                      labels = shape.label)
        
      }
      
      if(sum(idmodel.homogen)>0){
        
        if(is.null(linetype)){
          
          linetype = 1:length(levels(plotdf$linestyle))
          
        }else{
          
          linetype = linetype[match(levels(plotdf$linestyle),
                                    alllabel$style)]
          
        }
        names(linetype) = levels(plotdf$linestyle)
        
        linetype.label = alllabel$method.label[match(levels(plotdf$linestyle),
                                                     alllabel$style)]
        names(linetype.label) = names(linetype)
        
        plotdf.plot = plotdf.plot +
          ggplot2::geom_hline(data = subset(plotdf, pred.type=='homogen'),
                              ggplot2::aes(yintercept = value, linetype = linestyle),
                              color = 'black', linewidth = linewidth) +
          geom_ribbon(data = subset(plotdf, pred.type=='homogen'),
                      aes(x = as.numeric(site), ymin = llim, ymax = ulim),
                      fill = 'grey50', alpha = opacity) +
          ggplot2::scale_linetype_manual(values = linetype,
                                         labels = linetype.label)
        
      }
      
    }
    
    if(!missing(nrow.legend)){
      
      plotdf.plot = plotdf.plot +
        ggplot2::guides(color = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                                      # override.aes = list(
                                                      #   linetype = "solid",
                                                      #   shape = c(16, rep(NA, length(model.names)-1)))
        ),
        shape = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                      # override.aes = list(
                                      #   linetype = "solid",
                                      #   shape = c(16, rep(NA, length(model.names)-1)))
        ),
        linetype = ggplot2::guide_legend(nrow = nrow.legend, byrow=FALSE#,
                                         # override.aes = list(
                                         #   linetype = "solid",
                                         #   shape = c(16, rep(NA, length(model.names)-1)))
        )) #+
      # ggplot2::guides(linetype = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
      #                                               override.aes = list(
      #                                                 linetype = linetype_toplot
      #                                                 )),
      #                 shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
      #                                               override.aes = list(
      #                                                 shape = shape_toplot
      #                                               )),
      #                 color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
      #                                               override.aes = list(
      #                                                 color = color_toplot
      #                                               ))
      #                 )
      
    }
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='spe_heatmap'){
    
    ## spe as heatmap ----
    
    # temporary model names for each output
    model.names = paste('model', 1:length(output.list), sep = '')
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    plotdf = do.call('rbind.data.frame',
                     lapply(X = 1:length(output.list),
                            FUN = function(X){
                              
                              predsumm = apply(output.list[[X]]$muij_pred_by_mits, 
                                               2:3, mean)
                              
                              sqpredloss = (predsumm - output.empirical$musij_by_mits[1,,])
                              
                              rownames(sqpredloss) = paste0(rep(output.list[[X]]$datalist$causes,
                                                                each = length(output.list[[X]]$datalist$sites)), '_',
                                                            rep(output.list[[X]]$datalist$sites,
                                                                length(output.list[[X]]$datalist$causes)))
                              colnames(sqpredloss) = output.list[[X]]$datalist$causes
                              
                              sqpredloss_melt = reshape2::melt(sqpredloss)
                              
                              cbind.data.frame(sqpredloss_melt,
                                               'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                                length(output.list[[X]]$datalist$causes)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'mits' = rep(rep(output.list[[X]]$datalist$causes,
                                                                each = length(output.list[[X]]$datalist$sites)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'method' = model.names[X])
                              
                            }))
    
    plotdf$Var1 = factor(x = plotdf$Var1,
                         levels = rev(paste0(rep(output.list[[1]]$datalist$causes,
                                                 each = length(output.list[[1]]$datalist$sites)), '_',
                                             rep(output.list[[1]]$datalist$sites,
                                                 length(output.list[[1]]$datalist$causes)))))
    plotdf$Var2 = factor(x = plotdf$Var2,
                         levels = output.list[[1]]$datalist$causes)
    plotdf$site = factor(x = plotdf$site,
                         levels = output.list[[1]]$datalist$sites)
    plotdf$mits = factor(x = plotdf$mits,
                         levels = rev(output.list[[1]]$datalist$causes))
    plotdf$method = factor(x = plotdf$method, levels = model.names)
    
    # for manual cause label
    if(missing(cause.label)) cause.label = levels(plotdf$mits)
    if(missing(site.label)) site.label = levels(plotdf$site)
    
    names(cause.label) = output.list[[1]]$datalist$causes
    names(site.label) = output.list[[1]]$datalist$sites
    
    plotdf.plot = ggplot2::ggplot(plotdf, aes(Var2, mits, fill = value)) + 
      ggplot2::geom_tile(color="white", linewidth=.5) +
      # ggplot2::scale_fill_viridis_c(option = 'viridis', direction = -1) +
      # ggplot2::scale_fill_gradient(low="white", high="blue4") +
      ggplot2::scale_fill_gradient2(low="red", mid = 'white', high="blue",
                                    limits = c(limits.colorgradient[1], limits.colorgradient[2]),
                                    n.breaks = n.breaks.colorgradient) +
      ggplot2::facet_grid(method~site, #site~method,
                          labeller = ggplot2::labeller(method = method.label,
                                                       site = site.label)
      ) +
      ggplot2::scale_x_discrete(labels = cause.label) +
      ggplot2::scale_y_discrete(labels = cause.label) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=60, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=50),
                     axis.title.x = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size,
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     axis.ticks.y = element_line(linewidth = 2),
                     axis.ticks.length.y = unit(.5, "cm"),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = 40, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = 40, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     # legend.title = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 40, face = "bold"),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"),
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position) #+
    # ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)),
    #                 shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)))
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title, fill = NULL)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }else if(toplot=='predintscore_heatmap'){
    
    ## interval score of 95% prediction interval as heatmap ----
    
    # temporary model names for each output
    model.names = paste('model', 1:length(output.list), sep = '')
    if(is.null(method.label)) method.label = model.names
    names(method.label) = model.names
    
    plotdf = do.call('rbind.data.frame',
                     lapply(X = 1:length(output.list),
                            FUN = function(X){
                              
                              predI = apply(output.list[[X]]$muij_pred_by_mits, 2:3,
                                            FUN = function(v){
                                              
                                              quantile(x = v, probs = c(.025, .975), na.rm = T)
                                              
                                            })
                              
                              predintscore = 
                                matrix(data = scoringutils::interval_score(true_values = as.numeric(output.empirical$musij_by_mits[1,,]),
                                                                           lower = as.numeric(predI[1,,]),
                                                                           upper = as.numeric(predI[2,,]),
                                                                           interval_range = 95),
                                       nrow = length(output.list[[X]]$datalist$causes)*
                                         length(output.list[[X]]$datalist$sites),
                                       ncol = length(output.list[[X]]$datalist$causes))
                              
                              rownames(predintscore) = paste0(rep(output.list[[X]]$datalist$causes,
                                                                  each = length(output.list[[X]]$datalist$sites)), '_',
                                                              rep(output.list[[X]]$datalist$sites,
                                                                  length(output.list[[X]]$datalist$causes)))
                              colnames(predintscore) = output.list[[X]]$datalist$causes
                              
                              predintscore_melt = reshape2::melt(predintscore)
                              
                              cbind.data.frame(predintscore_melt,
                                               'site' = rep(rep(output.list[[X]]$datalist$sites,
                                                                length(output.list[[X]]$datalist$causes)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'mits' = rep(rep(output.list[[X]]$datalist$causes,
                                                                each = length(output.list[[X]]$datalist$sites)),
                                                            length(output.list[[X]]$datalist$causes)),
                                               'method' = model.names[X])
                              
                            }))
    
    plotdf$Var1 = factor(x = plotdf$Var1,
                         levels = rev(paste0(rep(output.list[[1]]$datalist$causes,
                                                 each = length(output.list[[1]]$datalist$sites)), '_',
                                             rep(output.list[[1]]$datalist$sites,
                                                 length(output.list[[1]]$datalist$causes)))))
    plotdf$Var2 = factor(x = plotdf$Var2,
                         levels = output.list[[1]]$datalist$causes)
    plotdf$site = factor(x = plotdf$site,
                         levels = output.list[[1]]$datalist$sites)
    plotdf$mits = factor(x = plotdf$mits,
                         levels = rev(output.list[[1]]$datalist$causes))
    plotdf$method = factor(x = plotdf$method, levels = model.names)
    
    # for manual cause label
    if(missing(cause.label)) cause.label = levels(plotdf$mits)
    if(missing(site.label)) site.label = levels(plotdf$site)
    
    names(cause.label) = output.list[[1]]$datalist$causes
    names(site.label) = output.list[[1]]$datalist$sites
    
    plotdf.plot = ggplot2::ggplot(plotdf, aes(Var2, mits, fill = value)) + 
      ggplot2::geom_tile(color="white", linewidth=.5) +
      # ggplot2::scale_fill_viridis_c(option = 'viridis', direction = -1) +
      # ggplot2::scale_fill_gradient(low="white", high="blue4") +
      ggplot2::scale_fill_gradient2(low="red", mid = 'white', high="blue",
                                    limits = c(limits.colorgradient[1], limits.colorgradient[2]),
                                    n.breaks = n.breaks.colorgradient) +
      ggplot2::facet_grid(method~site, #site~method,
                          labeller = ggplot2::labeller(method = method.label,
                                                       site = site.label)
      ) +
      ggplot2::scale_x_discrete(labels = cause.label) +
      ggplot2::scale_y_discrete(labels = cause.label) +
      ggplot2::theme(plot.title = ggplot2::element_text(size=60, face="bold"),
                     plot.subtitle = ggplot2::element_text(size=50),
                     axis.title.x = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = ggplot2::element_text(size=50,
                                                          margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
                     axis.text.x = ggplot2::element_text(color = "black", size = axis.text.size,
                                                         angle = xtext.angle, hjust = xhjust, vjust = 1),
                     axis.text.y = ggplot2::element_text(color = "black", size = axis.text.size),
                     axis.ticks.x = element_line(linewidth = 2),
                     axis.ticks.length.x = unit(.5, "cm"),
                     axis.ticks.y = element_line(linewidth = 2),
                     axis.ticks.length.y = unit(.5, "cm"),
                     panel.border = element_rect(color='black', linetype = "solid",
                                                 fill = NA, linewidth = 2),
                     panel.background = ggplot2::element_blank(),
                     panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                                              colour = "grey90"),
                     strip.text.x = ggplot2::element_text(size = 40, face = "bold"),
                     strip.text.y = ggplot2::element_text(size = 40, face = "bold"),
                     strip.background = ggplot2::element_rect(color="black", linewidth=3),
                     # legend.title = ggplot2::element_blank(),
                     legend.title = ggplot2::element_text(size = 40, face = "bold"),
                     legend.key.width = ggplot2::unit(3, "cm"), legend.key.height = ggplot2::unit(2, "cm"),
                     legend.key.size = ggplot2::unit(30, "cm"),
                     legend.spacing.x = ggplot2::unit(2, 'cm'), legend.text=ggplot2::element_text(size=40),
                     legend.position = legend.position) #+
    # ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)),
    #                 shape = ggplot2::guide_legend(nrow=nrow.legend, byrow=FALSE,
    #                                               override.aes = list(size = 10)))
    
    if(missing(plot.title)) plot.title = ggplot2::waiver()
    if(missing(plot.subtitle)) plot.subtitle = ggplot2::waiver()
    if(missing(xlab.title)) xlab.title = ggplot2::waiver()
    if(missing(ylab.title)) ylab.title = ggplot2::waiver()
    
    plotdf.plot = plotdf.plot +
      ggplot2::labs(title = plot.title, subtitle = plot.subtitle,
                    x = xlab.title, y = ylab.title, fill = NULL)
    
    plotlist = c(plotlist, list(plotdf.plot))
    
  }
  
  return(plotlist)
  
}

