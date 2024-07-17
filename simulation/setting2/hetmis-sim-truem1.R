

rm(list = ls())


# sourcing sourcecodes, outputs ----
sourcecode.path = ...    # specifies path ".../sourcecode" to the sourcecode folder

# saving output ====
output_dir <- 'hetsim_output'
if(!dir.exists(output_dir)){
  dir.create(output_dir)
}


# sourcing misclassification modeling codes
source(file.path(sourcecode.path, 'hetmis-functions.R'))


# simulation specifics ----
nCause = 5  # number of causes in CHAMPS
allcauses = paste0('cause', 1:nCause)  # all cause names
nSites = 6  # number of sites in CHAMPS
allsites = paste0('site', 1:nSites)  # all site names
ai.true = seq(0.7-0.05, 0.7+0.05, length.out = nCause)  # intrinsic accuracy
alpha.true = c(0.1, 0.3, 0.2, 0.1, 0.3)  # pull

# true CSMF in each site
set.seed(1)
p.true = round(LaplacesDemon::rdirichlet(n = nSites, alpha = rep(1, nCause)), 2)
for(s in 1:nSites){
  
  p.true[s, which.max(p.true[s,])] = 1 - sum(p.true[s,-which.max(p.true[s,])])
  
}
rowSums(p.true)
rownames(p.true) = allsites
colnames(p.true) = allcauses
p.true

# other specifics
nMits = 100
nVA = 200
nReplicate = 50

# mcmc specifics
nBurn = 2000
nMCMC = 5000
nThin = 1

# true effect size
pss.true = cbind.data.frame('pull' = 100, 'hetsens' = 5)


# true missclassification matrix ----
musij_by_mits.true = NULL
for(i in 1:nCause){
  
  # country-specific sensitivity
  sensi_sp.true = ai.true[i] + (1-ai.true[i])*alpha.true[i]
  sensi.true = rbeta(n = 1,
                     shape1 = 0.5 + 2*pss.true$pull*sensi_sp.true, 
                     shape2 = 0.5 + 2*pss.true$pull*(1-sensi_sp.true))
  senssi.true = rbeta(n = nSites,
                      shape1 = 0.5 + 2*pss.true$hetsens*sensi.true, 
                      shape2 = 0.5 + 2*pss.true$hetsens*(1-sensi.true))
  
  # country-specific rfp
  qi.true = as.numeric(LaplacesDemon::rdirichlet(n = 1,
                                                 alpha = 0.5 + (nCause-1)*pss.true$pull*(alpha.true[-i]/(1-alpha.true[i]))))
  qsi.true = matrix(data = qi.true, nrow = nSites, 
                    ncol = nCause-1, byrow = T)
  
  # country-specific misclassification matrices
  if(i==1){
    
    musij_by_mits.true = rbind(musij_by_mits.true,
                               cbind(senssi.true, (1-senssi.true)*qsi.true))
    
  }else if(i==nCause){
    
    musij_by_mits.true = rbind(musij_by_mits.true,
                               cbind((1-senssi.true)*qsi.true, senssi.true))
    
  }else{
    
    musij_by_mits.true = rbind(musij_by_mits.true,
                               cbind((1-senssi.true)*qsi.true[,1:(i-1)],
                                     senssi.true, 
                                     (1-senssi.true)*qsi.true[,i:(nCause-1)]))
    
  }
  
}

musij_by_mits.true = round(musij_by_mits.true, 2)

for(l in 1:nrow(musij_by_mits.true)){
  
  musij_by_mits.true[l, which.max(musij_by_mits.true[l,])] =
    1 - sum(musij_by_mits.true[l, -which.max(musij_by_mits.true[l,])])
  
}

head(musij_by_mits.true)
range(rowSums(musij_by_mits.true))

siteid_by_mits = rep(1:length(allsites), nCause)

# arranging in an array
musij.true = array(dim = c(nSites, nCause, nCause),
                   dimnames = list(allsites, allcauses, allcauses))
for(s in 1:nSites){
  
  musij.true[s,,] = musij_by_mits.true[siteid_by_mits==s,]
  
  print(range(rowSums(musij.true[s,,])))
  
  print(allsites[s])
  
}

do.call("c",
        lapply(1:nSites,
               FUN = function(s){
                 
                 min(diag(musij.true[s,,]))
                 
               }))

do.call("rbind",
        lapply(1:nSites,
               FUN = function(s){
                 
                 sort(svd(musij.true[s,,])$d, decreasing = T)
                 
               }))


# true uncalibrated csmf ====
p.uncalib.true = do.call("rbind",
                         lapply(1:nSites,
                                FUN = function(s){
                                  
                                  as.numeric(crossprod(musij.true[s,,], p.true[s,]))
                                  
                                }))
rownames(p.uncalib.true) = allsites
colnames(p.uncalib.true) = allcauses
rowSums(p.uncalib.true)


# plotting true missclassification matrix ====
plotdf = data.frame('value' = as.numeric(musij_by_mits.true),
                    'mits' = rep(rep(allcauses, each = length(allsites)),
                                 nCause),
                    'va' = rep(allcauses, each = nCause*length(allsites)),
                    'site' = rep(rep(allsites, nCause), nCause))

plotdf$mits = factor(x = plotdf$mits, levels = allcauses)
plotdf$va = factor(x = plotdf$va, levels = allcauses)
plotdf$site = factor(x = plotdf$site, levels = allsites)

cause.label = paste0('Cause ', 1:nCause)
names(cause.label) = allcauses

site.label = paste0('Site ', 1:nSites)
names(site.label) = allsites

pointcolors = c('green4', 'black', 'chocolate',
                'blue', 'orchid', 'dodgerblue', 'red3')
names(pointcolors) = allsites

ggplot2::ggplot(data = plotdf) +
  ggplot2::facet_grid(
    mits~va,
    labeller = ggplot2::labeller(mits = cause.label,
                                 va = cause.label)
  ) +
  ggplot2::geom_rect(data = subset(plotdf, mits==va),
                     fill = NA, colour = "red3", linewidth = 1,
                     xmin = -Inf, xmax = Inf,
                     ymin = -Inf, ymax = Inf) +
  ggplot2::geom_rect(data = subset(plotdf, mits!=va),
                     fill = NA, colour = "black", linewidth = 1,
                     xmin = -Inf, xmax = Inf,
                     ymin = -Inf, ymax = Inf) +
  ggplot2::ylim(0,1) +
  ggplot2::geom_point(ggplot2::aes(x = site, y = value, color = site),
                      shape = 8, size = 3, stroke = .9) +
  ggplot2::scale_color_manual(values = pointcolors,
                              labels = site.label) +
  ggplot2::scale_x_discrete(labels = site.label) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(size=15,
                                         margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = ggplot2::element_text(size=15,
                                         margin = ggplot2::margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = ggplot2::element_text(color = "black", size = 10,
                                        angle = 45, hjust = 1, vjust = 1),
    axis.text.y = ggplot2::element_text(color = "black", size = 10),
    # axis.ticks.x = ggplot2::element_line(linewidth = 1),
    axis.ticks.length.x = ggplot2::unit(.15, "cm"),
    axis.ticks.length.y = ggplot2::unit(.15, "cm"),
    panel.border = ggplot2::element_rect(color='black', linetype = "solid",
                                         fill = NA, linewidth = 1),
    panel.background = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    panel.grid.minor = ggplot2::element_line(linewidth = 0.5, linetype = 'solid',
                                             colour = "grey90"),
    strip.text.x = ggplot2::element_text(size = 15, face = "bold"),
    strip.text.y = ggplot2::element_text(size = 15, face = "bold"),
    strip.background = ggplot2::element_rect(color="black", linewidth=1),
    legend.title = ggplot2::element_blank(),
    legend.key.width = ggplot2::unit(1, "cm"),
    legend.key.height = ggplot2::unit(1, "cm"),
    legend.key.spacing.x = ggplot2::unit(.5, 'cm'),
    legend.text=ggplot2::element_text(size=15),
    legend.position = 'bottom'
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(
    nrow = 1, byrow=FALSE#,
    # override.aes = list(
    #   linetype = "solid",
    #   shape = c(rep(16, length(levels(plotdf$method))-1), NA))
  )) +
  ggplot2::labs(y = 'Classification Probability',
                x = 'Countries')   # 8 X 10


# replication study ==== 
for(r in 1:nReplicate){
  
  # r = 1
  set.seed(r)
  
  
  ## observed missclassification matrix ====
  obserrormat_by_site = lapply(1:nSites,
                               FUN = function(s){
                                 
                                 errormat_s = t(apply(musij.true[s,,], 1,
                                                      FUN = function(v){
                                                        
                                                        as.numeric(rmultinom(n = 1, 
                                                                             size = nMits,
                                                                             prob = v))
                                                        
                                                      }))
                                 rownames(errormat_s) = colnames(errormat_s) = allcauses
                                 
                                 errormat_s
                                 
                               })
  names(obserrormat_by_site) = allsites[1:nSites]
  
  
  ## observed unlabeled data ====
  obsva = do.call("rbind",
                  lapply(1:nSites,
                         FUN = function(s){
                           
                           as.numeric(rmultinom(n = 1, size = nVA, prob = p.uncalib.true[s,]))
                           
                         }))
  rownames(obsva) = allsites
  colnames(obsva) = allcauses
  rowSums(obsva)
  
  
  ## in-sample  ====
  ### in-sample misclassification ====
  #### homogeneous ====
  fitout_m0 = hetmis(errormat_by_site = obserrormat_by_site,
                     model.choice = 'hom', cause.type = 'single',
                     pss.method = 'learn', shapeBase = .5, 
                     shape_pull = c(.5,.5),
                     sensprior_shape = c(1,1),
                     Dirprior_shape = rep(1, nCause),
                     nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                     sourcecode.path = sourcecode.path,
                     seed = 1, verbose = T, saveoutput = F)
  
  #### partly heterogeneous ====
  fitout_m1 = hetmis(errormat_by_site = obserrormat_by_site,
                     model.choice = 'het_part', cause.type = 'single',
                     pss.method = 'learn', shapeBase = .5, 
                     shape_pull = c(.5,.5), shape_hetsens = c(.5,.5),
                     sensprior_shape = c(1,1),
                     Dirprior_shape = rep(1, nCause),
                     nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                     sourcecode.path = sourcecode.path,
                     seed = 1, verbose = T, saveoutput = F)
  
  #### fully heterogeneous ====
  fitout_m2 = hetmis(errormat_by_site = obserrormat_by_site,
                     model.choice = 'het', cause.type = 'single',
                     pss.method = 'learn', shapeBase = .5, 
                     shape_pull = c(.5,.5), shape_hetsens = c(.5,.5), shape_hetrfp = c(.5,.5),
                     sensprior_shape = c(1,1),
                     Dirprior_shape = rep(1, nCause),
                     nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                     sourcecode.path = sourcecode.path,
                     seed = 1, verbose = T, saveoutput = F)
  
  
  # output array
  Mmat.insample = array(dim = c(nMCMC, nSites, nCause, nCause, 3),
                        dimnames = list(NULL, allsites, allcauses, allcauses,
                                        c(fitout_m0$input$model.choice,
                                          fitout_m1$input$model.choice,
                                          fitout_m2$input$model.choice)))
  Mmat.obs = array(dim = c(nSites, nCause, nCause),
                   dimnames = list(allsites, allcauses, allcauses))
  for(s in 1:nSites){
    
    Mmat.insample[,s,,,fitout_m0$input$model.choice] = fitout_m0$Mmat[[s]]
    Mmat.insample[,s,,,fitout_m1$input$model.choice] = fitout_m1$Mmat[[s]]
    Mmat.insample[,s,,,fitout_m2$input$model.choice] = fitout_m2$Mmat[[s]]
    
    Mmat.obs[s,,] = obserrormat_by_site[[s]]/rowSums(obserrormat_by_site[[s]])
    
    print(allsites[s])
    
  }
  
  ic.df = cbind(fitout_m0$MCMCout$ic.df, 
                fitout_m1$MCMCout$ic.df,
                fitout_m2$MCMCout$ic.df)
  colnames(ic.df) = c(fitout_m0$input$model.choice,
                      fitout_m1$input$model.choice,
                      fitout_m2$input$model.choice)
  
  
  ### calibration ====
  doParallel::registerDoParallel(cores = 6)
  p.insample = foreach::foreach(s = 1:nSites, .combine = "mysimcombine_insample", .multicombine = T) %dopar% {
    
    #### using homogeneous misclassification ====
    calibout_m0 = calibratedva_v2(va_unlabeled = obsva[s,],
                                  cause.type = 'single', model.choice = "seqpshrink_up",
                                  Mmat.asDirich = fitout_m0$Mmat.asDirich[[s]],
                                  nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                                  sourcecode.path = sourcecode.path,
                                  seed = 1, verbose = T, saveoutput = F)
    
    #### using partly heterogeneous misclassification ====
    calibout_m1 = calibratedva_v2(va_unlabeled = obsva[s,],
                                  cause.type = 'single', model.choice = "seqpshrink_up",
                                  Mmat.asDirich = fitout_m1$Mmat.asDirich[[s]],
                                  nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                                  sourcecode.path = sourcecode.path,
                                  seed = 1, verbose = T, saveoutput = F)
    
    #### using fully heterogeneous misclassification ====
    calibout_m2 = calibratedva_v2(va_unlabeled = obsva[s,],
                                  cause.type = 'single', model.choice = "seqpshrink_up",
                                  Mmat.asDirich = fitout_m2$Mmat.asDirich[[s]],
                                  nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                                  sourcecode.path = sourcecode.path,
                                  seed = 1, verbose = T, saveoutput = F)
    
    p_calib = array(dim = c(nMCMC, nCause, 3),
                    dimnames = list(NULL, allcauses, 
                                    c(fitout_m0$input$model.choice,
                                      fitout_m1$input$model.choice,
                                      fitout_m2$input$model.choice)))
    p_calib[,,fitout_m0$input$model.choice] = calibout_m0$MCMCout$p_calib
    p_calib[,,fitout_m1$input$model.choice] = calibout_m1$MCMCout$p_calib
    p_calib[,,fitout_m2$input$model.choice] = calibout_m2$MCMCout$p_calib
    
    print(s)
    
    list('uncalib' = calibout_m0$p_uncalib_obs, 'calib' = p_calib)
    
  }
  
  rownames(p.insample$uncalib) = dimnames(p.insample$calib)[[4]] = allsites
  
  # rearranging for convenience
  p_calib_temp = array(dim = c(nMCMC, nSites, nCause, 3),
                       dimnames = list(NULL, allsites, allcauses, 
                                       c(fitout_m0$input$model.choice,
                                         fitout_m1$input$model.choice,
                                         fitout_m2$input$model.choice)))
  for(s in 1:nSites){
    
    p_calib_temp[,s,,] = p.insample$calib[,,,s]
    
  }
  p.insample$calib = p_calib_temp
  rm(p_calib_temp)
  
  
  ## out-sample ====
  doParallel::registerDoParallel(cores = 6)
  outsample = foreach::foreach(s = 1:nSites, .combine = "mysimcombine_outsample", .multicombine = T) %dopar% {
    
    
    ### misclassification ====
    #### homogeneous ====
    fitout_m0_s = hetmis(errormat_by_site = obserrormat_by_site[-s],
                         model.choice = 'hom', cause.type = 'single',
                         pss.method = 'learn', shapeBase = .5, 
                         shape_pull = c(.5,.5),
                         sensprior_shape = c(1,1),
                         Dirprior_shape = rep(1, nCause),
                         nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                         sourcecode.path = sourcecode.path,
                         seed = 1, verbose = T, saveoutput = F)
    
    #### partly heterogeneous ====
    fitout_m1_s = hetmis(errormat_by_site = obserrormat_by_site[-s],
                         model.choice = 'het_part', cause.type = 'single',
                         pss.method = 'learn', shapeBase = .5, 
                         shape_pull = c(.5,.5), shape_hetsens = c(.5,.5),
                         sensprior_shape = c(1,1),
                         Dirprior_shape = rep(1, nCause),
                         nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                         sourcecode.path = sourcecode.path,
                         seed = 1, verbose = T, saveoutput = F)
    
    #### fully heterogeneous ====
    fitout_m2_s = hetmis(errormat_by_site = obserrormat_by_site[-s],
                         model.choice = 'het', cause.type = 'single',
                         pss.method = 'learn', shapeBase = .5, 
                         shape_pull = c(.5,.5), shape_hetsens = c(.5,.5), shape_hetrfp = c(.5,.5),
                         sensprior_shape = c(1,1),
                         Dirprior_shape = rep(1, nCause),
                         nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                         sourcecode.path = sourcecode.path,
                         seed = 1, verbose = T, saveoutput = F)
    
    
    ### calibration ====
    #### using homogeneous misclassification ====
    calibout_m0_s = calibratedva_v2(va_unlabeled = obsva[s,],
                                    cause.type = 'single', model.choice = "seqpshrink_up",
                                    Mmat.asDirich = fitout_m0_s$Mmat.asDirich$other,
                                    nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                                    sourcecode.path = sourcecode.path,
                                    seed = 1, verbose = T, saveoutput = F)
    
    #### using heterogeneous misclassification ====
    calibout_m1_s = calibratedva_v2(va_unlabeled = obsva[s,],
                                    cause.type = 'single', model.choice = "seqpshrink_up",
                                    Mmat.asDirich = fitout_m1_s$Mmat.asDirich$other,
                                    nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                                    sourcecode.path = sourcecode.path,
                                    seed = 1, verbose = T, saveoutput = F)
    
    #### using heterogeneous misclassification ====
    calibout_m2_s = calibratedva_v2(va_unlabeled = obsva[s,],
                                    cause.type = 'single', model.choice = "seqpshrink_up",
                                    Mmat.asDirich = fitout_m2_s$Mmat.asDirich$other,
                                    nMCMC = nMCMC, nBurn = nBurn, nThin = nThin, adapt_delta_stan = .9,
                                    sourcecode.path = sourcecode.path,
                                    seed = 1, verbose = T, saveoutput = F)
    
    
    # output array
    Mmat_s = array(dim = c(nMCMC, nCause, nCause, 3),
                   dimnames = list(NULL, allcauses, allcauses,
                                   c(fitout_m0_s$input$model.choice,
                                     fitout_m1_s$input$model.choice,
                                     fitout_m2_s$input$model.choice)))
    Mmat_s[,,,fitout_m0_s$input$model.choice] = fitout_m0_s$Mmat$other
    Mmat_s[,,,fitout_m1_s$input$model.choice] = fitout_m1_s$Mmat$other
    Mmat_s[,,,fitout_m2_s$input$model.choice] = fitout_m2_s$Mmat$other
    
    p_calib_s = array(dim = c(nMCMC, nCause, 3),
                      dimnames = list(NULL, allcauses, 
                                      c(fitout_m0_s$input$model.choice,
                                        fitout_m1_s$input$model.choice,
                                        fitout_m2_s$input$model.choice)))
    p_calib_s[,,fitout_m0_s$input$model.choice] = calibout_m0_s$MCMCout$p_calib
    p_calib_s[,,fitout_m1_s$input$model.choice] = calibout_m1_s$MCMCout$p_calib
    p_calib_s[,,fitout_m2_s$input$model.choice] = calibout_m2_s$MCMCout$p_calib
    
    print(s)
    
    list('Mmat' = Mmat_s, 'p_calib' = p_calib_s)
    
  }
  
  # rearranging for convenience
  Mmat_temp = array(dim = c(nMCMC, nSites, nCause, nCause, 3),
                    dimnames = list(NULL, allsites, allcauses, allcauses,
                                    c(fitout_m0$input$model.choice,
                                      fitout_m1$input$model.choice,
                                      fitout_m2$input$model.choice)))
  p_calib_temp = array(dim = c(nMCMC, nSites, nCause, 3),
                       dimnames = list(NULL, allsites, allcauses, 
                                       c(fitout_m0$input$model.choice,
                                         fitout_m1$input$model.choice,
                                         fitout_m2$input$model.choice)))
  for(s in 1:nSites){
    
    Mmat_temp[,s,,,] = outsample$Mmat[,,,,s]
    p_calib_temp[,s,,] = outsample$p_calib[,,,s]
    
  }
  outsample$Mmat = Mmat_temp
  outsample$p_calib = p_calib_temp
  rm(Mmat_temp, p_calib_temp)
  
  saveRDS(list('Mmat.obs' = Mmat.obs,
               'Mmat.insample' = Mmat.insample,
               'ic' = ic.df,
               'csmf.insample_uncalib' = p.insample$uncalib,
               'csmf.insample_calib' = p.insample$calib,
               'Mmat.outsample' = outsample$Mmat,
               'csmf.outsample_calib' = outsample$p_calib), file.path(output_dir, paste0('truem1_r', r)))
  
  
  print(paste0("Replication ", r))
  
}


saveRDS(list('nCause' = nCause, 'allcauses' = allcauses,
             'nSites' = nSites, 'allsites' = allsites,
             'ai.true' = ai.true, 'alpha.true' = alpha.true,
             'csmf.true' = p.true,
             'nMits' = nMits, 'nVA' = nVA, 'nReplicate' = nReplicate,
             'nBurn' = nBurn, 'nMCMC' = nMCMC, 'nThin' = nThin,
             'pss.true' = pss.true,
             'Mmat.true' = musij.true, 'csmf_uncalib.true' = p.uncalib.true),
        file.path(output_dir, 'truem1_datagen'))

