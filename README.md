# Modeling Structure and Cross-Country Heterogeneity in VA Misclassification

## Description

Codes to implement the country-specific VA misclassification modeling proposed in [Pramanik et al. (2024+)](https://arxiv.org/abs/2312.03192). 

## Credit

Sandipan Pramanik developed this repository and wrote the codes included here.

## Required Package Installation

Before sourcing the functions, we recommend installing latest versions of required packages `reshape2`, `tidyverse`, `doParallel`, `rstan`, `LaplacesDemon`, and `ggplot2`.

## Sourcing Functions

Specify the path to the `sourcecode` folder available from this `GitHub` repository. The functions can then be sourced as

``` r
sourcecode.path = ...    # specifies path ".../sourcecode" to the sourcecode folder
source(file.path(sourcecode.path, "hetmis-functions.R"))    # sources ".../sourcecode/hetmis-functions.R"
```

## Key Function

In the `sourcecode` folder, `hetmis()` function within `hetmis-functions.R` is the main function
 
 * Arguments:
   * `va_labeled`: Data frame. COD estimates from an algorithm (e.g. verbal autopsy algorithms such as EAVA, InsilicoVA, InterVA) 
     * Rows are individuals. Recommended to assign row names.
     * For $C$ causes there should be $C$ named columns for each cause and one column named `site` indicating which site each row (individual) belongs to.
     * `va_labeled[l,"cause1"]` is 1 if `"cause1"` is the cause of death for individual `l` and 0 otherwise. The entries under the cause columns `va_labeled[l,c("cause1", ... ,"causeC")]` must sum to 1.
     * `va_labeled[l,"site"]` or `va_labeled$site[l]` is a character and indicates the site/domain/country individual `l` is from.
     * Only used if `gold_standard` is also provided.
   * `gold_standard`: Data frame. COD estimates from a gold standard diagnosis.
     * Similarly structured as `va_labeled`.
     * Only used if `va_labeled` is also provided.
   * `errormat_by_site`: List of observed misclassification count matrices from each site/domain/country
     * The components must have names and should indicate the domains
     * `errormat_by_site[[l]]` is the observed misclassification count matrix in site/domain/country `l`. For $C$ causes it's a $C \times C$ matrix and it must have causes as row and column names.
   * `model.choice`: Character. Indicates which model to fit. Available options
     * `"hom"` for general homogeneous model
     * `"het_part"` for partly-heterogeneous model
     * `"het"` for fully-heterogeneous model
     * `"pooled_empirical"` for pooled empirical estimates
     * `"separate_empirical"` for country-specific empirical estimates
   * `pss.method`: Character. Available options
     * `"fixed"` for fixing effect size values
     * `"learn"` for learning effect sizes from data
   * `pss_fixed`: List of fixed effect sizes for `pss.method="fixed"`. Provide as `pss_fixed = list("pss_pull" = 2, "pss_hetsens" = 2, "pss_hetrfp" = 2)`.
   * `shape_pull`, `shape_hetsens`, and `shape_hetrfp` are shape parameters in Beta prior on transformed effect sizes. Default `c(.5, .5)`.
   * `sensprior_shape`: shape parameters in Beta prior on intrinsic accuracies. Default `c(1, 1)` indicating the Uniform prior.
   * `Dirprior_shape`: shape parameters in Dirichlet prior on pull. Default `rep(1, C)` indicating the Uniform prior on $C$-dimensional simplex.
   * `nMCMC`: Posterior sample size. Default 5000.
   * `nBurn`: Burn-in. Default 2000.
   * `nThin`: Thinning size. Default 1.
   * `adapt_delta_stan`: Same as `adapt_delta` in `rstan::sampling()`. Default 0.9.
   * `sourcecode.path`: File path (as ".../sourcecode") to the `sourcecode` folder
   * `seed`: Random seed
   * `verbose`: Logical. `T` for printing progress report (Default). `F` otherwise.
   * `saveoutput`: Logical. `T` for saving output (Default). `F` otherwise.
   * `output_dir`: Character. Name of the folder to save output when `saveoutput=T`. Default `"hetmis_output"`.
   * `output_filename`: Character. Name of the output when `saveoutput=T`. Default `paste0(model.choice, "_", cause.type)`.
 
 * Output: A list with the following components:
   * `MCMCout`: MCMC output from the stanfit. This contains posterior samples of parameters, and some MCMC diagnostics (`"max_Rhat"`, `"min_ess_bulk"`, `"num_divergent"`, and `"num_max_treedepth"`).
   * `Mmat`: List of posterior samples of country-specific misclassification matrices.
   * `Mmat.postmean`: List of posterior mean of country-specific misclassification matrices.
   * `Mmat.asDirich`: List. For each country and gold-standard cause (such as MITS), this provides a Dirichlet approximation of posterior distributions of misclassification rates. Has the same structure as `Mmat.postmean`. For country `l`, `Mmat.asDirich[[l]]` is a $C \times C$ matrix where `Mmat.asDirich[[l]][i,]` are Dirichlet shape parameters that best approximate the posterior distribution of algorithm in country `l` given gold standard cause `i`.

## Reproducing Simulation Results

* `mwe.R` within the `MWE` folder contains minimal working examples
* `R` scripts in the `simulation` folder reproduce simulation results. Scripts in `setting1`, `setting2`, and `setting3` correspond to the three simulation settings.
