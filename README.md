**Functions for Markov Chain Monte Carlo (MCMC) posterior inference on spatial Bayesian Empirical Likelihood (SBEL) models**
 
This repository contains the code to draw posterior samples using MCMC procedure via a Metropolis Hastings (MH) algorithm from SBEL_CAR models proposed by Jahan et al. (2020). 
 
Using the function, BEL_leroux_new, posterior samples for fixed regression effects, precision paramter and spatial random effects can be drawn for SBEL-Leroux (specifying appropriate value of spatial autocorrelation Rho), SBEL-BYM (specifying Rho =1). The function also can be used to draw posterior samples from SBEL-IG model (specifying Rho=0) proposed by Chaudhuri et al. (2011) applying Independent Gaussian priors for the spatial random effects. 

The function, BSHEL, can be used to draw posteriors of interest from the Bayesian semi-paramteric Hierarchical Empirical Likelihood Models proposed by Porter et al. (2015). 

**References:**
1. Jahan F, Kennedy DW, Duncan EW & Mengesen KL(2020). Evaluation of spatial Bayesian Empirical Likelihood models in analysis of small area. Working paper, submitted to Statistical Modelling.
2. Chaudhuri, S., & Ghosh, M. (2011). Empirical likelihood for small area estimation. Biometrika, 473-480. https://www.jstor.org/stable/23076164?seq=1#metadata_info_tab_contents
3. Porter, A. T., Holan, S. H., & Wikle, C. K. (2015). Bayesian semiparametric hierarchical empirical likelihood spatial models. Journal of Statistical Planning and Inference, 165, 78-90. https://www.sciencedirect.com/science/article/pii/S0378375815000749
