#' BEL spatial function for BYM and Leroux
#' @param y observed response from data
#' @param x design matrix including a column of 1 for intercept and the other from dataset
#' @param n number of observations/small areas in the dataset
#' @param p number of regression paramters or number of columns in design matrix x
#' @param var residual variance
#' @param rho fixed value for spatial dependence rho=1 gives BYM rho=0 gives non spatial model
#' @param niter number of iterations
#' @param beta_init initial values for beta
#' @param psi_init initial values of spatial random effect psi
#' @param tau_init initial value of precision paramter tau
#' @param R neighbourhood matrix from spatial data
#' @param wi current EL weights
#' @param sd_psi value of sd for proposal distribution for psi
#' @param sd_beta value of sd for proposal distribution for beta
#' @param sd_tau value of sd for proposal distribution for tau
#' @importFrom emplik el.test
#' @importFrom MASS mvrnorm
#' @import stats
#' @export


BEL_leroux_new<-function(y,x,n,p,var,rho,niter,beta_init, psi_init, tau_init,R, wi, sd_psi, sd_beta, sd_tau)
{
  n.psi<- 0
  n.beta<- 0
  n.tau<- 0
  psi_sample<- matrix(nrow= n, ncol=niter)
  beta_sample <- matrix(nrow= p, ncol=niter)
  tau_sample <- c()
  psi_sample[,1]<- psi_init
  beta_sample[,1]<-beta_init
  tau_sample[1]<-tau_init
  beta<-beta_sample[,1]
  psi<-psi_sample[,1]
  tau<-tau_sample[1]
  # setting a loop for iteration
  for(i in 2: niter){

    proposal.mean.psi<- psi_sample[,i-1]
    proposal.sd.psi<-sd_psi*diag(n)
    psi_star<-mvrnorm(1,proposal.mean.psi,proposal.sd.psi)
    beta<- beta_sample[,i-1]
    mu_new<-x%*%beta+ psi_star
    wi_star<-el.test(y-mu_new, 0)$wts
    wi_star<-wi_star/sum(wi_star)

    mu_old<-x%*%beta + psi
    wi_orig<-el.test(y-mu_old, 0)$wts
    wi_orig<-wi_orig/sum(wi_orig)

    # checking the constraint 14
    if(all(wi_star>0) & (sum(wi_star)-1) < 0.0001)
    {
      # perform a MH step for psi*
      D_g_inv<-(1-rho)*diag(n)+rho*R # specify rho and neighbourhood matrix R beforehand
      pdr_psi<- target_psi(w=wi_star, psi=psi_star, D=D_g_inv, tau=tau)-
        target_psi(w=wi_orig, psi=psi, D=D_g_inv, tau=tau)#posterior density ratio for step 2, sampling psi
      if(rexp(1) > -pdr_psi){
        psi<-psi_star
        wi<- wi_star
        n.psi<-n.psi+1
      }

      else{
        psi<-psi
        wi<-wi
      }
    }

    psi_sample[,i]<-psi

    # Step 3 : sampling fixed effect beta
    proposal.mean.beta<-unname(lm(y~x-1, weights = wi)$coefficients)
    proposal.var.beta<- sd_beta*diag(p) # specifying arbitrarily, large variance but not huge
    beta_proposed<- mvrnorm(1,proposal.mean.beta,proposal.var.beta)

    wi_beta<- el.test(y-x%*%beta_proposed-psi, 0)$wts
    wi_beta<- wi_beta/sum(wi_beta)
    wi_orig_2<-el.test(y-x%*%beta-psi,0)$wts
    wi_orig_2<-wi_orig_2/sum(wi_orig_2)


    if(all(wi_beta)>0 & (sum(wi_beta)-1)< 0.0001)
    {
      pdr_beta<-target_beta(beta_proposed,w=wi_beta,proposal.mean.beta,g=10,tau)-
        target_beta(beta,w=wi_orig_2,proposal.mean.beta,g=10,tau)# posterior density ratio for beta
      if(rexp(1) > -pdr_beta){
        wi<- wi_beta
        beta<-beta_proposed
        n.beta<-n.beta+1
      }
      else{
        wi<-wi
        beta<-beta
        n.beta<-n.beta}
    }
    beta_sample[,i]<-beta

    # step 4 : sampling precision parameter tau
    tau_proposed<- exp(rnorm(1,mean= log(tau),sd=sd_tau)) # random walk
    pdr_tau<- target_tau(tau_proposed,beta_sample[,i],proposal.mean.beta,1,1,psi_sample[,i],D_g_inv)-
      target_tau(tau,beta_sample[,i-1],proposal.mean.beta,1,1,psi_sample[,i-1],D_g_inv) # posterior density ratio

    if(rexp(1)> -(pdr_tau+log(tau_proposed)-log(tau))) # if(rexp(1)> pdr_tau)
    {
      tau<-tau_proposed
      n.tau<-n.tau+1
    }

    tau_sample[i]<- tau
  }
  acceptance<- data.frame(parameter= c("beta", "psi", "tau"),
                          acceptance_rate= c((n.beta/niter)*100, (n.psi/niter)*100,
                                             (n.tau/niter)*100))
  # effective sample size
  # autocorrelation
  acf_psi<-c()
  for(j in 1:n)
  {
    acf<-acf(psi_sample[j,],plot=F)[[1]]
    n<-length(acf)
    acf_psi[j]<-sum(acf[2:n])
  }
  acf_beta0<-acf(beta_sample[1,],plot=F)[[1]]
  acf_beta1<-acf(beta_sample[2,],plot=F)[[1]]
  acf_tau<-acf(tau_sample,plot=F)[[1]]
  # ess
  ess_psi<-c()
  for(j in 1:n)
  {
    ess_psi[j]<- niter/ (1+2*acf_psi[j])
  }
  ess_beta0<- n.beta/ (1+2*sum(acf_beta0))
  ess_beta1<- n.beta/ (1+2*sum(acf_beta1))
  ess_tau<-n.tau/(1+2*sum(acf_tau))
  ess<- list(ess_psi,ess_beta0, ess_beta1, ess_tau)
  output<- list(Beta= beta_sample,psi= psi_sample,
                tau= tau_sample, acceptance_rate=acceptance,
                effective_sample_size=ess)
  output
}
#' BEL spatial function for Moran Basis Prior (POrter et al. 2015)
#' @param y observed response from data
#' @param x design matrix including a column of 1 for intercept and the other from dataset
#' @param n number of observations/small areas in the dataset
#' @param p number of regression paramters or number of columns in design matrix x
#' @param q number of positive eigen values of MBM
#' @param var residual variance
#' @param niter number of iterations
#' @param beta_init initial values for beta
#' @param psi_init initial values of spatial random effect psi
#' @param tau_init initial value of precision paramter tau
#' @param M n × q matrix with the columns being the eigenvectors corresponding to the q largest nonzero eigenvalues of the matrix PBP,P be the latent process space and B be the adjacency matrix
#' @param MBM M(B_plus-B)M, where B_plus is a diagonal matrix having number of neighbours for each area in its diagonal
#' @param wi current EL weights
#' @param sd_psi value of sd for proposal distribution for psi
#' @param sd_beta value of sd for proposal distribution for beta
#' @param sd_tau value of sd for proposal distribution for tau
#' @importFrom emplik el.test
#' @importFrom MASS mvrnorm
#' @import stats
#' @export
BSHEL<-function(y,x,n,p,q,var,niter,beta_init, psi_init, tau_init,M,MBM, wi, sd_psi, sd_beta, sd_tau)
{
  # setting initial counter
  n.psi<- 0
  n.beta<- 0
  n.tau<- 0
  psi_sample<- matrix(0,nrow= q, ncol=niter)
  beta_sample<- matrix(0,nrow= p, ncol=niter)
  tau_sample<- c()
  psi_sample[,1]<- psi_init
  beta_sample[,1]<-beta_init
  tau_sample[1]<-tau_init
  beta<-beta_sample[,1]
  psi<-psi_sample[,1]
  tau<-tau_sample[1]
  # setting a loop for iteration
  for(i in 2: niter){

    proposal.mean.psi<- psi_sample[,i-1]
    proposal.sd.psi<-sd_psi*diag(q)
    psi_star<-mvrnorm(1,proposal.mean.psi,proposal.sd.psi)
    beta<- beta_sample[,i-1]
    mu_new<-x%*%beta+ M%*%psi_star
    wi_star<-el.test(y-mu_new, 0)$wts
    wi_star<-wi_star/sum(wi_star)
    mu_old<-x%*%beta + M%*%psi
    wi_orig<-el.test(y-mu_old, 0)$wts
    wi_orig<-wi_orig/sum(wi_orig)

    # checking the constraint 14
    if(all(wi_star>0) & (sum(wi_star)-1) < 0.0001)
    {
      # perform a MH step for psi*
      pdr_psi<- as.numeric(target_yq(w=wi_star, a=psi_star, MBM, tau=tau)-
                             target_yq(w=wi_orig, a=psi, MBM, tau=tau))
      if(rexp(1) > - pdr_psi){
        psi<-psi_star
        wi<- wi_star
        n.psi<-n.psi+1
      }

      else{
        psi<-psi
        wi<-wi
      }
    }

    psi_sample[,i]<-psi
    # Step 3 : sampling fixed effect beta
    proposal.mean.beta<-unname(lm(y~x-1, weights = wi)$coefficients)
    proposal.var.beta<- sd_beta*diag(p) # specifying arbitrarily, large variance but not huge
    beta_proposed<- mvrnorm(1,proposal.mean.beta,proposal.var.beta)
    wi_beta<- el.test(y-x%*%beta_proposed-M%*%psi, 0)$wts
    wi_beta<- wi_beta/sum(wi_beta)
    wi_orig_2<-el.test(y-x%*%beta-M%*%psi,0)$wts
    wi_orig_2<-wi_orig_2/sum(wi_orig_2)

    if(all(wi_beta>0) & (sum(wi_beta)-1)< 0.0001)
    {
      pdr_beta<-target_beta(beta_proposed,w=wi_beta,proposal.mean.beta,g=10,tau)-
        target_beta(beta,w=wi_orig_2,proposal.mean.beta,g=10,tau)# posterior density ratio for beta
      if(rexp(1) > -pdr_beta){
        wi<- wi_beta
        beta<-beta_proposed
        n.beta<-n.beta+1
        var<- as.numeric(var(y- x%*%beta))
      }
      else{
        wi<-wi
        beta<-beta
        n.beta<-n.beta
        var<- var}
    }
    beta_sample[,i]<-beta

    # step 4 : sampling precision parameter tau

    tau_proposed<-exp(rnorm(1,mean= log(tau),sd=sd_tau)) # random walk
    pdr_tau<- target_tau_porter(tau_proposed,beta_sample[,i],proposal.mean.beta,1,0.01,psi_sample[,i],MBM,p=p,q=q)-
      target_tau_porter(tau,beta_sample[,i],proposal.mean.beta,1,0.01,psi_sample[,i],MBM,p=p,q=q) # posterior density ratio

    if(rexp(1)> -(pdr_tau+log(tau_proposed)-log(tau))) # if(rexp(1)> pdr_tau)
    {
      tau<-tau_proposed
      n.tau<-n.tau+1
    }

    tau_sample[i]<- tau
  }

  acceptance<- data.frame(parameter= c("beta", "psi", "tau"),
                          acceptance_rate= c((n.beta/niter)*100, (n.psi/niter)*100,
                                             (n.tau/niter)*100))
  n.psi<-n.psi
  # effective sample size
  # autocorrelation
  acf_psi<-c()
  for(j in 1:q)
  {
    acf<-acf(psi_sample[j,],plot=F)[[1]]
    n<-length(acf)
    acf_psi[j]<-sum(acf[2:q])
  }
  acf_beta0<-acf(beta_sample[1,],plot=F)[[1]]
  acf_beta1<-acf(beta_sample[2,],plot=F)[[1]]
  acf_tau<-acf(tau_sample,plot=F)[[1]]
  # ess
  ess_psi<-c()
  for(j in 1:q)
  {
    ess_psi[j]<- niter/ (1+2*acf_psi[j])
  }
  ess_beta0<- n.beta/ (1+2*sum(acf_beta0))
  ess_beta1<- n.beta/ (1+2*sum(acf_beta1))
  ess_tau<-n.tau/(1+2*sum(acf_tau))
  ess<- list(ess_psi, ess_beta0, ess_beta1, ess_tau)
  output<- list(Beta= beta_sample, psi= psi_sample, tau= tau_sample, acceptance_rate=acceptance,
                effective_sample_size=ess)
  output
}
