#' Conditional Posterior for Beta
#' @param b current value of beta
#' @param w current EL weights
#' @param b_mean mean of the proposal distribution
#' @param g value of zellnor's g prior taken as 10
#' @param tau current value of precision  parameter
#' @export
target_beta<-function(b,w,b_mean,g,tau){
  sum(log(w))-0.5*(t(b-b_mean)%*%(b-b_mean))*g*tau
}

#' Conditional Posterior for psi
#' @param w current EL weights
#' @param psi current value of psi
#' @param D variance covariance matrix of spatial random effects
#' @param tau current value of precision parameter
#' @export
target_psi<-function(w,psi,D,tau){
  sum(log(w)) -0.5*t(psi)%*%D%*%psi*tau
}

#' Supporting function 1 to compute Conditional Posterior for tau
#' @param b current value of beta
#' @param b_mean mean of the proposal distribution
#' @param g value of zellnor's g prior taken as 10
#' @param tau current value of precision parameter
#' @export
beta_giv_tau<-function(b,b_mean,tau,g=10){
  btb<-as.vector(t(b-b_mean)%*%(b-b_mean))
  btb_tau<-btb*g*tau
  (-0.5*btb_tau)
}
#' Supporting function 2 to compute Conditional Posterior for tau
#' @param psi current value of psi
#' @param D variance covariance matrix of spatial random effects
#' @param tau current value of precision parameter
#' @export
psi_giv_tau<-function(psi,tau,D){-0.5*((psi%*%D)%*%psi)*tau}

#' Supporting function 3 to compute Conditional Posterior for tau
#' @param alpha_1 scale parameter of inverse gamma distribution
#' @param alpha_2 shape parameter of inverse gamma distribution
#' @param tau current value of precision parameter
#' @export
tau_prior<-function(tau,alpha_1,alpha_2){-(1+alpha_1)*log(tau)-(alpha_2/tau)}

#' Conditional Posterior for tau
#' @param alpha_1 shape parameter for tau inverse gamma prior.
#' @param alpha_2 scale parameter for tau inverse gamma prior.
#' @param beta current value of beta
#' @param beta_mean mean of the proposal distribution
#' @param g value of zellnor's g prior taken as 10
#' @param tau current value of precision parameter
#' @param psi current value of psi
#' @param D variance covariance matrix of spatial random effects
#' @export
target_tau<-function(tau,beta,beta_mean,alpha_1,alpha_2,psi,D,g=10){
  psi_giv_tau(psi,tau,D)+ beta_giv_tau(beta,beta_mean,tau,g=10)+tau_prior(tau,alpha_1,alpha_2)
}

#' For Porter (2015) creation of matrix M
#' @param y observed data for response variable
#' @param x design matrix
#' @param B adjacency matrix
#' @export
M_create=function(y,x,B)
{
  n<-length(y)
  xx<-t(x)%*%x
  p=diag(1,n)-x%*%solve(xx)%*%t(x)
  pbp=t(p)%*%B%*%p
  q<-sum(eigen(pbp)$values>0)
  Re(eigen(pbp)$vectors[,1:q])
}

#' For Porter (2015) creation of matrix M
#' @param M n × q matrix with the columns being the eigenvectors corresponding to the q largest nonzero eigenvalues of the matrix PBP,P be the latent process space and B be the adjacency matrix
#' @param B adjacency matrix
#' @param Bplus diagional matrix stating the number of neghbours of each small area
#' @export
MBM_create=function(M,B, Bplus)
{
  t(M)%*%(Bplus-B)%*%M
}

#' For Porter (2015) conditional posterior of yq
#' @param a current/proposed value of yq, dimension reduced psi, spatial random effects
#' @param w EL weights
#' @param MBM MBM M(B_plus-B)M, where B_plus is a diagonal matrix having number of neighbours for each area in its diagonal, a positive definite matrix, variance covariance for spatial random effects
#' @param tau current value of precision parameter
#' @export
target_yq<-function(a,w,MBM,tau){sum(log(w))-0.5*t(a)%*%MBM%*%a*tau}

#' For Porter (2015) supporting function for conditional posterior of tau
#' @param a current/proposed value of yq, dimension reduced psi, spatial random effects
#' @param tau current value of precision parameter
#' @param MBM M(B_plus-B)M, where B_plus is a diagonal matrix having number of neighbours for each area in its diagonal, a positive definite matrix,variance covariance for spatial random effects
#' @export
yq_giv_tau<-function(a,tau,MBM){as.numeric(-0.5*(a%*%MBM%*%a)*tau)}

#' For Porter (2015) conditional posterior of tau
#' @param tau current value of precision parameter
#' @param alpha_1 shape parameter for tau inverse gamma prior.
#' @param alpha_2 scale parameter for tau inverse gamma prior.
#' @param beta current value of beta
#' @param beta_mean mean of the proposal distribution
#' @param g value of zellnor's g prior taken as 10
#' @param psi current value of psi
#' @param MBM MBM M(B_plus-B)M, where B_plus is a diagonal matrix having number of neighbours for each area in its diagonal, a positive definite matrix, variance covariance for spatial random effects
#' @param p number of fixed effect parameters
#' @param q number of positive eigenvalues of PBP (Porter 2015)
#' @export
target_tau_porter<-function(tau,beta,beta_mean,alpha_1,alpha_2,psi,MBM,g,p,q=q){
  ((p+q)/2)*log(tau)+ yq_giv_tau(psi,tau,MBM)+ beta_giv_tau(beta,beta_mean,tau,g=10)+
    tau_prior(tau,alpha_1,alpha_2)
}


