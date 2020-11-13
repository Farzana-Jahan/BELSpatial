#' Estimating equation for GEL estimates of beta
#'
#' @param tet initial value of regression coefficients beta
#' @param x design matrix
#' @param var residual variance from the fitted model and data
#' @param y observed response variable

g1<-function(tet,x,var,y)
{
  n<-dim(x)[1]
  beta1=(1/n)* sum(y-x%*%tet)
  beta2=(1/n)* (sum((y-x%*%tet)^2/var) -1)
  beta = cbind(beta1,beta2)
  return(beta)
}

#' gradient function for estimating equations of beta
#'
#' @param tet initial value of regression coefficients beta
#' @param x design matrix
#' @param var residual variance from the fitted model and data
#' @param y observed response variable

dg<-function(tet,x,var,y)
{
  n<-dim(x)[1]
  xx<-t(x)%*%x
  yx<-t(x)%*%y
  G<-matrix(c((1/n)*sum(-x),
              (2/n)*sum((t(x)%*%((x%*%tet)/var)-t(x)%*%(y/var)))),nrow=2,ncol=1)
  return(G)
}

# See in this function how I have (1) made it exportable, so it can be used
# outside of the package, unlike the last ones, (2) it imports the gel function
# from the gmm package via the "importFrom" tag, and (3) I've used some inline functions
# in the code itself to manage the fact that you need y and var in g1 and dg but
# the inputs for gel are only allowed to have tet and x.

#' generalised Empirical likelihood estimation for of beta
#'
#' @param x design matrix
#' @param y observed response
#' @param tet Initial value of regression coefficients
#' @param var residual variance
#' @importFrom gmm gel
#' @export
mele<- function(x,tet,y,var)
{
  gel_fit <- gel(g = function(tet,x){g1(tet,x,var,y)},
      x = x,
      tet0 = tet,
      gradv = function(tet,x){dg(tet,x,var,y)})
  return(unname(gel_fit$coefficients))
}


