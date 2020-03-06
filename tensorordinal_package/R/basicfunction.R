epsilon=10^-4
theta_to_p=function(theta,omega){
  epsilon=10^-4
  k = length(omega)
  p = matrix(nrow = length(theta),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(theta + omega[i]))
  }
  p =  cbind(p,rep(1,length(theta)))-cbind(rep(0,length(theta)),p)
  p[p<epsilon]=epsilon ## regularize
  return(p)
}

logistic = function(x){
  return(1/(1+exp(-x)))
}


#' An ordinal tensor randomly sampled from the cumulative model
#'
#' Sample randomly from the cumulative logistic model with the parameter tensor and the cut-off points
#' @usage realization(theta,omega)
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points
#' @return an ordinal tensor randomly sampled from the cumulative logistic model
#' @references Lee and Wang (2020) <arXiv:2002.06524>.
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices)),dim = indices)
#' b <- qnorm((1:3)/4)
#' r_sample <- realization(arr,b);r_sample
#' @export
#' @import rTensor
realization = function(theta,omega){
  k=length(omega)
  theta=as.tensor(theta)
  theta_output <- c(theta@data)
  p = theta_to_p(theta_output,omega)
  for (j in 1:length(theta_output)) {
    theta_output[j] <-  sample(1:(k+1),1,prob = p[j,])
  }
  return(as.tensor(array(theta_output,dim =theta@modes)))
}


#' Estimation of tensor entries from the cumulative model
#'
#' Estimate the ordinal-valued tensor entries from latent parameters and a type of estimations
#' @usage estimation(theta,omega,type = c("mode","mean","median"))
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points
#' @param type type of estimations:
#'
#' \code{"mode"} specifies argmax based label estimation
#'
#' \code{"mean"} specifies mean based label estimation
#'
#' \code{"median"} specifies median based label estimation
#' @return an estimated ordinal tensor from the latent parameters and types of estimations
#' @references Lee and Wang (2020) <arXiv:2002.06524>.
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices)),dim = indices)
#' b <- qnorm((1:3)/4)
#' r_predict <- estimation(arr,b,type = "mode");r_predict
#' @export
estimation = function(theta,omega,type = c("mode","mean","median")){
  k=length(omega)
  if(is.matrix(theta)){
    theta_output=c(theta)
  }else{
    theta=as.tensor(theta)
    theta_output <- c(theta@data)
  }
  p = theta_to_p(theta_output,omega)
  if(type=="mode"){  # score prediction based on the mode
    for (j in 1:length(theta_output)) theta_output[j] <-  which.max(p[j,])
  }else if(type=="mean"){ # score prediction based on the mean
    for (j in 1:length(theta_output)) theta_output[j] <-  sum(p[j,]*(1:(k+1)))
  }else if(type=="median"){# score prediction based on the median
    for (j in 1:length(theta_output)) theta_output[j] <-  which(c(cumsum(p[j,]),1)>=0.5)[1] ## median
  }
  if(is.matrix(theta)) return(matrix(theta_output,nrow=dim(theta)[1],ncol=dim(theta)[2]))
  else return(as.tensor(array(theta_output,dim =theta@modes)))
}

#' Log-likelihood function (cost function)
#'
#' Return log-likelihood function (cost function) value evaluated at a given parameter tensor, an observed tensor and cut-off points
#' @usage likelihood(ttnsr,theta,omega,type = c("ordinal","Gaussian"))
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points
#' @param ttnsr an observed tensor data
#' @param type types of log-likelihood function
#'
#' \code{"ordinal"} specifies log-likelihood function with the cumulative logistic model
#'
#' \code{"Gaussian"} specifies log-likelihood function with the Gaussian model
#' @return log-likelihood value at given inputs
#' @export

likelihood = function(ttnsr,theta,omega=0,type = c("ordinal","Gaussian")){
  index=which(is.na(ttnsr)==F & is.na(theta)==F)
  ttnsr=ttnsr[index]
  theta=theta[index]

  if(type=="Gaussian") return(sqrt(sum((ttnsr-theta)^2)))
  k = length(omega)
  p=theta_to_p(theta,omega)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))

}

#' Bayesian Information Criterion (BIC) value
#'
#' Obtain Bayesian Information Criterion (BIC) evaluated at a given parameter tensor, an observed tensor, dimension and rank of the tensors
#' @usage bic(ttnsr,theta,omega,d,r)
#' @param ttnsr an observed tensor data
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points
#' @param d dimension of the tensor
#' @param r rank of the tensor
#' @return BIC value at given inputs
#' @export
bic = function(ttnsr,theta,omega=0,d,r){
  return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
}
