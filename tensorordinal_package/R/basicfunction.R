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
#' Randomly sample an ordinal tensor from the cumulative model
#'
#' Random sampling from the cumulative logistic model given theta and omega
#' @usage realization(theta,omega)
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points (k-levels)
#' @return ordinal tensor randomly sampled from the cumulative logistic model
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices)),dim = indices)
#' b <- runif(3)
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


#' Logistic function
#'
#' Logistic function
#' @usage logistic(x)
#' @param x numeric value
#' @return logistic function value evaluated at x
#' @export
logistic = function(x){
  return(1/(1+exp(-x)))
}

#' Predict the entries of a tensor from the cumulative model
#'
#' Predict the ordinal-valued entries of a tensor with given parameters and a type of predictions
#' @usage estimation(theta,omega,type = c("mode","mean","median"))
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points (k-levels)
#' @param type type of predictions:
#'
#' \code{"mode"} specifies argmax based label prediction
#'
#' \code{"mean"} specifies mean based label prediction
#'
#' \code{"median"} specifies median based label prediction
#' @return predicted ordinal tensor from the given parameters and types of prediction
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices)),dim = indices)
#' b <- runif(3)
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
#' Log-likelihood function (cost function)
#' @usage likelihood(ttnsr,theta,omega,type = c("ordinal","Gaussian"))
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points (k-levels)
#' @param ttnsr observed tensor
#' @param type types of log-likelihood function
#'
#' \code{"ordinal"} specifies log-likelihood function with the cumulative logistic model
#'
#' \code{"Gaussian"} specifies squared error function
#' @return a value evaluated on given type and points
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

#' Bayesian Information Criterion
#'
#' Bayesian Information Criterion evaluated at given parameters, dimension of tensor and rank of tensor
#' @usage bic(ttnsr,theta,omega,d,r)
#' @param ttnsr Observed tensor data
#' @param theta continuous parameter tensor (latent parameter)
#' @param omega the cut-off points (k-levels)
#' @param d Dimension of tensor
#' @param r Rank of tensor
#' @return BIC value given observed tensor, estimated parameter and rank
#' @export
bic = function(ttnsr,theta,omega=0,d,r){
  return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
}
