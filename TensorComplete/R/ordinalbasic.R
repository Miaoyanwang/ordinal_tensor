epsilon=10^-4

logistic = function(x){
  return(1/(1+exp(-x)))
}


k_unfold=function(tnsr,m=NULL){
  if(is.null(m)) stop("mode m must be specified")
  num_modes <- tnsr@num_modes
  rs <- m
  cs <- (1:num_modes)[-m]
  unfold(tnsr,row_idx=rs,col_idx=cs)
}

#' The probability matrix given latent parameters.
#'
#' Compute the probability matrix given latent parameters from the cumulative model.
#' @usage theta_to_p(theta,omega)
#' @param theta a continuous-valued tensor (latent parameters).
#' @param omega the cut-off points.
#' @return a probability matrix in which the number of columns is possible outcomes and each row vector is corresponding probabilities at an entry of the tensor.
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices),-2,2),dim = indices)
#' b <- c(-1.5,0,1.5)
#' probability <- theta_to_p(arr,b);probability
#' @export
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


#' An ordinal tensor randomly simulated from the cumulative model.
#'
#' Simulate an ordinal tensor from the cumulative logistic model with the parameter tensor and the cut-off points.
#' @usage realization(theta,omega)
#' @param theta a continuous-valued tensor (latent parameters).
#' @param omega the cut-off points.
#' @return an ordinal tensor randomly simulated from the cumulative logistic model.
#' @references C. Lee and M. Wang. Tensor denoising and completion based on ordinal observations. \emph{International Conference on Machine Learning (ICML)}, 2020.
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices)),dim = indices)
#' b <- qnorm((1:3)/4)
#' r_sample <- realization(arr,b);r_sample
#' @export
#' @import tensorregress
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


#' Estimation of tensor entries from the cumulative model.
#'
#' Estimate the ordinal-valued tensor entries given latent parameters and a type of estimations.
#' @usage estimation(theta,omega,type = c("mode","mean","median"))
#' @param theta a continuous-valued tensor (latent parameters).
#' @param omega the cut-off points.
#' @param type type of estimations:
#'
#' \code{"mode"} specifies argmax based label estimation.
#'
#' \code{"mean"} specifies mean based label estimation.
#'
#' \code{"median"} specifies median based label estimation.
#' @return an estimated ordinal tensor given latent parameters and a type of estimations.
#' @references C. Lee and M. Wang. Tensor denoising and completion based on ordinal observations. \emph{International Conference on Machine Learning (ICML)}, 2020.
#' @examples
#' indices <- c(10,20,30)
#' arr <- array(runif(prod(indices),-2,2),dim = indices)
#' b <- c(-1.5,0,1.5)
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

#' Log-likelihood function (cost function).
#'
#' Return log-likelihood function (cost function) value evaluated at a given parameter tensor, an observed tensor, and cut-off points.
#' @usage likelihood(ttnsr,theta,omega,type = c("ordinal","Gaussian"))
#' @param theta a continuous-valued tensor (latent parameters).
#' @param omega the cut-off points.
#' @param ttnsr an observed tensor data.
#' @param type types of log-likelihood function.
#'
#' \code{"ordinal"} specifies log-likelihood function based on the cumulative logistic model.
#'
#' \code{"Gaussian"} specifies log-likelihood function based on the Gaussian model.
#' @return log-likelihood value at given inputs.
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

#' Bayesian Information Criterion (BIC) value.
#'
#' Compute Bayesian Information Criterion (BIC) given a parameter tensor, an observed tensor, the dimension, and the rank.
#' @usage bic(ttnsr,theta,omega,d,r)
#' @param ttnsr an observed tensor.
#' @param theta a continuous-valued tensor (latent parameters).
#' @param omega the cut-off points.
#' @param d dimension of the tensor.
#' @param r rank of the tensor.
#' @return BIC value at given inputs.
#' @export
bic = function(ttnsr,theta,omega=0,d,r){
  return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
}
