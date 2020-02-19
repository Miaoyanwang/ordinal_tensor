library(MASS)
library(rTensor)
library(pracma)
epsilon=10^-4

################### some functions used for analysis #######################################
#### Output: randomly sampled ordinal tensor based on logistic models given input parameters
#### Input : 
#### theta: continuous parameter tensor (latent parameter)
#### omega: the cut-off points (k-levels)
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

#### Output: predicted ordinal tensor from the given input parameters and types of predction
#### Input : 
#### theta: continuous parameter tensor (latent parameter)
#### omega: the cut-off points (k-levels)
#### type:`mode' the mostly likely label,`mean' the mean label, `median' the median label
estimation = function(theta,omega,type){ 
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

#### Output: the probability tensor converted from the parameter tensor based on logistic cumulative model
#### Input : 
#### theta: continuous parameter tensor (latent parameter)
#### omega: the cut-off points (k-levels)
theta_to_p=function(theta,omega){
  k = length(omega)
  p = matrix(nrow = length(theta),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(theta + omega[i]))
  }
  p =  cbind(p,rep(1,length(theta)))-cbind(rep(0,length(theta)),p)
  p[p<epsilon]=epsilon ## regularize
  return(p)
}

#### logistic function
logistic = function(x){
  return(1/(1+exp(-x)))
}

#### Output: log-likelihood function (type = "ordinal"), squared error ( type = "Gaussian")
#### Input : 
#### theta: continuous parameter tensor (latent parameter)
#### omega: the cut-off points (k-levels)
#### ttnsr: observed tensor

likelihood = function(ttnsr,theta,omega=0,type="ordinal"){
  index=which(is.na(ttnsr)==F & is.na(theta)==F)
  ttnsr=ttnsr[index]
  theta=theta[index]
  
  if(type=="Gaussian") return(sqrt(sum((ttnsr-theta)^2)))
  k = length(omega)
  p=theta_to_p(theta,omega)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
  
}



#### Output: BIC value given observed tensor, estimated parameter and rank
#### Input :
#### ttnsr: Observed tensor data
#### d: Dimension of tensor
#### r: Rank of tensor
bic = function(ttnsr,theta,omega=0,d,r){
  return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
}


################### sub functions for the main function `fit_ordinal' ########################
#### cost functions: `h1',`hc'
#### gradient functions: `g1',`gradient_tensor',`gc'
#### update function for factor matrices: `comb'
#### update function for a core tensor: `corecomb'
#### type == "Gaussian" will be used as subfunctions for `fit_continuous'

#### cost function based on unfolded matrix 
h1 = function(A_1,W1,ttnsr,omega,type="ordinal"){
  k = length(omega)
  theta_output =W1%*%c(A_1)
  if(type=="Gaussian"){ 
    l=sqrt(sum((theta_output[is.na(ttnsr)==F]-ttnsr[is.na(ttnsr)==F])^2))
    return(l)  
  }
  p = theta_to_p(theta_output,omega)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
}

#### cost function based on Tucker representation 
hc = function(A_1,A_2,A_3,C,ttnsr,omega,type="ordinal"){
  k = length(omega)
  theta_output = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
  if(type=="Gaussian"){
    l=sqrt(sum((theta_output[is.na(ttnsr)==F]-ttnsr[is.na(ttnsr)==F])^2))
    return(l)
  }
  p = theta_to_p(theta_output,omega)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
}

#### gradient function based on unfolded matrix
g1 = function(A_1,W,ttnsr,omega,type="ordinal"){
  k = length(omega)
  theta_output =W%*%c(A_1)
  if(type=="Gaussian"){
    pretheta=theta_output-as.matrix(ttnsr)
    pretheta[is.na(pretheta)==T]=0
    l=2*t(pretheta)%*%W
    return(l)
  }
  p = matrix(nrow = length(theta_output),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(theta_output + omega[i]))
  }
  q = matrix(nrow = length(theta_output),ncol = k+1)
  q[,1] <- p[,1]-1
  if(k>=2){
    for (i in 2:k) {
      q[,i] <- p[,i]+p[,i-1]-1
    }
  }
  q[,k+1] <- p[,k]
  l <- Reduce("+",lapply(1:(k+1),function(i) apply(rbind(W[which(c(ttnsr)==i),,drop=F])*q[which(c(ttnsr)==i),i],2,sum)))
  return(l)
}

#### gradient function with respect to theta
gradient_tensor=function(A_1,A_2,A_3,C,ttnsr,omega,type="ordinal"){
  k = length(omega)
  theta_output = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
  if(type=="Gaussian"){
    pretheta=theta_output-ttnsr
    pretheta[is.na(pretheta)==T]=0
    l=2*pretheta
    return(l)
  }
  p = matrix(nrow = length(theta_output),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(theta_output + omega[i]))
  }
  q = matrix(nrow = length(theta_output),ncol = k+1)
  q[,1] <- p[,1]-1
  if(k>=2){
    for (i in 2:k) {
      q[,i] <- p[,i]+p[,i-1]-1
    }
  }
  q[,k+1] <- p[,k]
  output=ttnsr
  for(i in 1:(k+1)){
    output[which(c(ttnsr)==i)]=q[which(c(ttnsr)==i),i]
  }
  return(output)
}

#### gradient function with respect to the core tensor (gradient not using kronecker product at all)
gc = function(A_1,A_2,A_3,C,ttnsr,omega,type="ordinal"){
  g = gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega,type) ## first, take entrywise gradient w.r.t. theta
  
  g = ttl(as.tensor(g),list(t(A_1),t(A_2),t(A_3)),ms = 1:3)@data ## then, take entrywise gradient w.r.t. core tensor
  
  return(g)
}

#### update one factor matrix at a time while holding others fixed 
comb = function(A,W,ttnsr,k,omega,alpha=TRUE,type="ordinal"){
  nA = A
  tnsr1 <- k_unfold(as.tensor(ttnsr),k)@data
  if (alpha==TRUE) {
    l <- lapply(1:nrow(A),function(i){optim(A[i,],
                                            function(x) h1(x,W,tnsr1[i,],omega,type),
                                            function(x) g1(x,W,tnsr1[i,],omega,type),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }else{
    l <- lapply(1:nrow(A),function(i){constrOptim(A[i,],
                                                  function(x) h1(x,W,tnsr1[i,],omega,type),
                                                  function(x) g1(x,W,tnsr1[i,],omega,type),
                                                  ui = as.matrix(rbind(W,-W)),ci = rep(-alpha,2*nrow(W)),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }
  return(nA)
}

#### update core tensor holding other factors fixed
corecomb = function(A_1,A_2,A_3,C,ttnsr,omega,alpha=TRUE,type="ordinal"){
  h <- function(x) hc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega,type)
  g <- function(x) c(gc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega,type)) 
  d <- optim(c(C@data),h,g,method="BFGS") 
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  
  return(C)
}



######################### the main function `fit_ordinal' #############################################
#### Output:
#### C,A_1,A_2,A_3: Estimated core tensor and factor matrices from the input observed tensor.
#### theta: Estimated latent parameter tensor (theta = ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
#### iteration: The number of iterations 
#### cost: Evaluations of cost functions at each iterations
#### omega: Estimated cut-off points vector (k-level)

#### Input:
#### ttnsr: Observed tensor data (k-level)
#### C,A_1,A_2,A_3: Initial points for a core tensor and factor matrices
#### omega: Cut-off points if it is known, (by default, omega = TRUE)
#### alpha: Signal level (by default, alpha = TRUE)

fit_ordinal = function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alpha = TRUE){
  
  if(is.logical(alpha)) alpha_minus=alpha_minus2=TRUE
  else{
    alpha_minus=alpha-epsilon
    alpha_minus2=alpha-2*epsilon
    prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    if(max(abs(prevtheta))>alpha){
      print("Warning: the input tensor exceeds the magnitude bound. Perform rescaling on the core tensor...")
      C=C/max(abs(prevtheta))*alpha/10
    }
  }
  
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omg = omega
  k=length(unique(as.factor(c(ttnsr))))-is.element(NA,ttnsr) ## for NA not being included in length
  d=dim(ttnsr)
  ttnsr=array(as.numeric(as.factor(ttnsr)),dim=d) ## code labels from 1 to K 
  
  while ((error > 10^-4)&(iter<10) ) {
    iter = iter +1
    #update omega
    prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    
    currentomega=omega
    if(is.logical(omg)) {
      omega=tryCatch(polr(as.factor(c(ttnsr))~offset(-c(prevtheta)))$zeta,error=function(c)"omega cannot be reliably estimated",warning=function(c)"omega cannot be reliably estimated")
      if(inherits(omega,"numeric")==FALSE){
        print("Warning: omega cannot be estimated! Omega from previous step is used");
        if(iter==1) omega=logit(1:(k-1)/k)
        else omega=currentomega
      }
      currentomega=omega
    }
    
    
    prev <- likelihood(ttnsr,prevtheta,omega)
    
    
    #update A_1
    W <-kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
    A_1<- comb(A_1,W,ttnsr,1,omega,alpha_minus2)
    #orthognalize A_1
    qr_res=qr(A_1)
    A_1=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),1)
    
    # update A_2
    W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
    A_2 <- comb(A_2,W,ttnsr,2,omega,alpha_minus)
    #orthognalize A_2
    qr_res=qr(A_2)
    A_2=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),2)
    
    # update A_3
    W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
    A_3 <- comb(A_3,W,ttnsr,3,omega,alpha)
    #orthognalize A_3
    qr_res=qr(A_3)
    A_3=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),3)
    
    # update C
    C_update <- corecomb(A_1,A_2,A_3,C,ttnsr,omega)
    
    if(!is.logical(alpha) & max(abs(ttl(C_update,list(A_1,A_2,A_3),ms=1:3)@data))>=alpha_minus2){
      theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      new <- likelihood(ttnsr,theta,omega)
      cost = c(cost,new); break 
    }else{
      C=C_update
      theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      new <- likelihood(ttnsr,theta,omega)
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
    
    
    message(paste(iter,"-th  iteration -- cost value is",new," -----------------"))
  }
  
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$theta= theta
  result$iteration <- iter
  result$cost = cost
  result$omega=omega
  return(result)
}

######################## functions for alternative methods ##################################################
#### Continuous tensor decomposition: `fit_continuous'
#### Ordinal matrix completion: `fit_ordinal_matrix'
#### One-bit tensor completion: `M_to_one', `one_to_M' (we used parameter estimation code from the reference in the main paper)

#### Continous Tucker decomposition with possibly missing entries
fit_continuous=function(ttnsr,C,A_1,A_2,A_3,alpha = TRUE){ 
  if(is.logical(alpha)) alpha_minus=alpha_minus2=TRUE
  else{
    alpha_minus=alpha-epsilon
    alpha_minus2=alpha-2*epsilon
    prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    if(max(abs(prevtheta))>alpha){
      print("Warning: the input tensor exceeds the magnitude bound. Perform rescaling on the core tensor...")
      C=C/max(abs(prevtheta))*alpha/10
    }
  }
  
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omega=0
  
  while ((error > 10^-4)&(iter<10) ) {
    iter = iter +1
    
    prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    prev <- likelihood(ttnsr,prevtheta,omega,type="Gaussian") 
    
    
    
    W <-kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
    A_1 <- comb(A_1,W,ttnsr,1,omega,alpha= alpha_minus2,type="Gaussian")
    #orthognalize A_1
    qr_res=qr(A_1)
    A_1=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),1)
    
    # update A_2
    W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
    A_2<- comb(A_2,W,ttnsr,2,omega,alpha= alpha_minus,type="Gaussian")
    #orthognalize A_2
    qr_res=qr(A_2)
    A_2=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),2)
    
    # update A_3
    W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
    A_3 <- comb(A_3,W,ttnsr,3,omega,alpha= alpha,type="Gaussian")
    #orthognalize A_3
    qr_res=qr(A_3)
    A_3=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),3)
    
    C_update <- corecomb(A_1,A_2,A_3,C,ttnsr,omega,type="Gaussian")
    
    if(!is.logical(alpha) & max(abs(ttl(C_update,list(A_1,A_2,A_3),ms=1:3)@data))>=alpha_minus2){
      theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      new <- likelihood(ttnsr,theta,omega,type="Gaussian")
      cost = c(cost,new); break 
    }else{
      C=C_update
      theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      new <- likelihood(ttnsr,theta,omega,type="Gaussian")
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
    
    message(paste(iter,"-th  iteration -- cost value is",new," -----------------"))
  }
  
  
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$theta = theta
  result$iteration <- iter
  result$cost = cost
  return(result)
}


#### special version of the function `fit_ordinal()'; ordinal matrix decomposition
fit_ordinal_matrix = function(ttnsr,A_1,A_2,omega=TRUE,alpha = TRUE){
  
  if(is.logical(alpha)) alpha_minus=TRUE
  else alpha_minus=alpha-epsilon
  
  if(!is.logical(alpha) & max(A_1%*%t(A_2))>= alpha-2*epsilon){
    A_2=A_2*(alpha-2*epsilon)/max(A_1%*%t(A_2))
  }
  
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omg = omega
  k=length(unique(as.factor(c(ttnsr))))-is.element(NA,ttnsr) 
  d=dim(ttnsr)
  ttnsr=array(as.numeric(as.factor(ttnsr)),dim=d) ## code labels from 1 to K 
  
  while ((error > 10^-4)&(iter<10) ) {
    iter = iter +1
    #update omega
    prevtheta <- A_1%*%t(A_2)
    
    currentomega=omega
    if(is.logical(omg)) {
      omega=tryCatch(polr(as.factor(c(ttnsr))~offset(-c(prevtheta)))$zeta,error=function(c)"omega cannot be reliably estimated",warning=function(c)"omega cannot be reliably estimated")
      if(inherits(omega,"numeric")==FALSE){
        print("Warning: omega cannot be estimated! Omega from previous step is used");
        if(iter==1) omega=logit(1:(k-1)/k)
        else omega=currentomega
      }
      currentomega=omega
    }
    
    
    prev <- likelihood(ttnsr,prevtheta,omega)
    
    # update A_2
    A_2 <- comb(A_2,A_1,ttnsr,2,omega,alpha = alpha_minus)
    
    A_1_update <- comb(A_1,A_2,ttnsr,1,omega,alpha = alpha)
    
    
    
    if(!is.logical(alpha) & max(abs(A_1_update%*%t(A_2)))>= alpha_minus){
      #orthognalize A_1
      qr_res=qr(A_1)
      A_1=qr.Q(qr_res)
      A_2=A_2%*%t(qr.R(qr_res))
      theta <- A_1%*%t(A_2)
      new <- likelihood(ttnsr,theta,omega)
      cost = c(cost,new); break 
    }else{
      A_1=A_1_update
      #orthognalize A_1
      qr_res=qr(A_1)
      A_1=qr.Q(qr_res)
      A_2=A_2%*%t(qr.R(qr_res))
      
      theta <- A_1%*%t(A_2)
      new <- likelihood(ttnsr,theta,omega)
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
  }
  
  
  result$A_1=A_1; result$A_2=A_2; result$omega=omega;
  result$iteration <- iter
  result$cost = cost
  result$theta=theta
  return(result)
}



#### 1bit-category-T: convert an order-K multi-level tensor to an order-(K+1) binary tensor
#### 1bit-sign-T: convert to an order-K binary tensor by comparing to its average

M_to_one=function(tensor,K,type="category"){
  if(type=="sign"){
    ave=mean(tensor,na.rm=T) ## taking average
    scale=sqrt(sum((tensor-ave)^2,na.rm=T)) ## F norm
    binary=2*((tensor-ave>0)-0.5) ## 1 for positive and -1 for negative
    E=array(1,dim=dim(binary))
    E[is.na(binary)]=0
    binary[E==0]=0
    return(list("binary"=binary,"E"=E,"ave"=ave,"scale"=scale))
  }
  
  binary=array(0,dim=c(dim(tensor),K-1))
  if(type=="category"){
    for(k in 1:(K-1)){
      prepare=tensor   
      prepare[tensor!=k]=-1
      prepare[tensor==k]=1
      binary[,,,k]=prepare
    }}else if(type=="cumulative"){
      for(k in 1:(K-1)){
        prepare=tensor   
        prepare[tensor>k]=-1
        prepare[tensor<=k]=1
        binary[,,,k]=prepare
      }
    } 
  
  if(K==2) binary=binary[,,,1]
  E=array(1,dim=dim(binary))
  E[is.na(binary)]=0
  binary[E==0]=0 ## missingness should be encoded as zero
  return(list("binary"=binary,"E"=E))
}

#### convert output from 1-bit tensor decomposition to M-level probability tensor
one_to_M=function(prob,K,type="categorical"){
  output=array(0,dim=c(dim(prob)[1:3],K))
  if(type=="categorical"){
    output[,,,1:(K-1)]=prob
    output[,,,K]=1-apply(prob,c(1:3),sum)
    return(output)
  }else if(type=="cumulative"){
    output[,,,1]=prob[,,,1]
    for(k in 2:(K-1)){
      output[,,,k]=prob[,,,k]-prob[,,,(k-1)]
    }
    output[,,,K]=1-prob[,,,(K-1)]
    return(output)
  }
}
