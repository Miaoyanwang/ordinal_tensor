library(MASS)
library(rTensor)
library(pracma)
epsilon=10^-4

########### simulate ordinal tensors based on logist models with arbitrary k
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

########### predict ordinal observations based on the estimated parameter
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
  for (j in 1:length(theta_output)) theta_output[j] <-  sum(p[j,]*(1:(k+1))) ## why rounding to integer? : I thought prediction has to be an integer because all observations are integers.
  }else if(type=="median"){# score prediction based on the median 
  for (j in 1:length(theta_output)) theta_output[j] <-  which(c(cumsum(p[j,]),1)>=0.5)[1] ## median
  }
  if(is.matrix(theta)) return(as.matrix(theta_output,dim=dim(theta)))
  else return(as.tensor(array(theta_output,dim =theta@modes)))
}   

########### cost function with arbitrary k ###########
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

########### cost function based on Tucker representation ###########
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

########### gradient with arbitrary k ###########
# gradient with respect to factor matrices
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

# gradient with respect to theta
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

# gradient with respect to the core tensor (gradient not using kronecker product at all)
gc = function(A_1,A_2,A_3,C,ttnsr,omega,type="ordinal"){
    g = gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega,type) ## first, take entrywise gradient w.r.t. theta
  
    g = ttl(as.tensor(g),list(t(A_1),t(A_2),t(A_3)),ms = 1:3)@data ## then, take entrywise gradient w.r.t. core tensor
  
    return(g)
}


####### update one factor matrix at a time while holding others fixed ###########
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

####### update core tensor #######
corecomb = function(A_1,A_2,A_3,C,ttnsr,omega,alpha=TRUE,type="ordinal"){
  h <- function(x) hc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega,type)
  g <- function(x) c(gc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega,type)) 
  d <- optim(c(C@data),h,g,method="BFGS") 
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  
  return(C)
}
                                              
## you can use this when W size is small
# corecomb = function(C,W,ttnsr,omega,alpha=TRUE){
#   Cvec <- c(C@data)
#   h <- function(x) h1(x,W,ttnsr,omega)
#   g <- function(x) g1(x,W,ttnsr,omega)
#   d <- optim(Cvec,h,g,method="BFGS")  ## seems BFGS is faster? (there was another option; what was that?)
#   C <- new("Tensor",C@num_modes,C@modes,data =d$par)
#   return(C)
# }

### main function ##
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
    while ((error > 10^-4)&(iter<50) ) {
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
    }
  
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$theta= theta
  result$iteration <- iter
  result$cost = cost
  result$omega=omega
  return(result)
}

likelihood = function(ttnsr,theta,omega,type="ordinal"){
    index=which(is.na(ttnsr)==F & is.na(theta)==F)
    ttnsr=ttnsr[index]
    theta=theta[index]
    
    if(type=="Gaussian") return(sqrt(sum((ttnsr-theta)^2)))
    k = length(omega)
    p=theta_to_p(theta,omega)
    l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
    return(sum(unlist(l)))
    
}

#### convert parameter tensor to probability tensor, based on logistic cumulative model
theta_to_p=function(theta,omega){
    k = length(omega)
    p = matrix(nrow = length(theta),ncol = k)
    for (i in 1:k) {
        p[,i] = as.numeric(logistic(theta + omega[i]))
    }
    p =  cbind(p,rep(1,length(theta)))-cbind(rep(0,length(theta)),p)
    p[p<epsilon]=epsilon
    return(p)
}

logistic = function(x){
    output=1/(1+exp(-x))
    return(output)
}

## BIC: Inputs d and r are vectors
bic = function(ttnsr,theta,omega,d,r){
    return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
}

## continous Tucker decomposition with possibly missing entries
fit_continuous=function(ttnsr,C,A_1,A_2,A_3,alpha = TRUE){  ## To allow input without alpha, I changed alpha = TRUE
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
    
    while ((error > 10^-4)&(iter<50) ) {
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
    }
    
    
    result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
    result$theta = theta
    result$iteration <- iter
    result$cost = cost
    return(result)
}


### special version of the function "fit_ordinal()"; ordinal matrix decomposition
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
    
    while ((error > 10^-4)&(iter<10) ) {
        iter = iter +1
        #update omega
        prevtheta <- A_1%*%t(A_2)
        
        currentomega=omega
        if(is.logical(omg)) {
            omega=tryCatch(polr(as.factor(c(ttnsr))~offset(-c(prevtheta)))$zeta,error=function(c)"omega cannot be reliably estimated",warning=function(c)"omega cannot be reliably estimated")
            if(inherits(omega,"numeric")==FALSE){
                print("Warning: omega cannot be estimated!"); break
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

##########################End of scripts#######################################
## Unused functions are removed from this script. Those functions are now saved in old_code/unused_functions.R

