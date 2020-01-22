library(MASS)
library(rTensor)
library(pracma)
library(doParallel)
epsilon=10^-4

########### simulation ordinal tensors based on logist model with arbitrary k
realization = function(theta,omega){
  theta=as.tensor(theta)
  thet <- c(theta@data)
  p = theta_to_p(thet,omega)
  for (j in 1:length(thet)) {
    thet[j] <-  sample(1:(k+1),1,prob = p[j,])
  }
  return(as.tensor(array(thet,dim =theta@modes)))
}      


estimation = function(theta,omega,type){
  if(is.matrix(theta)){
  thet=c(theta)
  }else{
  theta=as.tensor(theta)
  thet <- c(theta@data)
  }
  p = theta_to_p(thet,omega)
  if(type=="max"){  # score prediction based on the mode 
  for (j in 1:length(thet)) thet[j] <-  which.max(p[j,])
  }else if(type=="mean"){ # score prediction based on the mean 
  for (j in 1:length(thet)) thet[j] <-  sum(p[j,]*(1:(k+1))) ## why rounding to integer?
  }else if(type=="median"){# score prediction based on the median 
  for (j in 1:length(thet)) thet[j] <-  which(c(cumsum(p[i,]),1)>=0.5)[1] ## median
  }
  if(is.matrix(theta)) return(as.matrix(thet,dim=dim(theta)))
  else return(as.tensor(array(thet,dim =theta@modes)))
}   

########### cost function with arbitrary k ###########
h1 = function(A_1,W1,ttnsr,omega,type="ordinal"){
  k = length(omega)
  thet =W1%*%c(A_1)
  if(type=="Gaussian"){
      l=sqrt(sum((thet[is.na(ttnsr)==F]-ttnsr[is.na(ttnsr)==F])^2))
      return(l)  
  }
  p = theta_to_p(thet,omega)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
}

########### cost function based on Tucker representation ###########
hc = function(A_1,A_2,A_3,C,ttnsr,omega,type="ordinal"){
  k = length(omega)
  thet = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
  if(type=="Gaussian"){
      l=sqrt(sum((thet[is.na(ttnsr)==F]-ttnsr[is.na(ttnsr)==F])^2))
      return(l)
  }
  p = theta_to_p(thet,omega)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
}

########### gradient with arbitrary k ###########
g1 = function(A_1,W,ttnsr,omega,type="ordinal"){
  k = length(omega)
  thet =W%*%c(A_1)
  if(type=="Gaussian"){
      pretheta=thet-as.matrix(ttnsr)
      pretheta[is.na(pretheta)==T]=0
      l=2*t(pretheta)%*%W
      return(l)
  }
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(thet + omega[i]))
  }
  q = matrix(nrow = length(thet),ncol = k+1)
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
    thet = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
    if(type=="Gaussian"){
        pretheta=thet-ttnsr
        pretheta[is.na(pretheta)==T]=0
        l=2*pretheta
        return(l)
    }
    p = matrix(nrow = length(thet),ncol = k)
    for (i in 1:k) {
        p[,i] = as.numeric(logistic(thet + omega[i]))
    }
    q = matrix(nrow = length(thet),ncol = k+1)
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

# gradient with respect to core (gradient not using kronecker product at all)
gc = function(A_1,A_2,A_3,C,ttnsr,omega){
    g = gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega) ## first, take entrywise gradient w.r.t. theta
  
    g = ttl(as.tensor(g),list(t(A_1),t(A_2),t(A_3)),ms = 1:3) ## then, take entrywise gradient w.r.t. core tensor
  
    return(g)
}
     
            

####### update a factor matrix at one time while holding others fixed ###########
comb = function(A,W,ttnsr,k,omega,alph=TRUE,type="ordinal"){
  nA = A
  tnsr1 <- k_unfold(as.tensor(ttnsr),k)@data
  if (alph==TRUE) {
l <- lapply(1:nrow(A),function(i){optim(A[i,],function(x) h1(x,W,tnsr1[i,],omega,type),function(x) g1(x,W,tnsr1[i,],omega,type),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }else{
l <- lapply(1:nrow(A),function(i){constrOptim(A[i,],function(x) h1(x,W,tnsr1[i,],omega,type),function(x) g1(x,W,tnsr1[i,],omega,type),ui = as.matrix(rbind(W,-W)),ci = rep(-alph,2*nrow(W)),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }
  return(nA)
}

####### update core tensor #######
corecomb = function(A_1,A_2,A_3,C,ttnsr,omega,alph=TRUE){
  h <- function(x) hc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega)
  g <- function(x) c(gc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega)@data) 
  d <- optim(c(C@data),h,g,method="BFGS") 
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  
  return(C)
}
                                              
## you can use this when W size is small
# corecomb = function(C,W,ttnsr,omega,alph=TRUE){
#   Cvec <- c(C@data)
#   h <- function(x) h1(x,W,ttnsr,omega)
#   g <- function(x) g1(x,W,ttnsr,omega)
#   d <- optim(Cvec,h,g,method="BFGS")  ## seems BFGS is faster? (there was another option; what was that?)
#   C <- new("Tensor",C@num_modes,C@modes,data =d$par)
#   return(C)
# }


fit_ordinal = function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph = TRUE){
    
  #alphbound <- alph+10^-4
  if(is.logical(alph)) alpha_minus=alpha_minus2=TRUE
  else{
        alpha_minus=alph-epsilon
        alpha_minus2=alph-2*epsilon
  }
    
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omg = omega
  
    while ((error > 10^-4)&(iter<10) ) {
      iter = iter +1
      #update omega
      prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      
      currentomega=omega
      if(is.logical(omg)) {
          omega=tryCatch(polr(as.factor(c(ttnsr))~offset(-c(prevtheta)))$zeta,error=function(c)"omega cannot be reliably estimated",warning=function(c)"omega cannot be reliably estimated")
          if(inherits(omega,"numeric")==FALSE){
              print("Warning: omega cannot be estimated!"); break
          }
          currentomega=omega
      }
   
      prev <- likelihood(ttnsr,prevtheta,omega)
      
    
      #update A_1
      W <-kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W,ttnsr,1,omega,alpha_minus2)
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
      A_3 <- comb(A_3,W,ttnsr,3,omega,alph)
      #orthognalize A_3
      qr_res=qr(A_3)
      A_3=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),3)
      
      # update C
      C_update <- corecomb(A_1,A_2,A_3,C,ttnsr,omega)
      
      if(!is.logical(alph) & max(abs(ttl(C_update,list(A_1,A_2,A_3),ms=1:3)@data))>=alpha_minus2){
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
  result$iteration <- iter
  result$cost = cost; result$omega=omega
  return(result)
}

likelihood = function(ttnsr,thet,alpha,type="ordinal"){
    index=which(is.na(ttnsr)==F & is.na(thet)==F)
    ttnsr=ttnsr[index]
    thet=thet[index]
    
    if(type=="Gaussian") return(sqrt(sum((ttnsr-thet)^2)))
    
    k = length(alpha)
    p = matrix(nrow = length(thet),ncol = k)
    for (i in 1:k) {
        p[,i] = as.numeric(logistic(thet + alpha[i]))
    }
    p =  cbind(p,rep(1,length(thet)))-cbind(rep(0,length(thet)),p)
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
    return(p)
}

logistic = function(x){
    return(1/(1+exp(-x)))
}

## BIC: Inputs d and r are vectors. 
bic = function(ttnsr,theta,omega,d,r){
    return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
}

## continous Tucker decomposition with possibly missing data
tucker_missing=function(ttnsr,r,alpha,ini=TRUE){
    if(is.logical(alpha)) alpha_minus=alpha_minus2=TRUE
    else{
        alpha_minus=alpha-epsilon
        alpha_minus2=alpha-2*epsilon
    }
    
    if(is.logical(ini)){
        d=dim(ttnsr)
        A_1 = randortho(d[1])[,1:r[1]]
        A_2 = randortho(d[2])[,1:r[2]]
        A_3 = randortho(d[3])[,1:r[3]]
        C = rand_tensor(modes = r)
    }else{
        A_1=ini[[1]];A_2=ini[[2]];A_3=ini[[3]];C=ini[[4]]
    }
    
    if(!is.logical(alpha) & max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=(alpha-3*epsilon)){
        C=C*(alpha-3*epsilon)/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))
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
        A_1 <- comb(A_1,W,ttnsr,1,omega,alph= alpha_minus2,type="Gaussian")
        #orthognalize A_1
        qr_res=qr(A_1)
        A_1=qr.Q(qr_res)
        C=ttm(C,qr.R(qr_res),1)
        
        # update A_2
        W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
        A_2 <- comb(A_2,W,ttnsr,2,omega,alph= alpha_minus,type="Gaussian")
        #orthognalize A_2
        qr_res=qr(A_2)
        A_2=qr.Q(qr_res)
        C=ttm(C,qr.R(qr_res),2)
        
        # update A_3
        W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
        A_3 <- comb(A_3,W,ttnsr,3,omega,alph= alpha,type="Gaussian")
        #orthognalize A_3
        qr_res=qr(A_3)
        A_3=qr.Q(qr_res)
        C=ttm(C,qr.R(qr_res),3)
        
        ## alternating
        #C_update <- corecomb(A_1,A_2,A_3,C,ttnsr,omega,type)
        W=kronecker(kronecker(A_3,A_2),A_1)
        fit=lm(c(ttnsr)~-1+W)
        C_update=as.tensor(array(fit$coef,dim=r))
        
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
    
    
    result$Z <- C; result$U=list(A_1,A_2,A_3)
    result$iteration <- iter
    result$cost = cost
    return(result)
}


### special version of the function "fit_ordinal()"; ordinal matrix decomposition
fit_ordinal_matrix = function(ttnsr,A_1,A_2,omega=TRUE,alph = TRUE){
    
    if(is.logical(alph)) alpha_minus=TRUE
    else alpha_minus=alpha-epsilon
    
    if(!is.logical(alph) & max(A_1%*%t(A_2))>= alph-2*epsilon){
        A_2=A_2*(alph-2*epsilon)/max(A_1%*%t(A_2))
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
        A_2 <- comb(A_2,A_1,ttnsr,2,omega,alph = alpha_minus)
        
        A_1_update <- comb(A_1,A_2,ttnsr,1,omega,alph = alpha)
        
        
      
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
    return(result)
}

##########################End of codes#######################################
#########Functions below are used in the current version. I did not revise them carefully..#############################


##### Minorize-Maximization scheme for updating factor matrices -- much faster than alternating minimization
fit_ordinal_MM=function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alpha=TRUE){
    
    result = list()
    error<- 3
    iter = 0
    cost=NULL
    omg = omega
    
    while ((error > epsilon)&(iter<50) ) {
        iter = iter +1
        prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
        currentomega=omega
        if(is.logical(omg)) {
            omega=tryCatch(polr(as.factor(c(ttnsr))~offset(-c(prevtheta)))$zeta,error=function(c)"omega cannot be reliably estimated",warning=function(c)"omega cannot be reliably estimated")
            if(inherits(omega,"numeric")==FALSE){
                print("Warning: omega cannot be estimated!"); break
            }
            currentomega=omega
        }
        
        prev <- likelihood(ttnsr,prevtheta,omega)
        
        newtheta=prevtheta-2*gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega)
        message=capture.output(decomp<-tucker_missing(newtheta,r,alpha,ini=list(A_1,A_2,A_3,C)))
        A_1=decomp$U[[1]]
        A_2=decomp$U[[2]]
        A_3=decomp$U[[3]]
        C=decomp$Z
        
        theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
        new <- likelihood(ttnsr,theta,omega)
        cost = c(cost,new)
        (error <- abs((new-prev)/prev))
        message(paste(iter,"-th  iteration -- cost function is ",new," -----------------"))
    }
    result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
    result$iteration <- iter
    result$cost = cost; result$omega=currentomega
    return(result)
}


####### ordinal tensor decomposition based on CP structure #######
fit_ordinal_cp=function(ttnsr,A_1,A_2,A_3,omega=TRUE,alpha = TRUE){
    
    if(is.logical(alpha)) alpha_minus=alpha_minus2=TRUE
    else{
        alpha_minus=alpha-epsilon
        alpha_minus2=alpha-2*epsilon
    }
    
    result = list()
    error<- 3
    iter = 0
    cost=NULL
    omg = omega
    
    while ((error > 10^-4)&(iter<50) ) {
        iter = iter +1
        
        prevtheta <- tensorize(A_1,A_2,A_3)
        #update omega
        if(is.logical(ome)) omega <- polr(as.factor(c(ttnsr))~offset(-c(prevtheta)))$zeta
        prev <- likelihood(ttnsr,prevtheta,omega)
        
        #update A_1
        W1 = KhatriRao(A_3,A_2)
        A_1 <- comb(A_1,W1,ttnsr,1,omega,alpha_minus2)
        
        
        # update A_2
        W2 <- KhatriRao(A_3,A_1)
        A_2 <- comb(A_2,W2,ttnsr,2,omega,alpha_minus)
        
        # update A_3
        W3 <-  KhatriRao(A_2,A_1)
        A_3_update <- comb(A_3,W3,ttnsr,3,omega,alpha)
        
        if(!is.logical(alpha) & max(abs(tensorize(A_1,A_2,A_3)))>alpha_minus2){
            pre=rescale(A_1,A_2,A_3)
            A_1=pre$A_1
            A_2=pre$A_2
            A_3=pre$A_3
            theta <- tensorize(A_1,A_2,A_3)
            new <- likelihood(ttnsr,theta,omega)
            cost = c(cost,new)
            error <- abs((new-prev)/prev); break 
        }else{
            A_3=A_3_update
            
            pre=rescale(A_1,A_2,A_3)
            A_1=pre$A_1
            A_2=pre$A_2
            A_3=pre$A_3
            theta <- tensorize(A_1,A_2,A_3)
            new <- likelihood(ttnsr,theta,omega)
            cost = c(cost,new)
            error <- abs((new-prev)/prev)
        }   
    }
    
    result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
    result$iteration <- iter
    result$cost = cost; result$omega=omega
    return(result)
}

##### construct tensor from CP factors 
tensorize=function(A_1,A_2,A_3){
  r=ncol(A_1)
  tensor=0
  for(i in 1:r){
    tensor=tensor+A_1[,i]%o%A_2[,i]%o%A_3[,i]
  }
  return(tensor)
}


rescale=function(A_1,A_2,A_3){
  r=ncol(A_1)
  for(i in 1:r){
    A_3[,i]=A_3[,i]*sqrt(sum(A_1[,i]^2))*sqrt(sum(A_2[,i]^2))
    A_2[,i]=A_2[,i]/sqrt(sum(A_2[,i]^2))
    A_1[,i]=A_1[,i]/sqrt(sum(A_1[,i]^2))     
  }
  return(list("A_1"=A_1,"A_2"=A_2,"A_3"=A_3)) 
}


########### Hessian w.r.t. A_1 = C with W1 = kronecker(A_1,A_2,A_3) with arbitrary k ###########
Hessi = function(A_1,W1,ttnsr,omega){
    k = length(omega)
    thet =W1%*%c(A_1)
    p = matrix(nrow = length(thet),ncol = k)
    for (i in 1:k) {
        p[,i] = as.numeric(logistic(thet + omega[i]))
    }
    q = matrix(nrow = length(thet),ncol = k+1)
    q[,1] = p[,1]*(1-p[,1])
    for (i in 2:k) {
        q[,i] =  p[,i]*(1-p[,i])+p[,i-1]*(1-p[,i-1])
    }
    q[,k+1] = p[,k]*(1-p[,k])
    l= lapply(1:(k+1),function(i) t(rbind(W1[which(c(ttnsr)==i),])*q[which(c(ttnsr)==i),i])%*%rbind(W1[which(c(ttnsr)==i),]))
    return(Reduce("+", l))
}


