library(MASS)
library(rTensor)
library(pracma)
library(doParallel)



########### simulation ordinal tensors based on logist model with arbitrary k
realization = function(theta,omega){
  k = length(omega)
  theta=as.tensor(theta)
  thet <- c(theta@data)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = logistic(thet + omega[i])
  }
  p =  cbind(p,rep(1,length(thet)))-cbind(rep(0,length(thet)),p)
  for (j in 1:length(thet)) {
    thet[j] <-  sample(1:(k+1),1,prob = p[j,])
  }
  return(as.tensor(array(thet,dim =theta@modes)))
}      




########### cost function  with arbitrary k ###########
h1 = function(A_1,W1,ttnsr,omega){
  k = length(omega)
  thet =W1%*%c(A_1)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(thet + omega[i]))
  }
  p =  cbind(p,rep(1,length(thet)))-cbind(rep(0,length(thet)),p)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
}

hc = function(A_1,A_2,A_3,C,ttnsr,omega){
  k = length(omega)
  thet = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(thet + omega[i]))
  }
  p =  cbind(p,rep(1,length(thet)))-cbind(rep(0,length(thet)),p)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
}

# hc(A_1,A_2,A_3,C,ttnsr,omega)
########### gradient with arbitrary k ###########
g1 = function(A_1,W,ttnsr,omega){
  k = length(omega)
  thet =W%*%c(A_1)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(thet + omega[i]))
  }
  q = matrix(nrow = length(thet),ncol = k+1)
  q[,1] <- p[,1]-1
   if(k>=2){
  for (i in 2:k) {
      #q[,i] <-  (p[,i]*(1-p[,i])-p[,i-1]*(1-p[,i-1]))/(p[,i-1]-p[,i])
      q[,i] <- p[,i]+p[,i-1]-1
  }
   }
  q[,k+1] <- p[,k]
  l <- Reduce("+",lapply(1:(k+1),function(i) apply(rbind(W[which(c(ttnsr)==i),])*q[which(c(ttnsr)==i),i],2,sum)))
  return(l)
}


gradient_tensor=function(A_1,A_2,A_3,C,ttnsr,omega){
    k = length(omega)
    thet = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
    p = matrix(nrow = length(thet),ncol = k)
    for (i in 1:k) {
        p[,i] = as.numeric(logistic(thet + omega[i]))
    }
    q = matrix(nrow = length(thet),ncol = k+1)
    q[,1] <- p[,1]-1
    if(k>=2){
    for (i in 2:k) {
        ##q[,i] <-  (p[,i]*(1-p[,i])-p[,i-1]*(1-p[,i-1]))/(p[,i-1]-p[,i])
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

gc = function(A_1,A_2,A_3,C,ttnsr,omega){
  k = length(omega)
  thet = c(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(thet + omega[i]))
  }
  q = matrix(nrow = length(thet),ncol = k+1)
  q[,1] <- p[,1]-1
  for (i in 2:k) {
      #q[,i] <-  (p[,i]*(1-p[,i])-p[,i-1]*(1-p[,i-1]))/(p[,i-1]-p[,i])
    q[,i] <- p[,i]+p[,i-1]-1
  }
  q[,k+1] <- p[,k]
  W = kronecker(A_2,A_1)
  d = c(nrow(A_1),nrow(A_2),nrow(A_3))
  cl <- makeCluster(20)
  registerDoParallel(cl)
  l <- foreach(j = 1:d[3],.combine = "+") %dopar% {
    Reduce("+",lapply(1:(k+1),function(i) apply(rbind(kronecker(A_3[j,,drop = F],W)[which(c(ttnsr)[(d[1]*d[2]*(j-1)+1):(d[1]*d[2]*(j-1)+d[1]*d[2])]==i),])*
                                                  q[which(c(ttnsr)[(d[1]*d[2]*(j-1)+1):(d[1]*d[2]*(j-1)+d[1]*d[2])]==i)+d[1]*d[2]*(j-1),i],2,sum)))
  }
  stopCluster(cl)
  return(l)
}



####### update a factor matrix at one time while holding others fixed ###########
comb = function(A,W,ttnsr,k,omega,alph=TRUE){
  nA = A
  tnsr1 <- k_unfold(as.tensor(ttnsr),k)@data
  if (alph==TRUE) {
l <- lapply(1:nrow(A),function(i){optim(A[i,],function(x) h1(x,W,tnsr1[i,],omega),function(x) g1(x,W,tnsr1[i,],omega),method = "BFGS")$par})
    
    
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }else{
l <- lapply(1:nrow(A),function(i){constrOptim(A[i,],function(x) h1(x,W,tnsr1[i,],omega),function(x) g1(x,W,tnsr1[i,],omega),ui = as.matrix(rbind(W,-W)),ci = rep(-alph,2*nrow(W)),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }
  return(nA)
}

####### update core tensor #######
corecomb = function(A_1,A_2,A_3,C,ttnsr,omega,alph=TRUE){
  h <- function(x) hc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega)
  #g <- function(x) gc(A_1,A_2,A_3,new("Tensor",C@num_modes,C@modes,data = x),ttnsr,omega) ## take longe
  #d <- optim(c(C@data),h,g,method="BFGS") 
  d <- optim(c(C@data),h,method="BFGS")  ## skip the gradient calculation
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  
  return(C)
}


# corecomb = function(C,W,ttnsr,omega,alph=TRUE){
#   Cvec <- c(C@data)
#   h <- function(x) h1(x,W,ttnsr,omega)
#   g <- function(x) g1(x,W,ttnsr,omega)
#   d <- optim(Cvec,h,g,method="BFGS")  ## seems BFGS is faster??
#   C <- new("Tensor",C@num_modes,C@modes,data =d$par)
#   return(C)
# }

##### Minorize-Maximization scheme for updating factor matrices -- much faster than alternating minimization
fit_ordinal_MM=function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph=TRUE){
    alphbound <- alph+10^-4
    result = list()
    error<- 3
    iter = 0
    cost=NULL
    omg = omega
    
    if (alph == TRUE) {
        while ((error > 10^-4)&(iter<200) ) {
            iter = iter +1
            prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
            if(sum(omg)==TRUE) {
                omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
            }
            prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
        
            ## update C, A_1, A_2, A_3 all together
            newtheta=prevtheta-4*gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega)/length(omega)
            message=capture.output(decomp<-tucker(as.tensor(newtheta),r))
            A_1=decomp$U[[1]]
            A_2=decomp$U[[2]]
            A_3=decomp$U[[3]]
            C=decomp$Z
        
        
            theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
            new <- likelihood(ttnsr[ttnsr>0],theta[ttnsr>0],omega)
            cost = c(cost,new)
            (error <- abs((new-prev)/prev))
            message(paste(iter,"-th  iteration -- cost function is ",new," -----------------"))
        }
    }### To-do: need to add the option (alph != TRUE)
    
        result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
        result$iteration <- iter
        result$cost = cost; result$omega=omega
        return(result)
        
}
        
fit_ordinal = function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph = TRUE){
  alphbound <- alph+10^-4
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omg = omega
  
  if (alph == TRUE) {
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      #update omega
      prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      if(sum(omg)==TRUE) {
        omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
      }
      
   
      prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
      
      # update C
      C <- corecomb(A_1,A_2,A_3,C,ttnsr,omega)
      #update A_1
    

      W <-kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W,ttnsr,1,omega)
      #orthognalize A_1
      qr_res=qr(A_1)
      A_1=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),1)
      
      # update A_2
      W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      A_2 <- comb(A_2,W,ttnsr,2,omega)
      #orthognalize A_2
      qr_res=qr(A_2)
      A_2=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),2)
      
      # update A_3
      W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      A_3 <- comb(A_3,W,ttnsr,3,omega)
      #orthognalize A_3
      qr_res=qr(A_3)
      A_3=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),3)
      
      theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      new <- likelihood(ttnsr[ttnsr>0],theta[ttnsr>0],omega)
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
  }else{
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      #update omega
      prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      if(sum(omg)==TRUE) {
        omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
      }
      
      
      prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
      
      # update C
      C <- corecomb(A_1,A_2,A_3,C,ttnsr,omega)
      #update A_1
      W =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W,ttnsr,1,omega,alphbound)
      if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
      #orthognalize A_1
      qr_res=qr(A_1)
      A_1=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),1)
      
      
      # update A_2
      W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      A_2 <- comb(A_2,W,ttnsr,2,omega,alphbound)
      if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
      #orthognalize A_2
      qr_res=qr(A_2)
      A_2=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),2)
      
      
      # update A_3
      W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      A_3 <- comb(A_3,W,ttnsr,3,omega,alphbound)
      if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
      #orthognalize A_3
      qr_res=qr(A_3)
      A_3=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),3)
      
      
      
      
      theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
      new <- likelihood(ttnsr[ttnsr>0],theta[ttnsr>0],omega)
      cost = c(cost,new)
      error <- abs((new-prev)/prev)
      if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
    }
  }
  
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$iteration <- iter
  result$cost = cost; result$omega=omega
  return(result)
}


####### ordinal tensor decomposition based on CP structure #######
fit_ordinal_cp=function(ttnsr,A_1,A_2,A_3,omega=TRUE,alph = TRUE){
  alphbound <- alph+10^-4
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omg = omega
   
  if (alph == TRUE) {
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      prevtheta <- tensorize(A_1,A_2,A_3)
      #update omega
      if(sum(ome)==TRUE) omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
      prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
      
      #update A_1
      W1 = KhatriRao(A_3,A_2)
      A_1 <- comb(A_1,W1,ttnsr,1,omega)
      A_1 = apply(A_1,2,function(x){x/norm(x,"2")})
      
      # update A_2
      W2 <- KhatriRao(A_3,A_1)
      A_2 <- comb(A_2,W2,ttnsr,2,omega)
      A_2 = apply(A_2,2,function(x){x/norm(x,"2")})
      
      # update A_3
      W3 <-  KhatriRao(A_2,A_1)
      A_3 <- comb(A_3,W3,ttnsr,3,omega)
      
      
      theta <- tensorize(A_1,A_2,A_3)
      new <- likelihood(ttnsr[ttnsr>0],theta[ttnsr>0],omega)
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
  }else{
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      prevtheta <- tensorize(A_1,A_2,A_3)
      #update omega
      if(sum(ome)==TRUE) omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
      prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
      
      #update A_1
      W1 = KhatriRao(A_3,A_2)
      A_1 <- comb(A_1,W1,ttnsr,1,omega,alphbound)
      if(max(abs(tensorize(A_1,A_2,A_3)))>=alph) break
      
      
      # update A_2
      W2 <- KhatriRao(A_3,A_1)
      A_2 <- comb(A_2,W2,ttnsr,2,omega,alphbound)
      if(max(abs(tensorize(A_1,A_2,A_3)))>=alph) break
      
      # update A_3
      W3 <-  KhatriRao(A_2,A_1)
      A_3 <- comb(A_3,W3,ttnsr,3,omega,alphbound)
      if(max(abs(tensorize(A_1,A_2,A_3)))>=alph) break
      
      pre=rescale(A_1,A_2,A_3)
      A_1=pre$A_1
      A_2=pre$A_2
      A_3=pre$A_3
      
      
      theta <- tensorize(A_1,A_2,A_3)
      new <- likelihood(ttnsr[ttnsr>0],theta[ttnsr>0],omega)
      cost = c(cost,new)
      error <- abs((new-prev)/prev)
      if(max(abs(tensorize(A_1,A_2,A_3)))>=alph) break
    }
  }
  
  result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$iteration <- iter
  result$cost = cost; result$omega=omega
  return(result)
}

likelihood = function(ttnsr,thet,alpha){
  k = length(alpha)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = as.numeric(logistic(thet + alpha[i]))
  }
  p =  cbind(p,rep(1,length(thet)))-cbind(rep(0,length(thet)),p)
  l = lapply(1:(k+1),function(i) -log(p[which(c(ttnsr)==i),i]))
  return(sum(unlist(l)))
  
}


logistic = function(x){
  return(1/(1+exp(-x)))
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


## BIC: Inputs d and r are vectors. 
bic = function(ttnsr,theta,omega,d,r){
  return(2*likelihood(ttnsr,theta,omega)+(prod(r)+sum(r*(d-r)))*log(prod(d)))
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


