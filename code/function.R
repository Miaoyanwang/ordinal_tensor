library(MASS)
library(rTensor)
library(pracma)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(Matrix)

########### simulate ordinal tensors based on logistic model
realization = function(tnsr,alpha){
  tnsr=as.tensor(tnsr)
  thet <- k_unfold(tnsr,1)@data
  theta1 <- thet + alpha[1]
  theta2 <- thet + alpha[2]
  result <- k_unfold(tnsr,1)@data
  p1 <- logistic(theta1)
  p2 <- logistic(theta2)-logistic(theta1)
  p3 <- matrix(1,nrow = nrow(thet),ncol = ncol(thet))-logistic(theta2)
  for (i in 1:nrow(thet)) {
    for(j in 1:ncol(thet)){
      result[i,j] <- sample(c(1,2,3),1,prob= c(p1[i,j],p2[i,j],p3[i,j]))
    }
  }
  return(k_fold(result,1,modes = tnsr@modes))
}      


########### Hessian w.r.t. A_1 ###########
Hessi = function(A_1,W4,ttnsr,omega){
  thet =W4%*%c(A_1)
  p1 = logistic(thet + omega[1])
  p2 = logistic(thet + omega[2])
  q1 = p1*(1-p1)
  q2 = p2*(1-p2)+p1*(1-p1)
  q3 = p2*(1-p2) 
  H = t(W4[which(c(ttnsr)==1),]*q1[which(c(ttnsr)==1)])%*%W4[which(c(ttnsr)==1),]+t(W4[which(c(ttnsr)==2),]*q2[which(c(ttnsr)==2)])%*%W4[which(c(ttnsr)==2),]+t(W4[which(c(ttnsr)==3),]*q3[which(c(ttnsr)==3)])%*%W4[which(c(ttnsr)==3),]
  return(H)
}

########### cost function ###########
h1 = function(A_1,W1,ttnsr,omega){
  thet =W1%*%c(A_1)
  p1 = logistic(thet + omega[1])
  p2 = logistic(thet + omega[2])
  p = cbind(p1,p2-p1,1-p2)
  return(-sum(log(c(p[which(c(ttnsr)==1),1],p[which(c(ttnsr)==2),2],p[which(c(ttnsr)==3),3]))))
}

########### gradient ###########
g1 = function(A_1,W1,ttnsr,omega){
  thet =W1%*%c(A_1)
  p1 = logistic(thet + omega[1])
  p2 = logistic(thet + omega[2])
  q1 <- p1-1
  q2 <- (p2*(1-p2)-p1*(1-p1))/(p1-p2)
  q3 <- p2
  gd = apply(as.matrix(W1[which(c(ttnsr)==1),])*q1[which(c(ttnsr)==1)],2,sum)+apply(as.matrix(W1[which(c(ttnsr)==2),])*q2[which(c(ttnsr)==2)],2,sum)+apply(as.matrix(W1[which(c(ttnsr)==3),])*q3[which(c(ttnsr)==3)],2,sum)
  
    return(gd)
}

####### update a factor matrix at one time while holding others fixed ###########
comb = function(A,W,ttnsr,k,omega,alph=TRUE){
  nA = A
  tnsr1 <- k_unfold(as.tensor(ttnsr),k)@data
  if (alph==TRUE) {
      l <- lapply(1:nrow(A),function(i){optim(A[i,],function(x) h1(x,W,tnsr1[i,],omega),function(x) g1(x,W,tnsr1[i,],omega),method = "BFGS")$par})
    
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }else{
    l <- lapply(1:nrow(A),function(i){constrOptim(A[i,],
                                                  function(x) h1(x,W,tnsr1[i,],omega),function(x) g1(x,W,tnsr1[i,],omega),
                                                  ui = as.matrix(rbind(W,-W)),ci = rep(-alph,2*nrow(W)),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }
  return(nA)
}

####### update core tensor #######
corecomb = function(C,W,ttnsr,omega,alph=TRUE){
  Cvec <- c(C@data)
  h <- function(x) h1(x,W,ttnsr,omega)
  g <- function(x) g1(x,W,ttnsr,omega)
  H <- function(x) Hessi(x,W,ttnsr,omega)
  d <- nlminb(Cvec,h,g,H) 
  ##d <- optim(Cvec,h,g,method="BFGS")  ## seems BFGS is faster??
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)

  return(C)
}

####### update core tensor #######
prevcorecomb = function(C,W,ttnsr,omega,alph=TRUE){
  Cvec <- c(C@data)
  h <- function(x) h1(x,W,ttnsr,omega)
  g <- function(x) g1(x,W,ttnsr,omega)
  H <- function(x) Hessi(x,W,ttnsr,omega)
  
  
  if (alph==TRUE) {
    d <- nlminb(Cvec,h,g,H)
    C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  }else{
    d <- constrOptim(Cvec,h,g,ui = rbind(W,-W),ci = rep(-alph,2*nrow(W)),method = "BFGS")
    C <- new("Tensor",C@num_modes,C@modes,data =d$par)
  }
  return(C)
}


########ordinal tensor decomposition based on Tucker structure #######
fit_ordinal = function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph = TRUE){
    alphbound <- alph+10^-4
    result = list()
    error<- 3
    iter = 0
    cost=NULL
    if (alph == TRUE) {
        while ((error > 10^-4)&(iter<50) ) {
            iter = iter +1
            
            #update omega
            prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
            if(sum(omega)==TRUE) omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
            prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
            
            # update C
            W4 <- kronecker(kronecker(A_3,A_2),A_1)
            C <- corecomb(C,W4,c(ttnsr),omega)
            
            #update A_1
            W1 = kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
            A_1 <- comb(A_1,W1,ttnsr,1,omega)
            #orthognalize A_1
            qr_res=qr(A_1)
            A_1=qr.Q(qr_res)
            C=ttm(C,qr.R(qr_res),1)
            
            # update A_2
            W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
            A_2 <- comb(A_2,W2,ttnsr,2,omega)
            #orthognalize A_2
            qr_res=qr(A_2)
            A_2=qr.Q(qr_res)
            C=ttm(C,qr.R(qr_res),2)
            
            # update A_3
            W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
            A_3 <- comb(A_3,W3,ttnsr,3,omega)
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
            if(sum(omega)==TRUE) omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
            prev <- likelihood(ttnsr[ttnsr>0],prevtheta[ttnsr>0],omega)
            
            #update A_1
            W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
            A_1 <- comb(A_1,W1,ttnsr,1,omega,alphbound)
            if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
            #orthognalize A_1
            qr_res=qr(A_1)
            A_1=qr.Q(qr_res)
            C=ttm(C,qr.R(qr_res),1)
            
            
            # update A_2
            W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
            A_2 <- comb(A_2,W2,ttnsr,2,omega,alphbound)
            if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
            #orthognalize A_2
            qr_res=qr(A_2)
            A_2=qr.Q(qr_res)
            C=ttm(C,qr.R(qr_res),2)


            # update A_3
            W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
            A_3 <- comb(A_3,W3,ttnsr,3,omega,alphbound)
            if(max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))>=alph) break
            #orthognalize A_3
            qr_res=qr(A_3)
            A_3=qr.Q(qr_res)
            C=ttm(C,qr.R(qr_res),3)

            
            # update C
            W4 <- kronecker(kronecker(A_3,A_2),A_1)
            C <- corecomb(C,W4,c(ttnsr),omega)
            
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
    if (alph == TRUE) {
        while ((error > 10^-4)&(iter<50) ) {
            iter = iter +1
            
            prevtheta <- tensorize(A_1,A_2,A_3)
            #update omega
            if(sum(omega)==TRUE) omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
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
            if(sum(omega)==TRUE) omega <- polr(as.factor(c(ttnsr[ttnsr>0]))~offset(-c(prevtheta[ttnsr>0])))$zeta
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
    p1 = logistic(c(thet) + alpha[1])
    p2 = logistic(c(thet) + alpha[2])
    p = cbind(p1,p2-p1,1-p2)
    return(-sum(log(c(p[which(c(ttnsr)==1),1],p[which(c(ttnsr)==2),2],p[which(c(ttnsr)==3),3]))))
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

########################################################################### previous code ################################################
########################
#fit_ordinal2 = function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph = TRUE){
#omega <- sort(rnorm(2))
  #alphbound <- alph+10^-4
  #result = list()
  #error<- 3
  #iter = 0
  #d1 <- nrow(A_1); d2 <- nrow(A_2); d3 <- nrow(A_3)
  #r1 <- ncol(A_1); r2 <- ncol(A_2); r3 <- ncol(A_3)
  #if (alph == TRUE) {
      #while ((error > 10^-4)&(iter<50) ) {
        #iter = iter +1
      
      #update A_1
      #prevtheta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      #prev <- likelihood(ttnsr,prevtheta,omega)
      #W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      #A_1 <- comb(A_1,W1,ttnsr,1,omega)
      
      
      # update A_2
      #W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      #A_2 <- comb(A_2,W2,ttnsr,2,omega)
      
      # update A_3
      #W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      #A_3 <- comb(A_3,W3,ttnsr,3,omega)
      
      # update C
      #W4 <- kronecker(kronecker(A_3,A_2),A_1)
      #C <- corecomb(C,W4,c(ttnsr),omega)

      #update omega
      #theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      #m <- polr(as.factor(c(ttnsr))~offset(-c(theta@data)))
      #omega <- m$zeta
      
      
      
      #theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      # new <- likelihood(ttnsr,theta,omega)
      # error <- abs((new-prev)/prev)
      # }
    #}else{
      #while ((error > 10^-4)&(iter<50) ) {
        #iter = iter +1
      
      #update A_1
      #prevtheta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      #prev <- likelihood(ttnsr,prevtheta,omega)
      #W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      #A_1 <- comb(A_1,W1,ttnsr,1,omega,alphbound)
      #if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      
      # update A_2
      #W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      #A_2 <- comb(A_2,W2,ttnsr,2,omega,alphbound)
      #if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      # update A_3
      #W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      #A_3 <- comb(A_3,W3,ttnsr,3,omega,alphbound)
      #if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      # update C
      #W4 <- kronecker(kronecker(A_3,A_2),A_1)
      #C <- corecomb(C,W4,c(ttnsr),omega)
      #if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      #update omega
      #theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      #m <- polr(as.factor(c(ttnsr))~offset(-c(theta@data)))
      #omega <- m$zeta
      
      
      #theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      #new <- likelihood(ttnsr,theta,omega)
      #error <- abs((new-prev)/prev)
      #}
    #}
  
  #result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  #result$iteration <- iter; result$omega <- omega
  #return(result)
  #}




