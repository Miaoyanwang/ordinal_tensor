library(MASS)
library(rTensor)
library(pracma)
library(doParallel)
epsilon=10^-4



## A_1, A_2, A_3 
g1_factor=function(A_1,A_2,A_3,C,ttnsr,k,omega,type="ordinal"){
    g = gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega,type) ## first, take entrywise gradient w.r.t. theta
    
    if(k==1){
        g = ttl(as.tensor(g),list(matrix(1),t(A_2),t(A_3)),ms = 1:3) ## take gradient w.r.t. row in factor 1
        g = k_unfold(g,1)@data%*%t(k_unfold(C,1)@data)
    }
    else if(k==2){
        g = ttl(as.tensor(g),list(t(A_1),matrix(1),t(A_3)),ms = 1:3) ## take gradient w.r.t. row in factor 2
        g = k_unfold(g,2)@data%*%t(k_unfold(C,2)@data)
    }
    else if(k==3){
        g = ttl(as.tensor(g),list(t(A_1),t(A_2),matrix(1)),ms = 1:3) ## take gradient w.r.t. row in factor 3
        g = k_unfold(g,3)@data%*%t(k_unfold(C,3)@data)
    }
    return(g)
}

####### update one factor matrix at a time while holding others fixed ###########
factorcomb = function(A_1,A_2,A_3,C,ttnsr,k,omega,alpha=TRUE,type="ordinal"){
    if(k==1){
        if (alpha==TRUE) {
            l <- lapply(1:nrow(A_1),function(i){optim(A_1[i,],
                function(x) hc(t(as.matrix(x)),A_2,A_3,C,ttnsr[i,,],omega,type),
                function(x) g1_factor(t(as.matrix(x)),A_2,A_3,C,ttnsr[i,,],1,omega,type),method = "BFGS")$par})
            nA <- matrix(unlist(l),nrow = nrow(A_1),byrow = T)
        }else{
            W <- kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
            l <- lapply(1:nrow(A_1),function(i){constrOptim(A_1[i,],
                function(x) hc(t(as.matrix(x)),A_2,A_3,C,ttnsr[i,,],omega,type),
                function(x) g1_factor(t(as.matrix(x)),A_2,A_3,C,ttnsr[i,,],1,omega,type),ui = as.matrix(rbind(W,-W)),ci = rep(-alpha,2*nrow(W)),method = "BFGS")$par})
            nA <- matrix(unlist(l),nrow = nrow(A_1),byrow = T)
        }
        return(nA)
    }
    else if(k==2){
        if (alpha==TRUE) {
            l <- lapply(1:nrow(A_2),function(i){optim(A_2[i,],
                function(x) hc(A_1,t(as.matrix(x)),A_3,C,ttnsr[,i,],omega,type),
                function(x) g1_factor(A_1,t(as.matrix(x)),A_3,C,ttnsr[,i,],2,omega,type),method = "BFGS")$par})
            nA <- matrix(unlist(l),nrow = nrow(A_2),byrow = T)
        }else{
            W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
            l <- lapply(1:nrow(A_2),function(i){constrOptim(A_2[i,],
                function(x) hc(A_1,t(as.matrix(x)),A_3,C,ttnsr[,i,],omega,type),
                function(x) g1_factor(A_1,t(as.matrix(x)),A_3,C,ttnsr[,i,],2,omega,type),ui = as.matrix(rbind(W,-W)),ci = rep(-alpha,2*nrow(W)),method = "BFGS")$par})
            nA <- matrix(unlist(l),nrow = nrow(A_2),byrow = T)
        }
        return(nA)
    }
    else if(k==3){
        if (alpha==TRUE) {
            l <- lapply(1:nrow(A_3),function(i){optim(A_3[i,],
                function(x) hc(A_1,A_2,t(as.matrix(x)),C,ttnsr[,,i],omega,type),
                function(x) g1_factor(A_1,A_2,t(as.matrix(x)),C,ttnsr[,,i],3,omega,type),method = "BFGS")$par})
            nA <- matrix(unlist(l),nrow = nrow(A_3),byrow = T)
        }else{
            W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
            l <- lapply(1:nrow(A_3),function(i){constrOptim(A_3[i,],
                function(x) hc(A_1,A_2,t(as.matrix(x)),C,ttnsr[,,i],omega,type),
                function(x) g1_factor(A_1,A_2,t(as.matrix(x)),C,ttnsr[,,i],3,omega,type),ui = as.matrix(rbind(W,-W)),ci = rep(-alpha,2*nrow(W)),method = "BFGS")$par})
            nA <- matrix(unlist(l),nrow = nrow(A_3),byrow = T)
        }
        return(nA)
    }
}

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
        
        message=capture.output(decomp<-tucker_missing(newtheta,C,A_1,A_2,A_3,alpha))
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
    result$theta = theta
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
    result$theta=theta
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
    theta =W1%*%c(A_1)
    p = matrix(nrow = length(theta),ncol = k)
    for (i in 1:k) {
        p[,i] = as.numeric(logistic(theta + omega[i]))
    }
    q = matrix(nrow = length(theta),ncol = k+1)
    q[,1] = p[,1]*(1-p[,1])
    for (i in 2:k) {
        q[,i] =  p[,i]*(1-p[,i])+p[,i-1]*(1-p[,i-1])
    }
    q[,k+1] = p[,k]*(1-p[,k])
    l= lapply(1:(k+1),function(i) t(rbind(W1[which(c(ttnsr)==i),])*q[which(c(ttnsr)==i),i])%*%rbind(W1[which(c(ttnsr)==i),]))
    return(Reduce("+", l))
}


