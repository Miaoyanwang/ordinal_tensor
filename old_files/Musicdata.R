source("functions.R")
load("InCar_Music.RData")

## Choosing argmax index
expectation = function(theta,omega){
  k = length(omega)
  theta=as.tensor(theta)
  thet <- c(theta@data)
  p = matrix(nrow = length(thet),ncol = k)
  for (i in 1:k) {
    p[,i] = logistic(thet + omega[i])
  }
  p =  cbind(p,rep(1,length(thet)))-cbind(rep(0,length(thet)),p)
  for (j in 1:length(thet)) {
    thet[j] <-  which.max(p[j,])
  }
  return(as.tensor(array(thet,dim =theta@modes)))
}      


fit_missing = function(tensor,ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph){
  alphbound <- alph+10^-4
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  omg = omega
  
  while ((error > 10^-4)&(iter<50) ) {
    iter = iter +1
    
    #update omega
    theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    if(sum(omg)==TRUE) {
      omega <- polr(as.factor(c(ttnsr[tensor>0]))~offset(-c(theta[tensor>0])))$zeta
    }
    (prev <- likelihood(ttnsr[tensor>0],theta[tensor>0],omega))
    
    
    #update A_1
    W =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
    A_1 <- comb(A_1,W,ttnsr,1,omega,alphbound)
    #orthognalize A_1
    qr_res=qr(A_1)
    A_1=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),1)
    theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    (new <- likelihood(ttnsr[tensor>0],theta[tensor>0],omega))
    
    if(max(abs(theta))>=alph) break
    
    # update A_2
    W <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
    A_2 <- comb(A_2,W,ttnsr,2,omega,alphbound)
    #orthognalize A_2
    qr_res=qr(A_2)
    A_2=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),2)
    theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    (new <- likelihood(ttnsr[tensor>0],theta[tensor>0],omega))
    if(max(abs(theta))>=alph) break
    
    
    # update A_3
    W <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
    A_3 <- comb(A_3,W,ttnsr,3,omega,alphbound)
    #orthognalize A_3
    qr_res=qr(A_3)
    A_3=qr.Q(qr_res)
    C=ttm(C,qr.R(qr_res),3)
    theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    (new <- likelihood(ttnsr[tensor>0],theta[tensor>0],omega))
    if(max(abs(theta))>=alph) break
    
    # update C
    W <- kronecker(kronecker(A_3,A_2),A_1)
    C <- corecomb(C,W,c(ttnsr),omega)
    theta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
    (new <- likelihood(ttnsr[tensor>0],theta[tensor>0],omega))
    cost = c(cost,new)
    error <- abs((new-prev)/prev)
    if(max(abs(theta))>=alph) break
    
  }
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$iteration <- iter
  result$cost = cost; result$omega=omega; result$lastcost = new
  return(result)
}





fit_ordinal_missing = function(tensor,C,A_1,A_2,A_3,omega = TRUE,alph = 50){
  #initialization
  ttnsr = tensor
  output = list()
  error<- 3
  iter = 0
  cost=NULL
  
  #initial update
  (result = fit_missing(tensor,ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph))
  pcost = result$lastcost
  theta = ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3)@data
  
  #update parameter
  C = result$C; A_1 = result$A_1; A_2 = result$A_2; A_3 = result$A_3
  
  #Data update
  ttnsr[which(tensor==-1)] = realization(theta,result$omega)@data[which(tensor==-1)]
  
  while ((error > 10^-4)&(iter<50)) {
    #update(maximization)
    (result = fit_missing(tensor,ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph))
    ncost = result$lastcost
    
    C = result$C; A_1 = result$A_1; A_2 = result$A_2; A_3 = result$A_3
    theta = ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3)@data
    
    #update data
    ttnsr[which(tensor==-1)] = realization(theta,result$omega)@data[which(tensor==-1)]
    (error = abs((ncost-pcost)/ncost))
    (cost = c(cost,ncost))
    pcost = ncost
    iter = iter+1
    
  }
  output$C = result$C; output$A_1 = result$A_1; output$A_2 = result$A_2; output$A_3 = result$A_3
  output$omega = result$omega
  output$iteration = iter; output$error= error; output$costvar = cost
  return(output)
}

### The first initial point
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)



### The second initial point
ip = c(length(which(tensor==1)),length(which(tensor==2)),length(which(tensor==3)), 
       length(which(tensor==4)), length(which(tensor==5)))/length(which(tensor>0))
ttnsr = tensor
ttnsr[which(tensor==-1)]=sample(1:5,length(which(tensor==-1)),replace = T,prob = ip)

initial = tucker(as.tensor(ttnsr),rank=c(4,4,4))
A_1 = initial$U[[1]]
A_2 = initial$U[[2]]
A_3 = initial$U[[3]]
C = initial$Z
