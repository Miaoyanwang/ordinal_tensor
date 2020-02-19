library(MASS)
library(rTensor)
library(pracma)


# initial point
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
prevtheta <- ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
theta=prevtheta-2*gradient_tensor(A_1,A_2,A_3,C,ttnsr,omega)

# constant term for the missing values
W = array(0, dim(tensor))
W[tensor>0] = 1

costf = function(C,A_1,A_2,A_3,theta,W){
  return(sum((W*(ttl(C,list(A_1,A_2,A_3),ms = 1:3)@data-theta))^2)/2)
}

cgm = function(C,A_1,A_2,A_3,theta,W,i){
  val = as.tensor(W*(ttl(C,list(A_1,A_2,A_3),ms = 1:3)@data-theta))
  if(i == 1){
    g = k_unfold(val,1)@data%*%t(k_unfold(ttl(C,list(A_2,A_3),ms = c(2,3)),1)@data)
  }else if(i == 2){
    g = k_unfold(val,2)@data%*%t(k_unfold(ttl(C,list(A_1,A_3),ms = c(1,3)),2)@data)
  }else if(i == 3){
    g = k_unfold(val,3)@data%*%t(k_unfold(ttl(C,list(A_1,A_2),ms = c(1,2)),3)@data)
  }else{
    g = ttl(val,list(t(A_1),t(A_2),t(A_3)),ms = 1:3)
  }
  return(g)
}

cgma = function(A,G,theta,W,i){
  g = (k_unfold(as.tensor(W),i)@data*(A%*%G-k_unfold(as.tensor(theta),i)@data))%*%t(G)
  return(g)
}


tucker_missing = function(A_1,A_2,A_3,C,theta,W,alph = TRUE){
  alphbound <- alph+10^-4
  result = list()
  error<- 3
  iter = 0
  cost=NULL
  if (alph == TRUE) {
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      (prev = costf(C,A_1,A_2,A_3,theta,W))
    
      # update A_1
      # tic()
      # l <- lapply(1:nrow(A_1),
      #             function(i){optim(A_1[i,],
      #                               function(x) costf(C,matrix(x,ncol =length(x)),A_2,A_3,theta[i,,,drop = F],W[i,,,drop= F]),
      #                               function(x) cgm(C,matrix(x,ncol = length(x)),A_2,A_3,theta[i,,,drop = F],W[i,,,drop = F],1),
      #                               method = "BFGS")$par})
      # #A_1 = matrix(unlist(l),nrow = nrow(A_1),byrow = T) #It is slower 
      # toc()
      # 
      # tic()
      # G = k_unfold(C,1)@data%*%t(kronecker(A_3,A_2))
      # f = function(x) return(costf(C,matrix(x,nrow = dim(A_1)[1]),A_2,A_3,theta,W))
      # g = function(x) return(cgma(matrix(x,nrow = dim(A_1)[1]),G,theta,W,1))
      # A_1 <- matrix(optim(c(A_1),f,g,method ="BFGS" )$par,nrow = dim(A_1)[1])
      
      # l <- lapply(1:nrow(A_1),
      #             function(i){optim(A_1[i,],
      #                               function(x) costf(C,matrix(x,ncol =length(x)),A_2,A_3,theta[i,,,drop = F],W[i,,,drop= F]),
      #                               function(x) cgma(x,G,theta[i,,,drop = F],W[i,,,drop = F],1),
      #                               method = "BFGS")$par}) # It is slower
      # toc()
      # 
      f = function(x) return(costf(C,matrix(x,nrow = dim(A_1)[1]),A_2,A_3,theta,W))
      g = function(x) return(cgm(C,matrix(x,nrow = dim(A_1)[1]),A_2,A_3,theta,W,1))
      A_1 <- matrix(optim(c(A_1),f,g,method ="BFGS" )$par,nrow = dim(A_1)[1])
      #orthognalize A_1
      qr_res=qr(A_1)
      A_1=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),1)
      
      # update A_2
      
      f = function(x) return(costf(C,A_1,matrix(x,nrow = dim(A_2)[1]),A_3,theta,W))
      g = function(x) return(cgm(C,A_1,matrix(x,nrow = dim(A_2)[1]),A_3,theta,W,2))
      A_2 <- matrix(optim(c(A_2),f,g,method = "BFGS")$par,nrow = dim(A_2)[1])
      
      #orthognalize A_2
      qr_res=qr(A_2)
      A_2=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),2)
      
      # update A_3
      f = function(x) return(costf(C,A_1,A_2,matrix(x,nrow = dim(A_3)[1]),theta,W))
      g = function(x) return(cgm(C,A_1,A_2,matrix(x,nrow = dim(A_3)[1]),theta,W,3))
      A_3 <- matrix(optim(c(A_3),f,g,method = "BFGS")$par,nrow = dim(A_3)[1])
      #orthognalize A_3
      qr_res=qr(A_3)
      A_3=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),3)
      
      # update C
      f = function(x) return(costf(new("Tensor",C@num_modes,C@modes,x),A_1,A_2,A_3,theta,W))
      g = function(x) return(c(cgm(new("Tensor",C@num_modes,C@modes,x),A_1,A_2,A_3,theta,W,4)@data))
      C <- new("Tensor",C@num_modes,C@modes,data =optim(c(C@data),f,g,method = "BFGS")$par)
      (new = costf(C,A_1,A_2,A_3,theta,W))
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
  }else{
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      (prev = costf(C,A_1,A_2,A_3,theta,W))
      #update A_1
      G = k_unfold(C,1)@data%*%t(kronecker(A_3,A_2))
      l <- lapply(1:nrow(A_1),
                  function(i){constrOptim(A_1[i,],
                                          function(x) costf(C,matrix(x,ncol =length(x)),A_2,A_3,theta[i,,,drop = F],W[i,,,drop= F]),
                                          function(x) cgma(x,G,theta[i,,,drop = F],W[i,,,drop = F],1),
                                          ui = as.matrix(rbind(t(G),-t(G))),ci = rep(-alph,2*nrow(t(G))),
                                          method = "BFGS")$par})
      A_1 = matrix(unlist(l),nrow = nrow(A_1),byrow = T)
      
      #orthognalize A_1
      qr_res=qr(A_1)
      A_1=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),1)
      if(max(round(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)),digits = 3)>=alph) break
      
      # update A_2
      G = k_unfold(C,2)@data%*%t(kronecker(A_3,A_1))
      l <- lapply(1:nrow(A_2),
                  function(i){constrOptim(A_2[i,],
                                          function(x) costf(C,A_1,matrix(x,ncol =length(x)),A_3,theta[,i,,drop = F],W[,i,,drop= F]),
                                          function(x) cgma(x,G,theta[,i,,drop = F],W[,i,,drop = F],2),
                                          ui = as.matrix(rbind(t(G),-t(G))),ci = rep(-alph,2*nrow(t(G))),
                                          method = "BFGS")$par})
      A_2 = matrix(unlist(l),nrow = nrow(A_2),byrow = T)
      #orthognalize A_2
      qr_res=qr(A_2)
      A_2=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),2)
      if(max(round(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)),digits = 3)>=alph) break
      
      # update A_3
      G = k_unfold(C,3)@data%*%t(kronecker(A_2,A_1))
      l <- lapply(1:nrow(A_3),
                  function(i){constrOptim(A_3[i,],
                                          function(x) costf(C,A_1,A_2,matrix(x,ncol =length(x)),theta[,,i,drop = F],W[,,i,drop= F]),
                                          function(x) cgma(x,G,theta[,,i,drop = F],W[,,i,drop = F],3),
                                          ui = as.matrix(rbind(t(G),-t(G))),ci = rep(-alph,2*nrow(t(G))),
                                          method = "BFGS")$par})
      A_3 = matrix(unlist(l),nrow = nrow(A_3),byrow = T)
      #orthognalize A_3
      qr_res=qr(A_3)
      A_3=qr.Q(qr_res)
      C=ttm(C,qr.R(qr_res),3)
      if(max(round(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data)),digits = 3)>=alph) break
      
      # update C
      f = function(x) return(costf(new("Tensor",C@num_modes,C@modes,x),A_1,A_2,A_3,theta,W))
      g = function(x) return(c(cgm(new("Tensor",C@num_modes,C@modes,x),A_1,A_2,A_3,theta,W,4)@data))
      C <- new("Tensor",C@num_modes,C@modes,data =optim(c(C@data),f,g,method = "BFGS")$par)
      (new = costf(C,A_1,A_2,A_3,theta,W))
      cost = c(cost,new)
      (error <- abs((new-prev)/prev))
    }
  }
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$iteration <- iter
  result$cost = cost
  return(result)
}
