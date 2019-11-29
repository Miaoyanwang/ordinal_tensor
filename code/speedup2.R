library(MASS)
library(rTensor)
library(pracma)
library(ggplot2)
library(ggthemes)
library(gridExtra)

realization = function(tnsr,alpha){
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



Hessi = function(A_1,W4,ttnsr,omega){
  thet =W4%*%c(A_1)
  p1 = logistic(thet + omega[1])
  p2 = logistic(thet + omega[2])
  q1 = p1*(1-p1)
  q2 = p2*(1-p2)+p1*(1-p1)
  q3 = p2*(1-p2)
  H = t(W4[which(c(ttnsr)==1),])%*%diag(q1[which(c(ttnsr)==1)])%*%W4[which(c(ttnsr)==1),]+
    t(W4[which(c(ttnsr)==2),])%*%diag(q2[which(c(ttnsr)==2)])%*%W4[which(c(ttnsr)==2),]+
    t(W4[which(c(ttnsr)==3),])%*%diag(q3[which(c(ttnsr)==3)])%*%W4[which(c(ttnsr)==3),]
  return(H)
}



h1 = function(A_1,W1,ttnsr,omega){
  thet =W1%*%c(A_1)
  p1 = logistic(thet + omega[1])
  p2 = logistic(thet + omega[2])
  p = cbind(p1,p2-p1,1-p2)
  return(-sum(log(c(p[which(c(ttnsr)==1),1],p[which(c(ttnsr)==2),2],p[which(c(ttnsr)==3),3]))))
}


g1 = function(A_1,W1,ttnsr,omega){
  thet =W1%*%c(A_1)
  p1 = logistic(thet + omega[1])
  p2 = logistic(thet + omega[2])
  q1 <- p1-1
  q2 <- (p2*(1-p2)-p1*(1-p1))/(p1-p2)
  q3 <- p2
  gd = apply(diag(q1[which(c(ttnsr)==1)])%*%W1[which(c(ttnsr)==1),],2,sum)+
    apply(diag(q2[which(c(ttnsr)==2)])%*%W1[which(c(ttnsr)==2),],2,sum)+
    apply(diag(q3[which(c(ttnsr)==3)])%*%W1[which(c(ttnsr)==3),],2,sum)
    return(gd)
}


comb = function(A,W,ttnsr,k,omega,alph=TRUE){
  nA = A
  tnsr1 <- k_unfold(as.tensor(ttnsr),k)@data
  if (alph==TRUE) {
    l <- lapply(1:nrow(A),function(i){optim(A[i,],
                                            function(x) h1(x,W,tnsr1[i,],omega),
                                            function(x) g1(x,W,tnsr1[i,],omega),
                                            method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }else{
    l <- lapply(1:nrow(A),function(i){constrOptim(A[i,],
                                                  function(x) h1(x,W,tnsr1[i,],omega),function(x) g1(x,W,tnsr1[i,],omega),
                                                  ui = rbind(W,-W),ci = rep(-alph,2*nrow(W)),method = "BFGS")$par})
    nA <- matrix(unlist(l),nrow = nrow(A),byrow = T)
  }
  return(nA)
}
tic()
optim(Cvec,h,g,method = "BFGS")
toc()
tic()
nlminb(Cvec,h,g,H)
toc()

corecomb = function(C,W,ttnsr,omega,alph=TRUE){
  Cvec <- c(C@data)
  h <- function(x) h1(x,W,ttnsr,omega)
  g <- function(x) g1(x,W,ttnsr,omega)
  H <- function(x) Hessi(x,W,ttnsr,omega)
  d <- nlminb(Cvec,h,g,H)
  C <- new("Tensor",C@num_modes,C@modes,data =d$par)

  return(C)
}

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




#####################################################################
fit_ordinal = function(ttnsr,C,A_1,A_2,A_3,omega,alph = TRUE){
  alphbound <- alph+10^-4
  result = list()
  error<- 3
  iter = 0
  d1 <- nrow(A_1); d2 <- nrow(A_2); d3 <- nrow(A_3)
  r1 <- ncol(A_1); r2 <- ncol(A_2); r3 <- ncol(A_3)
  if (alph == TRUE) {
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      #update A_1
      prevtheta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      prev <- likelihood(ttnsr,prevtheta,omega)
      W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W1,ttnsr,1,omega)
      
      
      # update A_2
      W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      A_2 <- comb(A_2,W2,ttnsr,2,omega)
      
      # update A_3
      W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      A_3 <- comb(A_3,W3,ttnsr,3,omega)
      
      # update C
      W4 <- kronecker(kronecker(A_3,A_2),A_1)
      C <- corecomb(C,W4,c(ttnsr),omega)
      theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      new <- likelihood(ttnsr,theta,omega)
      (error <- abs((new-prev)/prev))
    }
  }else{
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      #update A_1
      prevtheta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      prev <- likelihood(ttnsr,prevtheta,omega)
      W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W1,ttnsr,1,omega,alphbound)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      
      # update A_2
      W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      A_2 <- comb(A_2,W2,ttnsr,2,omega,alphbound)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      # update A_3
      W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      A_3 <- comb(A_3,W3,ttnsr,3,omega,alphbound)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      # update C
      W4 <- kronecker(kronecker(A_3,A_2),A_1)
      C <- corecomb(C,W4,c(ttnsr),omega)
      theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      new <- likelihood(ttnsr,theta,omega)
      error <- abs((new-prev)/prev)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
    }
  }
  
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$iteration <- iter
  return(result)
}


fit_ordinal2 = function(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph = TRUE){
  omega <- sort(rnorm(2))
  alphbound <- alph+10^-4
  result = list()
  error<- 3
  iter = 0
  d1 <- nrow(A_1); d2 <- nrow(A_2); d3 <- nrow(A_3)
  r1 <- ncol(A_1); r2 <- ncol(A_2); r3 <- ncol(A_3)
  if (alph == TRUE) {
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      #update A_1
      prevtheta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      prev <- likelihood(ttnsr,prevtheta,omega)
      W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W1,ttnsr,1,omega)
      
      
      # update A_2
      W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      A_2 <- comb(A_2,W2,ttnsr,2,omega)
      
      # update A_3
      W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      A_3 <- comb(A_3,W3,ttnsr,3,omega)
      
      # update C
      W4 <- kronecker(kronecker(A_3,A_2),A_1)
      C <- corecomb(C,W4,c(ttnsr),omega)

      #update omega
      theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      m <- polr(as.factor(c(ttnsr))~offset(-c(theta@data)))
      omega <- m$zeta
      
      
      
      theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      new <- likelihood(ttnsr,theta,omega)
      error <- abs((new-prev)/prev)
    }
  }else{
    while ((error > 10^-4)&(iter<50) ) {
      iter = iter +1
      
      #update A_1
      prevtheta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      prev <- likelihood(ttnsr,prevtheta,omega)
      W1 =kronecker(A_3,A_2)%*%t(k_unfold(C,1)@data)
      A_1 <- comb(A_1,W1,ttnsr,1,omega,alphbound)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      
      # update A_2
      W2 <- kronecker(A_3,A_1)%*%t(k_unfold(C,2)@data)
      A_2 <- comb(A_2,W2,ttnsr,2,omega,alphbound)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      # update A_3
      W3 <- kronecker(A_2,A_1)%*%t(k_unfold(C,3)@data)
      A_3 <- comb(A_3,W3,ttnsr,3,omega,alphbound)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      # update C
      W4 <- kronecker(kronecker(A_3,A_2),A_1)
      C <- corecomb(C,W4,c(ttnsr),omega)
      if(max(abs(ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)@data))>=alph) break
      
      #update omega
      theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      m <- polr(as.factor(c(ttnsr))~offset(-c(theta@data)))
      omega <- m$zeta
      
      
      theta <- ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)
      new <- likelihood(ttnsr,theta,omega)
      error <- abs((new-prev)/prev)
    }
  }
  
  result$C <- C; result$A_1 <- A_1; result$A_2 <- A_2; result$A_3 <- A_3
  result$iteration <- iter; result$omega <- omega
  return(result)
}




#construction when d=20
B_1 = matrix(runif(20*3,min=-1,max=1),nrow = 20)
B_2 = matrix(runif(20*3,min=-1,max=1),nrow = 20)
B_3 = matrix(runif(20*3,min=-1,max=1),nrow = 20)
D = as.tensor(array(runif(3^3,min=-1,max=1),dim = c(3,3,3)))
theta = ttm(ttm(ttm(D,B_1,1),B_2,2),B_3,3)*2
max(abs(theta@data))


omega = c(-0.2,0.2)

ttnsr <- realization(theta,omega)@data



#initial point
set.seed(18)
A_1 = randortho(20)[,1:3]
A_2 = randortho(20)[,1:3]
A_3 = randortho(20)[,1:3]
C = rand_tensor(modes = c(3,3,3))
thetainit = ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)



## with omega
tic()
result1 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega)
toc()

result2 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega,10)

thetahat1 <- ttm(ttm(ttm(result1$C,result1$A_1,1),result1$A_2,2),result1$A_3,3)
thetahat2 <- ttm(ttm(ttm(result2$C,result2$A_1,1),result2$A_2,2),result2$A_3,3)
likelihood(ttnsr,thetahat1,omega)
likelihood(ttnsr,thetahat2,omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)

md <- lm(c(thetahat1@data)~c(theta@data))
plot(c(theta@data),c(thetahat1@data),xlab = expression(theta),ylab = expression(hat(theta)),
     main = expression(paste("when d=20 with ",omega,"=(-0.2,0.2)")),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),
       col=c("red", "blue"), lty = c(1,1),cex=0.8)


d20womg <- as.data.frame(cbind(c(ttnsr),c(theta@data), c(thetahat1@data),c(thetainit@data)))
names(d20womg) <- c("data","theta","thetahat","thetainit") 
head(d20womg)
write.table(d20womg, file = "d20womg.csv")
head(read.table("d20womg.csv",head = T))
lm(theta~thetahat,a)

## without omega
tic()
result3 <- fit_ordinal2(ttnsr,C,A_1,A_2,A_3)
toc()
result4 <- fit_ordinal2(ttnsr,C,A_1,A_2,A_3,alph=10)


thetahat3 <- ttm(ttm(ttm(result3$C,result3$A_1,1),result3$A_2,2),result3$A_3,3)
thetahat4 <- ttm(ttm(ttm(result4$C,result4$A_1,1),result4$A_2,2),result4$A_3,3)
likelihood(ttnsr,thetahat3,omega)
likelihood(ttnsr,thetahat4,omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)


md <- lm(c(thetahat3@data)~c(theta@data))

plot(c(theta@data),c(thetahat3@data),xlab = expression(theta),ylab = expression(hat(theta)),
     main = expression(paste("when d=20 without ",omega)),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),
       col=c("red", "blue"), lty = c(1,1),cex=0.8)



d20woomg <- as.data.frame(cbind(c(ttnsr),c(theta@data), c(thetahat3@data),c(thetainit@data)))
names(d20woomg) <- c("data","theta","thetahat",'thetainit') 
head(d20woomg)
write.table(d20woomg, file = "d20woomg.csv")
head(read.table("d20womg.csv",head = T))


#construction when d=20
B_1 = matrix(runif(30*3,min=-1,max=1),nrow = 30)
B_2 = matrix(runif(30*3,min=-1,max=1),nrow = 30)
B_3 = matrix(runif(30*3,min=-1,max=1),nrow = 30)
D = as.tensor(array(runif(3^3,min=-1,max=1),dim = c(3,3,3)))
theta = ttm(ttm(ttm(D,B_1,1),B_2,2),B_3,3)*2
max(abs(theta@data))


omega = c(-0.2,0.2)

ttnsr <- realization(theta,omega)@data



#initial point
set.seed(18)
A_1 = randortho(30)[,1:3]
A_2 = randortho(30)[,1:3]
A_3 = randortho(30)[,1:3]
C = rand_tensor(modes = c(3,3,3))
thetainit = ttm(ttm(ttm(C,A_1,1),A_2,2),A_3,3)



## with omega
tic()
result1 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega)
toc()
result2 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega,10)

thetahat1 <- ttm(ttm(ttm(result1$C,result1$A_1,1),result1$A_2,2),result1$A_3,3)
thetahat2 <- ttm(ttm(ttm(result2$C,result2$A_1,1),result2$A_2,2),result2$A_3,3)
likelihood(ttnsr,thetahat1,omega)
likelihood(ttnsr,thetahat2,omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)

md <- lm(c(thetahat1@data)~c(theta@data))
plot(c(theta@data),c(thetahat1@data),xlab = expression(theta),ylab = expression(hat(theta)),
     main = expression(paste("when d=30 with ",omega,"=(-0.2,0.2)")),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),
       col=c("red", "blue"), lty = c(1,1),cex=0.8)


## without omega
tic()
result3 <- fit_ordinal2(ttnsr,C,A_1,A_2,A_3)
toc()
result4 <- fit_ordinal2(ttnsr,C,A_1,A_2,A_3,alph=10)


thetahat3 <- ttm(ttm(ttm(result3$C,result3$A_1,1),result3$A_2,2),result3$A_3,3)
thetahat4 <- ttm(ttm(ttm(result4$C,result4$A_1,1),result4$A_2,2),result4$A_3,3)
likelihood(ttnsr,thetahat3,omega)
likelihood(ttnsr,thetahat4,omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)


md <- lm(c(thetahat3@data)~c(theta@data))

plot(c(theta@data),c(thetahat3@data),xlab = expression(theta),ylab = expression(hat(theta)),
     main = expression(paste("when d=30 without ",omega)),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),
       col=c("red", "blue"), lty = c(1,1),cex=0.8)



###############################BIC###################################################
likelihood = function(ttnsr,thet,alpha){
  p1 = logistic(c(thet@data) + alpha[1])
  p2 = logistic(c(thet@data) + alpha[2])
  p = cbind(p1,p2-p1,1-p2)
  return(-sum(log(c(p[which(c(ttnsr)==1),1],p[which(c(ttnsr)==2),2],p[which(c(ttnsr)==3),3]))))
}

bic = function(ttnsr,theta,omega,d,r){
  return(2*likelihood(ttnsr,theta,omega)+(r^3+3*r*(d-r))*log(d^3))
}

#construction when d=20
B_1 = matrix(runif(20*3,min=-1,max=1),nrow = 20)
B_2 = matrix(runif(20*3,min=-1,max=1),nrow = 20)
B_3 = matrix(runif(20*3,min=-1,max=1),nrow = 20)
D = as.tensor(array(runif(3^3,min=-1,max=1),dim = c(3,3,3)))
theta = ttm(ttm(ttm(D,B_1,1),B_2,2),B_3,3)*2
max(abs(theta@data))

omega = c(-0.2,0.2)
ttnsr <- realization(theta,omega)@data


#initial point
#true rank =3 d=20
#true rank =6 d=20

set.seed(18)
bicr <- matrix(nrow = 5,ncol = 2)
colnames(bicr) <- c("rank","BIC")
for(i in 2:6){
  A_1 = randortho(20)[,1:i]
  A_2 = randortho(20)[,1:i]
  A_3 = randortho(20)[,1:i]
  C = rand_tensor(modes = c(i,i,i))
  result1 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega)
  thetahat1 <- ttm(ttm(ttm(result1$C,result1$A_1,1),result1$A_2,2),result1$A_3,3)
  bicr[i-1,1] <-  i
  bicr[i-1,2] <- bic(ttnsr,thetahat1,omega,20,i)
}

B_1 = matrix(runif(20*4,min=-1,max=1),nrow = 20)
B_2 = matrix(runif(20*4,min=-1,max=1),nrow = 20)
B_3 = matrix(runif(20*4,min=-1,max=1),nrow = 20)
D = as.tensor(array(runif(4^3,min=-1,max=1),dim = c(4,4,4)))
theta = ttm(ttm(ttm(D,B_1,1),B_2,2),B_3,3)*2
max(abs(theta@data))

omega = c(-0.2,0.2)
ttnsr <- realization(theta,omega)@data


#initial point
#true rank =3 d=20
#true rank =6 d=20

set.seed(18)
bicr <- matrix(nrow = 5,ncol = 2)
colnames(bicr) <- c("rank","BIC")
for(i in 2:6){
  A_1 = randortho(20)[,1:i]
  A_2 = randortho(20)[,1:i]
  A_3 = randortho(20)[,1:i]
  C = rand_tensor(modes = c(i,i,i))
  result1 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega)
  thetahat1 <- ttm(ttm(ttm(result1$C,result1$A_1,1),result1$A_2,2),result1$A_3,3)
  bicr[i-1,1] <-  i
  bicr[i-1,2] <- bic(ttnsr,thetahat1,omega,20,i)
}








