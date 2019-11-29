source("function.R")


########## simulate tucker tensor with d, r ##############################
set.seed(18)
d=c(20,20,20)
r=c(3,3,3)
B_1 = matrix(runif(d[1]*r[1],min=-1,max=1),nrow = d[1])
B_2 = matrix(runif(d[2]*r[2],min=-1,max=1),nrow = d[2])
B_3 = matrix(runif(d[3]*r[3],min=-1,max=1),nrow = d[3])
D = as.tensor(array(runif(prod(r),min=-1,max=1),dim = r))
theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
theta = theta/max(abs(theta))*7 ## specify signal = 7


omega = c(-0.2,0.2)

ttnsr <- realization(theta,omega)@data
missing = 1*(rand_tensor(modes = d)>0)@data ## uniform missingness with p=0.5
ttnsr =ttnsr *missing 


#initial point
set.seed(18)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data

## with omega
tic()
result1 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega)
toc()
tic()
result2 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega,10)
toc()

thetahat1 <- ttl(result1$C,list(result1$A_1,result1$A_2,result1$A_3),ms=1:3)@data
thetahat2 <- ttl(result2$C,list(result2$A_1,result2$A_2,result2$A_3),ms=1:3)@data

likelihood(ttnsr,thetahat1,omega)
likelihood(ttnsr,thetahat2,omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)

md <- lm(c(thetahat1)~c(theta))
plot(c(theta),c(thetahat1),xlab = expression(theta),ylab = expression(hat(theta)),main = expression(paste("when d=20 with ",omega,"=(-0.2,0.2)")),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),col=c("red", "blue"), lty = c(1,1),cex=0.8)


## without omega
tic()
result3 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE)
toc()
tic()
result4 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alph=10)
toc()

thetahat3 <- ttm(ttm(ttm(result3$C,result3$A_1,1),result3$A_2,2),result3$A_3,3)@data
thetahat4 <- ttm(ttm(ttm(result4$C,result4$A_1,1),result4$A_2,2),result4$A_3,3)@data
likelihood(ttnsr,thetahat3,result3$omega)
likelihood(ttnsr,thetahat4,result4$omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)


md <- lm(c(thetahat3)~c(theta))
plot(c(theta),c(thetahat3),xlab = expression(theta),ylab = expression(hat(theta)),main = expression(paste("when d=20 without ",omega)),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),col=c("red", "blue"), lty = c(1,1),cex=0.8)


########################### cp tensor with d, r ##############################
d=c(20,20,20); r=3
B_1 = matrix(runif(d[1]*r,min=-1,max=1),nrow = d[1])
B_2 = matrix(runif(d[2]*r,min=-1,max=1),nrow = d[2])
B_3 = matrix(runif(d[3]*r,min=-1,max=1),nrow = d[3])
theta = tensorize(B_1,B_2,B_3)
theta = theta/max(abs(theta))*7 ## specify signal = 7



omega = c(-0.2,0.2)
ttnsr <- realization(theta,omega)@data


#initial point
set.seed(18)
A_1 = as.matrix(randortho(d[1])[,1:r])
A_2 = as.matrix(randortho(d[1])[,1:r])
A_3 = as.matrix(randortho(d[1])[,1:r])
thetainit = tensorize(A_1,A_2,A_3)


## with omega
tic()
result1 <- fit_ordinal_cp(ttnsr,A_1,A_2,A_3,omega)
toc()
tic()
result2 <- fit_ordinal_cp(ttnsr,A_1,A_2,A_3,omega,10)
toc()
thetahat1 <- tensorize(result1$A_1,result1$A_2,result1$A_3)
thetahat2 <- tensorize(result2$A_1,result2$A_2,result2$A_3)
likelihood(ttnsr,thetahat1,omega) ## unconstrained est
likelihood(ttnsr,thetahat2,omega) ## constrained est
likelihood(ttnsr,theta,omega) ## truth
likelihood(ttnsr,thetainit,omega) ## initilization

md <- lm(c(thetahat1)~c(theta))
plot(c(theta),c(thetahat1),xlab = expression(theta),ylab=expression(hat(theta)),main=expression(paste("when d=20 with ",omega,"=(-0.2,0.2)")),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),col=c("red", "blue"), lty = c(1,1),cex=0.8)



## without omega
tic()
result3 <- fit_ordinal_cp(ttnsr,A_1,A_2,A_3,omega=TRUE)
toc()
result4 <- fit_ordinal_cp(ttnsr,A_1,A_2,A_3,alph=10,omega=TRUE)


thetahat3 <- tensorize(result3$A_1,result3$A_2,result3$A_3)
thetahat4 <- tensorize(result4$A_1,result4$A_2,result4$A_3)
likelihood(ttnsr,thetahat3,result3$omega)
likelihood(ttnsr,thetahat4,result4$omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)


md <- lm(c(thetahat3)~c(theta))
plot(c(theta),c(thetahat3),xlab = expression(theta),ylab=expression(hat(theta)),main=expression(paste("when d=20 without",omega)),xlim=c(-10,10),ylim=c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),col=c("red", "blue"), lty = c(1,1),cex=0.8)



######################## previous code ################
################### construction when d=30 ##############################
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



########################## BIC: construction when d=20, r=3 ##########################
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


########################## BIC: construction when d=20, r=4 ##########################
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








