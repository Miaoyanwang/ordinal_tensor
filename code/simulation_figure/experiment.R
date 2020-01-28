source("functions.R")
########## brain tensor analysis ##########
load("../data/dti_brain.RData")
##ttnsr=tensor[,,attr[,2]=="22-25"] ## subset of individuals
ttnsr=tensor
d=dim(ttnsr)

r=c(25,25,7)
### random initilization
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
########## Music analysis #######################
source("functions.R")###
load("../data/InCar_Music.RData")
ttnsr=tensor
ttnsr[ttnsr==-1]=NA
missing=which(is.na(ttnsr)==F)
d=dim(ttnsr)
set.seed(1)
r=c(1,1,1) ## need to select via BIC ..
alpha=100
A_1 = as.matrix(randortho(d[1])[,1:r[1]])
A_2 = as.matrix(randortho(d[2])[,1:r[2]])
A_3 = as.matrix(randortho(d[3])[,1:r[3]])
C = rand_tensor(modes = r)
C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10


############# experiement ##########################
mask=sample(1:2884,2884*0.2,replace=T)
index=which(is.na(ttnsr)==F)[mask]
true=ttnsr[index]
ttnsr[index]=NA
##### Method 1: continuous tucker #########
tic()
result = fit_continuous(ttnsr,C,A_1,A_2,A_3,alpha=alpha)
toc()
thetahat=result$theta
mean((true-thetahat[index])^2) #MSE= 18.96803
mean(abs(true-thetahat[index])) #MAD= 2.589925
mean(true!=round(thetahat[index])) ## MCR

###### Method 2: ordinal tucker #########
tic()
result2 = fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alpha =100)
toc()
thetahat2=result2$theta
bic(ttnsr,thetahat2,result2$omega,d,r)## 10860.46
est=estimation(thetahat2,result2$omega,"mean")
mean((true-est@data[index])^2) ## MSE = 2.381944
mean(abs(true-est@data[index])) ## MAD = 1.256944
mean(true!=round(est@data[index])) ## MCR

######### Method 3: ordinal matrix##################
data=k_unfold(as.tensor(ttnsr),1)@data
d=dim(data)
r=1
set.seed(18)
A_1=  as.matrix(randortho(d[1])[,1:r])
A_2 = matrix(rnorm(d[2]*r,0,1),nrow=d[2])
A_2 = A_2/max(abs(A_1%*%t(A_2)))*alpha/10

tic()
result3 = fit_ordinal_matrix(data,A_1,A_2,omega=TRUE,alpha = alpha)
toc()
thetahat3=result3$theta
est=estimation(thetahat3,c(0,0),"mean") ## warning: omega cannot be estimated!!
mean((true-est[index])^2) ## MSE = 4.173611
mean(abs(true-est[index])) ## MAD = 1.395833
mean(true!=round(est[index])) ## MCR
####### end of data analysis ###

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
ttnsr[missing==1]=NA


#initial point
set.seed(18)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
alpha=7
C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10 ## make sure theta is in the interior 
thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data



## with omega
tic()
result1 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE)
toc()

tic()
result2 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega,alpha=10)
toc()

thetahat1 <- result1$theta
thetahat2 <- result2$theta

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
result4 <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE,20)
toc()

thetahat3 <- result3$theta
thetahat4 <- result4$theta
likelihood(ttnsr,thetahat3,result3$omega)
likelihood(ttnsr,thetahat4,result4$omega)
likelihood(ttnsr,theta,omega)
likelihood(ttnsr,thetainit,omega)


md <- lm(c(thetahat3)~c(theta))
plot(c(theta),c(thetahat3),xlab = expression(theta),ylab = expression(hat(theta)),main = expression(paste("when d=20 without ",omega)),xlim = c(-10,10),ylim = c(-10,10))
abline(md,col = 'red')
abline(0,1,col = 'blue')
legend("topleft", legend=c(paste("slope = ",round(md$coefficients[2],3)), "slope = 1"),col=c("red", "blue"), lty = c(1,1),cex=0.8)





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







