source("functions.R")
########## simulate tucker tensor with d, r ##############################
set.seed(18)
dlist=c(20,30,40,50,60)
r=c(3,3,3)
RMSE=matrix(0,nrow=5,ncol=5)
set.seed(18)

for(i in 1:5){
for(nsim in 1:5){
d=rep(dlist[i],3)
B_1 = randortho(d[1])[,1:r[1]]
B_2 = randortho(d[2])[,1:r[2]]
B_3 = randortho(d[3])[,1:r[3]]
D = rand_tensor(modes = r)
theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
alpha=10
theta = theta/max(abs(theta))*alpha ## specify signal = 10
omega = sort(rnorm(2,0,1))

ttnsr <- realization(theta,omega)@data
#missing = 1*(rand_tensor(modes = d)>0)@data ## uniform missingness with p=0.5
#ttnsr[missing==1]=NA


#initial point
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10 ## make sure theta is in the interior 
thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data

tic()
result <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alpha=10)
toc()

RMSE[i,nsim]=sqrt(mean((result$theta-theta)^2))
}
}

###

ggplot(data=figure1, aes(x=rescalednpq,y=sqrtmse))+geom_line(aes(color=d1d2d3))+
scale_shape_manual(values=seq(0,15))+
geom_point(aes(shape=d1d2d3))+
labs(x='Rescaled sample size N', y="Root mean squared error (RMSE)", 
color="Number of clusters", shape="Number of clusters")+
geom_point(data=rs.order4,aes(x=rescalednpqs,y=sqrtmse,shape="(4,4,4,4)"))+
geom_line(data=rs.order4,aes(x=rescalednpqs,y=sqrtmse,color="(4,4,4,4)"))



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


#############
########## simulate tucker tensor with d, r ##############################
set.seed(18)
Output=NULL
dlist=c(20,30)
rlist=c(3,10)

for(i in dlist){
    for(j in rlist){
d=rep(i,3); r=rep(j,3);
B_1 = randortho(d[1])[,1:r[1]]
B_2 = randortho(d[2])[,1:r[2]]
B_3 = randortho(d[3])[,1:r[3]]
D = rand_tensor(modes = r)
theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
alpha=10
theta = theta/max(abs(theta))*alpha ## specify signal = 10
omega = c(-0.3,-0.1,0.1,0.3) ## equal spaced
        
ttnsr <- realization(theta,omega)@data
                
#initial point
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10 ## make sure theta is in the interior 
thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
        
time=system.time(result <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alpha=2*alpha))
time=time[2]/length(result$iteration)

#likelihood(ttnsr,result$theta,omega)
#likelihood(ttnsr,theta,omega)
#likelihood(ttnsr,thetainit,omega)

Output=rbind(Output,c(time,d[1],r[1],-likelihood(ttnsr,theta,omega)/prod(d),-result$cost/(prod(d))))
    }
}

colnames(Output)=c("dim","rank","true_parameter",paste(rep("iter",10),1:10,sep=""))
write.table(Output,"Figure2.txt",quote=F)
##################
table=read.table("Figure2.txt",header=T)
cost=c(unlist(table[,4:13]))
iteration=rep(1:10,rep(4,10))
setting=rep(c("d=50,r=5 (1.03 sec per iteration)","d=50,r=10 (1.25 sec per iteration)","d=60,r=5 (2.65 sec per iteration)","d=60,r=10 (2.80 sec per iteration)"),10)
data=data.frame(cost=cost,setting=setting,iteration=iteration)
p=ggplot(data=data,aes(x=iteration,y=cost))+geom_line(aes(color=setting))+geom_point(aes(shape=setting))+ labs(x='iteration', y='log-likelihood')

g <- ggplot_build(p)
col=unique(g$data[[1]]$colour)[c(2,1,4,3)]
for(r in 1:4){
    p=p+geom_hline(yintercept=table[r,3],lty=2,col=col[r])
}
pdf("algorithm.pdf",width=8,height=3)
p
dev.off()



