### Figure: MSE vs. dimension, under various ranks##
source("../../functions.R")
set.seed(1)
Output=NULL
dlist=seq(from=15,to=40,by=2)
rlist=c(3,5,8)
ntotal=30
MSE_raw=array(0,dim=c(length(rlist),length(dlist),ntotal))

for(j in 1:length(rlist)){
    for(i in 1:length(dlist)){
        for(nsim in 1:ntotal){
        d=rep(dlist[i],3); r=rep(rlist[j],3);
        B_1 = randortho(d[1])[,1:r[1]]
        B_2 = randortho(d[2])[,1:r[2]]
        B_3 = randortho(d[3])[,1:r[3]]
        D = rand_tensor(modes = r)
        theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
        alpha=10
        theta = theta/max(abs(theta))*alpha ## specify signal = 10
        omega = logit((1:4)/5) ## equal spaced on the probit scale
        
        ttnsr <- realization(theta,omega)@data
        
        #initial point
        A_1 = randortho(d[1])[,1:r[1]]
        A_2 = randortho(d[2])[,1:r[2]]
        A_3 = randortho(d[3])[,1:r[3]]
        C = rand_tensor(modes = r)
        C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10 ## initial theta is in the interior 
        thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
        
        result <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE,alpha=2*alpha) ## relaxed bound
        
        MSE_raw[j,i,nsim]=mean((theta-result$theta)^2)
    }
}
}

Output=data.frame(MSE=c(t(apply(MSE_raw,c(1:2),mean))),rank=as.factor(rep(rlist,rep(length(dlist),length(rlist)))),dim=rep(dlist,length(rlist)),sd=c(t(apply(MSE_raw,c(1:2),sd))))
save(Output,MSE_raw,file="Figure2.RData")
