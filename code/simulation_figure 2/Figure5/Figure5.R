### Figure 5: MSE vs. number of ordinal levels; half observations ##
source("../../functions.R")
set.seed(1)
Output=NULL
d=rep(20,3)
rlist=c(3,5,8)
rho=0.5
Klist=2:7
ntotal=30
MSE_raw=array(0,dim=c(length(rlist),length(Klist),ntotal))
for(j in 1:length(rlist)){      
        r=rep(rlist[j],3)
        B_1 = randortho(d[1])[,1:r[1]]
        B_2 = randortho(d[2])[,1:r[2]]
        B_3 = randortho(d[3])[,1:r[3]]
        D = rand_tensor(modes = r)
        theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
        alpha=10 #
        theta = theta/max(abs(theta))*alpha ## specify signal = 10
        
        for(i in 1:length(Klist)){
        K=Klist[i]
        omega = logit((1:(K-1))/K) ## equal space on the logit scale
        
        
        for(nsim in 1:ntotal){ 
        ttnsr <- realization(theta,omega)@data
        data=array(NA,dim=dim(ttnsr))
        Obs=sample(1:prod(d),prod(d)*rho,replace=F)
        data[Obs]=ttnsr[Obs]
        #initial point
        A_1 = randortho(d[1])[,1:r[1]]
        A_2 = randortho(d[2])[,1:r[2]]
        A_3 = randortho(d[3])[,1:r[3]]
        C = rand_tensor(modes = r)
        C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10 ## initial theta is in the interior 
        thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
        
        result <- fit_ordinal(data,C,A_1,A_2,A_3,omega=TRUE,alpha=2*alpha) ## a slightly relaxed boundary
        
        MSE_raw[j,i,nsim]=mean((theta-result$theta)^2)/mean(theta^2) ## more stable using rescaled MSE
            }
    }
}

Output=data.frame(MSE=c(t(apply(MSE_raw,c(1:2),mean))),rank=as.factor(rep(rlist,rep(length(Klist),3))),K=rep(Klist,length(rlist)),sd=c(t(apply(MSE_raw,c(1:2),sd))))
save(Output,MSE_raw,file="Figure5.RData")

