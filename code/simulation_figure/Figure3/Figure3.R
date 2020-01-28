### Figure 3: MSE vs. signal level (alpha), under various ranks##
source("../../functions.R")
set.seed(1)
Output=NULL
rlist=c(3,5,8)
##alist=seq(from=0.7,to=3.5,by=0.5)
##alist=seq(from=0.5,to=2.5,by=0.3)
alist=c(0.6,0.7,seq(from=0.8,to=2.5,by=0.3))
d=rep(20,3)
ntotal=30
MSE_raw=array(0,dim=c(length(rlist),length(alist),ntotal))

for(j in 1:length(rlist)){
    for(i in 1:length(alist)){
        for(nsim in 1:ntotal){
        alpha=10^(alist[i]); r=rep(rlist[j],3);
        B_1 = randortho(d[1])[,1:r[1]]
        B_2 = randortho(d[2])[,1:r[2]]
        B_3 = randortho(d[3])[,1:r[3]]
        D = rand_tensor(modes = r)
        theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
        theta = theta/max(abs(theta))*alpha ## specify signal = alpha
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
        #MSE=c(MSE,mean((theta-result$theta)^2)/mean(theta^2))
        MSE_raw[j,i,nsim]=mean((theta-result$theta)^2)/mean(theta^2)
        ##sprintf("%d\n%d\n%d",i,j,nsim)
    }
}
}

Output=data.frame(MSE=c(t(apply(MSE_raw,c(1:2),mean))),rank=as.factor(rep(rlist,rep(length(alist),length(rlist)))),alpha=rep(alist,length(rlist)),sd=c(t(apply(MSE_raw,c(1:2),sd))))
save(Output,MSE_raw,file="Figure3.RData")

################ plot ############
library(ggplot2)
load("Figure3.RData")
Output$rank=as.factor(Output$rank)
levels(Output$rank)=c("r=3","r=5","r=8")

p=ggplot(data=Output,aes(x=alpha,y=MSE))+geom_line(aes(color=rank),size=1.2)+ labs(x=expression("signal level (log10 scale)"), y='relative MSE',size=1.5)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

##p=ggplot(data=Output,aes(x=alpha,y=MSE))+geom_line(aes(color=rank),size=1.5)+ labs(x=expression("signal level"~(log[10]~alpha)), y='relative MSE',size=1.5)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(1,2,3)]


p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.01,
position=position_dodge(0.05),color=col[rep(rep(1:length(rlist)),length(alist))],size=1)+geom_point(aes(shape=rank),size=1.2)

pdf("error_alpha.pdf",width=5,height=4)
p
dev.off()







