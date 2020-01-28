### Figure: cost function vs. iteration ##
source("../../functions.R")
library(ggplot2)
set.seed(1)
Output=NULL
dlist=c(25,30)
rlist=c(5,10)

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
        C=C/max(abs(ttl(C,list(A_1,A_2,A_3),ms=1:3)@data))*alpha/10 ## initial theta is in the interior 
        thetainit=ttl(C,list(A_1,A_2,A_3),ms=1:3)@data
        
        time=system.time(result <- fit_ordinal(ttnsr,C,A_1,A_2,A_3,omega=TRUE))
        time=time[2]/length(result$iteration)
        
        Output=rbind(Output,c(time,d[1],r[1],-likelihood(ttnsr,theta,omega),-result$cost))
    }
}

colnames(Output)=c("time","dim","rank","true_parameter",paste(rep("iter",15),1:15,sep=""))
write.table(Output,"Figure2.txt",quote=F)


################ plot ############
Output=read.table("Figure2.txt",header=T)
cost=c(unlist(table[,5:19]))
iteration=rep(1:15,rep(4,15))
Setting=rep(c(sprintf("d=25,r=5 (%.2f sec/iter)",Output[1,1]),sprintf("d=25,r=10 (%.2f sec/iter)",Output[2,1]),sprintf("d=30,r=5 (%.2f sec/iter)",Output[3,1]),sprintf("d=30,r=10 (%.2f sec/iter)",Output[4,1])),15)
data=data.frame(cost=cost,Setting=Setting,iteration=iteration)
####
p=ggplot(data=data,aes(x=iteration,y=cost/(10^4)))+geom_line(aes(color=Setting))+geom_point(aes(shape=Setting))+ labs(x='iteration', y=expression("log-likelihood"%*%10^4))

g <- ggplot_build(p)
col=unique(g$data[[1]]$colour)[c(2,1,4,3)] ## match the line color to curve
for(r in 1:4){
    p=p+geom_hline(yintercept=table[r,4]/(10^4),lty=2,col=col[r])
}
pdf("algorithm.pdf",width=5,height=2.8)
p
dev.off()
