### Figure: MSE vs. sample size (|Omega|), under various rank##
library(ggplot2)
d=rep(20,3)
rlist=c(3,5,8)
rholist=seq(from=0.3,to=1,by=0.1)
load("Figure4.RData")

Output$rank=as.factor(Output$rank)
levels(Output$rank)=c("r=3","r=5","r=8")

p=ggplot(data=Output,aes(x=rho,y=MSE))+geom_line(aes(color=rank),size=1.5)+geom_point(aes(shape=rank))+ labs(x=expression('observation fraction'~(abs(Omega)/d^K)), y='relative MSE',size=1.2)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

#p=p+stat_function(fun=function(d)10^5/(d^2),lty=2)
g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(1,2,3)]

p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.01,
position=position_dodge(0.05),color=col[rep(rep(1:length(rlist)),length(rholist))],size=0.5)+geom_point(aes(shape=rank),size=2)+scale_shape_manual(values = c(0,2,16))

pdf("error_sample.pdf",width=5,height=4)
p
dev.off()








