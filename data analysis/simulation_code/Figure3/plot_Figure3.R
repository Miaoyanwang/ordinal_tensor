### Figure 3: MSE vs. signal level (alpha), under various ranks##
library(ggplot2)
rlist=c(3,5,8)
alist=c(0.6,0.7,seq(from=0.8,to=2.5,by=0.3))
d=rep(20,3)
################ plot ############
load("Figure3.RData")
Output$rank=as.factor(Output$rank)
levels(Output$rank)=c("r=3","r=5","r=8")
alist=alist[alist>0.6]
Output=Output[Output[,3]>0.6,]

#p=ggplot(data=Output,aes(x=alpha,y=MSE))+geom_line(aes(color=rank),size=1.2)+ labs(x=expression("signal level (log10 scale)"), y='relative MSE',size=1.5)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

p=ggplot(data=Output,aes(x=alpha,y=MSE))+geom_line(aes(color=rank),size=1.5)+ labs(x=expression("signal level"~(log[10]~alpha)), y='relative MSE',size=1.5)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(1,2,3)]


p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.01,
position=position_dodge(0.05),color=col[rep(rep(1:length(rlist)),length(alist))],size=0.5)+geom_point(aes(shape=rank),size=2)+scale_shape_manual(values = c(0,2,16))

pdf("error_alpha.pdf",width=5,height=4)
p
dev.off()







