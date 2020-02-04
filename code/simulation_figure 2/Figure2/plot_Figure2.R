### Figure: MSE vs. dimension, under various ranks##
library(ggplot2)
dlist=seq(from=15,to=40,by=2)
rlist=c(3,5,8)
################ plot ############
load("Figure2.RData")
Output$rank=as.factor(Output$rank)
levels(Output$rank)=c("r=3","r=5","r=8")

p=ggplot(data=Output,aes(x=dim,y=MSE))+geom_line(aes(color=rank),size=1.5)+ labs(x=expression("tensor dimension"~(d)), y='MSE',size=1.5)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

p=p+stat_function(fun=function(d) 10^2.96/(d^2),lty=2,size=0.9)
g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(1,2,3)]

p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.3,
position=position_dodge(0.05),color=col[rep(rep(1:3),length(dlist))],size=0.5)+geom_point(aes(shape=rank),size=1.5)

pdf("error_dimension.pdf",width=5,height=4)
p
dev.off()




