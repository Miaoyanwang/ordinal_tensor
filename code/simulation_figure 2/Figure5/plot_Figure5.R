### Figure 5: MSE vs. number of ordinal levels; half observations ##
library(ggplot2)
d=rep(20,3)
rlist=c(3,5,8)
Klist=2:7
################ plot ############
load("Figure5.RData")
Output$rank=as.factor(Output$rank)
levels(Output$rank)=c("r=3","r=5","r=8")

p=ggplot(data=Output,aes(x=K,y=MSE))+geom_line(aes(color=rank),size=1.2)+ labs(x=expression('Number of Ordinal Levels'~(L)), y='relative MSE')+xlim(2,7)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+ylim(0,2)

g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(1,2,3)]

p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.1,
position=position_dodge(0),color=col[rep(rep(1:3),length(Klist))],size=0.5)+geom_point(aes(shape=rank),size=1.5)

pdf("error_level.pdf",width=5,height=4)
p
dev.off()


