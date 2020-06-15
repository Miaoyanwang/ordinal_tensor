### Figure: cost function vs. iteration ##
library(ggplot2)
dlist=c(25,30)
rlist=c(5,10)

################ plot ############
Output=read.table("Figure1.txt",header=T)
cost=c(unlist(Output[,5:19]))
iteration=rep(1:15,rep(4,15))
Setting=rep(c(sprintf("d=25,r=5 (%.2f sec/iter)",Output[1,1]),sprintf("d=25,r=10 (%.2f sec/iter)",Output[2,1]),sprintf("d=30,r=5 (%.2f sec/iter)",Output[3,1]),sprintf("d=30,r=10 (%.2f sec/iter)",Output[4,1])),15)
data=data.frame(cost=cost,Setting=Setting,iteration=iteration)
####
p=ggplot(data=data,aes(x=iteration,y=cost/(10^4)))+geom_line(aes(color=Setting),size=1.2)+geom_point(aes(shape=Setting))+ labs(x='iteration', y=expression("log-likelihood"%*%10^4))

g <- ggplot_build(p)
col=unique(g$data[[1]]$colour)[c(2,1,4,3)] ## match the line color to curve
for(r in 1:4){
    p=p+geom_hline(yintercept=Output[r,4]/(10^4),lty=2,col=col[r])
}
pdf("algorithm.pdf",width=6,height=2.8)
p
dev.off()
