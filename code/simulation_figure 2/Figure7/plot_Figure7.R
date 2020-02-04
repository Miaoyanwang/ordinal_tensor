### Figure: MSE vs. number of levels (L), under different methods ##
source("../../functions.R")
library(ggplot2)
library("R.matlab")
d=rep(20,3)
r=rep(5,3)
rho=0.8
Klist=2:7
load("Figure8.RData")

Err1_mean=apply(Err1,c(1,2),mean)
Err1_sd=apply(Err1,c(1,2),sd)
Err2_mean=apply(Err2,c(1,2),mean)
Err2_sd=apply(Err2,c(1,2),sd)
Err3_mean=apply(Err3,c(1,2),mean)
Err3_sd=apply(Err3,c(1,2),sd)
Err4_mean=apply(Err4,c(1,2),mean)
Err4_sd=apply(Err4,c(1,2),sd)
Err5_mean=apply(Err5,c(1,2),mean)
Err5_sd=apply(Err5,c(1,2),sd)

index=3
MAD_for_mode=c(Err1_mean[,index],Err2_mean[,index],Err3_mean[,index],Err4_mean[,index],Err5_mean[,index])
MAD_for_mode_sd=c(Err1_sd[,index],Err2_sd[,index],Err3_sd[,index],Err4_sd[,index],Err5_sd[,index])
method=c("Ordinal-T (ours)","Continuous-T","Ordinal-M","1bit-caregory-T","1bit-sign-T")
method=rep(method,rep(length(Klist),5))
K=rep(Klist,5)
Output=data.frame(MAD=MAD_for_mode,sd=MAD_for_mode_sd,method=method,K=K)

################ plot ############
p=ggplot(data=Output,aes(x=K,y=MAD))+geom_line(aes(color=method),size=1.2)+geom_point(aes(shape=method))+ labs(x=expression('Number of Ordinal Levels'~(L)), y='MAD')+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(5,3,4,1,2)]

p=p+geom_errorbar(aes(ymin=MAD-sd, ymax=MAD+sd), width=0.01,
position=position_dodge(0.05),color=col[rep(rep(1:5),length(Klist))],size=0.5)+geom_point(aes(shape=method),size=1.5)

pdf("MAD_median.pdf",width=6.2,height=4)
p
dev.off()








