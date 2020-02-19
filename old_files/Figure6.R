### Figure: MSE vs. sample size (|Omega|), under various rank##
source("../../functions.R")
library("R.matlab")
set.seed(1)
Output=NULL
d=rep(20,3)
r=rep(5,3) ## r cannot 1, otherwise worse than other methods...
ntotal=30
rholist=seq(from=0.5,to=1,by=0.05)
Err=Err2=Err3=Err4=Err5=Err6=matrix(0,nrow=length(rholist),ncol=4)

#### start
B_1 = as.matrix(randortho(d[1])[,1:r[1]])
B_2 = as.matrix(randortho(d[2])[,1:r[2]])
B_3 = as.matrix(randortho(d[3])[,1:r[3]])
D = rand_tensor(modes = r)
theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
alpha=10
theta = theta/max(abs(theta))*alpha ## specify signal = 10
omega = logit((1:4)/5) ## equal spaced on the probit scale
#omega=0   
K=length(omega)+1
##omega=logit((1:2)/3)
ttnsr <- realization(theta,omega)@data
 
######## start############################
for(i in 1:length(rholist)){
rho=rholist[i]
## generate observation for fit_ordinal
data=array(NA,dim=dim(ttnsr))
Obs=sample(1:prod(d),prod(d)*rho,replace=F)
data[Obs]=ttnsr[Obs] 
save(data,theta,omega,file=sprintf("../Figure6/alpha_R/input_level%.2f.RData",rholist[i]))
#### generate observation for 1-bit tensor
input=M_to_one(data,K,"category")
writeMat(sprintf("../Figure6/alpha_matlab/input_level%.2f_categorical.mat",rholist[i]),data=input$binary,E=input$E)
input=M_to_one(data,K,"cumulative")
writeMat(sprintf("../Figure6/alpha_matlab/input_level%.2f_cumulative.mat",rholist[i]),data=input$binary,E=input$E)
input=M_to_one(data,K,"sign")
writeMat(sprintf("../Figure6/alpha_matlab/input_level%.2f_sign.mat",rholist[i]),data=input$binary,E=input$E,ave=input$ave,scale=input$scale)
}

######################## decomposition ################################
for(i in 1:length(rholist)){
    load(sprintf("../Figure6/alpha_R/input_level%.2f.RData",rholist[i]))
    ######### Method 1: ordinal tensor ##################
    d=dim(data)
    A_1 = as.matrix(randortho(d[1])[,1:r[1]])
    A_2 = as.matrix(randortho(d[2])[,1:r[2]])
    A_3 = as.matrix(randortho(d[3])[,1:r[3]])
    C = rand_tensor(modes = r)
    result <- fit_ordinal(data,C,A_1,A_2,A_3,omega=TRUE,alpha=TRUE)
    save(result,file=sprintf("../Figure6/alpha_R/output_level%.2f_ordinal.RData",rholist[i]))
    ######### Method 2: continous tensor##################
    result2 <- fit_continuous(data,C,A_1,A_2,A_3,alpha=TRUE)
    save(result2,file=sprintf("../Figure6/alpha_R/output_level%.2f_cont.RData",rholist[i]))
    ### Method 3: multi-level matrix decomposition
    matrix=k_unfold(as.tensor(data),1)@data
    d=dim(matrix)
    A_1=  as.matrix(randortho(d[1])[,1:r[1]])
    A_2 = matrix(rnorm(d[2]*min(r[1],r[2]*r[3]),0,1),nrow=d[2])
    result <- fit_ordinal_matrix(matrix,A_1,A_2,omega=TRUE,alpha=TRUE)
    save(result,file=sprintf("../Figure6/alpha_R/output_level%.2f_matrix.RData",rholist[i]))
}

######################## analysis ################################
for(i in 1:length(rholist)){
######### Method 1: ordinal tensor ##################
load(sprintf("../Figure6/alpha_R/output_level%.2f_ordinal.RData",rholist[i]))
pred=estimation(result$theta,result$omega,"mode")@data
true=estimation(theta,omega,"mode")@data
Err[i,1]=mean(abs(true-pred))
Err[i,2]=mean(pred!=true)
pred=estimation(result$theta,result$omega,"median")@data
true=estimation(theta,omega,"median")@data
Err[i,3]=mean(abs(true-pred))
Err[i,4]=mean(pred!=true)
######### Method 2: continous tensor##################
load(sprintf("../Figure6/alpha_R/output_level%.2f_cont.RData",rholist[i]))
true=estimation(theta,omega,"mode")@data
Err2[i,1]=mean(abs(true-result2$theta))
Err2[i,2]=mean(true!=round(result2$theta))
true=estimation(theta,omega,"median")@data
Err2[i,3]=mean(abs(true-result2$theta))
Err2[i,4]=mean(true!=round(result2$theta))
### Method 3: multi-level matrix decomposition
load(sprintf("../Figure6/alpha_R/output_level%.2f_matrix.RData",rholist[i]))
pred=estimation(result$theta,result$omega,"mode")
true=estimation(theta,omega,"mode")@data
Err3[i,1]=mean(abs(c(true)-c(pred)))
Err3[i,2]=mean(c(pred)!=c(true))
pred=estimation(result$theta,result$omega,"median")
true=estimation(theta,omega,"median")@data
Err3[i,3]=mean(abs(c(true)-c(pred)))
Err3[i,4]=mean(c(pred)!=c(true))
}


### Methods 4-6: 1-bit Tensor decomposition. 
### Step 1. run 1-bit_simulation.m in Matlab. The outputs are saved in alpha_matlab/....
### Step 2. analysis based on the output
for(i in 1:length(rholist)){
load(sprintf("../Figure6/alpha_R/input_level%.2f.RData",rholist[i])) ## load theta and omega
############################################################
############optioin 1: categorical encoding############
output=readMat(sprintf("../Figure6/alpha_matlab/output_level%.2f_categorical.mat",rholist[i]))
prob=one_to_M(logistic(output$T.recovered),K,type="categorical")
### mode estimation; accuracy
est1=apply(prob,1:3,function(x)which.max(x))
true=estimation(theta,omega,"mode")@data
Err4[i,1]=mean(abs(true-est1))
Err4[i,2]=mean(abs(true!=est1))
###median estimation; accuracy
est2=apply(prob,1:3,function(x)which(cumsum(x)>=0.5)[1])
true=estimation(theta,omega,"median")@data
Err4[i,3]=mean(abs(true-est2))
Err4[i,4]=mean(abs(true!=est2))
############################################################
############ optioin 2: cumulative encoding############
output=readMat(sprintf("../Figure6/alpha_matlab/output_level%.2f_cumulative.mat",rholist[i]))
prob=one_to_M(logistic(output$T.recovered),K,type="cumulative")
## mode estimation; accuracy
est1=apply(prob,1:3,function(x)which.max(x))
true=estimation(theta,omega,"mode")@data
Err5[i,1]=mean(abs(true-est1))
Err5[i,2]=mean(abs(true!=est1))
###median estimation; accuracy
est2=apply(prob,1:3,function(x)which(cumsum(x)>=0.5)[1])
true=estimation(theta,omega,"median")@data
Err5[i,3]=mean(abs(true-est2))
Err5[i,4]=mean(abs(true!=est2))
############################################################
############ optioin 3: sign encoding############
output=readMat(sprintf("../Figure6/alpha_matlab/output_level%.2f_sign.mat",rholist[i]))
est=output$T.recovered
true=estimation(theta,omega,"mode")@data
Err6[i,1]=mean(abs(true-est))
Err6[i,2]=mean(abs(true!=round(est)))
true=estimation(theta,omega,"median")@data
Err6[i,3]=mean(abs(true-est))
Err6[i,4]=mean(abs(true!=round(est)))

#true=estimation(theta,omega,"mean")@data
#est1=apply(prob_total,1:3,function(x){return(sum(x*c(1:length(x))))})
}
### end evaluation


Output=data.frame(MSE=c(t(apply(MSE_raw,c(1:2),mean))),rank=as.factor(rep(rlist,rep(length(rholist),length(rlist)))),rho=rep(rholist,length(rlist)),sd=c(t(apply(MSE_raw,c(1:2),sd))))
save(MSE_raw,Output,file="Figure4.RData")

################ plot ############
load("Figure4.RData")
library(ggplot2)
Output$rank=as.factor(Output$rank)
levels(Output$rank)=c("r=3","r=5","r=8")

p=ggplot(data=Output,aes(x=rho,y=MSE))+geom_line(aes(color=rank),size=1.5)+geom_point(aes(shape=rank))+ labs(x='Number of Observed Entries', y='relative MSE',size=1.2)+theme(text = element_text(size=rel(4)),legend.text = element_text(size = 15),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))

#p=p+stat_function(fun=function(d)10^5/(d^2),lty=2)
g <- ggplot_build(p)
col=unlist(unique(g$data[[1]]["colour"]))[c(1,2,3)]

p=p+geom_errorbar(aes(ymin=MSE-sd, ymax=MSE+sd), width=0.01,
position=position_dodge(0.05),color=col[rep(rep(1:length(rlist)),length(rholist))],size=1)+geom_point(aes(shape=rank),size=1.5)

pdf("error_sample.pdf",width=5,height=4)
p
dev.off()








