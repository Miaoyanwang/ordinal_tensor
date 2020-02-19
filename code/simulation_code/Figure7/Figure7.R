### Figure: MSE vs. number of levels (K) under different methods ##
source("../../functions.R")
library("R.matlab")
set.seed(1)
d=rep(20,3)
r=rep(5,3)
ntotal=30
Klist=2:7
rho=0.8
Err1=Err2=Err3=Err4=Err5=array(0,dim=c(length(Klist),4,ntotal))

for(nsim in 1:ntotal){
B_1 = as.matrix(randortho(d[1])[,1:r[1]])
B_2 = as.matrix(randortho(d[2])[,1:r[2]])
B_3 = as.matrix(randortho(d[3])[,1:r[3]])
D = rand_tensor(modes = r)
theta = ttl(D,list(B_1,B_2,B_3),ms=1:3)@data
alpha=10
theta = theta/max(abs(theta))*alpha ## specify signal = 10
for(i in 1:length(Klist)){
K=Klist[i]
omega = logit((1:(K-1))/K) ## equal spaced on the probit scale
#omega=0   
#K=length(omega)+1
##omega=logit((1:2)/3)
ttnsr <- realization(theta,omega)@data
 
######## prepare input files ############################
## generate observation for fit_ordinal
data=array(NA,dim=dim(ttnsr))
Obs=sample(1:prod(d),prod(d)*rho,replace=F)
data[Obs]=ttnsr[Obs] 
save(data,theta,omega,file=sprintf("../Figure8/K_R/input_level%.1f_sim%d.RData",Klist[i],nsim))
#### generate observation for 1-bit tensor
input=M_to_one(data,K,"category")
writeMat(sprintf("../Figure8/K_matlab/input_level%.1f_sim%d_categorical.mat",Klist[i],nsim),data=input$binary,E=input$E)
input=M_to_one(data,K,"sign")
writeMat(sprintf("../Figure8/K_matlab/input_level%.1f_sim%d_sign.mat",Klist[i],nsim),data=input$binary,E=input$E,ave=input$ave,scale=input$scale)
}
}

######################## decomposition ################################
for(nsim in 1:ntotal){
for(i in 1:length(Klist)){
    K=Klist[i]
    load(sprintf("../Figure8/K_R/input_level%.1f_sim%d.RData",Klist[i],nsim))
    ######### Method 1: ordinal tensor ##################
    d=dim(data)
    A_1 = as.matrix(randortho(d[1])[,1:r[1]])
    A_2 = as.matrix(randortho(d[2])[,1:r[2]])
    A_3 = as.matrix(randortho(d[3])[,1:r[3]])
    C = rand_tensor(modes = r)
    result <- fit_ordinal(data,C,A_1,A_2,A_3,omega=TRUE,alpha=TRUE)
    save(result,file=sprintf("../Figure8/K_R/output_level%.1f_sim%d_ordinal.RData",Klist[i],nsim))
    ######### Method 2: continous tensor##################
    result2 <- fit_continuous(data,C,A_1,A_2,A_3,alpha=TRUE)
    save(result2,file=sprintf("../Figure8/K_R/output_level%.1f_sim%d_cont.RData",Klist[i],nsim))
    ### Method 3: multi-level matrix decomposition
    matrix=k_unfold(as.tensor(data),1)@data
    d=dim(matrix)
    A_1=  as.matrix(randortho(d[1])[,1:r[1]])
    A_2 = matrix(rnorm(d[2]*min(r[1],r[2]*r[3]),0,1),nrow=d[2])
    result <- fit_ordinal_matrix(matrix,A_1,A_2,omega=TRUE,alpha=TRUE)
    save(result,file=sprintf("../Figure8/K_R/output_level%.1f_sim%d_matrix.RData",Klist[i],nsim))
}
}


######################## analysis ################################
for(nsim in 1:ntotal){
    for(i in 1:length(Klist)){
        load(sprintf("../Figure8/K_R/input_level%.1f_sim%d.RData",Klist[i],nsim))
        ######### Method 1: ordinal tensor ##################
        load(sprintf("../Figure8/K_R/output_level%.1f_sim%d_ordinal.RData",Klist[i],nsim))
        pred=estimation(result$theta,result$omega,"mode")@data
        true=estimation(theta,omega,"mode")@data
        Err1[i,1,nsim]=mean(abs(true-pred))
        Err1[i,2,nsim]=mean(pred!=true)
        pred=estimation(result$theta,result$omega,"median")@data
        true=estimation(theta,omega,"median")@data
        Err1[i,3,nsim]=mean(abs(true-pred))
        Err1[i,4,nsim]=mean(pred!=true)
        ######### Method 2: continous tensor##################
        load(sprintf("../Figure8/K_R/output_level%.1f_sim%d_cont.RData",Klist[i],nsim))
        true=estimation(theta,omega,"mode")@data
        Err2[i,1,nsim]=mean(abs(true-result2$theta))
        Err2[i,2,nsim]=mean(true!=round(result2$theta))
        true=estimation(theta,omega,"median")@data
        Err2[i,3,nsim]=mean(abs(true-result2$theta))
        Err2[i,4,nsim]=mean(true!=round(result2$theta))
        ### Method 3: multi-level matrix decomposition
        load(sprintf("../Figure8/K_R/output_level%.1f_sim%d_matrix.RData",Klist[i],nsim))
        pred=estimation(result$theta,result$omega,"mode")
        true=estimation(theta,omega,"mode")@data
        Err3[i,1,nsim]=mean(abs(c(true)-c(pred)))
        Err3[i,2,nsim]=mean(c(pred)!=c(true))
        pred=estimation(result$theta,result$omega,"median")
        true=estimation(theta,omega,"median")@data
        Err3[i,3,nsim]=mean(abs(c(true)-c(pred)))
        Err3[i,4,nsim]=mean(c(pred)!=c(true))
    }
}

### Methods 4-6: 1-bit Tensor decomposition. 
### Step 1. run 1-bit_simulation.m in Matlab. The outputs are saved in K_matlab/....
### Step 2. analysis based on the output
for(nsim in 1:ntotal){
    for(i in 1:length(Klist)){
        load(sprintf("../Figure8/K_R/input_level%.1f_sim%d.RData",Klist[i],nsim))
        ############################################################
        ############optioin 1: categorical encoding############
        output=readMat(sprintf("../Figure8/K_matlab/output_level%.1f_sim%d_categorical.mat",Klist[i],nsim))
        prob=one_to_M(logistic(output$T.recovered),Klist[i],type="categorical")
        ### mode estimation; accuracy
        est1=apply(prob,1:3,function(x)which.max(x))
        true=estimation(theta,omega,"mode")@data
        Err4[i,1,nsim]=mean(abs(true-est1))
        Err4[i,2,nsim]=mean(abs(true!=est1))
        ###median estimation; accuracy
        est2=apply(prob,1:3,function(x)which(cumsum(x)>=0.5)[1])
        true=estimation(theta,omega,"median")@data
        Err4[i,3,nsim]=mean(abs(true-est2))
        Err4[i,4,nsim]=mean(abs(true!=est2))
        ############################################################
        ############ optioin 2: sign encoding############
        output=readMat(sprintf("../Figure8/K_matlab/output_level%.1f_sim%d_sign.mat",Klist[i],nsim))
        est=output$T.recovered
        true=estimation(theta,omega,"mode")@data
        Err5[i,1,nsim]=mean(abs(true-est))
        Err5[i,2,nsim]=mean(abs(true!=round(est)))
        true=estimation(theta,omega,"median")@data
        Err5[i,3,nsim]=mean(abs(true-est))
        Err5[i,4,nsim]=mean(abs(true!=round(est)))
    }
}
### end evaluation
dimnames(Err1)[[2]]=dimnames(Err2)[[2]]=dimnames(Err3)[[2]]=dimnames(Err4)[[2]]=dimnames(Err5)[[2]]=c("MAD_for_mode","MCR_for_mode","MAD_for_median","MCR_for_median")
save(Err1,Err2,Err3,Err4,Err5,file="Figure8.RData")
