### Cross validation ######
source("functions.R")
load("../data/dti_brain.RData")
## I changed the number of iteration from 50 to 10 in fit_continous/fit_ordinal
set.seed(18)
i1=i2=i3=i4=i5=NULL

for(i in 1:3){ ### 5-fold sampling while preserving the lable proportion
    label=sample(which(tensor==i),length(which(tensor==i))) 
    l = length(label)%/%5+1
    i1 = c(i1,label[1:l])
    i2 = c(i2,label[(l+1):(2*l)])
    i3 = c(i3,label[(2*l+1):(3*l)])
    i4 = c(i4,label[(3*l+1):(4*l)])
    i5 = c(i5,label[(4*l+1):length(ind)])
}
cindex = list(i1,i2,i3,i4,i5)

# `index' is for 1:5  
# I used following lines when I put codes in the server in sh file 
# args index ==$SLURM_ARRAY_TASK_ID" filename.R
test_index = cindex[[index]]
train_index = setdiff(1:length(tensor),test_index)
train_tensor = tensor
train_tensor[test_index] = NA


#####################################################################################
################# Method 1: 1-bit tensor completion CV ################################
###### Step 1. prepare input for matlab
input=array(0,dim=c(dim(train_tensor),2))
for(k in 1:2){
    prepare=train_tensor  
    prepare[train_tensor!=k]=-1
    prepare[train_tensor==k]=1
    input[,,,k]=prepare
}
E=array(1,dim=dim(input))
E[is.na(input)]=0
input[E==0]=0 ## missingness encoded as zero in 1-bit tensor completion algorithm. 
writeMat(paste("HPC/input_CV23_23_8_",index,".mat",sep = ""),data=input,E=E)

######Step 2. open matlab, run the script HPC_1bit.m. 
######Step 3. analysis after getting the above output files
BCV = as.data.frame(matrix(nrow = 5, ncol = 4))
names(BCV) = c("mode_MAD","mode_MCR","median_MAD","median_MCR")
for (i in 1:5) {
test_index = cindex[[i]]
output=readMat(paste("HPC/output_CV23_23_8_",index,".mat",sep = ""))
p_sum=apply(logistic(output$T.recovered),c(1:3),sum) ## negative sign per convention of 1-bit tensor algorithm
prob_total=array(0,dim=c(dim(data),3)) ## Order-4 probability tensor 
prob_total[,,,1:2]=logistic(output$T.recovered)
prob_total[,,,3]=1-p_sum
est1=apply(prob_total,1:3,function(x)which.max(x)) ## prediction via mode
BCV[i,1] = mean(est1[test_index]-tensor[test_index])) ##MAD
BCV[i,2] = error_rate = mean(est1[test_index]!=tensor[test_index]) ##MCR
est2=apply(prob_total,1:3,function(x)which(cumsum(x)>=0.5)[1])## prediction via median
BCV[i,3] = mean(est1[test_index]-tensor[test_index])) ##MAD
BCV[i,4] = error_rate = mean(est1[test_index]!=tensor[test_index]) ##MCR
}

#####################################################################################
################# Method 2: Continuous Tucker decomposition CV ################################
d = dim(tensor)
set.seed(18)
r = c(23,23,8)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
result = fit_continuous(train_tensor,C,A_1,A_2,A_3,alpha=TRUE)
save(result,file = paste("CV23_23_8_",index,".RData",sep = ""))

########## Analysis after getting the above output files ################
CV = as.data.frame(matrix(nrow = 5, ncol = 4))
names(CV) = c("mode_MAD","mode_MCR","median_MAD","median_MCR")
for (i in 1:5) {
  test_index = cindex[[i]]
  load(paste("CV23_23_8_",i,".RData",sep = ""))
  ##theta =ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3)@data
  theta = result$theta # alternatively, direct output the theta from result object
  CV[i,1] = CV[i,3]=mean(abs(theta[test_index]-tensor[test_index]))
  CV[i,2] = CV[i,4]=error_rate = mean(round(theta)[test_index]!=tensor[test_index])
}
apply(CV,2,mean)

##################### Method 3: ordinal tensor decomposition CV#################
d = dim(tensor)
set.seed(18)
r = c(23,23,8)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
result = fit_ordinal(train_tensor,C,A_1,A_2,A_3)
save(result,file = paste("OCV23_23_8_",index,".RData",sep = ""))
########## Analysis after getting the above output files ################

OCV = as.data.frame(matrix(nrow = 5, ncol = 4))
names(OCV) = c("mode_MAD","mode_MCR","median_MAD","median_MCR")
for (i in 1:5) {
  test_index = cindex[[i]]
  load(paste("OCV23_23_8_",i,".RData",sep = ""))
  theta = result$theta # alternatively, direct output the theta from result object
  #theta =ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3)@data
  out = estimation(result$theta ,result$omega,type="mode")@data # we can also use other estimator via type="mode","mean", or "median"
  OCV[i,1] = mean(abs(out[test_index]-tensor[test_index]))
  OCV[i,2] = error_rate = mean(out[test_index]!=tensor[test_index])
  out2 = estimation(result$theta ,result$omega,type="median")@data # we can also use other estimator via type="mode","mean", or "median"
  OCV[i,3] = mean(abs(out2[test_index]-tensor[test_index]))
  OCV[i,4] = error_rate = mean(out2[test_index]!=tensor[test_index])
}
apply(OCV,2,mean)

library(kableExtra)
kable(OCV,'latex')
