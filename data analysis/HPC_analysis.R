### Cross validation ######
source("functions.R")
load("../data/dti_brain.RData")


### 5 fold ################ To Chanwoo: Please changes the sampling per instruction ####

# `index' is for 1:5  
# I used following lines when I put codes in the server in sh file 
# args index ==$SLURM_ARRAY_TASK_ID" filename.R
# i: repetition j: j-th testset

indset =matrix(nrow =50,ncol=2)
s = 0
for(i in 1:10){
  for(j in 1:5){
    # for(k in 1:2){
      s = s+1
      indset[s,] = c(i,j)
    # }
  }
}

i = indset[index,][1]
j = indset[index,][2]
# k = indset[index,][3]

set.seed(i)
l1 = split(sample(which(tensor==1),length(which(tensor==1))),as.factor(1:5))
l2 = split(sample(which(tensor==2),length(which(tensor==2))),as.factor(1:5))
l3 = split(sample(which(tensor==3),length(which(tensor==3))),as.factor(1:5))

cindex = list()
for (k in 1:5) {
  cindex[[k]] = c(l1[[k]],l2[[k]],l3[[k]])
}

test_index = cindex[[j]]
train_index = setdiff(1:length(tensor),test_index)
train_tensor = tensor
train_tensor[test_index] = NA

#####################################################################################
################# Continuous Tucker decomposition CV ################################
d = dim(tensor)
r = c(23,23,8)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = A_1
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)

result = fit_continuous(train_tensor,C,A_1,A_2,A_3)
save(result,file = paste("CV_conti_",i,"_",j,".RData",sep = ""))


########## Analysis after getting the above output files ################
   
CV = as.data.frame(matrix(nrow = 50, ncol = 2))
names(CV) = c("MAE","MCR")
s = 0
for(i in 1:10){
  for (j in 1:5) {
    s = s+1
    set.seed(i)
    l1 = split(sample(which(tensor==1),length(which(tensor==1))),as.factor(1:5))
    l2 = split(sample(which(tensor==2),length(which(tensor==2))),as.factor(1:5))
    l3 = split(sample(which(tensor==3),length(which(tensor==3))),as.factor(1:5))
    
    cindex = list()
    for (k in 1:5) {
      cindex[[k]] = c(l1[[k]],l2[[k]],l3[[k]])
    }
    
    test_index = cindex[[j]]
    test_index = cindex[[j]]
    load(paste("CV_conti_",i,"_",j,".RData",sep = ""))
    theta = result$theta
    CV[s,1] = mean(abs(theta[test_index]-tensor[test_index]))
    CV[s,2] =  error_rate = mean(round(theta)[test_index]!=tensor[test_index])
  }
}




##################### ordinal glm tucker decomposition CV#################
d = dim(tensor)
r = c(23,23,8)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = A_1
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
result = fit_ordinal(train_tensor,C,A_1,A_2,A_3)
save(result,file = paste("CV_ordinal_",i,"_",j,".RData",sep = ""))


########## Analysis after getting the above output files ################
### with median estimation
OCVmedian = as.data.frame(matrix(nrow = 50, ncol = 2))
names(OCVmedian) = c("MAE","MCR")
s = 0
for(i in 1:10){
  for (j in 1:5) {
    s = s+1
    set.seed(i)
    l1 = split(sample(which(tensor==1),length(which(tensor==1))),as.factor(1:5))
    l2 = split(sample(which(tensor==2),length(which(tensor==2))),as.factor(1:5))
    l3 = split(sample(which(tensor==3),length(which(tensor==3))),as.factor(1:5))
    
    cindex = list()
    for (k in 1:5) {
      cindex[[k]] = c(l1[[k]],l2[[k]],l3[[k]])
    }
    
    test_index = cindex[[j]]
    test_index = cindex[[j]]
    load(paste("CV_ordinal_",i,"_",j,".RData",sep = ""))
    theta = result$theta
    out = estimation(theta,result$omega,type="median")@data # we can also use other estimator via type="mode","mean", or "median"
    OCVmedian[s,1] = mean(abs(out[test_index]-tensor[test_index]))
    OCVmedian[s,2] = error_rate = mean(out[test_index]!=tensor[test_index])
  }
}

### with mode estimaiton
OCVmode = as.data.frame(matrix(nrow = 50, ncol = 2))
names(OCVmode) = c("MAE","MCR")
s = 0
for(i in 1:10){
  for (j in 1:5) {
    s = s+1
    set.seed(i)
    l1 = split(sample(which(tensor==1),length(which(tensor==1))),as.factor(1:5))
    l2 = split(sample(which(tensor==2),length(which(tensor==2))),as.factor(1:5))
    l3 = split(sample(which(tensor==3),length(which(tensor==3))),as.factor(1:5))
    
    cindex = list()
    for (k in 1:5) {
      cindex[[k]] = c(l1[[k]],l2[[k]],l3[[k]])
    }
    
    test_index = cindex[[j]]
    test_index = cindex[[j]]
    load(paste("CV_ordinal_",i,"_",j,".RData",sep = ""))
    theta = result$theta
    out = estimation(theta,result$omega,type="mode")@data # we can also use other estimator via type="mode","mean", or "median"
    OCVmode[s,1] = mean(abs(out[test_index]-tensor[test_index]))
    OCVmode[s,2] = error_rate = mean(out[test_index]!=tensor[test_index])
  }
}

###################### Getting a data summary##############################
################### continuous decomposition ###############
cv_rep = matrix(nrow =10,ncol = 2)
for(k in 1:10){
  cv_rep[k,] = apply(CV[((k-1)*5+1):(k*5),],2,mean)  
}

cv_summary = rbind(apply(cv_rep,2,mean),apply(cv_rep,2,sd))
colnames(cv_summary) = c("MAD","MCR")
rownames(cv_summary) = c("mean","sd")

################### ordinal decomposition with mode ######## 
ocvmode_rep = matrix(nrow =10,ncol = 2)
for(k in 1:10){
  ocvmode_rep[k,] = apply(OCVmode[((k-1)*5+1):(k*5),],2,mean)  
}

ocvmode_summary = rbind(apply(ocvmode_rep,2,mean),apply(ocvmode_rep,2,sd))
colnames(ocvmode_summary) = c("MAD","MCR")
rownames(ocvmode_summary) = c("mean","sd")

################### ordinal decomposition with median ######
ocvmedian_rep = matrix(nrow =10,ncol = 2)
for(k in 1:10){
  ocvmedian_rep[k,] = apply(OCVmedian[((k-1)*5+1):(k*5),],2,mean)  
}

ocvmedian_summary = rbind(apply(ocvmedian_rep,2,mean),apply(ocvmedian_rep,2,sd))
colnames(ocvmedian_summary) = c("MAD","MCR")
rownames(ocvmedian_summary) = c("mean","sd")


