### Cross validation ######
source("functions.R")
load("dti_brain.RData")
set.seed(18)
ind = sample(1:length(tensor),length(tensor))

### 5 fold ################ Please make changes to the sampling #####
l = length(tensor)%/%5+1
i1 = ind[1:l]
i2 = ind[(l+1):(2*l)]
i3 = ind[(2*l+1):(3*l)]
i4 = ind[(3*l+1):(4*l)]
i5 = ind[(4*l+1):length(tensor)]
cindex = list(i1,i2,i3,i4,i5)

# `index' is for 1:5  
# I used following lines when I put codes in the server in sh file 
# args index ==$SLURM_ARRAY_TASK_ID" filename.R
test_index = cindex[[index]]
train_index = setdiff(ind,test_index)
train_tensor = tensor
train_tensor[test_index] = NA

#####################################################################################
################# Continuous Tucker decomposition CV ################################
d = dim(tensor)
r = c(23,23,8)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
result = fit_continuous(train_tensor,C,A_1,A_2,A_3,alpha=TRUE)
save(result,file = paste("CV23_23_8_",index,".RData",sep = ""))

########## Analysis after getting the above output files ################
CV = as.data.frame(matrix(nrow = 5, ncol = 3))
names(CV) = c("MSE","MAE","Error_rate")
for (i in 1:5) {
  test_index = cindex[[i]]
  load(paste("CV23_23_8_",i,".RData",sep = ""))
  theta =ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3)@data
  ##theta = result$theta # alternatively, direct output the theta from result object
  CV[i,1] = mean((round(theta)[test_index]-tensor[test_index])^2)
  CV[i,2] = mean(abs(round(theta)[test_index]-tensor[test_index]))
  CV[i,3] = error_rate = mean(round(theta)[test_index]!=tensor[test_index])
}

##################### ordinal glm tucker decomposition CV#################
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

OCV = as.data.frame(matrix(nrow = 5, ncol = 3))
names(OCV) = c("MSE","MAE","Error_rate")
for (i in 1:5) {
  test_index = cindex[[i]]
  load(paste("OCV23_23_8_",i,".RData",sep = ""))
  #theta = result$theta # alternatively, direct output the theta from result object
  theta =ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3)@data
  out = estimation(theta,result$omega,type="mode")@data # we can also use other estimator via type="mode","mean", or "median"
  OCV[i,1] = mean((round(out)[test_index]-tensor[test_index])^2)
  OCV[i,2] = mean(abs(round(out)[test_index]-tensor[test_index]))
  OCV[i,3] = error_rate = mean(round(out[test_index])!=tensor[test_index])
}
apply(OCV,2,mean)

library(kableExtra)
kable(OCV,'latex')
