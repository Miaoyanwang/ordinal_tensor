

make_theta = function(result){
  return(ttl(result$C,list(result$A_1,result$A_2,result$A_3),ms=1:3))
}



### Cross validation ######
load("dti_brain.RData")
set.seed(18)
ind = sample(1:length(tensor),length(tensor))

### 5 fold ################
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
train_tensor[test_index] = -1

#####################################################################################
################# Continuous Tucker decomposition CV ################################
W = array(0, dim(tensor))
W[train_index] = 1

d = dim(tensor)
r = c(25,25,7)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
result = tucker_missing(A_1,A_2,A_3,C,train_tensor,W)
save(result,file = paste("CV25_25_7_",index,".RData",sep = ""))

########## Analysis after getting the above output files ################

CV = as.data.frame(matrix(nrow = 5, ncol = 3))
names(CV) = c("MSE","MAE","Error_rate")
for (i in 1:5) {
  
  load(paste("CV23_23_8_",i,".RData",sep = ""))
  theta = make_theta(result)@data
  CV[i,1] = mean((round(theta)[test_index]-tensor[test_index])^2)
  CV[i,2] = mean(abs(round(theta)[test_index]-tensor[test_index]))
  CV[i,3] = error_rate = 1-length(which(round(theta)[test_index]- 
                                          tensor[test_index]==0))/125773
}

##################### ordinal glm tucker decomposition CV#################
d = dim(tensor)
r = c(25,25,7)
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
result = fit_ordinal(train_tensor,C,A_1,A_2,A_3)
save(result,file = paste("OCV25_25_7_",index,".RData",sep = ""))

########## Analysis after getting the above output files ################

OCV = as.data.frame(matrix(nrow = 5, ncol = 3))
names(OCV) = c("MSE","MAE","Error_rate")
for (i in 1:5) {
  load(paste("OCV23_23_8_",i,".RData",sep = ""))
  theta = make_theta(result)@data
  out = estimation(theta,result$omega)@data # we can also use mestimation 
  OCV[i,1] = mean((out[test_index]-tensor[test_index])^2)
  OCV[i,2] = mean(abs(out[test_index]-tensor[test_index]))
  OCV[i,3] = error_rate = 1-length(which(out[test_index]- 
                                          tensor[test_index]==0))/125773
}


library(kableExtra)
kable(OCV,'latex')
