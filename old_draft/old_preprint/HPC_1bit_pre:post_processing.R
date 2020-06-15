############ Method 3: 1-bit tensor completion CV ########
##The code below is for a single instance of train_tensor. Please modify the code when multiple instances (either from cross-validation or from replicates) are needed. 

## j-th testset. Need to prepare cindex beforhand. lines 5-8 are copied from lines 38-41 in HPC.analysis.R
test_index = cindex[[j]]
train_index = setdiff(1:length(tensor),test_index)
train_tensor = tensor
train_tensor[test_index] = NA

###### Step 1. prepare input for matlab
input=array(0,dim=c(dim(train_tensor),2)) ## suppose the data tensor is a 3-level ordinal tensor. 
for(k in 1:2){
    prepare=train_tensor  
    prepare[train_tensor!=k]=-1
    prepare[train_tensor==k]=1 ## dummy variable encoding; convert multi-level encoding to binary encodding. 
    input[,,,k]=prepare
}
E=array(1,dim=dim(input))
E[is.na(input)]=0
input[E==0]=0 ## missingness encoded as zero in 1-bit tensor completion algorithm. 
writeMat(paste("HPC/input.mat",data=input,E=E)

######Step 2. open matlab, run the script HPC_1bit.m. 

######Step 3. analysis after getting the above output files
BCV = as.data.frame(matrix(nrow = 5, ncol = 4))
names(BCV) = c("mode_MAD","mode_MCR","median_MAD","median_MCR")

for (i in 1:5) { ##5-fold CV
    test_index = cindex[[i]]
    output=readMat(paste("HPC/output.mat"))
    p_sum=apply(logistic(output$T.recovered),c(1:3),sum) ## estimated probability for labels 1 to K-1
    prob_total=array(0,dim=c(dim(data),3)) ## estimated probability for labels 1 to K
    prob_total[,,,1:2]=logistic(output$T.recovered)
    prob_total[,,,3]=1-p_sum
    est1=apply(prob_total,1:3,function(x)which.max(x)) ## prediction via mode
    BCV[i,1] = mean(est1[test_index]-tensor[test_index])) ##MAD
    BCV[i,2] = error_rate = mean(est1[test_index]!=tensor[test_index]) ##MCR
    est2=apply(prob_total,1:3,function(x)which(cumsum(x)>=0.5)[1])## prediction via median
    BCV[i,3] = mean(est1[test_index]-tensor[test_index])) ##MAD
    BCV[i,4] = error_rate = mean(est1[test_index]!=tensor[test_index]) ##MCR
}


