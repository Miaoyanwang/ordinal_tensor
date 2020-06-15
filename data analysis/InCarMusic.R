### Cross validation ######
source("functions.R")
load("../data/InCar_Music.RData")
tensor[tensor==-1]=NA
###
#m=tensor
#m[is.na(m)]=0
#m[m!=0]=1
#count1=apply(m,1,mean)
#count2=apply(m,2,mean)
#count3=apply(m,3,mean)
#tensor=tensor[count1>=0.015,count2>=0.015,count3>=0.015]
###

tensor[tensor==-1]=NA
ind=which(is.na(tensor)==F)
d=dim(tensor)
set.seed(1)
library("R.matlab")
Err1=Err2=Err3=Err4=Err5=Err6=matrix(0,nrow=5,ncol=10)

### 5 fold ################
for(nsim in 1:10){
i1=i2=i3=i4=i5=NULL

for(i in 1:5){ ### 5-fold sampling while preserving the lable proportion
    label=sample(which(tensor==i),length(which(tensor==i))) 
    l = length(label)%/%5+1
    i1 = c(i1,label[1:l])
    i2 = c(i2,label[(l+1):(2*l)])
    i3 = c(i3,label[(2*l+1):(3*l)])
    i4 = c(i4,label[(3*l+1):(4*l)])
    i5 = c(i5,label[(4*l+1):length(label)])
}
cindex = list(i1,i2,i3,i4,i5)


# `index' is for 1:5  
# I used following lines when I put codes in the server in sh file 
# args index ==$SLURM_ARRAY_TASK_ID" filename.R
for(index in 1:5){
test_index = cindex[[index]]
train_index = setdiff(ind,test_index)
train_tensor = tensor
train_tensor[test_index] = NA

save(test_index,train_index,file=sprintf("InCarMusic_subset/input_%d_sim%d.RData",index,nsim))

### Method 1. 1-bit tensor completion
input=M_to_one(train_tensor,5,type="sign")
writeMat(sprintf("InCarMusic_subset/input_CV%d_sim%d.mat",index,nsim),data=input$binary,E=input$E,ave=input$ave,scale=input$scale)
}
}

#####################################################################################
################# Methods 2 and 3: Tucker decomposition CV ################################
for(nsim in 1:10){
    for(index in 1:5){
   load(sprintf("InCarMusic_subset/input_%d_sim%d.RData",index,nsim))     
   train_tensor = tensor
   train_tensor[test_index] = NA
   
d=dim(train_tensor)
r=c(2,2,2)
A_1 = as.matrix(randortho(d[1])[,1:r[1]])
A_2 = as.matrix(randortho(d[2])[,1:r[2]])
A_3 = as.matrix(randortho(d[3])[,1:r[3]])
C = rand_tensor(modes = r)
result <- fit_ordinal(train_tensor,C,A_1,A_2,A_3,omega=TRUE,alpha=100)
save(result,file=sprintf("InCarMusic_subset/new_output_CV%d_sim%d_ordinal.RData",index,nsim))

result2 <- fit_continuous(train_tensor,C,A_1,A_2,A_3,alpha=100)
save(result2,file=sprintf("InCarMusic_subset/new_output_CV%d_sim%d_cont.RData",index,nsim))
}
}


######## analysis 
for(nsim in 1:10){
        for(index in 1:5){
        
        load(sprintf("InCarMusic_subset/input_%d_sim%d.RData",index,nsim))     
        train_tensor = tensor
        train_tensor[test_index] = NA
        
        ############ methods 1-2: ############
        
        load(sprintf("InCarMusic_subset/new_output_CV%d_sim%d_ordinal.RData",index,nsim))
        pred=estimation(result$theta,result$omega,"median")@data
        Err1[index,nsim]=mean(abs(pred[test_index]-tensor[test_index]))
        Err2[index,nsim]=mean(round(pred[test_index])!=tensor[test_index])
        load(sprintf("InCarMusic_subset/new_output_CV%d_sim%d_cont.RData",index,nsim))
        Err3[index,nsim]=mean(abs(result2$theta[test_index]-tensor[test_index]))
        Err4[index,nsim]=mean(round(result2$theta[test_index])!=tensor[test_index])
            
        ############ methods 3: sign encoding############
        output=readMat(sprintf("InCarMusic_subset/output_CV%d_sim%d.mat",index,nsim))
        est=output$T.recovered
        
        load(sprintf("InCarMusic_subset/input_%d_sim%d.RData",index,nsim))    

        Err5[index,nsim]=mean(abs(tensor[test_index]-est[test_index]))
        Err6[index,nsim]=mean(abs(tensor[test_index]!=round(est[test_index])))
 }
}



library(kableExtra)
kable(OCV,'latex')
