load("dti_brain.RData")
load("rank_24_24_8.RData")
theta = make_theta(result)
thetan = tucker(theta,rank = result$C@modes)

# clust1 is clustering from mode 1, cluster 2 is clustering from mode 2
clust1 = kmeans(thetan$U[[1]]%*%diag(diag(k_unfold(thetan$Z,1)@data)),7,nstart = 20)
clust2 = kmeans(thetan$U[[2]]%*%diag(diag(k_unfold(thetan$Z,2)@data)),7,nstart = 20)
par(mfrow = c(1,2))
plot(clust1$cluster,ylab = "group")
plot(clust2$cluster,ylab = "group")

# making connection matrix for each cluster group 
clust1$size
set = clust1$cluster
braincluster = list()
for (i in 1:length(k)) {
  braincluster[[i]] = matrix(0,nrow =68,ncol = 68)
  ci = which(set==i)
  for(p in ci){
    for(q in ci){
      braincluster[[i]][p,q]=1
    }
  }
}
for(i in 1:length(k)){
  write(braincluster[[i]],file = paste("cluster",i,".edge",sep = ""))
}

write(braincluster,file = "cluster.txt")



