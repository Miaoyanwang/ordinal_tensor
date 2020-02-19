source("functions.R")

### This file contains the code of a simple demo for the tensor estimations and predictions.

### set tensor dimension and rank of the tensor
d=c(30,30,30) 
r=c(10,10,10)


### set true latent parameter
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)
theta = ttl(C,list(A_1,A_2,A_3),ms=1:3)@data*10

### set true cut-off points
omega = c(-1,1)

### Draw a randomly drawn tensor from above parameters with logistic model.
ttnsr = realization(theta,omega)


##########################################################################################
### Estimation of the latent parameter theta and cut-off points omega from the observed tensor.

### Initial points for `ordinal_fit' function ( we can use BIC to estimate rank )
A_1 = randortho(d[1])[,1:r[1]]
A_2 = randortho(d[2])[,1:r[2]]
A_3 = randortho(d[3])[,1:r[3]]
C = rand_tensor(modes = r)

### Estimation for parameters 
parameter_est = fit_ordinal(ttnsr@data,C,A_1,A_2,A_3,omega = TRUE,alpha = 100)

print(paste("MSE of theta: ",sum((parameter_est$theta-theta)^2)/prod(d),
            ",   MSE of omega: ",sum(parameter_est$omega-omega)^2/length(omega)))

### Predction from estimated parameters.
mode_prediction = estimation(parameter_est$theta,parameter_est$omega,"mode")

print(paste("MAD : ",mean(abs((mode_prediction-ttnsr)@data)),
            ",    MCR : ",mean(mode_prediction@data!=ttnsr@data)))

