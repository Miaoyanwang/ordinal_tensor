addpath('/workspace/miaoyan/TensorCompletion_1bit_noisy-master/tensorlab/') 
addpath('/workspace/miaoyan/TensorCompletion_1bit_noisy-master') 

p_cc=0.1; 
doingmax=1; 
doing1bit=1; 
f      = @(x) 1./(1+exp(-x));
fprime = @(x) exp(-x)./((1+exp(-x)).*(1+exp(-x)));
iter=10; 
alternations=10;
infbound=10; 

Klist=[2,3,4,5,6,7];
%%% start simulation
for nsim = 1:30
for i = 1:6
%% categorical encoding
input=importdata(sprintf("K_matlab/input_level%.1f_sim%d_categorical.mat",Klist(i),nsim));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
save(sprintf("K_matlab/output_level%.1f_sim%d_categorical.mat",Klist(i),nsim),'T_recovered');

%% sign encoding
input=importdata(sprintf("K_matlab/input_level%.1f_sim%d_sign.mat",Klist(i),nsim));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
T_recovered=T_recovered/norm(T_recovered(:))*input.scale;
T_recovered=T_recovered+input.ave;
save(sprintf("K_matlab/output_level%.1f_sim%d_sign.mat",Klist(i),nsim),'T_recovered');
end
end
%%% end simulation
