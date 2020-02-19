p_cc=0.1; 
doingmax=1; 
doing1bit=1; 
f      = @(x) 1./(1+exp(-x));
fprime = @(x) exp(-x)./((1+exp(-x)).*(1+exp(-x)));
iter=10; 
alternations=10;
infbound=10; 

%%% start simulation
for nsim = 1:10
for i = 1:5
%% categorical encoding
input=importdata(sprintf("InCarMusic_subset/input_CV%d_sim%d.mat",i,nsim));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
T_recovered=T_recovered/norm(T_recovered(:))*input.scale;
T_recovered=T_recovered+input.ave;
save(sprintf("InCarMusic_subset/output_CV%d_sim%d.mat",i,nsim),'T_recovered');
end
end
%%% end simulation

