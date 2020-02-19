p_cc=0.1; 
doingmax=1; 
doing1bit=1; 
f      = @(x) 1./(1+exp(-x));
fprime = @(x) exp(-x)./((1+exp(-x)).*(1+exp(-x)));
iter=10; 
alternations=10;
infbound=10; 

index=[1,2,3,4,5];

%%% begin HPC analysis
for i = 1:5
input=importdata(sprintf("HPC/input_CV23_23_8_%.0f.mat",index(i)));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
save(sprintf("HPC/output_CV23_23_8_%.0f.mat",index(i),'T_recovered');
end
%%% end HPC analysis

     
     


