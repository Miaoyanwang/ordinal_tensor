p_cc=0.1; 
doingmax=1; 
doing1bit=1; 
f      = @(x) 1./(1+exp(-x));
fprime = @(x) exp(-x)./((1+exp(-x)).*(1+exp(-x)));
iter=10; 
alternations=10;
infbound=10; 

rholist=[0.4,0.5,0.6,0.7,0.8,0.9,1.0];
%%% start simulation
for i = 1:7
%% categorical encoding
input=importdata(sprintf("simulation_figure/Figure6/alpha_matlab/input_level%.2f_categorical.mat",rholist(i)));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
save(sprintf("simulation_figure/Figure6/alpha_matlab/output_level%.2f_categorical.mat",rholist(i)),'T_recovered');

%% cumulative encoding
input=importdata(sprintf("simulation_figure/Figure6/alpha_matlab/input_level%.2f_cumulative.mat",rholist(i)));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
save(sprintf("simulation_figure/Figure6/alpha_matlab/output_level%.2f_cumulative.mat",rholist(i)),'T_recovered');

%% sign encoding
input=importdata(sprintf("simulation_figure/Figure6/alpha_matlab/input_level%.2f_sign.mat",rholist(i)));
T_recovered=main_noisyor1bit_TC(input.data,input.E,p_cc,doingmax,doing1bit,f,fprime,iter,alternations,infbound);
T_recovered=T_recovered/norm(T_recovered(:))*input.scale;
T_recovered=T_recovered+input.ave;
save(sprintf("simulation_figure/Figure6/alpha_matlab/output_level%.2f_sign.mat",rholist(i)),'T_recovered');
end
%%% end simulation
