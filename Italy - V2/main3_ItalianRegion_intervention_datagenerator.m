% * Codes for generating Monte Carlo sample using the fitted model
%  'Uncertainty in OR Epidemiological Modeling: A Global Sensitivity Approach', 
%   by Xuefei Lu and Emanuele Borgonovo, 2021
% 
% * Author: Xuefei Lu, xuefei.lu@ed.ac.uk
% * Date: July, 2021
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Cheynet, E. Generalized SEIR Epidemic Model (Fitting and Computation). Zenodo, 2020, doi:10.5281/ZENODO.3911854. 
% [2] Dong E, Du H, Gardner L. An interactive web-based dashboard to track COVID-19 in real time. Lancet Inf Dis. 
%     20(5):533-534. doi: 10.1016/S1473-3099(20)30120-1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;close all;clc;
load('ItalianRegion_intervention_fitting.mat')
%% Part 1. sample input parameters, independent inputs
% number of Monte Carlo run, set to 100 for a test, in the paper n = 1000000
n = 1000
% [x1,..., x6] = [\alpha \beta 1/\gamma, delta, I0, Intervention time]
x = nan(n,6); 
rng(1234)
x(:,1) = (alpha1*1.1 - alpha1*0.9).*rand(n,1) + alpha1*0.9; % alpha
x(:,2) = (beta1*1.1 - beta1*0.9).*rand(n,1) + beta1*0.9; % beta
x(:,3) = ((1/gamma1)*1.1 - (1/gamma1)*0.9).*rand(n,1) + (1/gamma1)*0.9; %1/gamma
x(:,4) = (delta1*1.3 - delta1*0.7).*rand(n,1) + delta1*0.7; %delta1
I0 = Confirmed(1);
x(:,5) =  randsample([round(I0*0.8):round(I0*1.2)],n,true); %(I0*1.2 - I0*0.8).*rand(n,1) + I0*0.8;
x(:,6) = randsample([nwi:nwi+7],n,true);  % intervention time
%% run simulation
y = nan(n,1);
tic
for i = 1:n
[y(i)] =  SEIQRDP_intervention_Italy(x(i,:));
end
toc

%% save data
%save(['SEIRdata20200420_Italy.mat'],'x','y')


%% Part 2. sample input parameters with correlations: corr(\alpha,\delta) = 0.5
clearvars
rho = 0.5;
% number of Monte Carlo run, set to 100 for a test, in the paper n = 1000000
n = 1000
% [x1,..., x6] = [\alpha \beta 1/\gamma, delta, I0, Intervention time]
x = nan(n,6);
rng(1234)
Z = mvnrnd([0 0],[1 rho; rho 1], n);
U = normcdf(Z);
alpha = U(:,1); delta = U(:,2);

x(:,1) = (alpha1*1.1 - alpha1*0.9).*alpha + alpha1*0.9; % alpha
x(:,2) = (beta1*1.1 - beta1*0.9).*rand(n,1) + beta1*0.9; % beta
x(:,3) = ((1/gamma1)*1.1 - (1/gamma1)*0.9).*rand(n,1) + (1/gamma1)*0.9; %1/gamma
x(:,4) = (delta1*1.3 - delta1*0.7).*delta + delta1*0.7; %delta1
I0 = Confirmed(1);
x(:,5) =  randsample([round(I0*0.8):round(I0*1.2)],n,true);
x(:,6) = randsample([nwi:nwi+7],n,true); 

%%
y = nan(n,1);

tic
for i = 1:n
[y(i)] =  SEIQRDP_intervention_Italy(x(i,:));
end
toc

%% save
%save(['SEIRdata20200420_Italy_correlation.mat'],'x','y')
