% * Codes for Fig.1 (a)-(b) generate Monte Carlo sample using the fitted model
%  'Are Interventions in the COVID-19 Outbreak Really Important? A Global Sensitivity Approach', 
%   by Xuefei Lu and Emanuele Borgonovo, 2020
% 
% * Author: Xuefei Lu, xuefei.lu@ed.ac.uk
% * Date: Dec, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Cheynet, E. Generalized SEIR Epidemic Model (Fitting and Computation). Zenodo, 2020, doi:10.5281/ZENODO.3911854. 
% [2] Dong E, Du H, Gardner L. An interactive web-based dashboard to track COVID-19 in real time. Lancet Inf Dis. 
%     20(5):533-534. doi: 10.1016/S1473-3099(20)30120-1
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;close all;clc;
load('ItalianRegion_intervention_fitting.mat')
%% sample input parameters
% number of Monte Carlo run, set to 100 for a test, in the paper n = 1000000
n = 100
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
