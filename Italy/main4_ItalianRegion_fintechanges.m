% * Codes for generate sample of locally finite change decomposition 
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
% number of Monte Carlo run, set to 10 for a test, in the paper n = 20000;
n = 10
% [x1,..., x6] = [\alpha \beta 1/\gamma, delta, I0, Intervention time]
x = nan(n,6);
rng(1234)
x(:,1) = (alpha1*1.1 - alpha1*0.9).*rand(n,1) + alpha1*0.9; % alpha
x(:,2) = (beta1*1.1 - beta1*0.9).*rand(n,1) + beta1*0.9; % beta
x(:,3) = ((1/gamma1)*1.1 - (1/gamma1)*0.9).*rand(n,1) + (1/gamma1)*0.9; %1/gamma
x(:,4) = (delta1*1.3 - delta1*0.7).*rand(n,1) + delta1*0.7; %delta
I0 = Confirmed(1);
% For discrete variables, avoid consecutive duplicates
pool = [round(I0*0.8):round(I0*1.2)]; %(I0*1.2 - I0*0.8).*rand(n,1) + I0*0.8;
x(1,5) =  randsample(pool,1,true);
for i = 1:(n-1)
    x(i+1,5) =  randsample(setdiff(pool,x(i,5)),1,true);
end

pool = [nwi:nwi+7]; % intervention time
x(1,6) =  randsample(pool,1,true);
for i = 1:(n-1)
    x(i+1,6) =  randsample(setdiff(pool,x(i,6)),1,true);
end


%% run simulation
k = size(x,2);
ymc = nan(n-1,2^k);ffmc = nan(n-1,2^k);phimc = nan(n-1,k);

tic
for i = 1:(n-1)
%input = [ alpha1, beta, 1/gamma,delta, I0, tinv]
x0 = x(i,:);
x1 = x(i+1,:); 
[k,U,DX,y,ff,phi]=finitechanges(x0,x1,@SEIQRDP_intervention_Italy);
ymc(i,:) = y;
ffmc(i,:) = ff;
phimc(i,:) = phi;
end
toc
%% save data
%save(['InteractionsItaly.mat'],'ffmc','ymc','phimc','U','x')
