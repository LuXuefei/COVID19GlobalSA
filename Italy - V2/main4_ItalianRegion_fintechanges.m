% * Codes for generate sample of locally finite change decomposition 
%  'Uncertainty in OR Epidemiological Modeling: A Global Sensitivity Approach', 
%   by Xuefei Lu and Emanuele Borgonovo, 2021
% 
% * Author: Xuefei Lu, xuefei.lu@ed.ac.uk
% * Date: July, 2021
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;close all;clc;
load('ItalianRegion_intervention_fitting.mat')
%% Part 1. sample input parameters, independent inputs
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
pool = [round(I0*0.8):round(I0*1.2)]; 
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

%% Part 2. sample input parameters with correlations: corr(\alpha,\delta) = 0.5
clearvars
% number of Monte Carlo run, set to 10 for a test, in the paper n = 20000;
n = 10
% [x1,..., x6] = [\alpha \beta 1/\gamma, delta, I0, Intervention time]
x = nan(n,6);

rho = 0.5
rng(1234)
Z = mvnrnd([0 0],[1 rho; rho 1], n);
U = normcdf(Z);
alpha = U(:,1); delta = U(:,2);
x(:,1) = (alpha1*1.1 - alpha1*0.9).*alpha + alpha1*0.9; % alpha
x(:,2) = (beta1*1.1 - beta1*0.9).*rand(n,1) + beta1*0.9; % beta
x(:,3) = ((1/gamma1)*1.1 - (1/gamma1)*0.9).*rand(n,1) + (1/gamma1)*0.9; %1/gamma
x(:,4) = (delta1*1.3 - delta1*0.7).*delta + delta1*0.7; %delta1
I0 = Confirmed(1);
pool = [round(I0*0.8):round(I0*1.2)]; 
x(1,5) =  randsample(pool,1,true);
for i = 1:(n-1)
    x(i+1,5) =  randsample(setdiff(pool,x(i,5)),1,true);
end

pool = [nwi:nwi+7]; 
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
%save(['InteractionsItaly_corr.mat'],'ffmc','ymc','phimc','U','x')