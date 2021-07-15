% * Codes for Fig.3 (a)-(d) generate Monte Carlo sample using the fitted model
%  'Uncertainty in OR Epidemiological Modeling: A Global Sensitivity Approach', 
%   by Xuefei Lu and Emanuele Borgonovo, 2021
% 
% * Author: Xuefei Lu, xuefei.lu@ed.ac.uk
% * Date: July, 2021
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clearvars
close all
load('ItalianRegion_intervention_fitting.mat')
%% Load data
% Independent inputs
load('SEIRdata20200420_Italy.mat','x','y');
% Correlated inputs
% load('SEIRdata20200420_Italy_correlation.mat','x','y');
sname = 'Italy'
%% Step 1. generate trajactories with input uncertainty (might take a few minutes)
n = 10000;
rng(1234)
ind = randsample(size(x,1),n,false);
x = x(ind,:);
%% 
i=1;
[yi,time1i,Totalnri,Qti,Rti,Dti] =  SEIQRDP_intervention_Italy(x(i,:));
len = length(time1i);
y = nan(n,1); %timelm = size(n,len); %same
Totalnrm = nan(n,len); Qtm = nan(n,len);Rtm = nan(n,len);Dtm = nan(n,len);
tic
for i = 1:n
[y(i),time1,Totalnr,Qt,Rt,Dt] =  SEIQRDP_intervention_Italy(x(i,:));
Totalnrm(i,:) = Totalnr;
Qtm(i,:) = Qt;
Rtm(i,:) = Rt;
Dtm(i,:) = Dt;
end
toc

%% Step 2. Plot
%% Plot Fig 3 (a) or (b): uncertainty quantification
yplot = prctile(Totalnrm, [2.5,97.5],1);
figure
col = get(gca, 'ColorOrder'); coli = 6;
plot(time,Confirmed,'-o','Color', col(1,:),'LineWidth',1.5);hold on
shade(time1i,yplot(2,:),time1i,yplot(1,:), 'Color', col(coli,:), ...
    'FillType',[1 2], 'FillColor',col(coli,:), 'FillAlpha',0.3)
plot(time,Confirmed,'-o','Color', col(1,:),'LineWidth',1.5)
xline(time(nwi));
text(time(floor(nwi*0.2)),0.6*max(get(gca,'Ylim')),'Before Interv.')
text(time(nwi + 10),0.6*max(get(gca,'Ylim')),'After Interv.')
pe = xline(time(end),'Color', col(coli,:),'LineWidth',4);
pe.Color(4) = 0.01;
ylim([0,430000])
leg = {'Real data','95%PI'};
legend(leg{:},'location','North')
title(['Cumulative total infections, ',sname])

%% Plot Fig (b) or (d): uncertainty quantification, histogram
%% Load data
clearvars;
load('ItalainRegionData0420.mat')
% Independent inputs
load('SEIRdata20200420_Italy.mat','x','y');
% Correlated inputs
% load('SEIRdata20200420_Italy_correlation.mat','x','y');
sname = 'Italy';
%%
yy = y;
figure
col = get(gca, 'ColorOrder'); 
h = histogram(y,100,'Facecolor',col(6,:),'EdgeColor','none','Normalization','pdf');
xline(prctile(yy,2.5),'Color', col(6,:),'LineWidth',3);
xline(prctile(yy,97.5),'Color', col(6,:),'LineWidth',3);
xline(Confirmed(end),'Color', col(1,:),'LineWidth',3);

anArrow = annotation('doublearrow','Linewidth',2, 'Color',col(6,:)) ;
anArrow.Parent = gca;  % or any other existing axes or figure
ypmax = max(get(gca,'Ylim'));
anArrow.Position = [prctile(yy,2.5), 0.6*ypmax,prctile(yy,97.5)-prctile(yy,2.5) , 0] ;
text(mean(yy), 0.63*ypmax, '95% PI','Fontsize',13,'Color',col(6,:))

text(Confirmed(end),0.8*ypmax,'Real','Color',col(1,:),'Fontsize',13)
xpmax = max(get(gca,'Xlim'));
textstat = {['Real: ', num2str(Confirmed(end),'%.2s')],...
    ['Pred. Median: ', num2str(median(yy),'%.2u')] ,...
    ['Pred. Mean: ', num2str(mean(yy),'%.2u')],...
    ['Pred. Std: ', num2str(std(yy),'%.2u')]};
text(0.55*xpmax, 0.73*ypmax, textstat,'Fontsize',13)

title(['Total number of infected people (Q+R+D) on ',datestr(time(end),'dd-mmm'),',',sname])