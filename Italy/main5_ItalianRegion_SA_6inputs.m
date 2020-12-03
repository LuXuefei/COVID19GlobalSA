% * Codes for Sensitivity Analysis of SEIR model
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
%% define input parameter names
ipnames = {'\alpha', '\beta', '\gamma^{-1}', '\delta', 'I_0', 'Interv. in March'};

%% Fig.2. importance measures
% run main3_ItalianRegion_intervention_datagenerator.m to obtain the sample below
load('SEIRdata20200420_Italy')
% may take few seconds
[b,d,t,e,w]=betaKS3(x,y,50,'CDFs');

% run main4_ItalianRegion_fintechanges.m to obtain the sample below
load('InteractionsItaly.mat')
%
xx = x; xx(end,:)=[];
ffmcc = ffmc;
n= size(x,1);
for i = 1:n-1
x0 = [x(i,1), x(i,2),x(i,3), x(i,4), x(i,5),x(i,6)];
x1 = [x(i+1,1), x(i+1,2),x(i+1,3), x(i+1,4), x(i+1,5),x(i+1,6)];
if sum(x0 == x1)
xx(i,:)= nan(1,6);
ffmcc(i,:) = nan(1,size(ffmc,2));
end
end
indvalid = find (~isnan(xx(:,1)) == 1);

AA=cov(ffmc(indvalid,1:6))/var(ymc(indvalid,64));
diagT=diag(AA);
T=diagT/2;
meandsimension=sum(T) %1.1683

%plot
%define colors
Col = [141,211,199; 255,255,179;
190,186,218;
251,128,114;
128,177,211;
179,222,105;
253,180,98]/255;

% Sensitivity measures used: [Total Sobol; first order Sobol, Kuiper index]
AAA=[T'; e; t];
figure
bp = bar(AAA); %,'EdgeColor','none'
for ii = 1:length(AAA)
bp(ii).FaceColor = Col(ii,:);%repmat(Col(ii,:),size(AAA,1),1);
end
xlabel('Sensitivity Measure','Fontsize',13)
set(gca,'XTickLabel',{'T_i','S_i','\beta_{i}^{KU}'},'FontSize',13);


%% Fig. 4. The portions of the output variance

figure
Col(7,:) = [0,0,0]
eplot = [e,1-sum(e)];
labels = {ipnames{:},'Interactions'};
ax = gca();
p = pie(ax,eplot,[zeros(1,6),1]);
ax.Colormap = Col; 
set(p,'EdgeColor','none')%,'LineStyle','none')

%% Fig. 5. Interaction effects
plotk=25; % number of bars in the first panel
absffmc=abs(ffmc(indvalid,:));
ordernr = sum(U);
ffmcplot = mean(absffmc); %meanabsffmc
ffmcplot2 = ffmcplot(7:end);
ordernr2 = ordernr(7:end);
U2 = U(:,7:end);
[~,ind] = sort(ffmcplot2,'descend');

figure
subplot(1,5,[1 2])
hB = bar(1:plotk,ffmcplot2(ind(1:plotk)),'FaceColor',[254 224 139]/255,'EdgeColor','none')
for ii = 1:3
text(ii,hB.YData(ii), ['{',strjoin(ipnames(find(U2(:,ind(ii)) == 1)),','),'}'], ...
          'VerticalAlignment','bottom', 'Fontsize',13)
end
yl = get(gca,'Ylim');
yl(2) = yl(2)*1.1
ylim(yl)
set(gca, 'XTick',  1:plotk)
set(gca,'XTickLabel',ordernr2(ind(1:plotk)) )
ylabel('Mean magnitudes of interactions','Fontsize',13)
xlabel('Order of interaction','Fontsize',13)

% sign of interactions
[~,ind] = sort(mean(absffmc),'descend');
indp = ind(ordernr(ind(1:20))>1);
meanabsffmc = mean(absffmc);
nn= length(indvalid);
ffphipartial = nan(nn,4);
for ii = 1:4
for i = 1:nn
%input = [ alpha1, beta, 1/gamma,delta, I0, tinv]
x0 = x(indvalid(i),:);
x1 = x(indvalid(i)+1,:);
ff = ffmc(indvalid(i),indp(ii));
xd = x1 - x0;
den = xd .* [U(:,indp(ii)) == 1]';
den(~den) =1;
ffphipartial(i,ii) = ff/prod(den); %--> partial derivartives % \frac{phi_ij}{\triangle(xi)*\triangel(xj)}
end
end

for ii = 1:3
subplot(1,5,ii+2)
if mean(ffphipartial(:,ii))<0
    histogram(ffphipartial(:,ii),30,'FaceColor',[102 189 99]/255)
else
   histogram(ffphipartial(:,ii),30,'FaceColor',[252 141 89]/255)
end
ylabel(['Interaction ', [strjoin(ipnames(find(U(:,indp(ii)) == 1)),' and ')]], 'Fontsize',13)
xlabel([strjoin(ipnames(find(U(:,indp(ii)) == 1)),' , ')],'Fontsize',13)
end