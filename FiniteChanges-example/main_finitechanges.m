% * Codes for compute the first order and higher order effects \phi_ijk at random
%  locations
%  Emanuele Borgonovo, Xuefei Lu(2020)
%  Is the COVID-19 Outbreak Sensitive to the Time of Intervention? Evidence from Italy
% 
% * Author:Emanuele Borgonovo and Xuefei Lu, xuefei.lu@unibocconi.it
% * Date: April, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  An example is made on 
% * The 21-input Additive Gassian simulator
%                             $Y=X1 * X2 + X3 $
%  where X_i ~ Uniform(1; 2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;clc;close all
%% generate dataset
n = 100;
x = nan(n,3); % alpha beta delta^-1 E0+I0
rng(1234)
x(:,1) = rand(n,1) + 1;
x(:,2) = rand(n,1) + 1;
x(:,3) = rand(n,1) + 1;
%% calculate effects at random sample
k = size(x,2);
ymc = nan(n-1,2^k);ffmc = nan(n-1,2^k);phimc = nan(n-1,k);

tic
for i = 1:(n-1)
x0 = x(i,:);
x1 = x(i+1,:);
[k,U,DX,y,ff,phi]=finitechanges(x0,x1,@funcg); % change the funcg as the desired one
ymc(i,:) = y;
ffmc(i,:) = ff;
phimc(i,:) = phi;
end
toc
%% delete replicates if necessary, especially for discrete variables
xx = x; xx(end,:)=[];
ffmcc = ffmc;
n= size(x,1);
for i = 1:n-1
x0 = x(i,:);
x1 = x(i+1,:);
if sum(x0 == x1)
xx(i,:)= nan(1,6);
ffmcc(i,:) = nan(1,size(ffmc,2));
end
end
indvalid = find (~isnan(xx(:,1)) == 1);
%% Mean effects
absffmc=abs(ffmc(indvalid,:)); 
figure
ordernr = sum(U);
bar(1:length(mean(absffmc)),mean(absffmc))
set(gca, 'XTick',  1:length(mean(absffmc)))
set(gca,'XTickLabel',ordernr )
ylabel('Mean of \phi_{i,j,k...}','Fontsize',16)
xlabel('Interactions in increasing order','Fontsize',16)
%% higher order effects direction
[~,ind] = sort(mean(absffmc),'descend');
indp = ind(ordernr(ind(1:8))>1); %ipnames(find(U(:,indp(1)) == 1))

% \frac{phi_ij}{\triangle(xi)*\triangel(xj)}
nn= length(indvalid);
ffphipartial = nan(nn,4);
for ii = 1:4
for i = 1:nn
x0 = x(indvalid(i),:);
x1 = x(indvalid(i)+1,:);
ff = ffmc(indvalid(i),indp(ii));
xd = x1 - x0;
den = xd .* [U(:,indp(ii)) == 1]';
den(~den) =1;
ffphipartial(i,ii) = ff/prod(den);
end
end

% top two higher order interactions
figure
for ii = 1:2
subplot(1,2,ii)
%hist(ffmc(:,indp(ii))) %\phi_ij
histogram(ffphipartial(:,ii),50)
ylabel(['\phi_{',num2str(find(U(:,indp(ii)) == 1)'),'}'],'Fontsize',16)
xlabel(num2str(find(U(:,indp(ii)) == 1)'),'Fontsize',16)
end

%%
AA=cov(ffmc(indvalid,1:size(x,2)))/var(ymc(indvalid,2^size(x,2)));
diagT=diag(AA);
% main effects
T=diagT/2;
% mean dimension
meandsimension=sum(T) %1.1679