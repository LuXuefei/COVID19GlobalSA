clearvars;close all;clc;
load('ItalainRegionData0420.mat')

%% Step 1. Estimate before-intervation parameters from 05-Mar-2020 to 15-Mar-2020
semilogy(time,Confirmed,'co',time,Quarantined,'ro',time,Recovered,'bo',time,Deaths,'ko');
%%
nwi=14;%11;
timewi = time(1:nwi);
% Initial conditions
Ew0 = Confirmed  (1); % Initial number of exposed cases. Unknown but unlikely to be zero.
Iw0 = Confirmed  (1); % Initial number of infectious cases. Unknown but unlikely to be zero.
Qw0 = Quarantined(1);
Rw0 = Recovered  (1);
Dw0 = Deaths     (1);

% Definition of the first estimates for the parameters
alpha_guess = 0.001;%0.06;
beta_guess = 1; %0.8; % Infection rate
LT_guess = 5; % latent time in days
Q_guess = 21; % rate at which infectious people enter in quarantine
lambda_guess = [0.1,0.03]; % recovery rate
kappa_guess = [0.07,0.03]; % death rate
% 
guess = [alpha_guess,...
    beta_guess,...
    1/LT_guess,...
    1/Q_guess,...
    lambda_guess,...
    kappa_guess];

rng(1234)
[alphaw,betaw,gammaw,deltaw,Lambdaw,Kappaw,residuals] = ...
    fit_SEIQRDP(Quarantined(1:nwi),Recovered(1:nwi),Deaths(1:nwi),Npop,Ew0,Iw0,timewi,guess,'Display','off');
r2 =@(y,r) 1- sum(r(:).^2)/sum((y(:) - mean(y)).^2);
rrsew(1) = r2(Quarantined(1:nwi),residuals(1,:));
rrsew(2) = r2(Recovered(1:nwi),residuals(2,:));
rrsew(3) = r2(Deaths(1:nwi),residuals(3,:));
mean(rrsew) %0.9542
%% Step 2. Estimate after-intervation parameters from 16-Mar-2020 to 07-Apr-2020

E0 = Confirmed  (nwi+1); % Initial number of exposed cases. Unknown but unlikely to be zero.
I0 = Confirmed  (nwi+1); % Initial number of infectious cases. Unknown but unlikely to be zero.
Q0 = Quarantined(nwi+1);
R0 = Recovered  (nwi+1);
D0 = Deaths     (nwi+1);

% Definition of the first estimates for the parameters
alpha_guess = 0.06; % protection rate
beta_guess = 1.0; % Infection rate
LT_guess = 5; % latent time in days
QT_guess = 21; % quarantine time in days
lambda_guess = [0.1,0.05]; % recovery rate
kappa_guess = [0.1,0.05]; % death rate
% 
guess = [alpha_guess,...
    beta_guess,...
    1/LT_guess,...
    1/Q_guess,...
    lambda_guess,...
    kappa_guess];

rng(1234)
[alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,residuals] = ...
    fit_SEIQRDP(Quarantined((nwi+1):end),Recovered((nwi+1):end),Deaths((nwi+1):end),Npop,E0,I0,time((nwi+1):end),guess,'Display','off');
r2 =@(y,r) 1- sum(r(:).^2)/sum((y(:) - mean(y)).^2);
rrse1(1) = r2(Quarantined((nwi+1):end),residuals(1,:));
rrse1(2) = r2(Recovered((nwi+1):end),residuals(2,:));
rrse1(3) = r2(Deaths((nwi+1):end),residuals(3,:));
mean(rrse1) % 0.9980
days(time(end) - time(nwi+1))
%% To simulate the cases after fitting
dt = 0.1; % time step
timew = time(1):dt:time(nwi); % before intervention
N = numel(timew);
tw = [0:N-1].*dt;
[Sw,Ew,Iw,Qw,Rw,Dw,Pw] = SEIQRDP(alphaw,betaw,gammaw,deltaw,Lambdaw,Kappaw,Npop,Ew0,Iw0,Qw0,Rw0,Dw0,tw);

%
figure
col = get(gca, 'ColorOrder'); hold on;
plot(time(1:nwi),Confirmed(1:nwi),'o','Color', col(6,:),'LineWidth',1.5)
plot(time(1:nwi),Quarantined(1:nwi),'o','Color', col(2,:),'LineWidth',1.5)
plot(time(1:nwi),Recovered(1:nwi),'o','Color', col(1,:),'LineWidth',1.5)
plot(time(1:nwi),Deaths(1:nwi),'o','Color', col(4,:),'LineWidth',1.5)
plot(timew,Qw+Rw+Dw,'-','Color', [0.5,0.5,0.5],'LineWidth',1.5)
plot(timew,Qw+Rw+Dw,'-','Color', col(6,:),'LineWidth',1.5)
plot(timew,Qw,'-','Color', col(2,:),'LineWidth',1.5)
plot(timew,Rw,'-','Color', col(1,:),'LineWidth',1.5)
plot(timew,Dw,'-','Color', col(4,:),'LineWidth',1.5)
legend({'Total','Quarantined','Recovered','Deceased','Fitted'},'Location', 'best')
title('Italy, fitting parameters before intervention');
ylabel('Number of cases')
xlabel('time(days)')
%set(gca,'yscale','lin')


% figure
% semilogy(timew,Qw+Rw+Dw,'c',timew,Qw,'r',timew,Rw,'b',timew,Dw,'k');hold on
% semilogy(time(1:nwi),Confirmed(1:nwi),'co',time(1:nwi),Quarantined(1:nwi),'ro',time(1:nwi),Recovered(1:nwi),'bo',time(1:nwi),Deaths(1:nwi),'ko');
% ylabel('Number of cases')
% xlabel('time (days)')
% title('Italy, fitting parameters before intervention');
% % leg = {'Total','Quarantined (confirmed infectious)','Recovered','Dead'};
% % legend(leg{:},'location','southoutside')
% set(gcf,'color','w')
% set(gca,'yscale','lin')
set(gcf, 'PaperPosition', [-0.5 0 15 8]); %[left bottom width height]
set(gcf, 'PaperSize', [14 8]);
saveas(gca,'Italy_intervention_fitting_beforeInt','pdf')

%%
dt = 0.1;
time1 = time(nwi+1):dt:time(end); % after intervention
N = numel(time1);
t = [0:N-1].*dt;
[S,E,I,Q,R,D,P] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,Npop,E0,I0,Q0,R0,D0,t);


figure
col = get(gca, 'ColorOrder'); hold on;
plot(time((nwi+1):end),Confirmed((nwi+1):end),'o','Color', col(6,:),'LineWidth',1.5)
plot(time((nwi+1):end),Quarantined((nwi+1):end),'o','Color', col(2,:),'LineWidth',1.5)
plot(time((nwi+1):end),Recovered((nwi+1):end),'o','Color', col(1,:),'LineWidth',1.5)
plot(time((nwi+1):end),Deaths((nwi+1):end),'o','Color', col(4,:),'LineWidth',1.5)
plot(time1,Q+R+D,'-','Color', [0.5,0.5,0.5],'LineWidth',1.5)
plot(time1,Q+R+D,'-','Color', col(6,:),'LineWidth',1.5)
plot(time1,Q,'-','Color', col(2,:),'LineWidth',1.5)
plot(time1,R,'-','Color', col(1,:),'LineWidth',1.5)
plot(time1,D,'-','Color', col(4,:),'LineWidth',1.5)
legend({'Total','Quarantined','Recovered','Deceased','Fitted'},'Location', 'best')
title('Italy, fitting parameters after intervention');
ylabel('Number of cases')
xlabel('time(days)')
%set(gca,'yscale','lin')

% figure
% semilogy(time1,Q+R+D,'c',time1,Q,'r',time1,R,'b',time1,D,'k'); hold on
% semilogy(time((nwi+1):end),Confirmed((nwi+1):end),'co',time((nwi+1):end),Quarantined((nwi+1):end),'ro',time((nwi+1):end),Recovered((nwi+1):end),'bo',time((nwi+1):end),Deaths((nwi+1):end),'ko');
% ylabel('Number of cases')
% xlabel('time (days)')
% title('Italy, fitting parameters after intervention');
% % leg = {'Total','Quarantined (confirmed infectious)','Recovered','Dead'};
% % legend(leg{:},'location','southoutside')
% set(gcf,'color','w')
% set(gca,'yscale','lin')
set(gcf, 'PaperPosition', [-0.5 0 15 8]); %[left bottom width height]
set(gcf, 'PaperSize', [14 8]);
saveas(gca,'Italy_intervention_fitting_AfterInt','pdf')

%%
save('ItalianRegion_intervention_fitting.mat',...
    'alphaw','betaw','gammaw','deltaw','Lambdaw','Kappaw','Npop',...
    'alpha1','beta1','gamma1','delta1','Lambda1','Kappa1',...
    'time','nwi',...
    'Confirmed','Quarantined','Recovered','Deaths')
