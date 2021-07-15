function [y,time1,Totalnr,Qt,Rt,Dt] = SEIQRDP_intervention_Italy(x,gfx)
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


% inputs should be x = [alpha1, beta1 , 1/gamma1, delta1, I0, tinv]
alpha1 =  x(1);
beta1 = x(2);
gamma1 = 1/x(3);
%delta1 = 1/x(4);
delta1 = x(4);
I0 = x(5);
tinv= x(6);



load('ItalianRegion_intervention_fitting.mat', 'alphaw', 'betaw', 'deltaw','gammaw',...
    'Npop','Lambdaw','Kappaw','Lambda1','Kappa1','time','Recovered','Quarantined','Deaths')
dt = 0.1;
time1 = time(1):dt:time(end); % after intervention
N = numel(time1);
tt = [0:N-1].*dt;

if(nargin<2)
    gfx='';
end

% default inputs
t0 = tt(tt <= tinv); % before intervention
t1 = tt(tt >= t0(end)); % after intervention

%% simulate days before intervention, parameters estimated from
% 05-Mar-2020 to 09-Mar-2020

% 
% Es0 = I0;Is0 = I0;
% Qs0 = 0;Rs0=0;Ds0=0;


Es0 = I0; % Initial number of exposed cases. Unknown but unlikely to be zero.
Is0 = I0; % Initial number of infectious cases. Unknown but unlikely to be zero.
Qs0 = Quarantined(1);
Rs0 = Recovered  (1);
Ds0 = Deaths     (1);


[Ss,Es,Is,Qs,Rs,Ds,Ps] = SEIQRDP(alphaw,betaw,gammaw,deltaw,Lambdaw,Kappaw,...
    Npop,Es0,Is0,Qs0,Rs0,Ds0,t0);
%% simulate days after intervention

E1 = Es(end); 
I1 = Is(end); 
Q1 = Qs(end);
R1 = Rs(end);
D1 = Ds(end);

[Sp,Ep,Ip,Qp,Rp,Dp,Pp] = SEIQRDP(alpha1,beta1,gamma1,delta1,Lambda1,Kappa1,...
    Npop,E1,I1,Q1,R1,D1,t1 - tinv); %start as last of tinv
%% Adjust tinv
t1 = t1(2:end);
Sp = Sp(2:end);
Ep = Ep(2:end);
Ip = Ip(2:end);
Qp = Qp(2:end);
Rp = Rp(2:end);
Dp = Dp(2:end);
Pp = Pp(2:end);
% entirer series
St = [Ss,Sp];
Et = [Es,Ep];
It = [Is,Ip];
Qt = [Qs,Qp];
Rt = [Rs,Rp];
Dt = [Ds,Dp];
Pt = [Ps,Pp];

%%
Totalnr = Qt + Rt + Dt;
y = Totalnr(end);

%%
if ~isempty(gfx)
figure
semilogy(time1,[Qt+Rt+Dt],'c',time1,Qt,'r',time1,Rt,'b',time1,Dt,'k'); hold on
xline(time(floor(tinv)) + (tinv -floor(tinv) ) );
ylabel('Number of cases')
xlabel('time (days)')
title('Italy, intervention');
% leg = {'Total','Quarantined (confirmed infectious)','Recovered','Dead'};
% legend(leg{:},'location','southoutside')
set(gcf,'color','w')
set(gca,'yscale','lin')
end

