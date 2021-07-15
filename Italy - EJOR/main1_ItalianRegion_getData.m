% * Codes for collecting COVID-2019 data for Italy
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


%% Initialisation
clearvars;close all;clc;

% Download the data from ref [1] and read them with the function
% getDataCOVID_ITA
tableCOVIDItaly = getDataCOVID_ITA();

time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));
fprintf(['Most recent update: ',datestr(time(end)),'\n'])


%% Cumulative data
% Perform nation-wide statistics by summing the data of all the regions.

% Merge regional data for each day
tableCOVIDItaly_Tot = varfun(@sum,tableCOVIDItaly, ...
    'InputVariables',tableCOVIDItaly.Properties.VariableNames(7:end), ...
    'GroupingVariables','Date');
% Remove the 'GroupCount' variable, should total to the number of Italian regions (19 + 2 autonomous provinces)
tableCOVIDItaly_Tot = removevars(tableCOVIDItaly_Tot,'GroupCount');
% Rename the accumulated variables with the original variable names
tableCOVIDItaly_Tot.Properties.VariableNames=[tableCOVIDItaly.Properties.VariableNames(1),tableCOVIDItaly.Properties.VariableNames(7:end)];

% population
Npop = 60.48e6; 

Recovered   = tableCOVIDItaly_Tot.Recovered'  ;
Deaths      = tableCOVIDItaly_Tot.Deaths'     ;
Confirmed   = tableCOVIDItaly_Tot.Confirmed'  ;
Quarantined = tableCOVIDItaly_Tot.Quarantined';
time        = tableCOVIDItaly_Tot.Date'       ;
time = unique(datetime(datestr(datenum(tableCOVIDItaly.Date,'yyyy-mm-DDThh:MM:ss'))));
time.Format = 'dd-MMM-yyyy'; 

%% focus on data from 24-Feb-2020 to 20-Apr-2020

ind = find(time >= '24-Feb-2020' & time < '21-Apr-2020' );
time = time(ind);
Recovered = Recovered(ind);
Deaths  = Deaths(ind);
Quarantined =Quarantined(ind);
Confirmed =Confirmed(ind) ;

%% save data
% save('ItalainRegionData0420.mat',...
%        'Npop','time',...
%     'Confirmed','Quarantined', 'Recovered','Deaths')

