close all; clear all; clc;
fsize=18; fname='times';

% cd 'C:\Users\Abhishek\Work\Lecturer_CU\Program\cooling_PM\Prishat_code\'
%cd '/home/abhishek/work/Program/cooling_PM/Prishat_code'
flpath =  'C:\Users\Prishat\Desktop\Project-1\Data\';   % alternatively, give full path

% prevent showing any figure
%set(0,'DefaultFigureVisible','off')


% read all file
% dat = read_dat_fun('ALL',flpath);



flnm = 'Q=8_Ti=600_Du=125_Dl=60_data.xlsx';
% Create a spreadsheet DataMod.xlsx the data with the following columns
% average plate temperature at current time step
% simulated coefficient of convective heat loss
% simulated ambient temperature
% volume flow rate (Q) of steam
% cooling rate


flnmpth = strcat(flpath, flnm);
num = xlsread(flnmpth);






%%%% Gaussian Process regression module %%%%
xdat = num(:,2:5);
ydat = num(:,6);

gprMdlfull_BM_G = fitrgp(xdat,ydat,'KernelFunction','squaredexponential',...
    'OptimizeHyperparameters','auto','HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus','Optimizer','randomsearch'));

%ypred1 = resubPredict(gprMdlfull_BM_G);
[ypred,~,yint] = resubPredict(gprMdlfull_BM_G,'Alpha',0.01);

figure;
h1 = area([yint(:,1) yint(:,2)-yint(:,1)],-8,...
'FaceColor',[0.85,0.85,0.85],'EdgeColor',[0.85,0.85,0.85]);
hold on;
h1(1).FaceColor = 'none'; % remove color from bottom area
h1(1).EdgeColor = 'none';
h2 = plot(ydat,'r');hold on;% Plot original response values
h3 = plot(ypred,'b--'); % Plot predicted response values
%axis([0 510 -7 65]);
hold off

