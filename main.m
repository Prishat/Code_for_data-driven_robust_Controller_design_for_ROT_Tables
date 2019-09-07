close all; clear all; clc;
fsize=18; fname='times';

% cd 'C:\Users\Abhishek\Work\Lecturer_CU\Program\cooling_PM\Prishat_code\'
%cd '/home/abhishek/work/Program/cooling_PM/Prishat_code'
%flpath = '../newDataWithAmbTemp/';   % alternatively, give full path
% prevent showing any figure
%set(0,'DefaultFigureVisible','off')


% read all file
% dat = read_dat_fun('ALL',flpath);

flpath = 'C:\Users\Prishat\Desktop\project code\matlab\ti=600Q(var)';

flnm = 'Q= 8 Ti=600 Du=125, Dl=60.xlsx';
flnmpth = strcat(flpath, flnm);
flnum = 'Q= 8,Ti=600,Du=125,Dl=60';
num = xlsread('C:\Users\Prishat\Desktop\project code\matlab\ti=600Q(var)\Q= 8,Ti=600,Du=125,Dl=60');
tsteps = num(:,1); Ta = num(:,4); T = num(:,8:12);
nsen = size(T,2);
Tavg = sum(T,2)./size(T,2);


%%% CHANGING Q PLOTS %%%
colline = {'.b','.r','.y','.g','.m'};
flseq = {'1','2','3','4','5'}; flseq{nsen+1} = 'average';
for ii=1:size(T,2)
    plot(tsteps, T(:,ii), colline{ii},'LineWidth',1.5); hold on;
end
plot(tsteps, Tavg, '-k', 'LineWidth', 4)
hold off;
ylabel('average temperature (^oC)','FontName',fname,'fontsize',fsize);
xlabel('time (s)','FontName',fname,'fontsize',fsize);
% xlim([0 700]); ylim([50 650])
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
legend(flseq);
%eval(['print -dpng ./figs/T_all_', num2str(flnum),'.png' ]);
%eval(['print -depsc ./figs/T_all_', num2str(flnum),'.eps' ]);


% CALCULATE COOLING RATE
njump = 100;
dt = tsteps(njump+1:njump:end) - tsteps(1:njump:end-njump);
dT = Tavg(njump+1:njump:end) - Tavg(1:njump:end-njump);
dTdt = dT./dt;


figure()
plot(Tavg(1:njump:end-njump),dTdt,'LineWidth',2.5); hold on;
ylabel('cooling rate (^oC/s)','FontName',fname,'fontsize',fsize);
xlabel('average temperature (^oC)','FontName',fname,'fontsize',fsize);
% xlim([0 400]); ylim([-30 30])
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
% create a folder named figs
%eval(['print -dpng ./figs/theta_all', num2str(flnum),'.png' ]);
%eval(['print -depsc ./figs/theta_all', num2str(flnum),'.eps' ]);


% CALCULATE CONVECTIVE HEAT LOSS COEFFICIENT
c = dTdt./(Ta(njump+1:njump:end)-Tavg(njump+1:njump:end));

figure()
scatter3(Tavg(1:njump:end-njump),Ta(njump+1:njump:end),c); hold on;
zlabel('convective coeff','FontName',fname,'fontsize',fsize);
ylabel('ambient T (^oC)','FontName',fname,'fontsize',fsize);
xlabel('avg T (^oC)','FontName',fname,'fontsize',fsize);
% xlim([0 400]); ylim([-30 30])
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
% create a folder named figs
%eval(['print -dpng ./figs/theta_all', num2str(flnum),'.png' ]);
%eval(['print -depsc ./figs/theta_all', num2str(flnum),'.eps' ]);




%% Using particle filter to compute heat trasnfer coefficient and ambient T

% Ta0 = Ta(1);  Ta1 = Ta(10);
% varC = 1e-1; varTa = 1e-1;  % variance of Ta of the same order as Ta
% varlik = 1e-1;
% Nparticlesamp = 1e2;
% [Ccal, TaCal] = heatTranStep(dTdt, Tavg(1:njump:end-njump), ...
%     Ta(1:njump:end-njump), Ta0, Ta1, dt, varC, varTa, varlik, Nparticlesamp);



% figure()
% scatter3(Tavg(1:njump:end-njump),TaCal,Ccal);
% 
% nstepred = length(TaCal);
% 
% figure()
% plot(1:nstepred, Ta(1:njump:end-njump),'-b', ...
%         1:nstepred, TaCal,'-r')
% 
% figure()
% plot(1:nstepred,c, '-b', 1:nstepred, Ccal,'-r')
% 
% figure()
% plot(1:nstepred, c,'-r')





%% Using the log-transformed data


Ta0 = 180;  Ta1 = 175;
Taexp = Ta(1:njump:end-njump);
varC = 1e-1; varTa = 1e3;  % variance of Ta of the same order as Ta
varlik = 1e-1;
Nparticlesamp = 1e2;
[Ccal, TaCal] = heatTranStep(dTdt, Tavg(1:njump:end-njump), Taexp, ...
    Ta0, Ta1, dt, varC, varTa, varlik, Nparticlesamp);



figure()
plot(1:size(Taexp,1),Taexp,'.b', 1:size(TaCal,1),TaCal,'.r')
%eval(['print -dpng  ./figs/TaCal.png' ]);
%eval(['print -depsc ./figs/TaCal.eps' ]);



figure()
scatter3(Tavg(1:njump:end-njump),Ta(1:njump:end-njump),c,'.b'); hold on;
scatter3(Tavg(1:njump:end-njump),TaCal,Ccal,'.r');
zlabel('convective coeff','FontName',fname,'fontsize',fsize);
ylabel('ambient T (^oC)','FontName',fname,'fontsize',fsize);
xlabel('avg T (^oC)','FontName',fname,'fontsize',fsize);
% xlim([0 400]); ylim([-30 30])
ha=gca;set(ha,'linewidth',1.5,'FontName',fname,'FontSize',fsize,'Box','on');
% create a folder named figs
%eval(['print -dpng ./figs/theta_particleFilter.png' ]);
%eval(['print -depsc ./figs/theta_particleFilter.eps' ]);





















% xiTa=Tasort(1500); xiC = Csort(1500);
% e = abs(q - xiC*(T-xiTa))
% 
% [Cten,Taten] = meshgrid(C_Ta(:,1),C_Ta(:,2));
% evec = abs(q - Csort.*(Tasort-T));
% figure()
% scatter3(Csort, Tasort, evec,4,evec)
% 
% figure()
% scatter3(Csort, Tasort, likval,4,likval)
% 
% 
% figure()
% plot(C_Ta(:,1),C_Ta(:,2),'.b',Csort,Tasort,'.r')
% 
% figure()
% plot(lhsnorm(0.006,1e-5,100,'off'))


