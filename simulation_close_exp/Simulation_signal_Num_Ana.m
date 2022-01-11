clear all
close all

path='/Volumes/CX/JHU/Beads/code/matlab/simulation_close_exp';
filename = 'e_sim_large_R22.mat';
file = fullfile(path, filename);
load(file);

filename = 'e_cal_large_R22.mat';
file = fullfile(path, filename);
load(file);

%%
figure(1);
plot([0:15],[0:15],'k--','linestyle','-','linewidth',4)
hold on

plot(e_sim,e_cal,'ko','markersize',12,'linewidth',2);
% title('Enhancement factor')
xlabel({'Numerical simulation';'rNOE (\alpha=1) (%)'});
ylabel({'Analytical calculation';'rNOE (\alpha=1) (%)'});

% set(gca,'fontsize',20,'xlim',[0,80],'ylim',[0,80]);

set(gca,'fontsize',34)
h=gcf;
set(gcf, 'Position',  [100, 100, 605, 505])
set(h,'PaperOrientation','landscape');


%%
kd = logspace(-2,2,50);

[~,ind]=max(e_sim);
kd(ind)

[~,ind]=max(e_cal);
kd(ind)

%%
figure(2)
plot(kd,e_cal,'ko-')