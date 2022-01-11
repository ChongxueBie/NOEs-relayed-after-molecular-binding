clear all
close all

path='/Volumes/CX/JHU/Beads/matlab_integrate_V2/simulation_close_exp';
filename = 'signal_Kd_Kon_log_20210721.mat';
file = fullfile(path, filename);
load(file);
signal=signal;

kbind_list = logspace(-2,8,50);
kd_list = logspace(-9,3,50);

% figure(1)
% plot(log10(kd), signal(:,20,:))
% 
% figure(2)
% plot(log10(kon), signal(30,:,:))

%%

figure(3)
imagesc(signal'*100, 'YData', log10(kbind_list), 'XData', log10(kd_list),[0, 20])

colormap(jet(256))
% altered_colormap()
colorbar()
set(gca,'fontsize',34)

set(gca,'XColor',[0 0 0]); % Set RGB value to what you want
set(gca,'YColor',[0 0 0]); % Set RGB value to what you want

ax=gca;
ax.LineWidth = 1.5;

xlabel('K_D (M)')
% ylabel(strcat('k_{on} (s^{-1})'))

ylabel('k''_{on} (s^{-1})')
% ylabel('k_{on} * [R] (M^{-1} s^{-1})')
set(gcf, 'Position',  [100, 100, 600, 500])
set(gcf,'PaperOrientation','landscape');


