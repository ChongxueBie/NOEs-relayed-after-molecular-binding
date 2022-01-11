% W---RelayOH---LB----L
% initial magnetization: M10,M20,M30,M40
clc;
clear all;
close all;

% colorList = {[0 0 0],[1 0 0],[0 0 1], [0.5 0.5 0.5], [0.47,0.67,0.19]};
colorList = {[0 0 0],[1 0 0],[0 0 1], [0,1,1], [0.5 0.5 0.5], [0.47,0.67,0.19]};

R11=1/2.8;   % R11: R1 for I; R12: R1 for S
R21=1/1.8;R22=2.5*5e5; R23=R22;  R24=1/0.6;% R21: R2 for I; R22: R2 for S
wa=0;  w2=1.2; w3=-20; w4=-20;

plist = 0.1;
pp = 0;
kon = 1e5;
na=9;

% kbind_list = [10,50,100,150,200];
Rt_list = [5e-4,1e-3,5e-3,1e-2,5e-2,1e-1];

% pR = 0.001;  % free receptor population:[R] M
% free receptor population:[R] M
wplist = 0:1:100; ww=0; 

for pR = Rt_list
    p = plist;
    pp = pp+1;
        
for kd = [1e-1] 
    
L = p; % free ligands;
LB = (L/(L+kd))*pR;
LT = L+LB;   % total ligands

% kon = koff/kd;
% kon = min(kbind/(pR-LB),1e9);
kbind = kon*(pR-LB); % kon*[R]
koff = kd*kon;   

% LB/(L+LB)

M40 = na*L; M30 = na*LB; M20 = M30; M10=106; 
k21 = 1000;
k12 = k21*M20/M10;

u23 = 2*5e5; u32 = u23*M20/M30;
n23 = -1*5e5; n32 = n23*M20/M30;
R12=-1*n23;R13=-1*n23;R14=1/1.8;

ww=0;
for wp=wplist
    ww = ww+1;

for t=50;
%     50; %mxing time
dwa=-25:0.1:5;
nn=length(dwa);
for j=1:nn
    dw1=(wa-dwa(j))*750*2*pi;
    dw2=(wa-dwa(j)+w2)*750*2*pi;  
    dw3=(wa-dwa(j)+w3)*750*2*pi;  
    dw4=(wa-dwa(j)+w4)*750*2*pi;

M0=[0.5 0 0 M10 0 0 M20 0 0 M30 0 0 M40];
% M0=[0.4 0 0 M40 0 0 M30 0 0 M20 0 0 M10]; 
% pulse period
%    wp=2*pi*2;
   wx=0;wy=wp*2*pi;
P =[0                       0       0       0       0       0       0       0        0         0        0        0         0;...
    0                       R21+k12 dw1     -wy     -k21    0       0       0        0         0        0        0         0;...
    0                       -dw1    R21+k12 wx      0       -k21    0       0        0         0        0        0         0;...
    -2*M10*R11              wy      -wx     R11+k12 0       0       -k21    0        0         0        0        0         0;...
    0                       -k12    0       0       R22+k21 dw2     -wy     u32      0         0        0        0         0;...
    0                       0       -k12    0       -dw2    R22+k21 wx      0        u32       0        0        0         0;...
    -2*M20*R12-2*n32*M30    0       0       -k12    wy      -wx     R12+k21 0        0         n32      0        0         0;...
    0                       0       0       0       u23     0       0       R23+koff dw3      -wy      -kbind     0         0;...
    0                       0       0       0       0       u23     0       -dw3     R23+koff  wx       0        -kbind     0;...
    -2*M30*R13-2*n23*M20    0       0       0       0       0       n23     wy       -wx       R13+koff 0        0         -kbind;...
    0                       0       0       0       0       0       0       -koff     0         0        R24+kbind  dw4       -wy;... 
    0                       0       0       0       0       0       0       0        -koff      0        -dw4     R24+kbind   wx;...
    -2*M40*R14              0       0       0       0       0       0       0        0         -koff      wy      -wx       R14+kbind;...
    ];

temp1=expm(-P*t)*M0';    
% M1(j)=temp1(4)/M10;

M4(j)=temp1(13)/M40;
end

figure(1);plot(dwa,M4,'o-'); hold on;
signal(pp,ww) =1- M4(dwa == w3);
set(gca,'XDir','reverse')

end
end
end

figure(2);
plot(wplist,(signal(pp,:)),'-','markersize',8,'linestyle','-','linewidth',4,'Color',colorList{pp}); hold on;
hold on;

% X = wplist;
% Y(:)=(signal(pp,:)/signal(pp,end));
% 
% options = optimset('MaxFunEvals',100000*8, 'MaxIter',1000000,'Display','off');
% x = [30, 1];
% lb=[0, 0];
% ub=[inf, inf];
% [zfit_lsq,fmin,residual] = lsqnonlin(@(x)KD_fitting_w_lsq(x,X', Y'),x,lb,ub,options);
% 
% zfit_lsq

% e_sim = zfit_lsq(2)
% e_cal = (1/((R11/k21) + (k12/k21) - (R11/n32) + (R11/koff))) * (kbind/koff) * (na*L/110) *100

% Xplot =wplist;
% zplot_lsq = KD_fitting_w(zfit_lsq, Xplot);
% plot(Xplot, zplot_lsq,'-','linewidth',2,'Color',colorList{pp})


end

%%

% colorList = {[0 0 0],[1 0 0],[0 0 1], [0.5 0.5 0.5], [0.47,0.67,0.19]};
colorList = {[0 0 0],[1 0 0],[0 0 1], [0,1,1], [0.5 0.5 0.5], [0.47,0.67,0.19]};

figure(3);
pp=0;
jetcustom = copper(6);

for pR = Rt_list
    
pp=pp+1;
plot(wplist*2*pi,(signal(pp,:)),'-','markersize',8,'linestyle','-','linewidth',4, 'Color', jetcustom(pp,:)); hold on;
hold on;

end
% colormap jet

xlabel('\omega_1 (rad/s)');
ylabel('\alpha')
set(gca,'ytick',[0:0.2:200],'ylim',[0,1.05],'xtick',[0:100:1000],'xlim',[0,600],'fontsize',20)

% lg = legend('k''_{on}=10 s^{-1}', 'k''_{on}=50 s^{-1}','k''_{on}=100 s^{-1}','k''_{on}=150 s^{-1}','k''_{on}=200 s^{-1}','location','southeast','fontsize',20);
lg = legend('R_{T}=5\times10^{-4} M', 'R_{T}=1\times10^{-3} M','R_{T}=5\times10^{-3} M','R_{T}=1\times10^{-2} M','R_{T}=5\times10^{-2} M','R_{T}=1\times10^{-1} M','location','southoutside','fontsize',20,'NumColumns',2,'Box','off');


set(gca,'fontsize',30)
% set(gcf, 'Position',  [100, 100, 580, 505])
set(gcf, 'Position',  [100, 100, 530, 705])
set(gcf,'PaperOrientation','landscape');

