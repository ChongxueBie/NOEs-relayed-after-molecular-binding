% W---RelayOH---LB----L
% initial magnetization: M10,M20,M30,M40
clc;
clear all;
close all;


R11=1/2.8;   % R11: R1 for I; R12: R1 for S
R21=1/1.8;R22=2.5*5e5; R23=R22;  R24=1/0.6;% R21: R2 for I; R22: R2 for S
wa=0;  w2=1.2; w3=-20; w4=-20;

plist = 0.1;
pp = 0;
kon = 1e5;
na=9;

% koff = 500;
pR = 0.001;  % receptor population:[R] M
wplist = 0:1:60; 

nindex = 0;
for p = plist
    pp = pp+1;
        
for kd = logspace(-4,2,50)
%     [1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,1,2,3,4,5,6,7,8,9,10]  
nindex = nindex+1;
%     
% LT = p;   % total ligands
% LB = (LT+pR+kd)/2-sqrt((LT+pR+kd)^2-4*LT*pR)/2; % bound ligands
% L = LT-LB; % free ligands;

L = p; % free ligands;
LB = (L/(L+kd))*pR;
LT = L+LB;   % total ligands

kbind = kon*(pR-LB); % kon*[R]
koff = kd*kon;   

M40 = na*L; M30 = na*LB; M20 = M30; M10=106; 
k21 = 1000;
k12 = k21*M20/M10;

% u23 = 4000000; u32 = u23*M20/M30;
% n23 = -4000000; n32 = n23*M20/M30;
% R12=-1*n23;R13=-1*n23;R14=1/2;

u23 = 2*5e5; u32 = u23*M20/M30;
n23 = -1*5e5; n32 = n23*M20/M30;
R12=-1*n23;R13=-1*n23;R14=1/1.8;

ww = 0;
for wp=wplist
    ww = ww+1;

for t=50; %mxing time
dwa=-25:0.1:1;
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
M1(j)=temp1(4)/M10;
end

figure(1);plot(dwa,M1,'o-'); hold on;
signal(ww) =1- M1(dwa == w3);
set(gca,'XDir','reverse')

end
end

figure(2); plot(wplist*2*pi,signal*100,'ko'); hold on
X = wplist*2*pi;
Y(:) = signal*100;
g = fittype(@(a, b, x) a*b*na*x.^2./(1+a*x.^2),'independent',{'x'}, 'dependent','y');

zfit = fit(X', Y',g, 'StartPoint', [20, 5]);
% zfit.a = 1/R13/R23;
% ci = confint(zfit,0.68)

Xplot =0:0.05:500;zplot = zfit(Xplot);
plot(Xplot, zplot,'r-')

xlabel('B_1 (\muT)');
xlabel('w_1 (rad/s)');
ylabel('glycoNOE intensity (%)')
set(gca,'ylim',[0,30],'fontsize',20)

set(gca,'fontsize',34)
set(gcf, 'Position',  [100, 100, 605, 505])
set(gcf,'PaperOrientation','landscape');

% e_sim(nindex) = (zfit.b/100*110) / LB
% e_cal(nindex) = (1/((R11/k21) + (k12/k21) - (R11/n32) + (R11/koff)))

e_sim(nindex) = na*(zfit.b)
e_cal(nindex) = (1/((R11/k21) + (k12/k21) - (R11/n32) + (R11/koff))) * (kbind/koff) * (na*L/106) *100

% e_sim(nindex) = (zfit.b/100*110) / L
% e_cal(nindex) = (1/((R11/k21) + (k12/k21) - (R11/n32) + (R11/koff))) * (kbind/koff)


end
end

%%
%%
% figure(3);
% kd = logspace(-4,-1,50);
% % [1e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,1,2,3,4,5,6,7,8,9,10];
% plot(kd,e_sim,'ko');
% hold on;
% plot(kd,e_cal,'r-');

%%
figure(4);
plot(e_sim,e_cal,'ko','markersize',12,'linewidth',2);
% title('Enhancement factor')
xlabel({'Numerical simulation';'rNOE (\alpha=1) (%)'});
ylabel({'Analytical calculation';'rNOE (\alpha=1) (%)'});
set(gca,'fontsize',20,'xlim',[0,15],'ylim',[0,15]);
hold on; 
plot([0:80],[0:80],'k--','linestyle','-','linewidth',2)
set(gca,'fontsize',34)
h=gcf;
set(gcf, 'Position',  [100, 100, 605, 505])
set(h,'PaperOrientation','landscape');

%%
path='/Volumes/CX/JHU/Beads/matlab_integrate_V2/simulation_close_exp';
filename = 'e_sim_large_R22.mat';
file = fullfile(path, filename);
save(file, 'e_sim');

filename = 'e_cal_large_R22.mat';
file = fullfile(path, filename);
save(file, 'e_cal');