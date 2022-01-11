% W---RelayOH---LB----L
% initial magnetization: M10,M20,M30,M40
clc;
clear all;
close all;

colorList = {[0 0 0],[1 0 0],[0 0 1]};

R11=1/2.8;   % R11: R1 for I; R12: R1 for S
R21=1/1.8;R22=2.5*5e5; R23=R22;  R24=1/0.6;% R21: R2 for I; R22: R2 for S
wa=0;  w2=1.2; w3=-4; w4=-4;

plist = 0.1;%0.4; % total ligands: M

pp = 0;
kon = 1*1e5;
% koff = 1e3;
na=9;

pR = 0.001;  % receptor population:[R] M

% pR = 0.21;
kd_list=[0.001, 0.1, 10];

for p = plist
        
for kd = kd_list
pp = pp+1; 

% LT = p;   % total ligands
% LB = (LT+pR+kd)/2-sqrt((LT+pR+kd)^2-4*LT*pR)/2; % bound ligands
% L = LT-LB; % free ligands;

L = p; % free ligands;
LB = (L/(L+kd))*pR;
LT = L+LB;   % total ligands

% kon = k off/kd;
kbind = kon*(pR-LB); % kon*[R]
% kon = kbind/(pR-LB);
koff = kd*kon;  

M40 = na*L; M30 = na*LB; M20 = M30; M10=106; 
k21 = 1000;
% k12 = 0.05/110
k12 =k21*M20/M10;
% 
% u23 = 4000000; u32 = u23*M20/M30;
% n23 = -4000000; n32 = n23*M20/M30;
% R12=-1*n23;R13=-1*n23;R14=1/2;

u23 = 2.5e5; u32 = u23*M20/M30;
n23 = -1*5e5; n32 = n23*M20/M30;
R12=-1*n23;R13=-1*n23;R14=1/1.8;

% for wp= 187/(2*pi)
for wp=600/(2*pi)

for t=50; %mxing time
dwa=-6:0.1:6;
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

figure(1);
plot(dwa,M1*100,'','markersize',5,'linestyle','-','linewidth',4,'Color',colorList{pp}); hold on;
hold on;

signal(pp) =1- M1(dwa == w3);
set(gca,'XDir','reverse')

end
end
end
end

%%
% kon = koff./kd
L = p; % free ligands;
% LB = (LT+pR+kd)./2-sqrt((LT+pR+kd).^2-4*LT*pR)./2 % bound ligands
LB = (L./(L+kd_list)).*pR
kbind = kon.*(pR-LB)
% kon = koff./kd_list
% kbind = kon.*(pR-LB) % kon*[R]


set(gca,'ylim', [70 100],'xlim',[-6 6],'xtick',[-10:2:10]);
% set(gca,'ylim', [0.6 1.01],'xlim',[-3 3]);
set(gca,'XDir','reverse')
% 
% lg = legend(['K_D=0.01 M' char(10) '(k_{on}*[R]=36 s^{-1})'], ['K_D=1 M' char(10) '(k_{on}*[R]=364 s^{-1})'],['K_D=100 M' char(10) '(k_{on}*[R]=400 s^{-1})']);

lg = legend(['K_D=10^{-3} M'], ['K_D=10^{-1} M'],['K_D=10^{1} M'],'location','southwest');
% 
% ax = gca;
% neworder = [3,1,2];
% ax.Children = ax.Children(neworder);
% 
% neworder = [3,1,2];
% ax.Children = ax.Children(neworder);

set(gca,'fontsize',34)
set(lg, 'fontsize',30);
xlabel('Offset from water (ppm)')
ylabel('S/S_0 (%)')
h=gcf;
set(gcf, 'Position',  [100, 100, 655, 505])
set(h,'PaperOrientation','landscape');
