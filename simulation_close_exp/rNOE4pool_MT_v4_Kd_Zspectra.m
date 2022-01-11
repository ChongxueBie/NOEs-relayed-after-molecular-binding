% W---RelayOH---LB----L
% initial magnetization: M10,M20,M30,M40
clc;
clear all;
close all;

colorList = {[0 0 0],[1 0 0],[0 0 1]};

R11=1/2.8;   % R11: R1 for I; R12: R1 for S
R21=1/1.8;R22=1/1e-6; R23=2.5*5e5; R24=R23; R25=1/0.6;% R21: R2 for I; R22: R2 for S
wa=0; w2=-2.3; w3=1.2; w4=-4; w5=-4;

plist = 1e-10;%0.4; % total ligands: M

pp = 0;
kon = 1*1e5;
na=9;

% kbind = 10^2.8;
% koff = 100;
pR = 0.001;  % receptor population:[R] M
kd_list=[0.001, 0.1, 10];

for p = plist
        
for kd = kd_list
%     [5e-2,1e-1,1.3e-1,1.5e-1,1.8e-1,2e-1,3e-1,4e-1,5e-1,6e-1,7e-1,8e-1,9e-1,1,1.2,1.4,1.6,1.8,2,3,4,5,6,7,8,9,10]  

pp = pp+1; 

% LT = p;   % total ligands
% LB = (LT+pR+kd)/2-sqrt((LT+pR+kd)^2-4*LT*pR)/2; % bound ligands
% L = LT-LB; % free ligands;

L = p; % free ligands;
LB = (L/(L+kd))*pR;
LT = L+LB;   % total ligands

% kon = koff/kd;
kbind = kon*(pR-LB); % kon*[R]
% kon = kbind/(pR-LB);
koff = kd*kon;   

M50 = na*L; M40 = na*LB; M30 = M40; M10=106; M20=0.05*M10; 
k31 = 1000;
k13 =k31*M30/M10;

k21 = 30;
k12 = k21*M20/M10;

u34 = 2*5e5; u43 = u34*M30/M40;
n34 = -1*5e5; n43 = n34*M30/M40;
R13=-1*n34;R14=-1*n34;R15=1/1.8; R12=1/1.;

for wp=600/(2*pi)

for t=50; %mxing time
dwa=-6:0.1:6;
nn=length(dwa);
for j=1:nn
    dw1=(wa-dwa(j))*750*2*pi;
    dw2=(wa-dwa(j)+w2)*750*2*pi;  
    dw3=(wa-dwa(j)+w3)*750*2*pi;  
    dw4=(wa-dwa(j)+w4)*750*2*pi;
    dw5=(wa-dwa(j)+w5)*750*2*pi;

M0=[0.5 0 0 M10 0 0 M20 0 0 M30 0 0 M30 0 0 M40];
% M0=[0.4 0 0 M40 0 0 M30 0 0 M20 0 0 M10]; 
% pulse period
%    wp=2*pi*2;
   wx=0;wy=wp*2*pi;
   
P= [0                       0                   0                   0                  0            0            0            0            0             0            0            0             0           0            0             0;...
    0                       R21+k12+k13        -dw1                 wy                -k21          0            0           -k31          0             0            0            0             0           0            0             0;...
    0                       dw1                 R21+k12+k13        -wx                 0           -k21          0            0           -k31           0            0            0             0           0            0             0;...
   -2*M10*R11               -wy                 wx                  R11+k12+k13        0            0           -k21          0            0           -k31           0            0            0             0           0             0;...
    0                      -k12                 0                   0                  R22+k21     -dw2          wy           0            0             0            0            0             0           0            0             0;...
    0                       0                  -k12                 0                  dw2          R22+k21 -wx               0            0             0            0            0             0           0            0             0;...
   -2*M20*R12               0                   0                  -k12               -wy           wx           R12+k21      0            0             0            0            0             0           0            0             0;...
    0                      -k13                 0                   0                  0            0            0            R23+k31     -dw3           wy           u43          0             0           0            0             0;...
    0                       0                  -k13                 0                  0            0            0            dw3          R23+k31      -wx           0            u43           0           0            0             0;...
   -2*M30*R13-2*n43*M40     0                   0                  -k13                0            0            0           -wy           wx            R13+k31      0            0             n43         0            0             0;...
    0                       0                   0                   0                  0            0            0            u34          0             0            R24+koff    -dw4            wy         -kbind        0             0;...
    0                       0                   0                   0                  0            0            0            0            u34           0            dw4          R24+koff     -wx          0           -kbind         0;...
   -2*M40*R14-2*n34*M30     0                   0                   0                  0            0            0            0            0             n34         -wy           wx            R14+koff    0            0            -kbind;...
    0                       0                   0                   0                  0            0            0            0            0             0           -koff         0             0           R25+kbind   -dw5           wy;...
    0                       0                   0                   0                  0            0            0            0            0             0            0           -koff          0           dw5          R25+kbind    -wx;...
    -2*M50*R15              0                   0                   0                  0            0            0            0            0             0            0            0            -koff        -wy           wx            R15+kbind];


temp1=expm(-P*t)*M0';    
M1(j)=temp1(4)/M10;
end

figure(1);
plot(dwa,M1*100,'','markersize',5,'linestyle','-','linewidth',4,'Color',colorList{pp}); hold on;
hold on;

signal(pp) =1- M1(dwa == w4);
set(gca,'XDir','reverse')

end
end
end
end

%%
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

