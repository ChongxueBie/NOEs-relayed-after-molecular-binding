

% initial magnetization Iz, Sz (M10,M20)

clc;
clear;

R11=1/3; R12=0.2;  % R11: R1 for I; R12: R1 for S
R21=1;R22=1e5;   % R21: R2 for I; R22: R2 for S
wa=0; %w1: chemical shift for I; 
w2=1; %w2: chemical shift for S in ppm
% u = 10;
% n = 40;
 % k12: exchange rate from I to S; K21 vice versa
k=[30];
for wp=150 %total exchange rate
t=10; %mxing time
p = 1;
dwa=-10:0.1:10;
nn=length(dwa);
for j=1:nn
    dw1=(wa-dwa(j))*750*2*pi;
    dw2=(wa-dwa(j)+w2)*750*2*pi;
    k21=k; k12=k*p/110;

M10=110;
M20=p;
%M0=[0.5 0 0 M10 0 0 M20];
M0=[0.5 0 0 M10 0 0 M20]; 
% pulse period
%    wp=2*pi*2;
   wx=0;wy=wp*2*pi;

P= [0 0 0 0 0 0 0;...
    0 R21+k12 -dw1 wy -k21 0 0;...
    0 dw1 R21+k12 -wx 0 -k21 0;...
    -2*M10*R11 -wy wx R11+k12 0 0 -k21;...
    0 -k12 0 0 R22+k21 -dw2 wy;...
    0 0 -k12 0 dw2 R22+k21 -wx;...
    -2*M20*R12 0 0 -k12 -wy wx R12+k21];


temp1=expm(-P*t)*M0';    

M1(j)=temp1(4)/M10;
M2(j)=temp1(7)/M20;
end

plot(dwa,M1,'o', dwa,M1); hold on;

% plot(dwa,M2,'s', dwa,M2); hold on;
        set(gca,'ylim', [0 1.05],'xlim',[-6 6],'xtick',[-6:1:6]);
%          set(gca,'ylim', [0.5 1.05],'xlim',[-4000 4000]);
        set(gca,'XDir','reverse')
end 

