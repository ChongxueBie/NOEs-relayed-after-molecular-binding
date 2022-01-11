% W---RelayOH---G
% initial magnetization: M10,M20,M30
clc;
clear;

R11=1/12;   % R11: R1 for I; R12: R1 for S
R21=1/1;R22=1/0.01; R23=1/0.01;  % R21: R2 for I; R22: R2 for S
wa=0;  w2=1.2; w3=-20;
plist = 0.1;
pp = 0; ww=0; wplist = 0:5:200;
p = plist
pp = pp+1;
k21 = [4000] 

k12 = k21*p/110;
u23 = 40; u32 = u23;

nindex = 0;
for n23 = [-200:10:-10]; 
    
    n32 = n23;
R12=-n23;R13=-n23;
nindex = nindex+1;
ww = 0;
for wp=wplist
    ww = ww+1;
    
t=100%mxing time
dwa=-25:0.1:5;
nn=length(dwa);
for j=1:nn
    dw1=(wa-dwa(j))*750*2*pi;
    dw2=(wa-dwa(j)+w2)*750*2*pi;  
    dw3=(wa-dwa(j)+w3)*750*2*pi;  
M10=110; M20 = p; M30 = p;
M0=[0.5 0 0 M10 0 0 M20 0 0 M30]; 
% pulse period
%    wp=2*pi*2;
   wx=0;wy=wp*2*pi;
P =[0                       0       0       0       0       0       0       0    0   0;...
    0                       R21+k12 dw1     -wy     -k21    0       0       0    0   0;...
    0                       -dw1    R21+k12 wx      0       -k21    0       0    0   0;...
    -2*M10*R11              wy      -wx     R11+k12 0       0       -k21    0    0   0;...
    0                       -k12    0       0       R22+k21 dw2     -wy     u32  0   0;...
    0                       0       -k12    0       -dw2    R22+k21 wx      0    u32 0;...
    -2*M20*R12-2*n32*M30    0       0       -k12    wy      -wx     R12+k21 0    0   n32;...
    0                       0       0       0       u23       0       0     R23  dw3 -wy;...
    0                       0       0       0       0       u23       0     -dw3 R23 wx;...
    -2*M30*R13-2*n23*M20    0       0       0       0       0       n23     wy   -wx R13];

temp1=expm(-P*t)*M0';    
M1(j)=temp1(4)/M10;
end
figure(1);plot(dwa,M1,'o-'); hold on;
signal(ww) =1- M1(dwa == w3);
set(gca,'XDir','reverse')
end 


figure(2); plot(wplist/42.577,signal*100,'ko'); hold on
X = wplist/42.577;
Y(:) = signal*100;
g = fittype(@(a, b, x) a*b*x.^2./(1+a*x.^2),'independent',{'x'}, 'dependent','y');

zfit = fit(X', Y',g, 'StartPoint', [0.01, 50])
ci = confint(zfit,0.68)

Xplot =0:0.05:5;zplot = zfit(Xplot);
plot(Xplot, zplot,'r-')

e_sim(nindex) = zfit.b/100/p*110;
e_cal(nindex) = -n23/((R13+k21)*(R11+k12)/k21-k12);

% e_sim2(nindex) = zfit.b;
% e_cal2(nindex) = (-n23/((R13+k21)*(R11+k12)/k21-k12)) * (M30/110)*100;
xlabel('B_1 (\muT)');
ylabel('glycoNOE intensity (%)')
set(gca, 'xlim',[0,5],'fontsize',20)

end
%%
figure(3);
plot([-200:10:-10],e_sim,'ko');
hold on;
plot([-200:10:-10],e_cal,'r-');

%%
figure(4);
plot(e_sim,e_cal,'ko');
title('Enhancement factor')
xlabel('Numerical simulation');
ylabel('Analytical calculation');
set(gca,'fontsize',20,'xlim',[0,800],'ylim',[0,800],'xtick',[100:200:1000],...
    'ytick',[100:200:1000]);
hold on; plot([0:800],[0:800],'k--')
h=gcf;
set(h,'Position',[50 50 400 400]);
set(h,'PaperOrientation','landscape');