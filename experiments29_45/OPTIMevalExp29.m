%%%
clear all; close all
parinit=[0.00304003,1.75,0.01849588/2];

[paropt,fval]=fmincon(@COSTfunExp29,parinit,[],[],[],[],...
    [9.7e-5, 0.175, 1.04/100],[1e-2, 2, 1.04]);

%%
Par = PARexp29frA02defaultfun;
Par(16) = paropt(1);
Par(18) = paropt(2);
Par(39) = paropt(3);

TCr = 0.03; %Total creatine pool

% initial conditions
PCri = 1e-9;
Cri = TCr - PCri;
NADi = 0.0005;
NADHi = 1e-9;
ATPi =  0.005;
ADPi =  1e-9;
AMPi = 1e-9;
Pii =  0.03;
glyi =  0.04;
G1Pi = 1e-9;
G6Pi = 1e-9;
F6Pi =  1e-9;
FBPi = 1e-9;
DHAPi = 1e-9;
G3Pi = 1e-9;
GAPi = 1e-9;
BPGi =  1e-9;
P3Gi = 1e-9;
P2Gi = 1e-9;
PEPi = 1e-9;
PYRi = 1e-9;
LACi = 1e-9;
mgi = 5.132658807e-4;
pH_calci = 7.8;
protonloadi =0;

%% 
initialcondition = [PCri;
Cri;
NADi;
NADHi;
ATPi;
ADPi;
AMPi;
Pii;
glyi;
G1Pi;
G6Pi;
F6Pi;
FBPi;
DHAPi;
G3Pi;
GAPi;
BPGi;
P3Gi;
P2Gi;
PEPi;
PYRi;
LACi;
mgi;
pH_calci;
protonloadi];

t_data = [ 0 0.0833 0.167 0.333 0.5 0.667 0.833 1 2 3 5 10 15 ];
[t, C]= ode15s(@ODEexp29,[0 15],initialcondition,[],Par);

%curve1
    PCr = C(:,1);
    %curve2
    LAC = C(:,22);
    %curve3
    ATP =  C(:,5);
    %curve4
    ADP =  C(:,6);
    %curve5
    AMP = C(:,7);
    %curve6
    G1P = C(:,10);
    G6P = C(:,11);
    F6P =  C(:,12);
    HMP = G1P+G6P+F6P;
    %curve7
    FBP = C(:,13);
    DHAP = C(:,14);
    GAP = C(:,16);
    FDP = FBP+(DHAP+GAP)./2;
    %curve8
    G3P = C(:,15);
    %curve9
    P3G = C(:,18);
    P2G = C(:,19);
    PYR = C(:,21);
    PG =  (P3G + P2G)+PYR;


PCr_data = [ 0 0.00278 0.00534 0.00845 0.00976 0.00957 0.01 ...
    0.0102 0.0116 0.0125 0.0134 0.0144 0.015 ];

LAC_data = [ 0 8.61E-6 0.00219 0.00443 0.00434 0.00509 0.00515 ...
    0.00518 0.00649 0.00671 0.00767 0.00895 0.00932 ];

ATP_data = [ 0.00497 0.0032 0.00443 0.0047 0.00481 0.00472 0.00483 ...
    0.00491 0.005 0.00481 0.00499 0.00487 0.00499 ];

ADP_data = [ 0 0.00184 8.59E-4 3.91E-4 3.91E-4 3.23E-4 2.94E-4 ...
    3.03E-4 2.55E-4 2.45E-4 2.64E-4 2.35E-4 1.96E-4 ];

AMP_data = [ 0 8.68E-4 1.87E-4 3.07E-5 2.1E-5 2.1E-5 3.07E-5 ... 
    2.1E-5 0 0 0 0 0 ];

HMP_data = [ 0 5.34E-4 2.42E-4 2.04E-6 3.98E-6 5.93E-6 2.93E-5 ...
    5.08E-5 3.91E-5 3.13E-5 2.35E-5 2.74E-5 3.13E-5 ];

FDP_data = [ 0 5.77E-4 0.00101 5.66E-5 0 0 0 0 0 0 0 0 0 ];

G3P_data = [ 0 2.72E-4 6.93E-4 6.18E-4 4.82E-4 3.21E-4 2.44E-4 ...
    2.13E-4 1.59E-4 1.49E-4 1.08E-4 1.31E-4 1.43E-4 ];

PG_data = [ 0 3.11E-4 7.31E-4 7.29E-4 5.4E-4 4.33E-4 3.3E-4 2.81E-4 ...
    1.99E-4 1.76E-4 1.44E-4 1.46E-4 1.58E-4 ];

%time in minutes
t_data = [ 0 0.0833 0.167 0.333 0.5 0.667 0.833 1 2 3 5 10 15 ];

% Figure 3
figure(3); clf; set(gca,'Fontsize',14)
semilogy(t,1000*ATP,'k--',t_data,1000*ATP_data,'ko','linewidth',1.5, ...
    'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('ATP, ADP, AMP (mM)'); box on; hold on; 
semilogy(t,1000*ADP,'k-',t_data,1000*ADP_data,'k^','linewidth',1.5, ...
    'Markerfacecolor',0.75*[1 1 1],'Markersize',8);
semilogy(t,1000*AMP,'k-.',t_data(1:4),1000*AMP_data(1:4),'ks', ...
    'linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
hold off;legend('ATP','ATP data','ADP','ADP data','AMP','AMP data','Location','SouthEast')
%print -f3 -deps -tiff 'figure3.eps'
fig3 = figure(3)
resp3 = fig2plotly(fig3, 'filename', 'fig3', 'strip', false)


%Figure 4
figure(4); clf; set(gca,'Fontsize',14)
plot(t,1000*PCr,'k--',t_data,1000*PCr_data,'ko','linewidth',1.5, ...
    'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('PCr, Lactate (mM)'); box on; hold on; 
plot(t,1000*LAC,'k-',t_data,1000.*LAC_data,'k^','linewidth',1.5, ...
    'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
hold off;legend('PCr','PCr data','Lactate', 'Lactate data','Location','Best')
%print -f4 -deps -tiff 'figure4.eps'
fig4 = figure(4)
resp4 = fig2plotly(fig4, 'filename', 'fig4', 'strip', false)


%Figure 5
figure(5); clf; set(gca,'Fontsize',14)
plot(t,1000*HMP,'k--',t_data,1000*HMP_data,'ko','linewidth',1.5, ...
    'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); 
ylabel('Hexose monophosphates, Fructose diphosphates (mM)'); 
box on; hold on; 
plot(t,1000*FDP,'k-',t_data(1:4),1000*FDP_data(1:4),'k^','linewidth', ...
    1.5,'Markerfacecolor',0*[1 1 1],'Markersize',8); 
hold off;legend('HMP', 'HMP data', 'FDP', 'FDP data','Location','Best')
%print -f5 -deps -tiff 'figure5.eps'

%Figure 6
figure(6); clf; set(gca,'Fontsize',14)
plot(t,1000*G3P,'k--',t_data,1000*G3P_data,'ko','linewidth',1.5, ...
    'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); 
ylabel('Glycerol-3-Phosphate, Phosphoglycerates (mM)'); box on; hold on; 
plot(t,1000*PG,'k-',t_data,1000*PG_data,'k^','linewidth',1.5, ...
    'Markerfacecolor',0*[1 1 1],'Markersize',8); 
hold off;legend('G3P','G3P data','PG','PG data','Location','Best')
%print -f6 -deps -tiff 'figure6.eps'
