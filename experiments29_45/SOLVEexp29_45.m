%%%
clear all;close all
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
%--------------------------------------------------------------------
%%
%PARexp29frA02default;
PARexp29frA02optim;
options = odeset('RelTol',1e-9);
[t,C]=ode15s(@ODEexp29, [0 15], initialcondition,options,Par);

EVALexp29;
%%

%plot results for this experiment
PLOTexp29frA02;
%% 
% 30 min PCr and Lactate concentrations for Experiment 29, fracA = 0.002
PARexp29frA02optim;
options = odeset('RelTol',1e-9);
[t,C]=ode15s(@ODEexp29, [0 30], initialcondition,options,Par);
display('PCr at 30 min for fracA 0.002 = ')
display(C(end,1))
display('Phosphate at 30 min for fracA 0.002 = ')
display(C(end,8))
%%
initialcondition(24)=7;
PARexp29frA02optimpH7;
options = odeset('RelTol',1e-9);
[tpH7,CpH7]=ode15s(@ODEexp29, [0 15], initialcondition,options,Par);
PCrpH7 = CpH7(:,1);
LACpH7 = CpH7(:,22);

%Figure 7
figure(7);clf; set(gca,'Fontsize',14)
plot(tpH7,1000.*PCrpH7,'k--',t_data,1000.*PCr_data,...
    'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('PCr, Lactate (mM)'); box on; hold on; 
plot(tpH7,1000.*LACpH7,'k-',t_data,1000.*LAC_data,...
    'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
set(gca,'YLim',[0 16])
hold off;legend('PCr','PCr data','Lactate', 'Lactate data','Location','SouthEast');
% print -f7 -dtiff -r1200 'figure7.tiff'
fig7 = figure(7)
resp7 = fig2plotly(fig7, 'filename', 'fig7', 'strip', false)

%%
initialcondition(24)=7.8;
PARexp29frA40optim;
options = odeset('RelTol',1e-9);
[t40,C40]=ode15s(@ODEexp29, [0 10], initialcondition,options,Par);
PCr40 = C40(:,1);
Cr40 = C40(:,2);
NAD40 = C40(:,3);
NADH40 = C40(:,4);
ATP40 =  C40(:,5);
ADP40 =  C40(:,6);
AMP40 = C40(:,7);
Pi40 =  C40(:,8);
gly40 =  C40(:,9);
G1P40 = C40(:,10);
G6P40 = C40(:,11);
F6P40 =  C40(:,12);
FBP40 = C40(:,13);
DHAP40 = C40(:,14);
G3P40 = C40(:,15);
GAP40 = C40(:,16);
BPG40 =  C40(:,17);
P3G40 = C40(:,18);
P2G40 = C40(:,19);
PEP40 = C40(:,20);
PYR40 = C40(:,21);
LAC40 = C40(:,22);
mg40 = C40(:,23);
pH_calc40 = C40(:,24); 
protonload40 =C40(:,25);

pH_cy40 = ((t40 <= 1) | (t40>1 & addbuffer==0)).*pH_calc40 + ...
    (~((t40 <= 1) | (t40>1 & addbuffer==0)))*pHstat;

%Lumped variables as defined by Scopes (1973)
HMP40 = G1P40+G6P40+F6P40;
FDP40 = FBP40+(DHAP40+GAP40)/2;
%plot results for this experiment
PLOTexp29frA40;

%%
PARexp29frA40optim;
options = odeset('RelTol',1e-9);
[t,C]=ode15s(@ODEexp29, [0 30], initialcondition,options,Par);
display('PCr 30 min for fracA 0.4 = ')
display(C(end,1))
display('Phosphate 30 min for fracA 0.4 = ')
display(C(end,8))
%%
%clear all;close all
TCr = 0.03; %Total creatine pool

% initial conditions
PCri = 1e-9;
Cri = 0.03;
NADi = 0.0005;
NADHi = 1e-7;
ATPi =  1e-7;
ADPi =  0.005;
AMPi = 1e-7;
Pii =  0.03;
glyi =  0.04;
G1Pi = 1e-7;
G6Pi = 1e-7;
F6Pi =  1e-7;
FBPi = 1e-7;
DHAPi = 1e-7;
G3Pi = 1e-7;
GAPi = 1e-7;
BPGi =  1e-7;
P3Gi = 1e-7;
P2Gi = 1e-7;
PEPi = 1e-7;
PYRi = 1e-7;
LACi = 1e-7;
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

opt.RelTol=1e-4; 
opt.AbsTol=1e-7;
opt.RecomputeJACFactor=0.001;
opt.StartNewtonWithZeros = 1;
[tckk,Cckk,stats]=radauMex(@ODEexp29radCKko, [0 15], initialcondition,opt);
PCrckk = Cckk(:,1);
Crckk = Cckk(:,2);
ATPckk =  Cckk(:,5);
ADPckk =  Cckk(:,6);
AMPckk = Cckk(:,7);
Pickk =  Cckk(:,8);
LACckk = Cckk(:,22);
pH_calcckk = Cckk(:,24); 
protonloadckk =Cckk(:,25);

pH_cyckk = ((tckk <= 1) | (tckk>1 & addbuffer==0)).*pH_calcckk + ...
    (~((tckk <= 1) | (tckk>1 & addbuffer==0)))*pHstat;

%plot results for this experiment
PLOTexp29pHCKko;

%%
PARexp45optim;

%%

TCr = 0.03; %Total creatine pool

% initial conditions
PCri = 1e-9;
Cri = TCr-PCri;
NADi = 0.0005;
NADHi = 1e-9;
ATPi =  0.005;
ADPi =  1e-9;
AMPi = 1e-9;
Pii =  0.035;
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
mgi = 4.913450725E-4;
pH_calci = 7.3;
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

options = odeset('RelTol',1e-9);
[t45,C45]=ode15s(@ODEexp29, [0 130], initialcondition,options,Par);
PCr45 = C45(:,1);
Pi45 =  C45(:,8);
pH_calc45 = C45(:,24); 
protonload45 =C45(:,25);


% PCr
PCr45data = [ 0 0.00445 0.00547 0.00745 0.00833 0.00956 0.0103 0.0112...
    0.0122 0.0144 0.0169 0.019 0.0204 0.0211 0.0219 0.0237 0.0275 0.0274...
     0.016 0.00902 0.00138 ];
t45data = [ 0 1 2 5 10 20 30 40 45 50 60 70 80 85 90 100 ...
    105 110 115 120 130 ];

% Pi
Pi45data = [ 0.035 0.0305 0.0296 0.0282 0.0268 0.0252 0.0249 0.0237 ...
    0.0226 0.0205 0.0176 0.0153 0.0137 0.0122 0.0103 0.00793 0.00322 ...
    0.00247 0.0154 0.0222 0.0304 ];

% Figure 15
figure(15); clf; set(gca,'Fontsize',14)
plot(t45,1000.*PCr45,'k--',t45data,1000.*PCr45data,...
    'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('PCr, Pi (mM)'); box on; hold on; 
plot(t45,1000.*Pi45,'k-',t45data,1000.*Pi45data,...
    'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
hold off;
legend('PCr','Pi','Location','Best')
% print -f15 -dtiff -r1200 'figure15.tiff'
fig15 = figure(15)
resp15 = fig2plotly(fig15, 'filename', 'fig15', 'strip', false)

