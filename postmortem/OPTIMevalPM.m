%%%
clear all;close all
TCr = 0.03; %Total creatine pool

% initial conditions
PCri = 0.024;
Cri = TCr - PCri;
NADi = 0.0005;
NADHi = 1e-9;
ATPi =  0.005;
ADPi =  1e-9;
AMPi = 1e-9;
IMPi = 1e-9;
Pii =  0.006;
glyi =  0.065;
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
mgi = 7.282174603E-4;
pH_calci = 7.25;
protonloadi =0;

%% 

initialcondition = [PCri;
Cri;
NADi;
NADHi;
ATPi;
ADPi;
AMPi;
IMPi;
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
options = optimset('TolX',1e-9,'TolFun',1e-9,'MaxFunEvals',10000);
paropt = fmincon(@COSTfunPM2,[5e-4,6e-3],[],[],[],[],...
     [1e-4,5e-3]./10,[1e-4,2e-3].*10,[],options);

%%

global Par
Par = PARexpPMfun();
Par(89) = paropt(1);
Par(102) = paropt(2);
opt.RelTol=1e-14; 
opt.AbsTol=1e-14;
opt.RecomputeJACFactor=0.001;
opt.StartNewtonWithZeros = 1;
[t,CPM,stats]=radauMex(@ODEexpPMrad, [0 302.3], initialcondition,opt);

pH_calcPM = CPM(:,25); 

pH_cyPM = pH_calcPM; 

% pH data: Scopes (1974)
pHPMdata = [ 7.21 6.94 6.84 6.72 6.43 6.24 6.08 5.94 5.74 5.53 5.5 5.51 5.49 ];
tPMdata = [ 11.5 33.5 51.9 66.9 99.2 120 137.3 156.9 178.8 205.4 236.5 267.7 302.3 ];

%Figure 12
figure(12);clf; set(gca,'Fontsize',12)
plot(t,pH_cyPM,'k-','linewidth',1.5);
xlabel('Time (min)'); ylabel('pH'); box on; hold on;
plot(tPMdata,pHPMdata,'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
hold off;
legend('Model', 'Data','Location','Best')
% print -f13 -dtiff -r1200 'figure12.tiff'
% print -f13 -deps -tiff 'figure12.eps'
%%

EVALexpPM;
figure(13);clf; set(gca,'Fontsize',12)
xlabel('Time (min)'); ylabel('Proton consumption fluxes (M/min)'); box on;hold on;
plot(t,CKprtflux,'k-','linewidth',1.5);
plot(t,ATPase.*deltaH_ATPase,'k--','linewidth',1.5);
plot(t,glycprtflux,'k:','linewidth',1.5);
plot(t, v_AMPDA.*deltaH_AMPDA,'k-.','linewidth',1.5);
plot(t,protons_consumed,'k','linewidth',3);
legend('CK','ATPase','Glycogenolysis','AMPDA','Total','Location','NorthEast')
set(gca,'YLim',[-0.001 0.001]);set(gca,'XLim',[0 302.3]);hold off;
% print -f14 -dtiff -r1200 'figure13.tiff'
% print -f14 -deps -tiff 'figure13.eps'
