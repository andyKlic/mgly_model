function wrsse = COSTfunPM(paropt)

TCr = 0.03; %Total creatine pool

%% initial conditions
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

% pH data: Scopes (1974)
pHPMdata = [ 7.21 6.94 6.84 6.72 6.43 6.24 6.08 5.94 5.74 5.53 5.5 5.51 5.49 ];
tPMdata = [ 11.5 33.5 51.9 66.9 99.2 120 137.3 156.9 178.8 205.4 236.5 267.7 302.3 ];
global Par
Par = PARexpPMfun();
Par(89) = paropt(1);
Par(102) = paropt(2);
opt.RelTol=1e-4; 
opt.AbsTol=1e-7;
opt.RecomputeJACFactor=0.001;
opt.StartNewtonWithZeros = 1;
[t,CPM,stats]=radau5Mex(@ODEexpPMrad, tPMdata, initialcondition,opt);
%ode15s(@ODEexpPMrad, [0 302.3], initialcondition);%
%pH_calcPM = CPM(:,25); 
%model = CPM(:,25);data = pHPMdata';
wrsse = norm(pHPMdata'-CPM(:,25))^2;


