function wrsse = COSTfunExp45(paropt)

%%
Par = PARexp45defaultfun;
Par(16) = paropt(1);
Par(18) = paropt(2);

%%
%clear all;close all
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

t45data = [ 0 1 2 5 10 20 30 40 45 50 60 70 80 85 90 100 ...
    105 110 115 120 130 ];
options = odeset('RelTol',1e-9);
[t45,C45]=ode15s(@ODEexp29, t45data, initialcondition,options,Par);
PCr45 = C45(:,1);
Pi45 =  C45(:,8);
pH_calc45 = C45(:,24); 
protonload45 =C45(:,25);

% PCr
PCr45data = [ 0 0.00445 0.00547 0.00745 0.00833 0.00956 0.0103 0.0112...
    0.0122 0.0144 0.0169 0.019 0.0204 0.0211 0.0219 0.0237 0.0275 0.0274...
     0.016 0.00902 0.00138 ]';

% Pi
Pi45data = [ 0.035 0.0305 0.0296 0.0282 0.0268 0.0252 0.0249 0.0237 ...
    0.0226 0.0205 0.0176 0.0153 0.0137 0.0122 0.0103 0.00793 0.00322 ...
    0.00247 0.0154 0.0222 0.0304 ]';

wrsse = norm([norm((PCr45-PCr45data)), norm((Pi45-Pi45data))])^2;
