
% Abbreviations used:
% Enzymes
%   GP = glycogen phosphorylase
%   PGLM = phosphoglucomutase
%   PGI = phosphoglucoisomerase
%   PFK = phosphofructokinase
%   ALD = aldolase
%   TPI = triose phosphate isomerase
%   GAPDH = glyceraldehyde 3-phosphate dehydrogenase
%   G3PDH = glycerol 3-phosphate dehydrogenase
%   PGK = phosphoglycerate kinase
%   PGM = phosphoglycerate mutase
%   EN = enolase
%   PK = pyruvate kinase
%   LDH = lactate dehydrogenase
%   CK = creatine kinase
%   ADK = adenylate kinase
%   Metabolites (variables)
%   gly = glycogen
%   G1P = glucose 1-phosphate
%   G6P = glucose 6-phosphate
%   F6P = fructose 6-phosphate
%   FBP = fructose bisphosphate
%   DHAP = dihydroxyacetone phosphate
%   G3P = glycerol 3-phosphate
%   GAP = glyceraldehyde 3-phosphate
%   BPG = 1,3 bisphosphoglycerate
%   P3G = 3-phosphoglycerate
%   P2G = 2-phosphoglycerate
%   PEP = phosphoenolpyruvate
%   PYR = pyruvate
%   LAC = lactate
%   PCr = phosphocreatine
%   Cr = creatine
%   Pi = inorganic phosphate
%   AMP, ADP, ATP, NAD, NADH = ditto

global Par
%Par = PARexpPMfun;

%GLOBAL PARAMETERS
 R = 8.314e-3;% KJ./(K.*mol); %Ideal gas constant
 T1 = 298.15;% K; %Temp at which physical constants are reported
%   T = 303.15;% K; %Experimental temperature   
T = Par(94);
%   I = 0.1;%0.17;% M; %Ionic strength
I = Par(90);
    
%Metabolite concentrations and pH
PCr = CPM(:,1);
Cr = CPM(:,2);
NAD = CPM(:,3);
NADH = CPM(:,4);
ATP =  CPM(:,5);
ADP =  CPM(:,6);
AMP = CPM(:,7);
IMP = CPM(:,8);
Pi =  CPM(:,9);
gly =  CPM(:,10);
G1P = CPM(:,11);
G6P = CPM(:,12);
F6P =  CPM(:,13);
FBP = CPM(:,14);
DHAP = CPM(:,15);
G3P = CPM(:,16);
GAP = CPM(:,17);
BPG =  CPM(:,18);
P3G = CPM(:,19);
P2G = CPM(:,20);
PEP = CPM(:,21);
PYR = CPM(:,22);
LAC = CPM(:,23);
mg = CPM(:,24);
pH_calc = CPM(:,25); 
protonload =CPM(:,26);
    
TCr = 0.03;

%  addbuffer = 1;% dimensionless;
addbuffer = Par(97);
%  pHstat = 7.4 ;%dimensionless; 
pHstat = Par(98);

%Experiment
expno = Par(100);
%Cytosolic pH
pH_cy = (expno==6)*pH_calc + ...
    (1-(expno==6))*(((t <= 1) | (t>1 & addbuffer==0)).*pH_calc + ...
    (~((t <= 1) | (t>1 & addbuffer==0))).*pHstat);

%--------------------------------------------------------------------------
%IONIC EQUILIBRIA, BUFFERING AND APPARENT EQUILIBRIUM CONSTANTS
%Ionic equilibria and Gibbs energy of formation of reference species
%Variable symbol definitions and units
%
%Symbol     Definition          Units
%
%P_metabolite       binding polynomial      dimensionless
%NH_species         No. of hydrogen atoms       dimensionless 
%           on a species
%Navg_metabolite    Average proton binding      dimensionless
%deltaGof_species   Gibbs energy of formation at    KJ./mol 
%           zero ionic strength
%deltaGpof_species  Gibbs energy of formation at    KJ./mol 
%           a given ionic strength and pH
%deltaGpo_reaction  Reference Gibbs energy of a KJ./mol 
%           reaction

%  mgT = 0.005;% total Mg++
mgT = Par(95);
%  k = 0.08;% total K+
k = Par(96);
c0 = 1; % reference concentration

%Correction factors

RT2dadT = 1.4775;% KJ.*L.^0.5./mol.^1.5;  
B = 1.6;% M.^-0.5;
Icorr=(RT2dadT.*I.^0.5./(1+B.*I.^0.5));
I1 = 0.1;% M, Ionic strength at which pKa's were reported 
alphadebye = 1.17582;% M.^-0.5; 
IcorrpKa = alphadebye.*(I1.^0.5./(1+B.*I1.^0.5)-I.^0.5./(1+B.*I.^0.5))./log(10);
TcorrpKa = ((1./T-1./T1)./(log(10).*R));
RTalpha = 2.91482; % KJ.*L.^0.5./mol.^1.5; 
IcorrdeltaGpof = RTalpha.*I.^0.5./(1+B.*I.^0.5);

pKak_Pi = 0.5; % dimensionless;
deltaH1o_Pi = 3;% KJ./mol,
deltaHmgo_Pi = -2.9;% KJ./mol;
deltaH1_Pi = deltaH1o_Pi + Icorr.*(2.^2+1.^2-1.^2);
deltaHmg_Pi = deltaHmgo_Pi + Icorr.*(2.^2+2.^2-0.^2);
pKa1_Pi = 6.75+IcorrpKa.*(2.^2+1.^2-1.^2) + TcorrpKa.*deltaH1_Pi;
pKamg_Pi = 1.65+IcorrpKa.*(2.^2+2.^2-0.^2) + TcorrpKa.*deltaHmg_Pi;
P_Pi = 1+10.^(-pH_cy+pKa1_Pi) + (mg./c0).*(10.^(pKamg_Pi))+ ...
    (k./c0).*(10.^(pKak_Pi));
HPi2 = 1./P_Pi;
H2Pi1 = 10.^(-pH_cy+pKa1_Pi).*HPi2;
kPi = (k./c0).*(10.^(pKak_Pi)).*HPi2;
mgPi =(mg./c0).*(10.^(pKamg_Pi)).*HPi2;
Navg_Pi = 1.*H2Pi1;

dNavgPidH = 10.^(pKa1_Pi).*(P_Pi-10.^(-pH_cy+pKa1_Pi))./(c0.*P_Pi.^2);
dNavgPidmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_Pi+pKamg_Pi)./P_Pi.^2;
 
dmgPidmg = (P_Pi.*(10.^pKamg_Pi)./c0-...
(mg./c0).*(10.^(2.*pKamg_Pi))./c0)./P_Pi.^2;
dmgPidpH = (mg./c0).*10.^(-pH_cy+pKa1_Pi+pKamg_Pi).*log(10)./P_Pi.^2;

NH_HPi2 =  1;% dimensionless 
deltaGof_HPi2 = -1096.10;% KJ./mol;
deltaGpof_HPi2 = deltaGof_HPi2 + NH_HPi2.*log(10).*R.*T.*pH_cy - ...
IcorrdeltaGpof.*(4-NH_HPi2);
    
deltaH1o_ATP = -5;% KJ./mol,
deltaHmgo_ATP = -18;% KJ./mol,
deltaHko_ATP = -1;% KJ./mol;
deltaH1_ATP = deltaH1o_ATP + Icorr.*(4.^2+1.^2-3.^2);
deltaHmg_ATP = deltaHmgo_ATP+Icorr.*(4.^2+2.^2-2.^2);
deltaHk_ATP = deltaHko_ATP+Icorr.*(4.^2+1.^2-3.^2);
pKa1_ATP = 6.48+IcorrpKa.*(4.^2+1.^2-3.^2) + TcorrpKa.*deltaH1_ATP;
pKamg_ATP = 4.19+IcorrpKa.*(4.^2+2.^2-2.^2)+TcorrpKa.*deltaHmg_ATP;
pKak_ATP = 1.17+IcorrpKa.*(4.^2+1.^2-3.^2)+TcorrpKa.*deltaHk_ATP;   
P_ATP = 1 + 10.^(-pH_cy+pKa1_ATP) + (mg./c0).*10.^pKamg_ATP + ...
    (k./c0).*10.^pKak_ATP;
ATP4 = 1./P_ATP;
HATP3 = 10.^(-pH_cy+pKa1_ATP).*ATP4;
mgATP2 = (mg./c0).*10.^pKamg_ATP.*ATP4;
kATP = (k./c0).*10.^pKak_ATP.*ATP4;
Navg_ATP = (0.*ATP4 + 1.*HATP3 + 0.*mgATP2 + 0.*kATP);
dNavgATPdH = 10.^(pKa1_ATP).*(P_ATP-10.^(-pH_cy+pKa1_ATP))./(c0.*P_ATP.^2);
dNavgATPdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_ATP+pKamg_ATP)./P_ATP.^2;

dmgATP2dmg = (P_ATP.*(10.^pKamg_ATP)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_ATP))./c0)./P_ATP.^2;
dmgATP2dpH = (mg./c0).*10.^(-pH_cy+pKa1_ATP+pKamg_ATP).*log(10)./P_ATP.^2;

NH_ATP4 = 12; deltaGof_ATP4 = -2768.10;% KJ./mol;
deltaGpof_ATP4 = deltaGof_ATP4 + NH_ATP4.*R.*T.*log(10).*pH_cy - ...
    IcorrdeltaGpof.*(16-NH_ATP4);       

pKak_ADP=1;% dimensionless;
deltaH1o_ADP = -3;% KJ./mol,
deltaHmgo_ADP = -15;% KJ./mol;
deltaH1_ADP = deltaH1o_ADP+Icorr.*(3.^2+1.^2-2.^2);
deltaHmg_ADP = deltaHmgo_ADP+Icorr.*(3.^2+2.^2-1.^2);
pKa1_ADP = 6.38+IcorrpKa.*(3.^2+1.^2-2.^2)+TcorrpKa.*deltaH1_ADP;
pKamg_ADP = 3.25+IcorrpKa.*(3.^2+2.^2-1.^2)+TcorrpKa.*deltaHmg_ADP;
P_ADP = (1 + 10.^(-pH_cy+pKa1_ADP) + (mg./c0).*10.^pKamg_ADP + ...
    (k./c0).*10.^pKak_ADP);
ADP3 = 1./P_ADP;
HADP2 = 10.^(-pH_cy+pKa1_ADP).*ADP3;
mgADP = ADP3.*(mg./c0).*10.^pKamg_ADP;
kADP = ADP3.*(k./c0).*10.^pKak_ADP;
Navg_ADP = (0.*ADP3 + 1.*HADP2 + 0.*mgADP + 0.*kADP);
dNavgADPdH = 10.^(pKa1_ADP).*(P_ADP-10.^(-pH_cy+pKa1_ADP))./(c0.*P_ADP.^2);
dNavgADPdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_ADP+pKamg_ADP)./P_ADP.^2;
dmgADPdmg = (P_ADP.*(10.^pKamg_ADP)./c0-...
        (mg./c0).*(10.^(2.*pKamg_ADP))./c0)./P_ADP.^2;
dmgADPdpH = (mg./c0).*10.^(-pH_cy+pKa1_ADP+pKamg_ADP).*log(10)./P_ADP.^2;
NH_ADP3 = 12; deltaGof_ADP3 = -1906.13;
deltaGpof_ADP3 = deltaGof_ADP3 + NH_ADP3.*R.*T.*log(10).*pH_cy -... 
    IcorrdeltaGpof.*(9-NH_ADP3);        

    
deltaH1o_AMP = -3;% KJ./mol;
deltaHmgo_AMP = -7.5;% KJ./mol;
deltaH1_AMP = deltaH1o_AMP + Icorr.*(2.^2+1.^2-1.^2);
deltaHmg_AMP = deltaHmgo_AMP+Icorr.*(2.^2+2.^2-0.^2);
pKa1_AMP = 6.29+IcorrpKa.*(2.^2+1.^2-1.^2)+TcorrpKa.*deltaH1_AMP;
pKamg_AMP = 1.92+IcorrpKa.*(2.^2+2.^2-0.^2)+TcorrpKa.*deltaHmg_AMP;
P_AMP = (1 + 10.^(pKa1_AMP-pH_cy)+ (mg./c0).*10.^pKamg_AMP);
AMP2 = 1./P_AMP;
HAMP1 = AMP2.*10.^(pKa1_AMP-pH_cy);
mgAMP = (mg./c0).*10.^pKamg_AMP.*AMP2;
Navg_AMP = 0.*AMP2 + HAMP1 + 0.*mgAMP;
dNavgAMPdH = 10.^(pKa1_AMP).*(P_AMP-10.^(-pH_cy+pKa1_AMP))./(c0.*P_AMP.^2);
dNavgAMPdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_AMP+pKamg_AMP)./P_AMP.^2;
dmgAMPdmg = (P_AMP.*(10.^pKamg_AMP)./c0-...
        (mg./c0).*(10.^(2.*pKamg_AMP))./c0)./P_AMP.^2;
dmgAMPdpH = (mg./c0).*10.^(-pH_cy+pKa1_AMP+pKamg_AMP).*log(10)./P_AMP.^2;

 NH_AMP2 = 12; deltaGof_AMP2 = -1040.45;% KJ./mol;
deltaGpof_AMP2 = deltaGof_AMP2 + NH_AMP2.*log(10).*R.*T.*pH_cy - ...
        IcorrdeltaGpof.*(4-NH_AMP2);
    
pKamg_IMP = 1.67;% dimensionless;
deltaH1o_IMP = -2;% KJ./mol;
deltaH1_IMP = deltaH1o_IMP + Icorr.*(2.^2+1.^2-1.^2);
pKa1_IMP = 6.34+IcorrpKa.*(2.^2+1.^2-1.^2)+TcorrpKa.*deltaH1_IMP;
P_IMP = (1 + 10.^(pKa1_IMP-pH_cy)+ (mg./c0).*10.^pKamg_IMP);
IMP2 = 1./P_IMP;
HIMP1 = IMP2.*10.^(pKa1_IMP-pH_cy);
mgIMP = (mg./c0).*10.^pKamg_IMP.*IMP2;
Navg_IMP = HIMP1;
dNavgIMPdH = 10.^(pKa1_IMP).*(P_IMP-10.^(-pH_cy+pKa1_IMP))./(c0.*P_IMP.^2);
dNavgIMPdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_IMP+pKamg_IMP)./P_IMP.^2;
dmgIMPdmg = (P_IMP.*(10.^pKamg_IMP)./c0-....
                (mg./c0).*(10.^(2.*pKamg_IMP))./c0)./P_IMP.^2;
dmgIMPdpH = (mg./c0).*10.^(-pH_cy+pKa1_IMP+pKamg_IMP).*log(10)./P_IMP.^2;
NH_IMP2 = 11;

pKak_PCR = 0.31;% dimensionless;
deltaH1o_PCR = 2.66;% KJ./mol,
deltaHmgo_PCR = 8.19;% KJ./mol;
deltaH1_PCR = deltaH1o_PCR+Icorr.*(2.^2+1.^2-1.^2);
deltaHmg_PCR = deltaHmgo_PCR+Icorr.*(2.^2+2.^2-0.^2);
pKa1_PCR = 4.5+IcorrpKa.*(2.^2+1.^2-1.^2)+TcorrpKa.*deltaH1_PCR;
pKamg_PCR = 1.6+IcorrpKa.*(2.^2+2.^2-0.^2)+TcorrpKa.*deltaHmg_PCR;
P_PCR = (1  + 10.^(-pH_cy+pKa1_PCR) + (mg./c0).*10.^pKamg_PCR+ ...
    (k./c0).*10.^pKak_PCR);
HPCR = 1./P_PCR; 
H2PCR = 10.^(-pH_cy+pKa1_PCR).*HPCR;
kPCR = (k./c0).*10.^pKak_PCR.*HPCR;
mgPCR = (mg./c0).*10.^pKamg_PCR.*HPCR;
Navg_PCR = H2PCR;
dNavgPCRdH = 10.^(pKa1_PCR).*(P_PCR-10.^(-pH_cy+pKa1_PCR))./(c0.*P_PCR.^2);
dNavgPCRdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_PCR+pKamg_PCR)./P_PCR.^2;
dmgPCRdmg = (P_PCR.*(10.^pKamg_PCR)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_PCR))./c0)./P_PCR.^2;
dmgPCRdpH = (mg./c0).*10.^(-pH_cy+pKa1_PCR+pKamg_PCR).*log(10)./P_PCR.^2;

 NH_HPCR = 8;% dimensionless;       

pKa1_CR = 2.3;% dimensionless;
P_CR = 1 + 10.^(-pH_cy+pKa1_CR);    
HCR = 1./P_CR;
H2CR = HCR.*10.^(-pH_cy+pKa1_CR);
Navg_CR  = (0.*HCR + 1.*H2CR);
dNavgCRdH = 10.^(pKa1_CR).*(P_CR-10.^(-pH_cy+pKa1_CR))./(c0.*P_CR.^2);
dNavgCRdmg = 0;     
NH_HCR = 9;% dimensionless;

deltaH1o_G1P = -1.7; %KJ./mol,
deltaHmgo_G1P = -12; %KJ./mol;
deltaH1_G1P = deltaH1o_G1P+Icorr.*(2.^2+1.^2-1.^2);
deltaHmg_G1P = deltaHmgo_G1P+Icorr.*(2.^2+2.^2-0.^2);
pKa1_G1P = 6.09+IcorrpKa.*(2.^2+1.^2-1.^2)+TcorrpKa.*deltaH1_G1P;
pKamg_G1P = 2.48+IcorrpKa.*(2.^2+2.^2-0.^2)+TcorrpKa.*deltaHmg_G1P;
P_G1P=(1 + 10.^(-pH_cy).*10.^pKa1_G1P + (mg./c0).*10.^pKamg_G1P);
UG1P = 1./P_G1P;
HG1P = UG1P.*10.^(-pH_cy+pKa1_G1P);
mgG1P = UG1P.*(mg./c0).*10.^pKamg_G1P;
Navg_G1P = HG1P;
dNavgG1PdH = 10.^(pKa1_G1P).*(P_G1P-10.^(-pH_cy+pKa1_G1P))./(c0.*P_G1P.^2);
dNavgG1Pdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_G1P+pKamg_G1P)./P_G1P.^2;

dmgG1Pdmg = (P_G1P.*(10.^pKamg_G1P)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_G1P))./c0)./P_G1P.^2;
dmgG1PdpH = (mg./c0).*10.^(-pH_cy+pKa1_G1P+pKamg_G1P).*log(10)./P_G1P.^2;

NH_UG1P = 11;
deltaGof_UG1P = -1756.87;% KJ./mol;
deltaGpof_UG1P = deltaGof_UG1P + NH_UG1P.*log(10).*R.*T.*pH_cy - ...
    IcorrdeltaGpof.*(4-NH_UG1P);

pKa1_G6P = 6.11;% dimensionless;        
P_G6P = (1 + 10.^(-pH_cy+pKa1_G6P)); 
UG6P = 1./P_G6P;
HG6P = UG6P.*10.^(-pH_cy+pKa1_G6P);
Navg_G6P = HG6P;
dNavgG6PdH = 10.^(pKa1_G6P).*(P_G6P-10.^(-pH_cy+pKa1_G6P))./(c0.*P_G6P.^2);
dNavgG6Pdmg = 0;
NH_UG6P = 11; 
deltaGof_UG6P = -1763.94;% KJ./mol;
deltaGpof_UG6P = deltaGof_UG6P + NH_UG6P.*log(10).*R.*T.*pH_cy - ...
    IcorrdeltaGpof.*(4-NH_UG6P);

pKa1_F6P = 5.89;% dimensionless;
P_F6P = (1 + 10.^(-pH_cy+pKa1_F6P));
UF6P = 1./P_F6P;
HF6P = UF6P.*10.^(-pH_cy+pKa1_F6P);
Navg_F6P = HF6P;
dNavgF6PdH = 10.^(pKa1_F6P).*(P_F6P-10.^(-pH_cy+pKa1_F6P))./(c0.*P_F6P.^2);
dNavgF6Pdmg = 0;
NH_UF6P = 11;
deltaGof_UF6P = -1760.80;% KJ./mol;
deltaGpof_UF6P = deltaGof_UF6P + NH_UF6P.*log(10).*R.*T.*pH_cy -... 
        IcorrdeltaGpof.*(4-NH_UF6P);

pKa1_FDP = 6.4; %dimensionless, 
pKa2_FDP = 5.92; %dimensionless,
pKamg_FDP = 2.7; %dimensionless;    
P_FDP = (1 + 10.^(-pH_cy+pKa1_FDP) + 10.^(-2.*pH_cy+pKa1_FDP+pKa2_FDP)+...
    (mg./c0).*10.^(pKamg_FDP));
UFDP = 1./P_FDP;
HFDP = UFDP.*10.^(-pH_cy+pKa1_FDP);
H2FDP = UFDP.*10.^(-2.*pH_cy+pKa1_FDP+pKa2_FDP);
%FruBP4- + Mg2+  --> FruBP2-
mgFDP = UFDP.*(mg./c0).*10.^(pKamg_FDP);
Navg_FDP = HFDP+2.*H2FDP;

dNavgFDPdH = (P_FDP.*(10.^(pKa1_FDP)+2.*10.^(-pH_cy+pKa1_FDP+pKa2_FDP))- ...
        (10.^(-pH_cy+pKa1_FDP)+2.*10.^(-2.*pH_cy+pKa1_FDP+pKa2_FDP)).* ...
    (10.^(pKa1_FDP)+2.*10.^(-pH_cy+pKa1_FDP+pKa2_FDP)))./(c0.*P_FDP.^2);
dNavgFDPdmg = -(10.^(-pH_cy+pKa1_FDP)+2.*10.^(-2.*pH_cy+pKa1_FDP+pKa2_FDP)).*...
    (10.^(pKamg_FDP))./(c0.*P_FDP.^2);
dmgFDPdmg = (P_FDP.*(10.^pKamg_FDP)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_FDP))./c0)./P_FDP.^2;
dmgFDPdpH = (mg./c0).*(10.^pKamg_FDP).* ... 
    (P_FDP.*10.^(-pH_cy+pKa1_FDP).*log(10)- ...
    10.^(-2.*pH_cy+pKa1_FDP+pKa2_FDP).*2.*log(10))./P_FDP.^2;

NH_UFDP = 10;
deltaGof_UFDP = -2601.40;% KJ./mol;
deltaGpof_UFDP = deltaGof_UFDP + NH_UFDP.*log(10).*R.*T.*pH_cy - ...
    IcorrdeltaGpof.*(16-NH_UFDP);

pKa1_GAP = 6.45;        
P_GAP = (1 + 10.^(-pH_cy+pKa1_GAP));
UGAP = 1./P_GAP;
HGAP = UGAP.*10.^(-pH_cy+pKa1_GAP);
Navg_GAP = HGAP;
dNavgGAPdH = 10.^(pKa1_GAP).*(P_GAP-10.^(-pH_cy+pKa1_GAP))./(c0.*P_GAP.^2);
dNavgGAPdmg = 0;
NH_UGAP = 5;
deltaGof_UGAP = -1288.60;% KJ./mol;
deltaGpof_UGAP = deltaGof_UGAP + NH_UGAP.*log(10).*R.*T.*pH_cy - ...
        IcorrdeltaGpof.*(4-NH_UGAP);
    
pKamg_G3P = 1.63; % dimensionless;
deltaH1o_G3P = -3.1; %KJ./mol;
deltaH1_G3P = deltaH1o_G3P + Icorr.*(2.^2+1.^2-1.^2);
pKa1_G3P = 6.22+IcorrpKa.*(2.^2+1.^2-1.^2)+TcorrpKa.*deltaH1_G3P;
P_G3P = (1+10.^(-pH_cy+pKa1_G3P)+(mg./c0).*10.^pKamg_G3P);
UG3P = 1./P_G3P;
HG3P = UG3P.*10.^(-pH_cy+pKa1_G3P);
mgG3P = UG3P.*(mg./c0).*10.^pKamg_G3P;
Navg_G3P = HG3P;
dNavgG3PdH = 10.^(pKa1_G3P).*(P_G3P-10.^(-pH_cy+pKa1_G3P))./(c0.*P_G3P.^2);
dNavgG3Pdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_G3P+pKamg_G3P)./P_G3P.^2;

dmgG3Pdmg = (P_G3P.*(10.^pKamg_G3P)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_G3P))./c0)./P_G3P.^2;
dmgG3PdpH = (mg./c0).*10.^(-pH_cy+pKa1_G3P+pKamg_G3P).*log(10)./P_G3P.^2;
NH_UG3P = 7;
deltaGof_UG3P = -1339.25; %KJ./mol;
deltaGpof_UG3P = deltaGof_UG3P + NH_UG3P.*log(10).*R.*T.*pH_cy - ...
    IcorrdeltaGpof.*(4-NH_UG3P);
                
pKa1_DHAP = 5.9; %dimensionless, 
pKamg_DHAP = 1.57; %dimensionless;
P_DHAP = (1 + 10.^(-pH_cy+pKa1_DHAP)+(mg./c0).*10.^pKamg_DHAP);
UDHAP = 1./P_DHAP;
HDHAP = UDHAP.*(10.^(-pH_cy+pKa1_DHAP));
mgDHAP = UDHAP.*(mg./c0).*(10.^pKamg_DHAP);
Navg_DHAP = HDHAP;
dNavgDHAPdH = 10.^(pKa1_DHAP).*(P_DHAP-10.^(-pH_cy+pKa1_DHAP))./(c0.*P_DHAP.^2);
dNavgDHAPdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_DHAP+pKamg_DHAP)./P_DHAP.^2;

dmgDHAPdmg = (P_DHAP.*(10.^pKamg_DHAP)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_DHAP))./c0)./P_DHAP.^2;
dmgDHAPdpH = (mg./c0).*10.^(-pH_cy+pKa1_DHAP+pKamg_DHAP).*log(10)./P_DHAP.^2;

NH_UDHAP = 5;
deltaGof_UDHAP = -1296.26;% KJ./mol;
deltaGpof_UDHAP = deltaGof_UDHAP + NH_UDHAP.*log(10).*R.*T.*pH_cy -... 
    IcorrdeltaGpof.*(4-NH_UDHAP);

pKa1_13DPG = 7.5;% dimensionless;
P_13DPG = (1 + 10.^(-pH_cy+pKa1_13DPG));
U13DPG = 1./P_13DPG;
H13DPG = U13DPG.*10.^(-pH_cy+pKa1_13DPG);
Navg_13DPG = H13DPG;
dNavg13DPGdH = 10.^(pKa1_13DPG).*(P_13DPG-10.^(-pH_cy+pKa1_13DPG))./...
    (c0.*P_13DPG.^2);
dNavg13DPGdmg = 0;
NH_U13DPG = 4;
deltaGof_U13DPG = -2356.14; %KJ./mol;
deltaGpof_U13DPG = deltaGof_U13DPG + NH_U13DPG.*log(10).*R.*T.*pH_cy -... 
    IcorrdeltaGpof.*(16-NH_U13DPG);

pKa1_3PG = 6.21; %dimensionless;
P_3PG = (1 + 10.^(-pH_cy+6.21));
U3PG = 1./P_3PG;
H3PG = U3PG.*10.^(-pH_cy+6.21);
Navg_3PG = H3PG;
dNavg3PGdH = 10.^(pKa1_3PG).*(P_3PG-10.^(-pH_cy+pKa1_3PG))./(c0.*P_3PG.^2);
dNavg3PGdmg = 0;

NH_U3PG = 4; 
deltaGof_U3PG = -1502.54;% KJ./mol;
deltaGpof_U3PG = deltaGof_U3PG + NH_U3PG.*log(10).*R.*T.*pH_cy - ...
    IcorrdeltaGpof.*(9-NH_U3PG);        

pKa1_2PG = 7;% dimensionless, 
pKamg_2PG = 2.45; %dimensionless,
pKak_2PG = 1.18; %dimensionless;
P_2PG = (1 + 10.^(-pH_cy+pKa1_2PG) + (mg./c0).*10.^pKamg_2PG + ...
    (k./c0).*10.^pKak_2PG);
U2PG = 1./P_2PG;
H2PG = U2PG.*10.^(-pH_cy+pKa1_2PG);
mg2PG = U2PG.*(mg./c0).*10.^pKamg_2PG;
k2PG = U2PG.*(k./c0).*10.^pKak_2PG;
Navg_2PG = H2PG;
dNavg2PGdH = 10.^(pKa1_2PG).*(P_2PG-10.^(-pH_cy+pKa1_2PG))./(c0.*P_2PG.^2);
dNavg2PGdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_2PG+pKamg_2PG)./P_2PG.^2;

dmg2PGdmg = (P_2PG.*(10.^pKamg_2PG)./c0-...
        (mg./c0).*(10.^(2.*pKamg_2PG))./c0)./P_2PG.^2;
dmg2PGdpH = (mg./c0).*10.^(-pH_cy+pKa1_2PG+pKamg_2PG).*log(10)./P_2PG.^2;

NH_U2PG = 4;
deltaGof_U2PG = -1496.38; %KJ./mol;
deltaGpof_U2PG = deltaGof_U2PG + NH_U2PG.*log(10).*R.*T.*pH_cy - ... 
        IcorrdeltaGpof.*(9-NH_U2PG);

pKa1_PEP = 6.35;% dimensionless, 
pKamg_PEP = 2.26;% dimensionless, 
pKak_PEP = 1.08;% dimensionless; 
P_PEP = (1 + 10.^(pKa1_PEP-pH_cy)+ (mg./c0).*10.^pKamg_PEP + ...
    (k./c0).*10.^pKak_PEP);
UPEP = 1./P_PEP;
HPEP = UPEP.*(10.^(pKa1_PEP-pH_cy));
kPEP = UPEP.*(k./c0).*10.^pKak_PEP;
mgPEP = UPEP.*(mg./c0).*(10.^pKamg_PEP);
Navg_PEP = HPEP;
dNavgPEPdH = 10.^(pKa1_PEP).*(P_PEP-10.^(-pH_cy+pKa1_PEP))./(c0.*P_PEP.^2);
dNavgPEPdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_PEP+pKamg_PEP)./P_PEP.^2;

dmgPEPdmg = (P_PEP.*(10.^pKamg_PEP)./c0- ...
        (mg./c0).*(10.^(2.*pKamg_PEP))./c0)./P_PEP.^2;
dmgPEPdpH = (mg./c0).*10.^(-pH_cy+pKa1_PEP+pKamg_PEP).*log(10)./P_PEP.^2;

 NH_UPEP = 2;
 deltaGof_UPEP = -1263.65; %KJ./mol;
deltaGpof_UPEP = deltaGof_UPEP + NH_UPEP.*log(10).*R.*T.*pH_cy -... 
    IcorrdeltaGpof.*(9-NH_UPEP);

pKa1_PYR = 2.49; %dimensionless;
P_PYR = (1 + 10.^(-pH_cy+pKa1_PYR));        
UPYR = 1./P_PYR;
HPYR = UPYR.*10.^(-pH_cy+pKa1_PYR);
Navg_PYR = HPYR;
dNavgPYRdH = 10.^(pKa1_PYR).*(P_PYR-10.^(-pH_cy+pKa1_PYR))./(c0.*P_PYR.^2);
dNavgPYRdmg = 0;
NH_UPYR = 3;
deltaGof_UPYR = -472.27;% KJ./mol;
deltaGpof_UPYR = deltaGof_UPYR + NH_UPYR.*log(10).*R.*T.*pH_cy - ...
    IcorrdeltaGpof.*(1-NH_UPYR);        

pKamg_LAC = 0.98;% dimensionless;
deltaH1o_LAC = -0.33;% KJ./mol;
deltaH1_LAC = deltaH1o_LAC+Icorr.*(1.^2+1.^2-0.^2);
pKa1_LAC = 3.67+IcorrpKa.*(1.^2+1.^2-0.^2)+TcorrpKa.*deltaH1_LAC;
P_LAC = (1 + 10.^(-pH_cy+pKa1_LAC)+ (mg./c0).*10.^pKamg_LAC);
ULAC = 1./P_LAC;
HLAC = ULAC.*10.^(-pH_cy+pKa1_LAC);
mgLAC = (mg./c0).*10.^pKamg_LAC.*ULAC;
Navg_LAC = HLAC;
dNavgLACdH = 10.^(pKa1_LAC).*(P_LAC-10.^(-pH_cy+pKa1_LAC))./(c0.*P_LAC.^2);
dNavgLACdmg = -(10.^(-pH_cy)./c0).*10.^(pKa1_LAC+pKamg_LAC)./P_LAC.^2;

dmgLACdmg = (P_LAC.*(10.^pKamg_LAC)./c0-...
        (mg./c0).*(10.^(2.*pKamg_LAC))./c0)./P_LAC.^2;
dmgLACdpH = (mg./c0).*10.^(-pH_cy+pKa1_LAC+pKamg_LAC).*log(10)./P_LAC.^2;

NH_ULAC = 5; 
deltaGof_ULAC = -516.72;% KJ./mol;
deltaGpof_ULAC = deltaGof_ULAC + NH_ULAC.*log(10).*R.*T.*pH_cy -... 
        IcorrdeltaGpof.*(1-NH_ULAC);

    
%Free energy change in glycosidic linkage being broken (Beard and Qian)
dNH_GLY = -10;% dimensionless; 
deltaGpo_GLY = 655.7 + dNH_GLY.*log(10).*R.*T.*pH_cy;


NH_NAD = 26;
deltaGof_NAD = 0;% KJ./mol;
deltaGpof_NAD = deltaGof_NAD + NH_NAD.*log(10).*R.*T.*pH_cy - ... 
    IcorrdeltaGpof.*(1-NH_NAD);

NH_NADH = 27;
deltaGof_NADH = 22.65;% KJ./mol;
deltaGpof_NADH = deltaGof_NADH + NH_NADH.*log(10).*R.*T.*pH_cy - ...
        IcorrdeltaGpof.*(4-NH_NADH);
        
NH_H2O = 2;
deltaGof_H2O = -237.19;% KJ./mol;
deltaGpof_H2O = deltaGof_H2O + NH_H2O.*log(10).*R.*T.*pH_cy - ...
        IcorrdeltaGpof.*(0-NH_H2O);

NH_NH4 = 4;% dimensionless;
    
NH_H = 1;% dimensionless, 
deltaGof_H = 0;% KJ./mol; 
deltaGpof_H = deltaGof_H+NH_H.*log(10).*R.*T.*pH_cy;

%Calculation of proton binding change and apparent equilibrium constants
 
%For each of the reactions, a reference reaction has to be defined which
%allows for the computation of net proton binding change. 
%The binding change
%comes close to integer values at pH >= 8. 

% creatine kinase: UPCR + ADP3 + H+= ATP4 + UCR
%net proton binding change
deltaH_CK = (Navg_ATP + Navg_CR - Navg_PCR - Navg_ADP) + ...
    (NH_ATP4+NH_HCR - NH_HPCR - NH_ADP3);

Kref_CK = 2.58e8; %CK reference equilibrium constant at T=298.15 K,I=0
deltaHo_CKo = -17.55;% KJ./mol;
deltaH1_CK = deltaHo_CKo + Icorr.*(2.^2+3.^2+1.^2-4.^2-0.^2);
Kref_CKI = exp(log(Kref_CK)+alphadebye.*I.^0.5.* ...
    (2.^2+3.^2+1.^2-4.^2-0.^2)./(1+B.*I.^0.5)); %Ionic strength correction
Kref_CKT = 10.^(log10(Kref_CKI)-TcorrpKa.*deltaH1_CK);
deltaGpo_CK = -R.*T.*log(Kref_CKT);
Kapp_CK = exp(-deltaGpo_CK./(R.*T)).*10.^(-pH_cy).*P_ATP.*P_CR./(P_PCR.*P_ADP);

% adenylate kinase: ATP4 + AMP2 = 2.*ADP3
%proton binding change
deltaH_ADK = 2.*Navg_ADP - Navg_ATP - Navg_AMP + (2.*NH_ADP3-NH_ATP4-NH_AMP2); 
    
deltaGpo_ADK = 2.*deltaGpof_ADP3-deltaGpof_ATP4-deltaGpof_AMP2;
Kapp_ADK = exp(-deltaGpo_ADK./(R.*T)).*P_ADP.^2./(P_ATP.*P_AMP);

% glycogen phosphorylase: glycogen(n) + HPi2 = glycogen(n-1) + UG1P
deltaH_GP = Navg_G1P - Navg_Pi + (1-NH_HPi2);
    
deltaGpo_GP = deltaGpo_GLY + deltaGpof_UG1P -deltaGpof_HPi2;
Kapp_GP = exp(-deltaGpo_GP./(R.*T)).*P_G1P./(P_Pi);
    
% PGLM: UG1P = UG6P
%proton binding change
deltaH_PGLM = Navg_G6P - Navg_G1P + (NH_UG6P-NH_UG1P);
deltaGpo_PGLM = deltaGpof_UG6P-deltaGpof_UG1P;
Kapp_PGLM = exp(-deltaGpo_PGLM./(R.*T)).*P_G6P./P_G1P;
    
% PGI: UG6P = UF6P
%proton binding change
deltaH_PGI= Navg_F6P - Navg_G6P + (NH_UF6P-NH_UG6P);
deltaGpo_PGI = deltaGpof_UF6P - deltaGpof_UG6P;
Kapp_PGI = exp(-deltaGpo_PGI./(R.*T)).*P_F6P./P_G6P;

% PFK: UF6P + ATP4 = UFDP + ADP3 + H+
%proton binding change
deltaH_PFK = Navg_ADP+Navg_FDP-Navg_F6P-Navg_ATP+...
    NH_ADP3+NH_UFDP-NH_UF6P-NH_ATP4;
    
deltaGpo_PFK = deltaGpof_UFDP + deltaGpof_ADP3 + deltaGpof_H ...
        - deltaGpof_UF6P - deltaGpof_ATP4;
Kapp_PFK = exp(-deltaGpo_PFK./(R.*T)).*P_FDP.*P_ADP./(P_F6P.*P_ATP.*10.^(-pH_cy)); 
        
% ALD: UFDP = UDHAP+UGAP
%proton binding change
deltaH_ALD = Navg_DHAP+Navg_GAP - Navg_FDP + ...
        (NH_UDHAP+NH_UGAP-NH_UFDP);
    
deltaGpo_ALD = deltaGpof_UDHAP + deltaGpof_UGAP - deltaGpof_UFDP;
Kapp_ALD = exp(-deltaGpo_ALD./(R.*T)).*P_GAP.*P_FDP./P_DHAP;  

% TPI: UGAP = UDHAP
%proton binding change
deltaH_TPI = Navg_DHAP - Navg_GAP + (NH_UDHAP-NH_UGAP);
    
deltaGpo_TPI = deltaGpof_UDHAP - deltaGpof_UGAP;
Kapp_TPI = exp(-deltaGpo_TPI./(R.*T)).*P_DHAP./P_GAP;

% GAPDH: UGAP + HPi2 + NAD = U13DPG + NADH + H+
%proton binding change
deltaH_GAPDH = Navg_13DPG - Navg_GAP - Navg_Pi +...
        (NH_U13DPG+NH_NADH-NH_UGAP-NH_HPi2-NH_NAD);
    
deltaGpo_GAPDH = deltaGpof_U13DPG+deltaGpof_NADH+deltaGpof_H ...
    -deltaGpof_HPi2-deltaGpof_UGAP-deltaGpof_NAD;
Kapp_GAPDH = exp(-deltaGpo_GAPDH./(R.*T)).*P_13DPG./(P_Pi.*P_GAP.*10.^(-pH_cy));

% G3PDH: UG3P + NAD = UDHAP + NADH + H+
deltaH_G3PDH = Navg_DHAP - Navg_G3P + ...
        (NH_UDHAP+NH_NADH-NH_NAD-NH_UG3P);
    
deltaGpo_G3PDH = deltaGpof_H+deltaGpof_NADH+deltaGpof_UDHAP ...
        -deltaGpof_NAD-deltaGpof_UG3P;
    
Kapp_G3PDH = exp(-deltaGpo_G3PDH./(R.*T)).*P_DHAP./(P_G3P.*10.^(-pH_cy));

% PGK: U13DPG+ADP3 = U3PG+ATP4
%proton binding change
deltaH_PGK = Navg_3PG+Navg_ATP-Navg_13DPG-Navg_ADP +...
        (NH_U3PG+NH_ATP4-NH_U13DPG-NH_ADP3);
    
deltaGpo_PGK = deltaGpof_ATP4+deltaGpof_U3PG ...
    -deltaGpof_U13DPG-deltaGpof_ADP3;
Kapp_PGK = exp(-deltaGpo_PGK./(R.*T)).*P_ATP.*P_3PG./(P_13DPG.*P_ADP);

% PGM: U3PG = U2PG
%proton binding change
deltaH_PGM = Navg_2PG - Navg_3PG + (NH_U2PG-NH_U3PG);
deltaGpo_PGM = deltaGpof_U2PG-deltaGpof_U3PG;
Kapp_PGM = exp(-deltaGpo_PGM./(R.*T)).*P_2PG./P_3PG;

% ENOL: U2PG = UPEP + H20
%proton binding change
deltaH_ENOL = Navg_PEP - Navg_2PG + (NH_H2O+NH_UPEP-NH_U2PG);
deltaGpo_ENOL = deltaGpof_H2O+deltaGpof_UPEP-deltaGpof_U2PG;
Kapp_ENOL = exp(-deltaGpo_ENOL./(R.*T)).*P_PEP./P_2PG;

% PK: UPEP+ADP3+H+ = UPYR+ATP4
%proton binding change
deltaH_PK = Navg_PYR+Navg_ATP-Navg_PEP-Navg_ADP+ ...
        (NH_UPYR+NH_ATP4-NH_UPEP-NH_ADP3);
deltaGpo_PK = deltaGpof_UPYR+deltaGpof_ATP4 ...
    -deltaGpof_H-deltaGpof_UPEP-deltaGpof_ADP3;
Kapp_PK = exp(-deltaGpo_PK./(R.*T)).*P_PYR.*P_ATP.*10.^(-pH_cy)./(P_PEP.*P_ADP);

% LDH: UPYR + NADH + H+= ULAC + NAD
%proton binding change
deltaH_LDH = Navg_LAC-Navg_PYR+(NH_ULAC+NH_NAD-NH_UPYR-NH_NADH);
deltaGpo_LDH = deltaGpof_ULAC+deltaGpof_NAD ...
        -deltaGpof_UPYR-deltaGpof_NADH-deltaGpof_H;
Kapp_LDH = exp(-deltaGpo_LDH./(R.*T)).*P_LAC.*10.^(-pH_cy)./P_PYR;

%ATPase: ATP4 + H2O = ADP3 + HPi2 + H+
%proton binding change
% dimensionless;
deltaH_ATPase = Navg_ADP+Navg_Pi-Navg_ATP + (NH_ADP3+NH_HPi2-NH_ATP4-NH_H2O);

deltaGpo_ATPase = deltaGpof_ADP3+deltaGpof_HPi2+deltaGpof_H ...
        -deltaGpof_H2O-deltaGpof_ATP4;

Kapp_ATPase = exp(-deltaGpo_ATPase./(R.*T)).*P_ADP.*P_Pi./(P_ATP.*10.^(-pH_cy));    

% AMP Deaminase
deltaH_AMPDA = Navg_IMP-Navg_AMP + (NH_IMP2+NH_NH4-NH_H2O-NH_AMP2);

%------------------------------------------------------------------------------------------------------
%KINETICS
%------------------------------------------------------------------------------------------------------

% FLUX expressions


% glycogen phosphorylase
%  Vfgly = 5.000000e-002 ;
Vfgly = Par(1);
%  fracA = 2.000000e-003 ;
expno= Par(100);
fracA = Par(2);
%  KgpA_glyf = 1.700000e-003 ;
KgpA_glyf = Par(3);
%  KgpA_pi = 4.000000e-003 ;
KgpA_pi = Par(4);
%  KgpA_igly = 2.000000e-003 ;
KgpA_igly = Par(5);
%  KgpA_ipi = 4.700000e-003 ;
KgpA_ipi = Par(6);
%  KgpA_glyb = 1.500000e-004 ;
KgpA_glyb = Par(7);
%  KgpA_g1p = 2.700000e-003 ;
KgpA_g1p = Par(8);
%  KgpA_ig1p = 1.010000e-002 ;
KgpA_ig1p = Par(9);

Dglya = 1 + gly./KgpA_glyf + Pi./KgpA_pi + gly.*Pi./...
    (KgpA_glyf.*KgpA_ipi) + gly./KgpA_glyb + G1P./KgpA_g1p... 
    + gly.*G1P./(KgpA_ig1p.*KgpA_glyb);
pa = 1.404./(1+10.^(5.94-pH_cy)+10.^(pH_cy-7.29));
VbglyA = pa.*Vfgly.*KgpA_glyb.*KgpA_ig1p./(KgpA_igly.*KgpA_pi.*Kapp_GP);
glyAF = (pa.*Vfgly.*Pi./(KgpA_igly.*KgpA_pi))./Dglya;
glyAR= (VbglyA.*gly./(KgpA_glyb.*KgpA_ig1p))./Dglya;
flux_GPa = fracA.*(gly.*glyAF - G1P.*glyAR);
        
%Glycogen Phosphorylase B
fracB = 1 - fracA;
%  KgpB_pi = 2.000000e-004 ;
KgpB_pi = Par(10);
%  KgpB_ipi = 4.600000e-003 ;
KgpB_ipi = Par(11);
%  KgpB_iglyf = 1.500000e-002 ;
KgpB_iglyf = Par(12);
%  KgpB_g1p = 1.500000e-003 ;
KgpB_g1p = Par(13);
%  KgpB_ig1p = 7.400000e-003 ;
KgpB_ig1p = Par(14);
%  KgpB_iglyb = 4.400000e-003 ;
KgpB_iglyb = Par(15);
%  Kgp_amp = 9.700000e-005 ;
Kgp_amp = Par(16);
%  interactioncoeff = 2.000000e-002 ;
interactioncoeff = Par(17);
%  nH = 1.750000e+000 ;
nH = Par(18);
% allosteric term for AMP activation of b
M = (((AMP./Kgp_amp).^nH)./interactioncoeff)./...
    (1 + ((AMP./Kgp_amp).^nH)./interactioncoeff);

Dglyb = (1 + gly./KgpB_iglyf + Pi./KgpB_ipi + gly./KgpB_iglyb ...
    + G1P./KgpB_ig1p + gly.*Pi./(KgpB_iglyf.*KgpB_pi) + gly.*G1P./ ...
     (KgpB_g1p.*KgpB_iglyb));
pb = 1.75./(1+10.^(6.12-pH_cy)+10.^(pH_cy-7.03));
VbglyB = pb.*Vfgly.*KgpB_g1p.*KgpB_iglyb./(KgpB_iglyf.*KgpB_pi.*Kapp_GP);
glyBF = pb.*M.*(Vfgly.*Pi./(KgpB_iglyf.*KgpB_pi))./Dglyb;
glyBR = M.*(VbglyB.*gly./(KgpB_g1p.*KgpB_iglyb))./Dglyb;
flux_GPb = fracB.*(gly.*glyBF - G1P.*glyBR);

% PGLM: glucose-1-phosphate = glucose-6-phosphate
%  Vffpglm  = 4.800000e-001 ;
Vffpglm = Par(19);
%  Kpglm_g1p = 6.300000e-005 ;
Kpglm_g1p = Par(20);
%  Kpglm_g6p = 3.000000e-005 ;
Kpglm_g6p = Par(21);
 
Vfpglm = Vffpglm.*(1.329)./(1+10.^(-pH_cy+6.64)+10.^(pH_cy-8.36));
Vbpglm = Vfpglm.*Kpglm_g6p./(Kpglm_g1p.*Kapp_PGLM);
v_PGLM = (Vfpglm.*G1P./Kpglm_g1p - Vbpglm.*G6P./Kpglm_g6p)./...
    (1 + G1P./Kpglm_g1p + G6P./Kpglm_g6p);
        
% PGI: glucose-6-phosphate = fructose-6-phosphate
%  Vbbpgi = 8.800000e-001 ;
Vbbpgi = Par(22);
%  Kpgi_g6p = 4.800000e-004 ;
Kpgi_g6p = Par(23);
%  Kpgi_f6p = 1.190000e-004 ;
Kpgi_f6p = Par(24);
Vbpgi = Vbbpgi./(1+10.^(-pH_cy+6.94)+10.^(pH_cy-9.35));
Vfpgi = Vbpgi.*Kpgi_g6p./Kpgi_f6p.*Kapp_PGI;
v_PGI = (Vfpgi.*G6P./Kpgi_g6p - Vbpgi.*F6P./Kpgi_f6p)./ ...
     (1+F6P./Kpgi_f6p + G6P./Kpgi_g6p);     
        
% PFK: Fructose-6-phosphate + ATP = fructose 1,6-bisphosphate + ADP
%  Vffpfk  = 5.600000e-002 ;
Vffpfk = Par(25);
%  Kpfk_f6p  = 1.800000e-004 ;
Kpfk_f6p = Par(26);
%  Kpfk_f6pT = 2.000000e-002 ;
Kpfk_f6pT = Par(27);
%  Kpfk_atp = 8.000000e-005 ;
Kpfk_atp = Par(28);
%  Kpfk_atpT = 2.500000e-004 ;
Kpfk_atpT = Par(29);
%  Kpfk_fbp = 4.020000e-003 ;
Kpfk_fbp = Par(30);
%  Kpfk_fbpT = 4.020000e-003 ;
Kpfk_fbpT = Par(31);
%  Kpfk_adp = 2.700000e-003 ;
Kpfk_adp = Par(32);
%  Kpfk_adpT = 2.700000e-003 ;
Kpfk_adpT = Par(33);
%  Kpfki = 8.700000e-004 ;
Kpfki = Par(34);
%  Kmpfk = 6.000000e-005 ;
Kmpfk = Par(35);
%  d = 1.000000e-002 ;
d = Par(36);
%  e = 1.000000e-002 ;
e = Par(37);
%  Lo = 13 ;
Lo = Par(38);

Vfpfk = Vffpfk./(1+(pH_cy./6.8).^(-30));
Vbpfk = Vfpfk.*Kpfk_fbp.*Kpfk_adp./(Kpfk_f6p.*Kpfk_atp.*Kapp_PFK);

L = Lo.*((1 + ATP./Kpfki)./(1 + d.*ATP./Kpfki).* ...
    (1 + e.*AMP./Kmpfk)./(1 + AMP./Kmpfk)).^4;
alpha = Kpfk_f6p.*Kpfk_atp./(Kpfk_f6pT.*Kpfk_atpT);
Delta = (1+F6P./Kpfk_f6p).*(1+ATP./Kpfk_atp) + FBP./Kpfk_fbp + ...
    (ADP./Kpfk_adp).*(1+FBP./Kpfk_fbp);
Deltap = (1+F6P./Kpfk_f6pT).*(1+ATP./Kpfk_atpT) + FBP./Kpfk_fbpT + ...
    (ADP./Kpfk_adpT).*(1+FBP./Kpfk_fbpT); 
v_PFK = (Vfpfk.*F6P.*ATP./(Kpfk_f6p.*Kpfk_atp)./Delta- ...
    Vbpfk.*ADP.*FBP./(Kpfk_adp.*Kpfk_fbp)./Delta).* ...
    (1+alpha.*L.*(Deltap./Delta).^3)./(1+L.*(Deltap./Delta).^4);        

% ALD: Fructose 1,6-bisphosphate = dihydroxyacetone phosphate+glyceraldehyde phosphate
%Aldolase
%  Vffald = 1.040000e-001 ;
Vffald = Par(39);
%  Kald_fbp = 5.000000e-005 ;
Kald_fbp = Par(40);
%  Kald_dhap = 2.000000e-003 ;
Kald_dhap = Par(41);
%  Kald_gap = 1.000000e-003 ;
Kald_gap = Par(42);

Vfald = Vffald.*(1.013)./(1+10.^(-pH_cy+5.32)+10.^(pH_cy-9.15));
Vbald = Vfald.*Kald_gap.*Kald_dhap./(Kald_fbp.*Kapp_ALD);
v_ALD = (Vfald.*FBP./Kald_fbp-Vbald.*GAP.*DHAP./(Kald_gap.*Kald_dhap))./ ...
(1+FBP./Kald_fbp + GAP./Kald_gap + DHAP./Kald_dhap);        
        
% TPI: glyceraldehyde phosphate = dihydoxyacetone phosphate
%  Vfftpi  = 12 ;
Vfftpi = Par(43);
%  Ktpi_gap = 3.200000e-004 ;
Ktpi_gap = Par(44);
%  Ktpi_dhap = 6.100000e-004 ;
Ktpi_dhap = Par(45);
Vftpi=  Vfftpi;
Vbtpi = Vftpi.*Ktpi_dhap./(Ktpi_gap.*(Kapp_TPI));
v_TPI = (Vftpi.*GAP./Ktpi_gap-Vbtpi.*DHAP./Ktpi_dhap)./ ...
    (1 + GAP./Ktpi_gap + DHAP./Ktpi_dhap);
        
%G3PDH
%  Vbbg3pdh = 8.250000e-002 ;
Vbbg3pdh = Par(46);
%  Kg3pdh_g3p = 1.800000e-004 ;
Kg3pdh_g3p = Par(47);
%  Kg3pdh_nad = 1.200000e-005 ;
Kg3pdh_nad = Par(48);
%  Kg3pdh_dhap = 2.200000e-004 ;
Kg3pdh_dhap = Par(49);
%  Kg3pdh_nadh = 8.000000e-006 ;
Kg3pdh_nadh = Par(50);

Dg3pdh = (1 + G3P./Kg3pdh_g3p + NADH./Kg3pdh_nadh).*...
    (1 + DHAP./Kg3pdh_dhap + NAD./Kg3pdh_nadh);
Vbg3pdh  = Vbbg3pdh;
Vfg3pdh = Vbg3pdh.*Kg3pdh_g3p.*Kg3pdh_nad.*Kapp_G3PDH./ ...
    (Kg3pdh_dhap.*Kg3pdh_nadh);
v_G3PDH = (Vfg3pdh.*G3P.*NAD./(Kg3pdh_g3p.*Kg3pdh_nad)- ...
    Vbg3pdh.*DHAP.*NADH./(Kg3pdh_dhap.*Kg3pdh_nadh))./Dg3pdh;

% GAPDH: glyceraldehyde phosphate + Pi + NAD = 1,3-bisphospholycerate + NADH
%  Vffgad = 1.265000e+000 ;
Vffgad = Par(51);
%  Kgapdh_gap = 2.500000e-006 ;
Kgapdh_gap = Par(52);
%  Kgapdh_nad = 9.000000e-005 ;
Kgapdh_nad = Par(53);
%  Kgapdh_pi = 2.900000e-004 ;
Kgapdh_pi = Par(54);
%  Kgapdh_bpg = 8.000000e-007 ;
Kgapdh_bpg = Par(55);
%  Kgapdh_nadh = 3.300000e-006 ;
Kgapdh_nadh = Par(56);

% denominator for GAPDH
Dgap = 1 + Pi./Kgapdh_pi + GAP./Kgapdh_gap + NAD./Kgapdh_nad + ... 
      GAP.*NAD./(Kgapdh_gap.*Kgapdh_nad) + GAP.*NAD.*Pi./(Kgapdh_gap.* ...
      Kgapdh_nad.*Kgapdh_pi) + BPG./Kgapdh_bpg + NADH./Kgapdh_nadh ...
      + BPG.*NADH./(Kgapdh_nadh.*Kgapdh_bpg);

Vfgad = Vffgad.*(0.0007.*exp(pH_cy.*0.8979));
Vbgad = Vfgad.*Kgapdh_bpg.*Kgapdh_nadh./ ...
    (Kgapdh_gap.*Kgapdh_pi.*Kgapdh_nad.*Kapp_GAPDH);
v_GAPDH= (Vfgad.*GAP.*NAD.*Pi./(Kgapdh_nad.*Kgapdh_gap.*Kgapdh_pi)- ...
    Vbgad.*BPG.*NADH./(Kgapdh_bpg.*Kgapdh_nadh))./Dgap;

    
% PGK: 1,3-Bisphosphoglycerate+ADP = 3-phosphoglycerate+ATP
%  Vbbpgk = 1.120000e+000 ;
Vbbpgk = Par(57);
%  Kpgk_bpg = 2.000000e-003 ;
Kpgk_bpg = Par(58);
%  Kpgk_adp = 8.000000e-006 ;
Kpgk_adp = Par(59);
%  Kpgk_3pg = 1.200000e-003 ;
Kpgk_3pg = Par(60);
%  Kpgk_atp = 3.500000e-004 ;
Kpgk_atp = Par(61);
Vbpgk  = Vbbpgk;
Vfpgk = Vbpgk.*Kpgk_bpg.*Kpgk_adp./(Kpgk_3pg.*Kpgk_atp).*Kapp_PGK;
D_PGK=(1 + ADP./Kpgk_adp + BPG./Kpgk_bpg + BPG.*ADP./(Kpgk_bpg.*Kpgk_adp)+ ...
    P3G./Kpgk_3pg + ATP./Kpgk_atp + P3G.*ATP./(Kpgk_3pg.*Kpgk_atp));
v_PGK = (Vfpgk.*BPG.*ADP./(Kpgk_adp.*Kpgk_bpg)- ...
    Vbpgk.*ATP.*P3G./(Kpgk_atp.*Kpgk_3pg))./D_PGK;
        
% PGM: 3-phophoglycerate = 2-phosphoglycerate
%  Vffpgm = 1.120000e+000 ;
Vffpgm = Par(62);
%  Kpgm_3pg = 2.000000e-004 ;
Kpgm_3pg = Par(63);
%  Kpgm_2pg = 1.400000e-005 ;
Kpgm_2pg = Par(64);
Vfpgm = Vffpgm.*(0.989)./(1+10.^(-pH_cy+5.62)+10.^(pH_cy-8.74));
Vbpgm = Vfpgm.*Kpgm_2pg./(Kpgm_3pg.*Kapp_PGM);
v_PGM = (Vfpgm.*P3G./Kpgm_3pg-Vbpgm.*P2G./Kpgm_2pg)./...
    (1+P3G./Kpgm_3pg + P2G./Kpgm_2pg);
        
% ENOL: 2-phosphoglycerate = phosphoenolpyruvate+H20
%Enolase
%  Vffen = 1.920000e-001 ;
Vffen = Par(65);
%  Ken_2pg = 1.000000e-004 ;
Ken_2pg = Par(66);
%  Ken_pep = 3.700000e-004 ;
Ken_pep = Par(67);
Vfen =  Vffen;
Vben = Vfen.*Ken_pep./(Ken_2pg.*Kapp_ENOL);
v_ENOL = (Vfen.*P2G./Ken_2pg-Vben.*PEP./Ken_pep)./ ...
    (1 +PEP./Ken_pep + P2G./Ken_2pg);

% PK: phosphoenolpyruvate+ADP = pyruvate+ATP
%  Vffpk = 1.440000e+000 ;
Vffpk = Par(68);
%  Kpk_pep = 8.000000e-005 ;
Kpk_pep = Par(69);
%  Kpk_adp = 3.000000e-004 ;
Kpk_adp = Par(70);
%  Kpk_pyr = 7.050000e-003 ;
Kpk_pyr = Par(71);
%  Kpk_atp = 1.130000e-003 ;
Kpk_atp = Par(72);

Vfpk = Vffpk.*(1.05)./(1+10.^(-pH_cy+5.58)+10.^(pH_cy-8.79));
Vbpk = Vfpk.*Kpk_pyr.*Kpk_atp./(Kpk_pep.*Kpk_adp.*Kapp_PK);
v_PK = (Vfpk.*PEP.*ADP./(Kpk_pep.*Kpk_adp)- ...
    Vbpk.*PYR.*ATP./(Kpk_pyr.*Kpk_atp))./ ...
(1+PEP./Kpk_pep+ADP./Kpk_adp + PEP.*ADP./(Kpk_pep.*Kpk_adp) ...
+ ATP./Kpk_atp + PYR./Kpk_pyr + PYR.*ATP./(Kpk_pyr.*Kpk_atp));
        
% LDH: Pyruvate + NADH = Lactate + NAD
%  Vffldh = 1.920000e+000 ;
Vffldh = Par(73);
%  Kldh_pyr = 3.350000e-004 ;
Kldh_pyr = Par(74);
%  Kldh_nadh  = 2.000000e-006 ;
Kldh_nadh = Par(75);
%  Kldh_lac = 1.700000e-002 ;
Kldh_lac = Par(76);
%  Kldh_nad = 8.490000e-004 ;
Kldh_nad = Par(77);
Vfldh = Vffldh.*(-0.1134.*pH_cy+1.6069);
Vbldh = Vfldh.*Kldh_lac.*Kldh_nad./(Kldh_pyr.*Kldh_nadh.*Kapp_LDH);
v_LDH = (Vfldh.*PYR.*NADH./(Kldh_pyr.*Kldh_nadh)- ...
    Vbldh.*LAC.*NAD./(Kldh_lac.*Kldh_nad))./...
(1 + PYR./Kldh_pyr +  NADH./Kldh_nadh + ...
PYR.*NADH./(Kldh_pyr.*Kldh_nadh) + LAC./Kldh_lac + NAD./Kldh_nad + ...
LAC.*NAD./(Kldh_lac.*Kldh_nad));        

%ATPase: ATP + H2O = ADP + Pi
%  VmaxATPase = 0 ;
VmaxATPase = Par(88);
%  Katp_ATPase = 1.000000e-004 ;
Katp_ATPase = Par(89);
ATPase = VmaxATPase.*ATP./(Katp_ATPase+ATP);

%AMP deaminase: AMP+H2O = IMP + NH3
%  Vmax_AMPDA = 4e-3;% M./min;
Vmax_AMPDA = Par(101);
%  Kamp_AMPDA = 1e-3;% M;
Kamp_AMPDA = Par(102);
v_AMPDA = Vmax_AMPDA.*AMP./(Kamp_AMPDA+AMP);
        
% creatine kinase: PCr + ADP = ATP + Cr
%  VforCK = 5.000000e-001 ;
VforCK = Par(78);
%  Kck_pcr = 1.110000e-003 ;
Kck_pcr = Par(79);
%  Kck_iatp = 3.500000e-003 ;
Kck_iatp = Par(80);
%  Kck_iadp = 1.350000e-004 ;
Kck_iadp = Par(81);
%  Kck_ipcr = 3.900000e-003 ;
Kck_ipcr = Par(82);
%  Kck_cr = 3.800000e-003 ;
Kck_cr = Par(83);

VrevCK = (VforCK./Kapp_CK).*(Kck_iatp.*Kck_cr./(Kck_iadp.*Kck_pcr));
CK = (VrevCK.*ATP.*Cr./(Kck_iatp.*Kck_cr)-(VforCK.*ADP.*PCr./...
  (Kck_iadp.*Kck_pcr)))./(1+ADP./Kck_iadp+PCr./Kck_ipcr+PCr.*ADP./...
  (Kck_iadp.*Kck_pcr)+ATP./Kck_iatp+Cr.*ATP./(Kck_cr.*Kck_iatp));

    
% adenylate kinase: ATP + AMP = 2.*ADP
%  Vfadk = 8.800000e-001 ;
Vfadk = Par(84);
%  Kadk_amp = 3.200000e-004 ;
Kadk_amp = Par(85);
%  Kadk_atp = 2.700000e-004 ;
Kadk_atp = Par(86);
%  Kadk_adp = 3.500000e-004 ;
Kadk_adp = Par(87);

Vbadk = Vfadk.*Kadk_adp.^2./(Kadk_amp.*Kadk_atp.*Kapp_ADK);
ADK=(Vfadk.*ATP.*AMP./(Kadk_atp.*Kadk_amp) - Vbadk.*(ADP./...
  Kadk_adp).^2)./(1 + ATP./Kadk_atp + AMP./Kadk_amp + ATP.*AMP./...
  (Kadk_atp.*Kadk_amp) + 2.*ADP./Kadk_adp + ADP.^2./Kadk_adp.^2);

%Net proton flux for the pathway, protons consumed = deltaH_rxn.*flux_rxn
%--------------------------------------------------------------------------


%Buffer capacity calculation
%   carnosine = 0.03; %M
carnosine = Par(92);
%   tris = 0.015; %M
tris = Par(91);
%   acetate = 0.01; %M
acetate = Par(93);
%   histidine = 0.03; %M
histidine = Par(103);
    
%Total buffer capacity including resting buffer capacity, inorganic phosphate
%and that due to phosphates on glycolytic intermediates, considering only those
%groups whose pKa's are close to the physiological pH 
bufcapfixed = log(10).*(carnosine).*(10.^(-pH_cy-6.87))./ ...
((10.^(-pH_cy)+ 10.^(-6.87)).^2)+ ...
    log(10).*(tris).*(10.^(-pH_cy-8.3))./((10.^(-pH_cy)+ 10.^(-8.3)).^2) +...
    log(10).*(histidine).*(10.^(-pH_cy-6.3))./((10.^(-pH_cy)+ 10.^(-6.3)).^2)+...
        log(10).*(acetate).*(10.^(-pH_cy-4.8))./((10.^(-pH_cy)+ 10.^(-4.8)).^2);

bufcapmetab = log(10).*10.^(-pH_cy).*c0.*( 1 + dNavgPidH.*Pi+ ...
dNavgATPdH.*ATP+dNavgADPdH.*ADP+dNavgAMPdH.*AMP+dNavgIMPdH.*IMP+dNavgPCRdH.*PCr+ ...
dNavgCRdH.*Cr+dNavgG1PdH.*G1P+dNavgG6PdH.*G6P+dNavgF6PdH.*F6P+ ...
dNavgFDPdH.*FBP+dNavgGAPdH.*GAP+dNavgDHAPdH.*DHAP+dNavgG3PdH.*G3P+ ...
dNavg13DPGdH.*BPG+dNavg3PGdH.*P3G+dNavg2PGdH.*P2G+dNavgPEPdH.*PEP+ ...
dNavgPYRdH.*PYR+dNavgLACdH.*LAC);

protons_consumed = deltaH_CK.*(-CK)+deltaH_ADK.*(ADK)+ ...
deltaH_ATPase.*(ATPase) + deltaH_AMPDA.*v_AMPDA+ ...
deltaH_PGLM.*v_PGLM+(deltaH_GP).*(flux_GPa+flux_GPb)+deltaH_PGI.*v_PGI+ ...
deltaH_PFK.*v_PFK+ deltaH_ALD.*v_ALD+ ...
deltaH_TPI.*v_TPI+deltaH_GAPDH.*v_GAPDH+ ...
deltaH_PGK.*v_PGK+deltaH_PGM.*v_PGM+ ...
deltaH_ENOL.*v_ENOL+deltaH_PK.*v_PK+ ...
deltaH_LDH.*v_LDH+ deltaH_G3PDH.*v_G3PDH;

CKprtflux = deltaH_CK.*(-CK);

glycprtflux = deltaH_PGLM.*v_PGLM+(deltaH_GP).*(flux_GPa+flux_GPb)+ ...
deltaH_PGI.*v_PGI+ ...
deltaH_PFK.*v_PFK+ deltaH_ALD.*v_ALD+ ...
deltaH_TPI.*v_TPI+deltaH_GAPDH.*v_GAPDH+ ...
deltaH_PGK.*v_PGK+deltaH_PGM.*v_PGM+ ...
deltaH_ENOL.*v_ENOL+deltaH_PK.*v_PK+ ...
deltaH_LDH.*v_LDH+ deltaH_G3PDH.*v_G3PDH;

pHODEterm1 = (protons_consumed)./(bufcapfixed+bufcapmetab);

pHODEterm2 = (dNavgPidmg.*Pi+ ...
dNavgATPdmg.*ATP+dNavgADPdmg.*ADP+dNavgAMPdmg.*AMP+dNavgIMPdmg.*IMP+...
dNavgPCRdmg.*PCr+ ...
dNavgCRdmg.*Cr+dNavgG1Pdmg.*G1P+dNavgG6Pdmg.*G6P+dNavgF6Pdmg.*F6P+ ...
dNavgFDPdmg.*FBP+dNavgGAPdmg.*GAP+dNavgDHAPdmg.*DHAP+dNavgG3Pdmg.*G3P+ ...
dNavg13DPGdmg.*BPG+dNavg3PGdmg.*P3G+dNavg2PGdmg.*P2G+dNavgPEPdmg.*PEP+ ...
dNavgPYRdmg.*PYR+dNavgLACdmg.*LAC)./(bufcapfixed+bufcapmetab);

denom_mgODE = -1 - (dmgATP2dmg.*ATP+dmgADPdmg.*ADP+dmgAMPdmg.*AMP+ ...
dmgPidmg.*Pi+dmgPCRdmg.*PCr+ dmgIMPdmg.*IMP+...
dmgG1Pdmg.*G1P+dmgFDPdmg.*FBP+dmgG3Pdmg.*G3P+dmgDHAPdmg.*DHAP+ ...
dmgPEPdmg.*PEP+dmg2PGdmg.*P2G+dmgLACdmg.*LAC);

RHSterm1_mgODE = (dmgATP2dpH.*ATP+dmgADPdpH.*ADP+dmgAMPdpH.*AMP+dmgIMPdpH.*IMP+ ...
dmgPidpH.*Pi+dmgPCRdpH.*PCr+ ...
dmgG1PdpH.*G1P+dmgFDPdpH.*FBP+dmgG3PdpH.*G3P+dmgDHAPdpH.*DHAP+ ...
dmgPEPdpH.*PEP+dmg2PGdpH.*P2G+dmgLACdpH.*LAC)./denom_mgODE;

denomMgpHODE = 1-RHSterm1_mgODE.*pHODEterm2;

%------------------------------------------------------------------------
%DIFFERENTIAL EQUATIONS
dPCrdt =  CK;
dCrdt = -CK;
dNADdt = -v_GAPDH-v_G3PDH+v_LDH;
dNADHdt =  v_GAPDH+v_G3PDH-v_LDH;
dATPdt =  -CK-ADK-v_PFK+v_PGK+v_PK-ATPase;
dADPdt =  CK+2.*ADK+v_PFK-v_PGK-v_PK+ATPase;
dAMPdt = -ADK-v_AMPDA;
dIMPdt = v_AMPDA;
dPidt =  -(flux_GPa+flux_GPb) - v_GAPDH + ATPase;
dglydt =  -(flux_GPa+flux_GPb);
dG1Pdt = (flux_GPa+flux_GPb) - v_PGLM;
dG6Pdt = v_PGLM-v_PGI;
dF6Pdt =  v_PGI-v_PFK;
dFBPdt = v_PFK-v_ALD;   
dDHAPdt = v_ALD+v_TPI+v_G3PDH;
dG3Pdt = -v_G3PDH;
dGAPdt = v_ALD-v_TPI-v_GAPDH;
dBPGdt =  v_GAPDH-v_PGK;
dP3Gdt = v_PGK-v_PGM;
dP2Gdt =  v_PGM-v_ENOL;
dPEPdt =  v_ENOL-v_PK;
dPYRdt = v_PK-v_LDH;    
dLACdt = v_LDH;


RHSterm2_mgODE = (mgATP2.*dATPdt+mgADP.*dADPdt+mgAMP.*dAMPdt+mgIMP.*dIMPdt+ ...
mgPi.*dPidt+mgPCR.*dPCrdt+ ...
mgG1P.*dG1Pdt+mgFDP.*dFBPdt+mgG3P.*dG3Pdt+mgDHAP.*dDHAPdt+ ...
mgPEP.*dPEPdt+mg2PG.*dP2Gdt+mgLAC.*dLACdt)./denom_mgODE;
fixmg = 1;
dmgdt = fixmg.*(RHSterm2_mgODE+RHSterm1_mgODE.*pHODEterm1)./ ...
        (1-RHSterm1_mgODE.*pHODEterm2);

%  fixpH = 1;
fixpH = Par(99);
dpH_calcdt = fixpH.*(pHODEterm1+RHSterm2_mgODE.*pHODEterm2)./ ...
    (1-RHSterm1_mgODE.*pHODEterm2);
dprotonloaddt = -protons_consumed;

dConcdt = [dPCrdt;
dCrdt;
dNADdt;
dNADHdt;
dATPdt;
dADPdt;
dAMPdt;
dIMPdt;
dPidt;
dglydt;
dG1Pdt;
dG6Pdt;
dF6Pdt;
dFBPdt; 
dDHAPdt;
dG3Pdt;
dGAPdt;
dBPGdt;
dP3Gdt;
dP2Gdt;
dPEPdt;
dPYRdt; 
dLACdt;
dmgdt;
dpH_calcdt;
dprotonloaddt];
