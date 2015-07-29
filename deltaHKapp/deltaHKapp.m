pH_cy = 5.5:0.01:8.5;

%GLOBAL PARAMETERS
 R = 8.314e-3;% KJ./(K.*mol); %Ideal gas constant
 T1 = 298.15;% K; %Temp at which physical constants are reported
%   T = 303.15;% K; %Experimental temperature   
T = 310.15;
%   I = 0.1;%0.17;% M; %Ionic strength
I = 0.1;
    
%------------------------------------------------------------------------------------------------
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
mgT = 0.005;
%  k = 0.08;% total K+
k = 0.08;
mg = 5.132658807e-4;
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

% AMP Deaminase: AMP2 + H2O + H+ = IMP2 + NH4
deltaH_AMPDA = Navg_IMP-Navg_AMP + (NH_IMP2+NH_NH4-NH_H2O-NH_AMP2);

%%
%Figure 1
figure(1); clf; set(gca,'Fontsize',14)
plot(pH_cy,deltaH_CK,'k-',pH_cy,deltaH_GAPDH,'k--',...
    pH_cy,deltaH_PGK,'k-.',pH_cy,deltaH_LDH,'k:','linewidth',1.5); 
xlabel('pH'); ylabel('Proton consumption stoichiometry'); box on;
legend('CK','GAPDH','PGK','LDH','Location','SouthEast')
% print -f1 -dtiff -r1200 'figure1.tif'

%%
%Figure 2
figure(2); clf; set(gca,'Fontsize',14)
semilogy(pH_cy,Kapp_CK,'k-',pH_cy,Kapp_GAPDH,'k--',...
    pH_cy,Kapp_PGK,'k-.',pH_cy,Kapp_LDH,'k:','linewidth',1.5); 
xlabel('pH'); ylabel('Apparent equilibrium constant'); box on;
legend('CK','GAPDH','PGK','LDH','Location','SouthEast')
% print -f2 -dtiff -r1200 'figure2.tiff'
