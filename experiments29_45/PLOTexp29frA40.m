
%Plots for GPa fraction 0.4
tdata40 = [ 0 0.041666667 0.083333333 0.125 0.166666667 0.333333333 ...
    0.5 0.666666667 0.833333333 1 1.5 2 3 5 10 ];
% PCr
PCrdata40 = [ 0 0.002686022 0.003622001 0.003855996 0.005259965 ...
    0.0127478 0.016725712 0.018363676 0.018539172 0.018480514 ...
    0.019299459 0.021405316 0.023335686 0.02386215 0.024271622 ];
% Pi
Pidata40 = [ 0.03 0.026728991 0.024623038 0.021288611 0.018714668 ...
    0.012396807 0.008301898 0.007073425 0.006488438 0.006488826 ...
    0.005377402 0.003154552 0.001458167 2.88246E-4 0 ];
% Lactate
LACdata40 = [ 0 5.35795E-5 5.35795E-5 9.3106E-4 0.002335029 0.007892407 ...
    0.009764365 0.010700345 0.010875841 0.010759037 0.011343998 ...
    0.014327296 0.015614208 0.016433153 0.016784129 ];
% HMP
HMPdata40 = [ 0 7.08194E-4 0.001433791 0.00146728 0.001277509 ...
    4.96097E-4 6.07727E-4 9.03547E-4 0.00121053 0.001405883 ...
    0.00184682 8.03084E-4 2.67262E-4 1.55632E-4 1.66795E-4 ];

% FDP
FDPdata40 = [ 0 1.89114E-4 8.86803E-4 0.001405883 0.002137061 ...
    0.001729611 4.29119E-4 8.86467E-5 3.28316E-5 1.60871E-5 ...
    1.48555E-6 6.03231E-5 7.8815E-5 9.56257E-5 9.73068E-5 ];

%Figure 8
figure(8); clf; set(gca,'Fontsize',14)
plot(t40,1000.*PCr40,'k--',tdata40,1000.*PCrdata40,...
    'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('PCr, Pi, Lactate (mM)'); box on; hold on; 
plot(t40,1000.*Pi40,'k-',tdata40,1000.*Pidata40,...
    'k^','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8);
plot(t40,1000.*LAC40,'k-.',tdata40,1000.*LACdata40,...
    'ks','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
hold off;
legend('PCr','PCr data','Pi','Pi data','Lactate', 'Lactate data','Location','Best')
% print -f8 -dtiff -r1200 'figure8.tiff'
fig8 = figure(8);
resp8 = fig2plotly(fig8, 'filename', 'fig8', 'strip', false);

%figure 9
figure(9); clf; set(gca,'Fontsize',14)
plot(t40,1000.*HMP40,'k--',tdata40,1000.*HMPdata40,...
    'ko','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); 
ylabel('Hexose monophosphates, Fructose diphosphates (mM)'); 
box on; hold on; 
plot(t40,1000.*FDP40,'k-',tdata40,1000.*FDPdata40,...
    'k^','linewidth',1.5,'Markerfacecolor',0*[1 1 1],'Markersize',8); 
hold off;legend('HMP', 'HMP data', 'FDP', 'FDP data','Location','Best')
% print -f9 -dtiff -r1200 'figure9.tiff'

%Figure 10
figure(10); hold on; set(gca,'Fontsize',14)
plot(t40,pH_cy40,'k:','linewidth',1.5);
hold off;

