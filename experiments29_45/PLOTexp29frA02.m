
%Data from Scopes (1973) experiment 29 with GPa fraction 0.002 or GPb
%fraction 0.998
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
% print -f3 -dtiff -r1200 'figure3.tiff'
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
% print -f4 -dtiff -r1200 'figure4.tiff'

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
% print -f5 -dtiff -r1200 'figure5.tiff'

%Figure 6
figure(6); clf; set(gca,'Fontsize',14)
plot(t,1000*G3P,'k--',t_data,1000*G3P_data,'ko','linewidth',1.5, ...
    'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); 
ylabel('Glycerol-3-Phosphate, Phosphoglycerates (mM)'); box on; hold on; 
plot(t,1000*PG,'k-',t_data,1000*PG_data,'k^','linewidth',1.5, ...
    'Markerfacecolor',0*[1 1 1],'Markersize',8); 
hold off;legend('G3P','G3P data','PG','PG data','Location','Best')
% print -f6 -dtiff -r1200 'figure6.tiff'

%Figure 10
figure(10); clf; set(gca,'Fontsize',14)
plot(t,pH_cy,'k-','linewidth',1.5); 
xlabel('Time (min)'); ylabel('pH'); box on; %hold on;
%%

figure(11); clf; set(gca,'Fontsize',14)
plot(t,CKprtflux,'k-','linewidth',1.5); 
xlabel('Time (min)'); ylabel('Proton Consumption Fluxes (M/min)'); box on; hold on;
plot(t,glycprtflux,'k:','linewidth',1.5); hold off;
set(gca,'YLim',[-0.04 0.02],'XLim',[0 1])
legend('CK proton flux','Other proton flux','Location','Best')
% Figure11 inset
axes('position',[0.5 0.2 0.35 0.35])
plot(t,CKprtflux,'k-','linewidth',1.5); 
xlabel('Time (min)'); ylabel('Proton Consumption Fluxes (M/min)'); 
box on; hold on;
plot(t,glycprtflux,'k:','linewidth',1.5); hold off;
set(gca,'YLim',[-4 0.1],'XLim',[0 0.05])
% print -f11 -dtiff -r1200 'figure11.tiff'
fig11 = figure(11)
resp11 = fig2plotly(fig11, 'filename', 'fig11', 'strip', false)
%%

% Figure 17
figure(18); clf; set(gca,'Fontsize',14)
plot(t,1000*protonload,'k--','linewidth',1.5,'Markerfacecolor', ...
    [1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('Proton load, Lactate production (mM)'); 
box on; hold on; 
plot(t,1000*LAC,'k-','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1], ...
    'Markersize',8); 
hold off;set(gca,'YLim',[0 15]);
legend('Proton load','Lactate production','Location','Best')
% print -f18 -dtiff -r1200 'figure17a.tiff'
fig17a = figure(18)
resp17a = fig2plotly(fig17a, 'filename', 'fig17a', 'strip', false)
