%%%
%Figure 10
figure(10); hold on; set(gca,'Fontsize',14)
plot(tckk,pH_cyckk,'k--','linewidth',1.5);
hold off;set(gca,'XLim',[0 2],'YLim',[7.1 8]);
legend('GPb 99.8%', 'GPb 60%', 'No CK','Location','NorthEast')
% Create textarrow
annotation1 = annotation(...
  10,'textarrow',...
  [0.5166 0.5166],[0.863 0.8208],...
  'String',{'pHstat starts'});
% % print -f10 -dtiff -r1200 'figure10arrow.tiff'
fig10 = figure(10);
resp10 = fig2plotly(fig10, 'filename', 'fig10', 'strip', false);
%

% Figure 16 panel a
figure(16); clf; set(gca,'Fontsize',14)
semilogy(tckk,1000.*PCrckk,'k--','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('PCr, Lactate (mM)'); box on; hold on; 
semilogy(tckk,1000.*LACckk,'k-','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8); 
hold off;set(gca,'XLim',[0 5],'YLim',[1e-7 20]);
legend('PCr','Lactate','Location',[0.6776 0.7166 0.1964 0.1095])
% print -f16 -dtiff -r1200 'figure16a.tiff'
fig16a = figure(16);
resp16a = fig2plotly(fig16a, 'filename', 'fig16a', 'strip', false);
%%

% Figure 16 panel b
figure(17); clf; set(gca,'Fontsize',14)
semilogy(tckk,1000.*ATPckk,'k--','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('ATP, ADP, AMP (mM)'); box on; hold on; 
semilogy(tckk,1000.*ADPckk,'k-','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1],'Markersize',8);
semilogy(tckk,1000.*AMPckk,'k-.','linewidth',1.5,'Markerfacecolor',[1 1 1],'Markersize',8); 
hold off; set(gca,'XLim',[0 5],'YLim',[0 10])
legend('ATP','ADP','AMP','Location','Best')
% % print -f17 -dtiff -r1200 'figure16b.tiff'
fig166 = figure(17);
resp166 = fig2plotly(fig166, 'filename', 'fig16b', 'strip', false);

% Figure 17 panel b
figure(19); clf; set(gca,'Fontsize',14)
plot(tckk,1000*protonloadckk,'k--','linewidth',1.5,'Markerfacecolor', ...
    [1 1 1],'Markersize',8); 
xlabel('Time (min)'); ylabel('Proton load, Lactate production (mM)'); 
box on; hold on; 
plot(tckk,1000*LACckk,'k-','linewidth',1.5,'Markerfacecolor',0.75*[1 1 1], ...
    'Markersize',8); 
hold off;set(gca,'YLim',[-2 13]);
legend('Proton load','Lactate production','Location','Best')
% % print -f19 -dtiff -r1200 'figure17b.tiff'
fig177 = figure(19);
resp177 = fig2plotly(fig177, 'filename', 'fig17b', 'strip', false);
